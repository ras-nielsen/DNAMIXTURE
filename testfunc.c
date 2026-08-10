#include<stdio.h>
#include<stdlib.h>
#include<ctype.h>
#include<math.h>
#include<string.h>
#include"neldermead.h"
#include"json_parser.h"
// Nucleotide encoding: A=0, C=1, G=2, T=3
#define NUC_A 0
#define NUC_C 1
#define NUC_G 2
#define NUC_T 3

// Global variables used by json_parser module
MultiSNPData *global_snp_data = NULL;
Options global_opts;

// Likelihood ratio types
typedef enum {
    LR_L1,           // Suspect vs No-suspect
    LR_L2,           // Suspect vs Sibling
    LR_L3,           // Suspect vs Parent
    LR_L4            // Suspect vs Cousin
} LikelihoodRatioType;

// Contaminant model types
typedef enum {
    CONTAM_POPULATION,       // General population contaminant
    CONTAM_SINGLE_INDIV      // Single individual contaminant
} ContaminantModel;

// Hypothesis type - determines which likelihood function to use
typedef enum {
    HYPO_SUSPECT,
    HYPO_NOSUSPECT,
    HYPO_SIBLING,
    HYPO_PARENT,
    HYPO_COUSIN
} HypothesisType;

// Analysis parameters - passed to likelihood functions
typedef struct {
    double error_adj;      // Error adjustment parameter (Equation 11)
    int cousin_k;          // Cousin degree for L4 calculations
} AnalysisParams;

// Optimization context - single static variable for amoeba callbacks
typedef struct {
    MultiSNPData *data;       // Pointer to SNP data
    AnalysisParams params;    // Analysis parameters (error_adj, cousin_k)
    double f1_init;           // Initial f1 value
    double f2_init;           // Initial f2 value
    HypothesisType hypo_type; // Which hypothesis to optimize
    ContaminantModel model;   // Population or single individual contaminant
} OptimizationContext;

// Static optimization context (required for amoeba callback compatibility)
static OptimizationContext *opt_context = NULL;

// ============================================================================
// MEMORY MANAGEMENT FUNCTIONS
// ============================================================================

// Allocate multi-SNP data structure
MultiSNPData* allocate_multisnp_data(int numsnps)
{
        MultiSNPData *data = (MultiSNPData*)malloc(sizeof(MultiSNPData));
        if (!data) {
                fprintf(stderr, "Error: Failed to allocate MultiSNPData\n");
                exit(1);
        }
        data->numsnps = numsnps;
        data->snps = (SNPData*)calloc(numsnps, sizeof(SNPData));
        if (!data->snps) {
                fprintf(stderr, "Error: Failed to allocate SNP array\n");
                exit(1);
        }
        // Initialize genotypes to -1 (missing) so we can detect unset values
        for (int i = 0; i < numsnps; i++) {
                data->snps[i].suspect_genotype[0] = -1;
                data->snps[i].suspect_genotype[1] = -1;
                data->snps[i].victim_genotype[0] = -1;
                data->snps[i].victim_genotype[1] = -1;
        }
        return data;
}

// Allocate reads array for a specific SNP
void allocate_snp_reads(SNPData *snp, int numreads)
{
        snp->numreads = numreads;
        snp->reads = (Read*)malloc(numreads * sizeof(Read));
        if (!snp->reads) {
                fprintf(stderr, "Error: Failed to allocate reads array\n");
                exit(1);
        }
}

// Free multi-SNP data structure
void free_multisnp_data(MultiSNPData *data)
{
        if (!data) return;
        if (data->snps) {
                for (int i = 0; i < data->numsnps; i++) {
                        if (data->snps[i].reads) {
                                free(data->snps[i].reads);
                        }
                }
                free(data->snps);
        }
        free(data);
}

// Set a single read at a SNP
void set_snp_read(SNPData *snp, int read_idx, int nucleotide, double quality)
{
        if (read_idx < 0 || read_idx >= snp->numreads) {
                fprintf(stderr, "Error: Read index out of bounds\n");
                exit(1);
        }
        snp->reads[read_idx].nucleotide = nucleotide;
        snp->reads[read_idx].quality = quality;
}

// Set allele frequencies for a SNP
void set_snp_af(SNPData *snp, double af[4])
{
        for (int i = 0; i < 4; i++) {
                snp->AF[i] = af[i];
        }
}

// Set genotypes for a SNP
void set_snp_genotypes(SNPData *snp, int victim_g1, int victim_g2, int suspect_g1, int suspect_g2)
{
        snp->victim_genotype[0] = victim_g1;
        snp->victim_genotype[1] = victim_g2;
        snp->suspect_genotype[0] = suspect_g1;
        snp->suspect_genotype[1] = suspect_g2;
}

int static z_rndu=137;
void SetSeed(int seed)
{
   z_rndu = 170*(seed%178) + 137;
}

double uniform()
{
/*

	U(0,1): AS 183: Appl. Stat. 31:188-190
   Wichmann BA & Hill ID.  1982.  An efficient and portable
   pseudo-random number generator.  Appl. Stat. 31:188-190

   x, y, z are any numbers in the range 1-30000.  Integer operation up
   to 30323 required.
   
   Suggested to me by Z. Yang who also provided me with
   the source code used here.
*/
   static int x_rndu=11, y_rndu=23;
   double r;

   x_rndu = 171*(x_rndu%177) -  2*(x_rndu/177);
   y_rndu = 172*(y_rndu%176) - 35*(y_rndu/176);
   z_rndu = 170*(z_rndu%178) - 63*(z_rndu/178);
   if (x_rndu<0) x_rndu+=30269;
   if (y_rndu<0) y_rndu+=30307;
   if (z_rndu<0) z_rndu+=30323;
   r = x_rndu/30269.0 + y_rndu/30307.0 + z_rndu/30323.0;
   return (r-(int)r);
}

// Global flag for population mismatch warning
static int global_population_mismatch = 0;
static double global_mismatch_zscore = 0.0;

// ============================================================================
// POPULATION-GENOTYPE MATCH TEST
// ============================================================================
// Tests if suspect genotype matches reference population allele frequencies
// Compares observed statistic to simulated distribution under HWE

void test_population_match(void)
{
	int i, j, s, idx;
	int num_valid_snps = 0;
	double t_observed = 0.0;
	double t_simulated[1000];
	double mean_t, var_t, sd_t, z_score;

	// First pass: count valid SNPs
        if (!global_random_snps || global_random_snps->count == 0) {
        fprintf(stderr, "Warning: No random SNPs for population match test\n");
        return;
        }

        num_valid_snps = global_random_snps->count;

	// Allocate arrays to store precomputed values for valid SNPs
	int *major_allele = (int *)malloc(num_valid_snps * sizeof(int));
	double *major_freq = (double *)malloc(num_valid_snps * sizeof(double));
	double (*cumsum)[4] = malloc(num_valid_snps * sizeof(*cumsum));

	if (!major_allele || !major_freq || !cumsum ) {
		fprintf(stderr, "Error: Failed to allocate memory for population match test\n");
		exit(1);
	}

	// Precompute major allele, frequency, and cumulative sums for each valid SNP
	idx = 0;
	for (s = 0; s < global_random_snps->count; s++) {
                idx = s;
		RandomSNP *rsnp = &global_random_snps->snps[s];

		// Find major allele (highest frequency)
		major_allele[idx] = 0;
		double max_freq = rsnp->AF[0];
		for (i = 1; i < 4; i++) {
			if (rsnp->AF[i] > max_freq) {
				max_freq = rsnp->AF[i];
				major_allele[idx] = i;
			}
		}
		major_freq[idx] = max_freq;

		// Precompute cumulative frequencies for sampling
		cumsum[idx][0] = rsnp->AF[0];
		for (i = 1; i < 4; i++) {
			cumsum[idx][i] = cumsum[idx][i-1] + rsnp->AF[i];
		}
	}

	// Calculate observed t for the suspect
	for (idx = 0; idx < num_valid_snps; idx++) {
		RandomSNP *rsnp = &global_random_snps->snps[idx];
		double f_i = major_freq[idx];

		// Code suspect genotype relative to major allele
		double g_i;
		int allele1 = rsnp->suspect_gt[0];
		int allele2 = rsnp->suspect_gt[1];
		if (allele1 == major_allele[idx] && allele2 == major_allele[idx]) {
			g_i = 1.0;
		} else if (allele1 == major_allele[idx] || allele2 == major_allele[idx]) {
			g_i = 0.5;
		} else {
			g_i = 0.0;
		}

		t_observed += (f_i - g_i) * (f_i - g_i);
	}

	// Simulate 1000 random genomes under HWE
	for (j = 0; j < 1000; j++) {
		t_simulated[j] = 0.0;

		for (idx = 0; idx < num_valid_snps; idx++) {
			double f_i = major_freq[idx];

			// Sample first allele
			double u1 = uniform();
			int sim_allele1 = 0;
			for (i = 0; i < 4; i++) {
				if (u1 < cumsum[idx][i]) {
					sim_allele1 = i;
					break;
				}
			}

			// Sample second allele
			double u2 = uniform();
			int sim_allele2 = 0;
			for (i = 0; i < 4; i++) {
				if (u2 < cumsum[idx][i]) {
					sim_allele2 = i;
					break;
				}
			}

			// Code simulated genotype relative to major allele
			double g_sim;
			if (sim_allele1 == major_allele[idx] && sim_allele2 == major_allele[idx]) {
				g_sim = 1.0;
			} else if (sim_allele1 == major_allele[idx] || sim_allele2 == major_allele[idx]) {
				g_sim = 0.5;
			} else {
				g_sim = 0.0;
			}

			t_simulated[j] += (f_i - g_sim) * (f_i - g_sim);
		}
	}

	// Free precomputed arrays
	free(major_allele);
	free(major_freq);
	free(cumsum);

	// Calculate mean and SD of simulated values
	mean_t = 0.0;
	for (j = 0; j < 1000; j++) {
		mean_t += t_simulated[j];
	}
	mean_t /= 1000.0;

	var_t = 0.0;
	for (j = 0; j < 1000; j++) {
		var_t += (t_simulated[j] - mean_t) * (t_simulated[j] - mean_t);
	}
	var_t /= 999.0;  // Sample variance
	sd_t = sqrt(var_t);

	// Calculate Z-score
	z_score = (t_observed - mean_t) / sd_t;
	global_mismatch_zscore = z_score;

	// Report results
	fprintf(stderr, "\n=== Population Match Test ===\n");
	fprintf(stderr, "Number of SNPs tested: %d\n", num_valid_snps);
	fprintf(stderr, "Observed t statistic: %.4f\n", t_observed);
	fprintf(stderr, "Simulated mean: %.4f, SD: %.4f\n", mean_t, sd_t);
	fprintf(stderr, "Z-score: %.2f standard deviations from mean\n", z_score);

	if (fabs(z_score) > 2.0) {
		global_population_mismatch = 1;
		fprintf(stderr, "\nMISMATCH BETWEEN SUSPECT AND REFERENCE POPULATION. RESULTS MAY NOT BE RELIABLE.\n");
	} else {
		global_population_mismatch = 0;
		fprintf(stderr, "Suspect genotype consistent with reference population.\n");
	}
	fprintf(stderr, "==============================\n");
}

// ============================================================================
// QUALITY SCORE CONVERSION FUNCTIONS
// ============================================================================

// Convert phred quality score to error probability
// quality: phred-scaled quality score (typically 0-60)
// Uses standard Phred+33 encoding: P_error = 10^(-Q/10)
double quality_to_error_prob(double quality)
{
        // Validate reasonable range
        if (quality < 0) {
                fprintf(stderr, "Error: Quality score %.2f is negative. Valid range is 0-60.\n", quality);
                fprintf(stderr, "This may indicate incorrect quality score encoding or data corruption.\n");
                exit(1);
        }
        if (quality > 60) {
                fprintf(stderr, "Error: Quality score %.2f exceeds maximum of 60. Valid range is 0-60.\n", quality);
                fprintf(stderr, "This may indicate incorrect quality score encoding or data corruption.\n");
                exit(1);
        }

        return pow(10.0, -quality / 10.0);
}

// Probability of observing 'observed' nucleotide given true nucleotide and quality
// observed: observed nucleotide (0-3)
// true_nt: true nucleotide (0-3)
// quality: phred quality score
// error_adj: error adjustment parameter for Equation 11
// Returns: p(observed | true_nt, quality)
// If error_adj > 0, applies Equation 11: p'(x|nt) = (p(x|nt) + e) / (1 + 4e)
// Note: denominator is 1+4e because e is added to all 4 nucleotides for proper normalization
double read_likelihood_from_quality(int observed, int true_nt, double quality, double error_adj)
{
        double p_error = quality_to_error_prob(quality);
        double p_read;

        if (observed == true_nt) {
                // Correct observation
                p_read = 1.0 - p_error;
        } else {
                // Error - uniform over 3 other bases
                p_read = p_error / 3.0;
        }

        // Apply Equation 11 if error_adj > 0
        // Corrected denominator: 1 + 4e (since e is added to all 4 nucleotides)
        if (error_adj > 0.0) {
                p_read = (p_read + error_adj) / (1.0 + 4.0 * error_adj);
        }

        return p_read;
}

// Read likelihood given diploid genotype and quality score
// read: the Read structure (nucleotide + quality)
// genotype: diploid genotype [allele1, allele2]
// error_adj: error adjustment parameter
// Returns: p(read | genotype) = average of p(read|allele1) and p(read|allele2)
double read_likelihood_given_genotype_q(Read read, int genotype[2], double error_adj)
{
        double p1 = read_likelihood_from_quality(read.nucleotide, genotype[0], read.quality, error_adj);
        double p2 = read_likelihood_from_quality(read.nucleotide, genotype[1], read.quality, error_adj);
        return (p1 + p2) / 2.0;
}

// Read likelihood given single allele and quality score
double read_likelihood_given_allele_q(Read read, int allele, double error_adj)
{
        return read_likelihood_from_quality(read.nucleotide, allele, read.quality, error_adj);
}


// Log-sum-exp trick for numerical stability
// Computes log(sum(exp(log_values[i]))) without overflow/underflow
double log_sum_exp(double* log_values, int n)
{
        if (n == 0) return -INFINITY;
        if (n == 1) return log_values[0];

        // Find maximum value
        double max_val = log_values[0];
        int i;
        for (i = 1; i < n; i++) {
                if (log_values[i] > max_val)
                        max_val = log_values[i];
        }

        // If max is -infinity, all values are -infinity
        if (isinf(max_val) && max_val < 0)
                return -INFINITY;

        // Compute sum of exp(log_values[i] - max_val)
        double sum = 0.0;
        for (i = 0; i < n; i++)
                sum += exp(log_values[i] - max_val);

        return max_val + log(sum);
}

// ============================================================================
// MULTI-SNP LIKELIHOOD FUNCTIONS
// ============================================================================
// These functions work with MultiSNPData and sum log-likelihoods across SNPs

// PARAMETER TRANSFORMATIONS for constrained optimization
// Softmax transformation: ensures f1>=0, f2>=0, f1+f2<=1
void transform_params_2d(const double params[], double *f1_out, double *f2_out)
{
        double a = params[0];
        double b = params[1];
        double denom = 1.0 + exp(a) + exp(b);
        *f1_out = exp(a) / denom;
        *f2_out = exp(b) / denom;
}

// Logit transformation for 1D optimization: ensures 0 <= f2 <= 1
void transform_params_1d(const double params[], double *f2_out)
{
        double b = params[0];
        *f2_out = exp(b) / (1.0 + exp(b));
}

// Inverse transformation (for initialization)
void inverse_transform_2d(double f1, double f2, double params[])
{
        // Given f1, f2, find a, b such that:
        // f1 = exp(a) / (1 + exp(a) + exp(b))
        // f2 = exp(b) / (1 + exp(a) + exp(b))
        // Solution: a = log(f1/(1-f1-f2)), b = log(f2/(1-f1-f2))
        double f3 = 1.0 - f1 - f2;
        if (f3 <= 0) f3 = 0.01;  // Safety check
        if (f1 <= 0) f1 = 0.01;
        if (f2 <= 0) f2 = 0.01;
        params[0] = log(f1 / f3);
        params[1] = log(f2 / f3);
}

// ============================================================================
// CORE MULTI-SNP LIKELIHOOD FUNCTIONS
// ============================================================================

// Contaminant likelihood at a single SNP: p(read | contaminant, AF)
// Sums over all possible alleles weighted by allele frequency
double log_read_given_contaminant(Read read, double AF[4], double error_adj)
{
        double log_values[4];
        int i;
        for (i = 0; i < 4; i++) {
                log_values[i] = log(read_likelihood_given_allele_q(read, i, error_adj)) + log(AF[i]);
        }
        return log_sum_exp(log_values, 4);
}

// Suspect likelihood for all reads at a single SNP
// p(reads_at_snp | suspect_genotype, victim_genotype, f1, f2, AF)
double log_suspect_likelihood_single_snp(SNPData *snp, double f1_val, double f2_val, double error_adj)
{
        double log_lik = 0.0;
        int i;

        for (i = 0; i < snp->numreads; i++) {
                Read read = snp->reads[i];
                double log_values[3];

                // Component 1: p(read|suspect) * f1
                log_values[0] = log(read_likelihood_given_genotype_q(read, snp->suspect_genotype, error_adj)) + log(f1_val);

                // Component 2: p(read|victim) * f2
                log_values[1] = log(read_likelihood_given_genotype_q(read, snp->victim_genotype, error_adj)) + log(f2_val);

                // Component 3: p(read|contaminant) * (1-f1-f2)
                log_values[2] = log_read_given_contaminant(read, snp->AF, error_adj) + log(1.0 - f1_val - f2_val);

                log_lik += log_sum_exp(log_values, 3);
        }

        return log_lik;
}

// No-suspect likelihood for all reads at a single SNP (f1=0)
double log_nosuspect_likelihood_single_snp(SNPData *snp, double f2_val, double error_adj)
{
        double log_lik = 0.0;
        int i;

        for (i = 0; i < snp->numreads; i++) {
                Read read = snp->reads[i];
                double log_values[2];

                // Component 1: p(read|victim) * f2
                log_values[0] = log(read_likelihood_given_genotype_q(read, snp->victim_genotype, error_adj)) + log(f2_val);

                // Component 2: p(read|contaminant) * (1-f2)
                log_values[1] = log_read_given_contaminant(read, snp->AF, error_adj) + log(1.0 - f2_val);

                log_lik += log_sum_exp(log_values, 2);
        }

        return log_lik;
}

// Parent likelihood for all reads at a single SNP
// One allele from suspect, one from population
double log_parent_likelihood_single_snp(SNPData *snp, double f1_val, double f2_val, double error_adj)
{
        double log_sum_products[8];  // 2 suspect alleles * 4 population alleles
        int idx = 0;
        int k, v, i;

        for (k = 0; k < 2; k++) {  // Loop over suspect's two alleles
                int suspect_allele = snp->suspect_genotype[k];

                for (v = 0; v < 4; v++) {  // Loop over population alleles
                        double log_product = 0.0;

                        for (i = 0; i < snp->numreads; i++) {
                                Read read = snp->reads[i];
                                double log_values[3];

                                // Parent genotype = (suspect_allele, v)
                                // p(read | parent genotype) = average of p(read|allele1) and p(read|allele2)
                                double p1 = read_likelihood_from_quality(read.nucleotide, suspect_allele, read.quality, error_adj);
                                double p2 = read_likelihood_from_quality(read.nucleotide, v, read.quality, error_adj);
                                double p_read_given_parent = (p1 + p2) / 2.0;

                                log_values[0] = log(p_read_given_parent) + log(f1_val);
                                log_values[1] = log(read_likelihood_given_genotype_q(read, snp->victim_genotype, error_adj)) + log(f2_val);
                                log_values[2] = log_read_given_contaminant(read, snp->AF, error_adj) + log(1.0 - f1_val - f2_val);

                                log_product += log_sum_exp(log_values, 3);
                        }

                        // Weight by allele frequency
                        log_sum_products[idx++] = log_product + log(snp->AF[v]);
                }
        }

        // Average over suspect's two alleles (factor of 1/2)
        return log_sum_exp(log_sum_products, 8) + log(0.5);
}

// Sibling likelihood for all reads at a single SNP
// Simplified by reusing functions for each IBD state
double log_sibling_likelihood_single_snp(SNPData *snp, double f1_val, double f2_val, double error_adj)
{
        double log_components[3];

        // IBD=1 (prob 0.50): Sibling shares one allele - same as parent likelihood
        // But need to add the 0.5 IBD probability weight
        log_components[0] = log_parent_likelihood_single_snp(snp, f1_val, f2_val, error_adj) + log(0.5);

        // IBD=2 (prob 0.25): Sibling identical to suspect
        log_components[1] = log_suspect_likelihood_single_snp(snp, f1_val, f2_val, error_adj) + log(0.25);

        // IBD=0 (prob 0.25): Sibling unrelated - same as no-suspect likelihood
        log_components[2] = log_nosuspect_likelihood_single_snp(snp, f2_val, error_adj) + log(0.25);

        return log_sum_exp(log_components, 3);
}

// Kth cousin likelihood for all reads at a single SNP
// For kth cousins:
//   IBD=0 with probability 1-(1/2)^(2k): unrelated
//   IBD=1 with probability (1/2)^(2k): share one allele
// k=1: first cousins, k=2: second cousins, etc.
// cousin_k: degree of cousinship
double log_cousin_likelihood_single_snp(SNPData *snp, double f1_val, double f2_val, double error_adj, int cousin_k)
{
        double log_components[2];
        double prob_ibd1 = pow(0.5, 2*cousin_k);  // Probability of IBD=1
        double prob_ibd0 = 1.0 - prob_ibd1;       // Probability of IBD=0

        // IBD=1: Cousin shares one allele - same as parent likelihood
        log_components[0] = log_parent_likelihood_single_snp(snp, f1_val, f2_val, error_adj) + log(prob_ibd1);

        // IBD=0: Cousin unrelated - same as no-suspect likelihood
        log_components[1] = log_nosuspect_likelihood_single_snp(snp, f2_val, error_adj) + log(prob_ibd0);

        return log_sum_exp(log_components, 2);
}

// ============================================================================
// SINGLE INDIVIDUAL CONTAMINANT MODEL (Equation 4 from paper)
// ============================================================================
// These functions assume contamination comes from a single individual
// rather than the general population. We sum over all possible genotypes
// of the contaminant individual.

// Helper: compute HWE genotype probability in log-space (UNORDERED)
double log_genotype_probability_hwe(int allele1, int allele2, double AF[4])
{
        if (allele1 == allele2) {
                // Homozygous: p(AA) = q_A^2
                return 2.0 * log(AF[allele1]);
        } else {
                // Heterozygous: p(AC) = 2*q_A*q_C
                return log(2.0) + log(AF[allele1]) + log(AF[allele2]);
        }
}

// Suspect likelihood with SINGLE INDIVIDUAL contaminant at a single SNP
// Implements Equation 4: p(x|f1,f2) = Σ_{g∈G} p(genotype=g) × Π_{reads} [...]
double log_suspect_likelihood_single_indiv_contam_single_snp(SNPData *snp, double f1_val, double f2_val, double error_adj)
{
        double log_genotype_terms[10];  // 10 possible genotypes (4 homozygous + 6 heterozygous)
        int idx = 0;
        int g1, g2;

        // Sum over all possible contaminant genotypes
        for (g1 = 0; g1 < 4; g1++) {
                for (g2 = g1; g2 < 4; g2++) {  // g2 >= g1 to avoid double-counting
                        int contaminant_genotype[2] = {g1, g2};

                        // Start with log probability of this genotype
                        double log_term = log_genotype_probability_hwe(g1, g2, snp->AF);

                        // Multiply across all reads
                        int i;
                        for (i = 0; i < snp->numreads; i++) {
                                Read read = snp->reads[i];
                                double log_values[3];

                                // Component 1: p(read|suspect) * f1
                                log_values[0] = log(read_likelihood_given_genotype_q(read, snp->suspect_genotype, error_adj)) + log(f1_val);

                                // Component 2: p(read|contaminant_genotype) * (1-f1-f2)
                                log_values[1] = log(read_likelihood_given_genotype_q(read, contaminant_genotype, error_adj)) + log(1.0 - f1_val - f2_val);

                                // Component 3: p(read|victim) * f2
                                log_values[2] = log(read_likelihood_given_genotype_q(read, snp->victim_genotype, error_adj)) + log(f2_val);

                                log_term += log_sum_exp(log_values, 3);
                        }

                        log_genotype_terms[idx++] = log_term;
                }
        }

        return log_sum_exp(log_genotype_terms, 10);
}

// No-suspect likelihood with SINGLE INDIVIDUAL contaminant (f1=0)
double log_nosuspect_likelihood_single_indiv_contam_single_snp(SNPData *snp, double f2_val, double error_adj)
{
        double log_genotype_terms[10];
        int idx = 0;
        int g1, g2;

        for (g1 = 0; g1 < 4; g1++) {
                for (g2 = g1; g2 < 4; g2++) {
                        int contaminant_genotype[2] = {g1, g2};

                        double log_term = log_genotype_probability_hwe(g1, g2, snp->AF);

                        int i;
                        for (i = 0; i < snp->numreads; i++) {
                                Read read = snp->reads[i];
                                double log_values[2];

                                // Component 1: p(read|contaminant_genotype) * (1-f2)
                                log_values[0] = log(read_likelihood_given_genotype_q(read, contaminant_genotype, error_adj)) + log(1.0 - f2_val);

                                // Component 2: p(read|victim) * f2
                                log_values[1] = log(read_likelihood_given_genotype_q(read, snp->victim_genotype, error_adj)) + log(f2_val);

                                log_term += log_sum_exp(log_values, 2);
                        }

                        log_genotype_terms[idx++] = log_term;
                }
        }

        return log_sum_exp(log_genotype_terms, 10);
}

// Parent likelihood with SINGLE INDIVIDUAL contaminant
double log_parent_likelihood_single_indiv_contam_single_snp(SNPData *snp, double f1_val, double f2_val, double error_adj)
{
        double log_parent_genotype_terms[8];  // 2 suspect alleles * 4 population alleles
        int parent_idx = 0;
        int k, v;

        // Sum over possible parent genotypes
        for (k = 0; k < 2; k++) {
                int suspect_allele = snp->suspect_genotype[k];

                for (v = 0; v < 4; v++) {
                        int parent_genotype[2] = {suspect_allele, v};
                        double log_parent_prob = log(0.5) + log(snp->AF[v]);

                        // Sum over contaminant genotypes IN PROBABILITY SPACE
                        double prob_sum_contaminants = 0.0;
                        int g1, g2;

                        for (g1 = 0; g1 < 4; g1++) {
                                for (g2 = g1; g2 < 4; g2++) {
                                        int contaminant_genotype[2] = {g1, g2};

                                        double log_term = log_genotype_probability_hwe(g1, g2, snp->AF);

                                        int i;
                                        for (i = 0; i < snp->numreads; i++) {
                                                Read read = snp->reads[i];
                                                double log_values[3];

                                                log_values[0] = log(read_likelihood_given_genotype_q(read, parent_genotype, error_adj)) + log(f1_val);
                                                log_values[1] = log(read_likelihood_given_genotype_q(read, contaminant_genotype, error_adj)) + log(1.0 - f1_val - f2_val);
                                                log_values[2] = log(read_likelihood_given_genotype_q(read, snp->victim_genotype, error_adj)) + log(f2_val);

                                                log_term += log_sum_exp(log_values, 3);
                                        }

                                        // Add this contaminant's contribution in probability space
                                        prob_sum_contaminants += exp(log_term);
                                }
                        }

                        // Multiply parent probability by contaminant sum, then convert to log
                        double log_parent_term = log(exp(log_parent_prob) * prob_sum_contaminants);
                        log_parent_genotype_terms[parent_idx++] = log_parent_term;
                }
        }

        return log_sum_exp(log_parent_genotype_terms, 8);
}

// Sibling likelihood with SINGLE INDIVIDUAL contaminant
// Simplified by reusing verified functions for each IBD state
double log_sibling_likelihood_single_indiv_contam_single_snp(SNPData *snp, double f1_val, double f2_val, double error_adj)
{
        double log_components[3];

        // IBD=2 (prob 0.25): Sibling identical to suspect
        log_components[0] = log_suspect_likelihood_single_indiv_contam_single_snp(snp, f1_val, f2_val, error_adj) + log(0.25);

        // IBD=1 (prob 0.50): Sibling shares one allele - same as parent likelihood
        log_components[1] = log_parent_likelihood_single_indiv_contam_single_snp(snp, f1_val, f2_val, error_adj) + log(0.5);

        // IBD=0 (prob 0.25): Sibling unrelated - same as no-suspect likelihood
        log_components[2] = log_nosuspect_likelihood_single_indiv_contam_single_snp(snp, f2_val, error_adj) + log(0.25);

        return log_sum_exp(log_components, 3);
}

// Kth cousin likelihood with SINGLE INDIVIDUAL contaminant
// For kth cousins:
//   IBD=0 with probability 1-(1/2)^(2k+2): unrelated
//   IBD=1 with probability (1/2)^(2k+1): share one allele
// k=1: first cousins, k=2: second cousins, etc.
// cousin_k: degree of cousinship
double log_cousin_likelihood_single_indiv_contam_single_snp(SNPData *snp, double f1_val, double f2_val, double error_adj, int cousin_k)
{
        double log_components[2];
        double prob_ibd1 = pow(0.5, 2*cousin_k + 1);  // Probability of IBD=1
        double prob_ibd0 = 1.0 - prob_ibd1;           // Probability of IBD=0

        // IBD=1: Cousin shares one allele - same as parent likelihood
        log_components[0] = log_parent_likelihood_single_indiv_contam_single_snp(snp, f1_val, f2_val, error_adj) + log(prob_ibd1);

        // IBD=0: Cousin unrelated - same as no-suspect likelihood
        log_components[1] = log_nosuspect_likelihood_single_indiv_contam_single_snp(snp, f2_val, error_adj) + log(prob_ibd0);

        return log_sum_exp(log_components, 2);
}

// ============================================================================
// MULTI-SNP LIKELIHOOD FUNCTIONS (sum across all SNPs)
// ============================================================================

// Function pointer types for single-SNP likelihood functions
typedef double (*likelihood_func_3param_t)(SNPData*, double, double, double);
typedef double (*likelihood_func_2param_t)(SNPData*, double, double);
typedef double (*likelihood_func_4param_t)(SNPData*, double, double, double, int);

// Generic function to sum log-likelihood across all SNPs (for functions with f1, f2, error_adj)
double sum_log_likelihood_3param(MultiSNPData *data, double f1_val, double f2_val, double error_adj,
                                  likelihood_func_3param_t single_snp_func)
{
        double total_log_lik = 0.0;
        int s;
        for (s = 0; s < data->numsnps; s++) {
                total_log_lik += single_snp_func(&data->snps[s], f1_val, f2_val, error_adj);
        }
        return total_log_lik;
}

// Generic function to sum log-likelihood across all SNPs (for no-suspect functions with f2, error_adj)
double sum_log_likelihood_2param(MultiSNPData *data, double f2_val, double error_adj,
                                  likelihood_func_2param_t single_snp_func)
{
        double total_log_lik = 0.0;
        int s;
        for (s = 0; s < data->numsnps; s++) {
                total_log_lik += single_snp_func(&data->snps[s], f2_val, error_adj);
        }
        return total_log_lik;
}

// Generic function to sum log-likelihood across all SNPs (for cousin functions with f1, f2, error_adj, cousin_k)
double sum_log_likelihood_4param(MultiSNPData *data, double f1_val, double f2_val, double error_adj, int cousin_k,
                                  likelihood_func_4param_t single_snp_func)
{
        double total_log_lik = 0.0;
        int s;
        for (s = 0; s < data->numsnps; s++) {
                total_log_lik += single_snp_func(&data->snps[s], f1_val, f2_val, error_adj, cousin_k);
        }
        return total_log_lik;
}

// ============================================================================
// UNIFIED OBJECTIVE FUNCTIONS FOR OPTIMIZATION
// ============================================================================
// These functions use the static opt_context to access data and parameters
// and select the appropriate likelihood function based on hypothesis type and model

// Unified objective function for 2D optimization (f1, f2)
double objective_2d(const double params[])
{
        double f1_opt, f2_opt;
        transform_params_2d(params, &f1_opt, &f2_opt);
        double loglik;

        // Select appropriate likelihood function based on hypothesis type and model
        if (opt_context->model == CONTAM_POPULATION) {
                switch (opt_context->hypo_type) {
                        case HYPO_SUSPECT:
                                loglik = sum_log_likelihood_3param(opt_context->data, f1_opt, f2_opt, opt_context->params.error_adj,
                                                                   log_suspect_likelihood_single_snp);
                                break;
                        case HYPO_SIBLING:
                                loglik = sum_log_likelihood_3param(opt_context->data, f1_opt, f2_opt, opt_context->params.error_adj,
                                                                   log_sibling_likelihood_single_snp);
                                break;
                        case HYPO_PARENT:
                                loglik = sum_log_likelihood_3param(opt_context->data, f1_opt, f2_opt, opt_context->params.error_adj,
                                                                   log_parent_likelihood_single_snp);
                                break;
                        case HYPO_COUSIN:
                                loglik = sum_log_likelihood_4param(opt_context->data, f1_opt, f2_opt, opt_context->params.error_adj,
                                                                   opt_context->params.cousin_k, log_cousin_likelihood_single_snp);
                                break;
                        default:
                                fprintf(stderr, "Error: Invalid hypothesis type for 2D optimization\n");
                                exit(1);
                }
        } else {  // CONTAM_SINGLE_INDIV
                switch (opt_context->hypo_type) {
                        case HYPO_SUSPECT:
                                loglik = sum_log_likelihood_3param(opt_context->data, f1_opt, f2_opt, opt_context->params.error_adj,
                                                                   log_suspect_likelihood_single_indiv_contam_single_snp);
                                break;
                        case HYPO_SIBLING:
                                loglik = sum_log_likelihood_3param(opt_context->data, f1_opt, f2_opt, opt_context->params.error_adj,
                                                                   log_sibling_likelihood_single_indiv_contam_single_snp);
                                break;
                        case HYPO_PARENT:
                                loglik = sum_log_likelihood_3param(opt_context->data, f1_opt, f2_opt, opt_context->params.error_adj,
                                                                   log_parent_likelihood_single_indiv_contam_single_snp);
                                break;
                        case HYPO_COUSIN:
                                loglik = sum_log_likelihood_4param(opt_context->data, f1_opt, f2_opt, opt_context->params.error_adj,
                                                                   opt_context->params.cousin_k, log_cousin_likelihood_single_indiv_contam_single_snp);
                                break;
                        default:
                                fprintf(stderr, "Error: Invalid hypothesis type for 2D optimization\n");
                                exit(1);
                }
        }

        return -loglik;  // Negative for minimization
}

// Unified objective function for 1D optimization (f2 only, for no-suspect)
double objective_1d(const double params[])
{
        double f2_opt;
        transform_params_1d(params, &f2_opt);
        double loglik;

        if (opt_context->model == CONTAM_POPULATION) {
                loglik = sum_log_likelihood_2param(opt_context->data, f2_opt, opt_context->params.error_adj,
                                                   log_nosuspect_likelihood_single_snp);
        } else {  // CONTAM_SINGLE_INDIV
                loglik = sum_log_likelihood_2param(opt_context->data, f2_opt, opt_context->params.error_adj,
                                                   log_nosuspect_likelihood_single_indiv_contam_single_snp);
        }

        return -loglik;
}

// ============================================================================
// GENERIC OPTIMIZATION HELPERS
// ============================================================================

// Generic 2D optimizer: optimize f1 and f2
// Uses opt_context to determine which hypothesis and model to optimize
double optimize_f1_f2(double *f1_hat, double *f2_hat)
{
        double x[2], fmin;

        // Nelder-Mead with parameter transformation
        // step 0.5, ftol 1e-10 (tight tolerance), max 5000 evaluations
        inverse_transform_2d(opt_context->f1_init, opt_context->f2_init, x);
        nelder_mead(x, 2, 0.5, 1e-10, 5000, objective_2d, &fmin);
        transform_params_2d(x, f1_hat, f2_hat);

        return -fmin;
}

// Generic 1D optimizer: optimize f2 only (for no-suspect)
// Uses opt_context to determine which model to optimize
double optimize_f2_only(double *f2_hat)
{
        double x[1], fmin;

        // Nelder-Mead with parameter transformation
        // step 0.5, ftol 1e-10 (tight tolerance), max 5000 evaluations
        x[0] = log(opt_context->f2_init / (1.0 - opt_context->f2_init + 1e-10));
        nelder_mead(x, 1, 0.5, 1e-10, 5000, objective_1d, &fmin);
        transform_params_1d(x, f2_hat);

        return -fmin;
}

// ============================================================================
// UNIFIED LIKELIHOOD RATIO FUNCTIONS
// ============================================================================
// NOTE: These functions assume opt_context is already set by the caller
// They use the unified objective functions and set hypothesis types as needed

// Generic likelihood ratio calculator
// Computes LR for suspect vs. alternative hypothesis
// Also returns individual log-likelihoods via output parameters
double compute_likelihood_ratio(HypothesisType alt_hypothesis,
                                 double *f1_hat_suspect, double *f2_hat_suspect,
                                 double *f1_hat_alt, double *f2_hat_alt,
                                 double *loglik_suspect_out, double *loglik_alt_out)
{
        double loglik_suspect, loglik_alt;

        // Optimize suspect hypothesis
        opt_context->hypo_type = HYPO_SUSPECT;
        loglik_suspect = optimize_f1_f2(f1_hat_suspect, f2_hat_suspect);

        // Optimize alternative hypothesis
        if (alt_hypothesis == HYPO_NOSUSPECT) {
                opt_context->hypo_type = HYPO_NOSUSPECT;
                loglik_alt = optimize_f2_only(f2_hat_alt);
                *f1_hat_alt = 0.0;  // No suspect means f1=0
        } else {
                opt_context->hypo_type = alt_hypothesis;
                loglik_alt = optimize_f1_f2(f1_hat_alt, f2_hat_alt);
        }

        // Return individual log-likelihoods if requested
        if (loglik_suspect_out) *loglik_suspect_out = loglik_suspect;
        if (loglik_alt_out) *loglik_alt_out = loglik_alt;

        return loglik_suspect - loglik_alt;
}

// ============================================================================
// COMMAND-LINE ARGUMENT PARSING
// ============================================================================

void print_usage(const char *progname)
{
        fprintf(stderr, "Usage: %s [options]\n\n", progname);
        fprintf(stderr, "Required options:\n");
        fprintf(stderr, "  -i, --infile <path>       Input JSON file\n");
        fprintf(stderr, "  -l, --lr <type>           Likelihood ratio type: L1, L2, L3, L4\n");
        fprintf(stderr, "                            L1: Suspect vs No-suspect\n");
        fprintf(stderr, "                            L2: Suspect vs Sibling\n");
        fprintf(stderr, "                            L3: Suspect vs Parent\n");
        fprintf(stderr, "                            L4: Suspect vs Cousin\n\n");
        fprintf(stderr, "Optional options:\n");
        fprintf(stderr, "  -o, --outfile <path>      Output file (default: stdout)\n");
        fprintf(stderr, "  -c, --contaminant <type>  Contaminant model: population, single_individual\n");
        fprintf(stderr, "                            (default: population)\n");
        fprintf(stderr, "  -f1 <value>               Initial f1 value (default: 0.2)\n");
        fprintf(stderr, "  -f2 <value>               Initial f2 value (default: 0.5)\n");
        fprintf(stderr, "  -k, --cousin_k <int>      Cousin degree for L4 (default: 1)\n");
        fprintf(stderr, "  -e, --error_adj <value>   Error adjustment parameter (default: 0.0)\n");
        fprintf(stderr, "  -h, --help                Show this help message\n");
}

Options parse_arguments(int argc, char *argv[])
{
        Options opts;
        // Set defaults
        opts.infile = NULL;
        opts.outfile = NULL;
        opts.lr_type = LR_L1;
        opts.contam_model = CONTAM_POPULATION;
        opts.f1_init = 0.2;
        opts.f2_init = 0.5;
        opts.cousin_k = 1;
        opts.error_adj = 0.0;

        int i;
        for (i = 1; i < argc; i++) {
                if (strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--infile") == 0) {
                        if (i + 1 < argc) {
                                opts.infile = argv[++i];
                        } else {
                                fprintf(stderr, "Error: -i/--infile requires an argument\n");
                                exit(1);
                        }
                } else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--outfile") == 0) {
                        if (i + 1 < argc) {
                                opts.outfile = argv[++i];
                        } else {
                                fprintf(stderr, "Error: -o/--outfile requires an argument\n");
                                exit(1);
                        }
                } else if (strcmp(argv[i], "-l") == 0 || strcmp(argv[i], "--lr") == 0) {
                        if (i + 1 < argc) {
                                i++;
                                if (strcmp(argv[i], "L1") == 0) {
                                        opts.lr_type = LR_L1;
                                } else if (strcmp(argv[i], "L2") == 0) {
                                        opts.lr_type = LR_L2;
                                } else if (strcmp(argv[i], "L3") == 0) {
                                        opts.lr_type = LR_L3;
                                } else if (strcmp(argv[i], "L4") == 0) {
                                        opts.lr_type = LR_L4;
                                } else {
                                        fprintf(stderr, "Error: Invalid LR type '%s'. Must be L1, L2, L3, or L4\n", argv[i]);
                                        exit(1);
                                }
                        } else {
                                fprintf(stderr, "Error: -l/--lr requires an argument\n");
                                exit(1);
                        }
                } else if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--contaminant") == 0) {
                        if (i + 1 < argc) {
                                i++;
                                if (strcmp(argv[i], "population") == 0) {
                                        opts.contam_model = CONTAM_POPULATION;
                                } else if (strcmp(argv[i], "single_individual") == 0) {
                                        opts.contam_model = CONTAM_SINGLE_INDIV;
                                } else {
                                        fprintf(stderr, "Error: Invalid contaminant model '%s'\n", argv[i]);
                                        exit(1);
                                }
                        } else {
                                fprintf(stderr, "Error: -c/--contaminant requires an argument\n");
                                exit(1);
                        }
                } else if (strcmp(argv[i], "-f1") == 0) {
                        if (i + 1 < argc) {
                                opts.f1_init = atof(argv[++i]);
                        } else {
                                fprintf(stderr, "Error: -f1 requires an argument\n");
                                exit(1);
                        }
                } else if (strcmp(argv[i], "-f2") == 0) {
                        if (i + 1 < argc) {
                                opts.f2_init = atof(argv[++i]);
                        } else {
                                fprintf(stderr, "Error: -f2 requires an argument\n");
                                exit(1);
                        }
                } else if (strcmp(argv[i], "-k") == 0 || strcmp(argv[i], "--cousin_k") == 0) {
                        if (i + 1 < argc) {
                                opts.cousin_k = atoi(argv[++i]);
                        } else {
                                fprintf(stderr, "Error: -k/--cousin_k requires an argument\n");
                                exit(1);
                        }
                } else if (strcmp(argv[i], "-e") == 0 || strcmp(argv[i], "--error_adj") == 0) {
                        if (i + 1 < argc) {
                                opts.error_adj = atof(argv[++i]);
                        } else {
                                fprintf(stderr, "Error: -e/--error_adj requires an argument\n");
                                exit(1);
                        }
                } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
                        print_usage(argv[0]);
                        exit(0);
                } else {
                        fprintf(stderr, "Error: Unknown option '%s'\n", argv[i]);
                        print_usage(argv[0]);
                        exit(1);
                }
        }

        // Validate required arguments
        if (opts.infile == NULL) {
                fprintf(stderr, "Error: Input file (-i/--infile) is required\n");
                print_usage(argv[0]);
                exit(1);
        }

        return opts;
}

// Helper function: convert nucleotide character to integer
// Skips leading/trailing whitespace
int nucleotide_to_int(char nuc)
{
        switch(nuc) {
                case 'A': case 'a': return NUC_A;
                case 'C': case 'c': return NUC_C;
                case 'G': case 'g': return NUC_G;
                case 'T': case 't': return NUC_T;
                case ' ': case '\t': case '\n': case '\r':
                        return -1;  // Whitespace - caller should skip
                default:
                        fprintf(stderr, "Error: Invalid nucleotide '%c' (ASCII %d)\n", nuc, (int)nuc);
                        exit(1);
        }
}

int calculate_likelihood_ratios(){
        fprintf(stderr, "\n=== Calculating Likelihood Ratios ===\n");

        FILE *outfp = stdout;
        if (global_opts.outfile) {
                outfp = fopen(global_opts.outfile, "w");
                if (!outfp) {
                        fprintf(stderr, "Error: Could not open output file %s\n", global_opts.outfile);
                        return 1;
                }
        }

        // Set up optimization context
        OptimizationContext context;
        context.data = global_snp_data;
        context.params.error_adj = global_opts.error_adj;
        context.params.cousin_k = global_opts.cousin_k;
        context.f1_init = global_opts.f1_init;
        context.f2_init = global_opts.f2_init;
        context.model = global_opts.contam_model;
        opt_context = &context;

        // Variables for results
        double log_lr;
        double f1_opt_suspect, f2_opt_suspect;
        double f1_opt_alt, f2_opt_alt;
        double loglik_suspect, loglik_alt;

        // Determine alternative hypothesis and labels
        HypothesisType alt_hypo;
        const char *lr_name, *alt_name;

        switch (global_opts.lr_type) {
                case LR_L1:
                        alt_hypo = HYPO_NOSUSPECT;
                        lr_name = "L1";
                        alt_name = "No-Suspect";
                        break;
                case LR_L2:
                        alt_hypo = HYPO_SIBLING;
                        lr_name = "L2";
                        alt_name = "Sibling";
                        break;
                case LR_L3:
                        alt_hypo = HYPO_PARENT;
                        lr_name = "L3";
                        alt_name = "Parent";
                        break;
                case LR_L4:
                        alt_hypo = HYPO_COUSIN;
                        lr_name = "L4";
                        alt_name = "Cousin";
                        break;
                default:
                        fprintf(stderr, "Error: Unknown likelihood ratio type\n");
                        if (outfp != stdout) fclose(outfp);
                        return 1;
        }

        // Compute likelihood ratio using unified function
        log_lr = compute_likelihood_ratio(alt_hypo, &f1_opt_suspect, &f2_opt_suspect, &f1_opt_alt, &f2_opt_alt,
                                          &loglik_suspect, &loglik_alt);

        // Print results
        const char *model_suffix = (global_opts.contam_model == CONTAM_SINGLE_INDIV) ? " (Single Individual Contaminant)" : "";

        fprintf(outfp, "\nLikelihood Ratio %s: Suspect vs. %s%s\n", lr_name, alt_name, model_suffix);
        if (global_opts.lr_type == LR_L4) {
                fprintf(outfp, "Cousin degree k=%d\n", global_opts.cousin_k);
        }
        fprintf(outfp, "Log likelihood ratio: %.6f\n\n", log_lr);

        fprintf(outfp, "Model 1 (Suspect): log_likelihood = %.6f, Parameter estimates: f1 = %.6f, f2 = %.6f\n",
                loglik_suspect, f1_opt_suspect, f2_opt_suspect);

        if (alt_hypo == HYPO_NOSUSPECT) {
                fprintf(outfp, "Model 2 (%s): log_likelihood = %.6f, Parameter estimates: f2 = %.6f\n",
                        alt_name, loglik_alt, f2_opt_alt);
        } else {
                fprintf(outfp, "Model 2 (%s): log_likelihood = %.6f, Parameter estimates: f1 = %.6f, f2 = %.6f\n",
                        alt_name, loglik_alt, f1_opt_alt, f2_opt_alt);
        }

        // Include population match test results
        fprintf(outfp, "\nPopulation match test Z-score: %.2f\n", global_mismatch_zscore);
        if (global_population_mismatch) {
                fprintf(outfp, "MISMATCH BETWEEN SUSPECT AND REFERENCE POPULATION. RESULTS MAY NOT BE RELIABLE.\n");
        }

        if (outfp != stdout) {
                fclose(outfp);
                fprintf(stderr, "=== Analysis complete ===\n");
                fprintf(stderr, "Results written to: %s\n", global_opts.outfile);
        } else {
                fprintf(stderr, "\n=== Analysis complete ===\n");
        }
        return 0;
}

int main(int argc, char *argv[])
{
        // Parse command-line arguments
        global_opts = parse_arguments(argc, argv);

        fprintf(stderr, "=== Forensic DNA Mixture Analysis ===\n");
        fprintf(stderr, "Input file: %s\n", global_opts.infile);
        fprintf(stderr, "Output file: %s\n", global_opts.outfile ? global_opts.outfile : "stdout");
        fprintf(stderr, "Likelihood ratio: L%d\n", global_opts.lr_type + 1);
        fprintf(stderr, "Contaminant model: %s\n",
                global_opts.contam_model == CONTAM_POPULATION ? "population" : "single_individual");
        fprintf(stderr, "Initial f1: %.3f, f2: %.3f\n", global_opts.f1_init, global_opts.f2_init);
        fprintf(stderr, "Error adjustment: %.3f\n", global_opts.error_adj);
        if (global_opts.lr_type == LR_L4) {
                fprintf(stderr, "Cousin degree k: %d\n", global_opts.cousin_k);
        }
        fprintf(stderr, "\n");

        // Parse the data file
        parse_data();

        // Test if suspect genotype matches reference population
        test_population_match();

        // Calculate the likelihood ratios
        calculate_likelihood_ratios();

        // Cleanup
        if (global_snp_data) {
                free_multisnp_data(global_snp_data);
        }

        return 0;
}        

