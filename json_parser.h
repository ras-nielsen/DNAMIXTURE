#ifndef JSON_PARSER_H
#define JSON_PARSER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Nucleotide encoding
#define NUC_A 0
#define NUC_C 1
#define NUC_G 2
#define NUC_T 3

// Forward declarations of types defined in testfunc.c
typedef struct {
    int nucleotide;      // observed nucleotide: 0=A, 1=C, 2=G, 3=T
    double quality;      // phred quality score (0-60 typical range)
} Read;

typedef struct {
    int victim_genotype[2];    // victim's diploid genotype at this SNP
    int suspect_genotype[2];   // suspect's diploid genotype at this SNP
    Read *reads;               // array of reads at this SNP
    int numreads;              // number of reads at this SNP
    double AF[4];              // allele frequencies [A, C, G, T]
    double *cache;             // precomputed f-independent log-likelihood terms
} SNPData;

typedef struct {
    SNPData *snps;       // array of SNP data
    int numsnps;         // total number of SNPs
} MultiSNPData;

typedef struct {
    char *infile;                     // Input JSON file path
    char *outfile;                    // Output file path (NULL = stdout)
    int lr_type;                      // Which likelihood ratio to calculate
    int lr_given;                     // 1 = -l/--lr was supplied (0 = lambda-only run)
    int contam_model;                 // Contaminant model
    double f1_init;                   // Initial f1 value
    double f2_init;                   // Initial f2 value
    int cousin_k;                     // Cousin degree (for L4)
    double error_adj;                 // Error adjustment parameter
    int compute_lambda;               // 1 = also compute the lambda evidence score (default)
    double lambda_threshold;          // gate threshold C on the LR scale (default 10.0)
    int no_victim;                    // 1 = no victim genome: f2 fixed at 0, 1D optimization of f1
} Options;

// Hash table for fast allele frequency lookup
typedef struct PopFreqNode {
        char *position_key;           // e.g., "10:1978952"
        double AF[4];                 // Allele frequencies [A, C, G, T]
        struct PopFreqNode *next;     // For chaining
} PopFreqNode;

typedef struct {
        PopFreqNode **buckets;
        int size;                     // Number of buckets
        int count;                    // Number of entries
} PopFreqHashTable;

typedef struct {
    double AF[4];
    int suspect_gt[2];
} RandomSNP;

typedef struct {
    int count;
    RandomSNP *snps;
} RandomSNPSet;

extern RandomSNPSet *global_random_snps;

// External global variables
extern PopFreqHashTable *global_popfreq_table;
extern MultiSNPData *global_snp_data;
extern Options global_opts;

// External functions from testfunc.c needed by parser
extern MultiSNPData* allocate_multisnp_data(int numsnps);
extern void allocate_snp_reads(SNPData *snp, int numreads);
extern void set_snp_read(SNPData *snp, int read_idx, int nucleotide, double quality);
extern void set_snp_af(SNPData *snp, double af[4]);
extern void set_snp_genotypes(SNPData *snp, int victim_g1, int victim_g2, int suspect_g1, int suspect_g2);
extern int nucleotide_to_int(char nuc);

// File loading
char* load_file(const char *filename, size_t *file_size);

// JSON parsing helper functions
char* find_string(char *haystack, const char *needle);
char* skip_whitespace(char *p);
char* parse_quoted_string(char **p);
double parse_number(char **p);
char* skip_json_value(char *p);
int string_to_nucleotide(const char *str);

// Hash table functions
unsigned int hash_position_key(const char *key, int table_size);
PopFreqHashTable* create_popfreq_table(int size);
void insert_popfreq(PopFreqHashTable *table, const char *position_key, double AF[4]);
int lookup_popfreq(PopFreqHashTable *table, const char *position_key, double AF[4]);
void free_popfreq_table(PopFreqHashTable *table);

// JSON parsing functions
void parse_population_freqs(char *json_start);
int lookup_population_freqs(char *json_start, const char *position_key, double AF[4]);
int count_positions(char *json_start);
char* parse_position_entry(char *p, int snp_idx, char *json_start);
int parse_data(void);

#endif // JSON_PARSER_H
