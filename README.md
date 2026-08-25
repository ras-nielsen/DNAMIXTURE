# Forensic DNA Mixture Analysis

A C-based tool for analyzing forensic DNA mixtures and calculating likelihood ratios to determine the probability that a suspect contributed DNA to a crime scene sample.

## QUICKSTART

### Install

Requires a C compiler and make:

```bash
make
```

This compiles all source files and creates the `DNAMIXTURE` executable.

### Run

```bash
./DNAMIXTURE -i <input_file> [options]
```

By default only the lambda evidence score is computed. Required arguments:

- `-i, --infile <file>`: Input JSON file containing DNA sequencing data

Optional arguments:

- `-l, --lr <type>`: Also compute a single likelihood ratio test: L1 (suspect vs. no-suspect), L2 (vs. sibling), L3 (vs. parent), or L4 (vs. cousin)
- `-o, --outfile <file>`: Output file path (default: stdout)
- `-c, --contaminant <model>`: Contaminant model: population (default) or single_individual
- `-f1 <value>`: Initial f1 value (suspect DNA proportion, default: 0.2)
- `-f2 <value>`: Initial f2 value (victim DNA proportion, default: 0.5)
- `-k, --cousin_k <degree>`: Cousin degree for L4 (default: 1 for first cousins)
- `-e, --error_adj <value>`: Error adjustment parameter (default: 0.0)
- `--no-lambda`: Skip the lambda evidence score and run only the requested single analysis (requires `-l`)
- `-C, --lambda-threshold <value>`: Gate threshold C for the lambda evidence score, on the likelihood-ratio scale (default: 10)
- `--no-victim`: No victim genome available; victim proportion fixed at 0 and only the suspect proportion is estimated

Examples:

```bash
./DNAMIXTURE -i data.json
./DNAMIXTURE -i data.json -l L1
./DNAMIXTURE -i data.json -l L4 -k 2
```

### Preparing input from BAM/VCF files

A converter that builds DNAMIXTURE's JSON input from standard formats — a BAM file for the crime-stain reads, VCFs for the suspect and (optionally) the victim, and an allele-frequency panel — is maintained in the companion pipeline repository:

```
git clone https://github.com/GeoGenetics/forensic-mixture-simulations
pip install pysam
python forensic-mixture-simulations/tools/dnamixture_prep.py \
    --bam stain.bam --suspect-vcf suspect.vcf.gz --victim-vcf victim.vcf.gz \
    --panel forensic-mixture-simulations/panels/1000g.phase3.maf01.sites500k.vcf.gz \
    --af-field EUR_AF -o case.json
DNAMIXTURE -i case.json
```

Omit `--victim-vcf` when no victim genome is available and add `--no-victim` to the DNAMIXTURE call. See the pipeline repository's README and `panels/README.md` for details. THE BUNDLED REFERENCE PANELS ARE SUPPLIED FOR TESTING PURPOSES ONLY. FOR FORENSICS APPLICATIONS PLEASE USE A PANEL TAILORED TO YOUR USE CASE.

## DETAILED INFORMATION

### Overview

This program analyzes DNA sequencing data from crime scene samples that contain mixtures of DNA from multiple individuals (suspect, victim, and possible contaminants). It calculates likelihood ratios to evaluate different hypotheses about the source of the DNA.

### Features

- **Multiple Likelihood Ratio Tests**:
  - **L1**: Suspect vs. No-Suspect
  - **L2**: Suspect vs. Sibling
  - **L3**: Suspect vs. Parent
  - **L4**: Suspect vs. Cousin (configurable degree)

- **Contaminant Models**:
  - Population-based contaminant (general population)
  - Single individual contaminant

- **Quality Score Handling**:
  - Incorporates numeric Phred quality scores (0-60) into likelihood calculations
  - Optional error adjustment parameter (Equation 14 of the paper)

- **Large-Scale Data Processing**:
  - Handles large JSON input files (tested with 186MB files)
  - Processes thousands of genomic positions
  - Memory-efficient parsing

- **Population Match Test**:
  - Automatically tests if suspect genotype is consistent with reference population
  - Compares observed statistic to simulated distribution under Hardy-Weinberg Equilibrium
  - Warns if suspect appears mismatched with reference population (results may be unreliable)

### Build details

**Other Make targets:**

```bash
make clean      # Remove compiled objects and executable
make debug      # Build with debug symbols
make test       # Build and run test with a.json
make install    # Install to /usr/local/bin (requires sudo)
make help       # Show all available targets
```

**Manual compilation** without Make:

```bash
gcc -c dnamixture.c -o DNAMIXTURE.o -Wall
gcc -c json_parser.c -o json_parser.o -Wall
gcc -c neldermead.c -o neldermead.o -Wall
gcc DNAMIXTURE.o json_parser.o neldermead.o -o DNAMIXTURE -lm
```

### Input File Format

The program expects a JSON file with the following structure:

```json
{
  "position_reads": {
    "chr:position": {
      "reads": [["A", 37], ["G", 37], ...],
      "suspect_gt": ["A", "G"],
      "victim_gt": ["G", "G"]
    },
    ...
  },
  "population_freqs": {
    "chr:position": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
    ...
  },
  "random_positions": {
    "chr:position": {
      "A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25,
      "suspect_gt": ["A", "G"]
    },
    ...
  }
}
```

Fields:

- **position_reads**: Dictionary of genomic positions
  - Key: "chromosome:position" (e.g., "10:1978952")
  - **reads**: Array of [nucleotide, quality_score] pairs
    - nucleotide: A, C, G, or T
    - quality_score: PHRED quality score (0-60)
  - **suspect_gt**: Suspect genotype [allele1, allele2]
  - **victim_gt**: Victim genotype [allele1, allele2]

- **population_freqs**: Allele frequencies for each position
  - Key: "chromosome:position"
  - Value: Dictionary of nucleotide frequencies (must sum to 1.0)

- **random_positions**: Independent SNPs used exclusively by the population match test
  - Key: "chromosome:position"
  - **A**, **C**, **G**, and **T**: Reference-population allele frequencies
  - **suspect_gt**: Suspect genotype [allele1, allele2]
  - These entries do not need corresponding reads or victim genotypes in **position_reads**

### Output Format

Example output (L1 with the default lambda evidence score):

```
Likelihood Ratio L1: Suspect vs. No-Suspect
Log likelihood ratio: 295.779987

Model 1 (Suspect): log_likelihood = -533.172644, Parameter estimates: f1 = 0.842945, f2 = 0.090513
Model 2 (No-Suspect): log_likelihood = -828.952632, Parameter estimates: f2 = 0.110461

Lambda evidence score (threshold C = 10):
  log L1 (population contaminant):        295.779987
  log L1 (single individual contaminant): 250.032791
  log L2 (population contaminant):        74.831178
  log L2 (single individual contaminant): 74.591649
  Gate passed: max(log L2) = 74.831178 > log(C) = 2.302585
  log Lambda = min(log L1) = 250.032791

Population match test Z-score: 0.80
```

Output fields:

- **Log likelihood ratio**: Natural logarithm of the likelihood ratio for the requested test
  - Positive values favor the suspect hypothesis
  - Negative values favor the alternative hypothesis
  - Magnitude indicates strength of evidence
- **Model 1 / Model 2**: Maximized log-likelihood and estimated mixture proportions (f1 = suspect or relative, f2 = victim) under each hypothesis
- **Lambda evidence score block**: The four component log likelihood ratios, the gate decision, and log Lambda (see below)
- **Population match test Z-score**: See below; reported as "not performed" if no background SNPs are supplied

### Lambda Evidence Score

The lambda evidence score is the program's default output: with no `-l`
option only lambda is computed, and with `-l` it is reported alongside the
requested likelihood ratio (disable with `--no-lambda`). Lambda combines the
suspect-vs-no-suspect and suspect-vs-sibling tests across both contaminant
models:

- Let L1 and L2 be the suspect-vs-no-suspect and suspect-vs-sibling log
  likelihood ratios computed under the population-contaminant model, and let
  L1~ and L2~ be the same quantities computed under the single-individual
  contaminant model.
- If max(L2, L2~) > log(C), report log Lambda = min(L1, L1~).
- Otherwise report log Lambda = 0 (no evidence).

The threshold C is on the likelihood-ratio scale and defaults to 10 (i.e., the
evidence must favor the suspect over a hypothesized sibling by at least a
factor of 10 under at least one contaminant model). It can be changed with
`-C/--lambda-threshold`. Gating on the maximum of the two L2 values makes the
score robust to misspecification of the contaminant model, and reporting the
minimum of the two L1 values is the conservative bound on the evidence. All
logarithms are natural logs, as elsewhere in the output.

### Interpreting Results

Positive log likelihood ratios favor the suspect hypothesis, negative values
the alternative hypothesis. To convert to a likelihood ratio: LR = exp(log_LR).

DNA mixture proportions:

- **f1**: Proportion of suspect (or relative) DNA in mixture
- **f2**: Proportion of victim DNA in mixture
- **1 - f1 - f2**: Proportion of contaminant DNA

Values are optimized to maximize likelihood across all genomic positions.

### Population Match Test

Before calculating likelihood ratios, the program tests whether the suspect's genotype is consistent with the reference population allele frequencies. This helps detect cases where the suspect may be from a different population than the reference panel, which could lead to unreliable results.

The test reads its suspect genotypes and allele frequencies from the top-level **random_positions** object. It does not use suspect genotypes from **position_reads** or frequencies from **population_freqs**. Those two sections continue to supply the DNA-mixture likelihood calculations.

**Test statistic:**

```
t = sum_{i=1}^{S} (f_i - g_i)^2
```

Where:

- S = number of entries loaded from **random_positions**
- f_i = frequency of the major allele in the reference population at SNP i
- g_i = suspect genotype coded as 0, 0.5, or 1 (homozygous minor, heterozygous, or homozygous major)

**Procedure:**
1. Calculate observed t for the suspect
2. Simulate 1000 random genotypes under Hardy-Weinberg Equilibrium using the reference allele frequencies
3. Calculate t for each simulated genotype
4. Compute Z-score: (observed t - mean of simulated t) / SD of simulated t

**Interpretation:**
- If |Z-score| > 2, a warning is printed: `MISMATCH BETWEEN SUSPECT AND REFERENCE POPULATION. RESULTS MAY NOT BE RELIABLE.`
- This warning appears both on stderr during execution and in the output file
- If no **random_positions** entries are supplied, the test is not performed and the output states this explicitly

### Likelihood Ratio Types

**L1: Suspect vs. No-Suspect.** Compares:
- H1: Mixture contains suspect + victim + contaminant
- H2: Mixture contains victim + contaminant (no suspect)

**L2: Suspect vs. Sibling.** Compares:
- H1: Mixture contains suspect + victim + contaminant
- H2: Mixture contains suspect's sibling + victim + contaminant

Uses IBD (Identity By Descent) probabilities:
- P(IBD=0) = 0.25
- P(IBD=1) = 0.50
- P(IBD=2) = 0.25

**L3: Suspect vs. Parent.** Compares:
- H1: Mixture contains suspect + victim + contaminant
- H2: Mixture contains suspect's parent + victim + contaminant

Parent and child always share exactly one allele (IBD=1).

**L4: Suspect vs. Cousin.** Compares:
- H1: Mixture contains suspect + victim + contaminant
- H2: Mixture contains suspect's k-th degree cousin + victim + contaminant

Uses cousin IBD probabilities:
- P(IBD=1) = (1/2)^(2k)
- P(IBD=0) = 1 - (1/2)^(2k)

Where k is the cousin degree:
- k=1: First cousins (share 1/8 of DNA)
- k=2: Second cousins (share 1/32 of DNA)
- k=3: Third cousins (share 1/128 of DNA)

### Algorithm Details

The values of f1 and f2 that maximize the likelihood across all SNP positions are found with the Nelder-Mead simplex algorithm in double precision. A softmax transformation maps the constrained parameters (f1, f2 ≥ 0, f1+f2 ≤ 1) to an unconstrained space.

If the maximum lies on a boundary of the parameter space (e.g., f1 = 1), the search stops after a fixed number of function evaluations and returns the best point found, printing a warning to stderr. The reported likelihood is still the maximum; the warning only signals a boundary optimum.

For each sequencing read, the likelihood incorporates:
1. Base call quality scores (PHRED-scaled error probabilities)
2. True nucleotide at the position
3. Observed nucleotide in the read
4. Optional error adjustment parameter (Equation 14 of the paper)

Error probability: P_error = 10^(-Q/10)

With error adjustment (if e > 0):
```
p'(x|nt) = (p(x|nt) + e) / (1 + 4e)
```

For each possible contributor genotype combination, the program calculates:
1. Probability of mixture given genotypes and proportions (f1, f2)
2. Probability of each read given the mixture
3. Product across all reads at the position
4. Sum across all SNP positions (in log space)

### Error Handling

The program performs validation on:
- **Quality scores**: Must be in range 0-60 (exits with error if outside range)
- **Input file**: Must exist and be valid JSON
- **Nucleotides**: Must be A, C, G, or T
- **Command-line arguments**: Required arguments must be provided

Error messages are written to stderr.

### Files

- **dnamixture.c**: Main program with likelihood calculations and optimization
- **json_parser.c**: JSON parsing functions and hash table implementation
- **json_parser.h**: Header file with type definitions and function declarations
- **neldermead.c**: Nelder-Mead optimization implementation
- **neldermead.h**: Header for Nelder-Mead functions
- **Makefile**: Build automation with compilation targets
- **README.md**: This file

## Citation

If you use this software in your research, please cite the associated paper:

[Citation information to be added]

## References

- Nelder, J. A. & Mead, R. (1965). A simplex method for function minimization. The Computer Journal 7:308-313
- PHRED quality scores: Ewing & Green (1998)
- Identity by descent calculations for relatives
- Forensic DNA mixture interpretation methods

## Contact

Rasmus Nielsen <rasmus_nielsen@berkeley.edu>

## License

Copyright (c) 2026 Rasmus Nielsen, Abigail Daisy Ramsø and Thorfinn Korneliussen.
Licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0
International License (CC BY-NC-SA 4.0). See the LICENSE file for details.

## Disclaimer

This program is provided "as is", without warranty of any kind, express or
implied, including any warranty of accuracy or fitness for a particular
purpose. The authors and copyright holders accept no liability for any claim
or damages arising from the program or its use. This is research software; it
is not certified for forensic casework, and its output should not be used as
the sole basis for any legal or investigative decision. Users are solely
responsible for validating results and for any decisions based on them.
