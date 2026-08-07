# Forensic DNA Mixture Analysis

A C-based tool for analyzing forensic DNA mixtures and calculating likelihood ratios to determine the probability that a suspect contributed DNA to a crime scene sample.

## Overview

This program analyzes DNA sequencing data from crime scene samples that contain mixtures of DNA from multiple individuals (suspect, victim, and possible contaminants). It calculates likelihood ratios to evaluate different hypotheses about the source of the DNA.

## Features

- **Multiple Likelihood Ratio Tests**:
  - **L1**: Suspect vs. No-Suspect
  - **L2**: Suspect vs. Sibling
  - **L3**: Suspect vs. Parent
  - **L4**: Suspect vs. Cousin (configurable degree)

- **Contaminant Models**:
  - Population-based contaminant (general population)
  - Single individual contaminant

- **Quality Score Handling**:
  - Supports PHRED33 and PHRED64 quality score encoding
  - Incorporates sequencing error rates into likelihood calculations
  - Optional error adjustment parameter (Equation 11)

- **Large-Scale Data Processing**:
  - Handles large JSON input files (tested with 186MB files)
  - Processes thousands of genomic positions
  - Memory-efficient parsing

- **Population Match Test**:
  - Automatically tests if suspect genotype is consistent with reference population
  - Compares observed statistic to simulated distribution under Hardy-Weinberg Equilibrium
  - Warns if suspect appears mismatched with reference population (results may be unreliable)

## Installation

### Prerequisites

- GCC compiler
- Standard C libraries (stdio, stdlib, string, math)
- Make utility (usually pre-installed on Unix/Linux/macOS)

### Compilation

#### Using Make (Recommended)

The project includes a Makefile for easy compilation:

```bash
make
```

This will compile all source files and create the `testfunc` executable.

**Other Make targets:**

```bash
make clean      # Remove compiled objects and executable
make debug      # Build with debug symbols
make test       # Build and run test with a.json
make install    # Install to /usr/local/bin (requires sudo)
make help       # Show all available targets
```

#### Manual Compilation

If you prefer to compile manually without Make:

```bash
gcc -c testfunc.c -o testfunc.o -Wall
gcc -c json_parser.c -o json_parser.o -Wall
gcc -c neldermead.c -o neldermead.o -Wall
gcc testfunc.o json_parser.o neldermead.o -o testfunc -lm
```

### Project Structure

The codebase is organized into the following files:

- **testfunc.c**: Main program with likelihood calculations and optimization
- **json_parser.c**: JSON parsing functions and hash table implementation
- **json_parser.h**: Header file with type definitions and function declarations
- **neldermead.c**: Nelder-Mead optimization implementation
- **neldermead.h**: Header for Nelder-Mead functions
- **Makefile**: Build automation

## Usage

### Basic Command Line

```bash
./testfunc -i <input_file> -l <likelihood_ratio> [options]
```

### Required Arguments

- `-i, --infile <file>`: Input JSON file containing DNA sequencing data
- `-l, --lr <type>`: Likelihood ratio to calculate: L1, L2, L3, or L4

### Optional Arguments

- `-o, --outfile <file>`: Output file path (default: stdout)
- `-c, --contaminant <model>`: Contaminant model: population (default) or single_individual
- `-f1 <value>`: Initial f1 value (suspect DNA proportion, default: 0.2)
- `-f2 <value>`: Initial f2 value (victim DNA proportion, default: 0.5)
- `-k, --cousin_k <degree>`: Cousin degree for L4 (default: 1 for first cousins)
- `-e, --error_adj <value>`: Error adjustment parameter (default: 0.0)

### Examples

#### Basic L1 Analysis (Suspect vs. No-Suspect)
```bash
./testfunc -i data.json -l L1
```

#### L2 Analysis with Output File
```bash
./testfunc -i data.json -l L2 -o results.txt
```

#### L4 Analysis for Second Cousins
```bash
./testfunc -i data.json -l L4 -k 2
```

#### Analysis with Single Individual Contaminant
```bash
./testfunc -i data.json -l L1 -c single
```

#### Analysis with Custom Initial Parameters
```bash
./testfunc -i data.json -l L1 -f1 0.3 -f2 0.6
```

#### Analysis with Error Adjustment
```bash
./testfunc -i data.json -l L1 -e 0.01
```

## Input File Format

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
  }
}
```

### Fields:

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

## Output Format

The program outputs tab-separated values with the following columns:

### L1 Output (Suspect vs. No-Suspect)
```
log_LR	f1_suspect	f2_suspect	f2_nosuspect
12884.852539	0.497724	0.502274	0.477554
```

### L2, L3, L4 Output (Suspect vs. Relative)
```
log_LR	f1_suspect	f2_suspect	f1_relative	f2_relative
1733.559570	0.497724	0.502274	0.539563	0.460437
```

### Output Fields:

- **log_LR**: Natural logarithm of the likelihood ratio
  - Positive values favor the suspect hypothesis
  - Negative values favor the alternative hypothesis
  - Magnitude indicates strength of evidence
- **f1_suspect/f1_relative**: Optimized proportion of primary contributor DNA
- **f2_suspect/f2_relative**: Optimized proportion of secondary contributor DNA
- **f2_nosuspect**: Optimized victim proportion under no-suspect hypothesis (L1 only)

## Interpreting Results

### Log Likelihood Ratio (log_LR)

The log likelihood ratio quantifies the strength of evidence:

- **log_LR > 0**: Evidence supports suspect hypothesis
- **log_LR < 0**: Evidence supports alternative hypothesis
- **|log_LR| > 10**: Very strong evidence (LR > 20,000)
- **|log_LR| > 5**: Strong evidence (LR > 148)
- **|log_LR| > 2**: Moderate evidence (LR > 7)
- **|log_LR| < 2**: Weak evidence

To convert to likelihood ratio: LR = exp(log_LR)

### DNA Mixture Proportions

- **f1**: Proportion of suspect (or relative) DNA in mixture
- **f2**: Proportion of victim DNA in mixture
- **1 - f1 - f2**: Proportion of contaminant DNA

Values are optimized to maximize likelihood across all genomic positions.

### Population Match Test

Before calculating likelihood ratios, the program tests whether the suspect's genotype is consistent with the reference population allele frequencies. This helps detect cases where the suspect may be from a different population than the reference panel, which could lead to unreliable results.

**Test statistic:**
```
t = sum_{i=1}^{S} (f_i - g_i)^2
```

Where:
- S = number of SNPs with non-missing suspect genotype
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

## Likelihood Ratio Types

### L1: Suspect vs. No-Suspect
Compares:
- H1: Mixture contains suspect + victim + contaminant
- H2: Mixture contains victim + contaminant (no suspect)

### L2: Suspect vs. Sibling
Compares:
- H1: Mixture contains suspect + victim + contaminant
- H2: Mixture contains suspect's sibling + victim + contaminant

Uses IBD (Identity By Descent) probabilities:
- P(IBD=0) = 0.25
- P(IBD=1) = 0.50
- P(IBD=2) = 0.25

### L3: Suspect vs. Parent
Compares:
- H1: Mixture contains suspect + victim + contaminant
- H2: Mixture contains suspect's parent + victim + contaminant

Parent and child always share exactly one allele (IBD=1).

### L4: Suspect vs. Cousin
Compares:
- H1: Mixture contains suspect + victim + contaminant
- H2: Mixture contains suspect's k-th degree cousin + victim + contaminant

Uses cousin IBD probabilities:
- P(IBD=1) = (1/2)^(2k+1)
- P(IBD=0) = 1 - (1/2)^(2k+1)

Where k is the cousin degree:
- k=1: First cousins (share 1/8 of DNA)
- k=2: Second cousins (share 1/32 of DNA)
- k=3: Third cousins (share 1/128 of DNA)

## Algorithm Details

### Optimization Method

The values of f1 and f2 that maximize the likelihood across all SNP positions are found with the Nelder-Mead simplex algorithm in double precision. A softmax transformation maps the constrained parameters (f1, f2 ≥ 0, f1+f2 ≤ 1) to an unconstrained space.

If the maximum lies on a boundary of the parameter space (e.g., f1 = 1), the search stops after a fixed number of function evaluations and returns the best point found, printing a warning to stderr. The reported likelihood is still the maximum; the warning only signals a boundary optimum.

### Read Likelihood Calculation

For each sequencing read, the likelihood incorporates:
1. Base call quality scores (PHRED-scaled error probabilities)
2. True nucleotide at the position
3. Observed nucleotide in the read
4. Optional error adjustment parameter (Equation 11)

Error probability: P_error = 10^(-Q/10)

With error adjustment (if e > 0):
```
p'(x|nt) = (p(x|nt) + e) / (1 + 4e)
```

### Genotype Likelihood

For each possible contributor genotype combination, the program calculates:
1. Probability of mixture given genotypes and proportions (f1, f2)
2. Probability of each read given the mixture
3. Product across all reads at the position
4. Sum across all SNP positions (in log space)

## Error Handling

The program performs validation on:
- **Quality scores**: Must be in range 0-60 (exits with error if outside range)
- **Input file**: Must exist and be valid JSON
- **Nucleotides**: Must be A, C, G, or T
- **Command-line arguments**: Required arguments must be provided

Error messages are written to stderr.

## Files

### Source Code
- **testfunc.c**: Main program with likelihood calculations and optimization
- **json_parser.c**: JSON parsing functions and hash table implementation
- **json_parser.h**: Header file with type definitions and function declarations
- **neldermead.c**: Nelder-Mead optimization implementation
- **neldermead.h**: Header for Nelder-Mead functions
- **Makefile**: Build automation with compilation targets

### Documentation
- **README.md**: This file - comprehensive documentation
- **paper_draft.tex**: LaTeX document with mathematical formulas (reference)

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
