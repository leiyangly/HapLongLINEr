# HapLongLINEr

**A modular pipeline for discovering and curating full-length young LINE-1 elements in haploid long-read human genome assemblies.**

## Overview

HapLongLINEr discovers and curates full-length (≥5 kb) young LINE-1 elements (L1HS, L1PA2, and potentially intact L1PA3) in haploid long-read assemblies.  
**Each L1 identified is marked with its intact ORF status.**  
The pipeline supports both RepeatMasker-based and RepeatMasker-free approaches and provides a curated pangenome-level L1 sequence repository.

## Features

- Full-length L1 discovery from haploid long-read assemblies
- Intact ORF detection for each L1
- RepeatMasker-based and SV-based workflows
- Curated pangenome-level L1 sequence repository
- Flexible input: BED, .out, gzipped formats
- Modern, modular CLI

## Installation

### Using Conda (recommended)

Enable the required channels:
```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Install HapLongLINEr:
```bash
conda install haplongliner
```

### Manual Installation

Clone the repository:
```bash
git clone https://github.com/yourusername/HapLongLINEr.gitcd HapLongLINEr
```

Install dependencies (if not using conda):
```bash
conda install -c bioconda seqtk minimap2 emboss ncbi-blast+
```

## System Requirements

seqtk
minimap2 (2.18+ if remote reference genome is needed)
EMBOSS (for getorf)
ncbi-blast+
Usage

### Module 1: RepeatMasker-based

haplongliner rm --in your.genome.fa --mask repeatmasker.bed --reference hs1 --out output_dir
Or with a custom reference:

haplongliner rm --in your.genome.fa --mask repeatmasker.bed --custom custom_reference.fa.gz --out output_dir

### Module 2: RepeatMasker-free, SV-based

Documentation and example command coming soon.

### Module 3: Sequence Repository

See the repository folder or FASTA files with curated pangenome-level L1s.

## Input & Output

### Module 1

Input:

Haploid assembly FASTA
RepeatMasker BED or .out file (plain or gzipped)
Reference genome FASTA (local or remote, e.g., hs1/hg38)

Output:

.HapLongLINEr.output.txt: Two-column, tab-delimited text files with L1 info from your assembly and corresponding hg38 coordinates
.L1HSPA2PA3AllORF.intact.fa: FASTA file containing all intact L1HS, L1PA2, and L1PA3 sequences from the input assembly
Additional curated FASTA files for the pangenome reference set
Each L1 record includes its intact ORF status

### Module 2

Documentation and example command coming soon.

### Module 3

Curated FASTA files for the pangenome reference set

## Authors

Lei Yang, Sara Nematbakhsh, Amanda Norseen, and Rick McLaughlin
