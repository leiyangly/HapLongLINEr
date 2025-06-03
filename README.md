# HapLongLINEr

**A modular pipeline for discovering and curating full-length young LINE-1 elements in haploid long-read human genome assemblies.**

**🚧 Early Alpha: This project is under active development. Features and interfaces may change without notice. Use with caution.**


## Overview

- HapLongLINEr discovers and curates full-length (≥5 kb) young LINE-1 elements (L1HS, L1PA2, and potentially intact L1PA3) in haploid long-read assemblies.  
- The pipeline supports both RepeatMasker-based and RepeatMasker-free approaches and provides a curated pangenome-level L1 sequence repository.
- Each L1 identified is marked with its intact ORF status.


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

Install dependencies (if not already installed):
```bash
conda install -c bioconda seqtk minimap2 emboss ncbi-blast+
```
Install HapLongLINEr with pip:
```bash
pip install -e .
```


## Usage

### Module 1: RepeatMasker-based

Input:
- Haploid assembly FASTA
- RepeatMasker BED or .out file (plain or gzipped)
- Reference genome FASTA (local or remote, e.g., hs1/hg38)

To run HapLongLINEr:
```bash
haplongliner rm --in your.genome.fa --mask repeatmasker.bed --reference hs1 --out output_dir
```
Or with a custom reference:
```bash
haplongliner rm --in your.genome.fa --mask repeatmasker.bed --custom custom_reference.fa.gz --out output_dir
```

Output:
- .HapLongLINEr.output.txt: Two-column, tab-delimited text files with L1 info from your assembly and corresponding hg38 coordinates
- .L1HSPA2PA3AllORF.intact.fa: FASTA file containing all intact L1HS, L1PA2, and L1PA3 sequences from the input assembly
- Additional curated FASTA files for the pangenome reference set
- Each L1 record includes its intact ORF status
- Each L1 record includes its lifted over coordinate on the used reference genome

### Module 2: RepeatMasker-free, SV-based

Documentation and example command coming soon.

### Module 3: Sequence Repository

See the repository folder or FASTA files with curated pangenome-level L1s.


## Authors

Lei Yang, Sara Nematbakhsh, Amanda Norseen, and Rick McLaughlin
