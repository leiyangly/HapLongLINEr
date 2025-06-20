# HapLongLINEr

**A modular pipeline for discovering and curating full-length young LINE-1 elements in haploid long-read human genome assemblies.**


## Overview

- HapLongLINEr discovers and curates full-length (≥5 kb by default) young LINE-1 elements (L1HS, L1PA2, and potentially intact L1PA3) in haploid long-read assemblies.  
- The pipeline supports both RepeatMasker-based and RepeatMasker-free approaches and provides a curated pangenome-level L1 sequence repository.
- Each L1 identified is marked with its intact ORF status.
- 🚧 Early Alpha: This project is under active development. Features and interfaces may change without notice. Use with caution.


## Installation

### Using Conda

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

### Using Git

Install dependencies if necessary:
```bash
conda install -c bioconda seqtk minimap2 bedtools emboss blast
```

Clone the repository:
```bash
git clone https://github.com/leiyangly/HapLongLINEr.git
```

Install HapLongLINEr with pip:
```bash
cd HapLongLINEr

pip install -e .
```


## Usage

### Module 1: RepeatMasker-based

Input:
- Haploid assembly FASTA
- RepeatMasker BED or .out file (plain or gzipped)
- Reference genome FASTA (local or remote, e.g., hs1/hg38)

Command:
```bash
haplongliner rm --in your.genome.fa --mask repeatmasker.bed --reference hs1 --out output_dir
```
Or:
```bash
haplongliner rm --in your.genome.fa --mask repeatmasker.bed --custom custom_reference.fa.gz --out output_dir
```
To troubleshoot malformed RepeatMasker entries, use the optional
`--log-skipped` parameter to record skipped lines:
```bash
haplongliner rm --in your.genome.fa --mask repeatmasker.bed \
  --reference hs1 --out output_dir --log-skipped skipped.log
```

Output:
- OUT.TXT file with L1 info from your assembly and corresponding refence genome (hs1/hg38) coordinates and ORF status
- LOG.TXT file that summarizes results of each step of the pipeline module
- FASTA file containing all full length (>=5kb by default) L1HS, L1PA2, and intact L1PA3 sequences from the input assembly


### Module 2: RepeatMasker-free, SV-based

Input:
- Haploid assembly FASTA
- Structural variant (SV) callset (e.g., VCF or BED)
- Pangenome-level L1 reference FASTA

Command:
```bash
haplongliner sv --in your.genome.fa --sv your.sv.vcf --l1ref pangenome_L1_reference.fa --out output_file.bed
```
Output:
- OUT.TXT file with L1 info from your assembly and corresponding refence genome (hs1/hg38) coordinates and ORF status
- LOG.TXT file that summarizes results of each step of the pipeline module
- FASTA file containing all full length (>=5kb by default) L1HS, L1PA2, and intact L1PA3 sequences from the input assembly

### Module 3: Sequence Repository

Module 3 builds a sequence repository from previously identified insertions.
Currently the only required argument is the output directory.

Command:
```bash
haplongliner db --out output_folder
```

Output:
- Annotated BED file that adds the following columns to each L1 record: Frequency of presence in HPRC haploids; Intactness status in HPRC; Liftover coordinate in the chosen reference (hs1 or hg38)
- FASTA file that Contains all L1 sequences at the insertion site from all HPRC haploids that carry that L1


## Authors

Lei Yang, Sara Nematbakhsh, Amanda Norseen, and Rick McLaughlin
