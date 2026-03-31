# HapLongLINEr

**A modular pipeline for discovering and curating full-length young LINE-1 elements in haploid long-read human genome assemblies.**

## Overview

- HapLongLINEr primarily discovers and curates full-length (≥5 kb) young LINE-1 elements (L1HS and L1PA2) in haploid long-read assemblies.
- The pipeline focuses on full-length L1s, but can be expanded to support all transposable elements at any length.
- HapLoneLINEr can detect intact ORFs when ran on L1s with the length cutoff of >=5kb.
- When only L1 families are selected, the minimum length defaults to 5 kb; otherwise it defaults to 0 bp.

## Installation

### Install using git and pip

Create the recommended conda environment:
```bash
conda env create -f environment.yml
conda activate haplongliner
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

This editable install provides the Python package plus the bundled small data
assets used by HapLongLINEr, including `L1rp.fa`, the ORF BLAST database, and
the HPRC sequence lookup archives. Large reference downloads such as `hs1.fa.gz`
and UCSC RepeatMasker annotations are cached under `~/.cache/haplongliner` by
default. Override these locations with `HAPLONGLINER_DATA_DIR` or
`HAPLONGLINER_CACHE_DIR` if needed.

### Install as a conda package

Enable the required channels:
```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Install HapLongLINEr once the conda package is published:
```bash
conda install haplongliner
```

## Usage

Run `haplongliner`, `haplongliner rm` (RM mode), `haplongliner mm` (MM mode),
`haplongliner sv` (SV mode) or `haplongliner seq` with no arguments to see the
available options for each command.

### RM mode: Based on RepeatMasker annotations, with optional auto-masking

Input:
- Haploid assembly FASTA
- RepeatMasker BED or .out file (plain or gzipped), or omit `--mask` to run RepeatMasker on the assembly first
- Reference genome FASTA (hs1, hg38, or custom; default hs1)
- Transposable element families to analyse (via `--te`, default `L1,L1PA3`)

Command with test genome:
```bash
haplongliner rm \
  --in tests/test.genome.HG00410.1.fa \
  --mask tests/test.genome.HG00410.1.out \
  --ref hs1 \
  --out test_output_rm
```

Or with user provided genome:
```bash
haplongliner rm \
  --in your.genome.fa \
  --mask your.repeatmasker.bed \
  --ref hs1 \
  --out your_output_dir
```

Or let HapLongLINEr mask the assembly for you first:
```bash
haplongliner rm \
  --in your.genome.fa \
  --ref hs1 \
  --out your_output_dir
```

To search additional TE families, include the `--te` option. For example,
`--te L1,SVA` searches for both L1 and SVA elements.

Intact ORF detection is performed only when `--length` is at least 5000 and the `--te` list includes L1HS, L1PA2, or L1PA3; otherwise these ORF-related BLASTP and sequence extraction steps are skipped.

Output:
- haplongliner_rm.bed file with TE info from your assembly and corresponding reference genome (hs1/hg38) coordinates and ORF status
- haplongliner_rm.fa file containing all full length (>=5kb by default) sequences for the selected TE families
- Log file (when using `-g/--log`) that summarizes results of each step of the pipeline mode

### MM mode: Based on minimap2 seeding to L1rp followed by candidate-only RepeatMasker

Input:
- Haploid assembly FASTA
- Reference genome FASTA (hs1, hg38, or custom; default hs1)
- L1 families to analyse (via `--te`, default `L1,L1PA3`; MM mode is L1-only)

Command with test genome:
```bash
haplongliner mm \
  --in tests/test.genome.HG00410.1.fa \
  --ref hs1 \
  --out test_output_mm
```

Command with user provided genome:
```bash
haplongliner mm \
  --in your.genome.fa \
  --te L1,L1PA3 \
  --ref hs1 \
  --out your_output_dir
```

MM mode first maps the assembly to the bundled `L1rp.fa` seed with minimap2 to nominate L1-like
intervals, then runs RepeatMasker only on those candidate sequences to assign
their subfamily before continuing with liftover and optional ORF validation.
Candidate-level RepeatMasker runs with `-pa 4` by default and can be adjusted
with `-pa/--pa`.

Output:
- haplongliner_mm.bed file with TE info from your assembly and corresponding reference genome (hs1/hg38) coordinates and ORF status
- haplongliner_mm.fa file containing all full length (>=5kb by default) sequences for the selected L1 families
- Log file (when using `-g/--log`) that summarizes results of each step of the pipeline mode

### SV mode: Based on structural variation calls, no RepeatMasker pre-masking of whole assembly needed

Input:
- Haploid assembly FASTA
- Structural variant (SV) callset (VCF or BED)
- For phased haplotype-resolved SV VCFs, `--hap auto` (default) will try to infer haplotype 1 or 2 from the assembly filename; use `--hap 1`, `--hap 2`, or `--hap both` to override
- Reference genome FASTA (hs1, hg38, or custom; default hs1)
- Transposable element families to analyse (via `--te`, default `L1,L1PA3`)
- Maximum insertion length for candidate extraction (via `-x/--xlength`; default 20kb or 3×`--length`, whichever is higher)
- The built-in pangenome L1 reference was generated by the RM mode using 94 haplotypes of HPRC release 1, ≥5kb cutoff, plus the L1s that were present in the reference genome (hs1 and hg38) but not in any of the 94 HPRC haploids

Command with test genome:
```bash
haplongliner sv \
  --in tests/test.genome.HG00410.1.fa \
  --sv tests/test.genome.HG00410.vcf.gz \
  --teref hprc \
  --ref hs1 \
  --out test_output_sv
```

Command with user provided genome:
```bash
haplongliner sv \
  --in your.genome.fa \
  --sv your.sv.vcf \
  --teref your.te.bed \
  --te L1,SVA \
  --ref you_choose_hs1_or_hg38 \
  --out your_output_dir
```

ORF validation is carried out only when `--length` is at least 5000 and the `--te` list includes L1HS, L1PA2, or L1PA3; otherwise the BLASTP and ORF sequence extraction steps are skipped.
RepeatMasker-based candidate validation runs with `-pa 4` by default and can be
adjusted with `-pa/--pa`. When the SV VCF contains phased genotypes such as
`1|0` or `0|1`, SV mode uses the selected haplotype to keep only the variants
present on the target haploid assembly.

Output:
- haplongliner_sv.bed file summarizing TE coordinates and ORF status. Names from insertion calls that do not overlap lifted TEs end with `;nr`.
- haplongliner_sv.fa file containing all full length (>=5kb by default) sequences for the selected TE families
- Log file (when using `-g/--log`) that summarizes results of each step of the pipeline mode


### Sequence retrieval function

This function retrieves the sequences of L1 insertions present in the HPRC
repository.  The user supplies a coordinate (``chr:start-end``) or a BED file of
coordinates and the function outputs the corresponding L1 sequences.  Additional
FASTA files from other individuals can be provided and will be appended to the
resulting FASTA.
Currently, coordinate lookup is only available for the **hs1** reference genome.
Queries on other references are not yet supported.

Example command with single coordinate:
```bash
haplongliner seq -q chr1:48841802-48841812 -o chr1:48841802-48841812.fa
```

Command using coordinates from BED file and user provided extra L1 sequences:
```bash
haplongliner seq -q sites.bed -e extra.fa -o out.fa
```

Output:
- FASTA file containing all L1 sequences found at the queried coordinates,
  optionally combined with sequences from ``-e/--extra``.


## Authors

Lei Yang, Sara Nematbakhsh, Amanda Norseen, and Rick McLaughlin
