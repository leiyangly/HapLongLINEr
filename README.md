*HapLongLINEr

**Description**

HapLongLINEr is a modular pipeline for discovering and curating full-length (≥5 kb) young LINE-1 elements (L1HS, L1PA2, and intact L1PA3) and marking their intact open reading frames (ORFs) in haploid long-read human genome assemblies. The pipeline is designed for maximum recovery and accuracy of full-length LINE-1s, supporting both RepeatMasker-based and RepeatMasker-free approaches. It also provides a sequence repository for pangenome-level reference L1 elements.

Modules

1. Find Full-Length Young L1s (RepeatMasker-based)
Identifies full-length, young L1 elements in haploid long-read assemblies using RepeatMasker-masked input.
Input:
Haploid assembly FASTA
RepeatMasker BED file
Reference genome FASTA

2. Find Full-Length Young L1s (RepeatMasker-free, SV-based)
Discovers full-length, young L1s in haploid long-read assemblies using structural variant (SV) calls and a pangenome-level L1 reference from the Human Pangenome Reference Consortium (HPRC), without requiring RepeatMasker masking.
Input:
Haploid assembly FASTA
Structural variant (SV) callset
Pangenome-level L1 reference FASTA

3. Sequence Repository of Pangenome-Level Reference L1s
Provides a curated FASTA repository of all young, full-length L1 sequences identified at the pangenome level, suitable for downstream analysis and benchmarking.

**Dependencies**

bash
perl
seqtk
minimap2
NCBI BLAST+
getorf (EMBOSS)
Authors

Lei Yang, Amanda Norseen and Rick McLaughlin

**Usage**

Module 1 (RepeatMasker-based):
HapLongLINEr.sh your.genome.fa repeatmasker.bed reference.genome.fa

Module 2 (RepeatMasker-free, SV-based):
Documentation and example command coming soon.

Module 3 (Sequence repository):
See the repository folder or FASTA files with curated pangenome-level L1s.

**Output**

Files ending with .HapLongLINEr.output.txt:
Two-column, tab-delimited text files with L1 info from your assembly and corresponding hg38 coordinates.
Files ending with .L1HSPA2PA3AllORF.intact.fa:
FASTA file containing all intact L1HS, L1PA2, and L1PA3 sequences from the input assembly.
Additional curated FASTA files for the pangenome reference set.
Running on an HPRC Individual

Edit the first 10 lines of LoopRunHapLongLINEr.sh to set directories and input files, then run:

sh LoopRunHapLongLINEr.sh
