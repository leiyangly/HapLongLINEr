HapLongLINEr

Description

HapLongLINEr is a pipeline for identifying young LINE-1 elements (L1HS, L1PA2, and L1PA3) with intact open reading frames (ORFs) in haploid long-read human genome assemblies. The pipeline also converts (lifts over) the coordinates of these LINE-1s to a reference genome. The name stands for Haploid Long read assembly-based LINE-1 Retriever.

HapLongLINEr is designed for use with haploid long-read assemblies, as these provide much better recovery of young LINE-1 elements compared to diploid or fragmented assemblies. The current version is compatible with the file and directory structures used by the HPRC and 1000G_ONT projects.

Dependencies

bash
perl
seqtk
minimap2
NCBI BLAST+
getorf (EMBOSS)
Authors

Lei Yang and Amanda Norseen

Usage

To run the main pipeline, use:

HapLongLINEr.sh your.genome.fa repeatmasker.bed reference.genome.fa
Output

Files ending with .HapLongLINEr.output.txt:
Two-column, tab-delimited text files.
First column: LINE-1 info from your assembly
Second column: Corresponding info in hg38
Each field contains: chromosome, coordinate, strand, intactness, and L1 family, separated by underscores.
Files ending with .L1HSPA2PA3AllORF.intact.fa:
FASTA file containing the sequences of all intact L1HS, L1PA2, and L1PA3 elements from the input assembly.
Running on an HPRC Individual

Edit the first 10 lines of LoopRunHapLongLINEr.sh to set directories and input files, then run:

sh LoopRunHapLongLINEr.sh
