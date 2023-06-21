# HaLoLIFe
Description:
HaLoLIFe is a pipeline that finds young (L1HS, L1PA2 and L1PA3) LINE-1s (L1s) with intact open reading frames based on haploid long read human genome assemblies and performs a "liftOver" by converting the coordinate of the LINE-1s to reference genomes. HaLoLIFe stands for "Haploid Long read assembly-based Intact LINE-1 Finder". We suggest running HaLoLIFe on haploid long read assemblies, because otherwise most of the young LINE-1s will be missing in the assemblies. At this moment, it is written to be compatible with the file name/directory structure compatible with the HPRC project.

Dependencies:
bash, 
perl, 
seqtk, 
minimap2, 
NCBI BLAST+, 
getorf (EMBOSS)

Authors:
Lei Yang

How to run the pipeline:
HaLoLIFe.sh your.genome.fa repeatmasker.bed reference.genome.fa

Output file:
File with "HaLoLIFe.output.txt" at the end contains two columns. The first column describes the L1 info in your assembly, and the second column contains the info in hg38. Within each column contains the information of the chromosome, coordinate, strand, intactness and L1 families, separated by "_". File with "L1HSPA2PA3AllORF.intact.fa" at the end is the fasta file that contains the sequence of all intact L1s from the input genome assembly.

Other utitilties:
...
