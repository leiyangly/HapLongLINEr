# HaLoLIFe
Description
HaLoLIFE is a pipeline that finds LINE-1s (L1s) with intact open reading frames based on haploid long read assemblies and performs a "liftOver" by converting the coordinate of the LINE-1s to reference genomes. It stands for "Haploid Long read assembly-based Intact LINE-1 Finder".

Dependencies
bash
perl
seqtk
minimap2
getorf (EMBOSS)

Authors
Lei Yang

How to run the pipeline
HaLoLIFe.sh your.genome.fa repeatmasker.bed reference.genome.fa

Output file
output.txt contains two column. The first column describes the L1 info in your assembly, and the second column contains the info in hg38. Within each column contains the information of the chromosome, coordinate, strand, intactness and L1 families, separated by "_".

Other utitilties
...
