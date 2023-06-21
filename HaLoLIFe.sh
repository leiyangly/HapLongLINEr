#Set variables
ind=$1
hap=$2
gd=$3
wd=$4
pd=$5
rg=$6


#Extract the coordinate of the full length young L1s (L1HS, L1PA2 and L1PA3)
perl ${pd}/scripts/ExtractL1.pl ${gd}/${ind}.${hap}.f1_assembly_v2_genbank_rm.bed > ${wd}/${ind}.${hap}.L1HSPA2PA3.bed

#Extract the sequence of the full length L1s
cat ${wd}/${ind}.${hap}.L1HSPA2PA3.bed |\
perl -lane 'print if $F[5] eq "\+"' |\
seqtk subseq ${gd}/${ind}.${hap}.f1_assembly_v2_genbank.fa.gz - |\
seqtk seq -U -l 0 - > ${wd}/${ind}.${hap}.L1HSPA2PA3.fa

cat ${wd}/${ind}.${hap}.L1HSPA2PA3.bed |\
perl -lane 'print if $F[5] eq "\-"' |\
seqtk subseq ${gd}/${ind}.${hap}.f1_assembly_v2_genbank.fa.gz - |\
seqtk seq -U -r -l 0 - >> ${wd}/${ind}.${hap}.L1HSPA2PA3.fa



#Extract the flanking 2kb regions of L1
cat ${wd}/${ind}.${hap}.L1HSPA2PA3.bed |\
perl -lane 'BEGIN{$"="\t"} $F[2]=$F[1]; $F[1]=$F[1]-2000; print "@F"' |\
seqtk subseq ${gd}/${ind}.${hap}.f1_assembly_v2_genbank.fa.gz - |\
seqtk seq -U -l 0 - > ${wd}/${ind}.${hap}.L1HSPA2PA3-2kb.fa

cat ${wd}/${ind}.${hap}.L1HSPA2PA3.bed |\
perl -lane 'BEGIN{$"="\t"} $F[1]=$F[2]; $F[2]=$F[2]+2000; print "@F"' |\
seqtk subseq ${gd}/${ind}.${hap}.f1_assembly_v2_genbank.fa.gz - |\
seqtk seq -U -l 0 - > ${wd}/${ind}.${hap}.L1HSPA2PA3+2kb.fa

#Use minimap2 to map the flanking regions of L1s to the reference genome
minimap2 -x asm5 ${gd}/${rg}.mmi ${wd}/${ind}.${hap}.L1HSPA2PA3-2kb.fa > ${wd}/${ind}.${hap}.L1HSPA2PA3-2kb.${rg}.minimap.txt

minimap2 -x asm5 ${gd}/${rg}.mmi ${wd}/${ind}.${hap}.L1HSPA2PA3+2kb.fa > ${wd}/${ind}.${hap}.L1HSPA2PA3+2kb.${rg}.minimap.txt

#Find all ORFs in all L1s that are identified by repeat masker
cat ${wd}/${ind}.${hap}.L1HSPA2PA3.fa |\
perl -plne 's/[\:\-]/\_/g if m/^\>/' > ${wd}/${ind}.${hap}.L1HSPA2PA3.rename.fa

getorf -sequence ${wd}/${ind}.${hap}.L1HSPA2PA3.rename.fa -find 1 -outseq ${wd}/${ind}.${hap}.L1HSPA2PA3AllORF.fa

rm ${wd}/${ind}.${hap}.L1HSPA2PA3.rename.fa

perl ${pd}/scripts/ProcessORF.pl ${wd}/${ind}.${hap}.L1HSPA2PA3AllORF.fa > ${wd}/${ind}.${hap}.L1HSPA2PA3AllORF.bed

#BLAST the ORFs found
blastp -db ${pd}/L1ref/L1rpORF12p.fa -query ${wd}/${ind}.${hap}.L1HSPA2PA3AllORF.fa -outfmt "6 std qlen slen sacc" > ${wd}/${ind}.${hap}.L1HSPA2PA3AllORF.blastp


#Find the longest alignment to ORF1 and ORF2 among all found ORFs
perl ${pd}/scripts/FindLongestORF.pl ${wd}/${ind}.${hap}.L1HSPA2PA3AllORF.blastp > ${wd}/${ind}.${hap}.L1HSPA2PA3AllORF.combine.blastp

#Find L1s with intact ORFs:
perl ${pd}/scripts/FindIntactORF.pl ${wd}/${ind}.${hap}.L1HSPA2PA3AllORF.combine.blastp > ${wd}/${ind}.${hap}.L1HSPA2PA3AllORF.intact.blastp


#Integrate the ORF finding and liftover information.
#(1)only the best hit in the minimap2 alignments for the flanking 2kb are taken
#(2)only the minimap2 alignments that are longer than 200bp are taken. I tried various lengths, but setting this too long will make some of the regions unable to align. I think for the purpose of finding the rough liftover position, using short cutoffs makes sense. Setting this number to even ~1.4kb will not change the results for intact L1s of HG02622, but we don't know what will happen when we investigate more genomes
#(3)only take the coordinates of minimap mapping when the plus and minus strand lifted over to the similar chromosomal locations.
#(4)prefer to take the alignments on chromosomes instead of the decoys.
perl ${pd}/scripts/CombineTable.pl ${wd}/${ind}.${hap}.L1HSPA2PA3+2kb.${rg}.minimap.txt ${wd}/${ind}.${hap}.L1HSPA2PA3-2kb.${rg}.minimap.txt ${wd}/${ind}.${hap}.L1HSPA2PA3AllORF.intact.blastp ${wd}/${ind}.${hap}.L1HSPA2PA3.bed > ${wd}/${ind}.${hap}.HaLoLIFe.output.txt


#Retrieve the sequence of all the L1s with intact ORFs:
cat ${wd}/${ind}.${hap}.HaLoLIFe.output.txt |\
perl -lane '@G=split/\_/,$F[0]; print "$G[0]\t$G[1]\t$G[2]" if $G[3] eq "\+" and $G[6] eq "intact"' |\
seqtk subseq ${gd}/${ind}.${hap}.f1_assembly_v2_genbank.fa.gz - |\
seqtk seq -U -l 0 - > ${wd}/${ind}.${hap}.L1HSPA2PA3AllORF.intact.fa

cat ${wd}/${ind}.${hap}.HaLoLIFe.output.txt |\
perl -lane '@G=split/\_/,$F[0]; print "$G[0]\t$G[1]\t$G[2]" if $G[3] eq "\-" and $G[6] eq "intact"' |\
seqtk subseq ${gd}/${ind}.${hap}.f1_assembly_v2_genbank.fa.gz - |\
seqtk seq -U -r -l 0 - >> ${wd}/${ind}.${hap}.L1HSPA2PA3AllORF.intact.fa



