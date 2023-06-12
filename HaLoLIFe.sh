#Extract the coordinate of the full length young L1s (L1HS, L1PA2 and L1PA3)
perl ExtractL1.pl $2 > L1HSPA2PA3.bed

#Extract the sequence of the full length L1s
cat L1HSPA2PA3.bed |\
perl -lane 'print if $F[5] eq "\+"' |\
seqtk subseq $1 - |\
seqtk seq -U -l 0 - > L1HSPA2PA3.fa

cat ~/HPRC_L1/HG02622.maternal.L1HSPA2PA3.bed |\
perl -lane 'print if $F[5] eq "\-"' |\
seqtk subseq $1 - |\
seqtk seq -U -r -l 0 - >> L1HSPA2PA3.fa

#Extract the flanking 2kb regions of L1
cat L1HSPA2PA3.bed |\
perl -lane 'BEGIN{$"="\t"} $F[2]=$F[1]; $F[1]=$F[1]-2000; print "@F"' |\
seqtk subseq $1 - |\
seqtk seq -U -l 0 - > L1HSPA2PA3-2kb.fa

cat L1HSPA2PA3.bed |\
perl -lane 'BEGIN{$"="\t"} $F[1]=$F[2]; $F[2]=$F[2]+2000; print "@F"' |\
seqtk subseq $1 - |\
seqtk seq -U -l 0 - > L1HSPA2PA3+2kb.fa

#Use minimap2 to map the flanking regions of L1s to the reference genome
minimap2 -x asm10 -d ${3}.mmi $3

minimap2 -x asm10 ${3}.mmi L1HSPA2PA3-2kb.fa > L1HSPA2PA3-2kb.minimap.txt

minimap2 -x asm10 ${3}.mmi L1HSPA2PA3+2kb.fa > L1HSPA2PA3+2kb.minimap.txt

#Find all ORFs in all L1s that are identified by repeat masker
perl -plne 's/[\:\-]/\_/g if m/^\>/' L1HSPA2PA3.fa > temp009

getorf -sequence temp009 -find 1 -outseq L1HSPA2PA3AllORF.fa

rm temp009

perl ProcessORF.pl L1HSPA2PA3AllORF.fa > L1HSPA2PA3AllORF.bed

#BLAST the ORFs found
blastp -db ./L1ref/L1rpORF12p.fa -query L1HSPA2PA3AllORF.fa -outfmt "6 std qlen slen sacc" > L1HSPA2PA3AllORF.blastp


#Find the longest alignment to ORF1 and ORF2 among all found ORFs
perl FindLongestORF.pl L1HSPA2PA3AllORF.blastp > L1HSPA2PA3AllORF.combine.blastp

#Find L1s with intact ORFs:
perl FindIntactORF.pl L1HSPA2PA3AllORF.combine.blastp > L1HSPA2PA3AllORF.intact.blastp


#Integrate the ORF finding and liftover information.
#(1)only the best hit in the minimap2 alignments for the flanking 2kb are taken
#(2)only the minimap2 alignments that are longer than 200bp are taken. I tried various lengths, but setting this too long will make some of the regions unable to align. I think for the purpose of finding the rough liftover position, using short cutoffs makes sense. Setting this number to even ~1.4kb will not change the results for intact L1s of HG02622, but we don't know what will happen when we investigate more genomes
#(3)only take the coordinates of minimap mapping when the plus and minus strand lifted over to the similar chromosomal locations.
#(4)prefer to take the alignments on chromosomes instead of the decoys.
perl CombineTables.pl L1HSPA2PA3+2kb.hg38.txt L1HSPA2PA3-2kb.hg38.txt L1HSPA2PA3AllORF.intact.blastp L1HSPA2PA3.bed > output.txt


