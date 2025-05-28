#!/bin/bash

# HapLongLINErSV.sh
# Module 2: Find Full-Length Young L1s (RepeatMasker-free, SV-based, HPRC anchor)
# Usage: sh HapLongLINErSV.sh <sample_id> <hap> <target_assembly.fa> <work_dir> <pipeline_dir> <hprc_l1_bed> <hprc_hs1_fa>

sample_id=$1
hap=$2
target_assembly=$3
work_dir=$4
pipeline_dir=$5
hprc_l1_bed=$6      # BED file: HPRC L1 coordinates on hs1
hprc_hs1_fa=$7      # hs1 reference FASTA

# Step 1: Extract ±2kb flanking sequences from HPRC L1 coordinates (on hs1)
cat ${hprc_l1_bed} |\
perl -lane 'BEGIN{$"="\t"} $F[2]=$F[1]; $F[1]=$F[1]-2000; print "@F"' > ${work_dir}/HPRC.L1.minus2kb.bed

cat ${hprc_l1_bed} |\
perl -lane 'BEGIN{$"="\t"} $F[1]=$F[2]; $F[2]=$F[2]+2000; print "@F"' > ${work_dir}/HPRC.L1.plus2kb.bed

seqtk subseq ${hprc_hs1_fa} ${work_dir}/HPRC.L1.minus2kb.bed | seqtk seq -U -l 0 - > ${work_dir}/HPRC.L1.minus2kb.fa
seqtk subseq ${hprc_hs1_fa} ${work_dir}/HPRC.L1.plus2kb.bed  | seqtk seq -U -l 0 - > ${work_dir}/HPRC.L1.plus2kb.fa

# Step 2: Map the flanking regions to the target assembly
minimap2 -x asm5 ${target_assembly} ${work_dir}/HPRC.L1.minus2kb.fa > ${work_dir}/${sample_id}.${hap}.L1.minus2kb.minimap.txt
minimap2 -x asm5 ${target_assembly} ${work_dir}/HPRC.L1.plus2kb.fa  > ${work_dir}/${sample_id}.${hap}.L1.plus2kb.minimap.txt

# Step 3: Integrate mapping results to infer L1 locations in the target assembly
perl ${pipeline_dir}/scripts/CombineTable_SV.pl \
    ${work_dir}/${sample_id}.${hap}.L1.plus2kb.minimap.txt \
    ${work_dir}/${sample_id}.${hap}.L1.minus2kb.minimap.txt \
    ${hprc_l1_bed} \
    > ${work_dir}/${sample_id}.${hap}.L1.SVanchors.bed

# Step 4: Extract the putative full-length L1 sequences from the target assembly
seqtk subseq ${target_assembly} ${work_dir}/${sample_id}.${hap}.L1.SVanchors.bed | seqtk seq -U -l 0 - > ${work_dir}/${sample_id}.${hap}.L1.SVanchors.fa

# Step 5: Find ORFs in the extracted L1s
getorf -sequence ${work_dir}/${sample_id}.${hap}.L1.SVanchors.fa -find 1 -outseq ${work_dir}/${sample_id}.${hap}.L1.SVanchors.AllORF.fa

# Step 6: Process ORFs and identify intact L1s
perl ${pipeline_dir}/scripts/ProcessORF.pl ${work_dir}/${sample_id}.${hap}.L1.SVanchors.AllORF.fa > ${work_dir}/${sample_id}.${hap}.L1.SVanchors.AllORF.bed

blastp -db ${pipeline_dir}/L1ref/L1rpORF12p.fa \
    -query ${work_dir}/${sample_id}.${hap}.L1.SVanchors.AllORF.fa \
    -outfmt "6 std qlen slen sacc" \
    > ${work_dir}/${sample_id}.${hap}.L1.SVanchors.AllORF.blastp

perl ${pipeline_dir}/scripts/FindLongestORF.pl ${work_dir}/${sample_id}.${hap}.L1.SVanchors.AllORF.blastp > ${work_dir}/${sample_id}.${hap}.L1.SVanchors.AllORF.combine.blastp

perl ${pipeline_dir}/scripts/FindIntactORF.pl ${work_dir}/${sample_id}.${hap}.L1.SVanchors.AllORF.combine.blastp > ${work_dir}/${sample_id}.${hap}.L1.SVanchors.AllORF.intact.blastp

# Step 7: Extract intact L1 sequences
cat ${work_dir}/${sample_id}.${hap}.L1.SVanchors.AllORF.intact.blastp |\
perl -lane '$"="\t"; print "$F[0]\t$F[1]\t$F[2]" if $F[-1] eq "intact"' |\
seqtk subseq ${target_assembly} - |\
seqtk seq -U -l 0 - > ${work_dir}/${sample_id}.${hap}.L1.SVanchors.intact.fa

echo "SV-based L1 discovery (HPRC anchor) complete for ${sample_id} ${hap}."
