#Set the genome directory
gd=~/Downloads/HPRC
#Set the working directory
wd=~/Downloads/HPRC_L1
#Set HaLoLIFe's program directory
pd=~/Programs/HaLoLIFe
#Set the reference genome name
rg=hg38.fa
#List of HPRC individuals to query
HPRC=HG02622\ HG02630

#Index the reference genome
minimap2 -x asm10 -d ${gd}/${rg}.mmi ${gd}/${rg}



#Download the genomes and repeatmasker annotations from HPRC
for ind in ${HPRC}
do
for hap in maternal paternal
do
wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/${ind}/assemblies/year1_f1_assembly_v2_genbank/${ind}.${hap}.f1_assembly_v2_genbank.fa.gz -P ${gd}
wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/${ind}/assemblies/year1_f1_assembly_v2_genbank/annotation/repeat_masker/${ind}.${hap}.f1_assembly_v2_genbank_rm.bed -P ${gd}
done
done



#Run the intact ORF calling script
for ind in ${HPRC}
do
for hap in maternal paternal
do
sh ~/Programs/HaLoLIFe/HaLoLIFe.sh ${ind} ${hap} ${gd} ${wd} ${pd} ${rg}
done
done


#Remove the genomes
for ind in ${HPRC}
do
for hap in maternal paternal
do
rm ${gd}/${ind}.${hap}.f1_assembly_v2_genbank.fa.gz
rm ${gd}/${ind}.${hap}.f1_assembly_v2_genbank_rm.bed
done
done



