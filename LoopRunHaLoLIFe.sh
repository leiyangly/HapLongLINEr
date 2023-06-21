#Set the genome directory
gd=~/Downloads/HPRC
#Set the working directory
wd=~/Downloads/HPRC_L1
#Set HaLoLIFe's program directory
pd=~/Programs/HaLoLIFe
#Set the reference genome name
rg=hg38.fa
#List of HPRC individuals to query
HPRC=HG01891\ HG02486\ HG02559\ HG02257\ HG01358\ HG01123\ HG01361\ HG01258\ HG03516\ HG02572\ HG02886\ HG02717\ HG02630\ HG02622\ HG03540\ HG03453\ HG03579\ HG01978\ HG01928\ HG02148\ HG01952\ HG01106\ HG01175\ HG00741\ HG00735\ HG01071\ HG00621\ HG00438\ HG00673

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



