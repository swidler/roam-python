#!/bin/bash
#$ -cwd
#$ -j y
#$ -S bin/bash 

sampname=${1}
vcf=${2}
outdir=${3}
logfile=${outdir}/${sampname}_log.txt

module load bcftools
date +"%d-%m-%Y %H:%M:%S" > ${logfile}

echo -e "\n################################################################################" >> ${logfile}
echo -e  "\tGetting positions to add and remove for " ${sampname} >> ${logfile}
echo -e "################################################################################\n" >> ${logfile}

# Get all CpG positions that are actually present in the individual
echo -e "\tGetting all CpG positions" >> ${logfile}
bcftools query -f '%CHROM\t%POS\t[%TGT]\n' ${vcf} | awk '$3 ~ "C" {print $1 "\t" $2+1}' | bcftools query -f '%CHROM\t%POS0\t[%TGT]\n' -T - ${vcf} | awk '$3 ~ "G" {print $1 "\t" $2}' > ${outdir}/${sampname}_CpGs.txt

# Get all reference CpG positions that are covered by VCF
echo -e "\tGetting all covered reference CpG positions" >> ${logfile}
bcftools query -i 'REF="C"' -f '%CHROM\t%POS\n' ${vcf} | awk '{print $1 "\t" $2+1}' | bcftools query -i 'REF="G"' -T - -f '%CHROM\t%POS0\n' ${vcf} > ${outdir}/${sampname}_ref_CpGs.txt

echo -e "\tNow checking for overlap...\n"
comm --output-delimiter=, <(sort ${outdir}/${sampname}_CpGs.txt) <(sort ${outdir}/${sampname}_ref_CpGs.txt) | tee >(awk -F',' '{print $1}' | grep -v '^$' | sort > ${outdir}/${sampname}_to_add.txt) >(awk -F',' '{print $2}' | grep -v '^$' | sort > ${outdir}/${sampname}_to_remove.txt) >(awk -F',' '{print $3}' | grep -v '^$' | echo -e "\tThe number of reference CpG positions that are present in ${sampname} is: `wc -l` " >> ${logfile}) > /dev/null

echo -e "\tThe number of actual CpG positions that are in the genotype is: `wc -l ${outdir}/${sampname}_CpGs.txt | awk '{print $1}'`" >> ${logfile}
echo -e "\tThe number of theoretical (reference) CpG positions that are covered by the VCF is: `wc -l ${outdir}/${sampname}_ref_CpGs.txt | awk '{print $1}'` \n" >> ${logfile}

# The positions that move on to next step
echo -e "The number of positions to add is: `wc -l ${outdir}/${sampname}_to_add.txt | awk '{print $1}'` " >> ${logfile}
echo -e "The number of positions to remove is: `wc -l ${outdir}/${sampname}_to_remove.txt | awk '{print $1}'` \n"  >> ${logfile}

echo "Finished" >> ${logfile}
date +"%d-%m-%Y %H:%M:%S" >> ${logfile}
