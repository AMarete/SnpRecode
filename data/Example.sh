#!/bin/bash
set -e
set -o pipefail

# andrew.marete@canada.ca, 2018

# Check pre-requisites
if  ! command -v snprecode &> /dev/null; then echo "SnpRecode not found"; exit;  fi
if  ! command -v bcftools &> /dev/null; then echo "bcftools not found"; exit;  fi
if  ! command -v plink &> /dev/null; then echo "plink not found"; exit;  fi

# Directories
path0=/home/$USER
path1=${path0}/results
tmp0=${path0}/tmp$RANDOM
tmp1=${tmp0}/tmp$RANDOM

if [[ ! -d $path1 ]]; then mkdir -p $path1 ; else rm -rf $path1 ; fi
if [[ ! -d $tmp0 ]]; then mkdir -p $tmp0 ; fi
if [[ ! -d $tmp1 ]]; then mkdir -p $tmp1 ; fi

cd $tmp0

# Process files from upto 10 chips
ln -s ref.vcf.gz .
ln -s val1.vcf.gz .
ln -s val2.ped .
ln -s val2.map .

# optional -- convert plink to vcf
plink —{species} --file val2 —recode vcf-iid —out val2
bgzip val2.vcf && bcftools index val2.vcf.gz

# Normalize the data
bcftools norm -d snps -cs -f $ref_genome -r $bta -Ou val1.vcf.gz | bcftools view -m2 -M2 -Oz -o ${tmp1}/dt1.vcf.gz &
bcftools norm -d snps -cs -f $ref_genome -r $bta -Ou val2.vcf.gz | bcftools view -m2 -M2 -Oz -o ${tmp1}/dt2.vcf.gz &
bcftools norm -d snps -cs -f $ref_genome -r $bta -Ou ref.vcf.gz | bcftools view -m2 -M2 -Oz -o REF.vcf.gz &
wait

# select some reference animals to mask
pc=0.3
to_mask=$(bcftools query -l REF.vcf.gz | wc -l | awk -v y=$pc '{print int($1 * y)}')
bcftools query -l REF.vcf.gz | shuf -n $to_mask > ref.masked.txt

# samples list
bcftools query -l ${tmp1}/dt1.vcf.gz > dt1.samples
bcftools query -l ${tmp1}/dt2.vcf.gz > dt2.samples
cat dt1.samples dt2.samples > val_samples.txt
cat val_samples.txt ref.masked.txt > val_ref_mask_samples.txt

# Compare SNP between chips
bcftools query -f '%CHROM\t%POS\n' REF.vcf.gz > ref.snp.txt &
bcftools query -f '%CHROM\t%POS\n' ${tmp1}/dt1.vcf.gz > dt12.snp.txt &
bcftools query -f '%CHROM\t%POS\n' ${tmp1}/dt2.vcf.gz >> dt12.snp.txt &
wait 

# remove duplicates
echo "$(awk '!seen[$x]++' dt12.snp.txt)" > dt12.snp.txt

# subset marker common to both val and ref
awk -F '\t' 'NR==FNR {id[$1"_"$2]; next} $1"_"$2 in id' dt12.snp.txt  ref.snp.txt | awk '{print $1"\t"$2}' >  ref.snp.masked.txt

# split the REF dataset:
# 1. does not have the samples whose genotypes is to be masked ie the new REF
# 2. contains the randomly selected subset of ref samples and marker available in the validation data
bcftools index REF.vcf.gz 
bcftools view -S ^ref.masked.txt -Oz -o ${tmp1}/ref.vcf.gz  REF.vcf.gz &
bcftools view -S ref.masked.txt  -R ref.snp.masked.txt -Oz -o ${tmp1}/ref.masked.vcf.gz REF.vcf.gz &
wait

# convert to fimpute format using SnpRecode
# $tmp1 has four files: ref.vcf.gz, ref.masked.vcf.gz, dt1.vcf.gz, and dt2.vcf.gz (or if not normalized val2.ped and val2.ped)
# SnpRecode generates three files: [tot.geno, tot.mark] for fimpute, and [tot.vcf_a1a2] to decode fimpute output
./snprecode -D $tmp1 -O tot

# impute
if [[ ! -d output ]]; then mkdir -p output ; else rm -f output/* ; fi
cat > par.txt <<EOF
title="Population Imputation";
genotype_file="tot.geno";
snp_info_file="tot.mark";
output_folder="output";
njob = 5;
EOF

./FImpute par.txt

# convert back to vcf
./snprecode \
-g output/genotypes_imp.txt \
-s output/snp_info.txt \
-n val_ref_mask_samples.txt \
-o bis \
-t 1 \
-a tot.allele

bcftools index bis.vcf.gz

# Imputation Accuracy
bcftools view -S ref.masked.txt -Oz -o true.vcf.gz  REF.vcf.gz &
bcftools view -S ref.masked.txt -Oz -o mask.vcf.gz  bis.vcf.gz &
wait

./snprecode --file true.vcf.gz mask.vcf.gz
mv genotype_R2.* $path1

# Extract study samples
bcftools view -S val_samples.txt -Oz -o ${path1}/imputed.vcf.gz bis.vcf.gz

cd $path0
rm -rf $tmp0

# $path1 should have three files: imputed_file, genotype_correlation files [.txt and .pdf]
