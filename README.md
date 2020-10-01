author: andrew.marete@canada.ca, (C), 2020

snprecode is a helper utility to convert between various genotype formats specifically for fimpute software<br/>
Fimpute is a powerful software that accepts and generates genotypes in specific format.<br/> 
Format conversion software does not come standard with the distribution of fimpute executable<br/>
This software aims to bridge this gap by:<br/> 

    simulatneously seemless conversion of upto 10 chips from vcf and/or ped/map to fimpute acceptable format
    Conversion from fimpute to vcf or ped/map format
    Calculation of alleclic correlation between two vcf files

snp recode is written in python 3.5 and tested on:

    Operating System: CentOS Linux 7 (Core)
    CPE OS Name: cpe:/o:centos:centos:7
    Kernel: Linux 3.10.0-957.21.3.el7.x86_64
    Architecture: x86-64

For snprecode to run properly, install Python v3.5+ and Biopython module

To use snprecode, download the standalone from:<br/>
```https://github.com/AMarete/FImpute-Imputation/blob/master/snprecode```<br/>
change to executable: ```chmod a+x snprecode```<br/>
for full argument list run ```./snprecode.py -h```

basic usage:  
    
    Convert to fimpute format: 
    ./snprecode \
    -F /path/to/vcf/plink/files \
    -O PREFIX 
    
    Convert from fimpute format:
    ./snprecode \
    -g genotypes_imp.txt \
    -s snps_info.txt \
    -n study_samples.txt \ 
    -t 1=vcf; 2=ped|map \ 
    -a alleles.txt \
    -o PREFIX


On occasion plink tends to flip alleles, its therefore best to use a normalized vcf file as follows:
convert plink to vcf:

	plink —{species} --file {file} —recode vcf-iid —out {file}

compress with bgzip:<br/> 
  ```bgzip {file}```

normalize the vcf with the reference fasta genome, bcftools is a great tool for this:
if working with sequence data (>1m snp/chrom) there is need to split data per chromosome, the -r flag can be used

	bcftools norm -d snps -cs -f {ref_genome.fa} -r {chrom} {file.vcf.gz} | bcftools view -m2 -M2 -Oz -o {file2.vcf.gz} 

However if one really needs to use plink, then the procedure below could be followed:<br/>
convert the plink file to binary format<br/>
change the snp name to chrom_pos_ref_alt<br/>
use option —recode 12 to generate plink files with alleles recoded as 1|2<br/>

	plink —{species} —file {dt} —make-bed --real-ref-alleles —out {dt1}
	echo "$(awk '{$2=$1"_"$4"_"$5"_"$6 ;print $0 }' dt1.bim)" > dt1.bim
	plink —bile {dt1} —recode 12 —out {dt2}


add the vcf files to a temp directory e.g. assuming {file2.vcf.gz}, {dt2.ped, dt2.map} are from two populations<br/>
	```mkdir tmp```<br/>
	```mv dt2* tmp/```<br/>
	```mv file2.vcf.gz tmp/```<br/>
	```mv {other_ready_files} tmp```/<br/>

run snprecode<br/>
	```./snprecode -F tmp -O {PREFIX} ```

the program will stop if:<br/>
    ```- duplicate snp(s) or sample(s) found within or between files```<br/>
    ```- differing number of chromosomes between files```

If such errors are found, a file with errors will be generated: ```Error.txt```<br/>

If no errors are found, snprecode generates 3 files: <br/>
	```{PREFIX}.geno, {PREFIX}.mark for fimpute input and``` <br/>
	```{PREFIX}.alleles to decode fimpute output```<br/>

run fimpute per fimpute guidelines<br/>
```http://animalbiosciences.uoguelp.ca/~msargol/fimpute/FImpute_documentation.pdf```<br/>For phased output, use option ```save haplotypes; diplotype```

recode from Fimpute format to vcf<br/>
	```./snprecode -g genotypes_impute.txt -s snp_info.txt -n {samples.txt} -t 1 -a  {PREFIX}.alleles -o {PREFIX_2}```

A final vcf with study samples will be produced for downstream analysis i.e. ```{PREFIX_2}.vcf.gz```<br/>
If fimpute was successful in phasing, the alleles will be phased with | separator, otherwise they’ll have a / separator as per vcf conventions 


VCF correlation :: coming soon!
