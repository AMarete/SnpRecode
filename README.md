author: gitahimart@gmail.com, (C), 2017

<b>SnpRecode:</b> is a helper utility to convert between various genotype formats specifically for fimpute software. Fimpute is a large scale genotype imputation software designed for use where hundreds of thousands of individuals are genotyped with different chips. Fimpute accepts and generates genotypes in specific format. Format conversion software does not come standard with the distribution of fimpute executable. SnpRecode software aims to bridge this gap by:<br/> 

    – simulatneously seemless conversion of upto 10 chips from vcf and/or ped/map to fimpute acceptable format
    – conversion from fimpute to vcf or ped/map format
    – calculation of genotype correlation between two vcf files, presumably (but not limited to) a real.vcf and masked.vcf
    - create an Fimpute acceptable snp info file from one or several files

SnpRecode is written in Python 3.6 and tested on Linux Architecture: ```x86-64```. It can probably work on other Linux-like machines including MacOS. Other modules required include Biopython and Matplotlib which can be installed using ```pip biopython matplotlib```. To use [SnpRecode](https://github.com/AMarete/fimpute-utils/raw/master/snprecode), click on the link to download the standalone then change to executable: ```chmod a+x snprecode```; for full argument list run ```./snprecode -h```

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
    
    Calculate genotype correlation between two (un)compressed vcf
    ./snprecode --file file1.vcf.gz file2.vcf.gz
    
    Create an fimpute acceptable snp_info file from one or more plink maps (bim and/or map), max=10 files
    ./snprecode --file [file_1.map, ..., file_n.bim]

<b>Example workflow:</b> On occasion plink tends to flip alleles, its therefore best to use a normalized vcf file as follows:<br/>
– convert plink to vcf: ```plink —{species} --file {file} —recode vcf-iid —out {file}```<br/>
– compress with bgzip and create a vcf index: ```bgzip {file} && bcftools index {file}```<br/>
– normalize the vcf with the reference fasta genome, [bcftools](https://samtools.github.io/bcftools/) is a great tool for this. If working with sequence data (>1million snp/chrom) there is need to split data per chromosome, the -r flag can be used: <br/>

	bcftools norm -d snps -cs -f {ref_genome.fa} -r {chrom} {file.vcf.gz} | bcftools view -m2 -M2 -Oz -o {file2.vcf.gz} 

However if one really needs to use plink, then the procedure below could be followed:<br/>
– convert the plink file to binary format<br/>
– change the snp name to chrom_pos_ref_alt<br/>
– use option —recode 12 to generate plink files with alleles recoded as 1|2<br/>

	plink --{species} --file {file} --make-bed --real-ref-alleles -—out {dt1}
	echo "$(awk '{$2=$1"_"$4"_"$5"_"$6 ;print $0 }' dt1.bim)" > dt1.bim
	plink --bfile {dt1} --recode 12 -—out {dt2}

add the files to a temp directory e.g. assuming {file2.vcf.gz}, {dt2.ped, dt2.map} are from two populations<br/>
	```mkdir tmp```<br/>
	```mv dt2* tmp/```<br/>
	```mv file2.vcf.gz tmp/```<br/>
	```mv {other_ready_files} tmp```/<br/>

run SnpRecode: 	```./snprecode -D tmp -O {prefix} ```<br/>

the program will stop if:<br/>
    ```– duplicate snp(s) or sample(s) found within or between files```<br/>
    ```– files have different chromosome count and/or different chromosome names```

If such errors are found, a file with errors will be generated: ```Error.txt```. If no errors are found, SnpRecode generates 3 files: <br/>
```{prefix}.geno, {prefix}.mark ```as inputs for fimpute and ```alleles.txt ``` to decode fimpute output.<br/>

run fimpute per [fimpute guidelines](https://animalbiosciences.uoguelph.ca/~msargol/fimpute/FImpute_documentation.pdf)<br/>

<b>Recode back to vcf</b>:```./snprecode -g geno_impute.txt -s snp_info.txt -n {samples.txt} -t 1 -a  alleles.txt -o {prefix2}```. <br/>If the ```-t 2``` switch is used, a plink ped/map file will be generated. A final file with study samples will be produced for downstream analysis e.g ```{prefix2}.vcf.gz```. If fimpute was successful in phasing, the alleles will be phased with ```|``` separator, otherwise they’ll have a ```/``` separator as per vcf conventions. Sorting alphanumeric chromosomes is currently unsupported and all chromosmes should be recoded to numeric prior to using SnpRecode.<br/>

<b>Genotype correlation</b>: To estimate the imputation accuracy, we have implemented a dosage correlation estimator between a real and masked genotype. One can mask a percentage of the reference population pre-imputation, impute it as if the alleles were missing, and compare the now imputed masked genotype to the original genotype. Such a comparison can be achieved by running ```./snprecode --file file1.vcf.gz file2.vcf.gz```. The two files have to have samples occuring in the same sequence otherwise SnpRecode will stop running. A successful run produce two files:<br/>

    – genotype_R2.txt: contains the R-square values per SNP
    – genotype_R2.pdf: a graph showing distribution of minor allele frequency by R-square as shown below:


<p align="center">
  <img width="474" height="421" src="https://github.com/AMarete/fimpute-utils/blob/master/data/genotype_R2.png">
</p>

<b>Plink snp info files</b>: If one already has fimpute acceptable genotype files from one or more files, and they want to generate an fimpute acceptable snp info file corresponding to the genotypes, the option ```./snprecode --file file1.map [file2.map ... filen.map]``` can be used to create such a snp info file. Upto 10 files are allowed since fimpute software can only impute 10 chips at a time. Both bim and ped files are acceptable.

[<b>Example implementation script</b>](https://raw.githubusercontent.com/AMarete/fimpute-utils/master/data/Example.sh)


