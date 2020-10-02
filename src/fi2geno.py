#!/usr/bin/env python
# -*- coding: utf-8 -*-

import timeit
from datetime import datetime
from Bio import bgzf
from funtools import flatten, bomb
from parse_args import my_parser
from recode_dict import fimpute_2_vcf

start = timeit.default_timer()

samples = open(my_parser().samples, "r")
snp_info = open(my_parser().snps, "r")
geno_info = open(my_parser().geno, "r")
allele_info = open(my_parser().allele, "r")

toto = my_parser().type_

if my_parser().out:
    pass
else:
    bomb('Missing argument, "-o PREFIX", "--out PREFIX"\n'
         'Error: run `./snprecode -h` for complete arguments list '
         'required to recode from FImpute\n')

fo = my_parser().out
# open headers
if toto == 1:
    geno_out = bgzf.BgzfWriter(fo + ".vcf.gz", "wb")
    # write header
    geno_out.write("".join(
        '''##fileformat=VCFv4.2
##filedate=%s
##source="snprecode v1.0.1"
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
''' % (datetime.today().strftime('%Y%m%d'))
    ))
elif toto == 2:
    geno_out = open(fo + ".ped", "w")
else:
    print(f'\nError!: Missing argument. Specify the recode type: `-t 1` (for VCF) or `-t 2` (for PED/MAP)]')
    print('run `./snprecode -h` for complete arguments list required to recode from FImpute\n')
    raise SystemExit

# get samples files
study_samples = []
with samples as file:
    for line in file:
        study_samples.append(line.strip())

# Fimpute snp info
snps = {}
with snp_info as file:
    next(file)
    for line in file:
        snpid, chrom, pos, chips = line.strip().split('\t', 3)
        snps[chrom + ":" + pos] = [chrom, pos, snpid]

# Alleles file
alleles = {}
with allele_info as file:
    for line in file:
        chrom = line.strip().split("_")[0]
        pos = line.strip().split("_")[1]
        ref = line.strip().split("_")[2]
        alt = line.strip().split("_")[3]
        alleles[chrom + ":" + pos] = [chrom, pos, line.strip(), ref, alt, ".", "PASS", ".", "GT"]

# For VCF only
# Add header for first 9 lines and write vcf
if set(snps).issubset(set(alleles)):
    vcf_snps = [v for v in alleles.values()]
vcf_snps = sorted(vcf_snps, key=lambda x: (int(x[0]), int(x[1])))
vcf_snps.insert(0, ["#CHROM", "POS", "sample_id", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"])
genotypes = []

# imputed geno file
with geno_info as gfile:
    next(gfile)
    for line in gfile:
        geno_tot = []
        sample_id, chip, geno = line.strip().split()
        if sample_id not in study_samples: continue
        if toto == 1:
            row = []
            for n, a in enumerate(geno):
                row.append(fimpute_2_vcf.get(str(a), "!"))
            row.insert(0, sample_id)
            genotypes.append(row)
        elif toto == 2:
            for n, a in enumerate(geno):
                try:
                    geno_ACGT = {0: '1 1', 1: '1 2', 2: '2 2', 5: '0 0'}
                    geno_tot.append(geno_ACGT[int(a)])
                except KeyError:
                    geno_ACGT = {0: '1 1', 1: '1 2', 2: '2 2', 3: '1 2', 4: '2 1', 5: '0 0',
                                 6: '0 0', 7: '0 0', 8: '0 0', 9: '0 0'}
                    geno_tot.append(geno_ACGT[int(a)])
            geno_out.write('%s %s %s \n' % (''.join(sample_id), ''.join(sample_id), ' '.join(geno_tot)))

# continue writing VCF
if toto == 1:
    for row in zip(vcf_snps, zip(*genotypes)):
        line = [list(flatten(i)) for i in list(row)]
        geno_out.write('\t'.join(flatten(line)) + '\n')

# Map file for Ped
if toto == 2:
    map_out = open(fo + ".map", "w")
    '''
    if my_parser().map:
        map_out = open(my_parser().out + ".map", "w")
    else:
        bomb('Missing argument, "-m PREFIX", "--map PREFIX"\n'
             'Error: run `./snprecode -h` for complete arguments list '
             'required to recode from FImpute to PED/MAP\n')
             '''
    map_ped = []
    for snp in snps.keys():
        try:
            ms = alleles[snp][0:3]
            ms.insert(0, '0')
            ms = [ms[i] for i in [1, 3, 0, 2]]
            map_ped.append(ms)
        except KeyError:
            ms = snps[snp]
            ms.insert(0, '0')
            ms = [ms[i] for i in [1, 3, 0, 2]]
            map_ped.append(ms)
    map_ped = sorted(map_ped, key=lambda x: (int(x[0]), int(x[3])))
    with map_out as file:
        for row in map_ped:
            file.write('\t'.join(row) + '\n')

geno_out.close()

if toto == 1:
    tt = fo+'.vcf.gz'
else:
    tt = ' and '.join([fo + ".ped", fo + ".map"])

stop = timeit.default_timer()
mins, secs = divmod(stop - start, 60)
hours, mins = divmod(mins, 60)

print(f"Successful conversion from FImpute format to {tt}")
print("Conversion Time (H:M:S): %d:%d:%d\n" % (hours, mins, secs))
raise SystemExit(0)
