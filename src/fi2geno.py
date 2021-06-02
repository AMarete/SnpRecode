#!/usr/bin/env python
# -*- coding: utf-8 -*-

import timeit
import re
from datetime import datetime
from Bio import bgzf
from funtools import flatten, bomb
from parse_args import my_parser
from recode_dict import fimpute_2_vcf

start = timeit.default_timer()

try:
    samples = open(my_parser().samples, "r")
except FileNotFoundError:
    bomb(f"Missing argument or '{my_parser().samples}' may be empty\n")
try:
    snp_info = open(my_parser().snps, "r")
except FileNotFoundError:
    bomb(f"Missing argument or '{my_parser().snps}' may be empty\n")
try:
    geno_info = open(my_parser().geno, "r")
except FileNotFoundError:
    bomb(f"Missing argument or '{my_parser().geno}' may be empty\n")
try:
    allele_info = open(my_parser().allele, "r")
except FileNotFoundError:
    bomb(f"Missing argument or '{my_parser().allele}' may be empty\n")

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
##source="snprecode v1.0.3"
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
plink_info = {}
count = 0
with snp_info as file:
    next(file)
    for line in file:
        snpid, chrom, pos, chips = re.sub('\s+','\t',line).split('\t', 3)
        snps[snpid.lower()] = [chrom, pos, snpid]
        # snps[chrom + ":" + pos] = [chrom, pos, snpid]
        # snps.append(chrom + ":" + pos)
        plink_info[count] = [snpid.lower()]
        count += 1

# Alleles file
vcf_a1a2 = {}
plink_ACGT = {}
count = 0
with allele_info as file:
    for line in file:
        chrom, snp, cm, pos, a1, a2 = line.strip().split()
        vcf_a1a2[snp.lower()] = [chrom, pos, snp, a1, a2, ".", "PASS", ".", "GT"]
        # vcf_a1a2[chrom + ":" + pos] = [chrom, pos, snp, a1, a2, ".", "PASS", ".", "GT"]
        if a1 == '.': a1 = '0'
        if a2 == '.': a2 = '0'
        if a1 == '0':  # monomorphic snp (A1=0)
            plink_ACGT[snp.lower()] = [a1 + ' ' + a1, a1 + ' ' + a1, a2 + ' ' + a2, a1 + ' ' + a1,
                                       a1 + ' ' + a1, a1 + ' ' + a1, a1 + ' ' + a1, a1 + ' ' + a1,
                                       a1 + ' ' + a1, a1 + ' ' + a1, snp]
        else:
            plink_ACGT[snp.lower()] = [a1 + ' ' + a1, a1 + ' ' + a2, a2 + ' ' + a2, a1 + ' ' + a2,
                                       a2 + ' ' + a1, '0 0', '0 0', '0 0', '0 0', '0 0', snp]
        count += 1

# For VCF only
# Add header for first 9 lines and write vcf
vcfsnps = {k: vcf_a1a2[k] for k in vcf_a1a2.keys() & set(snps.keys())}.values()
vcfsnps = sorted(vcfsnps, key=lambda x: (int(x[0]), int(x[1])))
vcfsnps.insert(0, ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"])
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
                geno_tot.append(plink_ACGT[plink_info[n][0]][int(a)])
            geno_out.write('%s %s \n' % (' '.join([sample_id, sample_id, '0','0','2','-9' ]), ' '.join(geno_tot)))

# continue writing VCF
if toto == 1:
    for row in zip(vcfsnps, [*zip(*genotypes)]):
        line = [list(flatten(i)) for i in list(row)]
        geno_out.write('\t'.join(flatten(line)) + '\n')

# Map file for Ped
if toto == 2:
    map_out = open(fo + ".map", "w")
    map_ped = []
    for snp in snps.keys(): # snps is a list need to ammend this
        try:
            ms = vcf_a1a2[snp][0:3]
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
    tt = fo + '.vcf.gz'
else:
    tt = ' and '.join([fo + ".ped", fo + ".map"])

stop = timeit.default_timer()
mins, secs = divmod(stop - start, 60)
hours, mins = divmod(mins, 60)

print(f"Successful conversion from FImpute format to {tt}")
print("Conversion Time (H:M:S): %d:%d:%d\n" % (hours, mins, secs))
raise SystemExit(0)
