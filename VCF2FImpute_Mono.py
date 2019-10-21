#!/usr/bin/env python
import sys
import gzip
import time
import datetime
from optparse import OptionParser

start = time.time()
usage = "usage: %prog [options] <args>"
parser = OptionParser()
#parser.add_option("-g", "--vcf", action="store", dest="vcf", help="VCf file")

parser.add_option("-o", "--geno", action="store", dest="geno", help="FImpute genotype file")
parser.add_option("-m", "--snp", action="store", dest="snp", help="FImpute snp_info file")

try:
    (options, args) = parser.parse_args()
except TypeError:
    print("Make sure input file is specified")
    sys.exit(2)

print("-----------------------------------------------------")
print("-   Convert VCF V4.* to FImpute version 3.0 format  -")
print("-   Author:  andrew.marete@canada.ca                -")
print("-----------------------------------------------------")


# Some Functions: convert the list into 2D list


def chunks(l, n):
    # For item i in a range that is a length of l,
    for i in range(0, len(l), n):
        # Create an index range for l of n items:
        yield l[i:i + n]


# create a recode dictionary for VCF to FImpute notation
# recode_genotyped = {'0/0': '0', '1/1': '2', '1/0': '1', '0/1': '1', './.': '5'}
# recode_phased = {'0|0': '0', '1|1': '2', '1|0': '1', '0|1': '1', './.': '5'}
recode = {'0/0': '0',
          '1/1': '2',
          '1/0': '1',
          '0/1': '1',
          '.|.': '5',
          '0|0': '0',
          '1|1': '2',
          '1|0': '1',
          '0|1': '1',
          './.': '5'}



# Open and read VCF file
#in_vcf = gzip.open(options.vcf, "rt")

# Open output files for saving
FoutGenot = open(options.geno, "w")
FoutMark = open(options.snp, "w")

# SNP file for FImpute
# has 4 fields: snp_name, bta, pos, hd_snp_index ld_snp_index (starting from 1)
# snp_index = 1
# FoutMark.write('SNP_ID Chr Pos Chip1\n')

# Treat the Reference VCF genotype file
geno = []
snp = []

#for line in in_vcf:
for line in sys.stdin:
    # get list of animals
    if "#CHR" in line:
        anim_ref = line.strip().split("\t")[9:]
    # get DNA as 2D list
    # get SNP_info list
    if "#" in line: continue
    geno.append(line.strip().split("\t")[9:])
    bta = line.strip().split("\t")[0].replace("Chr", "")
    pos = line.strip().split("\t")[1]
    snp_key = bta + ":" + pos + " " + bta + " " + pos
    snp.append(snp_key)

# recode to FImpute format
geno_recode = []
for row in geno:
    geno_recode.append([recode.get(row[x], row[x]) for x in range(0, len(row), 1)])

# Initialize an empty dictionary for appending genotypes
geno = {}
# geno.setdefault(key, []).append(value) # append multiple values to same key

# transpose using list comprehension
for snp_row in range(0, len(snp), 1):
    for anim_column in range(0, len(anim_ref), 1):
        key1 = anim_ref[anim_column]
        value1 = geno_recode[snp_row][anim_column]
        try:
            geno[key1].append(value1)
        except KeyError:
            geno[key1] = [value1]

# save the FImpute genotype file
FoutGenot.write('ID Chip Call...\n')
# Reference File
for key, value in geno.items():
    FoutGenot.write('{} {} {}\n'.format(key, 1, ''.join(value)))

# sort the snp lists by chromosome (index = 1) then position (index = 2)
# not necessary since by default, VCF files should be sorted by chrom then pos2
snp2 = []
for item in chunks(snp, 1):
    for ii in item:
        snp2.append(ii.split())

snp2 = sorted(snp2, key=lambda x: (int(x[1]), int(x[2])))

# save the FImpute  snp_info file adding an FImpute required snp index
index = 0
FoutMark.write('snp chrom pos chip\n')
for row in snp2:
    index += 1
    row.insert(3, str(index))
    FoutMark.write(" ".join(row) + '\n')

FoutGenot.close()
FoutMark.close()
