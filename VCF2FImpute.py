#!/usr/bin/env python
import sys, os, gzip
from optparse import OptionParser

usage = "usage: %prog [options] <args>"
parser = OptionParser()

parser.add_option("--g", action="store", dest="genot", help="Input VCf file")

parser.add_option("--og", action="store", dest="outputGenot", help="Output FImpute genotype file")
parser.add_option("--om", action="store", dest="outputMark", help="Output FImpute SNP info file")

try:
    (options, args) = parser.parse_args()
except TypeError:
    print("Make sure input file is specified")
    sys.exit(2)


print("-----------------------------------------------------")
print("-   Convert VCF V4.* to FImpute version 3.0 format  -")
print("-   Pipeline written by:  andrew.marete@canada.ca   -")
print("-----------------------------------------------------")
print(
'''
Assumptions:
This script assumes the input:
1. is gzip compressed 
1. is tab delimited and in VCF format
2. is sorted by chromosome then position
''')
# Some Functions: convert the list into 2D list


def chunks(l, n):
    # For item i in a range that is a length of l,
    for i in range(0, len(l), n):
        # Create an index range for l of n items:
        yield l[i:i+n]


# Open and read VCF file
Finvcf = gzip.open(options.genot, "rt")

# Open output files for saving
FoutMark = open(options.outputMark, "w")
FoutGenot = open(options.outputGenot, "w")


# SNP file for FImpute
# has 4 fields: snp_name, bta, pos, snp_index (- starting from 1)
snp_index = 1
FoutMark.write('SNP_ID Chr Pos Index\n')

# empty index fo the genotypes
geno = []

for line in Finvcf:
    # get list of animals
    if "#CHR" in line:
        # anim.append(line.strip().split("\t")[9:])
        anim = line.strip().split("\t")[9:]
    # get DNA as 2D list
    # write SNP_info list
    if "#" not in line:
        geno.append(line.strip().split("\t")[9:])
        bta = line.strip().split("\t")[0]
        pos = line.strip().split("\t")[1]
        snp = bta + ":" + pos
        # print (bta_pos,bta,pos,snp_index)
        FoutMark.write('%s %s %s %s\n' % (snp, bta, pos, snp_index))
        snp_index += 1
FoutMark.close()

# create a recode dictionary for VCF to FImpute notation
recode = {'1/1':'0','0/0':'2','1/0':'1','0/1':'1','./.':'5'}

# transpose the genotype, row == i, column == j
geno_transpose = [[geno[j][i] for j in range(len(geno))] for i in range(len(geno[0]))]

geno_recode = []
for row in geno_transpose:
    hh = ''.join([recode.get(row[x],row[x]) for x in range(0, len(row), 1)])
    geno_recode.append(hh)


# Convert anim list to to 2D list
anim = list(chunks(anim, 1))

# Add the chip number to the animal list
# Add the genotype of the animal
for j in range(0, len(anim)):
    anim[j].append(1)
    anim[j].append(geno_recode[j])

# save the FImpute genotype
FoutGenot.write('ID    Chip    Call...\n')
FoutGenot.writelines(' '.join(str(j) for j in i) + '\n' for i in anim)


FoutGenot.close()
Finvcf.close()

print("Success!, Convertion to Fimpute format completed")
print("\n")
