#!/usr/bin/env python
import sys
import gzip
import time
import datetime
import collections.abc
from optparse import OptionParser

start = time.time()
usage = "usage: %prog [options] <args>"
parser = OptionParser()

parser.add_option("-R", "--ref", action="store", dest="geno_HD", help="Reference VCf file")
parser.add_option("-V", "--val", action="store", dest="geno_LD", help="Validation VCf file")

parser.add_option("-g", "--geno", action="store", dest="outputGenot", help="FImpute genotype")
parser.add_option("-m", "--mark", action="store", dest="outputMark", help="FImpute markers")

try:
    (options, args) = parser.parse_args()
except TypeError:
    print("Make sure input file is specified")
    sys.exit(2)

print(
    '''
        ------------------------------------------------------------------------------
        -   Convert Reference and Validation VCF V4.* to FImpute version 3.0 format  -
        -               Pipeline written by:  andrew.marete@canada.ca 	             -
        ------------------------------------------------------------------------------

Assumptions:
This script assumes the inputs:
1. are gzip compressed
1. are tab delimited and in VCF format
2. are sorted by chromosome then position

If unsure, Normalize the VCF file(s) with BCFtools

''')

# Some Functions: convert the list into 2D list


def chunks(l, n):
    # For item i in a range that is a length of l,
    for i in range(0, len(l), n):
        # Create an index range for l of n items:
        yield l[i:i + n]


# Function to Index 2D list


def index_2d(myList, v):
    for i, row in enumerate(myList):
        if v in row:
            # return i, row.index(v)
            return str(row.index(v) + 1)


# Function to flatten an irregular list e.g. ['29:51484561', '29', '51484561', ['979066', '601']] to ['29:51484561', '29', '51484561', '979066', '601']
'''
def flatten(l):
    for el in l:
        if isinstance(el, collections.abc.Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el
'''
def flatten(container):
    for i in container:
        if isinstance(i, (list, tuple)):
            for j in flatten(i):
                yield j
        else:
            yield i


# create a recode dictionary for VCF to FImpute notation
# recode_unphased = {'0/0': '0', '1/1': '2', '1/0': '1', '0/1': '1', './.': '5'}
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
ref_vcf = gzip.open(options.geno_HD, "rt")
val_vcf = gzip.open(options.geno_LD, "rt")

# Open output files for saving
FoutGenot = open(options.outputGenot, "w")
FoutMark = open(options.outputMark, "w")

'''
# for testing
ref_vcf = gzip.open('test_ref.vcf.gz', 'rt')
val_vcf = gzip.open('test_val.vcf.gz', 'rt')

FoutGenot = open('genotypes', "w")
FoutMark = open('snp_info', "w")
'''
# SNP file for FImpute
# has 4 fields: snp_name, bta, pos, hd_snp_index ld_snp_index (starting from 1)
# snp_index = 1
# FoutMark.write('SNP_ID Chr Pos Chip1 Chip2\n')

# Treat the Reference VCF genotype file
geno_ref = []
snp_ref = []

for line in ref_vcf:
    # get list of animals
    if "#CHR" in line:
        anim_ref = line.strip().split("\t")[9:]
    # get DNA as 2D list
    # get SNP_info list
    if "#" in line: continue
    geno_ref.append(line.strip().split("\t")[9:])
    bta = line.strip().split("\t")[0]
    pos = line.strip().split("\t")[1]
    snp_key = bta + ":" + pos + " " + bta + " " + pos
    snp_ref.append(snp_key)

# recode to FImpute format
geno_ref_recode = []
for row in geno_ref:
    # hh = ' '.join([recode.get(row[x], row[x]) for x in range(0, len(row), 1)])
    # geno_ref_recode.append(hh)
    geno_ref_recode.append([recode.get(row[x], row[x]) for x in range(0, len(row), 1)])

# Initialize an empty dictionary for appending genotypes
geno_ref = {}
# geno.setdefault(key, []).append(value) # append multiple values to same key
# print(geno)

# transpose using list comprehension
for snp_row in range(0, len(snp_ref), 1):
    for anim_column in range(0, len(anim_ref), 1):
        key1 = anim_ref[anim_column]
        value1 = geno_ref_recode[snp_row][anim_column]
        try:
            geno_ref[key1].append(value1)
        except KeyError:
            geno_ref[key1] = [value1]

# Treat the Validation VCF genotype file
geno_val = []
snp_val = []

for line in val_vcf:
    # get list of animals
    if "#CHR" in line:
        anim_val = line.strip().split("\t")[9:]
    # get DNA as 2D list
    # get SNP_info list
    if "#" in line: continue
    geno_val.append(line.strip().split("\t")[9:])
    bta = line.strip().split("\t")[0]
    pos = line.strip().split("\t")[1]
    snp_key = bta + ":" + pos + " " + bta + " " + pos
    snp_val.append(snp_key)

# recode to FImpute format
geno_val_recode = []
for row in geno_val:
    # hh = ' '.join([recode.get(row[x], row[x]) for x in range(0, len(row), 1)])
    # geno_val_recode.append(hh)
    geno_val_recode.append([recode.get(row[x], row[x]) for x in range(0, len(row), 1)])

# Initialize an empty dictionary for appending genotypes
geno_val = {}

# transpose using list comprehension
for snp_row in range(0, len(snp_val), 1):
    for anim_column in range(0, len(anim_val), 1):
        key1 = anim_val[anim_column]
        value1 = geno_val_recode[snp_row][anim_column]
        try:
            geno_val[key1].append(value1)
        except KeyError:
            geno_val[key1] = [value1]

# save the FImpute genotype file
# header
FoutGenot.write('ID Chip Call...\n')
# Reference File
for key, value in geno_ref.items():
    FoutGenot.write('{} {} {}\n'.format(key, 1, ''.join(value)))
# Validation File
for key, value in geno_val.items():
    FoutGenot.write('{} {} {}\n'.format(key, 2, ''.join(value)))
FoutGenot.close()

# 2 - create a unique set of snp from the Ref snp and Val snp list
snp_uniq = list(set(snp_ref + snp_val))

# convert the all snp lists to a 2D lists
snp_tot = []
count = 0
for item in chunks(snp_uniq, 1):
    for ii in item:
        snp_tot.append(ii.split())

snp_ref2 = []
for item in chunks(snp_ref, 1):
    for ii in item:
        snp_ref2.append(ii.split())

snp_val2 = []
for item in chunks(snp_val, 1):
    for ii in item:
        snp_val2.append(ii.split())

# sort the snp lists by chromosome (index = 1) then position (index = 2)
snp_tot = sorted(snp_tot, key=lambda x: (int(x[1]), int(x[2])))
snp_ref2 = sorted(snp_ref2, key=lambda x: (int(x[1]), int(x[2])))
snp_val2 = sorted(snp_val2, key=lambda x: (int(x[1]), int(x[2])))

# add an FImpute required snp index to the larger list == snp_tot
snp_tot_index = []
count = 0
for row in snp_tot:
    count += 1
    row.insert(3, str(count))
    snp_tot_index.append([row[0], row[-1]])

# add an FImpute required snp index to the smaller list == snp_val
snp_val2_index = []
count = 0
for row in snp_val2:
    count += 1
    row.insert(3, str(count))
    snp_val2_index.append([row[0], row[-1]])

# convert the new list to dictionaries
snp_tot_dic = dict()
snp_val2_dic = dict()

for row in snp_tot_index:
    if row[0] in snp_tot_dic:
        # append the new number to the existing array at this slot
        snp_tot_dic[row[0]].append(row[1])
    else:
        # create a new array in this slot
        snp_tot_dic[row[0]] = row[1]

for row in snp_val2_index:
    if row[0] in snp_val2_dic:
        # append the new number to the existing array at this slot
        snp_val2_dic[row[0]].append(row[1])
    else:
        # create a new array in this slot
        snp_val2_dic[row[0]] = row[1]

# merge the two dictionaries with larger dictionary as base
dic_tot = snp_tot_dic
for k, v in snp_val2_dic.items():
    dic_tot[k] = [dic_tot[k], v] if k in dic_tot else v

# print({k: snp_tot_dic[k] for k in list(snp_tot_dic.keys())[:3]})
# print({k: snp_val2_dic[k] for k in list(snp_val2_dic.keys())[:3]})
# print({k: dic_tot[k] for k in list(dic_tot.keys())[:3]})

# convert the merged dictionary back to a 2d list
dic_tot_2d = [list(i) for i in list(dic_tot.items())]

# for the last column, add 0 to signify snp only present in Ref file and not in Val file
for row in dic_tot_2d:
    if not any(isinstance(el, list) for el in
               row):  # This line checks if there is an irregular list within a list (e.g. ['29:51484561', '29', '51484561', ['979066', '601']] ) and if present, ignores it
        row.append("0")

# save the snp_info file
FoutMark.write('snp bta pos chip1 chip2\n')
for row in dic_tot_2d:
    row = list(flatten(row))  # this line calls a function that flattens 'row' if it is irregular
    bta = row[0].split(":")[0]
    pos = row[0].split(":")[1]
    row.insert(1, bta)
    row.insert(2, pos)
    FoutMark.write(" ".join(row) + '\n')


# close all marker info files
FoutMark.close()
FoutGenot.close()

ref_vcf.close()
val_vcf.close()

end = time.time()
print("Success!, Conversion of Reference and Validation files to FImpute format completed")
print("Conversion Time = " + str(datetime.timedelta(seconds=(end - start))))
print("\n")

raise SystemExit(0)
