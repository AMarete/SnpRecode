#!/usr/bin/env python
import bz2
import datetime
import glob
import gzip
import sys
import time
from optparse import OptionParser

start = time.time()
usage = "usage: %prog [options] <args>"
parser = OptionParser()

parser.add_option("-g", "--geno", action="store", dest="geno", help="FImpute genotype")
parser.add_option("-m", "--mark", action="store", dest="snp", help="FImpute markers")

try:
    (options, args) = parser.parse_args()
except TypeError:
    print("Make sure input file is specified")
    sys.exit(2)

print(
    '''
        -------------------------------------------------------------------------
        -   Convert multiple VCF V4.* files to FImpute version 3.0 format       -
        -               Author:  andrew.marete@canada.ca                        -
        -------------------------------------------------------------------------

Assumptions : This script assumes the VCF files are in the same directory, symbolic links OK

''')
# Some Functions:


def bomb(message):
    print("ERROR: " + message)
    sys.exit()


# Function to read various file formats


def open_by_suffix(filename):
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')

# Function to flatten an irregular list e.g. ['29:51484561', '29', '51484561', ['979066', '601']] to ['29:51484561',
# '29', '51484561', '979066', '601']


def flatten(container):
    for i in container:
        if isinstance(i, (list, tuple)):
            for j in flatten(i):
                yield j
        else:
            yield i


# create a recode dictionary for VCF to FImpute notation
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


# Open output files
FoutGeno = open(options.geno, "w")
FoutMark = open(options.snp, "w")

# header for FImpute genotype
FoutGeno.write('ID Chip Call...\n')

# get a list of the VCF in study directory
vcf_list = sorted(glob.glob('*.vc*'))
if len(vcf_list) > 10:
    bomb('FImpute does not support more than 10 chips in the same genotype file')
elif len(vcf_list) < 2:
    bomb('Not all VCF files are present,check path')

# Initialize empty dictionary for the various markers
mark_tot = {}

# Read VCf files one at a time
for file in range(0, len(vcf_list), 1):
    geno = []  # Initialize empty list for each genotype file
    mark = []
    chip = file + 1  # Chip information
    snp_number = 0

    for line in open_by_suffix(vcf_list[file]):
        if "#CHROM" in line:  # get list of animals in each file
            anim = line.strip().split("\t")[9:]

        if "#" in line: continue  # get SNP_info list per chip

        geno.append(line.strip().split("\t")[9:])
        bta = line.strip().split("\t")[0].replace("Chr", "")
        pos = line.strip().split("\t")[1]

        key = bta + ":" + pos
        snp_number += 1
        mark.append(key)

        if key not in mark_tot.keys():
            mark_chip = list('0' * len(vcf_list))
            mark_chip[file] = snp_number
            mark_tot[key] = mark_chip

        elif key in mark_tot.keys():
            value = list(mark_tot[key])
            value[file] = snp_number
            mark_tot[key] = value

    # recode genotypes to FImpute format
    geno_recode = []
    for row in geno:
        geno_recode.append([recode.get(row[x], row[x]) for x in range(0, len(row), 1)])

    # Initialize an empty dictionary for appending recoded genotypes
    geno_ref = {}

    # transpose using list comprehension
    for snp_row in range(0, len(mark), 1):
        for anim_column in range(0, len(anim), 1):
            key = anim[anim_column]
            value = geno_recode[snp_row][anim_column]
            try:
                geno_ref[key].append(value)
            except KeyError:
                geno_ref[key] = [value]

    # save genotypes accepted by FImpute
    for key, value in geno_ref.items():
        FoutGeno.write('{} {} {}\n'.format(key, chip, ''.join(value)))

FoutGeno.close()

'''
dict_test = open('test_snp_info', 'w')
for k, v in mark_tot.items():
    dict_test.write(str(k) + ' >>> ' + str(v) + '\n')
'''

# convert the merged dictionary back to a 2d list
mark = []

for row in [list(flatten(i)) for i in list(mark_tot.items())]:
    row.insert(1, row[0].split(":")[0])
    row.insert(2, row[0].split(":")[1])
    mark.append(row)

# Sort by chrom and pos
mark = sorted(mark, key=lambda x: (int(x[1]), int(x[2])))

# save the snp_info file
FoutMark.write('snp bta pos chips ...\n')
for row in mark:
    FoutMark.write(' '.join(str(e) for e in row) + '\n')

# close all marker info files
FoutMark.close()


end = time.time()
print("Success!, Conversion of Reference and Validation files to FImpute format completed")
print("Conversion Time = " + str(datetime.timedelta(seconds=(end - start))))
print("\n")

raise SystemExit(0)
