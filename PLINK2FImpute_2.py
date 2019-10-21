#!/usr/bin/env python
import datetime
import sys
import time
import collections.abc
from optparse import OptionParser

start = time.time()

usage = "usage: %script [options] <args>"
parser = OptionParser()

parser.add_option("--g1", action="store", dest="Ref", help="Reference genotype file in Plink ped/map format")
parser.add_option("--g2", action="store", dest="Val", help="Validation genotype file in Plink ped/map format")
parser.add_option("--m1", action="store", dest="MapRef", help="Reference snp file in Plink ped/map format")
parser.add_option("--m2", action="store", dest="MapVal", help="Validation snp file in Plink ped/map format")

parser.add_option("--og", action="store", dest="Genot", help="snp file recoded to FImpute format")
parser.add_option("--om", action="store", dest="Mark", help="genotype file recoded to FImpute format")


try:
    (options, args) = parser.parse_args()
except TypeError:
    print("Make sure input file is specified")
    sys.exit(2)       

print(
    '''
        -----------------------------------------------------------------------------------
        -   Convert Reference and Validation Ped/Map files to FImpute version 3.0 format  -
        -               Pipeline written by:  andrew.marete@canada.ca 	                  -
        -----------------------------------------------------------------------------------

Assumptions:
This script assumes:
1. All SNP in the Validation File are present in the Reference file

If unsure, try "grep -vxf validation_snp reference_snp"; 
If there is an output, some Validation SNP are not included in the Reference panel and should be excluded

''')

FinRef = open(options.Ref, "r")
FinVal = open(options.Val, "r")

FinMapRef = open(options.MapRef, "r")
FinMapVal = open(options.MapVal, "r")

FoutGenot = open(options.Genot, "w")
FoutMark = open(options.Mark, "w")

# Some Functions: convert a list into 2D list


def chunks(l, n):
    # For item i in a range that is a length of l,
    for x in range(0, len(l), n):
        # Create an index range for l of n items:
        yield l[x:x + n]


# Function to Index 2D list


def index_2d(myList, v):
    for i, row in enumerate(myList):
        if v in row:
            # return i, row.index(v)
            return str(row.index(v) + 1)


# Function to flatten an irregular list e.g. ['29:51484561', '29', '51484561', ['979066', '601']] to ['29:51484561', '29', '51484561', '979066', '601']
def flatten(l):
    for el in l:
        if isinstance(el, collections.abc.Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el

'''
# Open input and output files
FinRef = open("test2.ped", "r")
FinVal = open("test1.ped", "r")

FinMapRef = open("ref.map", "r")
FinMapVal = open("val.map", "r")


FoutGenot = open("genotypes", "w")
FoutMark = open("snp_info", "w")
'''

recode = {'11': '0', '22': '2', '12': '1', '21': '1', '00': '5'}

FoutGenot.write('ID    Chip                   Call...\n')
for line in FinRef:
    breed, ID, sire, dam, sex, phe, geno = line.strip().split(' ', 6)
    geno = geno.split()
    genotype = [recode.get(geno[x] + geno[x + 1], '!') for x in range(0, len(geno) - 1, 2)]
    FoutGenot.write('%s %s %s\n' % (ID, 1, ''.join(genotype)))

for line in FinVal:
    breed, ID, sire, dam, sex, phe, geno = line.strip().split(' ', 6)
    geno = geno.split()
    genotype = [recode.get(geno[x] + geno[x + 1], '!') for x in range(0, len(geno) - 1, 2)]
    FoutGenot.write('%s %s %s\n' % (ID, 2, ''.join(genotype)))

FinRef.close()
FinVal.close()

snp_ref = []
for line in FinMapRef:
    bta = line.strip().split("\t")[0]
    pos = line.strip().split("\t")[-1]
    snp_ref.append(bta + ":" + pos + " " + bta + " " + pos)

snp_val = []
for line in FinMapVal:
    bta = line.strip().split("\t")[0]
    pos = line.strip().split("\t")[-1]
    snp_val.append(bta + ":" + pos + " " + bta + " " + pos)

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

'''
print({k: snp_tot_dic[k] for k in list(snp_tot_dic.keys())[:3]})
print({k: snp_val2_dic[k] for k in list(snp_val2_dic.keys())[:3]})
print({k: dic_tot[k] for k in list(dic_tot.keys())[:3]})
'''
# convert the merged dictionary back to a 2d list
dic_tot_2d = [list(i) for i in list(dic_tot.items())]

# for the last column, add 0 to signify snp only present in Ref file and not in Val file
for row in dic_tot_2d:
    if not any(isinstance(el, list) for el in row):  # This line checks if there is an irregular list within a list (e.g. ['29:51484561', '29', '51484561', ['979066', '601']] ) and if present, ignores it
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
FinMapRef.close()
FinMapVal.close()
FoutMark.close()

end = time.time()
print("Success!, Conversion of Reference and Validation files to FImpute format completed")
print("Conversion Time = " + str(datetime.timedelta(seconds=(end - start))))
print("\n")


raise SystemExit(0)
