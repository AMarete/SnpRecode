#!/usr/bin/env python
# -*- coding: utf-8 -*-


from funtools import flatten
from parse_args import my_parser
import timeit

start = timeit.default_timer()
'''
print("convert one or several bim and/or map files to fimpute format\n"
      "coded by: andrew marete (C) 2015\n")
'''

# container for all snps
markers = {}

# get the snp files in order, largest to smallest
in_files = []
for item in my_parser().file:
    if item.name.endswith(('bim', 'map')):
        in_files.append(item.name)

file_list = {}
for file in in_files:
    file_list[file] = len(open(file).readlines())
file_list = sorted(file_list, key=file_list.get, reverse=True)

# read file
for index, file in enumerate(file_list):
    snp_number = 1
    chip = index + 1
    for line in open(file, 'r'):
        try:
            chrom, snp, cm, pos, A1, A2 = line.strip().split('\t')
        except ValueError:
            chrom, snp, cm, pos = line.strip().split()
        if chrom + ":" + pos not in markers.keys():
            mark_chip = list('0' * len(file_list))
            mark_chip[index] = snp_number
            markers[snp + ":" + chrom + ":" + pos] = mark_chip
            snp_number += 1
        elif chrom + ":" + pos in markers.keys():
            value = list(markers[snp + ":" + chrom + ":" + pos])
            value[index] = snp_number
            markers[snp + ":" + chrom + ":" + pos] = value
            snp_number += 1

# Create a comprehensive snps file
# Sort the snps by chrom then pos and write index
mark = []
for row in [list(flatten(i)) for i in list(markers.items())]:
    row.insert(1, row[0].split(":")[0])
    row.insert(2, row[0].split(":")[1])
    row.insert(3, row[0].split(":")[2])
    mark.append(row[1:])

mark = sorted(mark, key=lambda x: (int(x[1]), int(x[2])))

f_out = "snpinfo.txt"
mark_out = open(f_out, "w")
mark_out.write('snp bta pos chips...\n')
for row in mark:
    mark_out.write(' '.join(str(e) for e in row) + '\n')

stop = timeit.default_timer()
mins, secs = divmod(stop - start, 60)
hours, mins = divmod(mins, 60)

print(f"\nSuccess! created {f_out}")
print("Runtime (H:M:S): %d:%d:%d\n" % (hours, mins, secs))
print("\n")
raise SystemExit
