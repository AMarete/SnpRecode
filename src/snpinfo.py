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

mark_chip = list('0' * len(file_list))

# read file
for chip, file in enumerate(file_list):
    for snp_number, line in enumerate(open(file, 'r')):
        try:
            chrom, snpname, cm, pos, A1, A2 = line.strip().split('\t')
        except ValueError:
            chrom, snpname, cm, pos = line.strip().split()

        key = chrom + ":" + pos

        if key not in markers:
            markers[key] = [snpname] + mark_chip
        if key in markers:
            value = list(markers[key])
            value[chip+1] = snp_number + 1
            markers[key] = value

# Create a comprehensive snps file
# Sort the snps by chrom then pos and write index
mark = []
for row in [list(flatten(i)) for i in list(markers.items())]:
    row.insert(0, row[1])
    row.insert(1, row[1].split(":")[0])
    row.insert(2, row[2].split(":")[1])
    row.pop(3)
    row.pop(3)
    mark.append(row)
'''
for row in mark:
    print(row)
'''
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

