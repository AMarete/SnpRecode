#!/usr/bin/env python
# -*- coding: utf-8 -*-

from parse_args import my_parser
from os.path import abspath, isdir
from funtools import bomb, file_by_size, to_mat, flatten

try:
    vcf_list = [file_ for file_ in file_by_size(my_parser().filePath, ['vcf']) if file_.endswith(('vcf', 'vcf.gz'))]
    plink_list = to_mat(sorted(file_by_size(my_parser().filePath, ["ped", "map"])), 2)
    file_list = vcf_list + plink_list
except FileNotFoundError:
    bomb(f'No files found at path = {abspath(my_parser().filePath)}\n')


# get a list of input files
def check_path():
    if not isdir(my_parser().filePath):
        return bomb(f'No files found at path = {my_parser().filePath}\n')
    elif isdir(my_parser().filePath) and \
            not [item for item in flatten(file_list) if item.endswith(('vcf', 'vcf.gz', 'ped', 'map'))]:
        return bomb(f'No vcf or plink files found at path = {abspath(my_parser().filePath)}\n')
