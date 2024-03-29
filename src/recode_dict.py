#!/usr/bin/env python
# -*- coding: utf-8 -*-
vcf_2_fimpute = {
    '0/0': '0', '1/1': '2',
    '0/1': '1', '1/0': '1',
    './.': '5', '0/.': '5',
    '1/.': '5', './0': '5',
    './1': '5'
}

plink_2_fimpute = {
    '11': '0',
    '22': '2',
    '21': '1',
    '12': '1',
    '00': '5',
    '10': '5',
    '20': '5',
    '01': '5',
    '02': '5'
}

fimpute_2_vcf = {'0': '0|0', '1': '0/1', '2': '1|1', '3': '0|1', '4': '1|0', '5': '.|.',
                 '6': '0|.', '7': '1|.', '8': '.|0', '9': '.|1'}

recode = {
    '0|0': 0, '1|1': 2, '0|1': 1, '1|0': 1, '.|.': 9, '0|.': 9, '1|.': 9, '.|0': 9, '.|1': 9,
    '0/0': 0, '1/1': 2, '0/1': 1, '1/0': 1, './.': 9, '0/.': 9, '1/.': 9, './0': 9, './1': 9
}