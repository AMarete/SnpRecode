#!/usr/bin/env python

import gzip
import os


# Welcome Message
def welcome_message():
    print(
        """----------------------------------------
FImpute helper utility (Version 1.0.1)
Copyright (C) 2018-2020 Andrew Marete
----------------------------------------""")


# function to display Error messages
def bomb(message):
    print(f"Error: {message}")
    raise SystemExit


# Function to read various index formats
def open_by_suffix(filename):
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')


# line count
'''
def line_count(filename):
    return int(subprocess.check_output(['wc', '-l', filename]).split()[0])
'''


# Function to create a list object of files sorted by size, largest to smallest
def file_by_size(dirname, filetype, reverse=True):
    filepaths = []
    for filename in os.listdir(dirname):
        if any(suffix in filename for suffix in filetype):
            filename = os.path.join(dirname, filename)
            if os.path.isfile(filename):
                filepaths.append(filename)

    # Re-populate list with filename, size tuples
    for i, x in enumerate(filepaths):
        filepaths[i] = (filepaths[i], os.path.getsize(filepaths[i]))

    # Sort list by index size
    # If reverse=True sort from largest to smallest

    # If reverse=False sort from smallest to largest
    def function_(filename_):
        return filename_[1]

    filepaths.sort(key=function_, reverse=reverse)

    # Re-populate list with just filenames
    for i, j in enumerate(filepaths):
        filepaths[i] = filepaths[i][0]

    return filepaths


def line_count(filename):
    with open_by_suffix(filename) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


# Function to flatten an irregular list e.g. ['29:51484561', '29', '51484561', ['979066', '601']] to ['29:51484561',
# '29', '51484561', '979066', '601']
def flatten(container):
    for i in container:
        if isinstance(i, (list, tuple)):
            for j in flatten(i):
                yield j
        else:
            yield i


# Function to creat 2D list from list
def to_mat(items, n):
    return [items[i:i + n] for i in range(0, len(items), n)]


# count dups
def counter(data):
    d = {}
    for item in data:
        d[item] = d.get(item, 0) + 1
    return d


# Find duplicated values in a dictionary
# value is in all dictionary values
# e.g. 2 is common in all values D = {'A': {1, 2, 3}, 'B': {2, 4, 5}, 'C': {1, 2, 7}}
# use : set.intersection(*D.values())
# if value is common between some of the dict.values
# e.g. 2 is only common twice  D = {'A': {1, 9, 3}, 'B': {2, 4, 5}, 'C': {1, 2, 7}}, use fxn below
def find_common(data):
    flipped = {}
    out = []
    for i in data:
        for value in data[i]:
            if value not in flipped:
                flipped[value] = [i]
            else:
                flipped[value].append(i)
    for i in flipped:
        if len(flipped[i]) > 1:
            out.append(i)
    return out


# Correlations
def mean_(list0):
    total = 0
    for a in list0:
        total += float(a)
    return total / len(list0)


def std(list0):
    listMean = mean_(list0)
    dev = 0.0
    for i in range(len(list0)):
        dev += (list0[i] - listMean) ** 2
    return dev ** (1 / 2.0)


def allelic_r2(list0, list1):
    # First establish the means and standard deviations for both lists.
    xMean = mean_(list0)
    yMean = mean_(list1)
    xStandDev = std(list0)
    yStandDev = std(list1)
    # r numerator
    rNum = 0.0
    for i in range(len(list0)):
        rNum += (list0[i] - xMean) * (list1[i] - yMean)

    # r denominator
    rDen = xStandDev * yStandDev
    try:
        return round((rNum / rDen) ** 2, ndigits=4)
    except ZeroDivisionError:
        return 0
