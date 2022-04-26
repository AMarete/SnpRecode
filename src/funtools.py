#!/usr/bin/env python

import gzip
import os
import subprocess


class MyCols:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''


# Welcome Message
def welcome_message():
    print("""
----------------------------------------
FImpute helper utility (Version 1.0.4)
Copyright (C) 2018-2021 Andrew Marete
----------------------------------------""")


# function to display Error messages
def bomb(message):
    print(f"{MyCols.FAIL}Error: {message}{MyCols.ENDC}")
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


# Function to flatten an irregular list e.g. ['29:51484561', '29', '51484561', ['979066', '601']] to ['29:51484561',
# '29', '51484561', '979066', '601']
def flatten(container):
    for i in container:
        if isinstance(i, (list, tuple)):
            yield from flatten(i)
        else:
            yield i


# count lines in file
'''
def line_count(filename):
    with open_by_suffix(filename) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
'''


def line_count(command0, filename):
    data = subprocess.Popen([command0], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,
                            universal_newlines=True)
    data = [i.strip() for i in data.stdout]
    return list(flatten([filename, data]))


# Function to creat 2D list from list
def to_mat(items, n):
    return [items[i:i + n] for i in range(0, len(items), n)]


# count dups
def counter(data):
    d = {}
    for item in data:
        d[item] = d.get(item, 0) + 1
    return d


# set counter
def set_counter(sets):
    """
    Symmetrical Difference Between List of Sets of Strings
    Find uncommon elements between sets

    set1 = {'A', 'B', 'C'}
    set2 = {'C', 'D', 'E'}
    set3 = {'D', 'E', 'F'}

    targets = [set1, set2, set3]
    result = [2, 0, 1]
    If all elements are same within set1, set2, set3, i.e. no unique element in any given set,
     then sum(result) = 0
    Explanation:
        in set1, A and B are not found in any of the other sets,
        in set2, there are no unique elements to the set,
        in set3, F is not found in any of the other sets
    """
    targets = [*sets]
    return [
        len(
            set_element.difference(
                set.union(*[x for x in targets if x is not set_element])
            )
        )
        for set_element in targets
    ]


# find duplicated values in some dictionary.values()
def find_common(data):
    """
    Find duplicated values in a dictionary value is in all dictionary values
    e.g. 2 is common in all values D = {'A': {1, 2, 3}, 'B': {2, 4, 5}, 'C': {1, 2, 7}}
    use : set.intersection(*D.values())

    if value is common between some of the dict.values
    e.g. 2 is only common twice  D = {'A': {1, 9, 3}, 'B': {2, 4, 5}, 'C': {1, 2, 7}}, use this fxn
    """
    flipped = {}
    for i in data:
        for value in data[i]:
            if value not in flipped:
                flipped[value] = [i]
            else:
                flipped[value].append(i)
    return [i for i, value_ in flipped.items() if len(value_) > 1]


# Correlations
def mean_(list0):
    total = sum(float(a) for a in list0)
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


# Capture subprocess
def std_capture(command):
    try:
        data = subprocess.Popen([command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,
                                universal_newlines=True)
        dat = [line.strip().split('\t') for line in data.stdout]
        return dat

    except FileNotFoundError:
        data = subprocess.Popen([command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,
                                universal_newlines=True)
        dat = [line.strip().split('\t') for line in data.stdout]
        return dat
