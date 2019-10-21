import os
import sys
import statistics
import time
import datetime
import gzip
import numpy as np
from statistics import mean
from statistics import stdev
from io import StringIO
from optparse import OptionParser

np.seterr(divide='ignore', invalid='ignore')

print('''
--------------------------------------------------------------------
--------------------------------------------------------------------
This program calculates allelic correlation (R^2) for two VCF files
    - For better results filter SNP with < 0.5% MAF 

Author: andrew.marete@canada.ca
--------------------------------------------------------------------
--------------------------------------------------------------------
''')
start = time.time()
usage = "usage: %prog [options] <args>"
parser = OptionParser()

parser.add_option("-t", "--true", action="store", dest="true", help="True genotype matrix")
parser.add_option("-m", "--mask", action="store", dest="mask", help="Mask genotype matrix")

parser.add_option("-o", "--out", action="store", dest="corr", help="Diagonal of correlation matrix")

try:
    (options, args) = parser.parse_args()
except TypeError:
    print("Make sure input files are specified")
    sys.exit(2)

# Some Functions:

# For outputting Error Messages
def bomb(message):
    print("ERROR: " + message)
    sys.exit()


# Function to read various file formats
def open_by_suffix(filename):
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')


# flatten list
def flatten(container):
    for i in container:
        if isinstance(i, (list, tuple)):
            for j in flatten(i):
                yield j
        else:
            yield i


# Correlate two  numpy arrays
def generate_correlation_map(x, y):
    """Correlate each n with each m.

    Parameters
    ----------
    x : np.array
      Shape N X T.

    y : np.array
      Shape M X T.

    Returns
    -------
    np.array
      N X M array in which each element is a correlation coefficient.
      The diagonal gives the correlation between same SNP but from different matrices.

    """
    mu_x = x.mean(1)
    mu_y = y.mean(1)
    n = x.shape[1]
    if n != y.shape[1]:
        raise ValueError('x and y must have the same number of timepoints.')
    s_x = x.std(1, ddof=n - 1)
    s_y = y.std(1, ddof=n - 1)
    cov = np.dot(x,
                 y.T) - n * np.dot(mu_x[:, np.newaxis],
                                   mu_y[np.newaxis, :])
    return cov / np.dot(s_x[:, np.newaxis], s_y[np.newaxis, :])


# Split list or array in 'n' chunks
def chunks(l, n):
    # For item i in a range that is a length of l,
    for i in range(0, len(l), n):
        # Create an index range for l of n items:
        yield l[i:i + n]


recode = {
    '1/1': '1',
    '0/0': '2',
    '1/0': '3',
    '0/1': '3',
    './.': '5',
    '1|1': '1',
    '0|0': '2',
    '1|0': '3',
    '0|1': '3',
    '.|.': '5'
}

# Input and Output Files
'''
t_mat = open('1_5samp.vcf', 'r')
m_mat = open('2_5samp.vcf', 'r')
o_cor = open('toto', "w")
'''
t_mat = open_by_suffix(options.true)
m_mat = open_by_suffix(options.mask)
o_cor = open(options.corr, "w")

geno_dict = {}

# 1st VCF file
for line in t_mat:
    if "#" in line: continue
    anims = len(line.split('\t')) - 9
    snp = line.strip().split("\t")[2]
    geno = line.strip().split("\t")[9:]
    try:
        geno_dict[snp].append(geno)
    except KeyError:
        geno_dict[snp] = [geno]

# second VCF file
for line in m_mat:
    if "#" in line: continue
    snp = line.strip().split("\t")[2]
    if snp not in geno_dict.keys(): continue
    geno = line.strip().split("\t")[9:]
    try:
        geno_dict[snp].append(geno)
    except KeyError:
        geno_dict[snp] = [geno]
        #


dict_geno = {k: geno_dict[k] for k in list(geno_dict)}

# Find the length of the dictionary values
# retain only values of length 2: i.e. contain genotypes for true and masked genotype
# Using these keys subset the key:value pairs in the larger dictionary
length_dict = {key: len(value) for key, value in dict_geno.items()}

snp_key = []
for key, value in length_dict.items():
    if 2 == value:
        snp_key.append(key)
    if 2 != value:
        print(key + ' is not common between the two VCF files and will be ignored')

dict_tot = dict((k, dict_geno[k]) for k in snp_key if k in dict_geno)

# get the genotypes and snps
genotypes = [i for i in list(dict_tot.values())]
snp = [i for i in list(dict_tot.keys())]


mat1 = []
mat2 = []
for item in genotypes:
    for geno1 in item[0]:
        mat1.append(int(recode.get(geno1)))
    for geno1 in item[1]:
        mat2.append(int(recode.get(geno1)))

# Split into 2D arrays with dimension SNP x Animals
mat1 = np.array(list(chunks(mat1, anims)))
mat2 = np.array(list(chunks(mat2, anims)))

# Concordance
# 1. find the dimensions of true - masked
# 2. If concordant, mat1(i,j) - mat2(i,j) = 0
# 3. Count all non zero values i.e. values not concordant
# 4. subtract 3 from 1 to get number of variants concordant
# 5. divide value from 4 by 1 to get ratio concordant

mat3 = mat1 - mat2
var_tot = len(mat3) * len(mat3[0])  # rows * columns
var_concordant = var_tot - np.count_nonzero(mat3)
concordance = round((var_concordant / var_tot) * 100, 2)

print('\n')
print('>>> Concordance is ', str(concordance) + '% ...')

# Split the 2D arrays into smaller arrays of size 500000 SNP x Animals
mat1 = np.array(list(chunks(mat1, 100000)))
mat2 = np.array(list(chunks(mat2, 100000)))
snp = list(chunks(snp, 100000))

# get the number of smaller 2D array chunks
# If the number of chuncks is not equal between the two matrices, through an error
if len(mat1) != len(mat2):
    bomb('The VCF files may not have the same number of SNP')

# get the correlation between mat1 and mat2
# save only the diagonal
# convert it to a dictionary as: snp_name : correlation
# save the dictionary as an output txt file

cor_to_print = []
for array in range(0, len(mat1), 1):
    hh = dict(zip(snp[array], np.diag(generate_correlation_map(mat1[array], mat2[array]))))
    for k, v in hh.items():
        o_cor.write(str(k) + ' ' + str(v) + '\n')
    cor_to_print.append([hh[key] for key in hh])

# close files
t_mat.close()
m_mat.close()
o_cor.close()

# filter nan and 0s
cor_to_print = [x for sublist in cor_to_print for x in sublist]
cor_to_print = [x for x in cor_to_print if str(x) != 'nan' and x != 0]


# function to get mean and standard deviation of R-square values
def Average(lst):
    mn = round(mean(lst), 3)
    std = round(stdev(lst), 3)
    return mn, std


print('>>> R-square is ' + str(Average(cor_to_print)[0]) + ' ± ' + str(Average(cor_to_print)[1]) + ' ...')

raise SystemExit(0)
