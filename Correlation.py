import os 
import statistics 
import time
import datetime
import gzip
import numpy as np
from io import StringIO
from optparse import OptionParser
np.seterr(divide='ignore', invalid='ignore')

print('''
--------------------------------------------------------------------
--------------------------------------------------------------------
This program calculates allelic correlation (R^2) for two VCF files
Assumptions for the two files:
        1. The files are (b)gzip compressed
        2. The column (sample) names are identical
        3. The row (SNP) names are identical 

Author: andrew.marete@canada.ca
--------------------------------------------------------------------
--------------------------------------------------------------------
''')


usage = "usage: %prog [options] <args>"
parser = OptionParser()

parser.add_option("-T", "--true", action="store", dest="geno_true", help="True genotype matrix")
parser.add_option("-M", "--mask", action="store", dest="geno_mask", help="Mask genotype matrix")

parser.add_option("-O", "--out", action="store", dest="corr_out", help="Diagonal of correlation matrix")

try:
    (options, args) = parser.parse_args()
except TypeError:
    print("Make sure input files are specified")
    sys.exit(2)

# Files
t_mat = gzip.open(options.geno_true, "rt")
m_mat = gzip.open(options.geno_mask, "rt")
o_cor = open(options.corr_out, "w")


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


def chunks(l, n):
    # For item i in a range that is a length of l,
    for i in range(0, len(l), n):
        # Create an index range for l of n items:
        yield l[i:i+n]



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

list1 = []
snp = []
for line in t_mat:
	if "#" in line: continue
	anims1 = len(line.split('\t')) - 9
	snp.append(line.strip().split("\t")[2])
	geno = line.strip().split("\t")[9:]
	for genotype in geno:
		list1.append(int(recode.get(genotype)))

list2 = []
for line in m_mat:
	if "#" in line: continue	
	anims2 = len(line.split('\t')) - 9
	#snp2 = line.strip().split("\t")[2]
	geno = line.strip().split("\t")[9:]
	for genotype in geno:
		list2.append(int(recode.get(genotype)))


mat1 = np.array(list(chunks(list1, anims1)))
mat2 = np.array(list(chunks(list2, anims2)))

# get the correlation between mat1 and mat2
# save only the diagonal
# convert it to a dictionary as: snp_name : correlation
# save the dictionary as an output txt file
hh = dict(zip(snp, np.diag(generate_correlation_map(mat1, mat2))))

for k, v in hh.items():
    o_cor.write(str(k) + ' ' + str(v) + '\n')

# close files
t_mat.close()
m_mat.close()
o_cor.close()

# correlation for printing in log file
numbers = [hh[key] for key in hh]
mean_ = round(statistics.mean(numbers), 3)
stdv_ = round(statistics.stdev(numbers, mean_), 3)

# Concordance
# 1. find the dimensions of true - masked
# 2. If concordant, mat1(i,j) - mat2(i,j) = 0
# 3. Count all non zero values
# 4. subtract 3 from 1 to get number of variants concordant
# 5. divide 4 by 1 to get ratio concordant

mat3 = mat1 - mat2
var_tot = len(mat3) * len(mat3[0])  # rows * columns
var_concordant = var_tot - np.count_nonzero(mat3)
concordance = round((var_concordant / var_tot) * 100, 2)


print('...Concordance is ', str(concordance) + '% ...')
print('...R-square is ' + str(mean_) + ' ± ' + str(stdv_) + '...')
