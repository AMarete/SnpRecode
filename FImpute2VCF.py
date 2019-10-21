#!/usr/bin/env python
import sys
import time
import datetime
import gzip
from Bio import bgzf
from optparse import OptionParser

start = time.time()
usage = "usage: %script [options] <args>"
parser = OptionParser()

'''define input files'''
parser.add_option("-r", "--RefVCF", action="store", dest="VCF", help="Reference VCF")
parser.add_option("-s", "--snpInfo", action="store", dest="snpInfo", help="FImpute snp file")
parser.add_option("-g", "--ImputedFile", action="store", dest="Impute", help="FImpute geno file")

'''define output files'''
parser.add_option("-o", "--VCFout", action="store", dest="geno", help="Imputed VCF")

'''
(options, args) = parser.parse_args()
'''
try:
    (options, args) = parser.parse_args()
except TypeError:
    print("Make sure input file is specified")
    sys.exit(2)

print(
    '''
         --------------------------------------------------------
         -       Convert FImpute genotype output to VCF         -
         -       Program written by:  andrew.marete@canada.ca   -
         --------------------------------------------------------

This script assumes:
	- The Reference VCF is gzip compressed
	- Genotypes as opposed to Haplotypes were written in the FImpute output

 ''')

# open files
FinRef = gzip.open(options.VCF, "rt")
FinSNPInfo = open(options.snpInfo, "r")
FinImpute = open(options.Impute, "r")

FoutGeno = bgzf.BgzfWriter(options.geno, "wb")
#FoutGeno = open(options.geno, "w")
'''
FinRef = gzip.open("ref.vcf.gz", "rt")
FinImpute = open("genotypes_imp.txt.org", "r")
FinSNPInfo = open("snp_info.txt", "r")

# FoutGeno = bgzf.BgzfWriter("test_out_vcf.vcf.gz", "wb")
FoutGeno = open("test_out_vcf.vcf", "w")
'''

# read the FImpute snp_info file
# save chr:pos as key in a dictionary
# check if its available in Reference VCF and keep alleles
SNP_info = {}
count = 0
for line in FinSNPInfo:
    if "SNPID" in line: continue
    SNPID, Chr, Pos, chip_1, chip_2 = line.strip().split()
    count += 1
    SNP_info[Chr + ":" + Pos] = count

# Read VCF Ref file and save the header
# Read in first nine columns of the Reference VCF file
# save only those SNP retained by FImpute  i.e on the snp_info file
header = []
ref_anims = []
allele_info = {}

for line in FinRef:
    if "##" in line:
        if "bcftools" in line: continue
        header.append(line)
    if "#CHROM" in line:
        ref_anims.append(line.strip().split("\t")[9:])
    if "#" in line: continue
    bta = line.strip().split("\t")[0]
    pos = line.strip().split("\t")[1]
    key = bta + ":" + pos
    if key in SNP_info.keys():
        allele_info[key] = line.strip().split("\t")[0:9]
# header.append("#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT")

# flatten the 2D list
ref_anims = [j for sub in ref_anims for j in sub]

# save the alleles as 2D list
snp_vcf = []
for row in allele_info.values():
    snp_vcf.append(row)

# sort by chrom then pos
snp_vcf = sorted(snp_vcf, key=lambda x: (int(x[0]), int(x[1])))
snp_vcf.insert(0, ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"])

# Read imputed genotype file
# keep animals not in ref list
# recode and transpose on the fly
recode = {'0': '0/0', '2': '1/1', '1': '1/0', '1': '0/1', '5': './.'}

with FinImpute:
    next(FinImpute)
    for line in FinImpute:
        sample, chip, geno = line.strip().split()
        if sample in ref_anims: continue
        snp_vcf[0].append(sample)
        for n, a in enumerate(geno):
            snp_vcf[n+1].append(recode.get(str(a), "!"))

# save the outputs
for row in header:
    FoutGeno.write(''.join(row))
for row in snp_vcf:
    FoutGeno.write('\t'.join(row) + '\n')

FinRef.close()
FinSNPInfo.close()
FinImpute.close()
FoutGeno.close()


end = time.time()
print("Success!, Conversion from FImpute format to VCF")
print("Conversion Time = " + str(datetime.timedelta(seconds=(end - start))))
print("\n")

raise SystemExit(0)
