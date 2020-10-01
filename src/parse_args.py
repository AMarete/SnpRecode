#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse


def msg():
    return '''
    
    Convert to fimpute format: 
    snprecode \\
    -F /path/to/vcf/plink/files \\
    -O PREFIX 
    
    Convert from fimpute format:
    snprecode \\
    -g genotypes_imp.txt \\
    -s snps.txt \\
    -n samples.txt \\ 
    -t 1=vcf; 2=ped|map \\ 
    -a alleles.txt \\
    -o PREFIX
    '''


def my_parser():
    parser = argparse.ArgumentParser(usage=msg())
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    # required.add_argument(' ', nargs='*')
    optional.add_argument("-F", "--FILEPATH", metavar="\b", action="store", dest="filePath",
                          help="dir with vcf and/or ped|map files for converting")
    optional.add_argument("-O", "--OUT", metavar="\b", action="store", dest="OUT",
                          help="PREFIX for fimpute input files")
    # optional.add_argument("-M", "--MARK", metavar="\b", action="store", dest="MARK",
    #                     help="fimpute snp input file")
    optional.add_argument("-g", "--geno", metavar="\b", action="store", dest="geno",
                          help="imputed genotype from fimpute")
    optional.add_argument("-s", "--snps", metavar="\b", action="store", dest="snps",
                          help="snp info file from fimpute")
    optional.add_argument("-n", "--samples", metavar="\b", action="store", dest="samples",
                          help="column file specifying the study samples")
    optional.add_argument("-o", "--out", metavar="\b", action="store", dest="out",
                          help="PREFIX for study genotypes to be saved as vcf or ped|map")
    # optional.add_argument("-m", "--map", metavar="\b", action="store", dest="map",
    #                     help="output map file in plink format")
    optional.add_argument("-t", "--type", metavar="\b", action="store", dest="type_", type=int,
                          help="link code: 1 for vcf, 2 for ped|map")
    optional.add_argument("-a", "--alleles", metavar="\b", action="store", dest="allele",
                          help="column file with alleles formatted as: chrom_pos_ref_alt")

    parser._action_groups.append(optional)
    args = parser.parse_args()
    return args
