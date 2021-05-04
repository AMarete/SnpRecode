#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re
import timeit
from recode_dict import vcf_2_fimpute, plink_2_fimpute
from funtools import bomb, open_by_suffix, flatten, file_by_size, to_mat, line_count
from parse_args import my_parser
from garbage import file_list

start = timeit.default_timer()

# from check_dups import vcf_dups, ped_dups
# Open output files
fi = my_parser().OUT + "_geno.txt"
fo = my_parser().OUT + "_snps.txt"
fe = my_parser().OUT + "_alleles"

geno_out = open(fi, "w")
mark_out = open(fo, "w")
allele_out = open(fe, 'w')

# header for FImpute genotype and snps files
geno_out.write('sample_id Chip Call...\n')
mark_out.write('snp bta pos chips...\n')

# Initialize empty dictionary for the various markers
mark_tot = {}
snps_list = {}

for index, file in enumerate(file_list):
    if str(file_list[index]).endswith(("vcf", "vcf.gz")):
        #  Read VCf files one at a time
        samp_ids = []
        geno = []
        chip = index + 1
        snp_number = 0

        with open_by_suffix(file_list[index]) as f:
            for line in f:
                if "##" in line: continue
                # CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,samples...
                bta, pos, snp, ref, alt, qual, filter_, info, format_, genotype = line.strip().split("\t", 9)
                bta = bta.replace('Chr', '')
                # print([ii.split(':')[0] for ii in genotype.split('\t')])
                geno.append(list(flatten([vcf_2_fimpute.get(i, i)
                                          for i in [ii.split(':')[0] for ii in genotype.split('\t')]])))

                if bta != '#CHROM ' and pos != 'POS':
                    if bta + ":" + pos not in snps_list.keys():
                        snps_list[bta + ":" + pos] = [bta, snp, "0", pos, ref, alt]
                    # elif bta + ":" + pos in snps_list.keys() and snps_list[bta + ":" + pos] != [ref + '_' + alt]:
                    elif bta + ":" + pos in snps_list.keys() and snps_list[bta + ":" + pos][4] != ref:
                        print('Some REF/ALT may be flipped \n'
                              'Normalize VCF files e.g. `bcftools norm -f [REF_GENOME] ... `\n')
                        raise SystemExit

                    key = bta + ":" + pos
                    snp_number += 1

                    if key not in mark_tot.keys():
                        mark_chip = list('0' * len(file_list))
                        mark_chip[index] = snp_number
                        mark_tot[key] = mark_chip

                    elif key in mark_tot.keys():
                        value = list(mark_tot[key])
                        value[index] = snp_number
                        mark_tot[key] = value

        for ii in [list(a) for a in zip(*geno)]:
            geno_out.write('{} {} {}\n'.format(ii[0], chip, ''.join(ii[1:])))

    elif str(file_list[index]).endswith(("ped", "map")):
        file_map = file_list[index].split(' ')[0]
        file_ped = file_list[index].split(' ')[1]
        # Parse and write PED files one at a time
        chip = index + 1  # Chip information
        snp_number = 0  # initiate snp number

        with open_by_suffix(file_ped) as f:
            for line in f:
                breed, sample_id, sire, dam, sex, phe, geno = line.strip().split(' ', 6)
                geno = geno.split()
                genotype = [plink_2_fimpute.get(geno[x] + geno[x + 1], '!') for x in range(0, len(geno) - 1, 2)]
                geno_out.write('%s %s %s\n' % (sample_id, str(chip), ''.join(genotype)))

        with open_by_suffix(file_map) as f:
            for line in f:
                snp_number += 1
                bta, snp, cm, pos = line.strip().split()
                key = bta + ":" + pos
                # key = snp
                if key not in mark_tot.keys():
                    value = list('0' * len(file_list))
                    value[index] = snp_number
                    mark_tot[key] = value

                elif key in mark_tot.keys():
                    value = list(mark_tot[key])
                    value[index] = snp_number
                    mark_tot[key] = value
                else:
                    pass
                try:
                    ref = snp.split('_')[2]
                    alt = snp.split('_')[3]
                    # if snp.split('_')[0] + ":" + snp.split('_')[1] not in snps_list.keys():
                    if bta + ":" + pos not in snps_list.keys():
                        # snps_list[snp.split('_')[0] + ":" + snp.split('_')[1]] = [ref + '_' + alt]
                        snps_list[bta + ":" + pos] = [bta, snp, cm, pos, ref, alt]
                    # elif snp.split('_')[0] + ":" + snp.split('_')[1] in snps_list.keys():
                    elif bta + ":" + pos in snps_list.keys() and snps_list[bta + ":" + pos] == [alt + '_' + ref]:
                        print(f'Warning: Allele flipped for SNP {snp} in PLINK file(s)')
                    elif bta + ":" + pos in snps_list.keys() and snps_list[bta + ":" + pos][4] != ref:
                        print(f'Warning: Possible erroneous allele for SNP {snp} in PLINK file(s), normalize with '
                              f'`bcftools norm`')
                except IndexError:
                    print(f'Warning: SNP {snp} in PLINK file(s) not indexed as `chrom_pos_ref_alt`')
                    # raise SystemExit
    else:
        bomb(f'Recheck files at {file_list}')
        raise SystemExit

# Create a comprehensive snps file
mark = []
for row in [list(flatten(i)) for i in list(mark_tot.items())]:
    row.insert(1, row[0].split(":")[0])
    row.insert(2, row[0].split(":")[1])
    mark.append(row)

# Sort the snps by chrom then pos and write index
mark = sorted(mark, key=lambda x: (int(x[1]), int(x[2])))

for row in mark:
    mark_out.write(' '.join(str(e) for e in row) + '\n')

# Write a snp list with ref/alt information
for row in [list(flatten(i)) for i in list(snps_list.values())]:
    #row = row[0].replace(':', '_') + '_' + row[1]
    allele_out.write('\t'.join(row) + '\n')

if len(mark) - len(snps_list) > 0:
    dd = len(mark) - len(snps_list)
    print(f'Warning: {dd} SNPs excluded when writing `{fe}`\n'
          f'These {dd} SNPs will not be available when converting FImpute haplotypes to VCF format\n'
          f'If SNPs already indexed as `chrom_pos_ref_alt`, the alleles may be flipped\n'
          f'Tip: Convert the file to VCF, Normalize, rename SNPs e.g. \n'
          f"     `bcftools norm -f [REF_GENOME] file.vcf.gz -Ou | bcftools annotate --set-sample_id '%CHROM\_%POS\_%REF\_%ALT\'` \n"
          f'Otherwise, converting Fimpute genotypes to PLINK format should be OK\n')
else:
    pass

# close files
geno_out.close()
mark_out.close()
allele_out.close()

stop = timeit.default_timer()
mins, secs = divmod(stop - start, 60)
hours, mins = divmod(mins, 60)

print(f"\nSuccess! files {fi}, {fo}, and {fe} created")
print("Conversion Time (H:M:S): %d:%d:%d\n" % (hours, mins, secs))
print("\n")
