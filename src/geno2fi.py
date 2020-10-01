#!/usr/bin/env python
# -*- coding: utf-8 -*-

import timeit
from recode_dict import vcf_2_fimpute, plink_2_fimpute
from funtools import bomb, open_by_suffix, flatten, file_by_size, to_mat, line_count
from parse_args import my_parser

start = timeit.default_timer()

# from check_dups import vcf_dups, ped_dups
# Open output files
fi = my_parser().OUT + ".geno"
fo = my_parser().OUT + ".mark"
fe = my_parser().OUT + ".allele"

geno_out = open(fi, "w")
mark_out = open(fo, "w")
allele_out = open(fe, 'w')

# header for FImpute genotype and snps files
geno_out.write('sample_id Chip Call...\n')
mark_out.write('snp bta pos chips...\n')

# get a list of input files
vcf_list = [file_ for file_ in file_by_size(my_parser().filePath, ['vcf']) if file_.endswith(('vcf', 'vcf.gz'))]
plink_list = to_mat(sorted(file_by_size(my_parser().filePath, ["ped", "map"])), 2)
unsorted_file_list = []
file_list = []

for file in vcf_list + plink_list:
    if "vcf" in file:
        unsorted_file_list.append([line_count(file), file])
    else:
        unsorted_file_list.append([line_count(file[0]), file])

for item in sorted(unsorted_file_list, key=lambda x: (int(x[0])), reverse=True):
    file_list.append(item[1])

if len(file_list) > 10:
    bomb('FImpute does not support more than 10 chips in the same genotype index')
elif len(file_list) < 1:
    bomb('Recheck path to genotype files')

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

        for line in open_by_suffix(file_list[index]):
            if "#CHROM" in line:  # get list of animals in each index
                samp_ids.append(line.strip().split("\t")[9:])
            if "#" in line:
                continue  # get SNP_info list per chip
            # SNP information
            bta = line.strip().split("\t")[0].replace("Chr", "")
            pos = line.strip().split("\t")[1]
            ref = line.strip().split("\t")[3]
            alt = line.strip().split("\t")[4]

            if bta + ":" + pos not in snps_list.keys():
                snps_list[bta + ":" + pos] = [ref + '_' + alt]
            # elif bta + ":" + pos in snps_list.keys() and snps_list[bta + ":" + pos] != [ref + '_' + alt]:
            elif bta + ":" + pos in snps_list.keys() and snps_list[bta + ":" + pos][0].split('_')[0] != ref:
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

            # genotype information
            genotype = line.strip().split("\t")[9:]
            geno.append([vcf_2_fimpute.get(genotype[x].split(":")[0], genotype[x].split(":")[0]) for x in
                         range(0, len(genotype), 1)])

        for x, row in enumerate(zip(*geno)):
            geno_out.write('{} {} {}\n'.format(samp_ids[0][x], chip, ''.join(row)))

    else:
        # Parse and write PED files one at a time
        chip = index + 1  # Chip information
        snp_number = 0  # initiate snp number

        for line in open(file_list[index][1], 'r'):
            breed, sample_id, sire, dam, sex, phe, geno = line.strip().split(' ', 6)
            geno = geno.split()
            genotype = [plink_2_fimpute.get(geno[x] + geno[x + 1], '!') for x in range(0, len(geno) - 1, 2)]
            geno_out.write('%s %s %s\n' % (sample_id, str(chip), ''.join(genotype)))

        for row in open(file_list[index][0], 'r'):
            snp_number += 1
            bta, snp, cm, pos = row.strip().split()
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
                    snps_list[bta + ":" + pos] = [ref + '_' + alt]
                # elif snp.split('_')[0] + ":" + snp.split('_')[1] in snps_list.keys():
                elif bta + ":" + pos in snps_list.keys() and snps_list[bta + ":" + pos] == [alt + '_' + ref]:
                    print(f'Warning: Allele flipped for SNP {snp} in PLINK file(s)')
                elif bta + ":" + pos in snps_list.keys() and snps_list[bta + ":" + pos][0].split('_')[0] != ref:
                    print(f'Warning: Possible erroneous allele for SNP {snp} in PLINK file(s), normalize with '
                          f'`bcftools norm`')
            except IndexError:
                print(f'Warning: SNP {snp} in PLINK file(s) not indexed as `chrom_pos_ref_alt`')
                # raise SystemExit

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
for row in [list(flatten(i)) for i in list(snps_list.items())]:
    row = row[0].replace(':', '_') + '_' + row[1]
    allele_out.write(''.join(row) + '\n')

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
