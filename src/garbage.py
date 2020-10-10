from statistics import mode
from funtools import open_by_suffix, file_by_size, to_mat, line_count, counter, find_common, bomb
from parse_args import my_parser
from os.path import isdir

# get a list of input files
if not isdir(my_parser().filePath):
    bomb(f'No files found at path = {my_parser().filePath}\n')
vcf_list = [file_ for file_ in file_by_size(my_parser().filePath, ['vcf']) if file_.endswith(('vcf', 'vcf.gz'))]
plink_list = to_mat(sorted(file_by_size(my_parser().filePath, ["ped", "map"])), 2)
unsorted_file_list = []
file_list = []

for file in vcf_list + plink_list:
    if str(file).endswith(("vcf", "vcf.gz")):
        unsorted_file_list.append([line_count(file), file])
    else:
        unsorted_file_list.append([line_count(file[0]), file])

for item in sorted(unsorted_file_list, key=lambda x: (int(x[0])), reverse=True):
    file_list.append(item[1])

# Error collectors
within_file_snps_dups = {}  # key=file, v = file_snps
within_file_sample_dups = {}  # key=file, v = duplicated sample
between_file_sample_dups = {}  # key=file, v = file_samples
between_file_chroms = {}  # key=file, v = set(file_chroms)

for index, file_ in enumerate(file_list):
    if str(file_list[index]).endswith(("vcf", "vcf.gz")):
        snps = []
        chroms = []
        for line in open_by_suffix(file_list[index]):
            if "#CHROM" in line:
                samples = line.strip().split("\t")[9:]
                if len(samples) > len(set(samples)):
                    within_file_sample_dups[file_list[index]] = set(
                        [x for x in samples if samples.count(x) > 1])  # multimode(samples)
                between_file_sample_dups[file_list[index]] = set(samples)

            elif "#" not in line:
                line = line.strip().split('\t')
                bta = line[0].replace('Chr', '')
                pos = line[1]
                chroms.append(bta)
                snps.append(bta + ":" + pos)
        within_file_snps_dups[file_list[index]] = snps
        between_file_chroms[file_list[index]] = set(chroms)

    if not str(file_list[index]).endswith(("vcf", "vcf.gz")):
        samples = []
        snps = []
        chroms = []
        for line in open(file_list[index][1], 'r'):
            breed, sample_id, sire, dam, sex, phe, geno = line.strip().split(' ', 6)
            samples.append(sample_id)
        if len(samples) > len(set(samples)):
            within_file_sample_dups[file_list[index][1]] = set(
                [x for x in samples if samples.count(x) > 1])  # multimode(samples)
        between_file_sample_dups[file_list[index][1]] = set(samples)

        for row in open(file_list[index][0], 'r'):
            bta, snp, cm, pos = row.strip().split()
            snps.append(bta.replace("Chr", '') + ":" + pos)
            chroms.append(bta.replace("Chr", ''))
        within_file_snps_dups[file_list[index][0]] = snps
        between_file_chroms[file_list[index][0]] = set(chroms)

'''
within_file_snps_dups = {}  # key=file, v = file_snps
within_file_sample_dups = {}  # key=file, v = duplicated sample
between_file_sample_dups = {}  # key=file, v = file_samples
between_file_chroms = {}  # key=file, v = set(file_chroms)
'''


def check_dups():
    err = open('Error.txt', 'w')

    # SNPs occurring more than once in same file
    dict0 = {}
    for k, v in within_file_snps_dups.items():
        if len(v) > len(set(v)):
            dict0[k] = set([x for x in v if v.count(x) > 1])  # multimode(v)  from statistics module

    if dict0:
        err.write(''.join('SNPs occurring more than once in same file'))
        err.write(''.join('\n------------------------------------------\n'))
        for k, v in dict0.items():
            err.write('{:<3} {:<5} {:<8}\n'.format(k, '>>', ' '.join(v)))

    # sample occurring more than once in same file
    if within_file_sample_dups:
        err.write(''.join('\n--------------------------------------------\n'))
        err.write(''.join('sample occurring more than once in same file'))
        err.write(''.join('\n--------------------------------------------\n'))
        for k, v in within_file_sample_dups.items():
            err.write('{:<3} {:<5} {:<8}\n'.format(k, '>>', ' '.join(v)))

    # sample occurring more than once two or more files
    if within_file_sample_dups:
        err.write(''.join('\n----------------------------------------------------\n'))
        err.write(''.join('sample occurring more than once in two or more files'))
        err.write(''.join('\n----------------------------------------------------\n'))
        err.write(' '.join(find_common(between_file_sample_dups)))

    '''
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
    '''
    targets = [*between_file_chroms.values()]
    result = []

    for set_element in targets:
        result.append(len(set_element.difference(set.union(*[x for x in targets if x is not set_element]))))

    if sum(result) > 0:
        bomb('> files have differing number of chromosomes or chromosome ids\n'
             '       > check `less -S Error.txt` for optional additional information\n')
