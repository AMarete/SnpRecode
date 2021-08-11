from os.path import abspath, isdir
from funtools import file_by_size, to_mat, line_count, find_common, bomb, std_capture, flatten, set_counter
from parse_args import my_parser

# get a list of input files
if not isdir(my_parser().filePath):
    bomb(f'No files found at path = {abspath(my_parser().filePath)}\n')
vcf_list = [file_ for file_ in file_by_size(my_parser().filePath, ['vcf']) if file_.endswith(('vcf', 'vcf.gz'))]
plink_list = to_mat(sorted(file_by_size(my_parser().filePath, ["ped", "map"])), 2)
files = vcf_list + plink_list
file_list = []

# first check
if len(files) > 10:
    bomb('FImpute does not support imputing more than 10 chips simultaneously')
elif len(files) < 1:
    bomb(f'No files found at path = {abspath(my_parser().filePath)}\n')

# sort file from largest smallest
'''
for file in vcf_list + plink_list:
    if str(file).endswith(("vcf", "vcf.gz")):
        unsorted_file_list.append([line_count(file), file])
    else:
        unsorted_file_list.append([line_count(file[0]), file])

for item in sorted(unsorted_file_list, key=lambda x: (int(x[0])), reverse=True):
    file_list.append(item[1])
'''
unsorted_file_list = []
for file in files:
    if str(file).endswith(("vcf", "vcf.gz")):
        command = "zgrep -Ec '$' " + file
        unsorted_file_list.append(line_count(command, file))
    else:
        command = "zgrep -Ec '$' " + file[0]
        unsorted_file_list.append(line_count(command, file))

sorted_l = sorted(unsorted_file_list, key=lambda x: x[-1], reverse=True)
file_list = []
for item in sorted_l:
    file_list.append(' '.join(item[:-1]))

# Error collectors
within_file_snps_dups = {}  # key=file, v = file_snps
within_file_sample_dups = {}  # key=file, v = duplicated sample
between_file_sample_dups = {}  # key=file, v = file_samples
between_file_chroms = {}  # key=file, v = set(file_chroms)

for index, file_ in enumerate(file_list):
    if str(file_list[index]).endswith("vcf.gz"):
        command0 = "zcat< " + file_list[index] + "|grep -m 1  '#CHR' |cut -f10-"
        command1 = "zcat< " + file_list[index] + "| grep -v '#'|cut -f1-2"

        chroms = [i[0] for i in std_capture(command1)]
        snps = [':'.join(i) for i in std_capture(command1)]

        within_file_snps_dups[file_list[index]] = snps
        between_file_chroms[file_list[index]] = set(chroms)

        samples = list(flatten(std_capture(command0)))
        if len(samples) > len(set(samples)):
            within_file_sample_dups[file_list[index]] = set(
                [x for x in samples if samples.count(x) > 1])  # multimode(samples)
        between_file_sample_dups[file_list[index]] = set(samples)

    if str(file_list[index]).endswith("vcf"):
        command0 = "cat " + file_list[index] + "|grep -m 1  '#CHR' |cut -f10-"
        command1 = "cat " + file_list[index] + "| grep -v '#'|cut -f1-2"

        chroms = [i[0] for i in std_capture(command1)]
        snps = [':'.join(i) for i in std_capture(command1)]

        within_file_snps_dups[file_list[index]] = snps
        between_file_chroms[file_list[index]] = set(chroms)

        samples = list(flatten(std_capture(command0)))
        if len(samples) > len(set(samples)):
            within_file_sample_dups[file_list[index]] = set(
                [x for x in samples if samples.count(x) > 1])  # multimode(samples)
        between_file_sample_dups[file_list[index]] = set(samples)

    if str(file_list[index]).endswith(("ped", "map")):
        file_map = file_list[index].split(' ')[0]
        file_ped = file_list[index].split(' ')[1]
        command0 = "awk '{print $2}' " + file_ped
        snps = []
        chroms = []

        samples = list(flatten(std_capture(command0)))
        if len(samples) > len(set(samples)):
            within_file_sample_dups[file_ped] = set(
                [x for x in samples if samples.count(x) > 1])  # multimode(samples)
        between_file_sample_dups[file_ped] = set(samples)

        with open(file_map, 'r') as f:
            for row in f:
                bta, snp, cm, pos = row.strip().split()
                snps.append(bta.replace("Chr", '') + ":" + pos)
                chroms.append(bta.replace("Chr", ''))
            within_file_snps_dups[file_map] = snps
            between_file_chroms[file_map] = set(chroms)

'''
within_file_snps_dups = {}  # key=file, v = file_snps
within_file_sample_dups = {}  # key=file, v = duplicated sample
between_file_sample_dups = {}  # key=file, v = file_samples
between_file_chroms = {}  # key=file, v = set(file_chroms)
'''


def check_dups():
    err = open('Error.txt', 'w')
    # SNPs occurring more than once in same file
    dict0 = {
        k: {x for x in v if v.count(x) > 1}
        for k, v in within_file_snps_dups.items()
        if len(v) > len(set(v))
    }

    if dict0:
        err.write(''.join('SNPs occurring more than once in same file'))
        err.write(''.join('\n------------------------------------------\n'))
        for k, v in dict0.items():
            err.write('{:<3} {:<5} {:<8}\n'.format(k, '>>', ' '.join(v)))
        bomb('> duplicate SNP found. SNP should be unique in each file\n'
             '       > check `less -S Error.txt` for optional additional information\n')

    # sample occurring more than once in same file
    if within_file_sample_dups:
        err.write(''.join('\n--------------------------------------------\n'))
        err.write(''.join('sample occurring more than once in same file'))
        err.write(''.join('\n--------------------------------------------\n'))
        for k, v in within_file_sample_dups.items():
            err.write('{:<3} {:<5} {:<8}\n'.format(k, '>>', ' '.join(v)))
        bomb('> duplicate sample(s) found. Sample(s) should be unique in each file\n'
             '       > check `less -S Error.txt` for optional additional information\n')

    # sample occurring more than once two or more files
    if len(file_list) > 1:
        if between_file_sample_dups:
            err.write(''.join('\n----------------------------------------------------\n'))
            err.write(''.join('sample occurring more than once in two or more files'))
            err.write(''.join('\n----------------------------------------------------\n'))
            err.write(' '.join(find_common(between_file_sample_dups)))

        sets0 = set_counter(between_file_sample_dups.values())
        sets1 = [len(i) for i in [*between_file_sample_dups.values()]]
        if sets0 != sets1:
            bomb('> sample(s) found in more than one file. Sample(s) should be unique to each file\n'
                 '       > check `less -S Error.txt` for optional additional information\n')

        result_chrom = set_counter(between_file_chroms.values())
        if sum(result_chrom) > 0:
            bomb('> files have differing number of chromosomes or chromosome ids\n'
                 '       > check `less -S Error.txt` for optional additional information\n')
