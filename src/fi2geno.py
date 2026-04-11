import timeit
import re
from datetime import datetime
from contextlib import suppress
from Bio import bgzf
from funtools import flatten, bomb
from parse_args import my_parser
from recode_dict import fimpute_2_vcf

start = timeit.default_timer()

def open_file(file_path):
    with suppress(FileNotFoundError):
        return open(file_path, "r")
    bomb(f"Missing argument or '{file_path}' may be empty\n")

samples = open_file(my_parser().samples)
snp_info = open_file(my_parser().snps)
geno_info = open_file(my_parser().geno)
allele_info = open_file(my_parser().allele)

toto = my_parser().type_

if not my_parser().out:
    bomb('Missing argument, "-o PREFIX", "--out PREFIX"\n'
         'Error: run `./snprecode -h` for complete arguments list '
         'required to recode from FImpute\n')

fo = my_parser().out
# open headers
if toto == 1:
    geno_out = bgzf.BgzfWriter(f"{fo}.vcf.gz", "wb")
    # write header
    geno_out.write(f'''##fileformat=VCFv4.2
##filedate={datetime.now().strftime('%Y%m%d')}
##source="snprecode v1.0.4"
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
''')
elif toto == 2:
    geno_out = open(f"{fo}.ped", "w")
else:
    print(f'\nError!: Missing argument. Specify the recode type: `-t 1` (for VCF) or `-t 2` (for PED/MAP)]')
    print('run `./snprecode -h` for complete arguments list required to recode from FImpute\n')
    raise SystemExit

# get samples files
study_samples = []
with samples as file:
    study_samples.extend(line.strip() for line in file)

# Fimpute marker info
snps = {}
plink_info = {}
count = 0
with snp_info as file:
    next(file)
    for line in file:
        # marker, chrom, pos, chips = re.sub('\s+', '\t', line).split('\t', 3)
        marker, chrom, pos, chips = line.split(maxsplit=3)
        snps[marker.lower()] = [chrom, pos, marker]
        plink_info[count] = [f"{chrom}:{pos}"]
        count += 1

# Alleles file
vcf_a1a2 = {}
plink_ACGT = {}
count = 0
with allele_info as file:
    for line in file:
        chrom, marker, cm, pos, a1, a2 = line.strip().split()
        vcf_a1a2[marker.lower()] = [chrom, pos, marker, a1, a2, ".", "PASS", ".", "GT"]
        # vcf_a1a2[f"{chrom}:{pos}"] = [chrom, pos, marker, a1, a2, ".", "PASS", ".", "GT"]
        if a1 == '.': a1 = '0'
        if a2 == '.': a2 = '0'
        if a1 == '0':  # monomorphic marker (A1=0)
            plink_ACGT[f"{chrom}:{pos}"] = [f"{a1} {a1}", f"{a1} {a1}", f"{a2} {a2}", f"{a1} {a1}",
                                             f"{a1} {a1}", f"{a1} {a1}", f"{a1} {a1}", f"{a1} {a1}",
                                             f"{a1} {a1}", f"{a1} {a1}", marker]
        else:
            plink_ACGT[f"{chrom}:{pos}"] = [f"{a1} {a1}", f"{a1} {a2}", f"{a2} {a2}", f"{a1} {a2}",
                                             f"{a2} {a1}", '0 0', '0 0', '0 0', '0 0', '0 0', marker]
        count += 1

# For VCF only
# Add header for first 9 lines and write vcf
vcfsnps = {k: vcf_a1a2[k] for k in vcf_a1a2.keys() & set(snps.keys())}.values()
vcfsnps = sorted(vcfsnps, key=lambda x: (int(x[0]), int(x[1])))
vcfsnps.insert(0, ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"])
genotypes = []

# imputed geno file
with geno_info as gfile:
    next(gfile)
    for line in gfile:
        geno_tot = []
        sample_id, chip, geno = line.strip().split()
        if sample_id not in study_samples: continue
        if toto == 1:
            row = [fimpute_2_vcf.get(str(a), "!") for a in geno]
            row.insert(0, sample_id)
            genotypes.append(row)
        elif toto == 2:
            for n, a in enumerate(geno):
                with suppress(KeyError):
                    geno_tot.append(plink_ACGT[plink_info[n][0]][int(a)])
            geno_out.write('%s %s \n' % (' '.join([sample_id, sample_id, '0', '0', '2', '-9']), ' '.join(geno_tot)))

# continue writing VCF
if toto == 1:
    for row in zip(vcfsnps, [*zip(*genotypes)]):
        line = [list(flatten(i)) for i in list(row)]
        geno_out.write('\t'.join(flatten(line)) + '\n')

# Map file for Ped
if toto == 2:
    map_out = open(f"{fo}.map", "w")
    map_ped = []
    for marker in snps:  # snps is a list need to ammend this
        try:
            ms = vcf_a1a2[marker][:3]
            ms.insert(0, '0')
            ms = [ms[i] for i in [1, 3, 0, 2]]
            map_ped.append(ms)
        except KeyError:
            ms = snps[marker]
            ms.insert(0, '0')
            ms = [ms[i] for i in [1, 3, 0, 2]]
            map_ped.append(ms)
    map_ped = sorted(map_ped, key=lambda x: (int(x[0]), int(x[3])))
    with map_out as file:
        for row in map_ped:
            file.write('\t'.join(row) + '\n')

geno_out.close()

tt = f'{fo}.vcf.gz' if toto == 1 else ' and '.join([f'{fo}.ped', f'{fo}.map'])
stop = timeit.default_timer()
mins, secs = divmod(stop - start, 60)
hours, mins = divmod(mins, 60)

print(f"Successful conversion from FImpute format to {tt}")
print("Conversion Time (H:M:S): %d:%d:%d\n" % (hours, mins, secs))
raise SystemExit(0)