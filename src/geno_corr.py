import timeit
from collections import defaultdict
from statistics import mean
import pylab as pl
from funtools import allelic_r2, open_by_suffix, bomb
from parse_args import my_parser
from recode_dict import recode

start = timeit.default_timer()

# inputs

in_files = [
    item.name
    for item in my_parser().file
    if item.name.endswith(('vcf', 'vcf.gz'))
]

if len(in_files) != 2:
    bomb('''two vcf files required to calculate genotype correlation\n\ttry: `snprecode --file vcf1 vcf2`''')

'''
file1 = my_parser().file[0].name
file2 = my_parser().file[1].name
'''
file1 = in_files[0]
file2 = in_files[1]
mat1 = open_by_suffix(file1)
mat2 = open_by_suffix(file2)

# outputs
file_txt = 'genotype_R2.txt'
file_plot = "genotype_R2.pdf"

with open(file_txt, 'w') as outfile:
    outfile.write('marker\tmaf\tR2\n')

    result = defaultdict(list)
    samples = []
    for line in mat1:
        if "CHR" in line:
            samples.append(line.strip().split('\t')[9:])
        if "#" in line:
            continue
        snp = line.strip().split('\t')[0] + ":" + line.strip().split('\t')[1]
        geno = line.strip().split('\t')[9:]
        if any("2" in s for s in geno):
            bomb(f'multi-allelic sites not supported\n\ttry bcftools norm -m- {file1}')
        af = round(((''.join(geno).count('1')) / (len(''.join(geno).replace('|', '').replace('/', '')))), ndigits=4)
        result[snp] = [[af]] if af <= 0.5 else [[round(1 - af, ndigits=4)]]
        dosage = [int(recode.get(i.split(":")[0])) for i in geno]
        try:
            result[snp].append(dosage)
        except KeyError:
            continue

    for line in mat2:
        if "CHR" in line:
            samples.append(line.strip().split('\t')[9:])
        if "#" in line:
            continue
        snp = line.strip().split('\t')[0] + ":" + line.strip().split('\t')[1]
        geno = line.strip().split('\t')[9:]
        if any("2" in s for s in geno):
            bomb(f'multi-allelic sites not supported\n\ttry bcftools norm -m- {file2}')
        dosage = [int(recode.get(i.split(":")[0])) for i in geno]
        try:
            result[snp].append(dosage)
        except KeyError:
            continue

    if samples[0] != samples[1]:
        print(f'''Error: different or unordered samples in {file1} and {file2}\n       
          Samples in the two vcf has to be same and ordered''')
        raise SystemExit

    result = {k: v for k, v in result.items() if len(v) == 3}

    r2 = {
        k: [float(str(v[0])[1:-1]), allelic_r2(v[1], v[2])]
        for k, v in result.items()
    }

    # save file
    dt_plot = defaultdict(list)
    for k, v in r2.items():
        outfile.write('{}\t{}\t{}\n'.format(k, str(v[0]), str(v[1])))
        if v[0] <= 0.01:
            try:
                dt_plot["0.005 - 0.01"].append(v[1])
            except KeyError:
                dt_plot["0.005 - 0.01"] = v[1]
        elif 0.01 < v[0] <= 0.05:
            try:
                dt_plot["0.01 - 0.05"].append(v[1])
            except KeyError:
                dt_plot["0.01 - 0.05"] = v[1]
        elif 0.05 < v[0] <= 0.10:
            try:
                dt_plot["0.05 - 0.10"].append(v[1])
            except KeyError:
                dt_plot["0.05 - 0.10"] = v[1]
        elif 0.1 < v[0] <= 0.2:
            try:
                dt_plot["0.1 - 0.2"].append(v[1])
            except KeyError:
                dt_plot["0.1 - 0.2"] = v[1]
        elif 0.2 < v[0] <= 0.3:
            try:
                dt_plot["0.2 - 0.3"].append(v[1])
            except KeyError:
                dt_plot["0.2 - 0.3"] = v[1]
        elif 0.3 < v[0] <= 0.4:
            try:
                dt_plot["0.3 - 0.4"].append(v[1])
            except KeyError:
                dt_plot["0.3 - 0.4"] = v[1]
        else:
            try:
                dt_plot["0.4 - 0.5"].append(v[1])
            except KeyError:
                dt_plot["0.4 - 0.5"] = v[1]

    # get the coordinates for plotting
    xy = {}
    for k, v in dt_plot.items():
        if 0 in v:
            v.remove(0)
        xy[k] = round(mean(v), ndigits=3)
    xy_order = ["0.005 - 0.01", "0.01 - 0.05", "0.05 - 0.10",
                "0.1 - 0.2", "0.2 - 0.3", "0.3 - 0.4", "0.4 - 0.5"]

    xy_ordered = {k: xy[k] for k in xy_order}

    # Mean correlation
    mu = list(xy_ordered.values())
    mu = round(mean(mu[1:]), ndigits=2)

    # plot
    f = pl.figure()
    x = [0, 1, 2, 3, 4, 5, 6]
    xTicks = list(xy_ordered.keys())
    y = list(xy_ordered.values())
    pl.xticks(x, xTicks)
    pl.xticks(range(7), xTicks, rotation=45)
    pl.plot(x, y, '-*b')
    pl.xlabel('MAF', fontsize=12, fontweight='bold')
    pl.ylabel('$R^2$', fontsize=12, fontweight='bold')
    pl.title(r'$summary\ genotype\ correlation,\ R^2={}$'.format(mu), fontsize=12, fontweight='bold')
    # pl.show()

    f.savefig(file_plot, bbox_inches='tight')

    mat1.close()
    mat2.close()
stop = timeit.default_timer()
mins, secs = divmod(stop - start, 60)
hours, mins = divmod(mins, 60)

print(f"\nSuccess! R-sq = {mu}\nfiles created: \n\t1. {file_txt} \n\t2. {file_plot}")
print("Runtime (H:M:S): %d:%d:%d\n" % (hours, mins, secs))
print("\n")
raise SystemExit
