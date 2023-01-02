import os

### snakemake parameters
regions_file = snakemake.input.regions_file
bin_size = snakemake.params.bin_size
bin_regions_file = snakemake.output.bin_regions_file
chrom_regions_files = snakemake.output.chrom_regions_files

### sort the input regions file to make sure chrom bin counter is accurate
regions_file_sorted = '/'.join(bin_regions_file.split('/')[:-1]) + "/regions.sorted.bed.tmp"
os.system("sort -k1V -k2n -k3n {} > {}".format(regions_file, regions_file_sorted))
###
current_chrom = ''
region_counter = 0
with open(regions_file_sorted, 'r') as in_f:
    in_regions = in_f.readlines()
for chrom_counter, chrom_regions_file in enumerate(chrom_regions_files):
    bin_counter = 1
    with open(chrom_regions_file, 'w') as chrom_regions_f:
        for line in in_regions[region_counter:]:
            columns = line.split("\t")
            chrom = columns[0]
            if (chrom != current_chrom) and (current_chrom != ''): # not first chrom
                current_chrom = chrom
                break
            elif chrom != current_chrom: # first chrom
                current_chrom = chrom
            region_counter += 1
            region_start = int(columns[1])
            region_end = int(columns[2])
            bin_start = region_start
            bin_end = bin_start + bin_size
            while bin_start < region_end:
                if bin_end > region_end:
                    bin_end = region_end
                bin_name = "{}{}".format(chrom_counter, bin_counter)
                print("{}\t{}\t{}\t{}".format(chrom, bin_start, bin_end, bin_name), file=chrom_regions_f)
                bin_counter += 1
                # move
                bin_start = bin_end
                bin_end += bin_size
### combine regions from multiple chroms into a single file
cmd = "cat"
for f in chrom_regions_files:
    cmd += " " + f
cmd += " > " + bin_regions_file
os.system(cmd)
### remove temp files
os.remove(regions_file_sorted)