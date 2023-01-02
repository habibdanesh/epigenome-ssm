import os

### snakemake parameters
regions_file = snakemake.input.regions_file
bin_size = snakemake.params.bin_size
bin_regions_file = snakemake.output.bin_regions_file

### sort the input regions file to make sure chrom bin counter is accurate
regions_file_sorted = '/'.join(bin_regions_file.split('/')[:-1]) + "/regions.sorted.bed.tmp"
os.system("sort -k1V -k2n -k3n {} > {}".format(regions_file, regions_file_sorted))
###
current_chrom = ''
chrom_counter = 0
bin_counter = 1
with open(regions_file_sorted, 'r') as in_f, open(bin_regions_file, 'w') as out_f:
    in_regions = in_f.readlines()
    for line in in_regions:
        columns = line.split("\t")
        chrom = columns[0]
        if chrom != current_chrom:
            current_chrom = chrom
            chrom_counter += 1
            bin_counter = 1
        region_start = int(columns[1])
        region_end = int(columns[2])
        bin_start = region_start
        bin_end = bin_start + bin_size
        while bin_start < region_end:
            if bin_end > region_end:
                bin_end = region_end
            bin_name = "{}_{}".format(chrom_counter, bin_counter)
            print("{}\t{}\t{}\t{}".format(chrom, bin_start, bin_end, bin_name), file=out_f)
            bin_counter += 1
            # move
            bin_start = bin_end
            bin_end += bin_size
###
os.remove(regions_file_sorted)