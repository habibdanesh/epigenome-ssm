import os

### snakemake parameters
in_files = snakemake.input.in_files
regions_file = snakemake.input.regions_file
bw_avg_over_bed_src = snakemake.input.bw_avg_over_bed_src
out_files = snakemake.output.out_files

for in_f, out_f in zip(in_files, out_files):
    out_tab_f = "{}.tab".format(out_f)
    out_bed_f = "{}.tmp".format(out_f)
    cmd = "{} {} {} {} -bedOut={}".format(
            bw_avg_over_bed_src, in_f, regions_file, out_tab_f, out_bed_f)
    print(cmd)
    os.system(cmd)
    print('')
    ### sort
    cmd = "sort -k1V -k2n -k3n {} > {}".format(out_bed_f, out_f)
    print(cmd)
    os.system(cmd)
    print('')
    ### remove temp files
    os.remove(out_tab_f)
    os.remove(out_bed_f)