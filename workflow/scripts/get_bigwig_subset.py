import os
import numpy as np

### snakemake parameters
in_files = snakemake.input.in_files
regions_files = snakemake.input.regions_files
bw_avg_over_bed_src = snakemake.input.bw_avg_over_bed_src
out_npz_files = snakemake.output.out_npz_files

for in_f, out_npz_f in zip(in_files, out_npz_files):
    temp_out_files = []
    for i, regions_file in enumerate(regions_files):
        out_bg_f = out_npz_f[:-len(".npz")] + ".bg"
        out_tab_f = "{}.tab.tmp".format(out_bg_f)
        out_bed_f = "{}.bed.tmp".format(out_bg_f)
        temp_out_f = "{}.{}.tmp".format(out_bg_f, i)
        temp_out_files.append(temp_out_f)
        cmd = "{} {} {} {} -bedOut={}".format(
                bw_avg_over_bed_src, in_f, regions_file, out_tab_f, out_bed_f)
        os.system(cmd)
        ### sort
        cmd = "sort -k1V -k2n -k3n {} > {}".format(out_bed_f, temp_out_f)
        os.system(cmd)
        ### remove temp files
        os.remove(out_tab_f)
        os.remove(out_bed_f)
    ### combine results from multiple out files
    cmd = "cat"
    for temp_out_f in temp_out_files:
        cmd += " " + temp_out_f
    cmd += " > " + out_bg_f
    os.system(cmd)
    ### remove temp files
    for temp_out_f in temp_out_files:
        os.remove(temp_out_f)
    ### convert bedgraph to npz
    a = np.loadtxt(out_bg_f, dtype=np.single, delimiter="\t", usecols=4)
    np.savez_compressed(out_npz_f[:-4], a) # [:-4] to remove .npz from end
    os.remove(out_bg_f)
### remove temp regions files
for regions_file in regions_files:
    os.remove(regions_file)