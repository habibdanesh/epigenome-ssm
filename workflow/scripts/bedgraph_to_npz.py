import os
import numpy as np

### snakemake parameters
in_files = snakemake.input.in_files
out_files = snakemake.output.out_files

for in_f, out_f in zip(in_files, out_files):
    a = np.loadtxt(in_f, dtype=np.single, delimiter="\t", usecols=4)
    np.savez_compressed(out_f[:-4], a) # [:-4] to remove .npz from end
    os.remove(in_f)