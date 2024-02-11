import json
import psutil
import numpy as np
import ssm
import time

### parameters
in_files_locator = snakemake.input.in_files_locator
model_json = snakemake.input.model_json
total_bins = snakemake.params.total_bins
max_bins = snakemake.params.max_bins
chunk_out_npy = snakemake.output.chunk_out_npy

def print_mem_usage():
    print("RAM Used (GB): {}\n".format(psutil.Process().memory_info().rss / (1024 * 1024 * 1024)))

### track input files
in_files = []
with open(in_files_locator, 'r') as in_f_locator:
    for line in in_f_locator.readlines():
        in_files.append(line.strip("\n").split("\t")[2])

### model parameters
with open(model_json, 'r') as model_f:
    model_params = json.load(model_f)
model_type = model_params["model_type"]
K = model_params['K']
E = model_params['E']
nonneg_state = model_params["nonneg_state"]
sumone_state = model_params["sumone_state"]
nonneg_em = model_params["nonneg_em"]
message_passing = model_params["message_passing"]
lambda_1_l2 = model_params["lambda_1_l2"]
lambda_2_l1 = model_params["lambda_2_l1"]
lambda_2_l2 = model_params["lambda_2_l2"]
lambda_3_l2 = model_params["lambda_3_l2"]
theta_m = np.array(model_params["theta_m"], dtype=np.single)
lambda_m = np.array(model_params["lambda_m"], dtype=np.single)
epigenomes = model_params["epigenomes"]
assays = model_params["assays"]
model_params = None

### get the region info
chunk_idx = int(chunk_out_npy[:-len(".npy")].split('_')[-1])
bin_start = chunk_idx * max_bins
bin_end = (chunk_idx+1) * max_bins
if bin_end > total_bins:
    bin_end = total_bins
chunk_nbins = bin_end - bin_start
print(f"### chunk {chunk_idx},  bin {bin_start}:{bin_end}")

print_mem_usage()
### read input data
X = np.empty((E, chunk_nbins), dtype=np.single)
## stacked model
if model_type == "stacked":
    for idx, in_f in enumerate(in_files):
        start_time = time.time()
        X[idx, :] = np.load(in_f, mmap_mode='r')[bin_start:bin_end]

## concatenated model
if model_type == "concatenated":
    rows = []
    for assay in assays:
        col_npy_objects = []
        for epigenome in epigenomes:
            track = "{}_{}".format(epigenome, assay)
            # find track input file
            track_in_file = None
            for in_f in in_files:
                if in_f.split('/')[-1].startswith(track):
                    track_in_file = in_f
                    break
            assert track_in_file != None
            col_npy_objects.append(np.load(track_in_file, mmap_mode='r'))
        rows.append(np.hstack(([npy_obj[bin_start:bin_end] for npy_obj in col_npy_objects])))
    X = np.vstack(([row for row in rows]))
##
print("X: {}".format(X.shape)) # shape is E x n_bins
num_tracks = X.shape[0]
num_positions = X.shape[1]

### data normalization
X = np.arcsinh(X)
### run ssm annotation
model = ssm.ssm(E=num_tracks, G=num_positions, K=K, \
            lambda_1_l2=lambda_1_l2, lambda_2_l1=lambda_2_l1, lambda_2_l2=lambda_2_l2, lambda_3_l2=lambda_3_l2, \
            positive_state=nonneg_state, sumone_state=sumone_state, positive_em=nonneg_em, message_passing=message_passing, \
            verbose=False)
model.set_x(X)
model.set_theta(theta_m)
model.set_lambda(lambda_m)
model.update_state()
print("y_m: {}\n".format(model.y_m.shape)) # shape is K x n_bins

### save the feature tracks
np.save(chunk_out_npy, model.y_m)