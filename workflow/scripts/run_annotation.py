import json
import numpy as np
import ssm

### parameters
in_files_locator = snakemake.input.in_files_locator
n_bins_file = snakemake.input.n_bins_file
model_json = snakemake.input.model_json
chunk_npz_file = snakemake.output.chunk_npz_file

def print_mem_usage():
    import psutil
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
theta_m = np.asmatrix(model_params["theta_m"])
lambda_m = np.asmatrix(model_params["lambda_m"])
epigenomes = model_params["epigenomes"]
assays = model_params["assays"]
model_params = None

### get the region info
chunk_idx = int(chunk_npz_file[:-len(".npz")].split('_')[-1])
with open(n_bins_file, 'r') as n_bins_f:
    lines = n_bins_f.readlines()
bin_start = 0
bin_end = int(lines[0].strip("\n").split("\t")[-1])
for i in range(1, chunk_idx+1):
    bin_start = bin_end
    bin_end += int(lines[i].strip("\n").split("\t")[-1])
print_mem_usage()
### read input data
X_df = np.matrix([], dtype=np.single)
## stacked model
if model_type == "stacked":
    npz_objects = []
    for in_f in in_files:
        npz_objects.append(np.load(in_f, mmap_mode='r'))
        npz_obj = np.load(in_f, mmap_mode='r')
        X_df = np.vstack([X_df, npz_obj["arr_0"][bin_start:bin_end]]) \
                            if X_df.size else npz_obj["arr_0"][bin_start:bin_end]
        #print("X_df size (MB): {}".format(X_df.size * X_df.itemsize / (1024 * 1024)))
        #print_mem_usage()
    print_mem_usage()
## concatenated model
if model_type == "concatenated":
    rows = []
    for assay in assays:
        col_npz_objects = []
        for epigenome in epigenomes:
            track = "{}_{}".format(epigenome, assay)
            # find track input file
            track_in_file = None
            for in_f in in_files:
                if in_f.split('/')[-1].startswith(track):
                    track_in_file = in_f
                    break
            assert track_in_file != None
            col_npz_objects.append(np.load(track_in_file, mmap_mode='r'))
        rows.append(np.hstack(([npz_obj["arr_0"][bin_start:bin_end] for npz_obj in col_npz_objects])))
    X_df = np.vstack(([row for row in rows]))
##
print("X_df: {}".format(X_df.shape)) # shape is E x n_bins
num_tracks = X_df.shape[0]
num_positions = X_df.shape[1]
### data normalization
X_df = np.arcsinh(X_df)
### run ssm annotation
model = ssm.ssm(E=num_tracks, G=num_positions, K=K, \
            lambda_1_l2=lambda_1_l2, lambda_2_l1=lambda_2_l1, lambda_2_l2=lambda_2_l2, lambda_3_l2=lambda_3_l2, \
            positive_state=nonneg_state, sumone_state=sumone_state, positive_em=nonneg_em, message_passing=message_passing, \
            verbose=False)
model.set_x(np.asmatrix(X_df))
model.set_theta(theta_m)
model.set_lambda(lambda_m)
model.update_state()
print("y_m: {}\n".format(model.y_m.shape)) # shape is K x n_bins

### save the feature tracks
np.savez_compressed(chunk_npz_file, model.y_m)