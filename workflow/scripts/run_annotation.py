import json
import numpy as np
import ssm

### parameters
in_files = snakemake.input.in_files
model_json = snakemake.input.model_json
feature_npz_files = snakemake.output.feature_npz_files

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

### read input data
if model_type == "stacked":
    npz_objects = []
    for in_f in in_files:
        npz_objects.append(np.load(in_f))
    X_df = np.vstack(([npz_obj["arr_0"] for npz_obj in npz_objects]))
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
            col_npz_objects.append(np.load(track_in_file))
        rows.append(np.hstack(([npz_obj["arr_0"] for npz_obj in col_npz_objects])))
    X_df = np.vstack(([row for row in rows]))
##
print("X_df: {}".format(X_df.shape)) # shape is E x G
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
print("y_m: {}".format(model.y_m.shape)) # shape is K x G

### save the feature tracks
## npz
for k in range(K):
    np.savez_compressed(feature_npz_files[k], model.y_m[k, :]) # [:-4] to remove .npz from end