import json
import numpy as np
import ssm

### parameters
in_files = snakemake.input.in_files
tracks = snakemake.params.tracks
model_type = snakemake.params.model_type
n_features = snakemake.params.n_features
model_json = snakemake.output.model_json
model_err_file = snakemake.output.model_err_file

### model parameters
LAMBDA_1_L2 = .1
LAMBDA_2_L1 = .0
LAMBDA_2_L2 = .1
LAMBDA_3_L2 = .1
NUM_ITERATIONS = 30
STOP_THRESHOLD = 0.001 # to stop training based on optimization performance

### read input data
if model_type == "stacked":
    npz_objects = []
    for in_f in in_files:
        npz_objects.append(np.load(in_f))
    X_df = np.vstack(([npz_obj["arr_0"] for npz_obj in npz_objects]))
    print("X_df: {}".format(X_df.shape)) # shape is E x G
    num_tracks = X_df.shape[0]
    num_positions = X_df.shape[1]

### data normalization
X_df = np.arcsinh(X_df)

### run ssm training
model = ssm.ssm(E=num_tracks, G=num_positions, K=n_features, \
            lambda_1_l2=LAMBDA_1_L2, lambda_2_l1=LAMBDA_2_L1, lambda_2_l2=LAMBDA_2_L2, lambda_3_l2=LAMBDA_3_L2, \
            positive_state=True, sumone_state=False, positive_em=True, message_passing=True, \
            verbose=False)
model.set_x(np.asmatrix(X_df))
model.optimization(iteration=NUM_ITERATIONS)
error_m = model.error_m

### save the model
model_params = {
    "model_type": model_type,
    'K': n_features,
    'E': num_tracks,
    "nonneg_state": True,
    "sumone_state": False,
    "nonneg_em": True,
    "message_passing": True,
    "lambda_1_l2": LAMBDA_1_L2,
    "lambda_2_l1": LAMBDA_2_L1,
    "lambda_2_l2": LAMBDA_2_L2,
    "lambda_3_l2": LAMBDA_3_L2,
    "theta_m": model.theta_m.tolist(),
    "lambda_m": model.lambda_m.tolist()
}
with open(model_json, 'w') as out_file:
    json.dump(model_params, out_file, indent=4)
## save the errors
np.array(error_m).tofile(model_err_file)