import json
import numpy as np
import ssm

### parameters
in_files = snakemake.input.in_files
tracks = snakemake.params.tracks
model_type = snakemake.params.model_type
n_features = snakemake.params.n_features
lambda_1_l2 = snakemake.params.lambda_1_l2
lambda_2_l1 = snakemake.params.lambda_2_l1
lambda_2_l2 = snakemake.params.lambda_2_l2
lambda_3_l2 = snakemake.params.lambda_3_l2
n_rand_init = snakemake.params.n_rand_init
max_iter = snakemake.params.max_iter
model_json = snakemake.output.model_json
model_err_file = snakemake.output.model_err_file

### read input data
## stacked model
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
            lambda_1_l2=lambda_1_l2, lambda_2_l1=lambda_2_l1, lambda_2_l2=lambda_2_l2, lambda_3_l2=lambda_3_l2, \
            positive_state=True, sumone_state=False, positive_em=True, message_passing=True, \
            verbose=False)
##
model.set_x(np.asmatrix(X_df))
## multiple random initialization
np.random.seed()
seeds = []
errors = []
for t in range(n_rand_init):
    s = np.random.randint(1000000)
    model.re_init(seed=s)
    seeds.append(s)
    errors.append(model.total_error())
best_seed = seeds[np.argmin(errors)]
model.re_init(seed=best_seed)
##
model.optimization(iteration=max_iter)
error_m = model.error_m

### save the model
model_params = {
    "seed": model.seed,
    "model_type": model_type,
    'K': n_features,
    'E': num_tracks,
    "nonneg_state": True,
    "sumone_state": False,
    "nonneg_em": True,
    "message_passing": True,
    "lambda_1_l2": lambda_1_l2,
    "lambda_2_l1": lambda_2_l1,
    "lambda_2_l2": lambda_2_l2,
    "lambda_3_l2": lambda_3_l2,
    "theta_m": model.theta_m.tolist(),
    "lambda_m": model.lambda_m.tolist()
}
with open(model_json, 'w') as out_file:
    json.dump(model_params, out_file, indent=4)
## save the errors
np.array(error_m).tofile(model_err_file)