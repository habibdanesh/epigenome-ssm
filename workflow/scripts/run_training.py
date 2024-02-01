import os
import json
import numpy as np
import ssm

### parameters
in_files = snakemake.input.in_files
epigenomes = snakemake.params.epigenomes
assays = snakemake.params.assays
tracks = snakemake.params.tracks
model_type = snakemake.params.model_type
n_features = snakemake.params.n_features
lambda_1_l2 = snakemake.params.lambda_1_l2
lambda_2_l1 = snakemake.params.lambda_2_l1
lambda_2_l2 = snakemake.params.lambda_2_l2
lambda_3_l2 = snakemake.params.lambda_3_l2
n_rand_init = snakemake.params.n_rand_init
max_iter = snakemake.params.max_iter
min_improvement = snakemake.params.min_improvement
existing_model = snakemake.params.existing_model
model_json = snakemake.output.model_json

def save_model(itr):
    ### save the model
    model_params = {
        "seed": model.seed,
        "model_type": model_type,
        'K': n_features,
        'E': num_tracks,
        "nonneg_state": True,
        "sumone_state": False,
        "nonneg_em": True,
        "message_passing": False,
        "lambda_1_l2": lambda_1_l2,
        "lambda_2_l1": lambda_2_l1,
        "lambda_2_l2": lambda_2_l2,
        "lambda_3_l2": lambda_3_l2,
        "theta_m": model.theta_m.tolist(),
        "lambda_m": model.lambda_m.tolist(),
        "error_m": model.error_m,
        "improve_m": model.improve_m,
        "opt_time_m": model.opt_time_m,
        "assays": assays,
        "epigenomes": epigenomes
    }
    with open(f"{model_json[:-5]}.itr{itr}.json", 'w') as out_file:
        json.dump(model_params, out_file, indent=4)
    with open(model_json, 'w') as out_file:
        json.dump(model_params, out_file, indent=4)

### read input data
## stacked model
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

### run ssm training
model = ssm.ssm(E=num_tracks, G=num_positions, K=n_features, \
            lambda_1_l2=lambda_1_l2, lambda_2_l1=lambda_2_l1, lambda_2_l2=lambda_2_l2, lambda_3_l2=lambda_3_l2, \
            positive_state=True, sumone_state=False, positive_em=True, message_passing=False, \
            verbose=False)
##
model.set_x(X_df)
if os.path.isfile(existing_model):
    # Load model parameters
    with open(existing_model, 'r') as model_f:
        model_params = json.load(model_f)
        model.set_theta(model_params["theta_m"])
        model.set_lambda(model_params["lambda_m"])
        model.set_error_m(model_params["error_m"])
        model.set_improve_m(model_params["improve_m"])
        model.set_opt_time(model_params["opt_time_m"])
else:
    # multiple random initialization
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
## optimization
iter_step = 5
model.optimization(iteration=iter_step)
iter_passed = iter_step
save_model(iter_passed)
print("# iterations passed: {}".format(iter_passed))
while iter_passed < max_iter:
    model.optimization(iteration=iter_step, carryon=True)
    iter_passed += iter_step
    save_model(iter_passed)
    print("# iterations passed: {}".format(iter_passed))
    if np.mean(model.improve_m[-iter_step:]) < min_improvement: # if training is not progressing
        break
print("# iterations: {}".format(iter_passed))
save_model(iter_passed)