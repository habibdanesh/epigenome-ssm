import argparse
import os
import sys


### parse command-line arguments
arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("task", help="Task")
arg_parser.add_argument("-i", "--indir", required=True, 
                        help="Path to the directory that contains prepared input data")
arg_parser.add_argument("-m", "--modeltype", required=False, default="stacked", choices=["stacked", "concatenated"], 
                        help="Type of model")
arg_parser.add_argument("-k", "--nfeatures", type=int, required=False, default=5, 
                        help="Number of chromatin state features")
arg_parser.add_argument("-o", "--outdir", required=False, 
                        help="Output directory to save the generated files")
arg_parser.add_argument("-p", "--cores", type=int, required=False, default=1, 
                        help="Number of CPU cores")
arg_parser.add_argument("--debug", required=False, action="store_true", 
                        help="Debug mode")
arg_parser.add_argument("-n", "--dryrun", required=False, action="store_true", 
                        help="Dry run")
args = arg_parser.parse_args()


def run_pipeline():
    ### parameters
    task = args.task
    root_path = '/'.join(os.path.abspath(sys.argv[0]).split('/')[:-1])
    in_dir = args.indir
    model_type = args.modeltype
    n_features = args.nfeatures
    out_dir = args.outdir
    n_cores = args.cores
    dry_run = args.dryrun
    debug = args.debug
    ### train
    if task == "train":
        cmd = "snakemake --snakefile {}/workflow/Snakefile_train --cores {}".format(root_path, n_cores)
        cmd += " --config root_path={} in_dir={} model_type={} n_features={} out_dir={} debug={}".format(
                            root_path, in_dir, model_type, n_features, out_dir, debug)
        if dry_run:
            cmd += " -n"
        print(cmd)
        print('')
        os.system(cmd)
    ### annotate
    elif task == "annotate":
        cmd = "snakemake --snakefile {}/workflow/Snakefile_annotate --cores {}".format(root_path, n_cores)
        cmd += " --config root_path={} in_dir={} out_dir={} debug={}".format(
                            root_path, in_dir, out_dir, debug)
        if dry_run:
            cmd += " -n"
        print(cmd)
        print('')
        os.system(cmd)
    ###
    else:
        print("Undefined task!")


if __name__ == "__main__":
    run_pipeline()