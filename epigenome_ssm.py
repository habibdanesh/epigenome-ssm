import argparse
import os
import sys


### parse command-line arguments
arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("task", help="Task")
arg_parser.add_argument("-i", "--indir", required=True, 
                        help="Path to the directory that contains prepared input data")
arg_parser.add_argument("--n_bins_file", required=False, 
                        help="Path to BED file containing number of bins in each region")
arg_parser.add_argument("-b", "--binsize", type=int, required=False, default=200, 
                        help="Bin size (resolution) to average the signals")
arg_parser.add_argument("--custom_element_file", required=False,
                        help="Path to tab-separated file containing elements coordinates and other info")
arg_parser.add_argument("-m", "--modeltype", required=False, default="stacked", choices=["stacked", "concatenated"], 
                        help="Type of model")
arg_parser.add_argument("-k", "--nfeatures", type=int, required=False, default=5, 
                        help="Number of chromatin state features")
arg_parser.add_argument("--n_chunks", type=int, required=False, default=10, 
                        help="Number of chunks to divide the genome into, for efficiency")
arg_parser.add_argument("-o", "--outdir", required=False, 
                        help="Output directory to save the generated files")
arg_parser.add_argument("--model", required=False, 
                        help="Path to an existing model")
arg_parser.add_argument("--ref_genome", required=False, default='hg38', choices=['hg19', 'hg38', 'mm10'],
                        help="Reference genome assembly for generating the final feature bigWig files")
arg_parser.add_argument("-x", "--x_array", required=False, 
                        help="Path to an .npy file containing the input array to the model.")
arg_parser.add_argument("-p", "--cores", type=int, required=False, default=4, 
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
    n_bins_file = args.n_bins_file
    bin_size = args.binsize
    custom_element_file = args.custom_element_file
    model_type = args.modeltype
    n_features = args.nfeatures
    n_chunks = args.n_chunks
    out_dir = args.outdir
    model = args.model
    ref_genome = args.ref_genome
    x_array = args.x_array
    n_cores = args.cores
    dry_run = args.dryrun
    debug = args.debug
    ### train
    if task == "train":
        cmd = "snakemake --snakefile {}/workflow/Snakefile_train".format(root_path)
        if n_cores == -1: # use all available cores
            cmd += " --cores"
        else:
            cmd += " --cores {}".format(n_cores)
        cmd += " --config root_path={} in_dir={} model_type={} n_features={} out_dir={} model={} x_array={} debug={}".format(
                            root_path, in_dir, model_type, n_features, out_dir, model, x_array, debug)
        if dry_run:
            cmd += " -n"
        print(cmd)
        print('')
        os.system(cmd)
    ### annotate
    elif task == "annotate":
        cmd = "snakemake --snakefile {}/workflow/Snakefile_annotate".format(root_path)
        if n_cores == -1: # use all available cores
            cmd += " --cores"
        else:
            cmd += " --cores {}".format(n_cores)
        cmd += " --config root_path={} in_dir={} out_dir={} model={} n_chunks={} ref_genome={} debug={}".format(
                            root_path, in_dir, out_dir, model, n_chunks, ref_genome, debug)
        if dry_run:
            cmd += " -n"
        print(cmd)
        print('')
        os.system(cmd)
    ### evaluate
    elif task == "evaluate":
        cmd = "snakemake --snakefile {}/workflow/Snakefile_evaluate".format(root_path)
        if n_cores == -1: # use all available cores
            cmd += " --cores"
        else:
            cmd += " --cores {}".format(n_cores)
        cmd += " --config root_path={} in_dir={} n_bins_file={} bin_size={} custom_element_file={} out_dir={} debug={}".format(
                            root_path, in_dir, n_bins_file, bin_size, custom_element_file, out_dir, debug)
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