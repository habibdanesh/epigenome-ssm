import argparse
import os
import sys


### parse command-line arguments
arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("-i", "--infiles", required=True, 
                        help="Path to a tab-delimited file containing input file information")
arg_parser.add_argument("-r", "--regions", required=True, 
                        help="Path to a BED file containing genomic regions to be included in the analysis")
arg_parser.add_argument("-b", "--binsize", type=int, required=False, default=200, 
                        help="Bin size (or resolution) to average the signals")
arg_parser.add_argument("-o", "--outdir", required=False, 
                        help="Output directory to save the generated files")
arg_parser.add_argument("-p", "--cores", type=int, required=False, default=1, 
                        help="Number of CPU cores")
arg_parser.add_argument("-n", "--dryrun", required=False, action="store_true", 
                        help="Dry run")
args = arg_parser.parse_args()


def run_pipeline():
    ### parameters
    root_path = '/'.join(os.path.abspath(sys.argv[0]).split('/')[:-1])
    in_files = args.infiles
    regions_file = args.regions
    bin_size = args.binsize
    out_dir = args.outdir
    n_cores = args.cores
    dry_run = args.dryrun
    ### run
    cmd = "snakemake --snakefile {}/workflow/Snakefile_prepare --cores {}".format(root_path, n_cores)
    cmd += " --config root_path={} in_files={} regions_file={} bin_size={} out_dir={}".format(
                        root_path, in_files, regions_file, bin_size, out_dir)
    if dry_run:
        cmd += " -n"
    print(cmd)
    print('')
    os.system(cmd)


if __name__ == "__main__":
    run_pipeline()