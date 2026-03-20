# Step-by-step tutorial: Running epigenome-ssm

This guide will walk you through setting up your environment, downloading a small set of sample ChIP-seq tracks, configuring the pipeline, and running the `epigenome-ssm` annotation workflow locally. 

This tutorial assumes you are using a Unix-based terminal (macOS or Linux).

### Step 1: Set up a working directory and clone the repository
We will start by creating a dedicated directory for this tutorial and cloning the `epigenome-ssm` repository into it.

```bash
# Move to your Desktop (or any other directory of your choice)
cd ~/Desktop

# Create a new working directory for the tutorial and enter it
mkdir ssm-annotation
cd ssm-annotation

# Clone the repository
git clone https://github.com/habibdanesh/epigenome-ssm.git
```

### Step 2: Install Miniforge and the conda environment
To manage the package dependencies properly, we recommend using Miniforge (which provides `conda` and `mamba`). 

1. **Install Miniforge:** If you don't have Miniforge installed, download and install it from the official [Miniforge GitHub repository](https://github.com/conda-forge/miniforge) by following the instructions for your operating system.
2. **Create the environment:** Navigate into the cloned repository to create the dedicated conda environment.

```bash
cd epigenome-ssm

# Create the environment using the provided yaml file
mamba env create -f environment.yaml

# Activate the environment
mamba activate epigenome-ssm

# Move back to the main tutorial directory to run the pipeline
cd ..
```

### Step 3: Prepare the mapping file and download the datasets
The pipeline needs to know where your `bigWig` files are via a tab-separated locator file.

1. **Copy the example configuration files** into your working directory:
```bash
cp epigenome-ssm/examples/in_files_locator.tsv .
cp epigenome-ssm/examples/config.json .
```

2. **Download the BigWig files:** Download the sample bigWig tracks from [here](https://www.icloud.com/iclouddrive/015MUy7zmmVqCoGA4IEpmA98A#input). 

3. **Extract the datasets:** Once downloaded, move `input.zip` to your `~/Desktop/ssm-annotation` folder and unzip it. Ensure the `.bw` files are extracted directly into your working directory (or update the paths in your locator file accordingly).

```bash
mv ~/Downloads/input.zip .
unzip input.zip
```

Your `in_files_locator.tsv` file should look like this (the first 3 columns are required):
```tsv
Epithelial	H3K27ac	input/Epithelial_H3K27ac.pval.bw
Epithelial	H3K36me3	input/Epithelial_H3K36me3.pval.bw
Epithelial	H3K9me3	input/Epithelial_H3K9me3.pval.bw
Myeloid	H3K27ac	input/Myeloid_H3K27ac.pval.bw
Myeloid	H3K36me3	input/Myeloid_H3K36me3.pval.bw
Myeloid	H3K9me3	input/Myeloid_H3K9me3.pval.bw
```
If you decide to place input bigWig files in a different directory, make sure you update the paths in the third column in `in_files_locator.tsv`. 

### Step 4: Verify the JSON config file
You can prepare your `config.json` file using the template below. Make sure the paths are correct.

```json
{
    "in_files_locator": "in_files_locator.tsv",
    
    "bin_size": 25,
    "n_features": 3,
    "max_iter": 10,
    "model_type": "stacked",
    "n_chunks": 10,
    "ref_genome": "hg38",
    
    "n_cores": 12,
    "debug": "False",
    "results_dir": "results/stacked-k3",

    "training": {
        "regions_file": "epigenome-ssm/data/misc/encodePilotRegions.hg38.liftover.bed",
        "chromosomes": "chr1 chr10 chr21",
        "out_dir": "results/processed-data/training"
    },
    
    "annotation": {
        "regions_file": "epigenome-ssm/data/misc/ENCFF356LFX_complement.bed",
        "chromosomes": "chr21",
        "out_dir": "results/processed-data/annotation"
    }
}
```

### Step 5: Run the Master Workflow
Everything is now set up! Execute the entire pipeline using the master `Snakefile`. 

*(Note: We pass `--cores 12` explicitly so Snakemake knows the maximum total parallel threads it is allowed to use globally).*

```bash
snakemake -s epigenome-ssm/Snakefile --configfile config.json --cores 12
```

Once the pipeline successfully finishes, all expected outputs will be stored within the `results` folder.