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
conda activate epigenome-ssm

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

2. **Download the BigWig files:** Open the `in_files_locator.tsv` file and open the links in the 4th column using your web browser. Download all 6 `.bw` files and place them directly in your `~/Desktop/ssm-annotation` folder. 
   
Your `in_files_locator.tsv` file should look like this:
```tsv
Epithelial	H3K27ac	Epithelial_H3K27ac.pval.bw	https://www.icloud.com/iclouddrive/058x9wp6wm_H7ClRsvNemC64A#Epithelial_H3K27ac
Epithelial	H3K36me3	Epithelial_H3K36me3.pval.bw	https://www.icloud.com/iclouddrive/028pkgp4TfND5DfLFIWdm5Veg#Epithelial_H3K36me3
Epithelial	H3K9me3	Epithelial_H3K9me3.pval.bw	https://www.icloud.com/iclouddrive/03fR9d2A5AwFl2uFRUjznihGg#Epithelial_H3K9me3
Myeloid	H3K27ac	Myeloid_H3K27ac.pval.bw	https://www.icloud.com/iclouddrive/0f2VObl2yMOVM_tneMOzn8o_Q#Myeloid_H3K27ac
Myeloid	H3K36me3	Myeloid_H3K36me3.pval.bw	https://www.icloud.com/iclouddrive/090GF5VpuI1QtqZaF7y4HGq2w#Myeloid_H3K36me3
Myeloid	H3K9me3	Myeloid_H3K9me3.pval.bw	https://www.icloud.com/iclouddrive/091DorxN9npaixpvX8R_RbbBA#Myeloid_H3K9me3
```
If you decide to place input bigWig files in a different directory, make sure you update the paths in the third column in `in_files_locator.tsv`. 

### Step 4: Verify the JSON config file
You can prepare your `config.json` file using the template below. Make sure the paths are correct.

```json
{
    "src_dir": "epigenome-ssm",
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