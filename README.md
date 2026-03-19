# epigenome-ssm

epigenome-ssm is a chromatin state annotation method that uses a non-negative state space model (SSM) and summarizes each genomic position as a vector of continuous chromatin state features, in contrast to the previous approach of assigning a discrete chromatin state label. This continuous approach can summarize complex high-dimensional datasets into a small number of interpretable chromatin state features. Unlike discrete labels, these continuous features preserve the underlying continuous nature of the epigenomic signal tracks, can easily represent varying strengths of a given genomic element, can represent combinatorial elements with multiple types of activities, and are shown to be useful for expressive visualizations because they map complex high-dimensional datasets onto a much smaller number of dimensions.

## Installation

1. **Clone the repository:**
   You need to clone the `epigenome-ssm` repository to your local machine to access the pipeline scripts:
   ```bash
   git clone https://github.com/your-username/epigenome-ssm.git
   cd epigenome-ssm
   ```

2. **Install Conda/Mamba:**
   The preferred way to manage dependencies is using Miniforge, which provides both `conda` and `mamba`. You can download and install it from [here](https://github.com/conda-forge/miniforge).

3. **Create the environment:**
   Once Miniforge is installed, create the dedicated environment using the provided `environment.yaml` file:
   ```bash
   mamba env create -f environment.yaml
   mamba activate epigenome-ssm
   ```

## Running the Workflow

The epigenome-ssm pipeline handles data preparation, model training, and annotation automatically through a single Snakemake workflow. It takes multiple epigenomic signal tracks in bigWig files as input, learns chromatin state features, and generates the final chromatin state feature tracks in both `numpy .npy` and `bigWig` files. See `examples` for step-by-step instructions on how to run the pipeline on a small set of ChIP-seq tracks.

### 1. Prepare your input file locator
Create a tab-separated file detailing the input bigWig files. An example can be found at `examples/in_files_locator.tsv`. The first three columns in this file must be:
1. **Epigenome or cell type ID** (e.g., Brain)
2. **Assay name** (e.g., H3K27ac)
3. **Path to the bigWig file**

### 2. Configure the run parameters
Define the pipeline arguments in a single `config.json` file. An example is provided at `examples/config.json`.

**Key Parameters:**
- `src_dir`: Path to epigenome-ssm repository
- `bin_size`: Resolution (in base pairs) to average the signals.
- `n_features`: Number of continuous chromatin state features to learn.
- `max_iter`: Maximum number of iterations for training the model.
- `n_chunks`: Number of chunks to divide the genome into for parallel efficiency during annotation.
- `ref_genome`: Reference genome assembly (options: `hg19`, `hg38`, `mm10`)
- `regions_file` (training): Path to a BED file containing the regions to train the model on. For faster training, we recommend using a subset of the genome such as the ENCODE pilot regions in `data/misc/encodePilotRegions.hg38.liftover.bed`, which include about 1% of the human genome.
- `regions_file` (annotation): Path to a BED file containing the regions to generate the final chromatin state features. For hg38, the file `data/misc/ENCFF356LFX_complement.bed` can be used, which covers the whole human genome except the `ENCODE Exclusion List Regions`.
- `chromosomes`: This parameter can be used to specify a subset of the chromsomes present in `regions_file` for `training` and `annotation`. 

### 3. Execute the pipeline
Run the whole pipeline by invoking the master `Snakefile`. Note that you must pass the `--cores` flag to let Snakemake know the maximum number of CPU cores it is allowed to use for the entire run:

```bash
snakemake -s path_to_repository/Snakefile --configfile config.json --cores 12
```

## Publications

```
Daneshpajouh, H., Moghul, I., Wiese, K. C., & Libbrecht, M. W. (2025). Pan-cell type continuous chromatin state annotation of all IHEC epigenomes. bioRxiv, 2025-02.

Daneshpajouh, H., Chen, B., Shokraneh, N., Masoumi, S., Wiese, K. C., & Libbrecht, M. W. (2022). Continuous chromatin state feature annotation of the human epigenome. Bioinformatics, 38(11), 3029-3036.
```