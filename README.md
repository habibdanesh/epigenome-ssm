# epigenome-ssm

epigenome-ssm is a chromatin state annotation method that uses a non-negative state space model (SSM) and summarizes each genomic position as a vector of continuous chromatin state features, in contrast to the previous approach of assigning a discrete chromatin state label. This continuous approach can summarize complex high-dimensional datasets into a small number of interpretable chromatin state features. Unlike discrete labels, these continuous features preserve the underlying continuous nature of the epigenomic signal tracks, can easily represent varying strengths of a given genomic element, can represent combinatorial elements with multiple types of activities, and are shown to be useful for expressive visualizations because they map complex high-dimensional datasets onto a much smaller number of dimensions.

## Dependencies
```
Snakemake
NumPy
```

## Overview of the epigenome-ssm workflow

The epigenome-ssm workflow has three main steps:

1. **Data preparation**
2. **Training**
3. **Annotation**

---

### Step 1: Data preparation

Take the input BigWig files and bin them (by averaging the signal) at the desired resolution (e.g., 25 bp).  
**Important:** Training and annotation data should be prepared separately using the `prepare_data.py` script.

```
python prepare_data.py -i infiles.tsv -r regions.bed
```

- `infiles.tsv` is a tab-delimited file in which each row represents an input BigWig file with three columns:
    1. a unique cell type name / ID
    2. assay name
    3. path to the BigWig file
- `regions.bed` is a BED file containing the genomic regions from which to extract data.

For faster training, we recommend using a subset of the genome such as the ENCODE pilot regions in `data/misc/encodePilotRegions.hg38.liftover.bed`.

To prepare data for annotating the full genome, you can use `data/misc/ENCFF356LFX_complement.bed`.

See `prepare_data.py` for a full list of available options and parameters.

---

### Step 2: Train an epigenome-ssm model

Once the training data have been prepared, you can train an epigenome-ssm model with:

```
python epigenome_ssm.py train
```

See `epigenome_ssm.py` for optional parameters (e.g., number of features).

---

### Step 3: Annotate the genome with a trained model

After training, you can generate continuous feature annotation tracks from a trained epigenome-ssm model with:

```
python epigenome_ssm.py annotate
```

This command applies the trained model to the prepared annotation data and outputs continuous feature tracks. Again, see `epigenome_ssm.py` for a full list of optional parameters (e.g., output paths, and other runtime settings).

## Publications

```
Daneshpajouh, H., Moghul, I., Wiese, K. C., & Libbrecht, M. W. (2025). Pan-cell type continuous chromatin state annotation of all IHEC epigenomes. bioRxiv, 2025-02.

Daneshpajouh, H., Chen, B., Shokraneh, N., Masoumi, S., Wiese, K. C., & Libbrecht, M. W. (2022). Continuous chromatin state feature annotation of the human epigenome. Bioinformatics, 38(11), 3029-3036.
```