# epigenome-ssm

epigenome-ssm is a chromatin state annotation method that uses a non-negative state space model (SSM) and summarizes each genomic position as a vector of continuous chromatin state features, in contrast to the previous approach of assigning a discrete chromatin state label. This continuous approach can summarize complex high-dimensional datasets into a small number of interpretable chromatin state features. Unlike discrete labels, these continuous features preserve the underlying continuous nature of the epigenomic signal tracks, can easily represent varying strengths of a given genomic element, can represent combinatorial elements with multiple types of activities, and are shown to be useful for expressive visualizations because they map complex high-dimensional datasets onto a much smaller number of dimensions.

### Dependencies
```
Snakemake
NumPy
```

## Publications

```
Daneshpajouh, Habib, et al. "Continuous chromatin state feature annotation of the human epigenome." Bioinformatics 38.11 (2022): 3029-3036.
```