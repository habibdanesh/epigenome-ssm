import sys
import numpy as np
import pandas as pd

def filter_annotations(args):
    gene_annotation_f = args[1]
    out_file = args[2]
    ### reads genes annotation
    gene_df = pd.read_csv(gene_annotation_f, sep="\t", header=None, skiprows=7,
                            names=["chrom", "type", "start", "end", "strand", "attributes"], 
                            usecols=[0, 2, 3, 4, 6, 8],
                            dtype={"chrom":np.string_, "type":np.string_, "start":np.double, "end":np.double, 
                            "strand":np.string_, "attributes":np.string_})
    ### to handle NA values, read the start/end cloumns as double first, then convert to int (because of a known issue in pandas)
    gene_df.dropna(inplace=True)
    gene_df["start"] = gene_df["start"].astype(pd.Int64Dtype())
    gene_df["end"] = gene_df["end"].astype(pd.Int64Dtype())
    ### filter
    gene_df = gene_df[(gene_df["type"] == "gene") & 
                    (gene_df["attributes"].str.contains("gene_type=protein_coding"))
                    ]
    ### convert from 1-based to 0-based index
    gene_df.start = gene_df.start.to_numpy() - 1
    gene_df.end = gene_df.end.to_numpy() - 1
    ### save
    gene_df.to_csv(out_file, sep="\t", index=False)

if __name__ == "__main__":
    filter_annotations(sys.argv)