"""
Get feature values within genebody per bin

Inupt:
    - gene positions
    - Annotation

Output:
    - Feature values within genebody per bin
"""

import os
import sys
import time
import subprocess

import numpy as np
import pandas as pd

k = None
epigenome = None
chromosome = None

region_col_names = ['epigenome', 'chrom', 'base', 'final', 'expression value', 'strand', 'expressed']
region_df = None
annotation_df = None
NUM_BINS = None
gene_region_file = None
TEMP_DIR = None
annotation_result_dir = None


def main(argvs):
    global k
    global epigenome
    global chromosome
    global region_df
    global annotation_df
    global NUM_BINS
    global gene_region_file
    global TEMP_DIR
    global annotation_result_dir

    k = int(argvs[1])
    epigenome = argvs[2]
    chromosome = argvs[3]
    NUM_BINS = int(argvs[4])
    gene_region_file = argvs[5]
    TEMP_DIR = argvs[6]
    annotation_result_dir = argvs[7]

    read_gene_info()
    read_annotation_data()

    feature_df = []
    for region in region_df.itertuples(index=False, name=None):
        feature_df.append(get_genebody_bin_features(region))
    feature_df = pd.DataFrame(feature_df, columns=['f'+str(i) for i in range(1, NUM_BINS*k+1)])
    region_df = region_df.join(feature_df)
    region_df.insert(0, 'epigenome', epigenome) # add epigenome column
    region_df.to_csv(os.path.join(TEMP_DIR, \
        epigenome+'-'+chromosome+'_genebody_bin_features-' + 'k'+str(k) + '.csv'), index=False)


def read_gene_info():
    global region_df
    global epigenome
    global chromosome
    global region_col_names
    global gene_region_file
    col_dtype = {'chrom':np.string_, 'base':np.uint32, 'final':np.uint32, 'expression value':np.float64, \
        'strand':np.int8, 'expressed':np.bool_}
    region_df = pd.read_csv(gene_region_file, \
        names=region_col_names[1:6], header=None, delim_whitespace=True, dtype=col_dtype)
    gexp_std = sorted(region_df['expression value'], reverse=True)[int(len(region_df['expression value'])/3)]
    region_df = region_df[region_df['chrom'] == chromosome]
    region_df['expressed'] = region_df['expression value'] >= gexp_std
    region_df.reset_index(drop=True, inplace=True) # reset the index column


def read_annotation_data():
    global annotation_df
    global k
    global epigenome
    global chromosome
    global annotation_result_dir
    col_names = ['start', 'end'] + ['f'+str(i) for i in range(1,k+1)]
    col_dtype = {'start':np.uint32, 'end':np.uint32}
    for i in range(k):
        col_dtype['f'+str(i+1)] = np.float64
    annotation_file = os.path.join(annotation_result_dir, epigenome, chromosome+'.bed')
    annotation_df = pd.read_csv(annotation_file, \
        names=col_names, header=None, delim_whitespace=True, dtype=col_dtype)


def get_genebody_bin_features(region):
    '''
    region: ['chrom', 'base', 'final', 'expression value', 'strand', 'expressed']
    '''
    global annotation_df
    global k
    global NUM_BINS
    result = np.zeros((NUM_BINS, k))
    search_start = region[1]
    search_end = region[2]
    num_pos_per_bin = int((search_end-search_start) / NUM_BINS)
    num_annotation_rows = annotation_df.shape[0]
    available_bins = 0
    bin_end = search_start # just for initialization
    idx_list = annotation_df.index[(annotation_df.start<=search_start) & (search_start<annotation_df.end)]
    if len(idx_list) == 0:
        return [None for i in range(NUM_BINS*k)]
    row_idx = idx_list[0]
    for bin_idx in range(NUM_BINS):
        bin_start = bin_end
        if bin_idx < (NUM_BINS-1):
            bin_end += num_pos_per_bin
        else: # the last iteration
            bin_end = search_end # because bin_end may not reach search_end in the last itr
        bin_result = np.zeros(k)
        pos_counter = 0
        if bin_idx > 0:
            row_idx -= 1
        while row_idx < num_annotation_rows:
            row = annotation_df.iloc[row_idx,:] # ['start', 'end', 'f1', ..., 'fk']
            start = row[0]
            end = row[1]
            if end <= bin_start:
                row_idx += 1
                continue
            if start >= bin_end:
                break
            if start < bin_start and end > bin_start:
                start = bin_start
            if start < bin_end and end > bin_end:
                end = bin_end
            bin_result += (end - start) * np.array(row[2:], dtype='float64')
            pos_counter += (end - start)
            row_idx += 1
        if pos_counter > 0:
            result[bin_idx,:] = bin_result / pos_counter
            available_bins += 1
    if available_bins < (0.8*NUM_BINS):
        return [None for i in range(NUM_BINS*k)]
    else:
        return result.flatten()


if __name__ == "__main__":
    sys.exit(main(sys.argv))




