'''
Example command to run this script:
python annotation_node.py E116 chr20 0 142566 3 200 DNase,H2A.Z,H3K27ac,H3K27me3,H3K36me3,H3K4me1,H3K4me2,H3K4me3,H3K79me2,H3K9ac,H3K9me3,H4K20me1 ../../../data/f-data/blacklist-rm-data/ temp output output/annotation-result-k3 0
'''
# import packages
import os
import sys
import time
import subprocess
import re
import json
import codecs

import numpy as np

import models.state_space_model.ssm as ssm


def main(argvs):
    # parse arguments
    epigenome = argvs[1]
    chromosome = argvs[2]
    start_idx = int(argvs[3])
    end_idx = int(argvs[4])
    k = int(argvs[5])
    resolution = int(argvs[6])
    assays = list((argvs[7]).split(','))
    annotation_data_dir = argvs[8]
    TEMP_DIR = argvs[9]
    OUTPUT_DIR = argvs[10]
    annotation_result_dir = argvs[11]
    pid = int(argvs[12])
    print(epigenome, chromosome, start_idx, end_idx, pid)
    # load model parameters
    obj_text = codecs.open(os.path.join(OUTPUT_DIR, 'model_k'+str(k)+'.json'), 'r', \
        encoding='utf-8').read()
    model_param = json.loads(obj_text)
    E = model_param['E']
    nonneg_state = model_param['nonneg_state']
    sumone_state = model_param['sumone_state']
    nonneg_em = model_param['nonneg_em']
    message_passing = model_param['message_passing']
    lambda_1_l2 = model_param['lambda_1_l2']
    lambda_2_l1 = model_param['lambda_2_l1']
    lambda_2_l2 = model_param['lambda_2_l2']
    lambda_3_l2 = model_param['lambda_3_l2']
    theta_m = np.asmatrix(model_param['theta_m'])
    lambda_m = np.asmatrix(model_param['lambda_m'])
    # read data
    resol_fs = []
    for assay in assays:
        resol_fs.append(open(os.path.join(annotation_data_dir, epigenome + '-' + assay + \
            '/resolution-' + str(resolution) + 'bp/', chromosome + '.bed'), 'r'))
    region_length = end_idx - start_idx + 1
    # go through each position in the region and see if data is available for all the assays
    data_array = np.array([[-1.0 for j in range(region_length)] \
        for i in range(E)], dtype=float) # data array, with shape: E x region_length
    for assay_idx in range(E):
        resol_f = resol_fs[assay_idx]
        while True:
            line = resol_f.readline()
            if not line:
                break
            inf = line.strip('\n').split()
            pos_idx = int(inf[0])
            if pos_idx < start_idx:
                continue
            if pos_idx > end_idx:
                break
            signal = float(inf[1])
            signal = np.arcsinh(signal) # data transformation
            data_array[assay_idx, pos_idx-start_idx] = signal
    # remove columns of the data matrix where there is missing values
    idx_list_remove = [] # list of positions to remove
    for pos_idx in range(region_length):
        for assay_idx in range(E):
            if data_array[assay_idx, pos_idx] == -1:
                idx_list_remove.append(pos_idx)
                break
    data_array = np.delete(data_array, idx_list_remove, axis=1)
    idx_list_keep = [idx+start_idx for idx in range(region_length) \
        if idx not in idx_list_remove] # list of positions to keep
    g = len(idx_list_keep)
    # apply ssm
    model = ssm.ssm(E=E, G=g, K=k, \
        lambda_1_l2=lambda_1_l2, lambda_2_l1=lambda_2_l1, lambda_2_l2=lambda_2_l2, lambda_3_l2=lambda_3_l2, \
        positive_state=nonneg_state, sumone_state=sumone_state, positive_em=nonneg_em, message_passing=message_passing, \
        verbose=False)
    model.set_x(np.asmatrix(np.asanyarray(data_array)))
    model.set_theta(theta_m)
    model.set_lambda(lambda_m)
    model.update_state()
    # save the annotation result
    with open(os.path.join(annotation_result_dir, \
        epigenome + chromosome + '_' + str(pid) + '.bed'), 'w') as annotation_f:
        index_counter = 0
        for t in range(model.G):
            print('{} {} '.format(idx_list_keep[index_counter] * resolution, \
                (idx_list_keep[index_counter]+1) * resolution), end='', file=annotation_f)
            for k in range(model.K):
                print('{} '.format(model.y_m[k, t]), end='', file=annotation_f)
            print('\n', end='', file=annotation_f)
            index_counter += 1


if __name__ == '__main__':
    sys.exit(main(sys.argv))
