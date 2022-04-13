'''
Example command to run this script:
python head.py train [3,5,8,10,12,15] 0 0 0 0 [E003,E017,E114,E116,E122,E123,E127,E128] all_assays all_chromosomes 200 f ../../../data/f-data/processed-pilot-data/ ../../../data/f-data/blacklist-rm-data/
'''
# import packages
import os
import sys
import time
import subprocess
import re
import json
import codecs
import platform, psutil

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression

import util.utility as util
import util.myconfig as myconfig
import models.state_space_model.ssm as ssm


# global parameters
task = None # the main task the script is asked to do, e.g. train, annotate, evaluate
configs = None
K = None # number of features/labels/states
E = None # number of assays

nonneg_state = None # non-negativity constraint on states
sumone_state = None # constraint on states to sum to one
nonneg_em = None # non-negativity constraint on emission parameters
message_passing = None

epigenomes = None
assays = None
chromosomes = None
resolution = None # bin size
data_type = None # fold-change (f) or p-value (p)
training_data_dir = None
annotation_data_dir = None

TEMP_DIR = 'temp'
temp_dir_k = None
OUTPUT_DIR = 'output'

# model parameters
LAMBDA_1_L2 = .1
LAMBDA_2_L1 = .0
LAMBDA_2_L2 = .1
LAMBDA_3_L2 = .1
TRAINING_ITR = 30
STOP_THRESHOLD = 0.001 # to stop training based on optimization performance
WINDOW_SIZE = 100 # chunks of chromosomes

# evaluation parameters
NUM_EVAL_BINS = 10
NUM_EVAL_TRIALS = 10


def main(argvs):
    global task
    global configs
    global K
    task = argvs[1]
    if task == 'train':
        init(argvs)
        print('\n### training..\n')
        run_training()
        #run_training_in_parallel()
    elif task == 'annotate':
        init(argvs)
        print('\n### training..\n')
        run_training()
        #run_training_in_parallel()
        print('\n### annotation..\n')
        run_annotation()
    elif task == 'evaluate':
        init(argvs)
        print('\n### training..\n')
        run_training()
        #run_training_in_parallel()
        print('\n### annotation..\n')
        run_annotation()
        print('\n### evaluation..\n')
        run_evaluation()
    elif task == 'visualize':
        viz(argvs)
    else:
        print('ERROR: Task is not defined.')
        return


def init(argvs):
    global task
    global util
    global configs
    global K
    global E
    global nonneg_state
    global sumone_state
    global nonneg_em
    global message_passing
    global epigenomes
    global assays
    global chromosomes
    global resolution
    global data_type
    global training_data_dir
    global annotation_data_dir
    K = [int(i) for i in argvs[2][1:len(argvs[2])-1].split(',')]
    nonneg_state = int(argvs[3])
    sumone_state = int(argvs[4])
    nonneg_em = int(argvs[5])
    message_passing = int(argvs[6])
    epigenomes = argvs[7][1:len(argvs[7])-1].split(',')
    assays = argvs[8]
    chromosomes = argvs[9]
    resolution = int(argvs[10])
    data_type = argvs[11]
    training_data_dir = argvs[12]
    annotation_data_dir = argvs[13]
    if assays == 'all_assays':
        assays = util.find_assay_list(epigenome_list=epigenomes, \
            data_dir=training_data_dir, for_concatenated_mode=True)
    else:
        assays = sorted(assays[1:len(assays)-1].split(','))
    E = len(assays)
    script_dir = os.path.split(os.path.realpath(__file__))[0]
    root_path = os.path.join(script_dir, "../../../")
    configs = myconfig.config(root_path, log_dir=script_dir, data_type=data_type)
    if chromosomes == 'all_chromosomes':
        chromosomes = configs.chromosome_list
    else:
        chromosomes = chromosomes[1:len(chromosomes)-1].split(',')
    # print
    print('')
    print('task: ', task)
    print('K: ', K)
    print('E: ', E)
    print('lambda_1_L2: ', LAMBDA_1_L2)
    print('lambda_2_L1: ', LAMBDA_2_L1)
    print('lambda_2_L2: ', LAMBDA_2_L2)
    print('lambda_3_L2: ', LAMBDA_3_L2)
    print('training_itr: ', TRAINING_ITR)
    print('stop_threshold: ', STOP_THRESHOLD)
    print('window_size: ', WINDOW_SIZE)
    print('nonneg_state: ', nonneg_state)
    print('sumone_state: ', sumone_state)
    print('nonneg_em: ', nonneg_em)
    print('message_passing: ', message_passing)
    print('epigenomes: ', epigenomes)
    print('assays: ', assays)
    print('chromosomes: ', chromosomes)
    print('resolution: ', resolution)
    print('data_type: ', data_type)
    print('training_data_dir: ', training_data_dir)
    print('annotation_data_dir: ', annotation_data_dir)
    print('')
    # save to log
    log_f_name = 'log-k'
    for k in K:
        log_f_name += '-' + str(k)
    configs.openLog(log_f_name)
    configs.toLog('task: {}'.format(task))
    configs.toLog('K: {}'.format(K))
    configs.toLog('E: {}'.format(E))
    configs.toLog('lambda_1_L2: {}'.format(LAMBDA_1_L2))
    configs.toLog('lambda_2_L1: {}'.format(LAMBDA_2_L1))
    configs.toLog('lambda_2_L2: {}'.format(LAMBDA_2_L2))
    configs.toLog('lambda_3_L2: {}'.format(LAMBDA_3_L2))
    configs.toLog('training_itr: {}'.format(TRAINING_ITR))
    configs.toLog('stop_threshold: {}'.format(STOP_THRESHOLD))
    configs.toLog('window_size: {}'.format(WINDOW_SIZE))
    configs.toLog('nonneg_state: {}'.format(nonneg_state))
    configs.toLog('sumone_state: {}'.format(sumone_state))
    configs.toLog('nonneg_em: {}'.format(nonneg_em))
    configs.toLog('message_passing: {}'.format(message_passing))
    configs.toLog('epigenomes: {}'.format(epigenomes))
    configs.toLog('assays: {}'.format(assays))
    configs.toLog('chromosomes: {}'.format(chromosomes))
    configs.toLog('resolution: {}'.format(resolution))
    configs.toLog('data_type: {}'.format(data_type))
    configs.toLog('training_data_dir: {}'.format(training_data_dir))
    configs.toLog('annotation_data_dir: {}'.format(annotation_data_dir))
    # system info
    configs.toLog('\n## system info')
    configs.toLog('platform: ' + str(platform.system()))
    configs.toLog('platform-release: ' + str(platform.release()))
    configs.toLog('platform-version: ' + str(platform.version()))
    configs.toLog('CPU: ' + str(platform.processor()))
    configs.toLog('Available CPU cores: ' + str(psutil.cpu_count(logical=False)))
    configs.toLog('RAM: ' + str(round(psutil.virtual_memory().total/(1024.0 **3))) + ' GB')
    # create necessary directories
    if not os.path.exists(TEMP_DIR):
        os.mkdir(TEMP_DIR)
    if not os.path.exists(OUTPUT_DIR):
        os.mkdir(OUTPUT_DIR)
    for k in K:
        temp_dir_k = os.path.join(TEMP_DIR, 'k' + str(k))
        if not os.path.exists(temp_dir_k):
            os.mkdir(temp_dir_k)
        annotation_result_dir = os.path.join(OUTPUT_DIR, 'annotation-result-k' + str(k))
        if not os.path.exists(annotation_result_dir):
            os.mkdir(annotation_result_dir)
    configs.closeLog()


def run_training():
    global configs
    global K
    global E
    global nonneg_state
    global sumone_state
    global nonneg_em
    global message_passing
    global epigenomes
    global assays
    global chromosomes
    global resolution
    global training_data_dir
    # get the training regions
    index_region = {}
    with open(configs.pilot_index_region(resolution), 'r') as index_region_f: 
        while True:
            line = index_region_f.readline()
            if not line:
                break
            inf = line.strip('\n').split()
            chromosome = inf[0]
            start = int(inf[1])
            end = int(inf[2])
            if inf[0] not in index_region:
                index_region[chromosome] = [(start, end)]
            else:
                index_region[chromosome].append((start, end))
    # run training for each k
    error_m = None # model error
    for k in K:
        temp_dir_k = os.path.join(TEMP_DIR, 'k' + str(k))
        # read input data
        data_array = [[] for i in range(E)]
        g = 0 # total length of training regions
        for epigenome in epigenomes:
            for chromosome in index_region:
                if chromosome not in chromosomes:
                    continue
                skip = 0
                for region in index_region[chromosome]:
                    start_idx = region[0]
                    end_idx = region[1]
                    region_data_array, region_length = get_data_array(\
                        epigenome, chromosome, start_idx, end_idx, skip, temp_dir_k)
                    # concatenate
                    for i in range(E):
                        data_array[i] += region_data_array[i]
                    g += region_length
                    skip += end_idx - start_idx + 1
        model_file_path = os.path.join(OUTPUT_DIR, 'model_k' + str(k) + '.json')
        # check if the models exist already
        if os.path.exists(model_file_path):
            print('\n### k=' + str(k) + ' model exists already!')
            continue
        else:
            print('\n### k=' + str(k) + ' model')
        # create a new model
        model = ssm.ssm(E=E, G=g, K=k, \
            lambda_1_l2=LAMBDA_1_L2, lambda_2_l1=LAMBDA_2_L1, lambda_2_l2=LAMBDA_2_L2, lambda_3_l2=LAMBDA_3_L2, \
            positive_state=nonneg_state, sumone_state=sumone_state, positive_em=nonneg_em, message_passing=message_passing, \
            verbose=False)
        model.set_x(np.asmatrix(data_array))
        model.optimization(iteration=TRAINING_ITR)
        error_m = model.error_m
        # save the model
        model_param = {
            'E': E,
            'nonneg_state': nonneg_state,
            'sumone_state': sumone_state,
            'nonneg_em': nonneg_em,
            'message_passing': message_passing,
            'lambda_1_l2': LAMBDA_1_L2,
            'lambda_2_l1': LAMBDA_2_L1,
            'lambda_2_l2': LAMBDA_2_L2,
            'lambda_3_l2': LAMBDA_3_L2,
            'theta_m': model.theta_m.tolist(),
            'lambda_m': model.lambda_m.tolist()
        }
        with open(model_file_path, 'w') as out_file:
            json.dump(model_param, out_file)
        # save the errors
        np.array(error_m).tofile(os.path.join(OUTPUT_DIR, 'error_m_k' + str(k) + '.txt'))


def get_data_array(epigenome, chromosome, start_idx, end_idx, skip, temp_dir_k):
    global resolution
    global assays
    global training_data_dir
    # prepare data files
    resolution_files = []
    skipped_files = []
    for assay in assays:
        resolution_file = os.path.join(training_data_dir, epigenome + "-" + assay + \
            "/resolution-" + str(resolution) + "bp/", chromosome + ".bedGraph")
        resolution_files.append(resolution_file)
        skipped_file = os.path.join(temp_dir_k, \
            epigenome + "-" + assay + chromosome + str(skip) + ".txt")
        skipped_files.append(skipped_file)
        subprocess.call(" ".join(["tail", "-n", \
            "+" + str(skip+1), resolution_file, ">", skipped_file]), shell=True)
    # read data
    skipped_fs = []
    for skipped_file in skipped_files:
        skipped_fs.append(open(skipped_file, 'r'))
    is_first_index = True
    g = 0 # region length
    data_array = [[] for i in range(E)]
    while True:
        for i in range(len(skipped_fs)):
            skipped_f = skipped_fs[i]
            line = skipped_f.readline().strip('\n')
            if not line:
                raise Exception("invalid line")
            inf = line.split()
            index = int(inf[0])
            signal = float(inf[1])
            if is_first_index:
                if index != start_idx:
                    raise Exception("start index not begin with 0")
                is_first_index = False
            data_array[i].append(np.arcsinh(signal)) # data transformation
        g += 1
        if index == end_idx:
            break
    # clean up the skipped_files
    for skipped_f in skipped_fs:
        skipped_f.close()
    for skipped_file in skipped_files:
        os.remove(skipped_file)
    return data_array, g


def run_training_in_parallel():
    global configs
    global K
    global E
    global nonneg_state
    global sumone_state
    global nonneg_em
    global message_passing
    global epigenomes
    global assays
    global chromosomes
    global resolution
    global training_data_dir
    # read data
    index_region = {}
    with open(configs.pilot_index_region(resolution), 'r') as index_region_f: 
        while True:
            line = index_region_f.readline()
            if not line:
                break
            inf = line.strip('\n').split()
            chromosome = inf[0]
            start = int(inf[1])
            end = int(inf[2])
            if inf[0] not in index_region:
                index_region[chromosome] = [(start, end)]
            else:
                index_region[chromosome].append((start, end))
    # run training for each k
    error_m = None # model error
    for k in K:
        temp_dir_k = os.path.join(TEMP_DIR, 'k' + str(k))
        model_file_path = os.path.join(OUTPUT_DIR, 'model_k' + str(k) + '.json')
        # check if the models exist already
        if os.path.exists(model_file_path):
            print('\n### k=' + str(k) + ' model exists already!')
            continue
        else:
            print('\n### k=' + str(k) + ' model')
        # create a new model
        model = ssm.ssm(E=E, G=WINDOW_SIZE, K=k, \
            lambda_1_l2=LAMBDA_1_L2, lambda_2_l1=LAMBDA_2_L1, lambda_2_l2=LAMBDA_2_L2, lambda_3_l2=LAMBDA_3_L2, \
            positive_state=nonneg_state, sumone_state=sumone_state, positive_em=nonneg_em, message_passing=message_passing, \
            verbose=False)
        # save the model
        model_param = {
            'E': E,
            'nonneg_state': nonneg_state,
            'sumone_state': sumone_state,
            'nonneg_em': nonneg_em,
            'message_passing': message_passing,
            'lambda_1_l2': LAMBDA_1_L2,
            'lambda_2_l1': LAMBDA_2_L1,
            'lambda_2_l2': LAMBDA_2_L2,
            'lambda_3_l2': LAMBDA_3_L2,
            'theta_m': model.theta_m.tolist(),
            'lambda_m': model.lambda_m.tolist()
        }
        with open(model_file_path, 'w') as out_file:
            json.dump(model_param, out_file)
        error_m = []
        for itr in range(TRAINING_ITR):
            #print('\n# itr ', itr+1)
            pid = 0 # process id
            p_list = [] # list of processes
            for epigenome in epigenomes:
                for chromosome in index_region:
                    if chromosome not in chromosomes:
                        continue
                    skip = 0
                    for region in index_region[chromosome]:
                        start_idx = region[0]
                        end_idx = region[1]
                        try:
                            p = subprocess.Popen(["python", "training_node.py", epigenome, chromosome, \
                                str(start_idx), str(end_idx), str(skip), str(k), str(resolution), \
                                ','.join(assays), training_data_dir, temp_dir_k, OUTPUT_DIR, str(pid)])
                        except Exception as e:
                            print(e)
                        p_list.append(p)
                        pid += 1
                        skip += end_idx - start_idx + 1
            # check if all the processes are done
            while True:
                all_done = True
                for p in p_list:
                    poll = p.poll()
                    if poll == None:
                        all_done = False
                if all_done:
                    break
            # take the sufficient statistics to update the model parameters
            theta_lhs = np.asmatrix(np.zeros((k, k)))
            theta_rhs = np.asmatrix(np.zeros((k, E)))
            lambda_lhs = np.asmatrix(np.zeros((k, k)))
            lambda_rhs = np.asmatrix(np.zeros((k, k)))
            err_m = .0
            g = .0
            for i in range(pid):
                file_path = os.path.join(temp_dir_k, str(i) + '.json')
                obj_text = codecs.open(file_path, 'r', \
                    encoding='utf-8').read()
                node_stats = json.loads(obj_text)
                theta_lhs += node_stats['theta_lhs']
                theta_rhs += node_stats['theta_rhs']
                lambda_lhs += node_stats['lambda_lhs']
                lambda_rhs += node_stats['lambda_rhs']
                err_m += node_stats['error_m']
                g += node_stats['g']
                os.remove(file_path)
            # normalize by total length of all regions
            theta_lhs /= g
            theta_rhs /= g
            lambda_lhs /= g
            lambda_rhs /= g
            err_m /= g
            # update the model parameters
            model.update_theta(theta_lhs, theta_rhs)
            model.update_lambda(lambda_lhs, lambda_rhs)
            # save the model
            model_param = {
                'E': E,
                'nonneg_state': nonneg_state,
                'sumone_state': sumone_state,
                'nonneg_em': nonneg_em,
                'message_passing': message_passing,
                'lambda_1_l2': LAMBDA_1_L2,
                'lambda_2_l1': LAMBDA_2_L1,
                'lambda_2_l2': LAMBDA_2_L2,
                'lambda_3_l2': LAMBDA_3_L2,
                'theta_m': model.theta_m.tolist(),
                'lambda_m': model.lambda_m.tolist()
            }
            with open(model_file_path, 'w') as out_file:
                json.dump(model_param, out_file)
            error_m.append(err_m)
            #print('## error: ', err_m)
        # save the errors
        np.array(error_m).tofile(os.path.join(OUTPUT_DIR, 'error_m_k' + str(k) + '.txt'))


def run_annotation():
    global util
    global configs
    global K
    global E
    global nonneg_state
    global nonneg_em
    global epigenomes
    global assays
    global chromosomes
    global resolution
    global training_data_dir
    global annotation_data_dir
    # read and prepare data
    index_region = {}
    with open(configs.assay_region_file, 'r') as assay_region_f:
        while True:
            line = assay_region_f.readline()
            if not line:
                break
            inf = line.strip('\n').split()
            chromosome = inf[0]
            start = int(inf[1])
            end = int(inf[2])
            start_index = util.resolution_index_locate(start, resolution)
            end_index = util.resolution_index_locate(end, resolution)
            if chromosome in index_region:
                raise(Exception("Multi chromosome region"))
            index_region[chromosome] = [(start_index, end_index)]
    # remove blacklist_index_region
    blacklist_index_region = {}
    with open(configs.blacklist_region_file, 'r') as blacklist_region_f:
        while True:
            line = blacklist_region_f.readline()
            if not line:
                break
            inf = line.strip('\n').split()
            chromosome = inf[0]
            start = int(inf[1])
            end = int(inf[2])
            start_index = util.resolution_index_locate(start, resolution)
            end_index = util.resolution_index_locate(end, resolution)
            if chromosome not in blacklist_index_region:
                blacklist_index_region[chromosome] = [(start_index, end_index)]
            else:
                blacklist_index_region[chromosome].append((start_index, end_index))
    # apply exclude intersection: assay_index_region - blacklist_index_region
    for chromosome in blacklist_index_region:
        for region_blacklist in blacklist_index_region[chromosome]:
            for temp in index_region[chromosome]:
                region_assay = index_region[chromosome][0]
                index_region[chromosome].remove(region_assay)
                index_region[chromosome] += util.rm_intersect(region_assay, region_blacklist, exclude=True)
    for chromosome in index_region:
        index_region[chromosome] = list(sorted(index_region[chromosome], key=lambda x:x[0]))
    # run annotation for each k
    assemble = False # whether or not to assemble annotation files at the end
    for k in K:
        print('\nk=' + str(k) + '\n')
        temp_dir_k = os.path.join(TEMP_DIR, 'k' + str(k))
        annotation_result_dir = os.path.join(OUTPUT_DIR, 'annotation-result-k' + str(k))
        # assign jobs to subprocesses
        pid = 0 # process id
        p_list = [] # list of processes
        for epigenome in epigenomes:
            # check if the annotation for this chromosome exists already
            epigenome_annotation_dir = os.path.join(annotation_result_dir, epigenome)
            if os.path.exists(epigenome_annotation_dir):
                print('Annotation for ' + epigenome + ' exists already!')
                continue
            for chromosome in index_region:
                if chromosome not in chromosomes:
                    continue
                for region in index_region[chromosome]:
                    start_idx = region[0]
                    end_idx = region[1]
                    assemble = True
                    try:
                        p = subprocess.Popen(["python", "annotation_node.py", epigenome, chromosome, \
                            str(start_idx), str(end_idx), str(k), str(resolution), \
                            ','.join(assays), annotation_data_dir, temp_dir_k, OUTPUT_DIR, \
                            annotation_result_dir, str(pid)])
                    except Exception as e:
                        print(e)
                    p_list.append(p)
                    pid += 1
        # check if all the processes are done
        while True:
            all_done = True
            for p in p_list:
                poll = p.poll()
                if poll == None:
                    all_done = False
            if all_done:
                break
        if assemble:
            # assemble annotation files generated by nodes
            sub_annotation = {}
            files = []
            for (_, _, filenames) in os.walk(annotation_result_dir):
                files.extend(filenames)
                break
            for myfile in files:
                file_id = int(re.search(r'_\d+', myfile).group()[1:])
                chromosome = re.search(r'.*_', myfile).group()[:-1]
                if chromosome not in sub_annotation:
                    sub_annotation[chromosome] = [(myfile, file_id)]
                else:
                    sub_annotation[chromosome].append((myfile, file_id))
            for chromosome in sub_annotation:
                sub_annotation[chromosome] = list(sorted(sub_annotation[chromosome], key=lambda x:x[1]))
                final_annotation_file = os.path.join(annotation_result_dir, chromosome + ".bed")
                open(final_annotation_file, 'w')
                for sub_annotation_file,_ in sub_annotation[chromosome]:
                    result = subprocess.call(" ".join(["cat", os.path.join(annotation_result_dir, sub_annotation_file), \
                        ">>", final_annotation_file]), shell=True)
                    os.remove(os.path.join(annotation_result_dir, sub_annotation_file))
            # move files into their epigenome's directory
            files = []
            for (_, _, filenames) in os.walk(annotation_result_dir):
                files.extend(filenames)
                break
            for file in files:
                # remove epigenome identifier from the file name
                epigenome_id = file[:file.find('chr')]
                epigenome_dir = os.path.join(annotation_result_dir, epigenome_id)
                if not os.path.exists(epigenome_dir):
                    os.makedirs(epigenome_dir)
                source_path = os.path.join(annotation_result_dir, file)
                target_path = os.path.join(epigenome_dir, file[file.find('chr'):])
                os.rename(source_path, target_path)


def run_evaluation():
    global K
    # run evaluation for each k
    for k in K:
        print('\nk=' + str(k) + '\n')
        eval_whole_gene_bin(k)
        eval_enhancer_bin(k)


def eval_whole_gene_bin(k):
    print('evaluation: whole gene bin..')
    element = 'genebody_bin' # the element to get the features from
    k_bin = k * NUM_EVAL_BINS
    # check if the feature file for the element already exists
    feature_file = os.path.join(OUTPUT_DIR, element+'_features-'+'k'+str(k)+'.csv')
    if not os.path.exists(feature_file):
        get_features(k, element)
    element_df = pd.read_csv(feature_file)
    element_df = element_df.dropna() # drop rows with Nan values
    element_df['expression value'] = np.arcsinh(element_df['expression value']) # transformation
    # prediction
    train_scores = []
    test_scores = []
    for t in range(NUM_EVAL_TRIALS):
        X_train, X_test, y_train, y_test = train_test_split(\
            element_df.iloc[:,-k_bin:], element_df['expression value'], test_size=0.2)
        # train
        num_points = X_train.shape[0]
        reg = LinearRegression(n_jobs=-1).fit(X_train, y_train)
        y_pred = reg.predict(X_train)
        adj_r2 = util.adj_r2_score(y_train, y_pred, num_points, k_bin)
        train_scores.append(adj_r2)
        # test
        num_points = X_test.shape[0]
        y_pred = reg.predict(X_test)
        adj_r2 = util.adj_r2_score(y_test, y_pred, num_points, k_bin)
        test_scores.append(adj_r2)
    train_score = np.mean(train_scores)
    test_score = np.mean(test_scores)
    print('train score: ' + str(train_score))
    print('test score: ' + str(test_score))
    print('')
    # save
    with open(os.path.join(OUTPUT_DIR, 'eval_' + element + '-' + 'k'+str(k) + '.txt'), 'w') \
        as out_file:
        print(str(NUM_EVAL_TRIALS) + ' runs', file=out_file)
        for score in test_scores:
            print(score, file=out_file)
        print('average:\n{}'.format(test_score), file=out_file)


def eval_enhancer_bin(k):
    print('evaluation: enhancer bin..')
    element = 'enhancer_bin' # the element to get the features from
    k_bin = k * NUM_EVAL_BINS
    # check if the feature file for the element already exists
    feature_file = os.path.join(OUTPUT_DIR, element + '_features-' + 'k'+str(k) + '.csv')
    if not os.path.exists(feature_file):
        get_features(k, element)
    element_df = pd.read_csv(feature_file)
    element_df = element_df.dropna() # drop rows with Nan values
    element_df['activity value'] = np.arcsinh(element_df['activity value']) # transformation
    # prediction
    train_scores = []
    test_scores = []
    for t in range(NUM_EVAL_TRIALS):
        X_train, X_test, y_train, y_test = train_test_split(\
            element_df.iloc[:,-k_bin:], element_df['activity value'], test_size=0.2)
        # train
        num_points = X_train.shape[0]
        reg = LinearRegression(n_jobs=-1).fit(X_train, y_train)
        y_pred = reg.predict(X_train)
        adj_r2 = util.adj_r2_score(y_train, y_pred, num_points, k_bin)
        train_scores.append(adj_r2)
        # test
        num_points = X_test.shape[0]
        y_pred = reg.predict(X_test)
        adj_r2 = util.adj_r2_score(y_test, y_pred, num_points, k_bin)
        test_scores.append(adj_r2)
    train_score = np.mean(train_scores)
    test_score = np.mean(test_scores)
    print('train score: ' + str(train_score))
    print('test score: ' + str(test_score))
    print('')
    # save
    with open(os.path.join(OUTPUT_DIR, 'eval_' + element + '-' + 'k'+str(k) + '.txt'), 'w') \
        as out_file:
        print(str(NUM_EVAL_TRIALS) + ' runs', file=out_file)
        for score in test_scores:
            print(score, file=out_file)
        print('average:\n{}'.format(test_score), file=out_file)


def get_features(k, element):
    global configs
    global epigenomes
    global chromosomes
    print('getting ' + element + ' features..')
    temp_dir_k = os.path.join(TEMP_DIR, 'k' + str(k))
    for epigenome in epigenomes:
        if epigenome == 'E017':
            continue # because gene exp. data is not available for E017
        print(epigenome)
        child_process_list = []
        region_file = ''
        if element == 'genebody_bin':
            region_file = os.path.join(configs.data_path, 'tss-data', epigenome + '-filtered_gexp1kb.bed')
        if element == 'enhancer_bin':
            region_file = os.path.join(configs.data_path, 'enhancer-data', epigenome + '-filtered_enhancer.bed')
        annotation_result_dir = os.path.join(OUTPUT_DIR, 'annotation-result-k' + str(k))
        for chromosome in chromosomes:
            success = False
            while not success:
                try:
                    node_script = element + '_features_node.py'
                    child_process = subprocess.Popen(['python', node_script, str(k), epigenome, chromosome, \
                        str(NUM_EVAL_BINS), region_file, temp_dir_k, annotation_result_dir], \
                        stdout=subprocess.PIPE)
                except Exception as e:
                    print(e)
                    time.sleep(60)
                    continue # rebatch after some sleep
                child_process_list.append(child_process)
                success = True
            # sleep after each job submission to avoid scoket time out error
            time.sleep(1)
        while True:
            all_done = True
            for child_process in child_process_list:
                poll = child_process.poll()
                if poll == None:
                    all_done = False
            if all_done:
                break
    print('combining the output files..')
    final_df = pd.DataFrame()
    for epigenome in epigenomes:
        if epigenome == 'E017':
            continue # because gene exp. data is not available for E017
        for chromosome in chromosomes:
            chr_df_file = os.path.join(temp_dir_k, \
                epigenome+'-'+chromosome+'_'+element+'_features-' + 'k'+str(k) + '.csv')
            chr_df = pd.read_csv(chr_df_file)
            final_df = pd.concat([final_df, chr_df], ignore_index=True)
            os.remove(chr_df_file)
    final_df.to_csv(os.path.join(OUTPUT_DIR, element+'_features-' + 'k'+str(k) + '.csv'), index=False)


def viz(argvs):
    global K
    K = [int(i) for i in argvs[2][1:len(argvs[2])-1].split(',')]
    plot_training_err()

def plot_training_err():
    global K
    # plot the training error
    err_plot = plt.figure().gca()
    for k in K:
        err_k = np.fromfile(os.path.join(OUTPUT_DIR, 'error_m_k' + str(k) + '.txt'))
        print('k'+str(k) + ' error: ', err_k)
        err_k = err_k[1:] # to get rid of the initial error which might be too high
        err_plot.plot(range(1, len(err_k)+1), err_k, label='k='+str(k))
    plt.ylabel('Negative log-likelihood')
    plt.xlabel('Iteration')
    err_plot.xaxis.set_major_locator(MaxNLocator(integer=True))
    err_plot.legend(loc='best')
    plt.savefig('error_plot.pdf')
    print('')

if __name__ == "__main__":
    sys.exit(main(sys.argv))
