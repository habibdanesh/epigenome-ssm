import os
import math
import numpy as np
import sklearn
from sklearn.metrics import r2_score

def intersect(a: tuple, b: tuple):
    if a[0] < b[0]:
        if a[1] <= b[0]:
            return None
        elif a[1] > b[1]:
            return b
        elif a[1] <= b[1]:
            return b[0], a[1]
    else:
        if a[0] >= b[1]:
            return None
        elif a[1] > b[1]:
            return a[0], b[1]
        elif a[1] <= b[1]:
            return a

def test_interset():
    print(intersect((1, 4), (0, 2)),
          intersect((1, 4), (0, 0.5)),
          intersect((1, 4), (2, 3)),
          intersect((1, 4), (3, 5)),
          intersect((1, 4), (5, 6)),
          intersect((1, 4), (0, 6)))

def rm_intersect(a:tuple, b:tuple, exclude=False):
    '''
    :param a: region A
    :param b: region B
    :return: an array of region from A - A intersect B
    '''
    result = []
    inter = intersect(a, b)
    if inter is not None:
        if inter[0] != a[0]:
            if exclude:
                if inter[0]-1 >= a[0]:
                    result.append((a[0], inter[0]-1))
            else:
                result.append((a[0], inter[0]))
        if inter[1] != a[1]:
            if exclude:
                if inter[1]+1 <= a[1]:
                    result.append((inter[1]+1, a[1]))
            else:
                result.append((inter[1], a[1]))
    else:
        result.append(a)
    return result

def test_rm_interset():
    print(rm_intersect((1, 4), (0, 2)),
          rm_intersect((1, 4), (0, 0.5)),
          rm_intersect((1, 4), (1, 4)),
          rm_intersect((1, 4), (2, 3)),
          rm_intersect((1, 4), (3, 5)),
          rm_intersect((1, 4), (5, 6)),
          rm_intersect((1, 4), (0, 6)))
    print(rm_intersect((1, 4), (0, 2), exclude=True),
          rm_intersect((1, 4), (0, 0.5), exclude=True),
          rm_intersect((1, 4), (1, 4), exclude=True),
          rm_intersect((1, 4), (2, 3), exclude=True),
          rm_intersect((1, 4), (3, 5), exclude=True),
          rm_intersect((1, 4), (5, 6), exclude=True),
          rm_intersect((1, 4), (0, 6), exclude=True))
'''
    probability  
'''


def log_add(p1, p2, base=math.e):
    return p1 + math.log(1 + math.exp(p2 - p1), base)


def log_multi(p1, p2, base=math.e):
    return math.log(p1, base) + math.log(p2, base)


def test_on_log_pro():
    p1 = math.exp(-30)
    p2 = 2 * math.exp(-30)
    pass


def create_resolution(in_file, out_file, base, final, resolution=100):
    '''
    make the chrs into resolution(window), default as 100bp
    base,final for start and end of frame
    :param in_file: input file
    :param out_file: output resolution file
    :param base: start position
    :param final: end position
    :param resolution: size of resolution
    :return:
    '''
    in_f = open(in_file, 'r')
    out_f = open(out_file, 'a')
    index = 0
    remain = 0
    cut = resolution
    while True:
        line = in_f.readline()
        if not line:
            break
        inf = line.strip('\n').split()
        start = int(inf[0])
        end = int(inf[1])

        # within the frame
        if end <= base:
            continue
        elif start >= final:
            break  # finish the resolution on this frame
        elif start < base and end > base:
            start = base
        elif start < final and end > final:
            end = final

        signal = float(inf[2])
        if cut >= end - start:
            remain += (end - start) * signal
            cut = (index + 1) * resolution + base - end
        else:
            while (index + 1) * resolution + base < end:  # careful here
                avg_signal = (remain + cut * signal) / resolution
                # print(index * resolution + base, base + (index + 1) * resolution, avg_signal, file=out_f)
                print(index, avg_signal, file=out_f)
                remain = 0
                cut = resolution
                index += 1
            remain = (end - index * resolution - base) * signal
            cut = (index + 1) * resolution + base - end
    # print(index * resolution + base, (index + 1) * resolution +base, remain / resolution, file=out_f)
    print(index, remain / resolution, file=out_f)
    out_f.close()
    in_f.close()

def create_resolution_plus(in_file, out_file, base, final, resolution=100, out_mod='w'):
    '''
    another version to create resolution
    allow to create resolution direct from raw data by specify chromosome
    :param in_file: input file
    :param out_file: output resolution file
    :param chr_resol: chromosome name
    :param base: start position
    :param final: end position
    :param resolution: size of resolution
    :return:
    '''
    in_f = open(in_file, 'r')
    out_f = open(out_file, out_mod)
    index = 0
    remain = 0
    cut = resolution
    while True:
        line = in_f.readline()
        if not line:
            break
        inf = line.strip('\n').split()
        start = int(inf[0])
        end = int(inf[1])
        signal = float(inf[2])
 
        same_chr = True
        # within the frame
        if end <= base:
            continue
        elif start >= final:
            break  # finish the resolution on this frame
        elif start < base and end > base:
            start = base
        elif start < final and end > final:
            end = final

        if cut >= end - start:
            remain += (end - start) * signal
            cut = (index + 1) * resolution + base - end
        else:
            while (index + 1) * resolution + base < end:  # careful here
                avg_signal = (remain + cut * signal) / resolution
                # print(index * resolution + base, base + (index + 1) * resolution, avg_signal, file=out_f)
                print(index, avg_signal, file=out_f)
                remain = 0
                cut = resolution
                index += 1
            remain = (end - index * resolution - base) * signal
            cut = (index + 1) * resolution + base - end
    # print(index * resolution + base, (index + 1) * resolution +base, remain / resolution, file=out_f)
    print(index, remain / resolution, file=out_f)
    out_f.close()
    in_f.close()


def output_matrix(matrix, filename):
    '''
    output a matrix data to a file
    :param filename:
    :return:
    '''
    with open(filename, 'w') as f:
        print(matrix_to_str(matrix), file=f)

def read_matrix(filename):
    '''
    read in a format
    :param filename:
    :return: numpy matrix
    '''
    in_f = open(filename, "r")
    inf = in_f.readline().strip('\n').split()
    m = int(inf[0])
    n = int(inf[1])
    matrix = np.asmatrix(np.zeros((m,n)))
    for i in range(m):
        inf = in_f.readline().strip('\n').split()
        for j in range(n):
            matrix[i,j] = float(inf[j])
    in_f.close()
    return matrix

def read_threshold_matrix(in_f, m, n):
    '''
    read in a format
    :param input file:
    :return: numpy matrix
    '''
    matrix = np.asmatrix(np.zeros((m,n)))
    for i in range(m):
        inf = in_f.readline().strip('\n').split()
        for j in range(n):
            matrix[i,j] = float(inf[j])
    in_f.close()
    return matrix

def matrix_to_str(m):
    """
    :param m: numpy matrix
    :return s: nicely formatted matrix data in string
    """
    s = ""
    row,col = m.shape
    s += ("{} {}\n".format(row, col))
    for i in range(row):
        for j in range(col):
            s += "{} ".format(m[i,j])
        s += "\n"
    return s

def test_matrix_readinout():
    a = np.asmatrix([[1, 2, 3], [3, 3, 4]])
    output_matrix(a, "out.txt")
    print(readin_matrix("out.txt"))

def read_trained_ssm(in_file):
    '''
    read in theta param for ssm
    :param in_file: input theta file
    :return positive_flag, penalize_flag, factor_count, lamda1, assay_count, theta_m:parameter of ssm model
    update 2018-July-01: currently abandoned as the parameter changes in trainning
    '''
    in_f = open(in_file)
    positive_flag = bool(in_f.readline().strip('\n').split()[1])
    penalize_flag = bool(in_f.readline().strip('\n').split()[1])
    lamda1 = float(in_f.readline().strip('\n').split()[1])
    factor_count, assay_count = [int(element) for element in in_f.readline().strip('\n').split()]
    theta_m = np.asmatrix(np.zeros((factor_count, assay_count)))
    for k in range(factor_count):
        inf = in_f.readline().strip('\n').split()
        for e in range(assay_count):
            theta_m[k, e] = inf[e]
    return positive_flag, penalize_flag, factor_count, lamda1, assay_count, theta_m

def resolution_index_locate(pos, resolution_size, start=0):
    """
    :param pos: bp position 
    :return index: index after resolution applied
    """
    return math.floor((pos - start)*1.0 / resolution_size)

def adj_r2_score(real_data, predict_data, n, p):
    """
    :param n: number of samples
    :param p: number of states
    """
    r2 = r2_score(real_data, predict_data)
    adj_r2 = 1 - (1-r2)*(n-1)/(n-p-1)
    return adj_r2

def find_assay_list(epigenome_list, data_dir, for_concatenated_mode=True):
    assay_list = []
    for (_, dirs, files) in os.walk(data_dir):
        for dir_name in dirs:
            if dir_name[:4] in epigenome_list:
                assay_list.append(dir_name)
    if for_concatenated_mode:
        #remove those assays which are not presented in ALL the epigenomes, so that all the epigneomes will have the same set of assays
        for assay in assay_list:
            presented_in_all_epigenomes = True
            assay = assay[4:] #remove the epigenome identifier
            for epigenome in epigenome_list:
                if epigenome+assay not in assay_list:
                    presented_in_all_epigenomes = False
                    break
            if not presented_in_all_epigenomes:
                assay_list = [a for a in assay_list if a[4:]!=assay]
        assay_list = [a[5:] for a in assay_list] #remove the epigenome identifier and '-'
        assay_list = sorted(list(set(assay_list))) #remove duplicates
    return assay_list

if __name__ == '__main__':
    test_rm_interset()

