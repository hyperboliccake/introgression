import sys
sys.path.append('..')
import global_params as gp
sys.path.append('../misc')
from mystats import *

def average_item(item):
    if item == '':
        return mean([])
    if item[0] == '[':
        l = []
        if item != '[]':            
            l = [float(x) for x in item[1:-1].split(',')]
        return mean(l)
    return float(item)

def read_summary_file(fn):
    f = open(fn, 'r')
    line = f.readline()
    labels = line[:-1].split('\t')
    ncols = len(labels)
    vals = [[] for l in labels]
    line = f.readline()
    while line != '':
        items = line[:-1].split('\t')
        items = [average_item(x) for x in items]
        for i in range(ncols):
            vals[i].append(items[i])
        line = f.readline()
    f.close()
    d_mean = dict(zip(labels, [mean(x) for x in vals]))
    d_std_err = dict(zip(labels, [std_err(x) for x in vals]))
    d_bootstrap = dict(zip(labels, [bootstrap(x) for x in vals]))
    return d_mean, d_std_err, d_bootstrap

def write_header(f, keys):
    
    f.write('line_type\ttag')
    for key in keys:
        f.write('\t' + key)
    f.write('\n')

def write_line_set(d_mean, d_std_err, d_bootstrap, keys, tag, f):

    # line type: mean
    f.write('mean\t' + tag)
    for key in keys:
        f.write('\t' + str(d_mean[key]))
    f.write('\n')

    # line type: std err
    f.write('std_err\t' + tag)
    for key in keys:
        f.write('\t' + str(d_std_err[key]))
    f.write('\n')

    # line type: bootstrap lower
    f.write('bs_lower\t' + tag)
    for key in keys:
        f.write('\t' + str(d_bootstrap[key][0]))
    f.write('\n')

    # line type: bootstrap upper
    f.write('bs_upper\t' + tag)
    for key in keys:
        f.write('\t' + str(d_bootstrap[key][1]))
    f.write('\n')

def aggregate_summary_files(fns, fn_all, tags):
    f_all = open(fn_all, 'w')
    header = True
    keys = None # for keeping order consistent
    for i in range(len(fns)):
        d_mean, d_std_err, d_bootstrap = read_summary_file(fns[i])
        if header:
            keys = d_mean.keys()
            write_header(f_all, keys)
            header = False
        write_line_set(d_mean, d_std_err, d_bootstrap, keys, tags[i], f_all)
    f_all.close()













