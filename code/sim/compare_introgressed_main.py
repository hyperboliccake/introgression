import sys
import sim_process
from compare_introgressed import *
sys.path.append('..')
import global_params as gp

# takes in two sets of introgression calls on same set of data and
# compares them (for example, actual and predicted regions, or calls
# from two different prediction methods)

tag = sys.argv[1]
suffix1 = sys.argv[2]
suffix2 = sys.argv[3]

gp_dir = '../'
fn1 = gp_dir + gp.sim_out_dir + gp.sim_out_prefix + tag + \
      '_introgressed_' + suffix1 + '.txt'
fn2 = gp_dir + gp.sim_out_dir + gp.sim_out_prefix + tag + \
      '_introgressed_' + suffix2 + '.txt'

f1 = open(fn1, 'r')
f2 = open(fn2, 'r')

line1 = f1.readline()
line2 = f2.readline()

while line1 != '':

    d1, rep1, line1 = sim_process.read_introgression(f1, line1)
    d2, rep2, line2 = sim_process.read_introgression(f2, line2)
    assert rep1 == rep2, str(rep1) + ' ' + str(rep2)
    print 'rep', rep1
    
    

assert line2 == '', line2
