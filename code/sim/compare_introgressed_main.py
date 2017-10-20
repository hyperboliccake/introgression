import sys
import process_args
import sim_process
from compare_introgressed import *
sys.path.append('..')
import global_params as gp

# takes in two sets of introgression calls on same set of data and
# compares them (for example, actual and predicted regions, or calls
# from two different prediction methods)

args, last_read = process_args.process_args(sys.argv)
suffix1 = sys.argv[last_read + 1]
suffix2 = sys.argv[last_read + 2]


gp_dir = '../'
fn1 = gp_dir + gp.sim_out_dir + gp.sim_out_prefix + args['tag'] + \
      '_introgressed_' + suffix1 + '.txt'
fn2 = gp_dir + gp.sim_out_dir + gp.sim_out_prefix + args['tag'] + \
      '_introgressed_' + suffix2 + '.txt'
f_out = gp_dir + gp.sim_out_dir + gp.sim_out_prefix + args['tag'] + \
      '_introgressed_compare_' + suffix1 + '_' + suffix2 + '.txt'

f1 = open(fn1, 'r')
f2 = open(fn2, 'r')
f_out = open(f_out, 'w')

line1 = f1.readline()
line2 = f2.readline()

write_compare_header(f_out, args['states'], suffix1, suffix2)
while line1 != '' and line2 != '':

    d1, rep1, line1 = sim_process.read_introgression_blocks(f1, line1, args['states'])
    d2, rep2, line2 = sim_process.read_introgression_blocks(f2, line2, args['states'])
    assert rep1 == rep2, str(rep1) + ' ' + str(rep2)
    print 'rep', rep1

    base_counts, avg_base_counts = count_bases(d1, d2, args)

    write_compare_line(avg_base_counts, f_out, args['states'], suffix1, suffix2)

if line1 != '' or line2 != '':
    print 'one of these files is incomplete and for some reason I\'m not bothering to tell you which!'

f_out.close()
