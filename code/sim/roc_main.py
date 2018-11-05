# somewhat hardcoded solution to making a roc curve to compare my method
# and phylonet-hmm

import sys
import os
import process_args
import sim_process
import sim_predict
import sim_predict_phylohmm
import roc
import gzip
sys.path.append('..')
import global_params as gp
sys.path.append('../misc')
import mystats

##======
# read in simulation and prediction parameters
##======

method = sys.argv[1]
sim_tag = sys.argv[3]
sim_args = process_args.process_args_by_tag(sys.argv[2], sim_tag)
predict_args = None
if method == 'predicted':
    predict_args, last_read = sim_predict.process_args(sys.argv, sim_args, i=3)
elif method == 'predicted_phylohmm':
    predict_args, last_read = sim_predict_phylohmm.process_args(sys.argv, sim_args, i=3)
else:
    print 'invalid method specified'
    sys.exit()

##======
# files 
##======

gp_dir = '../'
f_probs = open(gp.sim_out_dir_absolute + gp.sim_out_prefix + \
               sim_tag + '_introgressed_probs_' + method + '_' + \
               predict_args['predict_tag'] + '.txt','r')
actual_f = open(gp.sim_out_dir_absolute + gp.sim_out_prefix + sim_tag + \
                '_introgressed_actual.txt', 'r')
f = open(gp.sim_out_dir_absolute + gp.sim_out_prefix + sim_tag + \
         '_roc_' + method + '_' + predict_args['predict_tag'] + '.txt', 'w')

##======
# read in probabilties and positions
##======

print 'reading in posterior probabilities'
print 'rep 0/' + str(sim_args['num_reps']),

# list with one entry per rep
probs = []

line = f_probs.readline()
for rep in range(sim_args['num_reps']):
    print '\r' + 'rep ' + str(rep) + '/' + str(sim_args['num_reps']),
    sys.stdout.flush()
    # {1:{cer:.9,.9,..., par:.1,.1,...}}
    d, rep, line = sim_process.read_state_probs(f_probs, line)
    # {1:[{cer:.9, par:.1},{cer:.9, par:.1},...]}
    d = roc.reformat_probs(d)
    probs.append(d)
print
f_probs.close()


ps_fn = gp.sim_out_dir_absolute + gp.sim_out_prefix + sim_tag + \
        '_positions_' + method + '_' + predict_args['predict_tag'] + \
        '.txt.gz'
ps_f = gzip.open(ps_fn, 'rb')
ps = []
for rep in range(sim_args['num_reps']):
    x = [int(i) for i in ps_f.readline().strip().split('\t')[2:]]
    ps.append(x)
ps_f.close()

##======
# read in actual introgressed blocks
##======

print 'reading in actual introgressed blocks'

actual = []

line = actual_f.readline()
for rep in range(sim_args['num_reps']):
    d, rep, line = sim_process.read_introgression_blocks(actual_f, line, \
                                                         sim_args['species'])
    actual.append(d)

##======
# calculate statistics for different thresholds and write to file
##======

print 'generating ROC file'

#x = 20
#thresholds = [float(i) / x for i in range(0, x + 2)]
thresholds = [0, .00001, .00005, .0001, .0005, .001, .005, .01, .05, .1, .5, .6, .7, .8, .9, .99, 1, 1.1]
header = True
for threshold in thresholds:
    print 'threshold:', threshold
    all_stats = []
    for rep in range(sim_args['num_reps']):
        paths_t = {}
        for ind in probs[rep]:
            paths_t[ind] = roc.threshold_probs(probs[rep][ind], \
                                               threshold, \
                                               sim_args['species_to'], \
                                               ps[rep], 0, \
                                               sim_args['num_sites']-1)
        predicted = sim_process.convert_to_blocks(paths_t, predict_args['states'])
        stats = roc.get_stats(actual[rep], predicted, sim_args)
        all_stats.append(stats)
    avg_stats = {}
    for key in all_stats[0].keys():
        x = []
        for rep in range(sim_args['num_reps']):
            x.append(all_stats[rep][key])
        avg_stats[key] = mystats.mean(x)
    roc.write_roc_line(f, threshold, avg_stats, header)
    header = False
f.close()



