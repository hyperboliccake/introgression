import sys
import os
import gzip
import predict
sys.path.insert(0, '..')
import global_params as gp

args = predict.process_predict_args(sys.argv[1:])

header = open(gp.analysis_out_dir_absolute + args['tag'] + '/' + \
              'blocks_' + args['known_states'][0] + \
              '_' + args['tag'] + '_chr' + gp.chrms[0] + '_quality.txt', 'r').readline()

for species_from in args['known_states']:
    fn = gp.analysis_out_dir_absolute + args['tag'] + '/' + \
         'blocks_' + species_from + \
         '_' + args['tag'] + '_quality.txt'
    f = open(fn, 'w')
    f.write(header)
    for chrm in gp.chrms:
        fn_chrm = gp.analysis_out_dir_absolute + args['tag'] + '/' + \
                  'blocks_' + species_from + \
                  '_' + args['tag'] + '_chr' + chrm + '_quality.txt'
        try:
            fc = open(fn_chrm, 'r')
        except:
            continue
        fc.readline()
        for line in fc.readlines():
            f.write(line)
    f.close()

