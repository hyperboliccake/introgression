import os
import math
import numpy.random
import sys
sys.path.insert(0, '../misc/')
import mystats

def parse_list(l):
    assert len(l) >= 2, l
    if len(l) == 2:
        return []
    return [float(x) for x in l[1:-1].split(',')]

# power = true positives / all positives
# false positive rate = false positives / all negatives

params = [line.strip().split(' ') for line in open('sim_multi_model_args.txt', 'r').readlines()]
param_names = ['tag', 'model', 'N0', 'num_samples_par', 'num_samples_cer', 'par_cer_migration', 't_par_cer', 'num_sites', 'rho', 'outcross_rate', 'num_reps']
assert len(param_names) == len(params[0]), str(len(param_names)) + ' != ' + str(len(params[0]))
out_dir = '../../results/sim/run_3/'
prefix = 'sim_out_'
suffix = '_summary.txt'
ids = range(1, len(params) + 1)
f_all = open(out_dir + prefix + 'power_fpr' + suffix, 'w')
for i in ids:
    print i
    f = open(out_dir + prefix + str(i) + suffix, 'r')
    col_names = f.readline().strip().split('\t') # header
    assert col_names[0] == 'num_introgressed_cer'
    assert col_names[1] == 'num_introgressed_tracts_cer'
    assert col_names[2] == 'num_not_introgressed_tracts_cer'
    assert col_names[3] == 'num_predicted_introgressed_cer'
    assert col_names[4] == 'num_predicted_introgressed_tracts_cer'
    assert col_names[10] == 'num_introgressed_correct'    
    assert col_names[11] == 'num_predicted_tracts_actual'    
    assert col_names[12] == 'num_actual_tracts_predicted'    
    assert param_names[7] == 'num_sites'

    if i == ids[0]:
        for p in param_names:
            f_all.write(p + '\t')
        f_all.write('row_type\t')
        f_all.write('power\tfpr\tfdr\tpower_tracts\tfpr_tracts\tfdr_tracts\n')

    line = f.readline()
    power = []
    fpr = []
    fdr = []
    power_tracts = []
    fpr_tracts = []
    fdr_tracts = []

    while line != '':
        line = line.strip().split('\t')

        # bp
        power_temp = []
        tp = parse_list(line[10])[1:]
        p = parse_list(line[0])[1:]
        for x in range(len(tp)):
            if float(p[x]) != 0:
                power_temp.append(float(tp[x])/float(p[x]))
        m = mystats.mean(power_temp)
        if m != 'NA':
            power.append(m)

        fpr_temp = []
        all_pred_pos = parse_list(line[3])[1:]
        true_pred_pos  = parse_list(line[10])[1:]
        # false_pos = all_pred_pos - true_pred_pos
        nsites = float(params[int(i) - 1][7])
        pos_sites = parse_list(line[0])[1:]
        # actual negatives = nsites - pos_sites
        for x in range(len(all_pred_pos)):
            if nsites - float(pos_sites[x]) != 0:
                fpr_temp.append((float(all_pred_pos[x]) - float(true_pred_pos[x]))/\
                                    (nsites - float(pos_sites[x])))
        m = mystats.mean(fpr_temp)
        if m != 'NA':
            fpr.append(m)
        
        # fdr is false positives / total calls (i.e. proportion of
        # predictions that are correct)
        fdr_temp = []
        all_pred_pos = parse_list(line[3])[1:]
        true_pred_pos  = parse_list(line[10])[1:]
        # false_pos = all_pred_pos - true_pred_pos
        for x in range(len(all_pred_pos)):
            if float(all_pred_pos[x]) != 0:
                fdr_temp.append((float(all_pred_pos[x]) - float(true_pred_pos[x]))/\
                                    float(all_pred_pos[x]))
                assert(fdr_temp[-1]) <= 1, str(i) + '  '+ str(x) + ' ' + str(all_pred_pos[x]) + ' ' + str(true_pred_pos[x])
        m = mystats.mean(fdr_temp)
        if m != 'NA':
            fdr.append(m)
        
        
        # tracts
        power_tracts_temp = []
        tp = parse_list(line[12])
        p = parse_list(line[1])
        for x in range(len(tp)):
            if float(p[x]) != 0:
                power_tracts_temp.append(float(tp[x])/float(p[x]))
        m = mystats.mean(power_tracts_temp)
        if m != 'NA':
            power_tracts.append(m)
            
        fpr_tracts_temp = []
        all_pred_pos = parse_list(line[4])[1:]
        true_pred_pos = parse_list(line[11])[1:]
        neg_tracts = parse_list(line[2])[1:]
        for x in range(len(all_pred_pos)):
            if float(neg_tracts[x]) != 0:
                fpr_tracts_temp.append((float(all_pred_pos[x]) - float(true_pred_pos[x]))/\
                                           float(neg_tracts[x]))
        m = mystats.mean(fpr_tracts_temp)
        if m != 'NA':
            fpr_tracts.append(m)

        fdr_tracts_temp = []
        all_pred_pos = parse_list(line[4])[1:]
        true_pred_pos  = parse_list(line[11])[1:]
        # false_pos = all_pred_pos - true_pred_pos
        for x in range(len(all_pred_pos)):
            if float(all_pred_pos[x]) != 0:
                fdr_tracts_temp.append((float(all_pred_pos[x]) - float(true_pred_pos[x]))/\
                                    float(all_pred_pos[x]))
        m = mystats.mean(fdr_tracts_temp)
        if m != 'NA':
            fdr_tracts.append(m)

        line = f.readline()

    agg = [power, fpr, fdr, power_tracts, fpr_tracts, fdr_tracts]
    # mean row
    for p in params[int(i) - 1]:
        f_all.write(p + '\t')
    f_all.write('mean')
    for item in agg:
        m = mystats.mean(item)
        f_all.write('\t' + str(m))
    f_all.write('\n')
    # standard error row
    for p in params[int(i) - 1]:
        f_all.write(p + '\t')
    f_all.write('std_err')
    for item in agg:
        se = mystats.std_err(item)
        f_all.write('\t' + str(se))
    f_all.write('\n')
    # bootstrap
    bls = []
    bus = []
    for item in agg:
        bl, bu = mystats.bootstrap(item)
        print bl, bu
        bls.append(bl)
        bus.append(bu)
    print bls
    print bus
    # lower bootstrap row
    for p in params[int(i) - 1]:
        f_all.write(p + '\t')
    f_all.write('lower')
    for x in range(len(agg)):
        f_all.write('\t' + str(bls[x]))
    f_all.write('\n')
    # upper bootstrap row
    for p in params[int(i) - 1]:
        f_all.write(p + '\t')
    f_all.write('upper')
    for x in range(len(agg)):
        f_all.write('\t' + str(bus[x]))
    f_all.write('\n')
    f.close()

f_all.close()
