import os
import sys
import copy
import re
import numpy.random

sys.path.insert(0, '../hmm/')
from hmm_bw import *

resume = False

def get_seqs_chrm(rc, rp, x, chrm):
    # for now, just make each alignment block a separate sequence, but
    # maybe check at some point whether any of them should be
    # concatenated?

    # what to do about columns with gaps? don't really want to count
    # two gapped columns as matching...
    # ignore all columns where x is gapped, count gaps in refs as mismatches


    # coding:
    # 0 -> matches cer ref but not par ref
    # 1 -> matches par ref but not cer ref
    # 2 -> matches cer ref and par ref
    # 3 -> doesn't match cer ref or par ref


    seqs = []
    scores = []
    psx = []
    f = open('../../alignments/genbank/' + rc + '_' + rp + '_' + x + '_chr' + chrm + '.maf', 'r')
    f.readline()
    n = -1
    while True:
        line = f.readline()
        if line == '':
            print scores
            print len(scores)
            return seqs, psx

        if line[0] == '#' or line[0] == '\n':
            continue

        if line[0] == 'a':
            m = re.search('score=(?P<score>[0-9\.]+)', line)
            scores.append(float(m.group('score')))
            m = re.search('mult=(?P<n>[0-9\.]+)', line)
            n = int(m.group('n'))
            
        if line[0] == 's':
            if n == 1:
                continue
            
            seqc = ''
            seqp = ''
            seqx = ''
            srcx = ''
            startx = -1
            sizex = -1
            strandx = ''

            # first sequence in alignment block
            s, src, start, size, strand, src_size, text = line.split()
            if rc in src:
                seqc = text
            elif rp in src:
                seqp = text
            else:
                assert x in src, x + ' ' + src
                srcx = src
                seqx = text
                # note that mugsy indexes from 0 and that for the
                # reverse strand, it indexes from the other end
                startx = int(start)
                sizex = int(size)
                strandx = strand

            # second sequence in alignment block
            line = f.readline()
            assert line[0] == 's', line[:100]            
            s, src, start, size, strand, src_size, text = line.split()
            if rc in src:
                seqc = text
            elif rp in src:
                seqp = text
            else:
                assert x in src, x + ' ' + src
                srcx = src
                seqx = text
                startx = int(start)
                sizex = int(size)
                strandx = strand

            if n == 2:
                # don't care about cases where only the two references align
                if seqx == '':
                    continue
                seq = ''
                # if cer reference doesn't align
                if seqc == '':
                    for i in range(len(seqp)):
                        if seqx[i] == '-':
                            continue
                        elif seqp[i] == seqx[i]:
                            seq += '1'
                        else:
                            seq += '3'
                # if par reference doesn't align
                else:
                    for i in range(len(seqc)):
                        if seqx[i] == '-':
                            continue
                        elif seqc[i] == seqx[i]:
                            seq += '0'
                        else:
                            seq += '3'
                assert len(seq) == sizex, str(len(seq)) + ' ' + str(sizex)
                seqs.append(seq)
                psx.append((srcx, startx, strandx))
                continue

            # (optional) third sequence in alignment block
            line = f.readline()
            assert line[0] == 's', line[:100]            
            s, src, start, size, strand, src_size, text = line.split()
            if rc in src:
                seqc = text
            elif rp in src:
                seqp = text
            else:
                assert x in src, x + ' ' + src
                srcx = src
                seqx = text
                startx = int(start)
                sizex = int(size)
                strandx = strand

            seq = ''
            for i in range(len(seqc)):
                if seqx[i] == '-':
                    continue
                elif seqc[i] == seqx[i]:
                    if seqp[i] == seqx[i]:
                        seq += '2'
                    else:
                        seq += '0'
                else:
                    if seqp[i] == seqx[i]:
                        seq += '1'
                    else:
                        seq += '3'
            assert len(seq) == sizex, str(len(seq)) + ' ' + str(sizex)
            seqs.append(seq)
            psx.append((srcx, startx, strandx))

def get_seqs(rc, rp, x, chrms):
    all_seqs = []
    all_psx = []
    for chrm in chrms:
        seqs, psx = get_seqs_chrm(rc, rp, x, chrm)
        all_seqs += seqs
        all_psx += psx
    return all_seqs, all_psx

def predict_introgressed_hmm(seqs):

    # create a hidden markov model and determine which reference genome we
    # are most likely to be in at each variant site
    hmm = HMM()

    hmm.set_init([.85,.15])

    # 0,1 would also work here
    hmm.set_states(['cer', 'par'])

    hmm.set_trans([[.9997,.0003],[.05,.95]])

    # combined error in parent sequences and observed sequence -
    # in this case 'error' comes from amout of difference
    # reasonable to see between reference and other strains of
    # same species; remember we've coded 0 for cer (non
    # introgressed) and 1 for par (introgressed)
    # hmm.set_emis({'cer':{0:.95, 1:.05},'par':{0:.05, 1:.95}}) 
    #hmm.set_emis([{'0':.5, '1':.0001, '2':.4998, '3':.0001},{'0':.0001, '1':.5, '2':.4998, '3':.0001}])
    hmm.set_emis([{'0':.1, '1':.001, '2':.88, '3':.019},{'0':.0017, '1':.06, '2':.93, '3':.0083}])

    predicted = []
    print len(seqs), 'seqs'
    sys.stdout.flush()
    for i in range(len(seqs)):
        if i % 100 == 0:
            print 'viterbi on seq', i
            sys.stdout.flush()
        obs_seq = seqs[i]
        hmm.set_obs(obs_seq)
        # this returns a list of state indices, so 0 and 1 for not introgressed and introgressed
        predicted.append(hmm.viterbi())

    return predicted, hmm

def predict_introgressed_id(seqs):
    # ok so the aligned sequences can have gaps
    window_size = 1000
    window_shift = 500
    thresholds = [.7, .95, .96]
    predicted = []
    for seq in seqs:
        p = [0 for b in xrange(len(seq))]
        for i in range(0, len(seq) - window_size, window_shift):
            region = seq[i:i+window_size]
            count0 = region.count('0')
            count1 = region.count('1')
            count2 = region.count('2')
            #count3 = window_size - count0 - count1 - count2
            cer_match = float(count0 + count2) / len(region)
            par_match = float(count1 + count2) / len(region) 
            if cer_match > thresholds[0]:
                if cer_match < thresholds[1] and par_match > thresholds[2]:
                    # this automatically appends to the list if we've
                    # gone over the length, but we'll deal with that
                    # later
                    p[i:i+window_size] = [1] * window_size 
        predicted.append(p[:len(seq)])

    return predicted
 

ref_cer = 'S288c'
ref_par = 'CBS432'

c = ['Sigma1278b.fa', 'SK1.fa', 'yjm1078.fa', 'yjm1083.fa', 'yjm1129.fa', 'yjm1133.fa', 'yjm1190.fa', 'yjm1199.fa', 'yjm1202.fa', 'yjm1208.fa', 'yjm1242.fa', 'yjm1244.fa', 'yjm1248.fa', 'yjm1250.fa', 'yjm1252.fa', 'yjm1273.fa', 'yjm1304.fa', 'yjm1307.fa', 'yjm1311.fa', 'yjm1326.fa', 'yjm1332.fa', 'yjm1336.fa', 'yjm1338.fa', 'yjm1341.fa', 'yjm1342.fa', 'yjm1355.fa', 'yjm1356.fa', 'yjm1381.fa', 'yjm1383.fa', 'yjm1385.fa', 'yjm1386.fa', 'yjm1387.fa', 'yjm1388.fa', 'yjm1389.fa', 'yjm1399.fa', 'yjm1400.fa', 'yjm1401.fa', 'yjm1402.fa', 'yjm1415.fa', 'yjm1417.fa', 'yjm1418.fa', 'yjm1419.fa', 'yjm1433.fa', 'yjm1434.fa', 'yjm1439.fa', 'yjm1443.fa', 'yjm1444.fa', 'yjm1447.fa', 'yjm1450.fa', 'yjm1460.fa', 'yjm1463.fa', 'yjm1477.fa', 'yjm1478.fa', 'yjm1479.fa', 'yjm1526.fa', 'yjm1527.fa', 'yjm1549.fa', 'yjm1573.fa', 'yjm1574.fa', 'yjm1592.fa', 'yjm1615.fa', 'yjm189.fa', 'yjm193.fa', 'yjm195.fa', 'yjm244.fa', 'yjm248.fa', 'yjm270.fa', 'yjm271.fa', 'yjm320.fa', 'yjm326.fa', 'yjm428.fa', 'yjm450.fa', 'yjm451.fa', 'yjm453.fa', 'yjm456.fa', 'yjm470.fa', 'yjm541.fa', 'yjm554.fa', 'yjm555.fa', 'yjm681.fa', 'yjm682.fa', 'yjm683.fa', 'yjm689.fa', 'yjm693.fa', 'yjm969.fa', 'yjm972.fa', 'yjm975.fa', 'yjm978.fa', 'yjm981.fa', 'yjm984.fa', 'yjm987.fa', 'yjm990.fa', 'yjm993.fa', 'yjm996.fa']
c = [x[:-3] for x in c]

#wine/european should be 41
group1 = ['yjm1244', 'yjm189', 'yjm1356', 'yjm1242', 'yjm1477', 'yjm1332', 'yjm453', 'yjm244', 'yjm1526', 'yjm1129', 'yjm1336', 'yjm1417', 'yjm1415', 'yjm1341', 'yjm969', 'yjm993', 'yjm984', 'yjm990', 'yjm972', 'yjm987', 'yjm975', 'yjm978', 'yjm981', 'yjm996', \
              'yjm248', 'yjm1078', 'yjm1252', 'yjm270', 'yjm1574', 'yjm271', 'yjm193', 'yjm1387', 'yjm1463', 'yjm1383', 'yjm1250', 'yjm1355', 'yjm1450', 'yjm1386', 'yjm1326', 'yjm1385', 'yjm1549']

# other should be 53, or actually 52 because leaving out s228c
group2 = ['yjm693', 'yjm1478', 'yjm1433', 'yjm689', 'yjm1399', 'yjm1190', 'yjm1381', 'yjm1338', 'yjm682', \
              'yjm683', 'yjm326', 'yjm450', 'yjm1527', 'yjm1202', 'yjm1199', 'yjm1311', 'yjm428', 'yjm1419', 'yjm451', \
              'yjm1615', 'yjm1208', 'yjm1083', 'yjm681', 'yjm554', 'yjm541', 'yjm320', 'yjm555', 'yjm1304', 'yjm1133', \
              'yjm456', 'yjm470', 'yjm1307', 'yjm1444', 'yjm627', 'yjm195', 'yjm1439', 'yjm1248', 'yjm1418', 'yjm1447', \
              'yjm1342', 'yjm1479', 'yjm1400', 'yjm1401', 'yjm1460', 'yjm1592', 'yjm1389', 'yjm1388', 'yjm1443', 'yjm1573', \
              'yjm1402', 'yjm1434', 'yjm1273']


group_all = group1 + group2

all_seqs = []
all_ps = []

chrms_roman = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XIV']

print 'start'
for x in group_all:
    print x
    seqs, psx = get_seqs(ref_cer, ref_par, x, chrms_roman)
    all_seqs += seqs
    all_ps += psx
sys.stdout.flush()
fi = open('../../results/introgressed_hmm.txt', 'w')
predicted, hmm = predict_introgressed_hmm(all_seqs)
try:
    assert len(predicted) == len(all_seqs)
except:
    print 'problem'
for s in xrange(len(predicted)):
    print '====', all_ps[s]
    sys.stdout.flush()
    in_int = False
    start_int = 0
    end_int = 0
    for i in xrange(len(predicted[s])):
        if predicted[s][i] == 1:
            if in_int:
                end_int = i
            else:
                in_int = True
                start_int = i
                end_int = i
        elif in_int:
            # 
            fi.write(all_ps[s][0] + ', ' + all_ps[s][2] + ' strand, ' + str(all_ps[s][1] + start_int) + '-' + str(all_ps[s][1] + end_int) + ', ' + str(all_ps[s][1]) + ', ' + str(len(all_seqs[s])) + '\n')
            assert all_ps[s][1] + start_int >= all_ps[s][1] and all_ps[s][1] + end_int
            in_int = False
            fi.flush()
    if in_int:
        fi.write(all_ps[s][0] + ', ' + all_ps[s][2] + ' strand, ' + str(all_ps[s][1] + start_int) + '-' + str(all_ps[s][1] + end_int) + ', ' + str(all_ps[s][1]) + ', ' + str(len(all_seqs[s])) + '\n')
        fi.flush()
fi.close()


