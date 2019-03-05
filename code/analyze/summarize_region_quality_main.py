import sys
import os
import gzip
import predict
from collections import defaultdict
from summarize_region_quality import *
import global_params as gp
import misc.read_fasta as read_fasta
import misc.read_table as read_table
import misc.seq_functions as seq_functions

# in main function for profiling purposes
def main():

    args = predict.process_predict_args(sys.argv[2:])

    task_ind = int(sys.argv[1])
    species_ind = task_ind / len(gp.chrms)
    chrm_ind = task_ind % len(gp.chrms) 

    species_from = args['states'][species_ind]
    chrm = gp.chrms[chrm_ind]

    gp_dir = '../'

    fn = gp.analysis_out_dir_absolute + args['tag'] + '/' + \
         'blocks_' + species_from + \
         '_' + args['tag'] + '_labeled.txt'

    # region_id strain chromosome predicted_species start end number_non_gap
    regions, labels = read_table.read_table_columns(fn, '\t')
    regions_chrm = defaultdict(list)
    for i in range(len(regions['region_id'])):
        if regions['chromosome'][i] == chrm:
            for key in regions.keys():
                regions_chrm[key].append(regions[key][i])
    n = len(regions_chrm['region_id'])

    for s in args['known_states']:
        regions_chrm['match_nongap_' + s] = [0 for i in range(n)]
        regions_chrm['num_sites_nongap_' + s] = [0 for i in range(n)]
        regions_chrm['match_hmm_' + s] = [0 for i in range(n)]
        regions_chrm['match_nonmask_' + s] = [0 for i in range(n)]
        regions_chrm['num_sites_nonmask_' + s] = [0 for i in range(n)]

    info_string_symbols = list('.-_npbcxNPBCX')
    for s in info_string_symbols:
        regions_chrm['count_' + s] = [0 for i in range(n)]

    strains = set(regions_chrm['strain'])

    # get masked sites for all references, not just the current
    # species_from we're considering regions from
    masked_sites_refs = {}
    for s in range(len(args['known_states'])):
        species_from_prefix = gp.ref_fn_prefix[args['known_states'][s]]
        masked_sites_refs[s] = \
            convert_intervals_to_sites(read_masked_intervals(gp_dir + \
                                                             gp.mask_dir + \
                                                             species_from_prefix + \
                                                             '_chr' + chrm + \
                                                             '_intervals.txt'))


    # loop through chromosomes and strains, followed by species of
    # introgression so that we only have to read each alignment in once
    ps_fn = gp.analysis_out_dir_absolute + args['tag'] + \
            '/positions_' + args['tag'] + '.txt.gz'
    ps_f = gzip.open(ps_fn, 'rb')

    # get chromosomes and strains from position file
    #for line in ps_f:
    while True:
        line = None
        try:
            line = ps_f.readline()
        except:
            pass
        if line == None or line == '':
            break
        line = line.split('\t')
        strain = line[0]
        current_chrm = line[1]
        if current_chrm != chrm:
            continue
        print strain, chrm

        # indices of alignment columns used by HMM
        ps = [int(x) for x in line[2:]]

        fn = gp_dir + gp.alignments_dir + \
             '_'.join(gp.alignment_ref_order) + '_' + strain + \
             '_chr' + chrm + '_mafft' + gp.alignment_suffix
        headers, seqs = read_fasta.read_fasta(fn)

        # to go from index in reference seq to index in alignment
        ind_align = []
        for s in range(len(seqs)):
            ind_align.append(index_alignment_by_reference(seqs[s]))

        masked_sites = \
            convert_intervals_to_sites(read_masked_intervals(gp_dir + \
                                                             gp.mask_dir + \
                                                             strain + \
                                                             '_chr' + chrm + \
                                                             '_intervals.txt'))
        masked_sites_ind_align = []
        for s in range(len(args['known_states'])):
            masked_sites_ind_align.append(\
                [ind_align[s][x] for x in masked_sites_refs[s]])
        masked_sites_ind_align.append([ind_align[-1][x] for x in masked_sites])

        # convert position indices from indices in master reference to
        # indices in alignment
        ps_ind_align = [ind_align[0][x] for x in ps]

        # loop through all regions for the specified chromosome and the
        # current strain
        for i in range(n):

            if regions_chrm['strain'][i] != strain:
                continue

            regions_dir = gp.analysis_out_dir_absolute + args['tag'] + '/regions/'
            if not os.path.isdir(regions_dir):
                os.mkdir(regions_dir)
            fn_region = regions_dir + regions_chrm['region_id'][i] + \
                        gp.fasta_suffix + '.gz'
            f_region = gzip.open(fn_region, 'wb')

            # calculate:
            # - identity with each reference
            # - fraction of region that is gapped/masked

            # index of start and end of region in aligned sequences
            slice_start = ind_align[0][int(regions_chrm['start'][i])]
            slice_end = ind_align[0][int(regions_chrm['end'][i])]
            assert slice_start in ps_ind_align, str(slice_start) + ' ' + \
                regions_chrm['start'][i] + ' ' + regions_chrm['region_id'][i]
            assert slice_end in ps_ind_align, str(slice_end) + ' ' + \
                regions_chrm['end'][i] + ' ' + regions_chrm['region_id'][i]

            seqx = seqs[-1][slice_start:slice_end + 1]

            # . = all match
            # - = gap in one or more sequences
            # p = matches predicted reference
            # 
            #info_all = dict(zip(args['known_states'], \
            #                    [['.' for c in seqx] for sj in args['known_states']]))
            info = [{'gap_any_flag':False, \
                     'mask_any_flag':False, \
                     'unseq_any_flag':False, \
                     'hmm_flag':False, \
                     'gap_flag':[], \
                     'mask_flag':[], \
                     'unseq_flag':[], \
                     'match_flag':[]} 
                    for k in range(len(seqx))]

            for sj in range(len(args['known_states'])):
                seqj = seqs[sj][slice_start:slice_end+1]
                statej = args['known_states'][sj]

                # only alignment columns used by HMM (polymorphic, no
                # gaps in any strain)
                total_match_hmm, total_sites_hmm, infoj = \
                    seq_id_hmm(seqj, seqx, slice_start, ps_ind_align)
                if statej == species_from or species_ind >= len(args['known_states']):
                    regions_chrm['num_sites_hmm'][i] = total_sites_hmm

                for k in range(len(seqx)):
                    info[k]['gap_any_flag'] = info[k]['gap_any_flag'] or \
                                              infoj['gap_flag'][k]
                    info[k]['unseq_any_flag'] = info[k]['unseq_any_flag'] or \
                                                infoj['unseq_flag'][k]
                    info[k]['hmm_flag'] = infoj['hmm_flag'][k]
                    info[k]['gap_flag'].append(infoj['gap_flag'][k])
                    info[k]['unseq_flag'].append(infoj['unseq_flag'][k])
                    info[k]['match_flag'].append(infoj['match'][k])


                regions_chrm['match_hmm_' + statej][i] = total_match_hmm

                # all alignment columns, excluding ones with gaps in
                # these two sequences
                total_match_nongap, total_sites_nongap = \
                    seq_functions.seq_id(seqj, seqx)
 
                regions_chrm['match_nongap_' + statej][i] = total_match_nongap
                regions_chrm['num_sites_nongap_' + statej][i] = total_sites_nongap

                # all alignment columns, excluding ones with gaps or
                # masked bases or unsequenced in *these two sequences*
                total_match_nonmask, total_sites_nonmask, infoj = \
                    seq_id_unmasked(seqj, seqx, slice_start, \
                                    masked_sites_ind_align[sj],\
                                    masked_sites_ind_align[-1])
                for k in range(len(seqx)):
                    info[k]['mask_any_flag'] = info[k]['mask_any_flag'] or \
                                               infoj['mask_flag'][k]
                    info[k]['mask_flag'].append(infoj['mask_flag'][k])

                regions_chrm['match_nonmask_' + statej][i] = \
                    total_match_nonmask
                regions_chrm['num_sites_nonmask_' + statej][i] = \
                    total_sites_nonmask

            names = args['known_states'] + [strain]
            for sj in range(len(names)):

                # write sequence to region alignment file, along with
                # start and end coordinates
                startj = bisect.bisect_left(ind_align[sj], slice_start)
                endj = bisect.bisect_left(ind_align[sj], slice_end)
                f_region.write('> ' + names[sj] + ' ' + \
                               str(startj) + ' ' + str(endj) + '\n')
                f_region.write(seqs[sj][slice_start:slice_end+1] + '\n')

            # also write string with info about each site
            info_string = make_info_string(info, 0, species_ind)
            f_region.write('> info\n')
            f_region.write(info_string + '\n')
            # and keep track of each symbol count
            for sym in info_string_symbols:
                regions_chrm['count_' + sym][i] = info_string.count(sym)

            sys.stdout.flush()

    ps_f.close()

    labels = labels + ['match_nongap_' + x for x in args['known_states']]
    labels = labels + ['num_sites_nongap_' + x for x in args['known_states']]
    labels = labels + ['match_hmm_' + x for x in args['known_states']]
    labels = labels + ['match_nonmask_' + x for x in args['known_states']]
    labels = labels + ['num_sites_nonmask_' + x for x in args['known_states']]
    labels = labels + ['count_' + x for x in info_string_symbols]

    fn = gp.analysis_out_dir_absolute + args['tag'] + '/' + \
         'blocks_' + species_from + \
         '_' + args['tag'] + '_chr' + chrm + '_quality.txt'

    f = open(fn, 'w')
    f.write('\t'.join(labels) + '\n')

    for i in range(n):
        #if regions_chrm['chromosome'][i] != chrm:
        #    continue
        f.write('\t'.join([str(regions_chrm[label][i]) for label in labels]))
        f.write('\n')
    f.close()

if __name__ == '__main__':
    main()
