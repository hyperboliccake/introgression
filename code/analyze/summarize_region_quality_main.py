import sys
import os
import gzip
from analyze import predict
from analyze.summarize_region_quality import (convert_intervals_to_sites,
                                              read_masked_intervals,
                                              index_alignment_by_reference,
                                              seq_id_hmm,
                                              seq_id_unmasked,
                                              make_info_string)
import global_params as gp
from misc import read_fasta
from misc import read_table
from misc import seq_functions
import numpy as np
import bisect


def main():

    args = predict.process_predict_args(sys.argv[2:])

    task_ind = int(sys.argv[1])
    species_ind = task_ind // len(gp.chrms)
    chrm_ind = task_ind % len(gp.chrms)

    species_from = args['states'][species_ind]
    chrm = gp.chrms[chrm_ind]

    base_dir = gp.analysis_out_dir_absolute + args['tag']

    regions_dir = f'{base_dir}/regions/'
    if not os.path.isdir(regions_dir):
        os.mkdir(regions_dir)

    # region_id strain chromosome predicted_species start end number_non_gap
    regions_chrm, labels = read_table.read_table_columns(
        f'{base_dir}/blocks_{species_from}_{args["tag"]}_labeled.txt', '\t',
        group_by='strain',
        chromosome=chrm
    )

    for r_id in regions_chrm:
        n = len(regions_chrm[r_id]['region_id'])

        for s in args['known_states']:
            regions_chrm[r_id]['match_nongap_' + s] = [0] * n
            regions_chrm[r_id]['num_sites_nongap_' + s] = [0] * n
            regions_chrm[r_id]['match_hmm_' + s] = [0] * n
            regions_chrm[r_id]['match_nonmask_' + s] = [0] * n
            regions_chrm[r_id]['num_sites_nonmask_' + s] = [0] * n

        info_string_symbols = list('.-_npbcxNPBCX')
        for s in info_string_symbols:
            regions_chrm[r_id]['count_' + s] = [0] * n

    # get masked sites for all references, not just the current
    # species_from we're considering regions from
    masked_sites_refs = {}
    for s, state in enumerate(args['known_states']):
        species_from_prefix = gp.ref_fn_prefix[state]
        masked_sites_refs[s] = \
            convert_intervals_to_sites(
                read_masked_intervals(
                    f'{gp.mask_dir}{species_from_prefix}'
                    f'_chr{chrm}_intervals.txt'))

    # loop through chromosomes and strains, followed by species of
    # introgression so that we only have to read each alignment in once
    with gzip.open(f'{base_dir}/positions_{args["tag"]}.txt.gz', 'rt') \
            as positions:
        for line in positions:
            line = line.split('\t')
            strain = line[0]
            current_chrm = line[1]
            if current_chrm != chrm or strain not in regions_chrm:
                continue

            print(strain, chrm)

            # indices of alignment columns used by HMM
            ps = np.array([int(x) for x in line[2:]])

            headers, seqs = read_fasta.read_fasta(
                gp.alignments_dir + '_'.join(gp.alignment_ref_order)
                + f'_{strain}_chr{chrm}_mafft{gp.alignment_suffix}')

            # to go from index in reference seq to index in alignment
            ind_align = []
            for seq in seqs:
                ind_align.append(index_alignment_by_reference(seq))

            masked_sites = convert_intervals_to_sites(
                read_masked_intervals(
                    f'{gp.mask_dir}{strain}_chr{chrm}_intervals.txt'))

            masked_sites_ind_align = []
            for s in range(len(args['known_states'])):
                masked_sites_ind_align.append(
                    ind_align[s][masked_sites_refs[s]])

            # add in sequence of query strain
            masked_sites_ind_align.append(
                ind_align[-1][masked_sites])

            # convert position indices from indices in master reference to
            # indices in alignment
            ps_ind_align = ind_align[0][ps]

            # loop through all regions for the specified chromosome and the
            # current strain
            for i in range(len(regions_chrm[strain]['region_id'])):
                r_id = regions_chrm[strain]['region_id'][i]
                start = regions_chrm[strain]['start'][i]
                end = regions_chrm[strain]['end'][i]

                # calculate:
                # - identity with each reference
                # - fraction of region that is gapped/masked

                # index of start and end of region in aligned sequences
                slice_start = ind_align[0][int(start)]
                slice_end = ind_align[0][int(end)]
                assert slice_start in ps_ind_align, \
                    f'{slice_start} {start} {r_id}'
                assert slice_end in ps_ind_align, \
                    f'{slice_end} {end} {r_id}'

                seqx = seqs[-1][slice_start:slice_end + 1]
                len_seqx = slice_end - slice_start + 1
                len_states = len(args['known_states'])

                # . = all match
                # - = gap in one or more sequences
                # p = matches predicted reference

                info = {'gap_any_flag': np.zeros((len_seqx), bool),
                        'mask_any_flag': np.zeros((len_seqx), bool),
                        'unseq_any_flag': np.zeros((len_seqx), bool),
                        'hmm_flag': np.zeros((len_seqx), bool),
                        'gap_flag': np.zeros((len_seqx, len_states), bool),
                        'mask_flag': np.zeros((len_seqx, len_states), bool),
                        'unseq_flag': np.zeros((len_seqx, len_states), bool),
                        'match_flag': np.zeros((len_seqx, len_states), bool)}

                for sj, statej in enumerate(args['known_states']):
                    seqj = seqs[sj][slice_start:slice_end+1]

                    # only alignment columns used by HMM (polymorphic, no
                    # gaps in any strain)
                    total_match_hmm, total_sites_hmm, infoj = \
                        seq_id_hmm(seqj, seqx, slice_start, ps_ind_align)

                    if statej == species_from \
                            or species_ind >= len(args['known_states']):
                        regions_chrm[strain]['num_sites_hmm'][i] = \
                            total_sites_hmm

                    # only write once, the first index
                    if sj == 0:
                        info['hmm_flag'] = infoj['hmm_flag']

                    info['gap_any_flag'] = np.logical_or(
                        info['gap_any_flag'], infoj['gap_flag'])
                    info['unseq_any_flag'] = np.logical_or(
                        info['unseq_any_flag'], infoj['unseq_flag'])
                    info['gap_flag'][:, sj] = infoj['gap_flag']
                    info['unseq_flag'][:, sj] = infoj['unseq_flag']
                    info['match_flag'][:, sj] = infoj['match']

                    regions_chrm[strain][f'match_hmm_{statej}'][i] = \
                        total_match_hmm

                    # all alignment columns, excluding ones with gaps in
                    # these two sequences
                    total_match_nongap, total_sites_nongap = \
                        seq_functions.seq_id(seqj, seqx)

                    regions_chrm[strain][f'match_nongap_{statej}'][i] =\
                        total_match_nongap
                    regions_chrm[strain][f'num_sites_nongap_{statej}'][i] =\
                        total_sites_nongap

                    # all alignment columns, excluding ones with gaps or
                    # masked bases or unsequenced in *these two sequences*
                    total_match_nonmask, total_sites_nonmask, infoj = \
                        seq_id_unmasked(seqj, seqx, slice_start,
                                        masked_sites_ind_align[sj],
                                        masked_sites_ind_align[-1])

                    info['mask_any_flag'] = np.logical_or(
                        info['mask_any_flag'], infoj['mask_flag'])
                    info['mask_flag'][:, sj] = infoj['mask_flag']

                    regions_chrm[strain][f'match_nonmask_{statej}'][i] = \
                        total_match_nonmask
                    regions_chrm[strain][f'num_sites_nonmask_{statej}'][i] = \
                        total_sites_nonmask

                region_writer = open(
                    f'{regions_dir}{r_id}{gp.fasta_suffix}', 'w')

                names = args['known_states'] + [strain]
                for sj in range(len(names)):
                    # write sequence to region alignment file, along with
                    # start and end coordinates
                    startj = bisect.bisect_left(ind_align[sj], slice_start)
                    endj = bisect.bisect_left(ind_align[sj], slice_end)

                    region_writer.write(f'> {names[sj]} {startj} {endj}\n')
                    region_writer.write(
                        ''.join(seqs[sj][slice_start:slice_end+1]) + '\n')

                # also write string with info about each site
                info_string = make_info_string(info, 0, species_ind)
                region_writer.write('> info\n')
                region_writer.write(info_string + '\n')
                region_writer.close()

                # TODO this can be made faster with numpy
                # and keep track of each symbol count
                for sym in info_string_symbols:
                    regions_chrm[strain]['count_' + sym][i] = \
                        info_string.count(sym)

            sys.stdout.flush()

    labels = labels + ['match_nongap_' + x for x in args['known_states']]
    labels = labels + ['num_sites_nongap_' + x for x in args['known_states']]
    labels = labels + ['match_hmm_' + x for x in args['known_states']]
    labels = labels + ['match_nonmask_' + x for x in args['known_states']]
    labels = labels + ['num_sites_nonmask_' + x for x in args['known_states']]
    labels = labels + ['count_' + x for x in info_string_symbols]

    with open(f'{base_dir}/blocks_{species_from}'
              f'_{args["tag"]}_chr{chrm}_quality.txt', 'w') as writer:

        writer.write('\t'.join(labels) + '\n')

        assert labels[0] == 'region_id', 'Unexpected labeled format'

        # reorganize output as list of tuples ordered by label
        output = []
        strains = list(regions_chrm.keys())
        for strain in strains:
            # pop to limit memory usage
            d = regions_chrm.pop(strain)
            output += list(zip(*[d[l] for l in labels]))

        # sort by region id (index 0, remove r)
        for entry in sorted(output, key=lambda e: int(e[0][1:])):
            writer.write('\t'.join([str(e) for e in entry]) + '\n')


if __name__ == '__main__':
    main()
