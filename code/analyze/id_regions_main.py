import sys
from analyze import predict
from collections import defaultdict
from operator import itemgetter
import global_params as gp


def main():
    args = predict.process_predict_args(sys.argv[1:])

    # order regions by chromosome, start (break ties alphabetically by strain)
    all_regions_by_chrm = dict(zip(gp.chrms, [[] for chrm in gp.chrms]))
    fns = {}
    for species_from in args['states']:

        # strain chromosome predicted_species start end number_non_gap
        fn = gp.analysis_out_dir_absolute + args['tag'] + '/' + \
             'blocks_' + species_from + \
             '_' + args['tag'] + '.txt'

        # introgressed regions keyed by strain and then chromosome:
        # (start, end, number_non_gap)
        regions = predict.read_blocks(fn)

        for strain in regions.keys():
            for chrm in regions[strain].keys():
                for entry in regions[strain][chrm]:
                    start, end, number_non_gap = entry
                    all_regions_by_chrm[chrm].append((start, end, number_non_gap, \
                                                      strain, species_from))

        fns[species_from] = fn[:-4] + '_labeled.txt'

    fs = {}
    for species_from in args['states']:
        fs[species_from] = open(fns[species_from], 'w')
        fs[species_from].write('region_id\tstrain\tchromosome\tpredicted_species\tstart\tend\tnum_sites_hmm\n')

    idc = 1
    for chrm in gp.chrms:
        for entry in sorted(all_regions_by_chrm[chrm], key = itemgetter(0, 3)):
            species_from = entry[4]
            rid = 'r' + str(idc)
            fs[species_from].write(rid + '\t' + entry[3] + '\t' + \
                                   chrm + '\t' + species_from + '\t' + \
                                   str(entry[0]) + '\t' + str(entry[1]) + \
                                   '\t' + str(entry[2]) + '\n')
            idc += 1

    for species_from in args['states']:
        fs[species_from].close()


if __name__ == "__main__":
    main()
