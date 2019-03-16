import sys
from analyze import predict
from analyze import read_args
from operator import itemgetter
import global_params as gp


def main():
    args = read_args.process_predict_args(sys.argv[1:])

    # order regions by chromosome, start (break ties alphabetically by strain)
    all_regions_by_chrm = dict(zip(gp.chrms, [[] for chrm in gp.chrms]))
    output_files = {}
    base_dir = gp.analysis_out_dir_absolute + args['tag']
    for species_from in args['states']:

        # strain chromosome predicted_species start end number_non_gap
        fn = f'{base_dir}/blocks_{species_from}_{args["tag"]}.txt'

        # introgressed regions keyed by strain and then chromosome:
        # (start, end, number_non_gap)
        regions = predict.read_blocks(fn)

        for strain in regions:
            for chrm in regions[strain]:
                for entry in regions[strain][chrm]:
                    start, end, number_non_gap = entry
                    all_regions_by_chrm[chrm].append(
                        (start, end, number_non_gap, strain, species_from))

        output_files[species_from] = f'{fn[:-4]}_labeled.txt'

    writers = {}
    for species_from in args['states']:
        writers[species_from] = open(output_files[species_from], 'w')
        writers[species_from].write(
            'region_id\tstrain\tchromosome\tpredicted_species\t'
            'start\tend\tnum_sites_hmm\n')

    idc = 1
    for chrm in gp.chrms:
        for entry in sorted(all_regions_by_chrm[chrm], key=itemgetter(0, 3)):
            (start, end, number_non_gap, strain, species_from) = entry
            writers[species_from].write(
                f'r{idc}\t{strain}\t{chrm}\t{species_from}\t'
                f'{start}\t{end}\t{number_non_gap}\n')
            idc += 1

    for species_from in args['states']:
        writers[species_from].close()


if __name__ == "__main__":
    main()
