# this is for adding a few columns to introgressed_blocks_par_summary file:
# - number of genes it overlaps
# - longest stretch of gaps

import sys
import os
import gzip
import gene_predictions
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../misc/')
import read_fasta

cen_starts = [151465, 238207, 114385, 449711, 151987, 148510, 496920, 105586, 355629, 436307, 440129, 150828, 268031, 628758, 326584, 555957]
cen_starts = [x-1 for x in cen_starts]

cen_ends =[151582,238323,114501,449821,152104,148627,497038,105703,355745,436425,440246,150947,268149,628875,326702,556073]
cen_ends = [x-1 for x in cen_ends]

tel_coords = [1,801,229411,230218,1,6608,812379,813184,1,1098,315783,316620,1,904,1524625,1531933,1,6473,569599,576874,1,5530,269731,270161,1,781,1083635,1090940,1,5505,556105,562643,1,7784,439068,439888,1,7767,744902,745751,1,807,665904,666816,1,12085,1064281,1078177,1,6344,923541,924431,1,7428,783278,784333,1,847,1083922,1091291,1,7223,942396,948010]
tel_coords = [x-1 for x in tel_coords]
tel_left_starts = [tel_coords[i] for i in range(0, len(tel_coords), 4)]
tel_left_ends = [tel_coords[i] for i in range(1, len(tel_coords), 4)]
tel_right_starts = [tel_coords[i] for i in range(2, len(tel_coords), 4)]
tel_right_ends = [tel_coords[i] for i in range(3, len(tel_coords), 4)]

def distance_from_telomere(start, end, chrm):

    assert start <= end, str(start) + ' ' + str(end)

    i = gp.chrms.index(chrm)
    # region entirely on left arm
    if end <= cen_starts[i]:
        return start - tel_left_ends[i]
    # region entirely on right arm
    if start >= cen_ends[i]:
        return tel_right_starts[i] - end
    # region overlaps centromere: return minimum distance from either telomere
    return min(start - tel_left_ends[i], tel_right_starts[i] - end)

def distance_from_centromere(start, end, chrm):

    assert start <= end, str(start) + ' ' + str(end)

    i = gp.chrms.index(chrm)
    # region entirely on left arm
    if end <= cen_starts[i]:
        return cen_starts[i] - end
    # region entirely on right arm
    if start >= cen_ends[i]:
        return start - cen_ends[i]
    # region overlaps centromere: return 0
    return 0

def write_region_summary_plus(fn, regions, fields):
    f = open(fn, 'w')
    f.write('region_id\t' + '\t'.join(fields) + '\n')
    keys = sorted(regions.keys(), key=lambda x: int(x[1:]))
    for region_id in keys:
        f.write(region_id + '\t')
        f.write('\t'.join([str(regions[region_id][field]) for field in fields]))
        f.write('\n')
    f.close()

def gap_columns(seqs):
    g = 0
    for i in range(len(seqs[0])):
        for seq in seqs:
            if seq[i] == gp.gap_symbol:
                g += 1
                break
    return g

def longest_consecutive(s, c):
    max_consecutive = 0
    current_consecutive = 0
    in_segment = False
    for i in range(len(s)):
        if s[i] == c:
            current_consecutive += 1
            in_segment = True
        else:
            if in_segment:
                max_consecutive = max(max_consecutive, current_consecutive)
                current_consecutive = 0
            in_segment = False
    return max_consecutive

tag = sys.argv[1]
suffix = ''
if len(sys.argv) == 3:
    suffix = sys.argv[2]

fn = gp.analysis_out_dir_absolute + tag + '/' + \
     'introgressed_blocks' + suffix + '_par_' + tag + '_summary.txt'
# copy pasta :(
fields = ['strain', 'chromosome', 'predicted_species', 'start', 'end', \
          'number_non_gap', 'number_match_ref1', 'number_match_ref2', \
          'number_match_only_ref1', 'number_match_ref2_not_ref1', \
          'number_mismatch_all_ref']
regions = gene_predictions.read_region_summary(fn)

region_genes = {}
for chrm in gp.chrms:
    fn_genes = gp.analysis_out_dir_absolute + tag + '/' + \
               'genes_for_each_region_chr' + chrm + '_' + tag + '.txt'
    d = gene_predictions.read_genes_for_each_region_summary(fn_genes)
    region_genes.update(d)

fn_out = gp.analysis_out_dir_absolute + tag + '/' + \
        'introgressed_blocks' + suffix + '_par_' + tag + '_summary_plus.txt'
fields.append('aligned_length')
fields.append('number_genes')
fields.append('number_gaps')
fields.append('longest_gap')
fields.append('distance_from_telomere')
fields.append('distance_from_centromere')


i = 0
for region_id in regions:
    if i % 100 == 0:
        sys.stdout.write(str(i) + '/' + str(len(regions)) + '\r')
        sys.stdout.flush()
    i += 1

    regions[region_id]['number_genes'] = region_genes[region_id]['num_genes']

    fn_align = gp.analysis_out_dir_absolute + tag + '/' + \
               'regions/'  + region_id + '.maf.gz'
    headers, seqs = read_fasta.read_fasta(fn_align, gz=True)

    regions[region_id]['aligned_length'] = len(seqs[0])

    number_gaps = gap_columns(seqs)
    regions[region_id]['number_gaps'] = number_gaps

    longest_gaps = [longest_consecutive(seq, gp.gap_symbol) for seq in seqs]
    regions[region_id]['longest_gap'] = max(longest_gaps)

    regions[region_id]['distance_from_telomere'] = \
        distance_from_telomere(int(regions[region_id]['start']), \
                               int(regions[region_id]['end']), \
                               regions[region_id]['chromosome'])

    regions[region_id]['distance_from_centromere'] = \
        distance_from_centromere(int(regions[region_id]['start']), \
                                 int(regions[region_id]['end']), \
                                 regions[region_id]['chromosome'])
    
write_region_summary_plus(fn_out, regions, fields)
    
