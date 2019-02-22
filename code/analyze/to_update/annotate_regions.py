import gzip
import gene_predictions
import sys
import global_params as gp
sys.path.insert(0, '../misc/')



def get_block_by_site(all_regions, seq):

    introgressed_by_site = [' ' for site in seq]
    return introgressed_by_site
    for regions in all_regions:
        for entry in regions:
            start_ind = gene_predictions.index_ignoring_gaps(seq, entry[0], 0)
            end_ind = gene_predictions.index_ignoring_gaps(seq, entry[1], 0)
            for i in range(start_ind, end_ind+1):
                introgressed_by_site[i] = entry[-1]
    return introgressed_by_site


def write_predictions_annotated(alignment_headers, alignment_seqs, master, \
                                strain_labels, match_by_site, \
                                gene_by_site, block_by_site, masked, fn):

    f = gzip.open(fn, 'wb')
    num_seqs = len(alignment_seqs)
    num_align_cols = len(alignment_seqs[0])
    ref_inds = range(0, len(alignment_seqs) - 1)
    ind_ref = -1
    ind_ref_sub = 0
    ind_strain = -1
    ind_strain_sub = 0
    sep = '\t'

    individual_indices = [0] * num_seqs

    # header
    f.write('ps_ref' + sep + 'ps_strain' + sep + \
            sep.join(strain_labels) + sep + \
            'match' + sep + \
            'gene' + sep + 'block' + sep + \
            sep.join([lab + '_masked' for lab in strain_labels]) + '\n')

    lines = []
    # one line for each alignment column
    for i in range(num_align_cols):
        line = ''

        # index in reference
        ps_ref = None
        if alignment_seqs[master][i] == gp.gap_symbol:
            ind_ref_sub += 1
            ps_ref = str(ind_ref) + '.' + str(ind_ref_sub)
        else:
            ind_ref_sub = 0
            ind_ref += 1
            ps_ref = str(ind_ref)
        line += ps_ref + sep
        
        # index in strain
        ps_strain = None
        if alignment_seqs[-1][i] == gp.gap_symbol:
            ind_strain_sub += 1
            ps_strain = str(ind_strain) + '.' + str(ind_strain_sub)
        else:
            ind_strain_sub = 0
            ind_strain += 1
            ps_strain = str(ind_strain)
        line += ps_strain + sep

        for s in alignment_seqs:
            line += s[i] + sep

        for r in ref_inds:
            line += match_by_site[r][i]
        line += sep

        if gene_by_site[i] != None:
            line += gene_by_site[i]
        line += sep

        if block_by_site[i] != ' ':
            line += block_by_site[i]

        for si in range(num_seqs):
            line += sep
            if alignment_seqs[si][i] != gp.gap_symbol:
                # TODO update n to x
                if masked[si][individual_indices[si]] == 'n': #gp.masked_symbol:
                    line += gp.masked_symbol
                individual_indices[si] += 1
            
        line += '\n'
        
        lines.append(line)
    f.writelines(lines)

    f.close()

# TODO give this a more general name/place
def read_predictions_annotated(fn):
    sep = '\t'
    f = gzip.open(fn, 'rb')
    labels = f.readline()[:-1].split(sep)
    d = dict(zip(labels, [[] for l in labels]))
    line = f.readline()
    while line != '':
        line = line[:-1].split(sep)
        for i in range(len(labels)):
            d[labels[i]].append(line[i])
        #d[line[0]] = dict(zip(labels[1:], line[1:]))
        line = f.readline()
    f.close()
    return d



