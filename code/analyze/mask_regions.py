import sys
import os
import gene_predictions
import annotate_regions
sys.path.insert(0, '../align')
import mask_helpers
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../misc')
import read_fasta
import write_fasta

# generate masked region alignment and annotated alignment files
# (replace masked bases with x)


def get_masked(p, start, end, strains):
    
    seqs = dict(zip(strains, [[] for s in strains]))
    a_start = p['ps'].index(str(start))
    a_end = p['ps'].index(str(end))
    for i in range(a_start, a_end + 1):
        for s in strains:
            c = ''
            if p[s + '_masked'][i] == '':
                c = p[s][i]
            else:
                c = gp.masked_symbol
            seqs[s].append(c)
    return [''.join(seqs[s]) for s in strains]

tag = sys.argv[1]
chrm = sys.argv[2]
suffix = ''
if len(sys.argv) == 4:
    suffix = sys.argv[3]

fn = gp.analysis_out_dir_absolute + tag + '/' + \
     'introgressed_blocks' + suffix + '_par_' + tag + '_summary.txt'
regions = gene_predictions.read_region_summary(fn)

strains = set([])
for region in regions:
    if regions[region]['chromosome'] == chrm:
        strains.add(regions[region]['strain'])
strains = list(strains)

#i = 0
for strain in strains:
    print strain
    sys.stdout.flush()

    fn_p_ann = gp.analysis_out_dir_absolute + tag + '/site_summaries/predictions_' + \
               strain + '_chr' + chrm + '_site_summary.txt.gz'
    p = annotate_regions.read_predictions_annotated(fn_p_ann)

    for region_id in regions:
        #if i % 100 == 0:
        #    sys.stdout.write(str(i) + '/' + str(len(regions)) + '\r')
        #    sys.stdout.flush()
        #i += 1

        if regions[region_id]['strain'] != strain:
            continue
        if regions[region_id]['chromosome'] != chrm:
            continue
               
        fn_align = gp.analysis_out_dir_absolute + tag + '/' + \
                   'regions/'  + region_id + '.maf.gz'
        headers, seqs = read_fasta.read_fasta(fn_align, gz=True)

        #seqs_masked = [] # entire chromosome seqs
        #for h in headers:
        #    fn_masked = h[:-1].split()[-1]
        #    fn_masked = fn_masked[:-len(gp.fasta_suffix)] + '_masked' + gp.fasta_suffix
        #    h, seq = read_fasta.read_fasta(fn_masked)
        #    seqs_masked.append(seq)

        start = int(regions[region_id]['start'])
        end = int(regions[region_id]['end'])

        strains = [h.split()[0][1:] for h in headers]
        seqs_masked = get_masked(p, start, end, strains)

        fn_align = gp.analysis_out_dir_absolute + tag + '/' + \
                   'regions/'  + region_id + '_masked.maf'
        write_fasta.write_fasta(headers, seqs_masked, fn_align, gz=True)
    
        #fn_align_ann = gp.analysis_out_dir_absolute + tag + '/' + \
            #               'regions/'  + region_id + '_annotated.txt.gz'

    
