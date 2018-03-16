# format output files to be read easily and plotted in R

import re
import sys
import os
import copy
import gene_predictions
sys.path.insert(0, '..')
import global_params as gp
sys.path.insert(0, '../sim/')
import sim_analyze_hmm_bw as sim
sys.path.insert(0, '../misc/')
import mystats

##======
# read in analysis parameters
##======

suffix = ''
if len(sys.argv == 3):
    suffix = sys.argv[1]

all_predict_args = [x.strip().split() for x in open(sys.argv[2], 'r').readlines()]
all_predict_args = [{'tag':x[0], 'improvement_frac':x[1], 'threshold':x[2], \
                     'expected_length':x[-2], 'expected_frac':x[-1]} \
                    for x in all_predict_args]

l = range(0,36)
l.remove(19)
l.remove(25)
l = [0]
all_predict_args = [all_predict_args[i] for i in l]

'''
finished = range(1,36)
finished.remove(2)
finished.remove(8)
finished.remove(14)
finished.remove(20)
finished.remove(26)
finished.remove(28)
finished.remove(32)
all_predict_args = [all_predict_args[i] for i in finished]
'''

# for things we might want to compare across tags, put in one table;
# for other things (which tend to be larger anyway), make one table
# per tag

sep = '\t'

##======
# for plot: lengths of all introgressed regions
##======

# one table for each tag
# strain chrm region_length

# one table for all tags
# tag improvement_frac threshold expected_length expected_frac avg_length lower upper median min max total_num_regions

print 'working on region lengths'

f = open(gp.analysis_out_dir_absolute + 'plot_region_lengths.txt', 'w')
for i in range(len(all_predict_args)):
    print '-', i
    args = all_predict_args[i]
    f_tag = open(gp.analysis_out_dir_absolute + args['tag'] + '/' + \
                 'plot_region_lengths' + suffix + '_' + args['tag'] + '.txt', 'w')
    fn = gp.analysis_out_dir_absolute + args['tag'] + '/' + \
         'introgressed_blocks_par' + suffix + '_' + args['tag'] + '_summary_plus.txt'
    region_summary = gene_predictions.read_region_summary(fn)
    lengths_all = []
    for region in region_summary:
        length = int(region_summary[region]['end']) - \
                 int(region_summary[region]['start']) + 1
        if int(region_summary[region]['number_match_ref2_not_ref1']) >= 5:
            f_tag.write(region + sep + region_summary[region]['strain'] + sep + \
                        region_summary[region]['chromosome'] + sep + \
                        str(length) + '\n')
            lengths_all.append(length)
    f_tag.close()
    f.write(args['tag'] + sep + args['improvement_frac'] + sep + \
            args['threshold'] + sep + args['expected_length'] + sep + \
            args['expected_frac'] + sep)
    f.write(str(mystats.mean(lengths_all)) + sep)
    bs_lower, bs_upper = mystats.bootstrap(lengths_all)
    f.write(str(bs_lower) + sep + str(bs_upper) + sep)
    f.write(str(mystats.median(lengths_all)) + sep)
    f.write(str(min(lengths_all)) + sep)
    f.write(str(max(lengths_all)) + sep)
    f.write(str(len(lengths_all)) + '\n')
f.close()

print 'done'

sys.exit()

##======
# for plot: number of genes per introgressed region
##======

# one table for each tag
# strain chrm region number_genes

# one table for all tags
# tag improvement_frac threshold expected_length expected_frac avg_genes_per_region lower upper median min max

print 'working on number of genes for each region'

f = open(gp.analysis_out_dir_absolute + 'plot_number_genes_by_region.txt', 'w')
for i in range(len(all_predict_args)):
    print '-', i
    args = all_predict_args[i]
    f_tag = open(gp.analysis_out_dir_absolute + args['tag'] + '/' + \
                 'plot_number_genes_by_region_' + args['tag'] + '.txt', 'w')
    fn = gp.analysis_out_dir_absolute + args['tag'] + '/' + \
         'genes_for_each_region_' + args['tag'] + '.txt'
    genes = gene_predictions.read_genes_for_each_region_summary(fn)
    fn = gp.analysis_out_dir_absolute + args['tag'] + '/' + \
         'introgressed_blocks_par_' + args['tag'] + '_summary.txt'
    region_summary = gene_predictions.read_region_summary(fn)
    num_genes_all = []
    for region in genes:
        f_tag.write(region + sep + region_summary[region]['strain'] + sep + \
                    region_summary[region]['chromosome'] + sep + \
                    genes[region]['num_genes'] + '\n')
        num_genes_all.append(int(genes[region]['num_genes']))
    f_tag.close()
    f.write(args['tag'] + sep + args['improvement_frac'] + sep + \
            args['threshold'] + sep + args['expected_length'] + sep + \
            args['expected_frac'] + sep)
    f.write(str(mystats.mean(num_genes_all)) + sep)
    bs_lower, bs_upper = mystats.bootstrap(num_genes_all)
    f.write(str(bs_lower) + sep + str(bs_upper) + sep)
    f.write(str(mystats.median(num_genes_all)) + sep)
    f.write(str(min(num_genes_all)) + sep)
    f.write(str(max(num_genes_all)) + '\n')
f.close()

print 'done'


##======
# for plot: number of introgressed bases for each strain
##======

# one table for all tags
# tag improvement_frac threshold expected_length expected_frac strain number_bases

print 'working on number of bases for each strain'

f = open(gp.analysis_out_dir_absolute + \
         'plot_number_introgressed_bases_by_strain.txt', 'w')
for i in range(len(all_predict_args)):
    print '-', i
    args = all_predict_args[i]
    fn = gp.analysis_out_dir_absolute + '/' + args['tag'] + '/' + \
         'regions_for_each_strain_' + args['tag'] + '.txt'
    regions = gene_predictions.read_regions_for_each_strain(fn)
    for strain in regions:
        total = 0
        for r in regions[strain]['region_list']:
            total += int(r[1])
        f.write(args['tag'] + sep + args['improvement_frac'] + sep + \
                args['threshold'] + sep + args['expected_length'] + sep + \
                args['expected_frac'] + sep + strain + sep + str(total) + '\n')
f.close()

print 'done'    

##======
# for plot: number of introgressed genes for each strain
##======

# one table for all tags
# tag improvement_frac threshold expected_length expected_frac strain number_genes

print 'working on number of genes for each strain'

f = open(gp.analysis_out_dir_absolute + \
         'plot_number_introgressed_genes_by_strain.txt', 'w')
for i in range(len(all_predict_args)):
    print '-', i
    args = all_predict_args[i]
    fn = gp.analysis_out_dir_absolute + '/' + args['tag'] + '/' + \
         'genes_for_each_strain_' + args['tag'] + '.txt'
    genes = gene_predictions.read_genes_for_each_strain(fn)
    for strain in genes:
        f.write(args['tag'] + sep + args['improvement_frac'] + sep + \
                args['threshold'] + sep + args['expected_length'] + sep + \
                args['expected_frac'] + sep + strain + sep + \
                genes[strain]['num_genes'] + sep + '\n')
f.close()

print 'done'

##======
# for plot: number of strains each gene introgressed in 
##======

# one table for each tag
# gene num_strains

# one table for all tags
# tag improvement_frac threshold expected_length expected_frac avg_strains_per_gene lower upper median min max total_num_genes total_num_genes_1 total_num_genes_>1

print 'working on number of strains for each gene'

f = open(gp.analysis_out_dir_absolute + 'plot_number_strains_by_genes.txt', 'w')
for i in range(len(all_predict_args)):
    print '-', i
    args = all_predict_args[i]
    f_tag = open(gp.analysis_out_dir_absolute + args['tag'] + '/' + \
                 'plot_number_strains_by_genes_' + args['tag'] + '.txt', 'w')
    fn = gp.analysis_out_dir_absolute + args['tag'] + '/' + \
         'strains_for_each_gene_' + args['tag'] + '.txt'
    strains = gene_predictions.read_strains_for_each_gene(fn)
    num_strains_all = []
    for gene in strains:
        f_tag.write(gene + sep + strains[gene]['num_strains'] + '\n')
        num_strains_all.append(int(strains[gene]['num_strains']))
    f_tag.close()
    f.write(args['tag'] + sep + args['improvement_frac'] + sep + \
            args['threshold'] + sep + args['expected_length'] + sep + \
            args['expected_frac'] + sep)
    f.write(str(mystats.mean(num_strains_all)) + sep)
    bs_lower, bs_upper = mystats.bootstrap(num_strains_all)
    f.write(str(bs_lower) + sep + str(bs_upper) + sep)
    f.write(str(mystats.median(num_strains_all)) + sep)
    f.write(str(min(num_strains_all)) + sep)
    f.write(str(max(num_strains_all)) + sep)
    f.write(str(len(num_strains_all)) + sep)
    f.write(str(len(filter(lambda x: x == 1, num_strains_all))) + sep)
    f.write(str(len(filter(lambda x: x > 1, num_strains_all))) + '\n')
f.close()

print 'done'

##======
# for plot: average fraction of each (introgressed) gene that's introgressed 
##======

# one table for each tag
# gene avg_frac_introgressed lower upper median min max

print 'working on fraction of gene introgressed'

for i in range(len(all_predict_args)):
    print '-', i
    args = all_predict_args[i]
    f_tag = open(gp.analysis_out_dir_absolute + args['tag'] + '/' + \
                 'plot_frac_introgressed_by_genes_' + args['tag'] + '.txt', 'w')
    fn = gp.analysis_out_dir_absolute + args['tag'] + '/' + \
         'strains_for_each_gene_' + args['tag'] + '.txt'
    strains = gene_predictions.read_strains_for_each_gene(fn)
    for gene in strains:
        fracs = [float(x[1]) for x in strains[gene]['strain_list']]
        f_tag.write(gene + sep)
        f_tag.write(str(mystats.mean(fracs)) + sep)
        bs_lower, bs_upper = mystats.bootstrap(fracs)
        f_tag.write(str(bs_lower) + sep + str(bs_upper) + sep)
        f_tag.write(str(mystats.median(fracs)) + sep)
        f_tag.write(str(min(fracs)) + sep)
        f_tag.write(str(max(fracs)) + '\n')
    f_tag.close()

print 'done'
