# TODO migrate the sim args file into here and create a script that
# converts values in here to that file

#====
# biological parameters
#====

mu = 1.84 * 10 ** -10

#====
# file extensions
#====

# suffix for _all_ fasta files
fasta_suffix = '.fa'

# suffix for _all_ alignment files; this needs to match the suffix output by mugsy
alignment_suffix = '.maf'

#====
# sequence locations/names
#====

ref_dir = {'S288c':'/tigress/AKEY/akey_vol2/aclark4/nobackup/100_genomes/genomes/S288c_SGD-R64/', \
                'CBS432':'/tigress/anneec/projects/introgression/data/CBS432/'}


# will look at all *.fa files in these directories; expects filenames
# in the format strain_chrX.fa
non_ref_dirs = {'S288c':['/tigress/AKEY/akey_vol2/aclark4/nobackup/100_genomes/genomes_gb/'], \
                    'CBS432':[]}

# genbank file for _all_ species/strains
gb_all = '/tigress/AKEY/akey_vol2/aclark4/nobackup/100_genomes/sequence.gb'

#master_ref = 'S288c'

# gb master
gb_master_dir = '/tigress/AKEY/akey_vol2/aclark4/nobackup/'

#====
# alignment files
#====

# order that references appear in filename
alignment_ref_order = ['S288c', 'CBS432']

# this is kind of dumb, but the name of the fasta files to get seqs
# from won't necessarily just be the names we want to use
ref_fn_prefix = {'S288c':'S288c_SGD-R64', 'CBS432':'CBS432'}

# directory for mugsy to output alignments and for subsequent parts of
# the program to read the alignments; this directory needs to exist
# before running the program
alignments_dir = '../alignments/genbank/'

# should we leave the alignments already completed in the alignments
# directory alone?
resume_alignment = True

# overall indexing reference
master_ref = alignment_ref_order[0]

#====
# HMM
#====

match_symbol = '+'
mismatch_symbol = '-'
unknown_symbol = '?'

unsequenced_symbol = 'N'
gap_symbol = '-'
unaligned_symbol = '?'

#====
# simulations
#====

# output directory for simulpations
sim_out_dir = '../results/sim/'

# prefix for simulation output
sim_out_prefix = 'sim_out_'

# suffix for simulation output
sim_out_suffix = '.txt'

'''
# set of parameters for all simulations
sim_params = []
sim_params.append({'tag':'1', \
                       'topology':('cer', 'par', 375000000), \
                       'species_to':'cer', \
                       'species_from1':'par', \
                       'species_from2':None, \
                       'num_samples_species_to':100, \
                       'num_samples_species_from1':10, \
                       'num_samples_species_from2':0, \
                       'N0_samples_species_to':8000000, \
                       'N0_samples_species_from1':8000000, \
                       'N0_samples_species_from2':0, \
                       'migration_from1':0, \
                       'migration_from2':0, \
                       'expected_tract_lengths':{'par':100}, \
                       'expected_num_tracts':{'par':2}, \
                       'has_ref_from1':True, \
                       'has_ref_from1':False, \
                       'rho': 7.425e-06, \
                       'outcross_rate': 2e-05, \
                       'theta': None, \
                       'num_sites': 50000, \
                       'num_reps': 50\
                       })
'''

#====
# analysis
#====

analysis_out_dir_absolute = \
    '/tigress/anneec/projects/introgression/results/analysis/'

regions_out_dir_absolute = analysis_out_dir_absolute + '/regions/'

genes_out_dir_absolute = analysis_out_dir_absolute + '/genes/'

#====
# software install locations
#====

mugsy_install_path = '/tigress/anneec/software/mugsy/'

tcoffee_install_path = '/tigress/anneec/software/T-COFFEE_installer_Version_11.00.8cbe486_linux_x64/bin/'

mafft_install_path = '/tigress/anneec/software/mafft/bin/'

ms_install_path = '/tigress/anneec/software/msdir/'

#====
# other
#====

chrms = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI']

