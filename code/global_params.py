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

## now specified in setup_args file

#====
# alignment files
#====

## alignments directory now specified in setup_args file

mask_dir = '../alignments/masked/'
#mask_dir = '/tigress/tcomi/aclark4_temp/par4/masked/'

# should we leave the alignments already completed in the alignments
# directory alone?
resume_alignment = True

# master_ref now automatically assumed to be first
# reference specified in setup_args file

#====
# HMM
#====

match_symbol = '+'
mismatch_symbol = '-'
unknown_symbol = '?'

unsequenced_symbol = 'n'
gap_symbol = '-'
unaligned_symbol = '?'
masked_symbol = 'x'

#====
# simulations
#====

# output directory for simulpations
sim_out_dir_absolute = '/tigress/tcomi/aclark4_temp/results/sim'
#sim_out_dir_absolute = '/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/sim/'

# prefix for simulation output
sim_out_prefix = 'sim_out_'

# suffix for simulation output
sim_out_suffix = '.txt'

#====
# analysis
#====

analysis_out_dir_absolute = \
    '/tigress/AKEY/akey_vol2/aclark4/projects/introgression/results/analysis/'

regions_out_dir_absolute = analysis_out_dir_absolute + '/regions/'

genes_out_dir_absolute = analysis_out_dir_absolute + '/genes/'

#====
# software install locations
#====

mugsy_install_path = '/tigress/anneec/software/mugsy/'

tcoffee_install_path = '/tigress/anneec/software/T-COFFEE_installer_Version_11.00.8cbe486_linux_x64/bin/'

mafft_install_path = '/tigress/anneec/software/mafft/bin/'

ms_install_path = '/tigress/anneec/software/msdir/'

# including dustmasker
blast_install_path = '/tigress/anneec/software/ncbi-blast-2.7.1+-src/c++/ReleaseMT/bin/'

orffinder_install_path = '/tigress/anneec/software/'

ldselect_install_path = '/tigress/anneec/software/ldSelect/'

structure_install_path = '/tigress/anneec/software/structure/'

#====
# other
#====

chrms = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI']
#chrms = ['I']

chrms_ara = dict(zip(chrms, range(1, len(chrms)+1)))
