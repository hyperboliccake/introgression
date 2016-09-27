#####
# file extensions
#####

# suffix for _all_ fasta files
fasta_suffix = '.fa'

# suffix for _all_ alignment files; this needs to match the suffix output by mugsy
alignment_suffix = '.maf'

#####
# sequence locations/names
#####

# directory for cerevisiae reference
cer_ref_dir = '/net/akey/vol2/aclark4/nobackup/100_genomes/genomes/S288c_SGD-R64/'

# name of cerevisiae reference
cer_ref_strain = 'S288c_SGD-R64'

# directory for paradoxus reference
par_ref_dir = '/net/akey/vol2/aclark4/nobackup/100_genomes/paradoxus/strains/CBS432/assembly/CBS432/'

# name of paradoxus reference
par_ref_strain = 'CB432'

# all directories containing non-reference cerevisiae sequences;
# program will assume that _all_ fasta files in each of these
# directories should be used
dirs_cer = ['/net/akey/vol2/aclark4/nobackup/100_genomes/genomes_gb/']

# all directories containing non-reference paradoxus sequences;
# program will assume that _all_ fasta files in each of these
# directories should be used
dirs_par = []

#####
# alignment files
#####

# directory for mugsy to output alignments and for subsequent parts of
# the program to read the alignments; this directory needs to exist
# before running the program
alignments_dir = '../../alignments/genbank/'

# should we leave the alignments already completed in the alignments
# directory alone?
resume_alignment = True

#####
# software install locations
#####

mugsy_install_path = '~/software/mugsy'

#####
# other
#####

chrms_roman = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XIV']

