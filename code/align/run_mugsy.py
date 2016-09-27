import os
sys.path.insert(0, '..')
import global_params as gp

# get all non-reference strains of cerevisiae and paradoxus
s = []
for d in gp.dirs_cer + gp.dirs_par:
    fns = os.listdir(d)
    # only look at fasta files in the directory
    fns = filter(lambda x: x.endswith(gp.fasta_suffix), fns)
    # only look at files containing '_chr' which should be chromosome
    # sequence files
    fns = filter(lambda x: '_chr' in x, fns)    
    num_files = len(fns)
    if num_files == 0:
        print 'found no chromosome sequence files in', d, '(perhaps you should check the _chr naming convention?)'
    fns = list(set([x[:x.find('_chr')] for x in fns]))
    num_strains = len(fns)
    assert num_files == num_strains * len(gp.chrms_roman), \
        'some strains in ' + d + ' have the wrong number of chromosome sequence files'
    entries = [(x, d) for x in fns]
    s += entries

a = []
if gp.resume_alignment:
    a = os.listdir(gp.alignments_dir)

# need to add this on the start of each command because os.system()
# creates a new shell instance every time
cmd_string_start = 'export MUGSY_INSTALL=' + gp.mugsy_install_path + '; '
cmd_string_start += 'export PATH=$PATH:$MUGSY_INSTALL:$MUGSY_INSTALL/mapping; '
cmd_string_start += 'export PERL5LIB=$MUGSY_INSTALL/perllibs; '

for strain, d in s:

    cmd_string = cmd_string_start
        
    for chrm in gp.chrms_roman:
        if (gp.cer_ref_strain + '_' + gp.par_ref_strain + '_' + x + '_chr' + chrm + gp.alignment_suffix) not in a:
            cmd_string += gp.mugsy_instal_path + '/mugsy ' + \
                '--directory ' + gp.alignments_dir + \
                '--prefix ' + gp.cer_ref_strain + '_' + par_ref_strain + '_' + x + '_chr' + chrm + ' ' + \
                gp.cer_ref_dir + '/' + gp.cer_ref_strain + '_chr' + chrm + gp.fasta_suffix + ' ' + \
                gp.par_ref_dir + '/' + gp.par_ref_strain + '_chr' + chrm + gp.fasta_suffix + ' ' + \
                d + '/' + strain + '_chr' + chrm + gp.fasta_suffix + '; '
    print x

    # commands can only be up to a certain length so break it up this way
    os.system(cmd_string)
