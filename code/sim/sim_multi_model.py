# Run ms coalescent simulations under a variety of demographic models
# (actually, just continuous migration, after the last divergence, for
# now) 

# aiming for ~89% sequence id b/t cer and par and ~99.5% within cer,
# and ~70% between cer and bay (?)

import os
import sys
import sim_analyze_hmm_bw
sys.path.insert(0, '..')
import global_params as gp


tag, topology, species_to, species_from1, species_from2, \
    num_samples_species_to, num_samples_species_from1, num_samples_species_from2, \
    N0_species_to, N0_species_from1, N0_species_from2, \
    migration_from1, migration_from2, \
    expected_length_introgressed, \
    expected_num_introgressed_tracts, \
    has_ref_from1, has_ref_from2, \
    rho, outcross_rate, theta, num_sites, num_reps = \
    sim_analyze_hmm_bw.process_args(sys.argv)

num_samples = num_samples_species_to + num_samples_species_from1 + num_samples_species_from2

gp_dir = '../'
outfilename = gp.sim_out_prefix + tag + gp.sim_out_suffix

# start of ms command
# (in case you were thinking about it, DON'T subtract 1 from nsites for
# the -r option)
ms_command = \
    gp.ms_install_path + '/ms ' + str(num_samples) + ' ' + str(num_reps) + \
    ' -t ' + str(theta) + \
    ' -r ' + str(rho) + ' ' + str(num_sites)

# 2 species
if species_from2 == None:
    join_time = topology[2]
    ms_command += \
        ' -I 2 ' + str(num_samples_species_to) + ' ' + str(num_samples_species_from1)
    ms_command += \
        ' -m 1 2 ' + str(migration_from1) + \
        ' -em ' + str(join_time) + ' 1 2 0' # this is probably implied
    ms_command += \
        ' -ej ' + str(join_time) + ' 1 2'

# 3 species
else:
    if type(topology[0]) != type([]):
        left = topology[0]
        topology[0] = topology[1]
        topology[1] = left
    most_recent_join_time = topology[0][2]
    least_recent_join_time = topology[2]
    last_to_join = topology[1]
    first_to_join1 = topology[0][0]
    first_to_join2 = topology[0][1]

    label = {species_to:'1', species_from1:'2', species_from2:'3'}

    # note that we need to keep the species in the order to, from1,
    # from2 (because we're assuming this is true in the analysis)
    ms_command += ' -I 3 ' + str(num_samples_species_to) + ' ' + \
        str(num_samples_species_from1) + ' ' + str(num_samples_species_from2)

    ms_command += \
        ' -m 1 2 ' + str(migration_from1) + \
        ' -m 1 3 ' + str(migration_from2) + \
        ' -em ' + str(most_recent_join_time) + ' 1 2 0' + \
        ' -em ' + str(most_recent_join_time) + ' 1 3 0'

    ms_command += \
        ' -ej ' + str(most_recent_join_time) + ' ' + \
        label[first_to_join2] + ' ' + label[first_to_join1]
    ms_command += \
        ' -ej ' + str(least_recent_join_time) + ' ' + \
        label[last_to_join] + ' ' + label[first_to_join1]

ms_command += ' -T > ' + gp_dir + gp.sim_out_dir + '/ms/' + outfilename

print(ms_command)
os.system(ms_command)


# divergence of cerevisae and paradoxus ~5 million ya - Kellis et al. 2003
#
# mutation rate ~1.84 * 10^-10 - Fay & Benavides 2005
# mutation rate ~2.2 * 10^-10 - Tsai et al. 2008
#
# generations per year ~2920 - Fay & Benavides 2005
#
# divergence of vineyard and sake cerevisiae (wine/sake) ~12000 ya - Fay & Benavides 2005
# (could easily be an order of magnitude older)
#
# population size of paradoxus europe ~8600000 - Tsai et al. 2008
#
# population size of paradoxus far east ~7200000 - Tsai et al. 2008
#
# population size of cerevisiae europe 8000000 - guess from paradoxus
#
# population size of cerevisiae far east 8000000 - guess from paradoxus
#
# frequency of sex in paradoxus europe 1/1000 generations - Tsai et al. 2008
#
# frequency of sex in paradoxus far east 1/3000 generations - Tsai et al. 2008
#
# frequency of sex in cerevisiae 1/2000 - 1/9000 generations - Tsai et al. 2008
#
# divergence of paradoxus european and eastenr ~500000 ya - (guesstimate from phylogeny in) Liti et al. 2009
# 
