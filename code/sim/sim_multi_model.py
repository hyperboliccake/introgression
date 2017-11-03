# Run ms coalescent simulations under a variety of demographic models
# (actually, just continuous migration, after the last divergence, for
# now) 

# aiming for ~89% sequence id b/t cer and par and ~99.5% within cer,
# and ~70% between cer and bay (?)

import os
import sys
import process_args
sys.path.insert(0, '..')
import global_params as gp

args, last_read = process_args.process_args(sys.argv)

num_samples = args['num_samples_species_to'] + \
              args['num_samples_species_from1'] + \
              args['num_samples_species_from2']

gp_dir = '../'
outfilename = gp.sim_out_prefix + args['tag'] + gp.sim_out_suffix

migration1_time_scaling = 20

# start of ms command
# (in case you were thinking about it, DON'T subtract 1 from nsites for
# the -r option)
# the -p option is to insure there's enough granularity in segsite locations
ms_command = \
    gp.ms_install_path + '/ms ' + str(num_samples) + ' ' + str(args['num_reps']) + \
    ' -t ' + str(args['theta']) + \
    ' -r ' + str(args['rho']) + ' ' + str(args['num_sites']) + \
    ' -p 8'


# 2 species
if args['species_from2'] == None:
    join_time = args['topology'][2]
    ms_command += \
        ' -I 2 ' + str(args['num_samples_species_to']) + \
        ' ' + str(args['num_samples_species_from1'])
    ms_command += \
        ' -m 1 2 ' + str(args['migration_from1']) + \
        ' -em ' + str(join_time / float(migration1_time_scaling)) + ' 1 2 0' + \
        ' -em ' + str(join_time) + ' 1 2 0' # this is probably implied
    ms_command += \
        ' -ej ' + str(join_time) + ' 1 2'

# 3 species
else:
    if type(args['topology'][0]) != type([]):
        left = args['topology'][0]
        topology[0] = args['topology'][1]
        topology[1] = left
    most_recent_join_time = args['topology'][0][2]
    least_recent_join_time = args['topology'][2]
    last_to_join = args['topology'][1]
    first_to_join1 = args['topology'][0][0]
    first_to_join2 = args['topology'][0][1]

    label = {args['species_to']:'1', \
             args['species_from1']:'2', \
             args['species_from2']:'3'}

    # note that we need to keep the species in the order to, from1,
    # from2 (because we're assuming this is true in the analysis)
    ms_command += ' -I 3 ' + str(args['num_samples_species_to']) + ' ' + \
        str(args['num_samples_species_from1']) + ' ' + \
        str(args['num_samples_species_from2'])

    ms_command += \
        ' -m 1 2 ' + str(args['migration_from1'] * migration1_time_scaling) + \
        ' -m 1 3 ' + str(args['migration_from2'])  + \
        ' -em ' + str(most_recent_join_time / float(migration1_time_scaling)) + ' 1 2 0' + \
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
