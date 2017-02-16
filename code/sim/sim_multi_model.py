# Run ms coalescent simulations under a variety of demographic models

# aiming for ~89% sequence id b/t cer and par and ~99.5% within cer

import os
import sys
import sim_analyze_hmm_bw
sys.path.insert(0, '..')
import global_params as gp


tag, model, N0, include_bay, include_unk,\
        num_samples_cer, num_samples_par, num_samples_bay,\
        par_cer_migration, bay_cer_migration,\
        t_cer_par, t_cer_bar_bay,\
        num_sites, rho, theta, outcross_rate, num_reps = \
        sim_analyze_hmm_bw.process_args(sys.argv)

num_samples = num_samples_cer + num_samples_par + num_samples_bay

# TODO? make a script that generates param files based on args given
# in global param file

outfilename = gp.sim_out_prefix + tag + gp.sim_out_suffix

# start of ms command
# (in case you were thinking about it, DON'T subtract 1 from nsites for
# the -r option)
ms_command = \
    gp.ms_install_path + '/ms ' + str(num_samples) + ' ' + str(num_reps) + \
    ' -t ' + str(theta) + \
    ' -r ' + str(rho) + ' ' + str(num_sites) + \
    ' -I 2 ' + str(num_samples_cer) + ' ' + str(num_samples_par)


# introgression happens continuously from par to cer between the
# present and t_cer_par (joining of par and cer lineages)
if model == 'C':
    ms_command += \
        ' -m 1 2 ' + str(par_cer_migration) + \
        ' -em ' + str(t_cer_par) + ' 1 2 0' # this is probably implied

# introgression happens in one pulse, at some time before t_cer_par
# specified by the number after D
elif model[0] == 'I':
    # time of pulse of migration
    t_mig = float(model[1:])
    # pulse of migration for 1 generation of time
    par_cer_migration *= t_cer_par * 2 * N0 # need to scale by 2 N0
    ms_command += \
        ' -em ' + str(t_mig / (2.0 * N0)) + ' 1 2 ' + str(par_cer_migration) + \
        ' -em ' + str((t_mig + 1) / (2.0 * N0)) + ' 1 2 0'

# fail
else:
    print 'incorrect model selection'
    sys.exit()

ms_command += \
    ' -ej ' + str(t_cer_par) + ' 1 2' + \
    ' -T > ' + gp.sim_out_dir + '/' + outfilename

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
