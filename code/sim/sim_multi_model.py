# Run ms coalescent simulations under a variety of demographic models

import os
import sys

tag = sys.argv[1]
model = sys.argv[2]
outdir = '../../results/sim/'

# parameters specific to all of the models

# all population sizes are the same
N0 = int(sys.argv[3])

num_samples_par = int(sys.argv[4])
num_samples_cer_1 = int(sys.argv[5])
num_samples_cer_2 = int(sys.argv[6])
num_samples_cer = num_samples_cer_1 + num_samples_cer_2
num_samples = num_samples_par + num_samples_cer

# migration parameter is 2 * N0 * m, where mij is fraction of i made
# up of j each generation; need to figure out how to make migration
# rates equivalent for different models
par_cer_migration = 2 * N0 * float(sys.argv[7])

# in generations
t_cer = float(sys.argv[8]) / (2 * N0)
t_par_cer = float(sys.argv[9]) / (2 * N0)

# 13,500 sites to get about 10% with one recombination event, .3% with
# more than one (based on poisson(.1), 1 recombination per chromosome
# of average length 750,000)
num_sites = float(sys.argv[10])

# parameter is recombination rate between adjacent bp per generation
# should probably be 1/750000 + 6.1 * 10^-6 (where 750000 is average
# chr size)
rho = 2 * N0 * float(sys.argv[11]) * (num_sites - 1)

outcross_rate = float(sys.argv[12])

rho *= outcross_rate


mu = 1.84 * 10 ** -10
theta = mu * 2 * num_sites * N0

num_reps = int(sys.argv[13])

outfilename = 'sim_out_' + tag + '.txt'

# start of ms command
# (in case you were thinking about it, DON'T subtract 1 from nsites for
# the -r option)
ms_command = \
    '/net/gs/vol1/home/aclark4/software/msdir/ms ' + str(num_samples) + ' ' + str(num_reps) + \
    ' -t ' + str(theta) + \
    ' -r ' + str(rho) + ' ' + str(num_sites) + \
    ' -I 3 ' + str(num_samples_par) + ' ' + str(num_samples_cer_1) + ' ' + str(num_samples_cer_2)

# introgression happens continuously from par to cer1 between the
# present and t_cer (joining of the two cer lineages)
if model == 'A':
    ms_command += \
        ' -m 2 1 ' + str(par_cer_migration) + \
        ' -em ' + str(t_cer) + ' 2 1 0'

# introgression happens continuously from par to cer1 between t_cer
# (joining of the two cer lineages) and t_par_cer (joining of par and
# cer lineages)
elif model == 'B':
    # scale so that there's the same total amount of migration
    par_cer_migration *= t_cer/float(t_par_cer - t_cer)
    ms_command += \
        ' -em ' + str(t_cer) + ' 2 1 ' +  str(par_cer_migration) + \
        ' -em ' + str(t_par_cer) + ' 2 1 0' # this is probably implied

# introgression happens continuously from par to cer1 between the
# present and t_par_cer (joining of par and cer lineages)
elif model == 'C':
    # scale so that there's the same total amount of migration
    par_cer_migration *= t_cer/float(t_par_cer)
    ms_command += \
        ' -m 2 1 ' + str(par_cer_migration) + \
        ' -em ' + str(t_par_cer) + ' 2 1 0' # this is probably implied

# introgression happens in one pulse, at some time before t_par_cer
# specified by the number after D
elif model[0] == 'D':
    # time of pulse of migration
    t_mig = float(model[1:])
    # pulse of migration for 1 generation of time
    par_cer_migration *= t_cer * 2 * N0 # need to scale by 2 N0
    ms_command += \
        ' -em ' + str(t_mig / (2.0 * N0)) + ' 2 1 ' + str(par_cer_migration) + \
        ' -em ' + str((t_mig + 1) / (2.0 * N0)) + ' 2 1 0'

# fail
else:
    print 'incorrect model selection'
    sys.exit

ms_command += \
    ' -ej ' + str(t_cer) + ' 3 2' + \
    ' -ej ' + str(t_par_cer) + ' 2 1' + \
    ' -T > ' + outdir + outfilename

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
