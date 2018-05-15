import ils_rosenberg
import mpmath as mp

mp.dps = 200
for n in [8000000./2, 8000000., 8000000.*2]:
    for t in [375000000./2, 375000000., 375000000.*2]:
        print t, n
        print "%0.50f" % (1 - ils_rosenberg.monophyletic_concordance_2(1, 94, t/n))
        print "%0.50f" % (1 - ils_rosenberg.monophyletic_concordance_2(2, 94, t/n))
        print "%0.50f" % (1 - ils_rosenberg.monophyletic_concordance_2(10, 94, t/n))
        print '***'

"""
T = 375000000 / 4N

    320000000

theta = mu * 4 * num_sites * N

rho factor = 7.425e-06
"""
