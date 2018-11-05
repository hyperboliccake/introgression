import ils_rosenberg
import mpmath as mp
import sys

mp.dps = 200

# correct
# 11.71875 * 4 * n / n = 46.875
T = 46.875
for Ti in [T/10, T/5, T/2, T, T*2, T*5, T*10]:
    #print "%0.50f" % (1 - ils_rosenberg.monophyletic_concordance_2(1, 2, Ti))
    print Ti, "%0.50f" % (1 - ils_rosenberg.monophyletic_concordance_2(1, 94, Ti))

sys.exit()

# test 
print "%0.50f" % ils_rosenberg.monophyletic_concordance_2(2, 2, 5)



# old
T = 375000000 * 4 
print "%0.50f" % (1 - ils_rosenberg.monophyletic_concordance_2(1, 2, T))
print "%0.50f" % (1 - ils_rosenberg.monophyletic_concordance_2(1, 2, T/2))
print "%0.50f" % (1 - ils_rosenberg.monophyletic_concordance_2(1, 2, T/10))
print "%0.50f" % (1 - ils_rosenberg.monophyletic_concordance_2(1, 2, T/100))

sys.exit()


for n in [8000000./2, 8000000., 8000000.*2]:
    for t in [375000000./2, 375000000., 375000000.*2]:
        print t, n
        print "%0.50f" % (1 - ils_rosenberg.monophyletic_concordance_2(1, 2, t*4))
        print "%0.50f" % (1 - ils_rosenberg.monophyletic_concordance_2(1, 2, t/n))
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
