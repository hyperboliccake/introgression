import sys
import gzip
sys.path.insert(0, '..')
import global_params as gp


def convert(s1, s2):
    i2 = -1
    i2d = 0
    l = []
    for i in range(len(s1)):
        if s2[i] == gp.gap_symbol:
            i2d += 1
        else:
            i2 += 1
            i2d = 0
        if s1[i] != gp.gap_symbol:
            if i2d == 0:
                l.append(str(i2))
            else:
                l.append(str(i2) + '.' + str(i2d))
    return l

def write_coordinates(l, fn):
    f = gzip.open(fn, 'wb')
    f.write('\n'.join([str(x) for x in l]))
    f.write('\n')
    f.close()

