import math
import numpy.random

def mean(l):
    l = filter(lambda x: x != 'NA', l)
    if len(l) == 0:
        return 'NA'
    return float(sum(l)) / len(l)

def std_dev(l):
    l = filter(lambda x: x != 'NA', l)
    if len(l) == 0:
        return 'NA'
    if len(l) == 1:
        return 0
    m = mean(l)
    return math.sqrt(sum([(x - m)**2 for x in l]) / (len(l) - 1))

def std_err(l):
    l = filter(lambda x: x != 'NA', l)
    if len(l) == 0:
        return 'NA'
    return std_dev(l) / math.sqrt(len(l))

def bootstrap(l, n = 100, alpha = .05):
    l = filter(lambda x: x != 'NA', l)
    x = len(l)
    if x == 0:
        return 'NA', 'NA'
    a = []
    for i in range(n):
        a.append(mean(numpy.random.choice(l, size = x, replace = True)))
    a.sort()
    #print len(a), a.count(0)
    #print mean(a)
    return a[int(alpha * n * .5)], a[int((1 - alpha * .5) * n)]


