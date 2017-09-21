import sys
import math
import cPickle
import mpmath

def F(k, a, b, c, d):
    '''returns the probability that at least one interspecific coalescence
    occurs as the a lineages from species A, b lineages from species B,
    and c lineages from species C coalesce to a total of k lineages AND
    that the most recent interspecific coalescence joins a lineage from
    species A and a lineage from species B'''

    assert(k >= 1)

    if d.has_key((k, a, b, c)):
        return d[(k, a, b, c)]

    if a < 0 or b < 0 or c < 0 or a + b + c <= k:
        return 0

    sum_choose_2 = choose_2(a + b + c)
    result =  a * b / sum_choose_2 + \
        F(k, a - 1, b, c, d) * choose_2(a) / sum_choose_2 + \
        F(k, a, b - 1, c, d) * choose_2(b) / sum_choose_2 + \
        F(k, a, b, c - 1, d) * choose_2(c) / sum_choose_2        
    d[(k, a, b, c)] = result
    return result
        
def choose_2(x):
    return x * (x - 1) / 2

def falling_factorial(a, k):
    product = 1
    for i in range(0, k):
        product *= a - i
    return product

def rising_factorial(a, k):
    product = 1
    for i in range(0, k):
        product *= a + i
    return product

def g(i, j, T, d):
    '''returns probability that i lineages derive from j lineages that
    existed T coalescent time units in the _past_'''


    if d.has_key((i, j, T)):
        return d[(i, j, T)]

    # note that result will only be nonzero for i >= j
    total = 0
    for k in range(j, i + 1):
        k = mpmath.mpf(k)
        total +=  mpmath.mpf(math.e) ** (-k * (k - 1) * T / 2.) * \
                  (2 * k - 1) * (-1) ** (k - j) * \
                  rising_factorial(j, k - 1) * \
                  falling_factorial(i, k) / \
                  (mpmath.mpf(math.factorial(j)) * \
                       math.factorial(k - j) * \
                       rising_factorial(i, k))

    d[(i, j, T)] = total

    return total

def W(m, n, x, k, T, d, dg):
    ''' returns the probability that X1 = x, X2 = k - x | X1 + X2 = k,
    where X1 and X2 are the number of lineages ancestral to species A
    and B respectively at time T3 + T, and m and n are the numbers of
    ancestral lineages for the two species at time T3'''

    if d.has_key((m, n, x, k, T)):
        return d[(m, n, x, k, T)]

    denom = 0
    for i in range(1, k):
        denom += g(m, i, T, dg) * g(n, k - i, T, dg)

    result = g(m, x, T, dg) * g(n, k - x, T, dg) / denom

    d[(m, n, x, k, T)] = result

    return result

def delta(k, i):
    if k == i:
        return 1
    return 0
        
# THIS WAS SUPPOSED TO BE TOPOLOGICAL BUT ENDED UP PRODUCING VALUES
# THAT MATCH TAKAHATA... FML
def takahata_concordance(r, s, q, T3, T2, F_dic, g_dic, W_dic):
    '''returns probability of takahata (???) concordance'''

    total = 0
    for m in range(1, r + 1):
        m = mpmath.mpf(m)
        for n in range(1, s + 1):
            n = mpmath.mpf(n)
            for k in range(1, m + n + 1):
                k = mpmath.mpf(k)
                print m, n, k
                xsum = 0
                for x in range(1, k):
                    x = mpmath.mpf(x)
                    lsum = 0
                    for l in range(1, q): # MISSING +1
                        l = mpmath.mpf(l)
                        lsum += g(q, l, T3 + T2, g_dic) * F(1, x, k-x, l, F_dic)
                    xsum += W(m, n, x, k, T2, W_dic, g_dic) * lsum

                total += g(r, m, T3, g_dic) * g(s, n, T3, g_dic) * g(m + n, k, T2, g_dic) * \
                         (F(k, m, n, 0, F_dic) + \
                          (1 - F(k, m, n, 0, F_dic)) * \
                          xsum)

    return total

def topological_concordance(r, s, q, T3, T2, F_dic, g_dic, W_dic):
    '''returns probability of topological concordance'''

    total = 0
    for m in range(1, r + 1):
        m = mpmath.mpf(m)
        for n in range(1, s + 1):
            n = mpmath.mpf(n)
            for k in range(1, m + n + 1):
                k = mpmath.mpf(k)
                print m, n, k
                xsum = 0
                for x in range(1, k):
                    x = mpmath.mpf(x)
                    lsum = 0
                    for l in range(1, q + 1):
                        l = mpmath.mpf(l)
                        lsum += g(q, l, T3 + T2, g_dic) * F(1, x, k-x, l, F_dic)
                    xsum += W(m, n, x, k, T2, W_dic, g_dic) * lsum

                total += g(r, m, T3, g_dic) * g(s, n, T3, g_dic) * g(m + n, k, T2, g_dic) * \
                         (F(k, m, n, 0, F_dic) + \
                          (1 - F(k, m, n, 0, F_dic)) * \
                          xsum)

    return total

def monophyletic_concordance(r, s, q, T3, T2, F_dic, g_dic, W_dic):
    '''returns probability of monophyletic concordance'''

    total = 0
    for m in range(1, r + 1):
        m = mpmath.mpf(m)
        for n in range(1, s + 1):
            n = mpmath.mpf(n)
            for k in range(1, m + n + 1):
                k = mpmath.mpf(k)
                print m, n, k
                for l in range(1, q + 1):
                    l = mpmath.mpf(l)
                    xsum = 0
                    for x in range(1, k):
                        x = mpmath.mpf(x)
                        xsum += W(m, n, x, k, T2, W_dic, g_dic) * (1 - F(3, x, k - x, l, F_dic) - \
                            F(3, x, l, k - x, F_dic) - F(3, k - x, l, x, F_dic)) * 1 / 3.

                    total += g(r, m, T3, g_dic) * g(s, n, T3, g_dic) * g(m + n, k, T2, g_dic) * g(q, l, T3 + T2, g_dic) * \
                        (delta(k, 1) * (1 - F(2, m, n, 0, F_dic)) * 2 / (l * (l + 1)) + \
                             (1 - delta(k, 1)) * (1 - F(k, m, n, 0, F_dic)) * xsum)
                        
    return total

def monophyletic_concordance_2(r, s, T3, F_dic = {}, g_dic = {}):
    '''returns probability of monophyletic concordance for two species'''

    total = 0
    r = mpmath.mpf(r)
    s = mpmath.mpf(s)
    T3 = mpmath.mpf(T3)
    for m in range(1, r + 1):
        m = mpmath.mpf(m)
        for n in range(1, s + 1):
            n = mpmath.mpf(n)
            total += g(r, m, T3, g_dic) * \
                g(s, n, T3, g_dic) * \
                (1 - F(2, m, n, 0, F_dic))
    return total
    
def concordance(r, s, q, T3, T2):    

    # load dictionaries of calculated values from files if they exist
    F_dic = {}
    try:
        f = open('F_vals.txt', 'r')
        F_dic = cPickle.load(f)
        f.close()
    except:
        pass

    g_dic = {}
    try:
        f = open('g_vals.txt', 'r')
        g_dic = cPickle.load(f)
        f.close()
    except:
        pass

    W_dic = {}
    try:
        f = open('W_vals.txt', 'r')
        W_dic = cPickle.load(f)
        f.close()
    except:
        pass

    # calculate different types of concordance
    mc = monophyletic_concordance(r, s, q, T3, T2, F_dic, g_dic, W_dic)
    tc = topological_concordance(r, s, q, T3, T2, F_dic, g_dic, W_dic)

    # save new dictionaries of calculations
    f = open('F_vals_new.txt', 'w')
    cPickle.dump(F_dic, f)
    f.close()

    f = open('g_vals_new.txt', 'w')
    cPickle.dump(g_dic, f)
    f.close()

    f = open('W_vals_new.txt', 'w')
    cPickle.dump(W_dic, f)
    f.close()

    return {'topological concordance':tc, 'monophyletic concordance':mc}

#results =  concordance(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]))
#print results['topological concordance']
#print results['monophyletic concordance']

mpmath.mp.dps = 50
for i in range(11, 16):
    mpmath.nprint(monophyletic_concordance_2(10, 100, 375000000/8000000./i), 50)
#mpmath.nprint(monophyletic_concordance_2(10, 100, 375000000/8000000.), 50)
#for i in range(2, 11):
#    mpmath.nprint(monophyletic_concordance_2(10, 100, 375000000/8000000.*i), 50)
