import gzip

def read_table_rows(fn, sep):

    f = None
    if fn.endswith('.gz'):
        f = gzip.open(fn, 'rb')
    else:
        f = open(fn, 'r')
    labels = f.readline()[:-1].split(sep)
    regions = [line[:-1].split(sep) for line in f.readlines()]
    d = {}
    for region in regions:
        d[region[0]] = dict(zip(labels[1:], region[1:]))
    f.close()
    return d, labels

def read_table_columns(fn, sep):

    f = None
    if fn.endswith('.gz'):
        f = gzip.open(fn, 'rb')
    else:
        f = open(fn, 'r')
    labels = f.readline()[:-1].split(sep)
    d = dict(zip(labels, [[] for l in labels]))
    line = f.readline()
    while line != '':
        line = line[:-1].split(sep)
        for i in range(len(labels)):
            d[labels[i]].append(line[i])
        line = f.readline()
    f.close()
    return d, labels
