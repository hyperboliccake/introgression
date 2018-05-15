import gzip
def read_table_rows(fn, sep, header=True, key_ind=0):
    # returns dictionary of rows keyed by first item in row

    f = None
    if fn.endswith('.gz'):
        f = gzip.open(fn, 'rb')
    else:
        f = open(fn, 'r')
    labels = None
    if header:
        labels = f.readline()[:-1].split(sep)
    regions = [line[:-1].split(sep) for line in f.readlines()]
    d = {}
    for region in regions:
        if header:
            d[region[key_ind]] = \
                dict(zip(labels[:key_ind] + labels[key_ind + 1:], \
                         region[:key_ind] + region[key_ind + 1:]))
        else:
            d[region[key_ind]] = region[:key_ind] + region[key_ind + 1:]
    f.close()
    return d, labels

def read_table_columns(fn, sep):
    # returns dictionary of columns, keyed by labels

    f = None
    if fn.endswith('.gz'):
        f = gzip.open(fn, 'rb')
    else:
        f = open(fn, 'r')
    labels = f.readline()[:-1].split(sep)
    d = dict(zip(labels, [[] for l in labels]))
    lines = f.readlines()
    f.close()
    for line in lines:
        line = line[:-1].split(sep)
        for i in range(len(labels)):
            d[labels[i]].append(line[i])
    return d, labels
