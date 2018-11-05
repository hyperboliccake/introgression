def read_origin(fn):
    f = open(fn)
    line = f.readline()
    seq = ''
    while line != '':
        if line.startswith('ORIGIN'):
            line = f.readline()
            while not line.startswith('//'):
                seq += ''.join(line.strip().split()[1:])
                line = f.readline()
            return seq
        line = f.readline()
    return None
