import re

def read_mugsy(fn, required_mult = 1):
    f = open(fn, 'r')
    line = f.readline()
    while line[0] == '#':
        line = f.readline()
    blocks = {}
    while line != '':
        assert line[0] == 'a', line
        block = {}
        m = re.search('a score=(?P<score>[0-9]+) ' +\
                          'label=(?P<label>[a-zA-Z0-9]+) ' + \
                          'mult=(?P<mult>[0-9]+)', line)
        block['mult'] = int(m.group('mult'))
        if block['mult'] >= required_mult:
            block['score'] = int(m.group('score'))
            block['strains'] = {}
            for i in range(block['mult']):
                line = f.readline()
                assert line[0] == 's'
                line = line.strip().split() # splits on space and tab
                name = line[1][:line[1].find('_')]
                block['strains'][name] = {}
                block['strains'][name]['start'] = int(line[2])
                block['strains'][name]['length'] = int(line[3])
                block['strains'][name]['strand'] = line[4]
                block['strains'][name]['aligned_length'] = int(line[5])
                block['strains'][name]['sequence'] = line[6]
                block['strains'][name]['index'] = i # for sorting
            blocks[m.group('label')] = block
        line = f.readline()
        while line != '' and line[0] != 'a':
            line = f.readline()
    f.close()
    return blocks

def read_mugsy_block(label, fn):
    f = open(fn, 'r')
    line = f.readline()
    while line[0] == '#':
        line = f.readline()
    while line != '':
        assert line[0] == 'a', line
        block = {}
        m = re.search('a score=(?P<score>[0-9]+) ' +\
                          'label=(?P<label>[a-zA-Z0-9]+) ' + \
                          'mult=(?P<mult>[0-9]+)', line)
        if m.group('label') == label:
            block['mult'] = int(m.group('mult'))
            block['score'] = int(m.group('score'))
            block['strains'] = {}
            for i in range(block['mult']):
                line = f.readline()
                assert line[0] == 's'
                line = line.strip().split() # splits on space and tab
                name = line[1][:line[1].find('_')]
                block['strains'][name] = {}
                block['strains'][name]['start'] = int(line[2])
                block['strains'][name]['length'] = int(line[3])
                block['strains'][name]['strand'] = line[4]
                block['strains'][name]['aligned_length'] = int(line[5])
                block['strains'][name]['sequence'] = line[6]
                block['strains'][name]['index'] = i # for sorting
            f.close()
            return block
        line = f.readline()
        while line != '' and line[0] != 'a':
            line = f.readline()
    f.close()

def make_ungapped(blocks, strain):

    for block in blocks:
        if strain in block['strains']:
            gaps = []
            n = 0
            ungapped = ''
            for pos in xrange(len(block['strains'][strain]['sequence'])):
                if block['strains'][strain]['sequence'][pos] == '-':
                    gaps.append(n)
                else:
                    ungapped.append(block['strains'][strain]['sequence'][pos])
                    n += 1
            block['strains'][strain]['sequence'] = ungapped
            block['gaps'] = gaps # needs to be the same for all sequences

#def gapped(ungapped, gaps):
    
