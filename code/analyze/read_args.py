import sys
from align import align_helpers

def process_predict_args(args):
    
    d = {}
    
    i = 0
    d['tag'] = args[i]

    i += 1
    fn_setup = args[i]
    setup_args = read_setup_args(fn_setup)

    i += 1
    d['improvement_frac'] = float(args[i])
    
    i += 1
    d['threshold'] = args[i]
    if d['threshold'] != 'viterbi':
        d['threshold'] = float(d['threshold'])

    i += 1
    d['hmm_init_option'] = args[i]

    d['known_states'] = setup_args['references']
    d['unknown_states'] = []

    if d['hmm_init_option'] == 'estimate_hmm':
        d['expected_length'] = {}
        d['expected_frac'] = {}
        for state in d['known_states'][1:]:
            i += 1
            d['expected_length'][state] = float(args[i])
            i += 1
            d['expected_frac'][state] = float(args[i])
        d['expected_frac'][d['known_states'][0]] = 1 - sum(d['expected_frac'].values())
        d['expected_length'][d['known_states'][0]] = 0 # calculate later

        i += 1
        while i < len(args):
            state = args[i]
            d['unknown_states'].append(state)
            i += 1
            d['expected_length'][state] = float(args[i])
            i += 1
            d['expected_frac'][state] = float(args[i])
            i += 1

    else:
        assert d['hmm_init_option'] == 'provided_hmm'
        print('not implemented yet')
        sys.exit()

    d['states'] = d['known_states'] + d['unknown_states']

    d['setup_args'] = setup_args

    return d
    
def read_setup_args(fn):

    x = {}

    f = open(fn, 'r')
    line = f.readline()
    while line != '':
        line = line[:-1].split(' ')
        x[line[0]] = line[1:]
        line = f.readline()
    f.close()

    d = {}
    d['references'] = x['references']
    d['reference_directories'] =  dict(zip(x['references'], x['reference_directories']))
    d['alignments_directory'] = x['alignments_directory'][0]

    d['strain_dirs'] = \
        align_helpers.get_strains(x['test_strain_directories'])

    return d

def get_predict_args_by_tag(fn, tag):
    f = open(fn, 'r')
    line = f.readline()
    while line != '':
        line = line.split(' ')
        if line[0] == tag:
            return line
        line = f.readline()
    print(f'tag not found: {tag}')
    return None

