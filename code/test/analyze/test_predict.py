import analyze.predict as predict
import hmm.hmm_bw as hmm
import pytest
from pytest import approx
# TODO: with python 3 change to io.StringIO
from StringIO import StringIO
from collections import Counter, defaultdict
import random
import numpy as np


def test_gp_symbols():
    # because the following tests use symbols instead of the variables
    assert predict.gp.match_symbol == '+'
    assert predict.gp.mismatch_symbol == '-'
    assert predict.gp.gap_symbol == '-'
    assert predict.gp.unsequenced_symbol == 'n'


@pytest.fixture
def args():
    args = predict.process_predict_args('p4e2 .001 viterbi 10000 .025 10000\
                                         .025 10000 .025 10000 .025 unknown\
                                         1000 .01'.split())
    return args


def test_process_predict_args():
    # test with default args
    args = predict.process_predict_args('p4e2 .001 viterbi 10000 .025 10000\
                                         .025 10000 .025 10000 .025 unknown\
                                         1000 .01'.split())
    assert args['tag'] == 'p4e2'
    assert args['improvement_frac'] == 0.001
    assert args['threshold'] == 'viterbi'

    assert args['known_states'] == predict.gp.alignment_ref_order
    assert args['unknown_states'] == ['unknown']
    assert args['states'] == predict.gp.alignment_ref_order + ['unknown']

    assert args['expected_frac'] == {'DBVPG6304': 0.025,
                                     'UWOPS91_917_1': 0.025,
                                     'unknown': 0.01,
                                     'CBS432': 0.025,
                                     'N_45': 0.025,
                                     'S288c': 0.89}

    assert args['expected_tract_lengths'] == {'DBVPG6304': 10000.0,
                                              'UWOPS91_917_1': 10000.0,
                                              'unknown': 1000.0,
                                              'CBS432': 10000.0,
                                              'N_45': 10000.0,
                                              'S288c': 0}
    assert args['expected_num_tracts'] == {}
    assert args['expected_bases'] == {}

    assert len(args.keys()) == 10


def test_process_predict_args_threshold():
    args = predict.process_predict_args('p4e2 .001 test 10000 .025 10000\
                                         .025 10000 .025 10000 .025 unknown\
                                         1000 .01'.split())
    assert args['threshold'] == 'viterbi'

    args = predict.process_predict_args('p4e2 .001 0.1 10000 .025 10000\
                                         .025 10000 .025 10000 .025 unknown\
                                         1000 .01'.split())
    assert args['threshold'] == 0.1


def test_process_predict_args_exceptions():
    # not enough unknown values
    with pytest.raises(IndexError):
        predict.process_predict_args('p4e2 .001 0.1 10000 .025 10000\
                                      .025 10000 .025 10000 .025 unknown\
                                      1000'.split())

    # not enough arg values
    with pytest.raises(ValueError):
        predict.process_predict_args('p4e2 .001 0.1 10000 .025 10000\
                                      .025 10000 .025 .025 unknown\
                                      1000 0.01'.split())

    with pytest.raises(ValueError):
        predict.process_predict_args('p4e2 .001 0.1 10000 .025 10000\
                                      .025 10000 .025 10000 .025 unknown\
                                      1000 NotADouble'.split())


def test_write_blocks_header():
    writer = StringIO()
    predict.write_blocks_header(writer)

    assert writer.getvalue() == '\t'.join(['strain',
                                           'chromosome',
                                           'predicted_species',
                                           'start',
                                           'end',
                                           'num_sites_hmm']) + '\n'


def test_get_emis_symbols():
    assert predict.get_emis_symbols([1]*1) == ['+', '-']
    assert predict.get_emis_symbols([1]*3) == ['+++',
                                               '++-',
                                               '+-+',
                                               '+--',
                                               '-++',
                                               '-+-',
                                               '--+',
                                               '---',
                                               ]


def test_write_hmm_header():
    # TODO this has a return value which doesn't seem to be used
    writer = StringIO()
    predict.write_hmm_header([], [], [], writer)
    assert writer.getvalue() == 'strain\tchromosome\n'

    writer = StringIO()
    predict.write_hmm_header(['s1', 's2'], ['u1'], ['-', '+'], writer)

    header = 'strain\tchromosome\t'
    header += '\t'.join(
        ['init_{}'.format(s) for s in ['s1', 's2', 'u1']] +
        ['emis_{}_{}'.format(s, sym)
         for s in ['s1', 's2', 'u1']
         for sym in ['-', '+']] +
        ['trans_{}_{}'.format(s, s2)
         for s in ['s1', 's2', 'u1']
         for s2 in ['s1', 's2', 'u1']])

    assert writer.getvalue() == header + '\n'


def test_ungap_and_code_helper():
    # nothing in prediction
    sequence, positions = predict.ungap_and_code_helper(
        '---',  # predicted reference string
        ['abc', 'def', 'ghi'],  # several references
        0)  # reference index
    assert positions == []
    assert sequence == []

    # one match
    sequence, positions = predict.ungap_and_code_helper(
        'a--',
        ['abc', 'def', 'ghi'],
        0)
    assert positions == [0]
    assert sequence == ['+--']

    # no match from refs
    sequence, positions = predict.ungap_and_code_helper(
        'a--',
        ['abc', 'def', '-hi'],
        0)
    assert positions == []
    assert sequence == []

    # two matches
    sequence, positions = predict.ungap_and_code_helper(
        'ae-',
        ['abc', 'def', 'gei'],
        0)
    assert positions == [0, 1]
    assert sequence == ['+--', '-++']

    # mess with ref index
    sequence, positions = predict.ungap_and_code_helper(
        'a--e-',
        ['a--bc', 'deeef', 'geeei'],
        0)
    assert positions == [0, 1]
    assert sequence == ['+--', '-++']
    sequence, positions = predict.ungap_and_code_helper(
        'a--e-',
        ['a--bc', 'deeef', 'geeei'],
        1)
    assert positions == [0, 3]
    assert sequence == ['+--', '-++']


def test_ungap_and_code():
    sequence, ref_seqs, positions = predict.ungap_and_code(
        'a---ef--i',
        ['ab-dhfghi',
         'a-cceeg-i',
         'a-ceef-hh'],
        0)

    assert sequence == '+++ -++ +-+ ++-'.split()
    assert ref_seqs[0] == '+++ +-- +-- +-+ ++-'.split()
    assert ref_seqs[1] == '+++ -+- -++ -+- ++-'.split()
    assert ref_seqs[2] == '+++ --+ -++ +-+ --+'.split()
    assert len(ref_seqs) == 3
    assert positions == [0, 3, 4, 7]


def test_poly_sites():
    sequence, ref_seqs, positions = predict.poly_sites(
        '+++ -++ +-+ ++-'.split(),
        ['+++ +-- +-- +-+ ++-'.split(),
            '+++ -+- -++ -+- ++-'.split(),
            '+++ --+ -++ +-+ --+'.split()],
        [0, 3, 4, 7]
    )
    assert sequence == '-++ +-+ ++-'.split()
    # TODO check if this is intended, the last element is removed from ref
    assert ref_seqs[0] == '+-- +-- +-+'.split()
    assert ref_seqs[1] == '-+- -++ -+-'.split()
    assert ref_seqs[2] == '--+ -++ +-+'.split()
    assert len(ref_seqs) == 3
    assert positions == [3, 4, 7]


def test_set_expectations_default(args):
    prev_tract = dict(args['expected_tract_lengths'])
    assert args['expected_num_tracts'] == {}
    assert args['expected_bases'] == {}
    predict.set_expectations(args, 1e5)  # made number arbitrary
    print(args)
    assert args['expected_num_tracts'] == {'DBVPG6304': 0.025 * 10,
                                           'UWOPS91_917_1': 0.025 * 10,
                                           'CBS432': 0.025 * 10,
                                           'N_45': 0.025 * 10,
                                           'S288c': 1 + 1}

    assert args['expected_bases'] == {'DBVPG6304': 0.025 * 1e5,
                                      'UWOPS91_917_1': 0.025 * 1e5,
                                      'CBS432': 0.025 * 1e5,
                                      'N_45': 0.025 * 1e5,
                                      'S288c': 1e5 - 1e4}
    prev_tract['S288c'] = 45000
    assert args['expected_tract_lengths'] == prev_tract


def test_get_symbol_freqs():
    sequence = '-++ +-+ ++- ---'.split()
    symbol_test_helper(sequence)
    # get all len 10 symbols
    syms = predict.get_emis_symbols([1]*10)

    random.seed(0)
    for i in range(10):
        sequence = [random.choice(syms) for j in range(100)]
        symbol_test_helper(sequence)


def symbol_test_helper(sequence):
    # NOTE this function assumes all transitions are same length, test does not
    ind, symb, weigh = predict.get_symbol_freqs(sequence)

    individual = []
    weighted = []

    symbols = dict(Counter(sequence))
    total = float(len(sequence))
    for k in symbols:
        symbols[k] /= total
    symbols = defaultdict(int, symbols)

    sequence = np.array(sequence)
    sequence = sequence.view('S1').reshape((sequence.size, -1))

    # look along species
    for s in np.transpose(sequence):
        s = ''.join(s)
        counts = Counter(s)
        weighted.append(counts[predict.gp.match_symbol])
        # NOTE casting unncessary in python 3
        total = float(sum(counts.values()))
        for k in counts:
            counts[k] /= total
        individual.append(defaultdict(int, counts))

    total = float(sum(weighted))
    weighted = [w / total for w in weighted]

    assert ind == individual
    assert symb == symbols
    assert weigh == approx(weighted)


def test_norm_list():
    random.seed(0)
    for i in range(10):
        test_list = [random.randint(0, 100) for j in range(100)]

        mynorm = test_list / np.sum(test_list, dtype=np.float)
        assert predict.norm_list(test_list) == approx(mynorm)


def test_initial_probabilities(args):
    probs = predict.initial_probabilities(args['known_states'],
                                          args['unknown_states'],
                                          args['expected_frac'],
                                          [0.1, 0.2, 0.3, 0.4, 0.5])

    assert args['expected_frac'] == {'DBVPG6304': 0.025,
                                     'UWOPS91_917_1': 0.025,
                                     'unknown': 0.01,
                                     'CBS432': 0.025,
                                     'N_45': 0.025,
                                     'S288c': 0.89}
    p = [0.1 + (0.89 - 0.1) * 0.9,
         0.2 + (0.025 - 0.2) * 0.9,
         0.3 + (0.025 - 0.3) * 0.9,
         0.4 + (0.025 - 0.4) * 0.9,
         0.5 + (0.025 - 0.5) * 0.9,
         0.01]

    p = p / np.sum(p, dtype=np.float)

    assert probs == approx(p)


def test_emission_probabilities(args):
    # normal mode
    symbols = predict.get_emis_symbols([1]*5)

    # NOTE not sure why this takes the keys in predict.py or uses len-2
    emis = predict.emission_probabilities(args['known_states'],
                                          args['unknown_states'],
                                          symbols)

    np_emis = np_emission(args, symbols)
    is_approx_equal_list_dict(emis, np_emis)

    # too many symbols
    symbols = predict.get_emis_symbols([1]*6)
    emis = predict.emission_probabilities(args['known_states'],
                                          args['unknown_states'],
                                          symbols)
    np_emis = np_emission(args, symbols)
    is_approx_equal_list_dict(emis, np_emis)

    # more unknowns
    args['unknown_states'].append('test')
    symbols = predict.get_emis_symbols([1]*5)
    emis = predict.emission_probabilities(args['known_states'],
                                          args['unknown_states'],
                                          symbols)
    np_emis = np_emission(args, symbols)
    is_approx_equal_list_dict(emis, np_emis)

    # no unknowns
    args['unknown_states'] = []
    symbols = predict.get_emis_symbols([1]*5)
    emis = predict.emission_probabilities(args['known_states'],
                                          args['unknown_states'],
                                          symbols)
    np_emis = np_emission(args, symbols)
    is_approx_equal_list_dict(emis, np_emis)


def np_emission(args, symbols):
    probs = {'-+': 0.9,
             '++': 0.09,
             '--': 0.009,
             '+-': 0.001}
    mismatch_bias = 0.99

    for k in probs:
        probs[k] *= 2**(len(args['known_states']) - 2)

    # make char array
    symbols = np.array(symbols)
    syms = symbols.view('S1').reshape((symbols.size, -1))
    # add first column to rest
    first_col = np.tile(syms[:, 0:1], (1, len(args['known_states'])))
    syms = np.core.defchararray.add(first_col,
                                    syms[:, 0:len(args['known_states'])])
    # lookup in table
    emit = np.vectorize(probs.__getitem__)(syms)
    # normalize
    emit /= sum(emit)

    # for unknown
    num_match = np.vectorize(
        lambda x: x.count(predict.gp.match_symbol))(symbols)
    # convert to match * (1-bias) + mismatch * bias, simplified
    num_match = num_match + mismatch_bias * (len(symbols[0]) - 2 * num_match)
    # normalize
    num_match /= sum(num_match)
    # repeat for number of unknown states
    num_match = np.transpose(
        np.tile(num_match, (len(args['unknown_states']), 1)))

    # create result
    result = [defaultdict(float, {k: v for k, v in zip(symbols, emit[:, i])})
              for i in range(emit.shape[1])]
    result.extend([defaultdict(float, {k: v for k, v in
                                       zip(symbols, num_match[:, i])})
                   for i in range(num_match.shape[1])])

    return result


def is_approx_equal_list_dict(actual, expected):
    for i in range(len(actual)):
        for k in actual[i]:
            assert actual[i][k] == approx(expected[i][k]),\
                "failed at i={}, k={}".format(i, k)


def test_norm_dict():
    def mynorm(d):
        total = float(sum(d.values()))
        return {k: v/total for k, v in d.items()}

    random.seed(0)
    for i in range(10):
        d = {k: v for k, v in zip(range(i),
                                  [random.randint(0, 100) for j in range(i)])}
        is_approx_equal_list_dict([predict.norm_dict(d)], [mynorm(d)])


def test_transition_probabilities(args):
    args['expected_tract_lengths']['S288c'] = 45000
    trans = predict.transition_probabilities(args['known_states'],
                                             args['unknown_states'],
                                             args['expected_frac'],
                                             args['expected_tract_lengths'])

    np_trans = np_transition(args)
    for i in range(len(trans)):
        assert trans[i] == approx(np_trans[i])


def np_transition(args):
    fracs = np.array([args['expected_frac'][k]
                      for k in args['known_states'] + args['unknown_states']])
    lens = np.array([1./args['expected_tract_lengths'][k]
                     for k in args['known_states'] + args['unknown_states']])

    np_trans = np.outer(np.multiply(lens, 1/(1-fracs)),
                        fracs)
    for i in range(np_trans.shape[0]):
        np_trans[i, i] = 1 - lens[i]

    # NOTE due to some floating point error (I suspect), the norm_list doesn't
    # actually alter the list
    # np_trans /= sum(np_trans)
    np_trans = [[el for el in row] for row in np_trans]
    return np_trans


def test_initial_hmm_parameters(args):
    args['expected_tract_lengths']['S288c'] = 45000
    symbols = predict.get_emis_symbols([1]*5)
    probs, emis, trans = predict.initial_hmm_parameters(
        symbols,
        args['known_states'],
        args['unknown_states'],
        args['expected_frac'],
        args['expected_tract_lengths'])

    assert args['expected_frac'] == {'DBVPG6304': 0.025,
                                     'UWOPS91_917_1': 0.025,
                                     'unknown': 0.01,
                                     'CBS432': 0.025,
                                     'N_45': 0.025,
                                     'S288c': 0.89}
    p = [0.2 + (0.89 - 0.2) * 0.9,
         0.2 + (0.025 - 0.2) * 0.9,
         0.2 + (0.025 - 0.2) * 0.9,
         0.2 + (0.025 - 0.2) * 0.9,
         0.2 + (0.025 - 0.2) * 0.9,
         0.01]

    p = p / np.sum(p, dtype=np.float)
    assert probs == approx(p)

    np_emis = np_emission(args, symbols)
    is_approx_equal_list_dict(emis, np_emis)

    np_trans = np_transition(args)
    for i in range(len(trans)):
        assert trans[i] == approx(np_trans[i])


def test_predict_introgressed(args, capsys):
    seqs = [list('NNENNENNEN'),  # S2288c
            list('NNNENEENNN'),  # CBS432
            list('NN-NNEENNN'),  # N_45
            list('NEENN-ENEN'),  # DBVPG6304
            list('ENENNEENEN'),  # UWOPS..
            list('NNENNEENEN'),  # predicted
            ]
    ref = seqs[:-1]
    pred = seqs[-1]

    path, prob, hmm, hmm_init, ps = predict.predict_introgressed(
        ref, pred, args, train=True)

    # check hmm output
    captured = capsys.readouterr()
    out = captured.out.split('\n')
    assert 'finished in 10 iterations' in out[-2]

    # ps are locations of polymorphic sites, not counting missing '-'
    assert ps == [0, 1, 3, 6, 8]
    assert hmm.init == [1, 0, 0, 0, 0, 0]

    # check path
    assert path == ['S288c', 'S288c', 'UWOPS91_917_1',
                    'UWOPS91_917_1', 'UWOPS91_917_1']

    assert prob[0]['S288c'] == 1


def test_write_positions():
    output = StringIO()
    predict.write_positions([0, 1, 3, 5, 7], output, 'test', 'I')
    assert output.getvalue() == "{}\t{}\t{}\n".format(
        "test",
        "I",
        "\t".join([str(i) for i in (0, 1, 3, 5, 7)]))


def test_write_blocks():
    output = StringIO()
    block = [(0, 1), (4, 6), (10, 8)]
    pos = [i * 2 for i in range(20)]
    predict.write_blocks(block,
                         pos,
                         output, 'test', 'I', 'pred')

    result = "\n".join(
        ["\t".join(['test', 'I', 'pred',
                    str(pos[s]), str(pos[e]), str(e - s + 1)])
         for s, e in block]) + "\n"

    assert output.getvalue() == result


def test_read_blocks(mocker):
    block_in = StringIO('''
''')

    mocked_file = mocker.patch('analyze.predict.open',
                               return_value=block_in)
    output = predict.read_blocks('mocked')

    mocked_file.assert_called_with('mocked', 'r')
    assert output.keys() == []

    block_in = StringIO('''header
test\tI\tpred\t100\t200\t10
''')

    mocked_file = mocker.patch('analyze.predict.open',
                               return_value=block_in)
    output = predict.read_blocks('mocked')

    assert len(output) == 1
    assert output['test']['I'] == [(100, 200, 10)]

    block_in = StringIO('''header
test\tI\tpred\t100\t200\t10
test\tI\tpred\t200\t200\t30
test\tI\tpred\t300\t400\t40
test\tII\tpred\t300\t400\t40
test2\tIII\tpred\t300\t400\t47
''')

    mocked_file = mocker.patch('analyze.predict.open',
                               return_value=block_in)
    output = predict.read_blocks('mocked')

    assert len(output) == 2
    assert len(output['test']) == 2
    assert len(output['test2']) == 1
    assert output['test']['I'] == [
        (100, 200, 10),
        (200, 200, 30),
        (300, 400, 40),
    ]
    assert output['test']['II'] == [(300, 400, 40)]
    assert output['test2']['III'] == [(300, 400, 47)]


def test_read_blocks_labeled(mocker):
    block_in = StringIO('''
''')

    mocked_file = mocker.patch('analyze.predict.open',
                               return_value=block_in)
    output = predict.read_blocks('mocked', labeled=True)

    mocked_file.assert_called_with('mocked', 'r')
    assert output.keys() == []

    block_in = StringIO('''header
r1\ttest\tI\tpred\t100\t200\t10
''')

    mocked_file = mocker.patch('analyze.predict.open',
                               return_value=block_in)
    output = predict.read_blocks('mocked', labeled=True)

    assert len(output) == 1
    assert output['test']['I'] == [('r1', 100, 200, 10)]

    block_in = StringIO('''header
r1\ttest\tI\tpred\t100\t200\t10
r2\ttest\tI\tpred\t200\t200\t30
r3\ttest\tI\tpred\t300\t400\t40
r4\ttest\tII\tpred\t300\t400\t40
r5\ttest2\tIII\tpred\t300\t400\t47
''')

    mocked_file = mocker.patch('analyze.predict.open',
                               return_value=block_in)
    output = predict.read_blocks('mocked', labeled=True)

    assert len(output) == 2
    assert len(output['test']) == 2
    assert len(output['test2']) == 1
    assert output['test']['I'] == [
        ('r1', 100, 200, 10),
        ('r2', 200, 200, 30),
        ('r3', 300, 400, 40),
    ]
    assert output['test']['II'] == [('r4', 300, 400, 40)]
    assert output['test2']['III'] == [('r5', 300, 400, 47)]


def test_write_hmm():
    output = StringIO()

    hm = hmm.HMM()

    # empty hmm
    predict.write_hmm(hm, output, 'strain', 'I', list('abc'))
    assert output.getvalue() == 'strain\tI\n'

    hm.set_states(list('abc'))
    hm.set_init([0, 1, 0])
    hm.set_trans([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
    hm.set_emis([{'a': 1, 'b': 0, 'c': 0},
                 {'a': 0, 'b': 0, 'c': 1},
                 {'a': 0, 'b': 1, 'c': 0},
                 ])

    output = StringIO()
    predict.write_hmm(hm, output, 'strain', 'I', list('abc'))

    result = 'strain\tI\t'
    result += '\t'.join(list('010')) + '\t'  # init
    result += '\t'.join(list('100001010')) + '\t'  # emis
    result += '\t'.join(list('010100001')) + '\n'  # trans
    assert output.getvalue() == result


def test_write_state_probs():
    output = StringIO()
    predict.write_state_probs([{}], output, 'strain', 'I')

    assert output.getvalue() == 'strain\tI\n'

    output = StringIO()
    predict.write_state_probs([
        {'a': 0, 'b': 0, 'c': 1},
        {'a': 1, 'b': 0, 'c': 0},
        {'a': 0, 'b': 1, 'c': 1},
    ], output, 'strain', 'I')

    assert output.getvalue() == \
        ('strain\tI\t'
         'a:0.00000,1.00000,0.00000\t'
         'c:1.00000,0.00000,1.00000\t'
         'b:0.00000,0.00000,1.00000\n')
