from sim import sim_process
from hmm import hmm_bw as hmm
import pytest
import operator
import numpy as np
from collections import defaultdict
import random


@pytest.fixture
def hm():
    hm = hmm.HMM()
    hm.set_states(['N', 'E'])
    hm.set_emissions([{'N': 0.3, 'E': 0.7},
                      {'N': 0.8, 'E': 0.2}])
    hm.set_observations([list('NNENNENNEN'),
                         list('NNNNNEENNN'),
                         list('NNENNEENEN')])
    hm.set_transitions([[0.5, 0.5], [0.3, 0.7]])
    hm.set_initial_p([0.2, 0.8])

    return hm


def test_get_max_path(hm):
    post = hm.posterior_decoding()
    path, probs = sim_process.get_max_path(post[0])
    maxes = [max(d.items(), key=operator.itemgetter(1)) for d in post[0]]
    other_path, other_prob = map(list, zip(*maxes))
    assert path == other_path
    assert probs == other_prob


def test_get_threshold_predicted(hm):
    post = hm.posterior_decoding()
    path, probs = sim_process.get_max_path(post[0])
    for thresh in (0, 0.2, 0.5, 0.8, 1):
        pred = sim_process.threshold_predicted(path, probs, thresh, 'N')
        t = np.array(path)
        p = np.array(probs)
        t[p < thresh] = 'N'
        assert pred == list(t)


def test_convert_to_blocks_one():
    random.seed(0)
    states = [str(i) for i in range(10)]
    # TODO fix this as an error or throw another exception
    # help_test_convert_blocks(states, [])
    help_test_convert_blocks(states, list('1'))
    help_test_convert_blocks(states, list('12'))
    help_test_convert_blocks(states, list('1111'))

    for test in range(10):
        seq = [str(random.randint(0, 9)) for i in range(100)]
        help_test_convert_blocks(states, seq)


def help_test_convert_blocks(states, seq):
    blocks = sim_process.convert_to_blocks_one(seq, states)

    nseq = np.array(seq, int)
    # add element to the end to catch repeats on last index
    nseq = np.append(nseq, nseq[-1]+1)
    diff = np.diff(nseq)
    locs = np.nonzero(diff)[0]
    lens = np.diff(locs)
    lens = np.append(locs[0]+1, lens)

    current = 0
    result = defaultdict(list)
    for i, l in enumerate(locs):
        result[seq[l]].append((current, current + lens[i] - 1))
        current += lens[i]

    for k in blocks:
        assert blocks[k] == result[k]
