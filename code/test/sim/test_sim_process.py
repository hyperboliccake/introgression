from sim import sim_process
import numpy as np
from collections import defaultdict
import random


def test_get_max_path(hm):
    post = hm.posterior_decoding()[0]
    path, probs = sim_process.get_max_path(post, hm.hidden_states)

    max_path = []
    max_probs = []
    for site_probs in post:
        max_state = None
        max_prob = -1
        for i, prob in enumerate(site_probs):
            if prob > max_prob:
                max_prob = prob
                max_state = hm.hidden_states[i]
        max_path.append(max_state)
        max_probs.append(max_prob)

    assert path == max_path
    assert probs == max_probs


def test_get_threshold_predicted(hm):
    post = hm.posterior_decoding()
    path, probs = sim_process.get_max_path(post[0], hm.hidden_states)
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
