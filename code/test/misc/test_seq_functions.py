import misc.seq_functions as seq
import numpy as np


def test_seq_id():
    match, sites = seq.seq_id(np.array([]), np.array([]))
    assert match == 0
    assert sites == 0

    match, sites = seq.seq_id(np.array(list('atcg')), np.array(list('tcgafd')))
    assert match == 0
    assert sites == 4

    match, sites = seq.seq_id(np.array(list('axat')), np.array(list('txtt')))
    assert match == 1
    assert sites == 3

    match, sites = seq.seq_id(np.array(list('axxt')), np.array(list('xatx')))
    assert match == 0
    assert sites == 0
