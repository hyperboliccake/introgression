import misc.seq_functions as seq


def test_seq_id():
    match, sites = seq.seq_id([], [])
    assert match == 0
    assert sites == 0

    match, sites = seq.seq_id('atcg', 'tcgafd')
    assert match == 0
    assert sites == 4

    match, sites = seq.seq_id('axat', 'txtt')
    assert match == 1
    assert sites == 3

    match, sites = seq.seq_id('axxt', 'xatx')
    assert match == 0
    assert sites == 0
