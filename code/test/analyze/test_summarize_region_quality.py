import analyze.summarize_region_quality as summarize
from io import StringIO
import pytest
from pytest import approx
import numpy as np


def test_read_masked_intervals(mocker):
    lines = StringIO('')
    mocked_file = mocker.patch('analyze.summarize_region_quality.open',
                               return_value=lines)
    intervals = summarize.read_masked_intervals('mocked')
    mocked_file.assert_called_with('mocked', 'r')
    assert intervals == []

    lines = StringIO('I am a header')
    mocked_file = mocker.patch('analyze.summarize_region_quality.open',
                               return_value=lines)
    intervals = summarize.read_masked_intervals('mocked')
    mocked_file.assert_called_with('mocked', 'r')
    assert intervals == []

    lines = StringIO('I am a header\n'
                     'short and stout')
    mocked_file = mocker.patch('analyze.summarize_region_quality.open',
                               return_value=lines)
    with pytest.raises(ValueError):
        intervals = summarize.read_masked_intervals('mocked')
    mocked_file.assert_called_with('mocked', 'r')

    lines = StringIO('I am a header\n'
                     '1 and 2')
    mocked_file = mocker.patch('analyze.summarize_region_quality.open',
                               return_value=lines)
    intervals = summarize.read_masked_intervals('mocked')
    assert intervals == [(1, 2)]


def test_convert_intervals_to_sites():
    sites = summarize.convert_intervals_to_sites([])
    assert sites == approx([])

    sites = summarize.convert_intervals_to_sites([(1, 2)])
    assert sites == approx([1, 2])

    sites = summarize.convert_intervals_to_sites([(1, 2), (4, 6)])
    assert sites == approx([1, 2, 4, 5, 6])


def test_index_alignment_by_reference():
    assert summarize.gp.gap_symbol == '-'

    output = summarize.index_alignment_by_reference(np.array(list('abc')))
    assert output == approx([0, 1, 2])

    output = summarize.index_alignment_by_reference(np.array(list('a-b-c')))
    assert output == approx([0, 2, 4])


def test_seq_id_hmm():
    assert summarize.gp.gap_symbol == '-'
    assert summarize.gp.unsequenced_symbol == 'n'

    match, sites, d = summarize.seq_id_hmm(np.array(list('abd')),
                                           np.array(list('abc')),
                                           0, [1, 2, 5])
    assert match == 1  # only count matches in included sites
    assert sites == 2  # included, not matching
    assert d['gap_flag'] == approx([False] * 3)
    assert d['hmm_flag'] == approx([False, True, True])
    assert d['match'] == approx([True, True, False])
    assert d['unseq_flag'] == approx([False] * 3)
    assert len(d) == 4

    match, sites, d = summarize.seq_id_hmm(np.array(list('n-d')),
                                           np.array(list('--c')),
                                           1, [3, 5])
    assert match == 0
    assert sites == 1
    assert d['gap_flag'] == approx([True, True, False])
    assert d['hmm_flag'] == approx([False, False, True])
    assert d['match'] == approx([False, True, False])
    assert d['unseq_flag'] == approx([True, False, False])
    assert len(d) == 4

    with pytest.raises(AssertionError) as e:
        match, sites, d = summarize.seq_id_hmm(np.array(list('n-d')),
                                               np.array(list('--c')),
                                               1, [2, 5])
    assert '- - 1' in str(e)

    with pytest.raises(AssertionError) as e:
        match, sites, d = summarize.seq_id_hmm(np.array(list('n-d')),
                                               np.array(list('--c')),
                                               1, [1, 5])
    assert 'n - 0' in str(e)


def test_seq_id_unmasked():
    assert summarize.gp.gap_symbol == '-'
    assert summarize.gp.unsequenced_symbol == 'n'

    match, sites, d = summarize.seq_id_unmasked(np.array(list('abd')),
                                                np.array(list('abc')),
                                                0, [], [])
    assert match == 2
    assert sites == 3
    assert d['mask_flag'] == approx([False, False, False])

    match, sites, d = summarize.seq_id_unmasked(np.array(list('abd')),
                                                np.array(list('abc')),
                                                0, [0], [])
    assert match == 1
    assert sites == 2
    assert d['mask_flag'] == approx([True, False, False])

    match, sites, d = summarize.seq_id_unmasked(np.array(list('abd')),
                                                np.array(list('abc')),
                                                2, [0], [1])
    assert match == 2
    assert sites == 3
    assert d['mask_flag'] == approx([False, False, False])

    match, sites, d = summarize.seq_id_unmasked(np.array(list('abd')),
                                                np.array(list('abc')),
                                                0, [0], [1])
    assert match == 0
    assert sites == 1
    assert d['mask_flag'] == approx([True, True, False])


def test_make_info_string():
    len_seqx = 14
    len_states = 3
    info = {'hmm_flag': np.zeros((len_seqx), bool),
            'gap_flag': np.zeros((len_seqx, len_states), bool),
            'mask_flag': np.zeros((len_seqx, len_states), bool),
            'match_flag': np.zeros((len_seqx, len_states), bool)}

    info['gap_flag'][0, 0] = True  # -
    info['gap_flag'][11, 1] = True  # -
    info['mask_flag'][1, 0] = True  # _
    info['mask_flag'][12, 1] = True  # _
    info['match_flag'][2, :] = True  # .
    info['match_flag'][13, :] = True  # .
    info['match_flag'][(3, 4, 5, 6), 0] = True  # b and c
    info['match_flag'][(3, 4, 7, 8), 1] = True  # b and p
    # x is default
    info['hmm_flag'][[4, 6, 8, 10]] = True  # capitalize

    s = summarize.make_info_string(info, master_ind=0, predict_ind=1)
    assert s == '-_.bBcCpPxX-_.'
    #            01234567890123

    len_seqx = 0
    len_states = 3
    info = {'hmm_flag': np.zeros((len_seqx), bool),
            'gap_flag': np.zeros((len_seqx, len_states), bool),
            'mask_flag': np.zeros((len_seqx, len_states), bool),
            'match_flag': np.zeros((len_seqx, len_states), bool)}

    s = summarize.make_info_string(info, master_ind=0, predict_ind=1)
    assert s == ''


def test_info_string_unknown():
    len_seqx = 5
    len_states = 2
    info = {'gap_any_flag': np.zeros((len_seqx), bool),
            'mask_any_flag': np.zeros((len_seqx), bool),
            'match_flag': np.zeros((len_seqx, len_states), bool)}

    info['gap_any_flag'][0] = True  # -
    info['mask_any_flag'][1] = True  # _
    info['match_flag'][2, :] = True  # .
    info['match_flag'][3, 0] = True  # x
    info['match_flag'][4, 1] = True  # X

    s = summarize.make_info_string_unknown(info, master_ind=0)
    assert s == '-_.xX'
    s = summarize.make_info_string(info, master_ind=0, predict_ind=3)
    assert s == '-_.xX'

    len_seqx = 0
    len_states = 2
    info = {'gap_any_flag': np.zeros((len_seqx), bool),
            'mask_any_flag': np.zeros((len_seqx), bool),
            'match_flag': np.zeros((len_seqx, len_states), bool)}

    s = summarize.make_info_string_unknown(info, master_ind=0)
    assert s == ''
