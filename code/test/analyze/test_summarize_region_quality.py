import analyze.summarize_region_quality as summarize
from StringIO import StringIO
import pytest


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
    assert sites == []

    sites = summarize.convert_intervals_to_sites([(1, 2)])
    assert sites == [1, 2]

    sites = summarize.convert_intervals_to_sites([(1, 2), (4, 6)])
    assert sites == [1, 2, 4, 5, 6]


def test_index_alignment_by_reference():
    assert summarize.gp.gap_symbol == '-'

    output = summarize.index_alignment_by_reference('')
    assert output == []

    output = summarize.index_alignment_by_reference('abc')
    assert output == [0, 1, 2]

    output = summarize.index_alignment_by_reference('a-b-c')
    assert output == [0, 2, 4]


def test_seq_id_hmm():
    assert summarize.gp.gap_symbol == '-'
    assert summarize.gp.unsequenced_symbol == 'n'

    match, sites, d = summarize.seq_id_hmm('abd', 'abc', 0, [1, 2, 5])
    assert match == 1  # only count matches in included sites
    assert sites == 2  # included, not matching
    assert d['gap_flag'] == [False] * 3
    assert d['hmm_flag'] == [False, True, True]
    assert d['match'] == [True, True, False]
    assert d['unseq_flag'] == [False] * 3
    assert len(d) == 4

    match, sites, d = summarize.seq_id_hmm('n-d', '--c', 1, [3, 5])
    assert match == 0
    assert sites == 1
    assert d['gap_flag'] == [True, True, False]
    assert d['hmm_flag'] == [False, False, True]
    assert d['match'] == [False, True, False]
    assert d['unseq_flag'] == [True, False, False]
    assert len(d) == 4


def test_seq_id_unmasked():
    assert summarize.gp.gap_symbol == '-'
    assert summarize.gp.unsequenced_symbol == 'n'

    match, sites, d = summarize.seq_id_unmasked('abd', 'abc', 0, [], [])
    assert match == 2
    assert sites == 3

    match, sites, d = summarize.seq_id_unmasked('abd', 'abc', 0, [0], [])
    assert match == 1
    assert sites == 2
    assert d['mask_flag'] == [True, False, False]

    match, sites, d = summarize.seq_id_unmasked('abd', 'abc', 2, [0], [1])
    assert match == 2
    assert sites == 3
    assert d['mask_flag'] == [False, False, False]

    match, sites, d = summarize.seq_id_unmasked('abd', 'abc', 0, [0], [1])
    assert match == 0
    assert sites == 1
    assert d['mask_flag'] == [True, True, False]
