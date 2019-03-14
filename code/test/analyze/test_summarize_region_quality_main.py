import analyze.summarize_region_quality_main as main
from io import StringIO


def test_main(mocker):
    # setup global params to match expectations
    mocker.patch(
        'analyze.summarize_region_quality_main.gp.analysis_out_dir_absolute',
        'dir/')
    mocker.patch('analyze.summarize_region_quality_main.gp.chrms',
                 ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX',
                  'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI'])

    mocker.patch('sys.argv',
                 "test.py 5 tag .001 viterbi 10000 .025 10000 .025 \
                 10000 .025 10000 .025 unknown 1000 .01".split())
    mocker.patch('analyze.summarize_region_quality_main.os.path.isdir',
                 return_value=True)
    # TODO check call arguments
    mocker.patch('misc.read_table.read_table_columns',
                 return_value=({'s': {'region_id': []}}, ['region_id']))
    mocker.patch('analyze.summarize_region_quality_main.read_masked_intervals',
                 return_value=[(1, 2)])
    lines = StringIO('')
    mocked_file = mocker.patch(
        'analyze.summarize_region_quality_main.gzip.open',
        return_value=lines)

    mocked_file = mocker.patch('analyze.summarize_region_quality_main.open',
                               mocker.mock_open())

    main.main()

    assert mocked_file.call_count == 2
    mocked_file.assert_any_call(
        'dir/tag/blocks_unknown_tag_quality.txt', 'w')
    mocked_file.assert_any_call(
        'dir/tag/regions/unknown.pkl', 'wb')

    # just headers
    states = ['S288c', 'CBS432', 'N_45', 'DBVPG6304', 'UWOPS91_917_1']
    symbols = list('.-_npbcxNPBCX')
    mocked_file().write.assert_has_calls([
        mocker.call('\t'.join(['region_id'] +
                              ['match_nongap_' + x for x in states] +
                              ['num_sites_nongap_' + x for x in states] +
                              ['match_hmm_' + x for x in states] +
                              ['match_nonmask_' + x for x in states] +
                              ['num_sites_nonmask_' + x for x in states] +
                              ['count_' + x for x in symbols]
                              )
                    + '\n')
    ])
    assert True
