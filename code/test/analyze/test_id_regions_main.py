import analyze.id_regions_main as main
from operator import itemgetter


def test_main_blank(mocker):
    # setup global params to match expectations
    mocker.patch('analyze.predict.gp.alignment_ref_order',
                 ['ref', 'state1'])
    mocker.patch('analyze.id_regions_main.gp.chrms',
                 ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX',
                  'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI'])
    mocker.patch('analyze.id_regions_main.gp.analysis_out_dir_absolute',
                 'dir/')

    mocker.patch('sys.argv',
                 "test.py tag .001 viterbi 1000 .025 unknown 1000 .01".split())
    mocker.patch('analyze.predict.read_blocks',
                 return_value={})

    mocked_file = mocker.patch('analyze.id_regions_main.open',
                               mocker.mock_open())

    main.main()

    assert mocked_file.call_count == 3
    mocked_file.assert_any_call('dir/tag/blocks_ref_tag_labeled.txt', 'w')
    mocked_file.assert_any_call('dir/tag/blocks_state1_tag_labeled.txt', 'w')
    mocked_file.assert_any_call('dir/tag/blocks_unknown_tag_labeled.txt', 'w')

    # just headers
    mocked_file().write.assert_has_calls([
        mocker.call('region_id\tstrain\tchromosome\tpredicted_species'
                    '\tstart\tend\tnum_sites_hmm\n')
    ]*3)


def test_main(mocker):
    # setup global params to match expectations
    mocker.patch('analyze.predict.gp.alignment_ref_order',
                 ['ref', 'state1'])
    mocker.patch('analyze.id_regions_main.gp.chrms',
                 ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX',
                  'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI'])
    mocker.patch('analyze.id_regions_main.gp.analysis_out_dir_absolute',
                 'dir/')

    mocker.patch('sys.argv',
                 "test.py tag .001 viterbi 1000 .025 unknown 1000 .01".split())

    regions = {
        'ref': {
            'I': [(10, 100, 10), (10, 100, 1)],
            'IX': [(10, 100, 10), (10, 100, 1)],
            'VI': [(10, 100, 10), (10, 100, 1)],
        },
        'state1': {
            'II': [(10, 100, 10), (10, 100, 1)],
            'X': [(10, 100, 10), (10, 100, 1)],
            'V': [(10, 100, 10), (10, 100, 1)],
        },
        'unknown': {
            'II': [(10, 100, 10), (10, 100, 1)],
            'X': [(10, 100, 10), (10, 100, 1)],
            'V': [(10, 100, 10), (10, 100, 1)],
        }
    }
    mocker.patch('analyze.predict.read_blocks',
                 return_value=regions)

    mocked_file = mocker.patch('analyze.id_regions_main.open',
                               mocker.mock_open())

    main.main()

    assert mocked_file.call_count == 3
    mocked_file.assert_any_call('dir/tag/blocks_ref_tag_labeled.txt', 'w')
    mocked_file.assert_any_call('dir/tag/blocks_state1_tag_labeled.txt', 'w')
    mocked_file.assert_any_call('dir/tag/blocks_unknown_tag_labeled.txt', 'w')

    # headers
    calls = [
        mocker.call('region_id\tstrain\tchromosome\tpredicted_species'
                    '\tstart\tend\tnum_sites_hmm\n')
    ]*3

    rid = 1
    by_chrom = dict(zip(main.gp.chrms, [[] for chrm in main.gp.chrms]))
    for spec in ('ref', 'state1', 'unknown'):
        for s in sorted(regions):
            for c in regions[s]:
                for e in regions[s][c]:
                    start, end, num = e
                    by_chrom[c].append((start, end, num, s, spec))

    for c in main.gp.chrms:
        for e in sorted(by_chrom[c], key=itemgetter(0, 3)):
            start, end, num, s, spec = e
            calls.append(mocker.call(
                'r{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    rid, s, c, spec, start, end, num)))
            rid += 1

    mocked_file().write.assert_has_calls(calls)
