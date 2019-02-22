import align.align_helpers as helper
import pytest


def test_flatten():
    assert helper.flatten([]) == []
    assert helper.flatten([[]]) == []
    assert helper.flatten([[1]]) == [1]
    assert helper.flatten([[2], [1]]) == [2, 1]
    assert helper.flatten([[1, 2], [1]]) == [1, 2, 1]


def test_get_strains(mocker, capsys):
    # ensure params are correct
    mocker.patch('align.align_helpers.gp.fasta_suffix',
                 '.fa')
    mocker.patch('align.align_helpers.gp.chrms',
                 ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX',
                  'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI'])

    mock_dir = mocker.patch('os.listdir')

    mock_dir.return_value = []
    helper.get_strains(['mock'])
    assert "found no chromosome sequence files in mock" in\
        capsys.readouterr().out

    mock_dir.return_value = ['nothing']
    helper.get_strains(['invalid'])
    assert "found no chromosome sequence files in invalid" in\
        capsys.readouterr().out

    with pytest.raises(AssertionError) as e:
        mock_dir.return_value = ['missing_chr.fa']
        helper.get_strains(['noStrains'])

    assert "some strains in noStrains are missing" in str(e)

    mock_dir.return_value = ['strain_chr{}.fa'.format(chrm)
                             for chrm in helper.gp.chrms]
    strains = helper.get_strains(['one_strain'])
    assert strains == [('strain', 'one_strain')]

    mock_dir.return_value = ['strain{}_chr{}.fa'.format(strain, chrm)
                             for chrm in helper.gp.chrms
                             for strain in (1, 2)]
    strains = helper.get_strains(['two/strains'])
    assert strains == [('strain1', 'two/strains'),
                       ('strain2', 'two/strains')]

    mock_dir.side_effect = [['strain{}_chr{}.fa'.format(strain, chrm)
                             for chrm in helper.gp.chrms
                             for strain in (1, 2)],
                            ['strain_chr{}.fa'.format(chrm)
                             for chrm in helper.gp.chrms]]
    strains = helper.get_strains(['two_strains', 'one_strain'])
    assert strains == [('strain', 'one_strain'),
                       ('strain1', 'two_strains'),
                       ('strain2', 'two_strains')]
