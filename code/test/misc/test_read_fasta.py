from misc.read_fasta import read_fasta
from io import StringIO
import pytest
from pytest import approx
import numpy as np


def test_read_fasta_empty(mocker):
    fasta = StringIO('>\n')
    mocked_file = mocker.patch('misc.read_fasta.open', return_value=fasta)
    headers, seqs = read_fasta('mocked')

    assert headers == ['>']
    assert seqs.tolist() == [[]]
    mocked_file.assert_called_with('mocked', 'r')

    fasta = StringIO('')
    mocked_file = mocker.patch('misc.read_fasta.open', return_value=fasta)
    # TODO probably handle empty files better
    with pytest.raises(IndexError):
        headers, seqs = read_fasta('mocked')


def test_read_fasta_single(mocker):
    fasta = StringIO('''
not read
> headseq headfname.fa
actg
---
atcg
''')
    mocked_file = mocker.patch('misc.read_fasta.open', return_value=fasta)
    headers, seqs = read_fasta('mocked')

    assert headers == ['> headseq headfname.fa']
    assert seqs == approx(np.asarray([list('actg---atcg')]))
    mocked_file.assert_called_with('mocked', 'r')


def test_read_fasta_multi(mocker):
    fasta = StringIO('''
not read
> headseq headfname.fa
actg
---
atcg
> headseq headfname.fa
actg
actg
---
> headseq2 headfname.fa
actg-
cataaa
''')
    mocked_file = mocker.patch('misc.read_fasta.open', return_value=fasta)
    headers, seqs = read_fasta('mocked')

    assert headers == ['> headseq headfname.fa',
                       '> headseq headfname.fa',
                       '> headseq2 headfname.fa']
    print(seqs)
    assert seqs == approx(np.array([list('actg---atcg'),
                                    list('actgactg---'),
                                    list('actg-cataaa')]))
    mocked_file.assert_called_with('mocked', 'r')
