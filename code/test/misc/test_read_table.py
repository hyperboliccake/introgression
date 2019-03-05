import misc.read_table as read_table
from io import StringIO
import pytest


def test_read_table_rows_empty(mocker):
    table = StringIO('')
    mocked_file = mocker.patch('misc.read_table.open', return_value=table)
    mocked_gz = mocker.patch('misc.read_table.gzip.open', return_value=table)

    d, labels = read_table.read_table_rows('mocked', '\t', False)
    mocked_file.assert_called_with('mocked', 'r')
    mocked_gz.assert_not_called()
    assert d == {}
    assert labels is None

    table = StringIO('')
    mocked_file = mocker.patch('misc.read_table.open', return_value=table)
    mocked_gz = mocker.patch('misc.read_table.gzip.open', return_value=table)

    d, labels = read_table.read_table_rows('mocked.gz', '\t', False)
    mocked_gz.assert_called_with('mocked.gz', 'rb')
    mocked_file.assert_not_called()
    assert d == {}
    assert labels is None


def test_read_table_rows_no_header(mocker):
    table = StringIO('a,b,c,d\ne,f,g\nh,i,j,k\nl,m,n,o,p\n')
    mocked_file = mocker.patch('misc.read_table.open', return_value=table)
    d, labels = read_table.read_table_rows('mocked', ',', False, 0)

    mocked_file.assert_called_with('mocked', 'r')
    assert d == {'a': list('bcd'), 'e': list('fg'),
                 'h': list('ijk'), 'l': list('mnop')}
    assert labels is None

    table = StringIO('a,b,c,d\ne,f,g\nh,i,j,k\nl,m,n,o,p\n')
    mocked_file = mocker.patch('misc.read_table.open', return_value=table)
    d, labels = read_table.read_table_rows('mocked', ',', False, 2)

    mocked_file.assert_called_with('mocked', 'r')
    assert d == {'c': list('abd'), 'g': list('ef'),
                 'j': list('hik'), 'n': list('lmop')}
    assert labels is None

    table = StringIO('a\n\nh,i,j,k\nl,m,n,o,p\n')
    mocked_file = mocker.patch('misc.read_table.open', return_value=table)
    d, labels = read_table.read_table_rows('mocked', ',', False, 0)

    mocked_file.assert_called_with('mocked', 'r')
    assert d == {'a': [], '': [],
                 'h': list('ijk'), 'l': list('mnop')}
    assert labels is None


def test_read_table_rows_with_header(mocker):
    table = StringIO('a,b,c,d\ne,f,g,h\ni,j,k,l\nm,n,o,p\n')
    mocked_file = mocker.patch('misc.read_table.open', return_value=table)
    d, labels = read_table.read_table_rows('mocked', ',', True, 0)

    mocked_file.assert_called_with('mocked', 'r')
    assert d == {'e': {k: v for k, v in zip(list('bcd'), list('fgh'))},
                 'i': {k: v for k, v in zip(list('bcd'), list('jkl'))},
                 'm': {k: v for k, v in zip(list('bcd'), list('nop'))},
                 }
    assert labels == list('abcd')

    # fewer headers than rest
    table = StringIO('a,b,c\ne,f,g,h\ni,j,k,l\nm,n,o,p\n')
    mocked_file = mocker.patch('misc.read_table.open', return_value=table)
    d, labels = read_table.read_table_rows('mocked', ',', True, 1)

    mocked_file.assert_called_with('mocked', 'r')
    assert d == {'f': {k: v for k, v in zip(list('ac'), list('eg'))},
                 'j': {k: v for k, v in zip(list('ac'), list('ik'))},
                 'n': {k: v for k, v in zip(list('ac'), list('mo'))},
                 }
    assert labels == list('abc')

    # fewer columns than headers
    table = StringIO('a,b,c,d\ne,f,g\ni,j,k\nm,n,o\n')
    mocked_file = mocker.patch('misc.read_table.open', return_value=table)
    d, labels = read_table.read_table_rows('mocked', ',', True, 0)

    mocked_file.assert_called_with('mocked', 'r')
    assert d == {'e': {k: v for k, v in zip(list('bc'), list('fg'))},
                 'i': {k: v for k, v in zip(list('bc'), list('jk'))},
                 'm': {k: v for k, v in zip(list('bc'), list('no'))},
                 }
    assert labels == list('abcd')


def test_read_table_columns_empty(mocker):
    table = StringIO('')
    mocked_file = mocker.patch('misc.read_table.open', return_value=table)
    mocked_gz = mocker.patch('misc.read_table.gzip.open', return_value=table)

    d, labels = read_table.read_table_columns('mocked', '\t')
    mocked_file.assert_called_with('mocked', 'r')
    mocked_gz.assert_not_called()
    assert d == {'': []}
    assert labels == ['']

    table = StringIO('')
    mocked_file = mocker.patch('misc.read_table.open', return_value=table)
    mocked_gz = mocker.patch('misc.read_table.gzip.open', return_value=table)

    d, labels = read_table.read_table_columns('mocked.gz', '\t')
    mocked_gz.assert_called_with('mocked.gz', 'rb')
    mocked_file.assert_not_called()
    assert d == {'': []}
    assert labels == ['']


def test_read_table_columns(mocker):
    table = StringIO('a,b,c,d\n'
                     'e,f,g,h\n'
                     'i,j,k,l\n'
                     'm,n,o,p\n')
    mocked_file = mocker.patch('misc.read_table.open', return_value=table)
    d, labels = read_table.read_table_columns('mocked', ',')

    mocked_file.assert_called_with('mocked', 'r')
    assert d == {'a': list('eim'),
                 'b': list('fjn'),
                 'c': list('gko'),
                 'd': list('hlp')
                 }
    assert labels == list('abcd')

    # fewer headers than rest
    table = StringIO('a,b,c\ne,f,g,h\ni,j,k,l\nm,n,o,p\n')
    mocked_file = mocker.patch('misc.read_table.open', return_value=table)
    d, labels = read_table.read_table_columns('mocked', ',')

    mocked_file.assert_called_with('mocked', 'r')
    assert d == {'a': list('eim'),
                 'b': list('fjn'),
                 'c': list('gko'),
                 }
    assert labels == list('abc')

    # fewer columns than headers
    table = StringIO('a,b,c,d\ne,f,g\ni,j,k\nm,n,o\n')
    mocked_file = mocker.patch('misc.read_table.open', return_value=table)

    with pytest.raises(IndexError):
        d, labels = read_table.read_table_columns('mocked', ',')
    mocked_file.assert_called_with('mocked', 'r')
