from misc.binary_search import present


def test_present():
    assert present([], 1) is False
    assert present([1], 1) is True
    assert present([1, 2], 1) is True
    assert present([1, 2], 2) is True
    assert present([1, 2], 3) is False
