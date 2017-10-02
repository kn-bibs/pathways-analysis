import pytest

from command_line.types import Slice
from command_line.types import Indices
from command_line.types import positive_int


def test_slice():
    items = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    cases = {
        '2:3': [2],
        '2:5': [2, 3, 4],
        ':5': [0, 1, 2, 3, 4],
        ':': items,
        ':-2': items[:-2]
    }

    for constructor, result in cases.items():
        tested_slice = Slice(constructor)
        assert tested_slice.get(items) == result


def test_indices():
    items = ['a', 'b', 'c', 'd', 'e']

    cases = {
        '0': ['a'],
        '0,1': ['a', 'b'],
        '1,3': ['b', 'd']
    }

    for constructor, result in cases.items():
        testes_indices = Indices(constructor)
        assert testes_indices.get(items) == result


def test_positive_int():

    with pytest.raises(ValueError):
        positive_int('-5')

    assert positive_int('5') == 5
