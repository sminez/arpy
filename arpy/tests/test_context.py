import pytest
from arpy import ARContext


oi_allowed = [
    'p', '23', '31', '12', '0', '023', '031', '012',
    '123', '1', '2', '3', '0123', '01', '02', '03'
]
io_allowed = [
    'p', '23', '31', '12', '0', '023', '031', '012',
    '123', '1', '2', '3', '0123', '10', '20', '30'
]

ctx1 = ARContext(oi_allowed, '+---', 'into')
ctx2 = ARContext(io_allowed, '+---', 'into')


def test_builtins():
    '''The builtins should update with the context they are created in'''
    with ctx1 as ar:
        E1 = ar('E')
        assert E1 == ar('{ 01 02 03 }')

    with ctx2 as ar:
        E2 = ar('E')
        assert E2 == ar('{ 10 20 30 }')

    assert E1 != E2


def test_allowed():
    '''
    The threading of the config through to operations should cause them to
    fail when they encounter a value from another config.
    '''
    with ctx1 as ar:
        A1 = ar('a01')
        A2 = ar('a02')

    with pytest.raises(ValueError):
        with ctx2 as ar:
            ar('A1 ^ A2')
