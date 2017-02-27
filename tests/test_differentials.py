# import pytest
# from utils import ap, neg_ap, imetric, txyz_ijmetric, txyz_ijkmetric
from arpy import Alpha, Pair
from arpy.algebra.differential import component_partial


def test_new_component_partial():
    '''
    component_partial returns a new object rather than modifying
    the original and passing it back
    '''
    original = Pair('012')
    wrt = Alpha('2')
    differentiated = component_partial(original, wrt, 'by', (-1, 1, 1, 1))
    assert differentiated is not original
