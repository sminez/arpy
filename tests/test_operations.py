import pytest
from utils import imetric
from arpy import Alpha, find_prod, inverse


@pytest.mark.parametrize('index,metric', imetric)
def test_inverse(index, metric):
    '''
    Inverting an α and multiplying by the original should give αp
    '''
    alpha = Alpha(index)
    inverse_alpha = inverse(alpha, metric=metric)
    assert find_prod(alpha, inverse_alpha, metric=metric) == Alpha('p')
