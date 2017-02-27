import pytest
from utils import ap, neg_ap, imetric, ijmetric, ijkmetric
from arpy import Alpha, find_prod, ALLOWED


@pytest.mark.parametrize('index,metric', imetric)
def test_alpha_p(index, metric):
    '''
    Left and right multiplication by αp is idempotent.
    Left and right multiplication by -αp negates.
    '''
    test_alpha = Alpha(index)
    assert find_prod(ap, test_alpha, metric=metric) == test_alpha
    assert find_prod(test_alpha, ap, metric=metric) == test_alpha

    negative_test_alpha = Alpha(index, -1)
    assert find_prod(neg_ap, test_alpha, metric=metric) == negative_test_alpha
    assert find_prod(test_alpha, neg_ap, metric=metric) == negative_test_alpha


@pytest.mark.parametrize('i,j,metric', ijmetric)
def test_swapped_indices_ij(i, j, metric):
    '''
    Swapping two indices negates the product
    '''
    i, j = Alpha(i), Alpha(j)
    ij, ji = find_prod(i, j, metric=metric), find_prod(j, i, metric=metric)
    assert ij.index == ji.index
    assert ij.sign == ji.sign * -1


@pytest.mark.parametrize('i,j,k,metric', ijkmetric)
def test_swapped_indices_ijk(i, j, k, metric):
    '''
    Swapping two indices negates the product
    '''
    i, j, k = Alpha(i), Alpha(j), Alpha(k)
    ijk = find_prod(find_prod(i, j, metric=metric), k, metric=metric)
    jik = find_prod(find_prod(j, i, metric=metric), k, metric=metric)
    assert ijk.index == jik.index
    assert ijk.sign == jik.sign * -1


@pytest.mark.parametrize('index,metric', imetric)
def test_squaring(index, metric):
    '''
    The product of an element with itself is +-αp
    '''
    pos, neg = Alpha(index), Alpha(index, -1)

    assert find_prod(pos, pos, metric) in [ap, neg_ap]
    assert find_prod(neg, neg, metric) in [ap, neg_ap]
    assert find_prod(pos, neg, metric) in [ap, neg_ap]
    assert find_prod(pos, neg, metric) in [ap, neg_ap]
