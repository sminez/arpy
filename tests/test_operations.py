import pytest
from utils import metrics, ap, neg_ap
from arpy import Alpha, Pair, MultiVector, find_prod, inverse, dagger, \
    commutator, full, ALLOWED


@pytest.mark.parametrize('metric', metrics)
def test_inverse(metric):
    '''
    Inverting an α and multiplying by the original should give αp
    '''
    for index in ALLOWED:
        alpha = Alpha(index)
        inverse_alpha = inverse(alpha, metric=metric)
        assert find_prod(alpha, inverse_alpha, metric=metric) == Alpha('p')


@pytest.mark.parametrize('metric', metrics)
def test_dagger(metric):
    '''Dagger negates elements that square to -1'''
    alphas = [
        Alpha(a, find_prod(Alpha(a), Alpha(a), metric=metric).sign)
        for a in ALLOWED
    ]
    negated = MultiVector([Pair(a) for a in alphas])
    assert negated == dagger(MultiVector(ALLOWED), metric=metric)


@pytest.mark.parametrize('metric', metrics)
def test_commutator(metric):
    '''Commutator results should always be +-αp'''
    for i in ALLOWED:
        ai = Alpha(i)
        for j in ALLOWED:
            aj = Alpha(j)
            assert commutator(ai, aj, metric=metric) in [ap, neg_ap]
