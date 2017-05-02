from utils import metrics
from arpy import Alpha, Pair, MultiVector, find_prod, inverse, dagger, \
    commutator, project, ALLOWED

ap = Alpha('p')
neg_ap = Alpha('-p')


def test_inverse():
    '''
    Inverting an α and multiplying by the original should give αp
    '''
    for metric in metrics:
        for index in ALLOWED:
            alpha = Alpha(index)
            inverse_alpha = inverse(alpha, metric=metric)
            assert find_prod(alpha, inverse_alpha, metric=metric) == ap


def test_dagger():
    '''Dagger negates elements that square to -1'''
    for metric in metrics:
        alphas = [
            Alpha(a, find_prod(Alpha(a), Alpha(a), metric=metric).sign)
            for a in ALLOWED
        ]
        negated = MultiVector([Pair(a) for a in alphas])
        assert negated == dagger(MultiVector(ALLOWED), metric=metric)


def test_commutator():
    '''Commutator results should always be +-αp'''
    for metric in metrics:
        for i in ALLOWED:
            ai = Alpha(i)
            for j in ALLOWED:
                aj = Alpha(j)
                assert commutator(ai, aj, metric=metric) in [ap, neg_ap]


def test_project_alpha():
    '''
    Projecting an Alpha should return the value if and only if the
    number of indices on the Alpha is the grade we are projecting
    '''
    # NOTE:: (deggen) The project function is in arpy/algebra/operations.py
    pass


def test_project_pair():
    '''
    Projecting a Pair should return the value if and only if the
    number of indices on the Pair's Alpha is the grade we are projecting
    '''
    pass


def test_project_multivector():
    '''
    Projecting a MultiVector should return a new MultiVector that contains
    only pairs that are of the correct grade
    '''
    # NOTE:: (Deggen) The `new MultiVector` part is important! ;)
    #        Have a look at test_new_component_partial in the
    #        test_differentials file.
    pass
