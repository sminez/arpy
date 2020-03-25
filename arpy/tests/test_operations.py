from .. import (
    Alpha,
    ARConfig,
    MultiVector,
    Term,
    commutator,
    config,
    dagger,
    find_prod,
    inverse,
    project,
)
from .utils import metrics

ap = Alpha("p")
neg_ap = Alpha("-p")


def test_inverse():
    """
    Inverting an α and multiplying by the original should give αp
    """
    for metric in metrics:
        for index in config.allowed:
            new_config = ARConfig(config.allowed, metric, config.division_type)
            alpha = Alpha(index)
            inverse_alpha = inverse(alpha, cfg=new_config)
            assert find_prod(alpha, inverse_alpha, cfg=new_config) == ap


def test_dagger_alpha():
    """dagger works correctly for alphas"""
    for metric in metrics:
        new_config = ARConfig(config.allowed, metric, config.division_type)

        for ix in new_config.allowed:
            sign = find_prod(Alpha(ix), Alpha(ix), cfg=new_config)._sign
            a1 = Alpha(index=ix, sign=sign, cfg=new_config)
            a2 = dagger(Alpha(ix, cfg=new_config), cfg=new_config)
            assert a1 == a2


def test_dagger():
    """Dagger negates elements that square to -1"""
    for metric in metrics:
        new_config = ARConfig(config.allowed, metric, config.division_type)
        alphas = [
            Alpha(index=a, sign=find_prod(Alpha(a), Alpha(a), cfg=new_config)._sign, cfg=new_config)
            for a in config.allowed
        ]
        negated = MultiVector([Term(a) for a in alphas], cfg=new_config)
        assert negated == dagger(MultiVector(config.allowed), cfg=new_config)


def test_commutator():
    """Commutator results should always be +-αp"""
    for metric in metrics:
        new_config = ARConfig(config.allowed, metric, config.division_type)
        for i in config.allowed:
            ai = Alpha(i)
            for j in config.allowed:
                aj = Alpha(j)
                assert commutator(ai, aj, cfg=new_config) in [ap, neg_ap]


def test_project_alpha():
    """
    Projecting an Alpha should return the value if and only if the
    number of indices on the Alpha is the grade we are projecting
    """
    a = Alpha("01")
    assert project(a, 2) == a
    assert project(a, 0) is None
    assert project(ap, 0) == ap
    assert project(ap, 1) is None


def test_project_term():
    """
    Projecting a Term should return the value if and only if the
    number of indices on the Term's Alpha is the grade we are projecting
    """
    t1 = Term("01")
    t2 = Term("p")
    assert project(t1, 2) == t1
    assert project(t1, 0) is None
    assert project(t2, 0) == t2
    assert project(t2, 1) is None


def test_project_multivector():
    """
    Projecting a MultiVector should return a new MultiVector that contains
    only terms that are of the correct grade
    """
    m1 = MultiVector([Alpha("01")])
    m2 = MultiVector([Alpha("p")])
    assert project(m1, 2) == m1
    assert project(m1, 0) == MultiVector()
    assert project(m2, 0) == m2
    assert project(m2, 2) == MultiVector()
