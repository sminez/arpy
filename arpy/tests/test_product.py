from arpy import Alpha, ARConfig, Term, config, find_prod, full
from utils import ij_pairs, ijk_triplets, metrics

ap = Alpha("p")
neg_ap = Alpha("-p")


def test_alpha_p():
    """
    Left and right multiplication by αp is idempotent.
    Left and right multiplication by -αp negates.
    """
    for index in config.allowed:
        for metric in metrics:
            new_config = ARConfig(config.allowed, metric, config.division_type)
            test_alpha = Alpha(index)
            assert find_prod(ap, test_alpha, cfg=new_config) == test_alpha
            assert find_prod(test_alpha, ap, cfg=new_config) == test_alpha

            negative_test_alpha = Alpha(index, -1)
            left_neg_ap = find_prod(neg_ap, test_alpha, cfg=new_config)
            right_neg_ap = find_prod(test_alpha, neg_ap, cfg=new_config)
            assert left_neg_ap == negative_test_alpha
            assert right_neg_ap == negative_test_alpha


def test_swapped_indices_ij():
    """
    Swapping two indices negates the product
    """
    for i, j in ij_pairs:
        i, j = Alpha(i), Alpha(j)
        ij, ji = find_prod(i, j), find_prod(j, i)
        assert ij._index == ji._index
        assert ij._sign == ji._sign * -1


def test_swapped_indices_ijk():
    """
    Swapping two indices negates the product
    """
    for i, j, k in ijk_triplets:
        i, j, k = Alpha(i), Alpha(j), Alpha(k)
        ijk = find_prod(find_prod(i, j), k)
        jik = find_prod(find_prod(j, i), k)
        assert ijk._index == jik._index
        assert ijk._sign == jik._sign * -1


def test_squaring():
    """
    The product of an element with itself is +-αp
    """
    for index in config.allowed:
        for metric in metrics:
            new_config = ARConfig(config.allowed, metric, config.division_type)
            pos, neg = Alpha(index), Alpha(index, -1)
            assert find_prod(pos, pos, new_config) in [ap, neg_ap]
            assert find_prod(neg, neg, new_config) in [ap, neg_ap]
            assert find_prod(pos, neg, new_config) in [ap, neg_ap]
            assert find_prod(pos, neg, new_config) in [ap, neg_ap]


def test_term_products():
    """
    The alpha generated by forming the product of two terms is the
    same as the alpha generated by forming the product of the term alphas
    on their own.
    """
    for istr in config.allowed:
        for jstr in config.allowed:
            ai, aj = Alpha(istr), Alpha(jstr)
            ti, tj = Term(istr), Term(jstr)
            assert full(ti, tj).alpha == full(ai, aj)


def test_full_alpha_term():
    """
    Should get the same result for Alpha^Term and Term^Alpha if this is just
    adding the same xi to either Alpha.
    """
    for i in config.allowed:
        for j in config.allowed:
            ai, pi = Alpha(i), Term(i, "test")
            aj, pj = Alpha(j), Term(j, "test")
            assert full(ai, pj) == full(pi, aj)
