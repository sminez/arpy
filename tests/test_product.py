from itertools import combinations
from arpy import Alpha, find_prod, ALLOWED


ap = Alpha('p')
neg_ap = Alpha('p', -1)
txyz = ['0', '1', '2', '3']


def test_alpha_p():
    '''
    Left and right multiplication by αp is idempotent.
    Left and right multiplication by -αp negates.
    '''
    for index in ALLOWED:
        test_alpha = Alpha(index)
        assert find_prod(ap, test_alpha) == test_alpha
        assert find_prod(test_alpha, ap) == test_alpha

        negative_test_alpha = Alpha(index, -1)
        assert find_prod(neg_ap, test_alpha) == negative_test_alpha
        assert find_prod(test_alpha, neg_ap) == negative_test_alpha


def test_swapped_indices():
    '''
    Swapping two indices negates the product
    '''
    two_cases = [c for c in combinations(txyz, 2)]
    three_cases = [c for c in combinations(txyz, 3)]

    for case in two_cases:
        i, j = Alpha(case[0]), Alpha(case[1])
        ij, ji = find_prod(i, j), find_prod(j, i)
        assert ij.index == ji.index
        assert ij.sign == ji.sign * -1

    for case in three_cases:
        i, j, k = Alpha(case[0]), Alpha(case[1]), Alpha(case[2])
        ijk = find_prod(find_prod(i, j), k)
        jik = find_prod(find_prod(j, i), k)
        assert ijk.index == jik.index
        assert ijk.sign == jik.sign * -1


def test_squaring():
    '''
    The product of an element with itself is +-αp
    '''
    for index in ALLOWED:
        pos = Alpha(index)
        neg = Alpha(index, -1)

        assert find_prod(pos, pos) in [ap, neg_ap]
        assert find_prod(neg, neg) in [ap, neg_ap]
        assert find_prod(pos, neg) in [ap, neg_ap]
        assert find_prod(pos, neg) in [ap, neg_ap]
