from arpy import Alpha, find_prod, ALLOWED


def test_alpha_p():
    '''αp.αμ = αμ ∀ μ ∈ allowed'''
    ap = Alpha('p')
    neg_ap = Alpha('p', -1)

    for index in ALLOWED:
        test_alpha = Alpha(index)
        assert find_prod(ap, test_alpha) == test_alpha
        assert find_prod(test_alpha, ap) == test_alpha

        negative_test_alpha = Alpha(index, -1)
        assert find_prod(neg_ap, test_alpha) == negative_test_alpha
        assert find_prod(test_alpha, neg_ap) == negative_test_alpha
