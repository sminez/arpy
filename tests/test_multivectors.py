from arpy import MultiVector, Alpha


def test_simplification():
    '''MultiVector auto simplification works'''
    m1 = MultiVector('1 2 3')
    m2 = MultiVector('1 2')
    m3 = MultiVector('1 2 12')

    # Like MultiVectors cancel entirely
    assert m1 - m1 == MultiVector()
    # Unmatched terms are unaffected
    assert m1 - m2 == MultiVector('3')
    # Subtraction works along with simplification
    assert m1 - m3 == MultiVector([Alpha('3'), Alpha('-12')])
