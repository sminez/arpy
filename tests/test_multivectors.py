from arpy import MultiVector, Alpha, Pair


m1 = MultiVector('1 2 3')
m2 = MultiVector('1 2')
m3 = MultiVector('1 2 12')


def test_multivector_construction():
    '''Alternate construction methods are equal'''
    assert MultiVector('3 1 2') == m1
    assert MultiVector(['1', '2', '3']) == m1
    assert MultiVector([Alpha(i) for i in '123']) == m1
    assert MultiVector(Alpha(i) for i in '123') == m1
    assert MultiVector(Pair(i) for i in '123') == m1


def test_simplification():
    '''MultiVector auto simplification works'''
    # Like MultiVectors cancel entirely
    assert m1 - m1 == MultiVector()
    # Unmatched terms are unaffected
    assert m1 - m2 == MultiVector('3')
    # Subtraction works along with simplification
    assert m1 - m3 == MultiVector([Alpha('3'), Alpha('-12')])


def test_addition():
    '''Adding multivectors should combine their components _and_ simplify'''
    assert m1 + m1 == MultiVector('1 2 3 1 2 3')
    assert m1 + m2 == MultiVector('1 1 2 2 3')
    assert m1 + m2 + m3 == MultiVector('1 1 1 2 2 2 3 12')
