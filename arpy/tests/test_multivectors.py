import pytest
from arpy import MultiVector, Alpha, Pair


m1 = MultiVector("1 2 3")
m2 = MultiVector("1 2")
m3 = MultiVector("1 2 12")


def test_multivector_construction():
    """Alternate construction methods are equal"""
    assert MultiVector("3 1 2") == m1
    assert MultiVector("-1 -2 -3") == -m1
    assert MultiVector(["1", "2", "3"]) == m1
    assert MultiVector(("1", "2", "3")) == m1
    assert MultiVector({"1", "2", "3"}) == m1
    assert MultiVector([Alpha(i) for i in "123"]) == m1
    assert MultiVector(Alpha(i) for i in "123") == m1
    assert MultiVector(Pair(i) for i in "123") == m1
    assert MultiVector(m1) == m1

    with pytest.raises(ValueError):
        MultiVector([1, 2, 3, 4])
        MultiVector("foo")


def test_simplification():
    """MultiVector auto simplification works"""
    # Like MultiVectors cancel entirely
    res1 = m1 - m1
    res1.cancel_terms()
    assert res1 == MultiVector()
    # Unmatched terms are unaffected
    res2 = m1 - m2
    res2.cancel_terms()
    assert res2 == MultiVector("3")
    # Subtraction works along with simplification
    res3 = m1 - m3
    res3.cancel_terms()
    assert res3 == MultiVector([Alpha("3"), Alpha("-12")])


def test_addition():
    """Adding multivectors should combine their components _and_ simplify"""
    assert m1 + m1 == MultiVector("1 2 3 1 2 3")
    assert m1 + m2 == MultiVector("1 1 2 2 3")
    assert m1 + m2 + m3 == MultiVector("1 1 1 2 2 2 3 12")


def test_contains():
    """We can test for membership in a MultiVector"""
    assert Alpha("2") in m1
    assert Pair("2") in m1
    assert Alpha("12") not in m1
    assert Pair("12") not in m1
    assert "Maxwell" not in m1


def test_iteration_pairs():
    """We can iterate over a MultiVector and get pairs"""
    for p in m3:
        assert isinstance(p, Pair)


def test_get_item():
    """Dict syntax should return a list of pairs"""
    assert m1["1"] == [Pair("1")]
    assert m1["2"] == [Pair("2")]
    assert m1["3"] == [Pair("3")]


def test_del_item():
    """We can delete Xis on a given index"""
    m = MultiVector("1 2 3 1 2 3")
    del m["1"]
    assert m == MultiVector("2 2 3 3")
    del m["2"]
    assert m == MultiVector("3 3")
    del m["3"]
    assert m == MultiVector()


def test_copy():
    """Copy returns a new object"""
    new_mvec = m1.copy()
    assert new_mvec is not m1


def test_relabel():
    """Relabelling gives a new MultiVector with correct components"""
    relabelled = m1.relabel("1", "foo")
    assert relabelled is not m1
    assert relabelled == MultiVector([Pair("1", "foo"), "2", "3"])
    relabelled_again = relabelled.relabel_many([("2", "bar"), ("3", "baz")])
    assert relabelled_again is not relabelled
    assert relabelled_again is not m1
    assert relabelled_again == MultiVector([Pair("1", "foo"), Pair("2", "bar"), Pair("3", "baz")])


def test_len():
    """The length of a Multivector is the number of Xi components"""
    assert len(m1) == 3
    assert len(m2) == 2
    assert len(m3) == 3
    large = MultiVector("1 1 2 3 12 12 31 p p p")
    assert len(large) == 10
