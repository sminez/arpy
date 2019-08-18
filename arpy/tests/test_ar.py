import pytest
from arpy import Alpha, Term, ar, full


def test_creation():
    """ar creates equivalent types to the standard classes"""
    assert ar("a23") == Alpha("23")
    assert ar("p1") == Term("1")


def test_product():
    """ar full product work correctly"""
    assert ar("a1 ^ a2") == Alpha("12")
    assert ar("a2 ^ a1") == Alpha("-12")
    assert ar("a2 ^ -a1") == Alpha("12")
    assert ar("-p1 ^ a2") == Term("-12", "1")
    assert ar("p2 ^ a1") == Term("-12", "2")
    assert ar("a1 ^ p2") == Term("12", "2")
    assert ar("a2 ^ p1") == Term("-12", "1")


def test_division():
    """ar division works correctly"""
    assert ar("a1 / a2") == Alpha("-12")
    assert ar("a2 / a1") == Alpha("12")
    assert ar("p1 / a2") == Term("-12", "1")
    assert ar("p2 / a1") == Term("12", "2")

    with pytest.raises(NotImplementedError):
        ar("a1 / p2")
    with pytest.raises(NotImplementedError):
        ar("a2 / p1")


def test_differential():
    """Operations with differentials work as expected"""
    from arpy import d, A

    assert ar("d A") == d(A)
    assert ar("A d") == d(A, div="by")
    assert ar("A (d A)") == full(A, d(A))
    # TODO: look at pytest's capsys feature
    assert ar("foo ") is None
    foo = Alpha("1")
    assert ar("foo ^ foo") == Alpha("-p")
