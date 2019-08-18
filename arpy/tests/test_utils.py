from arpy import reorder_allowed

allowed = [
    "p",
    "23",
    "31",
    "12",
    "0",
    "023",
    "031",
    "012",
    "123",
    "1",
    "2",
    "3",
    "0123",
    "01",
    "02",
    "03",
]


def test_reorder_allowed_idempotent():
    """Reordering to the same order is idempotent"""
    new = reorder_allowed(allowed, order="pBtThAqE")
    assert new == allowed


def test_reorder_allowed_same_elements():
    """Reodering should just be a permutation of the current allowed"""
    new = reorder_allowed(allowed, order="pthqEATB")
    assert set(new) == set(allowed)
