import pytest

from .. import ARConfig, config

original_allowed = set(config.allowed)
new_allowed = [
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
    "10",
    "20",
    "30",
]


def test_set_reset():
    """Resetting correctly resets the config"""
    cfg = ARConfig(metric=config.metric, allowed=config.allowed, div=config.division_type)
    cfg.allowed = new_allowed
    assert set(cfg.allowed) == set(new_allowed)
    cfg.reset()
    assert set(cfg.allowed) == set(original_allowed)


def test_invalid_allowed_length():
    """Changing allowed gets rejected if it is the wrong length"""
    cfg = ARConfig(metric=config.metric, allowed=config.allowed, div=config.division_type)

    with pytest.raises(ValueError):
        # valid indices but too short
        cfg.allowed = new_allowed[:-1]

    with pytest.raises(ValueError):
        # valid indices but too long
        cfg.allowed = new_allowed + ["01"]


def test_invalid_allowed_indices():
    """Changing allowed gets rejected if indices aren't p0123"""
    cfg = ARConfig(metric=config.metric, allowed=config.allowed, div=config.division_type)

    with pytest.raises(ValueError):
        # invalid indices
        cfg.allowed = ["not", "valid", "at", "all!"]

    with pytest.raises(ValueError):
        # invalid indices, correct length
        cfg.allowed = [i + "foo" for i in new_allowed]


def test_invalid_metric_length_string():
    """Changing the metric gets rejected if it is too short"""
    cfg = ARConfig(metric=config.metric, allowed=config.allowed, div=config.division_type)

    with pytest.raises(ValueError):
        cfg.metric = "+--"  # short

    with pytest.raises(ValueError):
        cfg.metric = "+---+"  # long


def test_invalid_metric_invalid_string():
    """Changing the metric gets rejected if it is too short"""
    cfg = ARConfig(metric=config.metric, allowed=config.allowed, div=config.division_type)

    with pytest.raises(ValueError):
        cfg.metric = "foo"  # short

    with pytest.raises(ValueError):
        cfg.metric = "foop"  # correct len

    with pytest.raises(ValueError):
        cfg.metric = "foobar"  # long


def test_invalid_metric_length_tuple():
    """Changing the metric gets rejected if it is too short"""
    cfg = ARConfig(metric=config.metric, allowed=config.allowed, div=config.division_type)

    with pytest.raises(ValueError):
        cfg.metric = (1, -1, -1)  # short

    with pytest.raises(ValueError):
        cfg.metric = (1, 1, 1, -1, -1)  # long


def test_invalid_metric_invalid_tuple():
    """Changing the metric gets rejected if it is too short"""
    cfg = ARConfig(metric=config.metric, allowed=config.allowed, div=config.division_type)

    with pytest.raises(ValueError):
        cfg.metric = ("foo", "bar")  # short

    with pytest.raises(ValueError):
        cfg.metric = ("not", "valid", "because", "reasons")  # correct len

    with pytest.raises(ValueError):
        cfg.metric = ("foo", "bar", "baz", "spam", "grok")  # long
