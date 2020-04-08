from enum import Enum
from functools import lru_cache


class DivisionType(Enum):
    """
    Given that flipping two elements in a computation results in a change of sign,
    we need to define whether division acts from the right (by) or from the left
    (into) as a convention to use for Python's "/" operator
    """

    BY = "by"
    INTO = "into"


class Zet(Enum):
    """
    Zets form a partition of the algebra and each consists of a single `time
    like' element paired with three `space like' elements. While this sounds
    similar to traditional Relativistic four-vectors, Zets do not transform as
    four-vectors; though they do posess a number of very useful properties that
    shed light on the internal structure of AR and the full product.

    Zet partitions are independent of the metric and internal ordering of their
    elements, only depending on the set of basis alphas making up each element.
    """

    B = "B"
    T = "T"
    A = "A"
    E = "E"

    def superscript(self) -> str:
        _superscripts = {"B": "ᴮ", "A": "ᴬ", "T": "ᵀ", "E": "ᴱ"}
        return _superscripts[self.name]

    @classmethod
    @lru_cache(maxsize=16)
    def from_index(self, ix: str) -> "Zet":
        if ix in map(frozenset, "p 23 31 12".split()):
            return Zet.B
        elif ix in map(frozenset, "0 023 031 012".split()):
            return Zet.T
        elif ix in map(frozenset, "123 1 2 3".split()):
            return Zet.A
        elif ix in map(frozenset, "0123 01 02 03".split()):
            return Zet.E

        raise ValueError(f"{ix} is an invalid index")


class Orientation(Enum):
    """
    Each of the cardinal basis dimensions within the algebra that an element
    or value can take as a primary orientation.
    """

    T = "time"
    X = "space-x"
    Y = "space-y"
    Z = "space-z"

    @classmethod
    @lru_cache(maxsize=16)
    def from_index(self, ix: str) -> "Orientation":
        if ix in map(frozenset, "p 0 123 0123".split()):
            return Orientation.T
        elif ix in map(frozenset, "23 023 1 01".split()):
            return Orientation.X
        elif ix in map(frozenset, "31 031 2 02".split()):
            return Orientation.Y
        elif ix in map(frozenset, "12 012 3 03".split()):
            return Orientation.Z

        raise ValueError(f"{ix} is an invalid index")
