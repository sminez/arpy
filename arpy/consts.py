from dataclasses import dataclass
from enum import Enum

from .config import ARConfig


@dataclass
class ZetElements:
    e: str
    x: str
    y: str
    z: str

    @property
    def all(self):
        return [self.e, self.x, self.y, self.z]

    def __iter__(self):
        yield from self.all


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

    @classmethod
    def from_index(self, ix_str: str) -> "Zet":
        ix = frozenset(ix_str)
        if ix in [frozenset(s) for s in "p 23 31 12".split()]:
            return Zet.B
        elif ix in [frozenset(s) for s in "0 023 031 012".split()]:
            return Zet.T
        elif ix in [frozenset(s) for s in "123 1 2 3".split()]:
            return Zet.A
        elif ix in [frozenset(s) for s in "0123 01 02 03".split()]:
            return Zet.E

        raise ValueError(f"{ix} is an invalid index")

    def elements(self, config: ARConfig) -> ZetElements:
        comps = config.zet_comps[self.name]
        return ZetElements(**comps)

    @property
    def superscript(self) -> str:
        # _superscripts = {"B": "ᴮ", "A": "ᴬ", "T": "ᵀ", "E": "ᴱ"}
        _superscripts = {"B": "ᴮ", "A": "", "T": "ᵀ", "E": "ᴱ"}
        return _superscripts[self.name]

    @property
    def time_like(self) -> str:
        if self.name == "B":
            return "p"
        elif self.name == "T":
            return "t"
        elif self.name == "A":
            return "h"
        elif self.name == "E":
            return "q"

        raise ValueError("Invalid zet")


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
    def from_index(self, ix_str: str) -> "Orientation":
        ix = frozenset(ix_str)
        if ix in [frozenset(s) for s in "p 0 123 0123".split()]:
            return Orientation.T
        elif ix in [frozenset(s) for s in "23 023 1 01".split()]:
            return Orientation.X
        elif ix in [frozenset(s) for s in "31 031 2 02".split()]:
            return Orientation.Y
        elif ix in [frozenset(s) for s in "12 012 3 03".split()]:
            return Orientation.Z

        raise ValueError(f"{ix} is an invalid index")
