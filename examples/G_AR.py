"""
Looking at G_AR - AR without sign.

Formally, this is the abelian group of the Power set on n elements using set
symmetric difference as the group operator. This gives a surprising amount of
the structure within AR and provides an easier introduction to understanding
how the full product works.

NOTE:
  This script expects to be run in a terminal that understands ANSI color
  escape sequences in order to be able to generate its output.
"""
from functools import reduce
from itertools import combinations
from typing import Set


class Element:
    def __init__(self, raw: str) -> "Element":
        self.raw = raw
        self.zet = self.__zet()
        self.orientation = self.__orientation()

    def zet_repr(self):
        return f"{self.zet}{self.orientation}"

    def color_repr(self, zet_only=False, orientation_only=False):
        base = {"B": 24, "T": 91, "A": 160, "E": 28}
        offset = {"t": 0, "x": 6, "y": 12, "z": 18}

        if zet_only:
            s, i = self.orientation, "BTAE".index(self.zet) + 1
        elif orientation_only:
            s, i = self.zet, "txyz".index(self.orientation) + 1
        else:
            s, i = self.zet_repr(), base[self.zet] + offset[self.orientation]

        return f"\x1b[48;5;{str(i).rjust(3, '0')}m{s} "

    def __zet(self):
        if self.raw in "p 23 13 12".split():
            return "B"
        elif self.raw in "0 023 013 012".split():
            return "T"
        elif self.raw in "123 1 2 3".split():
            return "A"
        elif self.raw in "0123 01 02 03".split():
            return "E"

    def __orientation(self):
        if self.raw in "p 0 123 0123".split():
            return "t"
        elif self.raw in "23 023 1 01".split():
            return "x"
        elif self.raw in "13 013 2 02".split():
            return "y"
        elif self.raw in "12 012 3 03".split():
            return "z"


def compose(a: Set, b: Set) -> str:
    s = a.symmetric_difference(b)
    product = "".join(sorted(s))
    return product if product else "p"


def flatten(col):
    if isinstance(col, list):
        return reduce(lambda a, b: a + b, map(flatten, col))
    else:
        return [col]


class G_AR:
    def __init__(self, generated=False):
        self.order = 4
        self.indices = [str(n) for n in range(self.order)]
        if generated:
            self.elements = flatten(
                [[set(c) for c in combinations(self.indices, k)] for k in range(self.order + 1)]
            )
        else:
            self.elements = [set()] + [
                set(s) for s in "0 1 2 3 01 02 03 23 31 12 023 031 012 123 0123".split()
            ]

        self.cayley = []

        for a in self.elements:
            row = []
            for b in self.elements:
                product = compose(a, b)
                product = "p" if product == "" else product
                row.append(Element(product))
            self.cayley.append(row)

    def cayley(self):
        for row in self.cayley:
            print(" ".join(e.raw.ljust(self.order) for e in row))

    def zet_cayley(self, zet_only=False, orientation_only=False):
        for row in self.cayley:
            s = "".join(e.color_repr(zet_only, orientation_only) for e in row)
            print(f"{s}\x1b[0m")
        print()


if __name__ == "__main__":
    g = G_AR()
    g.zet_cayley()
    g.zet_cayley(zet_only=True)
    g.zet_cayley(orientation_only=True)
