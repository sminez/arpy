from copy import copy

from ...config import ARConfig
from ...config import config as cfg
from ...utils.utils import SUB_SCRIPTS


class Alpha:
    """
    An Alpha represents a pure element of the algebra without magnitude.
    It is composed of 0-4 Dimensions with the number of dimensions determining
    its nature: i.e. scalar, vector, bivector, trivector, quadrivector
    """

    def __init__(self, index: str, sign: int = 1, cfg: ARConfig = cfg):
        if sign not in [1, -1]:
            raise ValueError("Invalid α sign: {}".format(sign))

        if index.startswith("-"):
            index = index[1:]
            sign *= -1

        if index not in cfg.allowed + cfg.allowed_groups:
            raise ValueError("Invalid α index: {}".format(index))

        self._index = index
        self._sign = sign
        self.allowed = cfg.allowed
        self.allowed_groups = cfg.allowed_groups

    @property
    def sign(self):
        return self._sign

    def __repr__(self):
        neg = "-" if self._sign == -1 else ""
        try:
            return "{}α{}".format(neg, "".join(SUB_SCRIPTS[i] for i in self._index))
        except KeyError:
            return "{}α{}".format(neg, self._index)

    def __tex__(self):
        neg = "-" if self.sign == -1 else ""
        return neg + "\\alpha_{" + self._index + "}"

    def __eq__(self, other):
        if not isinstance(other, Alpha):
            return False

        return all([(self._index == other._index), (self._sign == other._sign)])

    def __lt__(self, other):
        try:
            allowed = self.allowed + self.allowed_groups
            return allowed.index(self._index) < allowed.index(other._index)
        except ValueError:
            raise TypeError(
                f"Inconsistant config detected:\n{self} -> {self.cfg}\n{other} -> {other.cfg}"
            )

    def __neg__(self):
        neg = copy(self)
        neg._sign *= -1
        return neg

    def __hash__(self):
        return hash((self._index, self._sign))
