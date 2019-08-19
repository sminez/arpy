from copy import deepcopy

from ...utils.utils import SUB_SCRIPTS
from ..config import config as cfg
from .alpha import Alpha


class Xi:
    """
    A Xi represents a magnitude within the algebra. The principle of Absolute Relativity
    states that no quantity may appear without its proper Space-time Nature (Alpha).
    As such, raw Xi values have no meaning within AR. That said, we need a method of
    constructing them
    """

    def __init__(self, val, partials=None, sign=1, tex=None, cfg=cfg):
        if isinstance(val, Alpha):
            val = val._index

        if isinstance(val, str) and val.startswith("-"):
            sign = sign * -1
            val = val[1:]

        self._val = val
        self._sign = sign
        self._partials = partials if partials else []
        self._tex_val = tex
        self.cfg = cfg

    @property
    def val(self):
        return self._val

    @property
    def sign(self):
        return self._sign

    def __hash__(self):
        return hash((self._val, self._sign, tuple(self._partials)))

    def __eq__(self, other):
        return all(
            [
                isinstance(other, Xi),
                (self._val == other._val),
                (self._partials == other._partials),
                (self._sign == other._sign),
            ]
        )

    def __lt__(self, other):
        if not isinstance(other, Xi):
            raise TypeError()

        if not all(v in self.cfg.allowed for v in [self._val, other._val]):
            return self._val < other._val

        if self._val != other._val:
            allowed = self.cfg.allowed
            return allowed.index(self._val) < allowed.index(other._val)

        # Comparison for lists is short circuiting & element-wise
        if self._partials != other._partials:
            return self._partials < other._partials

        return self._sign < other._sign

    def __neg__(self):
        neg = deepcopy(self)
        neg._sign *= -1
        return neg

    def __repr__(self):
        sign = "" if self._sign == 1 else "-"
        partials = "".join(
            "∂{}".format("".join(SUB_SCRIPTS[i] for i in p._index))
            for p in reversed(self._partials)
        )
        try:
            display_val = "".join(SUB_SCRIPTS[i] for i in self._val)
            return "{}{}ξ{}".format(sign, partials, display_val)
        except KeyError:
            return "{}{}{}".format(sign, partials, self._val)

    def __tex__(self):
        sign = "+" if self._sign == 1 else "-"
        partials = "".join("\\partial_{" + p._index + "}" for p in reversed(self._partials))
        if self._tex_val is not None:
            return sign + partials + self._tex_val
        elif self._val in self.cfg.allowed:
            return sign + partials + "\\xi_{" + self._val + "}"
        else:
            return sign + partials + self._val
