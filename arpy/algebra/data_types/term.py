import re
from copy import copy
from typing import List, Union

from ...utils.utils import SUB_SCRIPTS, power_notation
from ..config import ARConfig
from ..config import config as cfg
from .alpha import Alpha
from .xi import Xi

AlphaOrString = Union[Alpha, str]

# Regex for parsing mathematical Xi expressions such as the following:
#   a23[sin kt]
explicit_xi = r"[p0123]*\[.*\]"


class Term:
    """
    A Term represents an AR compliant value: i.e. some value bound to its
    Space-time Nature (Alpha). When constructing a Term from existing Alpha
    and Xi values, the resultant sign of the Term is extracted and the signs
    of all constituant elements are set to 1.
    """

    def __init__(
        self, alpha: AlphaOrString, components: List[Xi] = None, sign: int = 1, cfg: ARConfig = cfg
    ):
        if isinstance(alpha, Alpha):
            if alpha._sign == -1:
                alpha._sign = 1
                sign *= -1
        elif isinstance(alpha, str):
            if len(alpha) > 0 and alpha[0] == "-":
                alpha = alpha[1:]
                sign *= -1

            if components is None and re.match(explicit_xi, alpha):
                # This is something of the form 012[Sin(kx-ωt)]
                # -> split into `012` and `Sin(kx-ωt)`
                alpha, _comp = alpha.split("[")
                components = [Xi(_comp[:-1], cfg=cfg)]  # Remove the trailing ']'

            alpha = Alpha(alpha, cfg=cfg)
        else:
            raise ValueError("Terms are constructed using an Alpha or string")

        if components is None:
            components = [Xi(alpha, cfg=cfg)]

        if isinstance(components, str):
            components = [Xi(components, cfg=cfg)]

        if not all(isinstance(c, Xi) for c in components):
            raise ValueError(f"Term components must be instances of Xi: {components}")

        for comp in components:
            if comp._sign == -1:
                comp._sign = 1
                sign *= -1

        self._components = components
        self._component_partials = []
        self._alpha = alpha
        self._sign = sign
        self.cfg = cfg

    @property
    def alpha(self):
        """The Alpha value of this term corrected for sign"""
        alpha = copy(self._alpha)
        alpha._sign = self._sign
        return alpha

    @alpha.setter
    def alpha(self, a):
        self._sign = a._sign
        a._sign = 1
        self._alpha = a

    @property
    def sign(self):
        return self._sign

    @property
    def index(self):
        """The Alpha index for this Term. i.e. Alpha without sign"""
        return self._alpha._index

    @property
    def component_partials(self):
        return self._component_partials

    @component_partials.setter
    def component_partials(self, val):
        if not isinstance(val, list) and all(isinstance(p, Alpha) for p in val):
            raise ValueError(f"Invalid partials for Term: {val}")

        self._component_partials = sorted(val)

    def __eq__(self, other):
        if not isinstance(other, Term):
            return False

        return all(
            [
                self.cfg == other.cfg,
                self._alpha == other._alpha,
                self._sign == other._sign,
                sorted(self._components) == sorted(other._components),
                sorted(self._component_partials) == sorted(other._component_partials),
            ]
        )

    def __repr__(self):
        sgn = "-" if self._sign == -1 else ""
        comps = ".".join(power_notation(sorted(self._components)))
        partial_strs = [
            "∂{}".format("".join(SUB_SCRIPTS[i] for i in p._index))
            for p in sorted(self._component_partials)
        ]
        partials = "".join(power_notation(partial_strs))

        return f"({sgn}{self._alpha}, {partials}{comps})"

    def _repr_no_alpha(self, ix=0, count=1):
        """
        Provides a string representation of the Xi part of a Term, optionally
        interpolating a count in the case of repeated terms in an iterable.

        The 'ix' parameter allows for ignoring leading terms in Xi products if
        they have been factored out elsewhere within a given string representation.
        """
        sgn = "-" if self._sign == -1 else "+"
        comps = ".".join(power_notation(sorted(self._components[ix:])))
        count = "" if count == 1 else count
        partial_strs = [
            "∂{}".format("".join(SUB_SCRIPTS[i] for i in p._index))
            for p in sorted(self._component_partials)
        ]
        partials = "".join(power_notation(partial_strs))

        return f"{sgn} {count}{partials}{comps}"

    def __hash__(self):
        return hash((self._sign, self._alpha, tuple(sorted(self._components))))

    def __neg__(self):
        neg = copy(self)
        neg._sign *= -1
        return neg

    def __lt__(self, other):
        if not isinstance(other, Term):
            raise TypeError()

        if self._alpha != other._alpha:
            return self._alpha < other._alpha

        return sorted(self._components) < sorted(other._components)
