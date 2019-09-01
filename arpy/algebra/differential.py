"""
A selection of different implementations for symbolically computing
the 4-vector 4-differential Dμ and other Differential operators.

In Cartesian coordinates Dμ is:
    Dμ = ∂μ/αμ = ∂ / (αμ∂xμ) = α0∂0 - αi∂i = α0∂0 - ∇i

All other differential operators follow the same restrictions of
Absolute Relativity and should only operate on MultiVectors.

NOTE:: Specific operators (such as Dmu) are defined in the __init__ file.
"""
from copy import deepcopy

from ..utils.utils import SUB_SCRIPTS
from .config import config as cfg
from .data_types import Alpha, MultiVector
from .operations import div_by, div_into, full, inverse


class AR_differential:
    """Differential operator: can be used inside of ar()"""

    def __init__(self, wrt, cfg=cfg):
        if isinstance(wrt, MultiVector):
            self.wrt = [term._alpha for term in wrt]
        else:
            if isinstance(wrt, str):
                wrt = wrt.split()

            if not isinstance(wrt, list):
                raise ValueError(
                    "Differential operators must be initialised with either"
                    " a MultiVector, list or string of alpha indices"
                )

            # Conversion to Alpha catches invalid indices
            self.wrt = [Alpha(comp, cfg=cfg) for comp in wrt]

        self.cfg = cfg

        alphas = ", ".join([str(a) for a in self.wrt])
        self.__doc__ = "Differnetiate with respect to: {}".format(alphas)

    def __call__(self, mvec, cfg=None, div=None):
        """
        Compute the result of Differentiating a each component of a MultiVector
        with respect to a given list of unit elements under the algebra.
        """
        comps = []
        if cfg is None:
            cfg = self.cfg

        for term in mvec:
            for element in self.wrt:
                result = term_partial(term, element, cfg, div)
                comps.append(result)

        return MultiVector(comps, cfg=cfg)

    def __repr__(self):
        elements = [
            "{}∂{}".format(str(inverse(a, cfg=self.cfg)), "".join(SUB_SCRIPTS[i] for i in a._index))
            for a in self.wrt
        ]
        return "{ " + " ".join(elements) + " }"

    def __tex__(self):
        elements = []
        for a in self.wrt:
            inv = inverse(a, cfg=self.cfg)
            sign = "" if inv._sign == 1 else "-"

            elements.append("%s\\alpha_{%s}\\partial_{%s}" % (sign, a._index, a._index))

        return r"\{ " + " ".join(elements) + r" \}"


def _div(alpha, wrt, cfg, div=None):
    """Divide an alpha component based on the set division type"""
    div = div if div else cfg.division_type
    if div == "by":
        return div_by(alpha, wrt, cfg)
    elif div == "into":
        return div_into(wrt, alpha, cfg)
    else:
        raise ValueError("Invalid division specification: {}".format(cfg.division_type))


def term_partial(term, wrt, cfg, div):
    """
    Symbolically differentiate a term by storing the partials and
    converting the alpha value using the correct division type.
    """
    new_term = deepcopy(term)
    new_term.alpha = _div(new_term.alpha, wrt, cfg, div)

    if len(new_term._components) == 1:
        new_term._components[0].partials = [wrt] + new_term._components[0]._partials
    else:
        new_term.component_partials = [wrt] + new_term._component_partials

    return new_term


@full.add((AR_differential, MultiVector))
def _full_differential_mvec(diff, mvec, cfg=cfg):
    res = diff(mvec, cfg=cfg)
    return res


@full.add((MultiVector, AR_differential))
def _full_mvec_differential_mvec(mvec, diff, cfg=cfg):
    res = diff(mvec, cfg=cfg, div="by")
    return res
