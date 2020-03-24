from itertools import groupby

from ..algebra.data_types import Term
from ..utils.utils import Nat, Zet


def replace_grad(mvec):
    """
    Parse out any gradients in a given multivector:

        Grad f = αx[df/dx] + αy[df/dy] + αz[df/dz]

    NOTE: This includes equivalent term patterns generated by
          differentials with respect to other elements of the
          algebra.
    """
    parsed_terms = []

    for xi, terms in groupby(mvec._terms, lambda t: t._components):
        if len(terms) != 3:
            continue


def replace_curl(mvec):
    """
    Parse out any curls in a given multivector:

        Curl F = αx[dFz/dy-dFy/dz] + αy[dFx/dz-dFz/dx] + αz[dFy/dx-dFx/dy]

    NOTE: This includes equivalent term patterns generated by
          differentials with respect to other elements of the
          algebra.
    """
    parsed_terms = []


def replace_div(mvec):
    """
    Parse out any divergences in a given multivector:

        Div F = dFx/dx + dFy/dy + dFz/dz

    NOTE: This includes equivalent term patterns generated by
          differentials with respect to other elements of the
          algebra.
    """
    parsed_terms = []


def replace_3vector_partials(mvec):
    """
    Parse out any instances of differentiating all of the components
    of a 3vector by the same index in a given multivector:
    """
    parsed_terms = []
