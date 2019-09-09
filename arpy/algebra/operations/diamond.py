from ..data_types import MultiVector
from .project import project


def diamond(mvec):
    """
    The diamond conjugate of a MultiVector is defined as
        'M_diamond = 2<M>0 - M

    It negates everything with a 'direction' (e.g. not ap)
    """
    if not isinstance(mvec, MultiVector):
        raise ValueError("Can only compute the diamond of a MultiVector")

    scalar = project(mvec, 0) * 2
    return scalar - mvec
