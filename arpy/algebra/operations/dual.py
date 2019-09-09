from ..data_types import MultiVector
from .full import full


def dual(mvec):
    """
    The dual of a MultiVector is defined as 'M_bar = -a0123 ^ M'
    """
    if not isinstance(mvec, MultiVector):
        raise ValueError("Can only compute the dual of a MultiVector")

    # Ensure that we take q = a0123 as the correct Alpha to be in alignment
    # with the MultiVector that was passed
    q = mvec.cfg.q._terms[0].alpha
    return full(-q, mvec)


def MM_bar(mvec, no_cancel=False):
    """
    The product M ^ M_bar
    """
    result = full(mvec, dual(mvec))

    if not no_cancel:
        result.cancel_terms()

    return result
