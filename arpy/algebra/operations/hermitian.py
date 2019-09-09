from copy import copy

from ...utils.concepts.dispatch import dispatch_on
from ..config import config as cfg
from ..data_types import Alpha, MultiVector, Term
from .full import full
from .rev import rev


@dispatch_on(index=0)
def hermitian(obj, cfg=cfg):
    """
    Compute the Hermitian conjugate of the argument.

    The Hermitian Conjugate is defined to by 'a0 ^ rev(M) ^ a0' for
    MultiVectors with the notation implying that the product is formed
    individually for each term within the MultiVector.
    """
    raise NotImplementedError


# Alias given that the Hermitian Conjugate is more ofter referred to as
# "Mvec dagger" when discussing calculations as a short hand.
dagger = hermitian


def _a0Ma0(elem, cfg):
    a0 = Alpha("0", cfg=cfg)
    Ma0 = full(rev(elem), a0, cfg=cfg)
    return full(a0, Ma0, cfg=cfg)


@hermitian.add(Alpha)
def _hermitian_alpha(alpha, cfg=cfg):
    # return _a0Ma0(alpha, cfg)
    res = copy(alpha)
    if full(alpha, alpha, cfg)._sign == -1:
        res._sign = -1

    return res


@hermitian.add(Term)
def _hermitian_term(term, cfg=cfg):
    # return _a0Ma0(term, cfg)
    res = copy(term)
    if full(term._alpha, term._alpha, cfg)._sign == -1:
        res.sign = -1

    return res


@hermitian.add(MultiVector)
def _hermitian_mvec(mvec, cfg=cfg):
    # return MultiVector([_a0Ma0(t, mvec.cfg) for t in mvec._terms], cfg=mvec.cfg)
    _neg = [
        Alpha(a, cfg=cfg)
        for a in cfg.allowed
        if full(Alpha(a, cfg=cfg), Alpha(a, cfg=cfg), cfg)._sign == -1
    ]
    mvec = copy(mvec)
    new_vec = []
    for term in mvec:
        if term._alpha in _neg:
            term._sign *= -1
        new_vec.append(term)
    res = MultiVector(new_vec, cfg=cfg)

    return res
