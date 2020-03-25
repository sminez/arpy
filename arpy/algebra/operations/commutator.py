from ...config import config as cfg
from ...utils.concepts.dispatch import dispatch_on
from ..data_types import Alpha
from .full import full, inverse


@dispatch_on((0, 1))
def commutator(a, b, cfg=cfg):
    """
    Computes the group commutator [a, b] = (a . b . a^-1 . b^-1) for Alphas.
    """
    raise NotImplementedError


@commutator.add((Alpha, Alpha))
def _group_commutator(a, b, cfg=cfg):
    product = full(a, b, cfg)
    product = full(product, inverse(a, cfg=cfg), cfg)
    product = full(product, inverse(b, cfg=cfg), cfg)
    return product


# @commutator.add((Term, Term))
# def _ring_commutator(a, b):
#     '''Computes the ring commutator [a, b] = ab - ba for Pairs'''
#     return full(a, b) - full(b, a)
