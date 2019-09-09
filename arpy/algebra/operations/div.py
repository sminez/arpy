from ...utils.concepts.dispatch import dispatch_on
from ..config import config as cfg
from ..data_types import Alpha, Term
from .full import find_prod, inverse


@dispatch_on((0, 1))
def div_by(a, b, cfg=cfg):
    """Divide one element by another"""
    raise NotImplementedError


@div_by.add((Alpha, Alpha))
def _div_by_alpha_alpha(a, b, cfg=cfg):
    return find_prod(a, inverse(b, cfg), cfg)


@div_by.add((Term, Alpha))
def _div_by_term_alpha(a, b, cfg=cfg):
    alpha = find_prod(a.alpha, inverse(b, cfg), cfg)
    return Term(alpha, a._components, cfg=cfg)


@dispatch_on((0, 1))
def div_into(a, b, cfg=cfg):
    """Divide one element into another"""
    raise NotImplementedError


@div_into.add((Alpha, Alpha))
def _div_into_alpha_alpha(a, b, cfg=cfg):
    return find_prod(inverse(a, cfg), b, cfg)


@div_into.add((Alpha, Term))
def _div_into_alpha_term(a, b, cfg=cfg):
    alpha = find_prod(inverse(a, cfg), b.alpha, cfg)
    return Term(alpha, b._components, cfg=cfg)
