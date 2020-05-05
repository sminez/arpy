"""
By inspection we can show that for alphas with a single index,
or the Quadrivector, the net sign following pops is unchanged.
For Bivectors and Trivectors the sign is reversed:

a -> a       (no pops: no sign change)
ab -> ba     (1 pop:   sign change)
abc -> cba   (3 pops:  sign change)
  -acb
   cab
  -cba
abcd -> dcba (6 pops: no sign change)
  -abdc
   adbc
  -dabc
   dacb
  -dcab
   dcba

Even though we only need to carry this operation out for objects of
grade 0 -> 4, we can show that the number of pops required for reversing
an Alpha of grade n is the (n-1)th triangular number.
"""
from copy import copy

from ...utils.concepts.dispatch import dispatch_on
from ..data_types import Alpha, MultiVector, Term


@dispatch_on(index=0)
def rev(obj):
    """
    Reverse the order basis elements within an object and then resolve
    back into permitted Alpha values.

    In notation, this is denoted with an over tilde (~).
    """
    raise NotImplementedError


@rev.add(Alpha)
def _rev_alpha(alpha):
    res = copy(alpha)
    if len(res._index) in [1, 4]:
        return res
    return -res


@rev.add(Term)
def _rev_term(term):
    res = copy(term)
    if len(res.index) in [1, 4] and res.index[0] != "p":
        return res
    return -res


@rev.add(MultiVector)
def _rev_multivector(mvec):
    return MultiVector([rev(t) for t in mvec._terms], cfg=mvec.cfg)
