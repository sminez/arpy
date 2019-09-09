"""
Multiplying αs
==============
This is based on a set of simplification rules based on allowed manipulations
of elements in the algebra. (NOTE: In all notation, αμ.αν is simplified to αμν)

(1)   αpμ == αμp == αμ
    'Multiplication by αp (r-point) is idempotent. (αp is the identity)'
(2i)  α0^2 == αp
    'Repeated α0 indices can just be removed.'
(2ii) αi^2 == -αp
    'Repeated αi indices can be removed by negating'
(2iii) α^2 == +-αp
    'All elements square to either +αp or -αp'
(3)   αμν == -ανμ
    'Adjacent indices can be popped by negating.'


Counting pops
=============
I am converting the current product into an array of integers in order to allow
for the different orderings of each final product in a flexible way. Ordering
is a mapping of index (0,1,2,3) to position in the final product. This should
be stable regardless of how we define the 16 elements of the algebra.

The algorithm makes use of the fact that for any ordering we can dermine the
whether the total number of pops is odd or even by looking at the first element
alone and then recursing on the rest of the ordering as a sub-problem.
If the final position of the first element is even then it will take an odd
number of pops to correctly position it. We can then look only at the remaining
elements and re-label them with indices 1->(n-1) and repeat the process until
we are done.
"""
from copy import copy
from functools import wraps

from ...utils.concepts.dispatch import dispatch_on
from ..config import config as cfg
from ..data_types import Alpha, MultiVector, Term

POINT = "p"


def product_cache(func):
    cache = dict()

    @wraps(func)
    def wrapped(i, j, cfg=cfg):
        args = (i, j, tuple(cfg.metric), tuple(cfg.allowed))
        result = cache.get(args)

        if result:
            return copy(result)

        result = func(i, j, cfg=cfg)
        cache[args] = result
        return copy(result)

    wrapped.cache = cache
    return wrapped


@product_cache
def find_prod(i, j, cfg=cfg):
    """
    Compute the product of two alpha values in the algebra. This uses some
    optimisations and observations that I've made in order to speed up the
    computation.

    NOTE: find_prod ALWAYS returns a new alpha as we don't want to mutate
          the values passed in as that will mess up any future calculations!
    """
    if not (i.allowed == j.allowed == cfg.allowed):
        err = "Inconsistant allowed values detected when computing a product.\n"
        err += 'Config allowed: {}\nPassed values: "{}" "{}"'.format(cfg.allowed, i, j)
        raise ValueError(err)

    metric = {k: v for k, v in zip("0123", cfg.metric)}
    targets = {frozenset(a): a for a in cfg.allowed}
    sign = i._sign * j._sign
    components = i._index + j._index

    # Multiplication by αp is idempotent
    if POINT in components:
        index = components.replace(POINT, "", 1)
        return Alpha(index, sign, cfg=cfg)

    # Pop and cancel matching components
    for repeated in set(i._index).intersection(set(j._index)):
        first = components.find(repeated)
        second = components.find(repeated, first + 1)
        n_pops = second - first - 1
        sign *= -1 if (n_pops % 2 == 1) else 1
        sign *= metric[repeated]
        components = "".join(c for c in components if c != repeated)

    if len(components) == 0:
        return Alpha(POINT, sign, cfg=cfg)

    target = targets[frozenset(components)]

    if target == components:
        return Alpha(target, sign, cfg=cfg)

    ordering = {c: i + 1 for i, c in enumerate(target)}
    current = [ordering[c] for c in components]

    while len(current) > 1:
        if current[0] % 2 == 0:
            sign *= -1
        current = current[1:]
        new_order = {j: i + 1 for i, j in enumerate(sorted(current))}
        current = [new_order[k] for k in current]

    return Alpha(target, sign, cfg=cfg)


def inverse(a, cfg=cfg):
    """Find the inverse of an Alpha element"""
    return Alpha(a._index, (find_prod(a, a, cfg)._sign * a._sign), cfg=cfg)


@dispatch_on((0, 1))
def full(a, b, cfg=cfg):
    """Compute the Full product of two elements"""
    raise NotImplementedError


@full.add((Alpha, Alpha))
def _full_alpha_alpha(a, b, cfg=cfg):
    return find_prod(a, b, cfg)


@full.add((Alpha, Term))
def _full_alpha_term(a, b, cfg=cfg):
    alpha = find_prod(a, b.alpha, cfg)
    return Term(alpha, b._components, cfg=cfg)


@full.add((Term, Alpha))
def _full_term_alpha(a, b, cfg=cfg):
    alpha = find_prod(a.alpha, b, cfg)
    return Term(alpha, a._components, cfg=cfg)


@full.add((Term, Term))
def _full_term_term(a, b, cfg=cfg):
    a, b = copy(a), copy(b)
    alpha = find_prod(a.alpha, b.alpha, cfg)

    # TODO: This will be incorrect for products of terms that already contain
    #       products of terms that have their own partials. We need to sort out
    #       the correct behaviour in this case in order to be able to perform more
    #       complicated products.
    return Term(alpha, a._components + b._components, cfg=cfg)


@full.add((MultiVector, MultiVector))
def _full_mvec_mvec(mv1, mv2, cfg=cfg):
    prod = MultiVector((full(i, j, cfg) for i in mv1 for j in mv2), cfg=cfg)
    return prod


@full.add((Alpha, MultiVector))
def _full_alpha_mvec(a, m, cfg=cfg):
    prod = MultiVector((full(a, comp, cfg) for comp in m), cfg=cfg)
    return prod


@full.add((MultiVector, Alpha))
def _full_mvec_alpha(m, a, cfg=cfg):
    prod = MultiVector((full(comp, a, cfg) for comp in m), cfg=cfg)
    return prod


# NOTE:: Definitions of the full product involving differnetials are found in
#        differential.py due to import conflicts.
