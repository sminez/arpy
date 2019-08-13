"""
arpy (Absolute Relativity in Python)
Copyright (C) 2016-2018 Innes D. Anderson-Morrison All rights reserved.

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
from copy import deepcopy
from functools import wraps
from types import FunctionType
from .config import config as cfg
from .multivector import MultiVector
from ..utils.concepts.dispatch import dispatch_on
from .ar_types import Alpha, Pair, XiProduct


POINT = "p"


def product_cache(func):
    cache = dict()

    @wraps(func)
    def wrapped(i, j, cfg=cfg):
        args = (i, j, tuple(cfg.metric), tuple(cfg.allowed))
        result = cache.get(args)

        if result:
            return result

        result = func(i, j, cfg=cfg)
        cache[args] = result
        return result

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
    sign = i.sign * j.sign
    components = i.index + j.index

    # Multiplication by αp is idempotent
    if POINT in components:
        index = components.replace(POINT, "", 1)
        return Alpha(index, sign, cfg=cfg)

    # Pop and cancel matching components
    for repeated in set(i.index).intersection(set(j.index)):
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
    return Alpha(a.index, (find_prod(a, a, cfg).sign * a.sign), cfg=cfg)


##############################################################################
@dispatch_on((0, 1))
def full(a, b, cfg=cfg):
    """Compute the Full product of two elements"""
    raise NotImplementedError


@full.add((Alpha, Alpha))
def _full_alpha_alpha(a, b, cfg=cfg):
    return find_prod(a, b, cfg)


@full.add((Alpha, Pair))
def _full_alpha_pair(a, b, cfg=cfg):
    alpha = find_prod(a, b.alpha, cfg)
    return Pair(alpha, b.xi, cfg=cfg)


@full.add((Pair, Alpha))
def _full_pair_alpha(a, b, cfg=cfg):
    alpha = find_prod(a.alpha, b, cfg)
    return Pair(alpha, a.xi, cfg=cfg)


@full.add((Pair, Pair))
def _full_pair_pair(a, b, cfg=cfg):
    a, b = deepcopy(a), deepcopy(b)
    alpha = find_prod(a.alpha, b.alpha, cfg)
    return Pair(alpha, XiProduct([a.xi, b.xi]), cfg=cfg)
    # NOTE:: Not sure which method is correct as this is mixing
    #        usign with magsign...

    # if alpha.sign == -1:
    #     a.xi.sign *= -1
    #     b.xi.sign *= -1
    #     alpha.sign = 1
    #     return Pair(alpha, XiProduct([a.xi, b.xi]))
    # else:
    #     return Pair(alpha, XiProduct([a.xi, b.xi]))


@full.add((MultiVector, MultiVector))
def _full_mvec_mvec(mv1, mv2, cfg=cfg):
    prod = MultiVector((full(i, j, cfg) for i in mv1 for j in mv2), cfg=cfg)
    prod.replacements.extend(mv1.replacements + mv2.replacements)
    return prod


@full.add((Alpha, MultiVector))
def _full_alpha_mvec(a, m, cfg=cfg):
    prod = MultiVector((full(a, comp, cfg) for comp in m), cfg=cfg)
    prod.replacements = m.replacements
    return prod


@full.add((MultiVector, Alpha))
def _full_mvec_alpha(m, a, cfg=cfg):
    prod = MultiVector((full(comp, a, cfg) for comp in m), cfg=cfg)
    prod.replacements = m.replacements
    return prod


# NOTE:: Definitions of the full product involving differnetials are found in
#        differential.py due to import conflicts.


##############################################################################
@dispatch_on((0, 1))
def div_by(a, b, cfg=cfg):
    """Divide one element by another"""
    raise NotImplementedError


@div_by.add((Alpha, Alpha))
def _div_by_alpha_alpha(a, b, cfg=cfg):
    return find_prod(a, inverse(b, cfg), cfg)


@div_by.add((Pair, Alpha))
def _div_by_pair_alpha(a, b, cfg=cfg):
    alpha = find_prod(a.alpha, inverse(b, cfg), cfg)
    return Pair(alpha, a.xi, cfg=cfg)


##############################################################################
@dispatch_on((0, 1))
def div_into(a, b, cfg=cfg):
    """Divide one element into another"""
    raise NotImplementedError


@div_into.add((Alpha, Alpha))
def _div_into_Alpha_Alpha(a, b, cfg=cfg):
    return find_prod(inverse(a, cfg), b, cfg)


@div_into.add((Alpha, Pair))
def _div_into_Alpha_Pair(a, b, cfg=cfg):
    alpha = find_prod(inverse(a, cfg), b.alpha, cfg)
    return Pair(alpha, b.xi, cfg=cfg)


##############################################################################
@dispatch_on(index=0)
def project(element, grade, cfg=cfg):
    """
    Implementation of the grade-projection operator <A>n.
    Return only the elements of A that are of grade n.
    NOTE:: αp is a grade-0 scalar element.
    """
    raise NotImplementedError


@project.add(Alpha)
def _project_alpha(element, grade, cfg=cfg):
    if element.index == POINT:
        if grade == 0:
            return element
        else:
            return None
    elif len(element.index) == grade:
        return element
    else:
        return None


@project.add(Pair)
def _project_pair(element, grade, cfg=cfg):
    if element.alpha.index == POINT:
        if grade == 0:
            return element
        else:
            return None
    elif len(element.alpha.index) == grade:
        return element
    else:
        return None


@project.add(MultiVector)
def _project_multivector(element, grade, cfg=cfg):
    correct_grade = []
    if grade == 0:
        for component in element:
            if component.alpha.index == POINT:
                correct_grade.append(component)
    else:
        for component in element:
            ix = component.alpha.index
            if len(ix) == grade and ix != POINT:
                correct_grade.append(component)
    res = MultiVector(correct_grade, cfg=cfg)
    res.replacements.extend(element.replacements)
    return res


##############################################################################


@dispatch_on((0, 1, 2))
def prod_apply(arg1, arg2, arg3=None, cfg=cfg):
    """
    Apply a function to the cartesian product of two multivectors
    NOTE:: The function must act on a single Multivector.
    """
    raise NotImplementedError


@prod_apply.add((FunctionType, MultiVector, MultiVector))
def _prod_apply_dmm(func, mv1, mv2, cfg=cfg):
    if not isinstance(mv1, MultiVector) and isinstance(mv2, MultiVector):
        raise TypeError("Arguments must be a MultiVectors")
    return MultiVector((full(i, j, cfg) for i in func(mv1, cfg=cfg) for j in mv2), cfg=cfg)


@prod_apply.add((MultiVector, FunctionType, MultiVector))
def _prod_apply_mdm(mv1, func, mv2, cfg=cfg):
    if not isinstance(mv1, MultiVector) and isinstance(mv2, MultiVector):
        raise TypeError("Arguments must be a MultiVectors")
    return MultiVector((full(i, j, cfg) for i in mv1 for j in func(mv2, cfg=cfg)), cfg=cfg)


@prod_apply.add((FunctionType, tuple, type(None)))
def _prod_apply_d_mm(func, mvecs=(None, None), _=None, cfg=cfg):
    if not isinstance(mvecs[0], MultiVector) and isinstance(mvecs[1], MultiVector):
        raise TypeError("Arguments must be a MultiVectors")
    return func(full(mvecs[0], mvecs[1], cfg))


##############################################################################


@dispatch_on(index=0)
def dagger(obj, cfg=cfg):
    """
    Compute the Hermitian conjugate of the argument.
    """
    raise NotImplementedError


@dagger.add(MultiVector)
def _dagger_mvec(mvec, cfg=cfg):
    _neg = [
        Alpha(a, cfg=cfg)
        for a in cfg.allowed
        if full(Alpha(a, cfg=cfg), Alpha(a, cfg=cfg), cfg).sign == -1
    ]
    mvec = deepcopy(mvec)
    new_vec = []
    for pair in mvec:
        if pair.alpha in _neg:
            pair.alpha.sign *= -1
        new_vec.append(pair)
    res = MultiVector(new_vec, cfg=cfg)
    res.replacements.extend(mvec.replacements)
    return res


@dagger.add(Alpha)
def _dagger_alpha(alpha, cfg=cfg):
    res = deepcopy(alpha)
    if full(alpha, alpha, cfg).sign == -1:
        res.sign *= -1

    return res


@dagger.add(Pair)
def _dagger_pair(pair, cfg=cfg):
    res = deepcopy(pair)
    if full(pair.extract_alpha(), pair.extract_alpha(), cfg).sign == -1:
        res.xi.sign *= -1

    return res


##############################################################################


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


# @commutator.add((Pair, Pair))
# def _ring_commutator(a, b):
#     '''Computes the ring commutator [a, b] = ab - ba for Pairs'''
#     return full(a, b) - full(b, a)
