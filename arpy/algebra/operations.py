'''
arpy (Absolute Relativity in Python)
Copyright (C) 2016-2017 Innes D. Anderson-Morrison All rights reserved.

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
'''
from copy import deepcopy
from types import FunctionType
from .config import ALLOWED, METRIC
from .multivector import MultiVector
from ..utils.concepts.dispatch import dispatch_on
from .ar_types import Alpha, Pair, XiProduct


def find_prod(i, j, metric=METRIC, allowed=ALLOWED):
    '''
    Compute the product of two alpha values in the algebra. This uses some
    optimisations and observations that I've made in order to speed up the
    computation.
    '''
    # NOTE:: These are here so that they the user can change metric and allowed
    #        without having to reload the module.
    metric = {k: v for k, v in zip('0123', metric)}
    targets = {frozenset(a): a for a in allowed}

    # Rule (1) :: Multiplication by αp is idempotent
    if i.index == 'p':
        return Alpha(j.index, (i.sign * j.sign), allowed=allowed)
    elif j.index == 'p':
        return Alpha(i.index, (i.sign * j.sign), allowed=allowed)

    # Rule (2) :: Squaring and popping
    sign = i.sign * j.sign
    components = i.index + j.index
    intersection = set(i.index).intersection(set(j.index))

    for repeated in intersection:
        first = components.find(repeated)
        second = components.find(repeated, first + 1)
        n_pops = second - first - 1
        sign *= (-1 if (n_pops % 2 == 1) else 1)
        sign *= metric[repeated]
        components = ''.join(c for c in components if c != repeated)

    # If everything cancelled then i == j and we are left with αp (r-point)
    if len(components) == 0:
        return Alpha('p', sign)

    # Rule (3) :: Popping to the correct order
    target = targets[frozenset(components)]

    if target == components:
        return Alpha(target, sign, allowed=allowed)

    ordering = {c: i+1 for i, c in enumerate(target)}
    current = [ordering[c] for c in components]

    while len(current) > 1:
        if current[0] % 2 == 0:
            sign *= -1
        current = current[1:]
        new_order = {j: i+1 for i, j in enumerate(sorted(current))}
        current = [new_order[k] for k in current]

    return Alpha(target, sign, allowed=allowed)


def inverse(a, metric=METRIC, allowed=ALLOWED):
    '''Find the inverse of an Alpha element'''
    return Alpha(
        a.index,
        (find_prod(a, a, metric, allowed).sign * a.sign),
        allowed=allowed
    )


##############################################################################
@dispatch_on((0, 1))
def full(a, b, metric=METRIC, allowed=ALLOWED):
    '''Compute the Full product of two elements'''
    raise NotImplementedError


@full.add((Alpha, Alpha))
def _full_alpha_alpha(a, b, metric=METRIC, allowed=ALLOWED):
    return find_prod(a, b, metric, allowed)


@full.add((Alpha, Pair))
def _full_alpha_pair(a, b, metric=METRIC, allowed=ALLOWED):
    alpha = find_prod(a, b.alpha, metric, allowed)
    return Pair(alpha, b.xi)


@full.add((Pair, Alpha))
def _full_pair_alpha(a, b, metric=METRIC, allowed=ALLOWED):
    alpha = find_prod(a.alpha, b, metric, allowed)
    return Pair(alpha, a.xi)


@full.add((Pair, Pair))
def _full_pair_pair(a, b, metric=METRIC, allowed=ALLOWED):
    a, b = deepcopy(a), deepcopy(b)
    alpha = find_prod(a.alpha, b.alpha, metric, allowed)
    return Pair(alpha, XiProduct([a.xi, b.xi]))
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
def _full_mvec_mvec(mv1, mv2, metric=METRIC, allowed=ALLOWED):
    prod = MultiVector(full(i, j, metric, allowed) for i in mv1 for j in mv2)
    prod._simplify()
    return prod

# NOTE:: Definitions of the full product involving differnetials are found in
#        differential.py due to import conflicts.


##############################################################################
@dispatch_on((0, 1))
def div_by(a, b, metric=METRIC, allowed=ALLOWED):
    '''Divide one element by another'''
    raise NotImplementedError


@div_by.add((Alpha, Alpha))
def _div_by_alpha_alpha(a, b, metric=METRIC, allowed=ALLOWED):
    return find_prod(a, inverse(b, metric, allowed), metric, allowed)


@div_by.add((Pair, Alpha))
def _div_by_pair_alpha(a, b, metric=METRIC, allowed=ALLOWED):
    alpha = find_prod(a.alpha, inverse(b, metric, allowed), metric, allowed)
    return Pair(alpha, a.xi)


##############################################################################
@dispatch_on((0, 1))
def div_into(a, b, metric=METRIC, allowed=ALLOWED):
    '''Divide one element into another'''
    raise NotImplementedError


@div_into.add((Alpha, Alpha))
def _div_into_Alpha_Alpha(a, b, metric=METRIC, allowed=ALLOWED):
    return find_prod(inverse(a, metric, allowed), b, metric, allowed)


@div_into.add((Alpha, Pair))
def _div_into_Alpha_Pair(a, b, metric=METRIC, allowed=ALLOWED):
    alpha = find_prod(inverse(a, metric, allowed), b.alpha, metric, allowed)
    return Pair(alpha, b.xi)


##############################################################################
@dispatch_on(index=0)
def project(element, grade):
    '''
    Implementation of the grade-projection operator <A>n.
    Return only the elements of A that are of grade n.
    NOTE:: αp is a grade-0 scalar element.
    '''
    raise NotImplementedError


@project.add(Alpha)
def _project_alpha(element, grade):
    if grade == 0 and element.index == 'p':
        return element
    elif len(element.index) == grade:
        return element
    else:
        return None


@project.add(Pair)
def _project_pair(element, grade):
    if grade == 0 and element.alpha.index == 'p':
        return element
    elif len(element.alpha.index) == grade:
        return element
    else:
        return None


@project.add(MultiVector)
def _project_multivector(element, grade):
    correct_grade = []
    if grade == 0:
        for component in element:
            if component.alpha.index == 'p':
                correct_grade.append(component)
    else:
        for component in element:
            ix = component.alpha.index
            if len(ix) == grade and ix != 'p':
                correct_grade.append(component)
    return MultiVector(correct_grade)


##############################################################################

@dispatch_on((0, 1, 2))
def prod_apply(arg1, arg2, arg3=None, metric=METRIC, allowed=ALLOWED):
    '''
    Apply a function to the cartesian product of two multivectors
    NOTE:: The function must act on a single Multivector.
    '''
    raise NotImplementedError


@prod_apply.add((FunctionType, MultiVector, MultiVector))
def _prod_apply_dmm(func, mv1, mv2, metric=METRIC, allowed=ALLOWED):
    if not isinstance(mv1, MultiVector) and isinstance(mv2, MultiVector):
        raise TypeError('Arguments must be a MultiVectors')
    return MultiVector(
        full(i, j, metric, allowed)
        for i in func(mv1, metric=metric, allowed=allowed)
        for j in mv2
    )


@prod_apply.add((MultiVector, FunctionType, MultiVector))
def _prod_apply_mdm(mv1, func, mv2, metric=METRIC, allowed=ALLOWED):
    if not isinstance(mv1, MultiVector) and isinstance(mv2, MultiVector):
        raise TypeError('Arguments must be a MultiVectors')
    return MultiVector(
        full(i, j, metric, allowed)
        for i in mv1
        for j in func(mv2, metric=metric, allowed=allowed)
    )


@prod_apply.add((FunctionType, tuple, type(None)))
def _prod_apply_d_mm(func, mvecs=(None, None), _=None,
                     metric=METRIC, allowed=ALLOWED):
    if not isinstance(mvecs[0], MultiVector) and \
            isinstance(mvecs[1], MultiVector):
        raise TypeError('Arguments must be a MultiVectors')
    return func(full(mvecs[0], mvecs[1], metric, allowed))


##############################################################################


def dagger(mvec, metric=METRIC, allowed=ALLOWED):
    '''return ther Hermitian conjugate of the Multivector'''
    _neg = [
        Alpha(a) for a in allowed
        if full(Alpha(a), Alpha(a), metric, allowed).sign == -1
    ]
    mvec = deepcopy(mvec)
    new_vec = []
    for pair in mvec:
        if pair.alpha in _neg:
            pair.alpha.sign *= -1
        new_vec.append(pair)
    return MultiVector(new_vec)


##############################################################################

@dispatch_on((0, 1))
def commutator(a, b, metric=METRIC, allowed=ALLOWED):
    '''
    Computes the group commutator [a, b] = (a . b . a^-1 . b^-1) for Alphas.
    '''
    raise NotImplementedError


@commutator.add((Alpha, Alpha))
def _group_commutator(a, b, metric=METRIC, allowed=ALLOWED):
    product = full(a, b, metric, allowed)
    product = full(product, inverse(a), metric, allowed)
    product = full(product, inverse(b), metric, allowed)
    return product


# @commutator.add((Pair, Pair))
# def _ring_commutator(a, b):
#     '''Computes the ring commutator [a, b] = ab - ba for Pairs'''
#     return full(a, b) - full(b, a)
