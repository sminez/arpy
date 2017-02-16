'''
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
from .config import ALLOWED, METRIC
from ..utils.concepts.dispatch import dispatch_on
from .ar_types import Alpha, Pair, MultiVector, XiProduct


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
        return Alpha(j.index, (i.sign * j.sign))
    elif j.index == 'p':
        return Alpha(i.index, (i.sign * j.sign))

    # Rule (2) :: Squaring and popping
    sign = i.sign * j.sign
    components = i.index + j.index
    intersection = set(i.index).intersection(set(j.index))

    for repeated in intersection:
        first = components.find(repeated)
        second = components.find(repeated, first + 1)
        n_pops = second - first - 1
        sign *= -1 if (n_pops % 2 == 1) else 1
        sign *= metric[repeated]
        components = ''.join(c for c in components if c != repeated)

    # If everything cancelled then i == j and we are left with αp (r-point)
    if len(components) == 0:
        return Alpha('p', sign)

    # Rule (3) :: Popping to the correct order
    target = targets[frozenset(components)]

    if target == components:
        return Alpha(target, sign)

    ordering = {c: i for i, c in enumerate(target)}
    current = [ordering[c] for c in components]

    while len(current) > 0:
        sign *= -1 if (current[0] % 2 == 0) else 1
        current = current[1:]
        new_order = {j: i for i, j in enumerate(sorted(current))}
        current = [new_order[k] for k in current]

    return Alpha(target, sign)


def inverse(a, metric=METRIC):
    '''Find the inverse of an Alpha element'''
    return Alpha(a.index, (find_prod(a, a, metric).sign * a.sign))


##############################################################################
@dispatch_on((0, 1))
def wedge(a, b, metric=METRIC):
    '''Compute the Wedge Product of two elements'''
    raise NotImplementedError


@wedge.add((Alpha, Alpha))
def _wedge_alpha_alpha(a, b, metric=METRIC):
    return find_prod(a, b, metric)


@wedge.add((Alpha, Pair))
def _wedge_alpha_pair(a, b, metric=METRIC):
    alpha = find_prod(a, b.alpha, metric)
    return Pair(alpha, b.xi)


@wedge.add((Pair, Alpha))
def _wedge_pair_alpha(a, b, metric=METRIC):
    alpha = find_prod(a.alpha, b, metric)
    return Pair(alpha, a.xi)


@wedge.add((Pair, Pair))
def _wedge_pair_pair(a, b, metric=METRIC):
    a, b = deepcopy(a), deepcopy(b)
    alpha = find_prod(a.alpha, b.alpha, metric)
    if alpha.sign == -1:
        axi, bxi = a.xi, b.xi
        axi.sign, bxi.sign = -1, -1
        alpha.sign = 1
        return Pair(alpha, XiProduct([axi, bxi]))
    else:
        return Pair(alpha, XiProduct([a.xi, b.xi]))


##############################################################################
@dispatch_on((0, 1))
def dot(a, b, metric=METRIC):
    '''Compute the dot product of two elements'''
    raise NotImplementedError


@dot.add((Alpha, Alpha))
def dot_alpha_alpha(a, b, metric=METRIC):
    raise NotImplementedError


@dot.add((Alpha, Pair))
def dot_alpha_pair(a, b, metric=METRIC):
    raise NotImplementedError


@dot.add((Pair, Alpha))
def dot_pair_alpha(a, b, metric=METRIC):
    raise NotImplementedError


##############################################################################
def full(a, b, metric=METRIC):
    '''Compute the Geometric product of two elements'''
    return MultiVector([dot(a, b), wedge(a, b)])


##############################################################################
@dispatch_on((0, 1))
def div_by(a, b, metric=METRIC):
    '''Divide one element by another'''
    raise NotImplementedError


@div_by.add((Alpha, Alpha))
def _div_by_alpha_alpha(a, b, metric=METRIC):
    return find_prod(a, inverse(b, metric), metric)


@div_by.add((Pair, Alpha))
def _div_by_pair_alpha(a, b, metric=METRIC):
    alpha = find_prod(inverse(a.alpha, metric), b, metric)
    return Pair(alpha, a.xi)


##############################################################################
@dispatch_on((0, 1))
def div_into(a, b, metric=METRIC):
    '''Divide one element into another'''
    raise NotImplementedError


@div_into.add((Alpha, Alpha))
def _div_into_Alpha_Alpha(a, b, metric=METRIC):
    return find_prod(inverse(a, metric), b, metric)


@div_into.add((Alpha, Pair))
def _div_into_Alpha_Pair(a, b, metric=METRIC):
    alpha = find_prod(inverse(a, metric), b.alpha, metric)
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
