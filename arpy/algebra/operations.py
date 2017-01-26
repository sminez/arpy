'''
Multiplying αs
==============
This is based on a set of simplification rules based on allowed
manipulations of elements in the algebra.
(NOTE:: In all notation, αμ.αν is simplified to αμν)

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
I am converting the current product into an array of integers in order
to allow for the different orderings of each final product in a flexible
way. Ordering is a mapping of index (0,1,2,3) to position in the final
product. This should be stable regardless of how we define the 16
elements of the algebra.

The algorithm makes use of the fact that for any ordering we can
dermine the whether the total number of pops is odd or even by looking
at the first element alone and then recursing on the rest of the
ordering as a sub-problem.
If the final position of the first element is even then it will take an
odd number of pops to correctly position it. We can then look only at
the remaining elements and re-label them with indices 1->(n-1) and
repeat the process until we are done.
NOTE:: I proved this by brute force. (i.e. listing all options and
showing that the proposition holds...!)
'''
from .ar_types import Alpha, Pair
from .config import ALLOWED, METRIC


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


def inverse(a):
    '''Find the inverse of an Alpha element'''
    return Alpha(a.index, (find_prod(a, a).sign * a.sign))


def wedge(a, b):
    '''Compute the Wedge Product of two values'''
    if isinstance(a, Alpha):
        if isinstance(b, Alpha):
            return find_prod(a, b)
        if isinstance(b, Pair):
            alpha = find_prod(a, b.alpha)
            return Pair(alpha, b.xi)
    elif isinstance(a, Pair):
        if isinstance(b, Alpha):
            alpha = find_prod(a.alpha, b)
            return Pair(alpha, a.xi)
        elif isinstance(b, Pair):
            raise ValueError("Unable to wedge Pairs")
    else:
        raise ValueError("Invalid wedge product: {} {}".format(a, b))


def div_by(a, b):
    '''Divide one element by another'''
    if isinstance(a, Alpha) and isinstance(b, Alpha):
        return find_prod(a, inverse(b))
    else:
        raise ValueError("Invalid division: {} / {}".format(a, b))


def div_into(a, b):
    '''Divide one element into another'''
    if isinstance(a, Alpha):
        if isinstance(b, Alpha):
            return find_prod(inverse(a), b)
        elif isinstance(b, Pair):
            alpha = find_prod(inverse(a), b.alpha)
            return Pair(alpha, b.xi)
    else:
        raise ValueError("Invalid division: {} \ {}".format(a, b))
