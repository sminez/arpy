'''
Classes, functions and default representations of Ξ vectors in the algebra.

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
from ..config import ALLOWED, ALLOWED_GROUPS, METRIC, DIVISION_TYPE


class Alpha():
    '''Unit elements in the algebra'''
    def __init__(self, index, sign=None):
        '''
        Handle multiple constructor methods for αs
        '''
        if sign is None:
            if index.startswith('-'):
                index, sign = index[1:], -1
            else:
                sign = 1

        if index not in ALLOWED and index not in ALLOWED_GROUPS:
            raise ValueError('Invalid α index: {}'.format(index))

        if sign not in [1, -1]:
            raise ValueError('Invalid α sign: {}'.format(sign))

        self.index = index
        self.sign = sign

    def __repr__(self):
        neg = '-' if self.sign == -1 else ''
        return '{}α{}'.format(neg, self.index)

    def __eq__(self, other):
        return (self.index == other.index) and (self.sign == other.sign)

    def __mul__(self, other):
        if isinstance(other, Alpha):
            return find_prod(self, other)
        if isinstance(other, Pair):
            alpha = self * other.alpha
            return Pair(alpha, other.xi)

    def __neg__(self):
        self.sign *= -1
        return self

    def inverse(self):
        return Alpha(self.index, ((self * self).sign * self.sign))

    def __truediv__(self, other):
        if DIVISION_TYPE == 'by':
            if isinstance(other, Alpha):
                return self * other.inverse()
            elif isinstance(other, Pair):
                raise TypeError('undefined division')
        elif DIVISION_TYPE == 'into':
            if isinstance(other, Alpha):
                return self.inverse() * other
            elif isinstance(other, Pair):
                alpha = self.inverse() * other.alpha
                return Pair(alpha, other.xi)
        else:
            raise ValueError('Invalid division type: ' + DIVISION_TYPE)


class Xi():
    '''A symbolic Real value'''
    def __init__(self, val, unit=None, partials=[]):
        if isinstance(val, Alpha):
            val = val.index

        if unit is None:
            unit = val

        self.val = val
        self.unit = unit
        self.partials = partials

    def __eq__(self, other):
        return (self.val == other.val) and (self.partials == other.partials)

    def __repr__(self):
        partials = ('∂{}'.format(p.index) for p in reversed(self.partials))
        return '{}ξ{}'.format(''.join(partials), self.val)


class Pair:
    def __init__(self, a, xi=None):
        if xi is None:
            xi = Xi(a)

        if isinstance(a, Alpha):
            self.alpha = a
        else:
            self.alpha = Alpha(a)
        self.xi = xi

    def __eq__(self, other):
        return (self.alpha == other.alpha) and (self.xi == other.xi)

    def __repr__(self):
        return '{} {}'.format(self.alpha, self.xi)

    def __mul__(self, other):
        if isinstance(other, Alpha):
            alpha = self.alpha * other
            return Pair(alpha, self.xi)
        elif isinstance(other, Pair):
            raise NotImplemented
            # alpha = self.alpha * other.alpha
            # xi = self.xi * other.xi
            # return ξα(alpha, xi)

    def __truediv__(self, other):
        if DIVISION_TYPE == 'by':
            if isinstance(other, Alpha):
                alpha = self.alpha * other.invert()
                return Pair(alpha, self.xi)
            elif isinstance(other, Pair):
                raise TypeError('undefined division')
        elif DIVISION_TYPE == 'into':
            raise TypeError('undefined division')
        else:
            raise ValueError('Invalid division type: ' + DIVISION_TYPE)


##############################################################################
def a(index, sign=None):
    return Alpha(index, sign)


def xi(val, unit=None, partials=[]):
    return Xi(val, unit, partials)


# Prebuild vectors to work with based on the 4-vector components
XiM = [Pair('p')] + [Pair(a) for a in ALLOWED if len(a) == 2 and '0' not in a]
XiT = [Pair('0')] + [Pair(a) for a in ALLOWED if len(a) == 3 and '0' in a]
XiA = [Pair('123')] + [Pair(a) for a in ALLOWED if len(a) == 1 and a != '0']
XiE = [Pair('0123')] + [Pair(a) for a in ALLOWED if len(a) == 2 and '0' in a]
XiG = [Pair(a) for a in ALLOWED]
# Prebuild vectors based on length of index
Xi1 = [Pair(a) for a in ALLOWED if len(a) == 1]
Xi2 = [Pair(a) for a in ALLOWED if len(a) == 2]
Xi3 = [Pair(a) for a in ALLOWED if len(a) == 3]


##############################################################################


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
        return a(j.index, (i.sign * j.sign))
    elif j.index == 'p':
        return a(i.index, (i.sign * j.sign))

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
        return a('p', sign)

    # Rule (3) :: Popping to the correct order
    target = targets[frozenset(components)]

    if target == components:
        return a(target, sign)

    ordering = {c: i for i, c in enumerate(target)}
    current = [ordering[c] for c in components]

    while len(current) > 0:
        sign *= -1 if (current[0] % 2 == 0) else 1
        current = current[1:]
        new_order = {j: i for i, j in enumerate(sorted(current))}
        current = [new_order[k] for k in current]

    return a(target, sign)
