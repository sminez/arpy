'''
Classes, functions and default representations of Ξ vectors in the algebra.
'''
from .config import ALLOWED, ALLOWED_GROUPS


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

    def __neg__(self):
        self.sign *= -1
        return self


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
        return '({},{})'.format(self.alpha, self.xi)


##############################################################################
# vectors to work with based on the 4-vector components
XiM = [Pair('p')] + [Pair(a) for a in ALLOWED if len(a) == 2 and '0' not in a]
XiT = [Pair('0')] + [Pair(a) for a in ALLOWED if len(a) == 3 and '0' in a]
XiA = [Pair('123')] + [Pair(a) for a in ALLOWED if len(a) == 1 and a != '0']
XiE = [Pair('0123')] + [Pair(a) for a in ALLOWED if len(a) == 2 and '0' in a]
XiG = [Pair(a) for a in ALLOWED]
# Prebuild vectors based on length of index
Xi1 = [Pair(a) for a in ALLOWED if len(a) == 1]
Xi2 = [Pair(a) for a in ALLOWED if len(a) == 2]
Xi3 = [Pair(a) for a in ALLOWED if len(a) == 3]

# For use in the REPL
Xi_vecs = {'XiM': XiM, 'XiT': XiT, 'XiA': XiA, 'XiE': XiE, 'XiG': XiG,
           'Xi1': Xi1, 'Xi2': Xi2, 'Xi3': Xi3}
