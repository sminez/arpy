'''
Classes, functions and default representations of Ξ vectors in the algebra.
'''
import collections.abc
from .config import ALLOWED, ALLOWED_GROUPS, XI_GROUPS


class Alpha:
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

    def __hash__(self):
        return hash((self.index, self.sign))


class Xi:
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


class MultiVector(collections.abc.Set):
    '''A custom container type for working efficiently with multivectors'''
    def __init__(self, components=[]):
        '''
        Given a list of pairs, build the mulitvector by binding the ξ values
        '''
        self.basis_blades = {Alpha(a): None for a in ALLOWED}

        if not all([isinstance(comp, Pair) for comp in components]):
            raise ValueError("Multivectors can only contain Xi/Alpha pairs")

        for comp in components:
            self.basis_blades[comp.alpha] = comp.xi

    def __len__(self):
        '''
        Only initialised blades count towards the length of a multivector
        '''
        return len([blade for blade in self.basis_blades if blade is not None])

    def __contains__(self, other):
        '''Return True if the requested alpha value has been initialised'''
        if isinstance(other, Alpha):
            return self.basis_blades[other] is not None
        elif isinstance(other, Pair):
            return self.basis_blades[other.alpha] is not None
        else:
            return False

    def __getitem__(self, key):
        '''mvec[alpha] returns a pair'''
        if isinstance(key, str):
            # Allow retreval by bare string as well as Alpha
            key = Alpha(key)
        if not isinstance(key, Alpha):
            raise KeyError

        xi = self.basis_blades[key]
        if xi is not None:
            return Pair(key, xi)
        else:
            raise KeyError

    def __iter__(self):
        # '''Iteration over a multivector always gives all 16 elements'''
        for a in ALLOWED:
            blade = Alpha(a)
            if self.basis_blades[blade] is not None:
                yield Pair(blade, self.basis_blades[blade])

    def __repr__(self):
        comps = ['α{}{}'.format(a, self.basis_blades[Alpha(a)])
                 for a in ALLOWED if self.basis_blades[Alpha(a)]]
        return '{' + ', '.join(comps) + '}'
