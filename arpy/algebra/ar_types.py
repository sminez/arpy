'''
arpy (Absolute Relativity in Python)
Copyright (C) 2016-2017 Innes D. Anderson-Morrison All rights reserved.

Classes, functions and default representations of Ξ vectors in the algebra.
'''
from .config import ALLOWED, ALLOWED_GROUPS, SUB_SCRIPTS


class Alpha:
    '''Unit elements in the algebra'''
    def __init__(self, index, sign=None,
                 allowed=ALLOWED, allowed_groups=ALLOWED_GROUPS):
        '''
        Handle multiple constructor methods for αs
        '''
        self.allowed = ALLOWED
        self.allowed_groups = ALLOWED_GROUPS

        if sign is None:
            if index.startswith('-'):
                index, sign = index[1:], -1
            else:
                sign = 1

        if index not in allowed and index not in allowed_groups:
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

    def __lt__(self, other):
        return self.allowed.index(self.index) < self.allowed.index(other.index)

    def __neg__(self):
        self.sign *= -1
        return self

    def __hash__(self):
        return hash((self.index, self.sign))


class Xi:
    '''A symbolic Real value'''
    def __init__(self, val, partials=None, sign=1):
        if isinstance(val, Alpha):
            val = val.index

        self.val = val
        self.sign = sign
        self.partials = partials if partials else []

    @property
    def components(self):
        return [self]

    def __eq__(self, other):
        return (self.val == other.val) and (self.partials == other.partials)

    def __lt__(self, other):
        try:
            return ALLOWED.index(self.val) < ALLOWED.index(other.val)
        except:
            return self.val < other.val

    def __repr__(self):
        sign = '+' if self.sign == 1 else '-'
        partials = (
            '∂{}'.format(''.join(SUB_SCRIPTS[i] for i in p.index))
            for p in reversed(self.partials)
        )
        if self.val in ALLOWED + ALLOWED_GROUPS:
            display_val = ''.join(SUB_SCRIPTS[i] for i in self.val)
            return '{}{}ξ{}'.format(sign, ''.join(partials), display_val)
        else:
            return '{}{}{}'.format(sign, ''.join(partials), self.val)


class XiProduct:
    '''Symbolic Xi valued products with a single sign'''
    def __init__(self, components):
        self.components = tuple(components)
        self.partials = []
        self.sign_base = 1

    @property
    def sign(self):
        s = self.sign_base
        for comp in self.components:
            s *= comp.sign
        return s

    @sign.setter
    def sign(self, val):
        self.sign_base = val

    @property
    def val(self):
        # Expressing the product values as a dotted list of indices
        return '.'.join(c.val for c in self.components)

    def __eq__(self, other):
        same_sign = (self.sign == other.sign)
        same_components = (self.components == other.components)
        return same_sign and same_components

    def __repr__(self):
        sign = '+' if self.sign == 1 else '-'
        # Stripping component signs as we have taken care of the overall
        # product sign at initialisation.
        partials = (
            '∂{}'.format(''.join(SUB_SCRIPTS[i] for i in p.index))
            for p in reversed(self.partials)
        )
        comps = ''.join(str(c)[1:] for c in self.components)
        return sign + ''.join(partials) + comps


class Pair:
    '''A Pair may be any object along with an Alpha value'''
    def __init__(self, a, x=None):
        if x is None:
            x = Xi(a)

        if isinstance(a, Alpha):
            self.alpha = a
        else:
            self.alpha = Alpha(a)

        if not isinstance(x, (Xi, XiProduct)):
            x = Xi(x)

        self.xi = x

    def __eq__(self, other):
        return (self.alpha == other.alpha) and (self.xi == other.xi)

    def __repr__(self):
        return '({}, {})'.format(self.alpha, self.xi)
