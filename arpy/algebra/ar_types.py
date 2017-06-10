'''
arpy (Absolute Relativity in Python)
Copyright (C) 2016-2017 Innes D. Anderson-Morrison All rights reserved.
'''
from copy import deepcopy
from .config import ALLOWED, ALLOWED_GROUPS, SUB_SCRIPTS, BXYZ_LIKE


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
        try:
            return '{}α{}'.format(
                neg, ''.join(SUB_SCRIPTS[i] for i in self.index)
            )
        except:
            return '{}α{}'.format(neg, self.index)

    def __tex__(self):
        neg = '-' if self.sign == -1 else ''
        return neg + '\\alpha_{' + self.index + '}'

    def __eq__(self, other):
        if not isinstance(other, Alpha):
            return False
        return (self.index == other.index) and (self.sign == other.sign)

    def __lt__(self, other):
        return self.allowed.index(self.index) < self.allowed.index(other.index)

    def __neg__(self):
        neg = deepcopy(self)
        neg.sign *= -1
        return neg

    def __hash__(self):
        return hash((self.index, self.sign))


class Xi:
    '''A symbolic Real value'''
    def __init__(self, val, partials=None, sign=1):
        if isinstance(val, Alpha):
            val = val.index

        if isinstance(val, str) and val.startswith('-'):
            self.sign = sign * -1
            self.val = val[1:]
        else:
            self.val = val
            self.sign = sign

        self.partials = partials if partials else []

    def __hash__(self):
        return hash((self.val, self.sign, tuple(self.partials)))

    @property
    def components(self):
        return [self]

    def __eq__(self, other):
        if not isinstance(other, Xi):
            return False

        return (self.val == other.val) and \
               (self.partials == other.partials) and \
               (self.sign == other.sign)

    def __lt__(self, other):
        try:
            return ALLOWED.index(self.val) < ALLOWED.index(other.val)
        except:
            return self.val < other.val

    def __neg__(self):
        neg = deepcopy(self)
        neg.sign *= -1
        return neg

    def __repr__(self):
        sign = '' if self.sign == 1 else '-'
        partials = (
            '∂{}'.format(''.join(SUB_SCRIPTS[i] for i in p.index))
            for p in reversed(self.partials)
        )
        if self.val in ALLOWED + ALLOWED_GROUPS:
            display_val = ''.join(SUB_SCRIPTS[i] for i in self.val)
            return '{}{}ξ{}'.format(sign, ''.join(partials), display_val)
        else:
            return '{}{}{}'.format(sign, ''.join(partials), self.val)

    def __tex__(self):
        sign = '+' if self.sign == 1 else '-'
        partials = ''.join(
            '\\partial_{' + p.index + '}' for p in reversed(self.partials)
        )
        if self.val in ALLOWED + ALLOWED_GROUPS:
            return sign + partials + '\\xi_{' + self.val + '}'
        else:
            return sign + partials + self.val

    def bxyz(self):
        '''Return a string representing only {b,x,y,z} information'''
        sign = '+' if self.sign == 1 else '-'
        partials = [BXYZ_LIKE[p.index] for p in self.partials]
        partial_str = ''.join(['∂{}'.format(p) for p in reversed(partials)])
        val = BXYZ_LIKE[self.val]
        return sign + partial_str + '[' + val + ']'


class XiProduct:
    '''Symbolic Xi valued products with a single sign'''
    def __init__(self, components):
        self.components = tuple(components)
        self.partials = []
        self.sign_base = 1

    def __hash__(self):
        return hash((self.components, tuple(self.partials), self.sign_base))

    @property
    def sign(self):
        s = 1
        for comp in self.components:
            s *= comp.sign
        return s * self.sign_base

    @sign.setter
    def sign(self, val):
        if val not in [1, -1]:
            raise ValueError()
        if self.sign == -1:
            val = -val
        self.sign_base = val

    @property
    def val(self):
        # ordered tuple of indices
        return tuple(c.val for c in self.components)

    def __eq__(self, other):
        if not isinstance(other, XiProduct):
            return False

        same_sign = (self.sign == other.sign)
        same_partials = self.partials == other.partials
        self_comps = [c for c in deepcopy(self.components)]
        for c in self_comps:
            c.sign = 1
        other_comps = [c for c in deepcopy(other.components)]
        for c in other_comps:
            c.sign = 1
        same_components = (set(self_comps) == set(other_comps))
        return same_sign and same_partials and same_components

    def __neg__(self):
        neg = deepcopy(self)
        neg.sign_base *= -1
        return neg

    def __repr__(self):
        partials = (
            '∂{}'.format(''.join(SUB_SCRIPTS[i] for i in p.index))
            for p in reversed(self.partials)
        )
        comps = [str(c) for c in self.components]
        comps = [c[1:] if c[0] == '-' else c for c in comps]
        sign = '' if self.sign == 1 else '-'
        return sign + ''.join(partials) + '.'.join(comps)

    def __tex__(self):
        partials = ''.join(
            '\\partial{}'.format(p.index) for p in reversed(self.partials)
        )
        comps = [str(c) for c in self.components]
        comps = '.'.join(c[1:] if c[0] == '-' else c for c in comps)
        sign = '' if self.sign == 1 else '-'
        return sign + partials + comps


class Pair:
    '''A Pair may be any object along with an Alpha value'''
    def __init__(self, a, x=None):
        if x is None:
            if isinstance(a, str) and a.startswith('-'):
                x = Xi(a[1:])
            else:
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
