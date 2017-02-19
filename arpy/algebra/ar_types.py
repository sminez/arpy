'''
Classes, functions and default representations of Ξ vectors in the algebra.
'''
import collections.abc
from itertools import groupby
from .config import ALLOWED, ALLOWED_GROUPS, ALPHA_TO_GROUP


CW = {1: 2, 2: 3, 3: 1}
ACW = {1: 3, 2: 1, 3: 2}


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

    def __lt__(self, other):
        return ALLOWED.index(self.index) < ALLOWED.index(other.index)

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
        partials = ('∂{}'.format(p.index) for p in reversed(self.partials))
        return '{}{}ξ{}'.format(sign, ''.join(partials), self.val)


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
        partials = ('∂{}'.format(p.index) for p in reversed(self.partials))
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


class MultiVector(collections.abc.Set):
    '''A custom container type for working efficiently with multivectors'''
    def __init__(self, components=[]):
        # Given a list of pairs, build the mulitvector by binding the ξ values
        self.components = {Alpha(a): [] for a in ALLOWED}

        for comp in components:
            if isinstance(comp, (str, Alpha)):
                comp = Pair(comp)
            if not isinstance(comp, Pair):
                raise ValueError('Arguments must be Alphas, Pairs or Strings')
            try:
                self.components[comp.alpha].append(comp.xi)
            except KeyError:
                # Negative Alpha value
                alpha, xi = comp.alpha, comp.xi
                alpha.sign = 1
                xi.sign *= -1
                self.components[alpha].append(xi)

    def __repr__(self):
        comps = ['  α{} {}'.format(a, self._nice_xi(Alpha(a), False, True))
                 for a in ALLOWED if self.components[Alpha(a)]]
        return '{\n' + '\n'.join(comps) + '\n}'

    def __len__(self):
        # All allowed values are initialised with [] so we are
        # only counting componets that have a Xi value set.
        return len([v for v in self.components.values() if v != []])

    def __add__(self, other):
        # Allow for the addition of Multivectors and Pairs.
        # This will always return a new MultiVector.
        if not isinstance(other, (Pair, MultiVector)):
            raise TypeError()

        comps = [p for p in self]
        if isinstance(other, Pair):
            comps.append(other)
        elif isinstance(other, MultiVector):
            comps.extend(p for p in other)

        return MultiVector(comps)

    def __contains__(self, other):
        if isinstance(other, Alpha):
            return self.components[other] != []
        elif isinstance(other, Pair):
            return other.xi in self.components[other.alpha]
        else:
            return False

    def __getitem__(self, key):
        if isinstance(key, str):
            # Allow retreval by string as well as Alpha
            key = Alpha(key)
        if not isinstance(key, Alpha):
            raise KeyError

        xis = self.components[key]
        return [Pair(key, x) for x in xis]

    def __iter__(self):
        for alpha in ALLOWED:
            xi = self._nice_xi(Alpha(alpha), False)
            if xi:
                if isinstance(xi, list):
                    for x in xi:
                        yield Pair(alpha, x)
                else:
                    yield Pair(alpha, xi)

    def _nice_xi(self, alpha, raise_key_error=True, for_print=False):
        '''Single element xi lists return their value raw'''
        try:
            xi = self.components[alpha]
        except KeyError:
            if raise_key_error:
                raise KeyError
        if len(xi) == 1:
            return xi[0]
        else:
            if for_print:
                return '(' + ', '.join(str(x) for x in xi) + ')'
            else:
                return xi

    def MTAE_grouped(self):
        '''
        Print an MTAE grouped representation of the MultiVector
        NOTE:: This deliberately does not return a new MultiVector as we
               should always be working with strict alpha values not grouped.
        '''
        by_alpha = groupby(self, key=lambda x: ALPHA_TO_GROUP[x.alpha.index])
        MTAE = [
            (group, tuple(c.xi for c in components))
            for (group, components) in by_alpha
        ]
        print('{')
        for group, comps in MTAE:
            print('  α{}'.format(group).ljust(7), comps)
        print('}')

    def del_notation(self):
        '''
        Print an MTAE grouped representation of the multivector with del
        vector derivative notation if possible.
        NOTE:: This deliberately does not return a new MultiVector as we
               should always be working with strict alpha values not grouped.
        '''
        print(del_notation(self))

    def collected_terms(self):
        '''Display the multivector with factorised Xi values'''
        print('{')

        sorted(self, key=lambda x: (x.alpha, x.xi.components[0]))
        g = groupby(self, lambda x: (x.alpha, x.xi.components[0]))
        alphas_seen = []

        for key, full in g:
            current_alpha = key[0]
            comps = [f.xi.components[1:] for f in full]
            s = []
            for c in comps:
                s.extend(c)
            if s:
                if current_alpha not in alphas_seen:
                    print('  {}:'.format(key[0]))
                    alphas_seen.append(current_alpha)
                print('    {}({})'.format(key[1], ''.join(str(x) for x in s)))
        print('}')


def to_del(group, components, replacement, sign=None):
    '''Replace a list of components with an alternative Pair'''
    sign = sign if sign else components[0].alpha.sign
    return Pair(group, replacement)


def del_notation(mvec):
    '''
    Return a formatted string representation of a MultiVector that is
    MTAE grouped and expressed in del notation.
    '''
    def curl(MTAE_group):
        return MTAE_group

    def grad(MTAE_group):
        return MTAE_group

    def div(MTAE_group):
        return MTAE_group

    def partials(MTAE_group):
        return MTAE_group

    def terms(MTAE_group):
        '''Ensure that all terms have MTAE alphas'''
        return [comp.xi for _, comp in MTAE_group]

    del_replaced = []
    grouped = groupby(mvec, key=lambda x: ALPHA_TO_GROUP[x.alpha.index])

    for group, components in grouped:
        MTAE_group = [(group, comp) for comp in components]
        del_terms = terms(partials(div(grad(curl(MTAE_group)))))
        del_replaced.append((group, del_terms))

    formatted = [
        '  α{}'.format(group).ljust(7) + '{}'.format(comps)
        for (group, comps) in del_replaced
    ]
    return '{\n' + '\n'.join(formatted) + '\n}'
