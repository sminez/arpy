import collections.abc
from itertools import groupby
from .ar_types import Alpha, Pair
from .del_grouping import del_grouped
from .config import ALLOWED, ALLOWED_GROUPS, ALPHA_TO_GROUP


class MultiVector(collections.abc.Set):
    '''A custom container type for working efficiently with multivectors'''
    _allowed_alphas = ALLOWED

    def __init__(self, components=[]):
        # Given a list of pairs, build the mulitvector by binding the ξ values
        self.components = {Alpha(a): [] for a in self._allowed_alphas}

        for comp in components:
            if isinstance(comp, (str, Alpha)):
                comp = Pair(comp)
            if not isinstance(comp, Pair):
                raise ValueError('Arguments must be Alphas, Pairs or Strings')
            if comp.alpha.index in self._allowed_alphas:
                try:
                    self.components[comp.alpha].append(comp.xi)
                except KeyError:
                    # Negative Alpha value
                    alpha, xi = comp.alpha, comp.xi
                    alpha.sign = 1
                    xi.sign *= -1
                    self.components[alpha].append(xi)

    def __repr__(self):
        comps = ['  α{}{}'.format(str(a).ljust(5), self._nice_xi(Alpha(a)))
                 for a in self._allowed_alphas if self.components[Alpha(a)]]
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
        for alpha in self._allowed_alphas:
            try:
                for xi in self.components[Alpha(alpha)]:
                    yield Pair(alpha, xi)
            except KeyError:
                pass

    def _nice_xi(self, alpha, raise_key_error=False, for_print=True):
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

    def simplified(self, on='xi', ix=0):
        '''
        Display the multivector with simplified Xi values.
        '''
        ix2 = 1 if ix == 0 else 1
        if on == 'xi':
            def key(x):
                return (x.xi.components[ix], ALPHA_TO_GROUP[x.alpha.index])
        elif on == 'alpha':
            def key(x):
                return (ALPHA_TO_GROUP[x.alpha.index], x.xi.components[ix])
        else:
            raise ValueError('"on" must be one of "xi" or "alpha"')

        g = groupby(sorted(self, key=key), key)
        seen = []

        print('{')
        for key, full in g:
            comps = [Pair(f.alpha, f.xi.components[ix2]) for f in full]
            comps = sorted(del_grouped(comps), key=lambda x: x.alpha.index)
            agrouped = groupby(comps, lambda x: x.alpha.index)
            if key[0] not in seen:
                print('  {}'.format(key[0]))
                seen.append(key[0])
            for alpha, components in agrouped:
                components = [c for c in components]
                for component in components:
                    if component.alpha.sign == -1:
                        component.alpha.sign = 1
                        component.xi.sign *= -1
                comp_str = ' '.join([str(c.xi) for c in components])
                print('    {}: {}'.format(Alpha(alpha), comp_str))
        print('}')


class DelMultiVector(MultiVector):
    _allowed_alphas = ALLOWED_GROUPS  # + ALLOWED
    # _allowed_alphas = '0 123 i 0jk p 0123 i0 jk'.split()
