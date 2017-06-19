'''
arpy (Absolute Relativity in Python)
Copyright (C) 2016-2017 Innes D. Anderson-Morrison All rights reserved.
'''


class ARConfig:
    '''The arpy paramater configuration object'''
    def __init__(self, allowed, metric, div):
        '''Bind in all parameters'''
        self.allowed = allowed
        self.metric = metric
        self.division_type = div

        # Generate the config
        self.update_config()

    def update_config(self):
        '''
        Define algebra level data and mappings.
        NOTE: The ARConfig class is extended to include an `update_env` method
        in the main __init__.py.
        '''
        self._h = [a for a in self.allowed if len(a) == 3 and '0' not in a][0]
        self._q = [a for a in self.allowed if len(a) == 4][0]
        self._B = [a for a in self.allowed if len(a) == 2 and '0' not in a]
        self._T = [a for a in self.allowed if len(a) == 3 and '0' in a]
        self._A = [a for a in self.allowed if len(a) == 1 and a not in 'p0']
        self._E = [a for a in self.allowed if len(a) == 2 and '0' in a]

        # Map α to 4set membership
        self.four_sets = {comp: 'B' for comp in ['p'] + self._B}
        self.four_sets.update({comp: 'T' for comp in ['0'] + self._T})
        self.four_sets.update({comp: 'A' for comp in [self._h] + self._A})
        self.four_sets.update({comp: 'E' for comp in [self._q] + self._E})

        # Fast lookup of 4set components in {t,x,y,z} order
        _dims = 'b x y z'.split()
        self.four_set_comps = {
            'B': dict(zip(_dims, ['p'] + self._B)),
            'T': dict(zip(_dims, ['0'] + self._T)),
            'A': dict(zip(_dims, [self._h] + self._A)),
            'E': dict(zip(_dims, [self._q] + self._E))
        }

        _groups = [['p'] + self._B, ['0'] + self._T,
                   [self._h] + self._A, [self._q] + self._E]
        bxyz_pairings = [(s[0], 'b') for s in _groups]
        bxyz_pairings.extend([(s[1], 'x') for s in _groups])
        bxyz_pairings.extend([(s[2], 'y') for s in _groups])
        bxyz_pairings.extend([(s[3], 'z') for s in _groups])

        self.bxyz_like = dict(bxyz_pairings)

        # How the 3-vector components are grouped and under what names
        # TODO: Have a way to dynamically alter these names?
        self.xi_groups = {
            'i': ['1', '2', '3'],
            'i0': [a for a in self.allowed if len(a) == 2 and '0' in a],
            'jk': [a for a in self.allowed if len(a) == 2 and '0' not in a],
            '0jk': [a for a in self.allowed if len(a) == 3 and '0' in a]
        }

        # Names to group the results of calculations under: scalars & 3-vectors
        self.allowed_groups = ['p', '0', '123', '0123'] + \
            [g for g in self.xi_groups.keys()]

        # For a given alpha, find the group it should be assigned to
        self.alpha_to_group = self._build_alpha_to_group()

    def _build_alpha_to_group(self):
        # Scalars are grouped individually
        _pairs = [
            ('p', 'p'), ('0', '0'), (self._h, self._h), (self._q, self._q)]
        # 3-vector components are grouped under the vector name
        flipped = [[(v, group) for v in vals]
                   for group, vals in self.xi_groups.items()]
        for group in flipped:
            _pairs.extend(group)
        return dict(_pairs)


# The labelling and ordering of the 16 elements of the algebra.
# NOTE:: The order will affect the visualisation of the Cayley Table
#       but not the results of finding products.
_B = ['p', '23', '31', '12']     # ΞB :: Magnetic Field and rest mass
_T = ['0', '023', '031', '012']  # ΞΤ :: Angular-Momentum/Charge density
_A = ['123', '1', '2', '3']      # ΞΑ :: Current Density and hedgehog
_E = ['0123', '01', '02', '03']  # ΞE :: Electric Field and dual rest mass
allowed = _B + _T + _A + _E

# The space-time metric that will be used
metric = (1, -1, -1, -1)

config = ARConfig(allowed, metric, 'into')
