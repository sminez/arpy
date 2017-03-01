'''
arpy (Absolute Relativity in Python)
Copyright (C) 2016-2017 Innes D. Anderson-Morrison All rights reserved.
'''
# The labelling and ordering of the 16 elements of the algebra.
# NOTE:: The order will affect the visualisation of the Cayley Table
#       but not the results of finding products.
_B = ['p', '23', '31', '12']     # ΞB :: Magnetic Field and rest mass
_T = ['0', '023', '031', '012']  # ΞΤ :: Angular Momentum Density and time
_A = ['123', '1', '2', '3']      # ΞΑ :: Charge Current Density and hedgehog
_E = ['0123', '10', '20', '30']  # ΞE :: Electric Field and dual rest mass
ALLOWED = _B + _T + _A + _E


SUPER_SCRIPTS = {'B': 'ᴮ', 'A': 'ᴬ', 'T': 'ᵀ', 'E': 'ᴱ'}
SUB_SCRIPTS = {
    '0': '₀', '1': '₁', '2': '₂', '3': '₃',
    'p': 'ₚ', '123': '₁₂₃', '0123': '₀₁₂₃'
}
GROUP_TO_4SET = {'jk': 'B', 'i': 'A', '0jk': 'T', 'i0': 'E'}


# Map α to 4set membership
FOUR_SETS = {comp: 'B' for comp in _B}
FOUR_SETS.update({comp: 'T' for comp in _T})
FOUR_SETS.update({comp: 'A' for comp in _A})
FOUR_SETS.update({comp: 'E' for comp in _E})

# Fast lookup of 4set components in {t,x,y,z} order
_dims = 'b x y z'.split()
FOUR_SET_COMPS = {
    'B': dict(zip(_dims, _B)),
    'T': dict(zip(_dims, _T)),
    'A': dict(zip(_dims, _A)),
    'E': dict(zip(_dims, _E))
}

bxyz_pairings = [(s[0], 'b') for s in [_B, _T, _A, _E]]
bxyz_pairings.extend([(s[1], 'x') for s in [_B, _T, _A, _E]])
bxyz_pairings.extend([(s[2], 'y') for s in [_B, _T, _A, _E]])
bxyz_pairings.extend([(s[3], 'z') for s in [_B, _T, _A, _E]])
BXYZ_LIKE = dict(bxyz_pairings)


# How the 3-vector components are grouped and under what names
XI_GROUPS = {
    'i': ['1', '2', '3'],
    'i0': [a for a in ALLOWED if len(a) == 2 and '0' in a],
    'jk': [a for a in ALLOWED if len(a) == 2 and '0' not in a],
    '0jk': [a for a in ALLOWED if len(a) == 3 and '0' in a]
}

# Names to group the results of calculations under: scalars & 3-vectors
ALLOWED_GROUPS = ['p', '0', '123', '0123'] + [g for g in XI_GROUPS.keys()]

# The space-time metric that will be used
METRIC = (1, -1, -1, -1)

# Whether division is defined as 'by' or 'into'
# NOTE:: Any other values will raise exceptions within the rest of the code!
DIVISION_TYPE = 'into'


##############################################################################
def _build_alpha_to_group():
    # Scalars are grouped individually
    _pairs = [('p', 'p'), ('0', '0'), ('123', '123'), ('0123', '0123')]
    # 3-vector components are grouped under the vector name
    flipped = [[(v, group) for v in vals] for group, vals in XI_GROUPS.items()]
    for group in flipped:
        _pairs.extend(group)
    return dict(_pairs)
##############################################################################

# For a given alpha, find the group it should be assigned to
ALPHA_TO_GROUP = _build_alpha_to_group()
