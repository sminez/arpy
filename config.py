# The labelling and ordering of the 16 elements of the algebra.
# NOTE:: The order will affect the visualisation of the Cayley Table
#        but not the results of finding products.
ALLOWED = [
    'p', '23', '31', '12',          # ΞM :: Magnetic and rest mass
    '0', '023', '031', '012',       # ΞE :: Electric and dual rest mass
    '123', '1', '2', '3',           # ΞΑ :: Vector Potential and dual time(?)
    '0123', '10', '20', '30'        # ΞΤ :: Spin and time
]

# How the 3-vector components are grouped and under what names
ξ_GROUPS = [
    ('i', ['1', '2', '3']),
    ('i0', [a for a in ALLOWED if len(a) == 2 and '0' in a]),
    ('jk', [a for a in ALLOWED if len(a) == 2 and '0' not in a]),
    ('0jk', [a for a in ALLOWED if len(a) == 3 and '0' in a])
]

# Names to group the results of calculations under: scalars & 3-vectors
ALLOWED_GROUPS = ['p', '0', '123', '0123'] + [g[0] for g in ξ_GROUPS]

# The space-time metric that will be used
METRIC = [1, -1, -1, -1]

# Whether division is defined as 'by' or 'into'
# NOTE:: Any other values will raise exceptions within the rest of the code!
DIVISION_TYPE = 'by'


##############################################################################
def build_alpha_to_group():
    # Scalars are grouped individually
    _pairs = [('p', 'p'), ('0', '0'), ('123', '123'), ('0123', '0123')]
    # 3-vector components are grouped under the vector name
    for g in [[(μ, group) for μ in αs] for (group, αs) in ξ_GROUPS]:
        _pairs += g
    return dict(_pairs)
##############################################################################

# For a given alpha, find the group it should be assigned to
α_TO_GROUP = build_alpha_to_group()
