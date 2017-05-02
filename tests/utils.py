from itertools import combinations
from arpy import Alpha


# Data used in multiple tests
metrics = [                             # All possible +- metrics
    (t, x, y, z)                        # NOTE:: This is not just +---/-+++
    for t in [1, -1] for x in [1, -1]   # > it covers any set of 4
    for y in [1, -1] for z in [1, -1]
]

# Pre-zipped combinations of indices for paramatizing tests
txyz = ['0', '1', '2', '3']
ij_pairs = [(c[0], c[1]) for c in combinations(txyz, 2)]
ijk_triplets = [(c[0], c[1], c[2]) for c in combinations(txyz, 3)]
