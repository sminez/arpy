from itertools import combinations
from arpy import Alpha, ALLOWED


# Data used in multiple tests
ap = Alpha('p')                         # αp
neg_ap = Alpha('p', -1)                 # -αp
txyz = ['0', '1', '2', '3']             # space-time indices
metrics = [                             # All possible +- metrics
    (t, x, y, z)                        # NOTE:: This is not just +---/-+++
    for t in [1, -1] for x in [1, -1]   # > it covers any set of 4
    for y in [1, -1] for z in [1, -1]
]

# Pre-zipped combinations of indices for paramatizing tests
imetric = [i_metric for i_metric in zip(ALLOWED, metrics)]
ijmetric = [
    (c[0], c[1], m)
    for c in combinations(txyz, 2) for m in metrics
]
ijkmetric = [
    (c[0], c[1], c[2], m)
    for c in combinations(txyz, 3) for m in metrics
]
