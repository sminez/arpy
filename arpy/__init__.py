from .algebra.config import ALLOWED, XI_GROUPS, METRIC, DIVISION_TYPE, ALPHA_TO_GROUP
from .algebra.ar_types import Alpha, Xi, Pair, MultiVector
from .algebra.operations import find_prod, inverse, wedge, dot, full, \
        div_by, div_into, project
from .algebra.differential import AR_differential, Dmu, DG, differential_operator
from .utils.lexparse import ARContext


##############################################################################
# Multi-vectors to work with based on the 4-vector components
XiG = MultiVector([Pair(a) for a in ALLOWED])
XiM = MultiVector([Pair('p')] + [Pair(a) for a in XI_GROUPS['jk']])
XiT = MultiVector([Pair('0')] + [Pair(a) for a in XI_GROUPS['0jk']])
XiA = MultiVector([Pair('123')] + [Pair(a) for a in XI_GROUPS['i']])
XiE = MultiVector([Pair('0123')] + [Pair(a) for a in XI_GROUPS['i0']])
# Multi-vectors based on grade
Xi1 = project(XiG, 1)
Xi2 = project(XiG, 2)
Xi3 = project(XiG, 3)
Xi4 = project(XiG, 4)

# Build the default context for computation
# NOTE:: The user can create a new context in the same way or modify the
#        properties of the original context using .metric and .division
ar = ARContext(METRIC, DIVISION_TYPE)
