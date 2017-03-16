# arpy (Absolute Relativity in Python)
# Copyright (C) 2016-2017 Innes D. Anderson-Morrison All rights reserved.

__version__ = '0.1.4'

from copy import deepcopy
from .algebra.config import ALLOWED, XI_GROUPS, METRIC, DIVISION_TYPE, \
        ALPHA_TO_GROUP, ALLOWED_GROUPS, FOUR_SET_COMPS, FOUR_SETS, \
        BXYZ_LIKE
from .algebra.ar_types import Alpha, Xi, Pair
from .algebra.multivector import MultiVector, DelMultiVector
from .algebra.operations import find_prod, inverse, full, div_by, div_into, \
        project, prod_apply, dagger, commutator
from .algebra.differential import AR_differential, differential_operator
from .algebra.del_grouping import del_grouped
from .utils.lexparse import ARContext
from .utils.visualisation import cayley, sign_cayley, sign_distribution


##############################################################################
# Horrible hack to get arround cyclic imports #
###############################################
def invert_multivector(self):
    # ~mvec as a shortcut for the Hermitian conjugate
    inverted = deepcopy(self)
    for alpha, xis in inverted.components.items():
        if full(alpha, alpha).sign == -1:
            for xi in xis:
                xi.sign *= -1
    return inverted

MultiVector.__invert__ = invert_multivector


##############################################################################
# Multi-vectors to work with based on the 3/4-vectors #
#######################################################
A = MultiVector('0 1 2 3')                            # The potentials
B = MultiVector(XI_GROUPS['jk'])                      # The Magnetic field
E = MultiVector(XI_GROUPS['i0'])                      # The Electric field
F = E + (-B)                                          # The Farady tensor
T = MultiVector([a for a in ALLOWED if len(a) == 3])  # The trivectors
G = MultiVector(ALLOWED)                              # The general multivector
##############################################################################
# Multi-vectors to work with based on the 4Set components #
###########################################################
B4 = MultiVector([Pair('p')] + [Pair(a) for a in XI_GROUPS['jk']])
T4 = MultiVector([Pair('0')] + [Pair(a) for a in XI_GROUPS['0jk']])
A4 = MultiVector([Pair('123')] + [Pair(a) for a in XI_GROUPS['i']])
E4 = MultiVector([Pair('0123')] + [Pair(a) for a in XI_GROUPS['i0']])
Fp = F + MultiVector('p')
F4 = Fp + MultiVector('0123')


##############################################################################
# Sepcific Differnetial operators #
###################################
Dmu = d = differential_operator(['0', '1', '2', '3'])
DG = differential_operator(ALLOWED)
DF = differential_operator(F)

DB = differential_operator(B4)
DT = differential_operator(T4)
DA = differential_operator(A4)
DE = differential_operator(E4)

# Build the default context for computation
# NOTE:: The user can create a new context in the same way or modify the
#        properties of the original context using .metric and .division
ar = ARContext(METRIC, DIVISION_TYPE)


# All values that will be imported when the user does `from arpy import *`
__all__ = [
    # Config initialised values
    'ALLOWED', 'XI_GROUPS', 'METRIC', 'DIVISION_TYPE', 'ALPHA_TO_GROUP',
    'ALLOWED_GROUPS', 'FOUR_SET_COMPS', 'FOUR_SETS', 'BXYZ_LIKE',
    # Data structures
    'Alpha', 'Xi', 'Pair', 'MultiVector', 'DelMultiVector',
    # Non differential operators
    'find_prod', 'inverse', 'full', 'div_by', 'div_into',
    'project', 'prod_apply', 'dagger', 'commutator',
    # Differential operators
    'Dmu', 'd', 'DG', 'DF', 'DB', 'DT', 'DA', 'DE',
    # Differential operator helpers
    'AR_differential', 'differential_operator', 'del_grouped',
    # Visulaisation functions
    'cayley', 'sign_cayley', 'sign_distribution',
    # Pre-defined MultiVectors
    'G', 'F', 'Fp', 'B', 'T', 'A', 'E',
    'B4', 'T4', 'A4', 'E4', 'F4',
    # The a pre-defined ar() context function
    'ar'
]
