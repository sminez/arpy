# arpy (Absolute Relativity in Python)
# Copyright (C) 2016-2017 Innes D. Anderson-Morrison All rights reserved.

__version__ = '0.1.2'

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
# Multi-vectors to work with based on the 4Set components #
###########################################################
B = MultiVector([Pair('p')] + [Pair(a) for a in XI_GROUPS['jk']])
T = MultiVector([Pair('0')] + [Pair(a) for a in XI_GROUPS['0jk']])
A = MultiVector([Pair('123')] + [Pair(a) for a in XI_GROUPS['i']])
E = MultiVector([Pair('0123')] + [Pair(a) for a in XI_GROUPS['i0']])
G = MultiVector([Pair(a) for a in ALLOWED])
F = B + E


##############################################################################
# Sepcific Differnetial operators #
###################################
Dmu = Dμ = differential_operator(['0', '1', '2', '3'])
DG = differential_operator(ALLOWED)
DF = differential_operator(XI_GROUPS['jk'] + XI_GROUPS['i0'])

DB = differential_operator(XI_GROUPS['jk'] + ['p'])
DT = differential_operator(XI_GROUPS['0jk'] + ['0'])
DA = differential_operator(XI_GROUPS['i'] + ['123'])
DE = differential_operator(XI_GROUPS['i0'] + ['0123'])

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
    'Dmu', 'Dμ', 'DG', 'DF', 'DB', 'DT', 'DA', 'DE',
    # Differential operator helpers
    'AR_differential', 'differential_operator', 'del_grouped',
    # Visulaisation functions
    'cayley', 'sign_cayley', 'sign_distribution',
    # Pre-defined MultiVectors
    'XiG', 'XiB', 'XiT', 'XiA', 'XiE',
    # The a pre-defined ar() context function
    'ar'
]
