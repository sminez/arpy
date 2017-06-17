# arpy (Absolute Relativity in Python)
# Copyright (C) 2016-2017 Innes D. Anderson-Morrison All rights reserved.

__version__ = '0.1.7'

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
from .utils.utils import Tex
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

# Indices for alphas
_h = [a for a in ALLOWED if len(a) == 3 and '0' not in a][0]
_q = [a for a in ALLOWED if len(a) == 4][0]
_B = [a for a in ALLOWED if len(a) == 2 and '0' not in a]
_T = [a for a in ALLOWED if len(a) == 3 and '0' in a]
_A = ['0', '1', '2', '3']
_E = [a for a in ALLOWED if len(a) == 2 and '0' in a]

##############################################################################
# Multi-vectors to work with based on the 3/4-vectors #
#######################################################

p = MultiVector('p')
h = MultiVector(_h)
q = MultiVector(_q)
t = MultiVector('0')

A = MultiVector(_A)
B = MultiVector(_B)
E = MultiVector(_E)
F = E + B
T = MultiVector(_T)
G = MultiVector(ALLOWED)

##############################################################################
# Multi-vectors to work with based on the 4Set components #
###########################################################
B4 = MultiVector(['p'] + _B)
T4 = MultiVector(['0'] + _T)
A4 = MultiVector([_h] + _A)
E4 = MultiVector([_q] + _E)
Fp = F + p
F4 = F + p + q


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
ar = ARContext(METRIC, ALLOWED)


def arpy_info():
    '''Display some information about arpy'''
    print('\nNow running arpy version:\t', __version__)
    print('=======================================')
    print('Allowed αs:\t', ', '.join([str(Alpha(a)) for a in ALLOWED]))
    print('Division:\t', DIVISION_TYPE)
    print('Metric:\t\t', ''.join(['+' if i == 1 else '-' for i in METRIC]))


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
    # The a pre-defined ar() context function and Tex output
    'ar', 'Tex', 'arpy_info'
]
