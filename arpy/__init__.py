# arpy (Absolute Relativity in Python)
# Copyright (C) 2016-2017 Innes D. Anderson-Morrison All rights reserved.

__version__ = '0.1.7'

import types
from copy import deepcopy
from .algebra.config import config, ARConfig
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
# Horrible hacks to get arround cyclic imports #
################################################
def invert_multivector(self):
    '''~mvec as a shortcut for the Hermitian conjugate'''
    inverted = deepcopy(self)
    for alpha, xis in inverted.components.items():
        if full(alpha, alpha).sign == -1:
            for xi in xis:
                xi.sign *= -1
    return inverted


MultiVector.__invert__ = invert_multivector


def update_env(self):
    '''Update the list of predefined operators and multivectors'''
    # Multi-vectors to work with based on the 3-vectors
    self.p = MultiVector('p')
    self.h = MultiVector(self._h)
    self.q = MultiVector(self._q)
    self.t = MultiVector('0')

    self.A = MultiVector(self._A)
    self.B = MultiVector(self._B)
    self.E = MultiVector(self._E)
    self.F = self.E + self.B
    self.T = MultiVector(self._T)
    self.G = MultiVector(self.allowed)

    self.B4 = MultiVector(['p'] + self._B)
    self.T4 = MultiVector(['0'] + self._T)
    self.A4 = MultiVector([self._h] + self._A)
    self.E4 = MultiVector([self._q] + self._E)
    self.Fp = self.F + self.p
    self.F4 = self.F + self.p + self.q

    # Differential operators
    self.Dmu = self.d = differential_operator(['0', '1', '2', '3'])
    self.DG = differential_operator(self.allowed)
    self.DF = differential_operator(self.F)

    self.DB = differential_operator(self.B4)
    self.DT = differential_operator(self.T4)
    self.DA = differential_operator(self.A4)
    self.DE = differential_operator(self.E4)


# Add the update_env method to ARConfig _and_ the config instance
ARConfig.update_env = update_env
config.update_env = types.MethodType(update_env, config)

##############################################################################

# Bring the config definitions into scope
config.update_config()
config.update_env()

p = config.p
h = config.h
q = config.q
t = config.t

A = config.A
B = config.B
E = config.E
F = config.F
T = config.T
G = config.G

B4 = config.B4
T4 = config.T4
A4 = config.A4
E4 = config.E4
Fp = config.Fp
F4 = config.F4

Dmu = d = config.Dmu
DG = config.DG
DF = config.DF
DB = config.DB
DT = config.DT
DA = config.DA
DE = config.DE

# Build the default context for computation
# NOTE:: The user can create a new context in the same way or modify the
#        properties of the original context using .metric and .division
ar = ARContext(config)


def arpy_info():
    '''Display some information about arpy'''
    print('\nNow running arpy version:\t', __version__)
    print('=======================================')
    print('Allowed αs:\t', ', '.join([str(Alpha(a)) for a in config.allowed]))
    print('Division:\t', config.division_type)
    metric = ['+' if i == 1 else '-' for i in config.metric]
    print('Metric:\t\t', ''.join(metric))


# All values that will be imported when the user does `from arpy import *`
__all__ = [
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
    'ar', 'Tex', 'arpy_info', 'config', 'ARConfig'
]
