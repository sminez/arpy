# arpy (Absolute Relativity in Python)
# Copyright (C) 2016-2018 Innes D. Anderson-Morrison All rights reserved.

__version__ = "0.3.5"

import types
from copy import copy
from ctypes import c_int, py_object, pythonapi
from sys import _getframe

from .algebra.data_types import Alpha, MultiVector, Term, Xi
from .algebra.differential import AR_differential
from .algebra.operations import (
    MM_bar,
    commutator,
    dagger,
    diamond,
    div_by,
    div_into,
    dual,
    find_prod,
    full,
    hermitian,
    inverse,
    project,
    rev,
)
from .config import ARConfig, config
from .reductions.reducers import del_grouped, replace_all
from .utils.lexparse import ARContext
from .utils.utils import Nat, Tex, Zet, reorder_allowed
from .utils.visualisation import cayley, js_cayley, op_block, sign_cayley, sign_distribution


##############################################################################
# Horrible hacks to get arround cyclic imports #
################################################
def invert_multivector(self):
    """~mvec as a shortcut for the Hermitian conjugate"""
    inverted = copy(self)
    for alpha, xis in inverted.components.items():
        if full(alpha, alpha).sign == -1:
            for xi in xis:
                xi.sign *= -1
    return inverted


MultiVector.__invert__ = invert_multivector


def update_env(self, lvl=2):
    """Update the list of predefined operators and multivectors"""

    def _bind_to_calling_scope(defs, lvl):
        """
        Inject the default Multivectors and operators into the main scope
        of the repl. (THIS IS HORRIFYING!!!)
        NOTE: This uses some not-so-nice abuse of stack frames and the
              ctypes API to make this work and as such it will almost
              certainly not run under anything other than cPython.
        """
        # Grab the stack frame that the caller's code is running in
        frame = _getframe(lvl)
        # Dump the matched variables and their values into the frame
        frame.f_locals.update(defs)
        # Force an update of the frame locals from the locals dict
        pythonapi.PyFrame_LocalsToFast(py_object(frame), c_int(0))

    # Multi-vectors to work with based on the 3-vectors
    self.p = MultiVector("p", cfg=self)
    self.h = MultiVector(self._h, cfg=self)
    self.q = MultiVector(self._q, cfg=self)
    self.t = MultiVector("0", cfg=self)

    self.A = MultiVector(self._A, cfg=self)
    self.B = MultiVector(self._B, cfg=self)
    self.E = MultiVector(self._E, cfg=self)
    self.F = self.E + self.B
    self.T = MultiVector(self._T, cfg=self)
    self.G = MultiVector(self.allowed, cfg=self)

    self.zet_B = MultiVector(["p"] + self._B, cfg=self)
    self.zet_T = MultiVector(["0"] + self._T, cfg=self)
    self.zet_A = MultiVector([self._h] + self._A, cfg=self)
    self.zet_E = MultiVector([self._q] + self._E, cfg=self)
    self.Fp = MultiVector(["p"] + self._B + self._E, cfg=self)
    self.zet_F = MultiVector(["p", self._q] + self._B + self._E, cfg=self)
    self.Fpq = self.zet_F

    # Differential operators
    self.Dmu = self.d = AR_differential(["0", "1", "2", "3"], cfg=self)
    self.DG = AR_differential(self.allowed, cfg=self)
    self.DF = AR_differential(self.F, cfg=self)

    self.DB = AR_differential(self.zet_B, cfg=self)
    self.DT = AR_differential(self.zet_T, cfg=self)
    self.DA = AR_differential(self.zet_A, cfg=self)
    self.DE = AR_differential(self.zet_E, cfg=self)

    _vars = [
        "p",
        "h",
        "q",
        "t",
        "A",
        "B",
        "E",
        "F",
        "T",
        "G",
        "zet_B",
        "zet_T",
        "zet_A",
        "zet_E",
        "Fp",
        "zet_F",
        "Dmu",
        "d",
        "DG",
        "DF",
        "DB",
        "DT",
        "DA",
        "DE",
    ]
    defs = dict(zip(_vars, (getattr(self, var) for var in _vars)))
    _bind_to_calling_scope(defs, lvl)


# Add the update_env method to ARConfig _and_ the config instance
ARConfig.update_env = update_env
config.update_env = types.MethodType(update_env, config)

##############################################################################

# Bring the config definitions into scope
config.update_config()
config.update_env()

# Build the default context for computation
# NOTE:: The user can create a new context in the same way or modify the
#        properties of the original context using .metric and .division
ar = ARContext(cfg=config)


def arpy_info():
    """Display some information about arpy"""
    print("\nNow running arpy version:\t", __version__)
    print("=======================================")
    print("Allowed αs:\t", ", ".join([str(Alpha(a)) for a in config.allowed]))
    print("Division:\t", config.division_type)
    metric = ["+" if i == 1 else "-" for i in config.metric]
    print("Metric:\t\t", "".join(metric))


# All values that will be imported when the user does `from arpy import *`
__all__ = [
    # Data structures
    "Alpha",
    "Xi",
    "Term",
    "MultiVector",
    # Non differential operators
    "find_prod",
    "inverse",
    "full",
    "div_by",
    "div_into",
    "project",
    "dagger",
    "hermitian",
    "commutator",
    "diamond",
    "rev",
    "dual",
    "MM_bar",
    # Differential operators
    "Dmu",
    "d",
    "DG",
    "DF",
    "DB",
    "DT",
    "DA",
    "DE",
    # Differential operator helpers
    "AR_differential",
    "del_grouped",
    # Visulaisation functions
    "cayley",
    "sign_cayley",
    "sign_distribution",
    "js_cayley",
    "op_block",
    # Pre-defined MultiVectors
    "G",
    "F",
    "Fp",
    "B",
    "T",
    "A",
    "E",
    "zet_B",
    "zet_T",
    "zet_A",
    "zet_E",
    "zet_F",
    # Util functions
    "ar",
    "Tex",
    "arpy_info",
    "config",
    "ARConfig",
    "ARContext",
    "reorder_allowed",
    "Zet",
    "Nat",
    "replace_all",
]
