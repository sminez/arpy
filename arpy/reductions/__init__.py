"""
A utility module for performing term level simplification and notation
changes on MultiVectors.
"""
from .reducers import SUB_SCRIPTS, SUPER_SCRIPTS, Template, Term, Replacement, chain_reducers


__all__ = ["SUB_SCRIPTS", "SUPER_SCRIPTS", "Template", "Term", "Replacement", "chain_reducers"]
