'''
Tools for visualising results from calculations and investigations.
'''
from ..algebra.config import ALLOWED
from ..algebra.ar_types import Alpha
from ..algebra.operations import full


def cayley(op=full, padding=6):
    '''
    Print current Cayley table to the terminal allowing for specification
    of the operation used to compute the table.

    Any function that accepts two Alphas can be passed as op.
    '''
    comps = (
        ' '.join([str(op(Alpha(a), Alpha(b))).rjust(padding) for a in ALLOWED])
        for b in ALLOWED
    )
    for comp in comps:
        print(comp)
