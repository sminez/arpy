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
        ' '.join([str(op(Alpha(a), Alpha(b))).rjust(padding) for b in ALLOWED])
        for a in ALLOWED
    )
    for comp in comps:
        print(comp)


def sign_cayley(op=full):
    '''
    Print +-1 signs for the current Cayley table to the terminal allowing
    for specification of the operation used to compute the table.

    Any function that accepts two Alphas can be passed as op.
    '''
    divider = '      ' + ''.join('+---------' for _ in range(4)) + '+'
    comps = (
        ' '.join([
            '■' if op(Alpha(a), Alpha(b)).sign == -1 else '□'
            for b in ALLOWED
        ])
        for a in ALLOWED
    )

    print('          ', '         '.join(['B', 'A', 'T', 'E']))
    print(divider)

    for i, comp in enumerate(comps):
        comp = '| '.join(comp[n:n+8] for n in range(0, len(comp), 8))
        print(str(Alpha(ALLOWED[i])).ljust(5), '|', comp, '|')
        # Divide after each 4-Set
        if (i + 1) % 4 == 0:
            print(divider)
