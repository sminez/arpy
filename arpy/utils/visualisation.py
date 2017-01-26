'''
Tools for visualising results from calculations and investigations.
http://seaborn.pydata.org/generated/seaborn.heatmap.html
'''
from ..config import ALLOWED
from .algebra import a


def print_cayley(op="full"):
    '''
    Dump the current Cayley table to the terminal allowing for specification
    of the operation used to compute the table.
    '''
    if op == "full":
        cayley = [[a(i) * a(j) for j in ALLOWED] for i in ALLOWED]
    elif op == "into":
        cayley = [[a(i).inverse() * a(j) for j in ALLOWED] for i in ALLOWED]
    elif op == "by":
        cayley = [[a(i) * a(j).inverse() for j in ALLOWED] for i in ALLOWED]

    for k in cayley:
        line = [str(a).rjust(6) for a in k]
        print(" ".join(line))
