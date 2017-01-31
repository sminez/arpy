from .algebra.config import *
from .algebra.ar_types import *
from .algebra.operations import *
from .algebra.differential import *
from .utils.lexparse import ar


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
