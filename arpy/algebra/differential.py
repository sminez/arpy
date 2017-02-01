'''
A selection of different implementations for symbolically computing
the 4-vector 4-differential Dμ.

In Cartesian coordinates this is:
    Dμ = ∂μ/αμ = ∂ / (αμ∂xμ) = α0∂0 - αi∂i = α0∂0 - ∇i
'''
from itertools import groupby

from .config import ALLOWED, ALPHA_TO_GROUP, DIVISION_TYPE
from .ar_types import Alpha, Xi, Pair, MultiVector
from .operations import div_by, div_into


def component_partial(component, wrt, div):
    '''
    Symbolically differentiate a component by storing the partials and
    converting the alpha value using the correct division type.
    '''
    if div == 'by':
        alpha = div_by(component.alpha, wrt)
    elif div == 'into':
        alpha = div_into(component.alpha, wrt)
    else:
        raise ValueError('Invalid division specification: %s' % div)
    partials = component.xi.partials + [wrt]
    return Pair(alpha, Xi(component.xi.val, component.xi.unit, partials))


def AR_differential(mvec, wrt, div=DIVISION_TYPE):
    '''
    Compute the result of Differentiating a each component of a MultiVector
    with respect to a given list of unit elements under the algebra.
    '''
    return MultiVector([
        component_partial(component, Alpha(element), div)
        for element in wrt for component in mvec
    ])


def Dmu(mvec):
    '''The main operator from the paper'''
    return AR_differential(mvec, ['0', '1', '2', '3'])


def DG(mvec):
    '''A full derivative with respect to all components'''
    return AR_differential(mvec, ALLOWED)

##############################################################################


# TODO:: Move these to being methods on MultiVectors
def by_alpha(vec, vector_groups=False):
    '''
    Group results of a derivative operation together in order to form the
    resulting coupled differential equations.
    '''
    vec = sorted(vec, key=lambda a: ALLOWED.index(a.alpha.index))

    if vector_groups:
        return groupby(vec, lambda x: ALPHA_TO_GROUP[x.alpha.index])
    else:
        return groupby(vec, key=lambda x: x.alpha.index)


def show(vec, vector_groups=False):
    '''Display a result and optionally group to scalars and 3-vectors'''
    for alpha, group in by_alpha(vec, vector_groups):
        xis = [('-' if g.alpha.sign == -1 else '+') + str(g.xi) for g in group]
        print('α{} ['.format(alpha), ' '.join(xis), ']')
