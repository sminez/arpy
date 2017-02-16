'''
A selection of different implementations for symbolically computing
the 4-vector 4-differential Dμ.

In Cartesian coordinates this is:
    Dμ = ∂μ/αμ = ∂ / (αμ∂xμ) = α0∂0 - αi∂i = α0∂0 - ∇i
'''
from copy import deepcopy
from .config import ALLOWED, DIVISION_TYPE
from .ar_types import Alpha, MultiVector
from .operations import div_by, div_into


def component_partial(component, wrt, div):
    '''
    Symbolically differentiate a component by storing the partials and
    converting the alpha value using the correct division type.
    '''
    # NOTE:: using deep copy so that all of the objects inside of the
    #        component get copied as well.
    new_component = deepcopy(component)
    if div == 'by':
        new_component.alpha = div_by(new_component.alpha, wrt)
    elif div == 'into':
        new_component.alpha = div_into(new_component.alpha, wrt)
    else:
        raise ValueError('Invalid division specification: %s' % div)
    new_component.xi.add_partial(wrt)
    return new_component


def AR_differential(mvec, wrt, div=DIVISION_TYPE):
    '''
    Compute the result of Differentiating a each component of a MultiVector
    with respect to a given list of unit elements under the algebra.
    '''
    comps = []
    for component in mvec:
        for element in wrt:
            comps.append(component_partial(component, Alpha(element), div))
    return MultiVector(comps)


def Dmu(mvec):
    '''The main operator from the paper'''
    return AR_differential(mvec, ['0', '1', '2', '3'])


def DG(mvec):
    '''A full derivative with respect to all components'''
    return AR_differential(mvec, ALLOWED)
