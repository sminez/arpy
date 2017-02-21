'''
A selection of different implementations for symbolically computing
the 4-vector 4-differential Dμ and other Differential operators.

In Cartesian coordinates Dμ is:
    Dμ = ∂μ/αμ = ∂ / (αμ∂xμ) = α0∂0 - αi∂i = α0∂0 - ∇i

All other differential operators follow the same restrictions of
Absolute Relativity and should only operate on MultiVectors.
'''
from copy import deepcopy
from .config import ALLOWED, DIVISION_TYPE, METRIC
from .ar_types import Alpha, MultiVector
from .operations import div_by, div_into


def component_partial(component, wrt, div, metric):
    '''
    Symbolically differentiate a component by storing the partials and
    converting the alpha value using the correct division type.
    '''
    # NOTE:: using deep copy so that all of the objects inside of the
    #        component get copied as well.
    new_component = deepcopy(component)
    if div == 'by':
        new_component.alpha = div_by(new_component.alpha, wrt, metric)
    elif div == 'into':
        new_component.alpha = div_into(new_component.alpha, wrt, metric)
    else:
        raise ValueError('Invalid division specification: %s' % div)
    # new_component.xi.add_partial(wrt)
    new_component.xi.partials = [wrt] + new_component.xi.partials

    return new_component


def AR_differential(mvec, wrt, div=DIVISION_TYPE, metric=METRIC):
    '''
    Compute the result of Differentiating a each component of a MultiVector
    with respect to a given list of unit elements under the algebra.
    '''
    comps = []
    for component in mvec:
        for element in wrt:
            comp = component_partial(component, Alpha(element), div, metric)
            comps.append(comp)
    return MultiVector(comps)


def differential_operator(wrt):
    '''Define a new operator as a function for later use'''
    def operator(mvec, div=DIVISION_TYPE, metric=METRIC):
        return AR_differential(mvec, wrt, div, metric)
    return operator


def Dmu(mvec, div=DIVISION_TYPE, metric=METRIC):
    '''The main operator from the paper'''
    return AR_differential(mvec, ['0', '1', '2', '3'], div, metric)


def DG(mvec, div=DIVISION_TYPE, metric=METRIC):
    '''A full derivative with respect to all components'''
    return AR_differential(mvec, ALLOWED, div, metric)
