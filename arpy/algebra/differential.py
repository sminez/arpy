'''
arpy (Absolute Relativity in Python)
Copyright (C) 2016-2017 Innes D. Anderson-Morrison All rights reserved.

A selection of different implementations for symbolically computing
the 4-vector 4-differential Dμ and other Differential operators.

In Cartesian coordinates Dμ is:
    Dμ = ∂μ/αμ = ∂ / (αμ∂xμ) = α0∂0 - αi∂i = α0∂0 - ∇i

All other differential operators follow the same restrictions of
Absolute Relativity and should only operate on MultiVectors.

NOTE:: Specific operators (such as Dmu) are defined in the __init__ file.
'''
from copy import deepcopy
from .config import ALLOWED, DIVISION_TYPE, METRIC
from .ar_types import Alpha
from .multivector import MultiVector
from .operations import div_by, div_into
from .del_grouping import del_grouped


def _div(alpha, wrt, metric, allowed, div):
    '''Divide an alpha component based on the set division type'''
    if div == 'by':
        return div_by(alpha, wrt, metric, allowed)
    elif div == 'into':
        return div_into(wrt, alpha, metric, allowed)
    # NOTE:: Reversed division types divide the differential alpha by that
    #        of the multivector component.
    elif div == 'revby':
        return div_by(wrt, alpha, metric, allowed)
    elif div == 'revinto':
        return div_into(alpha, wrt, metric, allowed)
    else:
        raise ValueError('Invalid division specification: %s' % div)


def component_partial(component, wrt, div, metric, allowed):
    '''
    Symbolically differentiate a component by storing the partials and
    converting the alpha value using the correct division type.
    '''
    # NOTE:: using deep copy so that all of the objects inside of the
    #        component get copied as well.
    new_component = deepcopy(component)
    new_component.alpha = _div(new_component.alpha, wrt, metric, allowed, div)
    new_component.xi.partials = [wrt] + new_component.xi.partials
    return new_component


def AR_differential(mvec, wrt, div=DIVISION_TYPE, metric=METRIC,
                    allowed=ALLOWED, as_del=False):
    '''
    Compute the result of Differentiating a each component of a MultiVector
    with respect to a given list of unit elements under the algebra.
    '''
    comps = []
    for component in mvec:
        for element in wrt:
            element = Alpha(element)
            comp = component_partial(component, element, div, metric, allowed)
            comps.append(comp)
    result = MultiVector(comps)

    if as_del:
        return MultiVector(del_grouped(result))
    else:
        return result


def differential_operator(wrt):
    '''Define a new operator as a function for later use'''
    def operator(mvec, div=DIVISION_TYPE, metric=METRIC,
                 allowed=ALLOWED, as_del=False):
        return AR_differential(mvec, wrt, div, metric, allowed, as_del)
    return operator
