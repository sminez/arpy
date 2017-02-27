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
from .del_grouping import del_grouped


def _div(alpha, wrt, metric, div):
    '''Divide an alpha component based on the set division type'''
    if div == 'by':
        return div_by(alpha, wrt, metric)
    elif div == 'into':
        return div_into(wrt, alpha, metric)
    else:
        raise ValueError('Invalid division specification: %s' % div)


def component_partial(component, wrt, div, metric):
    '''
    Symbolically differentiate a component by storing the partials and
    converting the alpha value using the correct division type.
    '''
    # NOTE:: using deep copy so that all of the objects inside of the
    #        component get copied as well.
    new_component = deepcopy(component)
    new_component.alpha = _div(new_component.alpha, wrt, metric, div)
    new_component.xi.partials = [wrt] + new_component.xi.partials
    return new_component


def AR_differential(mvec, wrt, div=DIVISION_TYPE, metric=METRIC, as_del=False):
    '''
    Compute the result of Differentiating a each component of a MultiVector
    with respect to a given list of unit elements under the algebra.
    '''
    comps = []
    for component in mvec:
        for element in wrt:
            comp = component_partial(component, Alpha(element), div, metric)
            comps.append(comp)

    result = MultiVector(comps)
    if as_del:
        return del_grouped(result)
    else:
        return result


def differential_operator(wrt):
    '''Define a new operator as a function for later use'''
    def operator(mvec, div=DIVISION_TYPE, metric=METRIC, as_del=False):
        return AR_differential(mvec, wrt, div, metric, as_del)
    return operator


##############################################################################
# Sepcific operators #
######################
def Dmu(mvec, div=DIVISION_TYPE, metric=METRIC, as_del=False):
    '''The main operator from the paper'''
    return AR_differential(mvec, ['0', '1', '2', '3'], div, metric, as_del)


def DG(mvec, div=DIVISION_TYPE, metric=METRIC, as_del=False):
    '''A full derivative with respect to all components'''
    return AR_differential(mvec, ALLOWED, div, metric, as_del)
