'''
A selection of different implementations for symbolically computing
the 4-vector 4-differential Dμ.

In Cartesian coordinates this is:
    Dμ = ∂μ/αμ = ∂ / (αμ∂xμ) = α0∂0 - αi∂i = α0∂0 - ∇i
'''
from itertools import groupby

from ..config import ALLOWED, α_TO_GROUP
from .algebra import α, ξ, ξα


def _partial(component, wrt):
    '''
    Symbolically differentiate a component by storing the partials and
    converting the alpha value using the correct division type.
    '''
    alpha = component.alpha / wrt
    partials = component.xi.partials + [wrt]
    return ξα(alpha, ξ(component.xi.val, component.xi.unit, partials))


def partial(vec, wrt):
    '''
    Symbolically differentiate a whole vector of ξα pairs
    '''
    return [_partial(comp, wrt) for comp in vec]


def Dμ(vec):
    '''
    Compute the result of [Dμ vec] under the algebra. this result is a
    single list of terms which almost always will get passed through
    `by_α` or `by_vec_notation` to make it more readable.
    '''
    grouped = [partial(vec, α(a)) for a in ['0', '1', '2', '3']]
    output = []
    for g in grouped:
        output += g
    return output


##############################################################################


def by_α(vec, vector_groups=False):
    '''
    Group results of a derivative operation together in order to form the
    resulting coupled differential equations.
    '''
    vec = sorted(vec, key=lambda a: ALLOWED.index(a.alpha.index))

    if vector_groups:
        return groupby(vec, lambda x: α_TO_GROUP[x.alpha.index])
    else:
        return groupby(vec, key=lambda x: x.alpha.index)
