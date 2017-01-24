'''
A selection of different implementations for symbolically computing
the 4-vector 4-differential Dμ.

In Cartesian coordinates this is:
    Dμ = ∂μ/αμ = ∂ / (αμ∂xμ) = α0∂0 - αi∂i = α0∂0 - ∇i
'''
from itertools import groupby

from ..config import ALLOWED, ALPHA_TO_GROUP
from .algebra import a, xi, Pair


def _partial(component, wrt):
    '''
    Symbolically differentiate a component by storing the partials and
    converting the alpha value using the correct division type.
    '''
    alpha = component.alpha / wrt
    partials = component.xi.partials + [wrt]
    return Pair(alpha, xi(component.xi.val, component.xi.unit, partials))


def partial(vec, wrt):
    '''
    Symbolically differentiate a whole vector of ξα pairs
    '''
    return [_partial(comp, wrt) for comp in vec]


def _differential(vec, wrt=None):
    '''
    Compute the result of Differentiating a xi vector with respect to a given
    list of unit elements under the algebra. This result is a
    single list of terms which almost always will get passed through
    `by_α` or `by_vec_notation` to make it more readable.

    This is used as base to build other differential operators.
    '''
    grouped = [partial(vec, a(comp)) for comp in wrt]
    output = []
    for g in grouped:
        output += g
    return output


def del_i(vec):
    '''The 3D gradient of a vector'''
    return _differential(vec, ['1', '2', '3'])


def Dmu(vec):
    '''The main operator from the paper'''
    return _differential(vec, ['0', '1', '2', '3'])


def DG(vec):
    '''A full derivative with respect to all components'''
    return _differential(vec, ALLOWED)


def DM(vec):
    '''Differentiate with respect to the magnetic field'''
    M = ['p'] + [a for a in ALLOWED if len(a) == 2 and '0' not in a]
    return _differential(vec, M)


def DE(vec):
    '''Differentiate with respect to the electric field'''
    E = ['0123'] + [a for a in ALLOWED if len(a) == 2 and '0' in a]
    return _differential(vec, E)


##############################################################################


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
