'''
A selection of different implementations for symbolically computing
the 4-vector 4-differential Dμ and other Differential operators.

In Cartesian coordinates Dμ is:
    Dμ = ∂μ/αμ = ∂ / (αμ∂xμ) = α0∂0 - αi∂i = α0∂0 - ∇i

All other differential operators follow the same restrictions of
Absolute Relativity and should only operate on MultiVectors.
'''
from copy import deepcopy
from .config import ALLOWED, DIVISION_TYPE, METRIC, FOUR_SET_COMPS, \
        ALPHA_TO_GROUP
from .ar_types import Alpha, Pair, MultiVector, DelMultiVector
from .operations import div_by, div_into


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
    if as_del:
        return _4set_differential(mvec, wrt, div, metric)
    comps = []
    for component in mvec:
        for element in wrt:
            comp = component_partial(component, Alpha(element), div, metric)
            comps.append(comp)
    return MultiVector(comps)


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


##############################################################################
# 4set derivatives #
####################
def _4set_result(left, right):
    '''Return the 4set name that is produced by left acting on right'''
    if left == 'B':
        return right
    elif right == 'B':
        return left
    elif left == right:
        return 'B'
    else:
        # Helper dict for mapping compositions to resultant sets
        _mappings = [(('A', 'E'), 'T'), (('A', 'T'), 'E'), (('E', 'T'), 'A')]
        _4setmap = {frozenset(lr): res for lr, res in _mappings}
        return _4setmap[frozenset([left, right])]


def _find_4sets(alphas):
    '''Retuns a 3-tuple of full, timelike and 3vector 4set component names'''
    paired_blades = []
    three_vecs = []
    matched = []

    for s, comps in FOUR_SET_COMPS.items():
        if comps['b'] in alphas:
            paired_blades.append(s)
            matched.append(comps['b'])

        _3vec = [comps[i] for i in ['x', 'y', 'z']]
        if all(c in alphas for c in _3vec):
            three_vecs.append(s)
            matched.extend(_3vec)

    if sorted(matched) != sorted(alphas):
        raise ValueError("MultiVector is not composed of 4Set elements")

    return paired_blades, three_vecs


def replace_partial(diff_4set, mvec_4set, component_blade, metric, div):
    diff_index = FOUR_SET_COMPS[diff_4set]['b']
    component_index = FOUR_SET_COMPS[mvec_4set][component_blade]
    alpha = _div(Alpha(component_index), Alpha(diff_index), metric, div)
    group_index = ALPHA_TO_GROUP[alpha.index]
    if component_blade == 'b':
        replacement_xi = '∂{}Ξ{}'.format(diff_index, component_index)
    else:
        replacement_xi = '∂{}{}'.format(diff_index, mvec_4set)
    return Pair(Alpha(group_index, alpha.sign), replacement_xi)


def replace_grad(diff_4set, mvec_4set, metric, div):
    diff_index = FOUR_SET_COMPS[diff_4set]['y']
    component_index = FOUR_SET_COMPS[mvec_4set]['b']
    alpha = _div(Alpha(component_index), Alpha(diff_index), metric, div)
    group_index = ALPHA_TO_GROUP[alpha.index]
    return Pair(Alpha(group_index, alpha.sign), '∇Ξ{}'.format(component_index))


def replace_div(diff_4set, mvec_4set, metric, div):
    diff_index = FOUR_SET_COMPS[diff_4set]['y']
    div_index = FOUR_SET_COMPS[mvec_4set]['y']
    alpha = _div(Alpha(div_index), Alpha(diff_index), metric, div)
    group_index = ALPHA_TO_GROUP[alpha.index]
    return Pair(Alpha(group_index, alpha.sign), '∇•{}'.format(mvec_4set))


def replace_curl(diff_4set, mvec_4set, metric, div):
    diff_index = FOUR_SET_COMPS[diff_4set]['y']
    curl_index = FOUR_SET_COMPS[mvec_4set]['x']
    alpha = _div(Alpha(curl_index), Alpha(diff_index), metric, div)
    group_index = ALPHA_TO_GROUP[alpha.index]
    return Pair(Alpha(group_index, alpha.sign), '∇x{}'.format(mvec_4set))


def _4set_differential(mvec, wrt, div, metric):
    '''Shortcut differential that prints the result in del notation'''
    components = []

    alphas = [a for a in ALLOWED if Alpha(a) in mvec]
    mvec_blade_4sets, mvec_3vec_4sets = _find_4sets(alphas)
    wrt_blade_4sets, wrt_3vec_4sets = _find_4sets(wrt)

    # Compute the action of any paired blades on 4-Sets
    for diff_4set in wrt_blade_4sets:
        # blade|blade: ∂{paired blade}{paired blade}
        for mvec_4set in mvec_blade_4sets:
            components.append(
                replace_partial(diff_4set, mvec_4set, 'b', metric, div)
            )
        # blade|vec: ∂{paired blade}{3vector}
        for mvec_4set in mvec_3vec_4sets:
            components.append(
                replace_partial(diff_4set, mvec_4set, 'x', metric, div)
            )
    # Compute the action of any 3vectors on the 4set
    for diff_4set in wrt_3vec_4sets:
        # vec|blade: ∇Ξ{paired blade}
        for mvec_4set in mvec_blade_4sets:
            components.append(
                replace_grad(diff_4set, mvec_4set, metric, div)
            )
        # vec|vec: ∇•{3vector} & ∇x{3vector}
        for mvec_4set in mvec_3vec_4sets:
            components.extend([
                replace_div(diff_4set, mvec_4set, metric, div),
                replace_curl(diff_4set, mvec_4set, metric, div)
            ])

    return DelMultiVector(components)
