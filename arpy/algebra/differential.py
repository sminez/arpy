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
from .ar_types import Alpha, MultiVector
from .operations import div_by, div_into


def _div(alpha, wrt, metric, div):
    '''Divide an alpha component based on the set division type'''
    if div == 'by':
        return div_by(alpha, wrt, metric)
    elif div == 'into':
        return div_into(alpha, wrt, metric)
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


def _4set_differential(mvec, wrt, div, metric):
    '''Shortcut differential that prints the result in del notation'''
    def index_and_sign(comp_ix, diff_ix):
        alpha = _div(Alpha(comp_ix), Alpha(diff_ix), metric, div)
        return alpha.index, alpha.sign

    PARTIAL = '∂{}{}'
    GRAD = '∇Ξ{}'
    DIV = '∇•{}'
    CURL = '∇x{}'

    differential_components = []

    alphas = [a for a in ALLOWED if Alpha(a) in mvec]
    mvec_blade_4sets, mvec_3vec_4sets = _find_4sets(alphas)
    wrt_blade_4sets, wrt_3vec_4sets = _find_4sets(wrt)

    # Compute the action of any paired blades on 4-Sets
    for diff_4set in wrt_blade_4sets:
        diff_index = FOUR_SET_COMPS[diff_4set]['b']

        # blade|blade: ∂{paired blade}{paired blade}
        for mvec_4set in mvec_blade_4sets:
            component_index = FOUR_SET_COMPS[mvec_4set]['b']
            index, sign = index_and_sign(component_index, diff_index)
            group_index = ALPHA_TO_GROUP[index]
            replacement_xi = PARTIAL.format(diff_index, component_index)
            differential_components.append(
                (Alpha(group_index, sign), replacement_xi)
            )

        # blade vec: ∂{paired blade}{3vector}
        for mvec_4set in mvec_3vec_4sets:
            result_4set = _4set_result(diff_4set, mvec_4set)
            component_index = FOUR_SET_COMPS[mvec_4set]['x']
            index, sign = index_and_sign(component_index, diff_index)
            group_index = ALPHA_TO_GROUP[index]
            replacement_xi = PARTIAL.format(diff_index, result_4set)
            differential_components.append(
                (Alpha(group_index, sign), replacement_xi)
            )

    # vec blade: ∇Ξ{paired blade}
    for diff_4set in wrt_3vec_4sets:
        diff_index = FOUR_SET_COMPS[diff_4set]['y']
        for mvec_4set in mvec_blade_4sets:
            component_index = FOUR_SET_COMPS[mvec_4set]['b']
            index, sign = index_and_sign(component_index, diff_index)
            group_index = ALPHA_TO_GROUP[index]
            differential_components.append(
                (Alpha(group_index, sign), GRAD.format(component_index))
            )

        # vec vec: ∇•{3vector} & ∇x{3vector}
        for mvec_4set in mvec_3vec_4sets:
            div_index = FOUR_SET_COMPS[mvec_4set]['y']
            index, sign = index_and_sign(div_index, diff_index)
            group_index = ALPHA_TO_GROUP[index]
            differential_components.append(
                (Alpha(group_index, sign), DIV.format(mvec_4set))
            )

            curl_index = FOUR_SET_COMPS[mvec_4set]['x']
            index, sign = index_and_sign(curl_index, diff_index)
            group_index = ALPHA_TO_GROUP[index]
            differential_components.append(
                (Alpha(group_index, sign), CURL.format(mvec_4set))
            )

    return differential_components
