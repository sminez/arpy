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
from .config import config as cfg
from ..utils.utils import SUB_SCRIPTS
from .ar_types import Alpha
from .multivector import MultiVector, DelMultiVector
from .operations import div_by, div_into, inverse, full


class AR_differential:
    '''Differential operator: can be used inside of ar()'''
    def __init__(self, wrt, cfg=cfg):
        if isinstance(wrt, MultiVector):
            self.wrt = [pair.alpha for pair in wrt]
        else:
            if isinstance(wrt, str):
                wrt = wrt.split()

            if isinstance(wrt, list):
                # Conversion to Alpha catches invalid indices
                self.wrt = [Alpha(comp, cfg=cfg) for comp in wrt]
            else:
                raise ValueError(
                    'Differential operators must be initialised with either'
                    ' a MultiVector, list or string of alpha indices')

        alphas = ', '.join([str(a) for a in self.wrt])
        self.__doc__ = 'Differnetiate with respect to: {}'.format(alphas)

    def __call__(self, mvec, cfg=cfg, div=None, as_del=False):
        '''
        Compute the result of Differentiating a each component of a MultiVector
        with respect to a given list of unit elements under the algebra.
        '''
        comps = []
        for comp in mvec:
            for element in self.wrt:
                result = component_partial(comp, element, cfg, div)
                comps.append(result)
        derivative = MultiVector(comps, cfg=cfg)
        derivative.replacements.extend(mvec.replacements)

        if as_del:
            return DelMultiVector(derivative)
        else:
            return derivative

    def __repr__(self):
        elements = [
            '{}∂{}'.format(
                str(inverse(a)),
                ''.join(SUB_SCRIPTS[i] for i in a.index)
            )
            for a in self.wrt
        ]
        return '{ ' + ' '.join(elements) + ' }'


def _div(alpha, wrt, cfg, div=None):
    '''Divide an alpha component based on the set division type'''
    div = div if div else cfg.division_type
    if div == 'by':
        return div_by(alpha, wrt, cfg)
    elif div == 'into':
        return div_into(wrt, alpha, cfg)
    else:
        raise ValueError(
            'Invalid division specification: {}'.format(cfg.division_type))


def component_partial(component, wrt, cfg, div):
    '''
    Symbolically differentiate a component by storing the partials and
    converting the alpha value using the correct division type.
    '''
    # NOTE:: using deep copy so that all of the objects inside of the
    #        component get copied as well.
    new_component = deepcopy(component)
    new_component.alpha = _div(new_component.alpha, wrt, cfg, div)
    new_component.xi.partials = [wrt] + new_component.xi.partials
    return new_component


def differential_operator(wrt, cfg=cfg):
    '''Define a new operator as a function for later use'''
    return AR_differential(wrt, cfg=cfg)


@full.add((AR_differential, MultiVector))
def _full_differential_mvec(diff, mvec, cfg=cfg):
    res = diff(mvec, cfg=cfg)
    return res


@full.add((MultiVector, AR_differential))
def _full_mvec_differential_mvec(mvec, diff, cfg=cfg):
    res = diff(mvec, cfg=cfg, div='by')
    return res
