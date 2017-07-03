'''
arpy (Absolute Relativity in Python)
Copyright (C) 2016-2017 Innes D. Anderson-Morrison All rights reserved.

This module provides helpr functions for writing reductions on multivectors.


TODO:
    Simplified Div/Grad/Curl/Partials
    Second Derivatives
    Cross, Dot, Wedge products
'''
from collections import namedtuple

from ..algebra.ar_types import Alpha, Xi, XiProduct, Pair
from ..algebra.multivector import MultiVector
from ..algebra.config import config as cfg
from ..utils.utils import SUPER_SCRIPTS, SUB_SCRIPTS
from ..utils.concepts.prelude import groupby, tee

term = namedtuple('term', ['d', 'xi', 'alpha', 'sign', 'pair'])
patrn = namedtuple('pattern', ['d', 'xi', 'alpha', 'sign'])
subs = namedtuple('subs', ['aix', 'asign', 'xval', 'tex'])


def make_term(pair, neg=False):
    '''
    Convert a Pair to a term in order to have easier access to attributes.
    Passing neg=True allows you to make the negative of this term which is
    useful for checking for cancelling terms in multivectors.
    NOTE: In an Mvec, all alpha signs are +ve which is why we are using the
          sign of the Xi.
    '''
    if not isinstance(pair.xi, (Xi, XiProduct)):
        raise ValueError('Non Xi xi value: {}'.format(pair.xi))

    sign = pair.xi.sign
    if neg:
        sign *= -1

    return term(d=tuple(pair.xi.partials), xi=pair.xi,
                alpha=pair.alpha.index, sign=sign, pair=pair)


def cancel_like_terms(mvec):
    '''
    For each alpha in the multivector, cancel terms that match their
    negative and return a new multivector of the remaining terms.
    '''
    filtered_pairs = []
    for g in groupby(mvec, lambda p: p.alpha):
        alpha, pairs = g
        seen = {}

        for p in pairs:
            if -p.xi in seen:
                # We already have this term's negative so cancel both
                del seen[-p.xi]
            else:
                # Store this term for now and use it to check for
                # future cancelling terms
                seen[p.xi] = p

        filtered_pairs.extend(list(seen.values()))
    return MultiVector(filtered_pairs)


def substitute(pattern, replacement=None):
    '''
    Search for a given pattern in an input stream and replace
    the terms with an alternative. If the replacement is None
    then the matched terms will just be dropped.

    The pattern
    '''
    def new_filter(stream):
        pass

    return new_filter


def filter_on_partials(terms, partials, cfg=cfg, n=0):
    '''
    For a stream of terms, return only those that contain a partial
    derivative from partials at position n.
    '''
    def npartial(t, n):
        return t.xi.partials[n].index

    filtered_terms = (
        term(npartial(t, n), t.xi.val, t.alpha.index, t.xi.sign, t)
        for t in terms
        if t.xi.partials != [] and npartial(t, n) in partials
        and not isinstance(t.xi.val, tuple)
    )

    return sorted(filtered_terms, key=lambda t: cfg.alpha_to_group[t.xi])


def present_4sets_derivatives(pairs, cfg, n=0):
    '''
    Return the set of all four sets with at least one partial derivative
    component in `pairs`
    '''
    return {cfg.four_sets[p.xi.partials[n].index]
            for p in pairs if p.xi.partials != []}
