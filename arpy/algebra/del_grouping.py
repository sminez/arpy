from itertools import groupby
from collections import namedtuple
from .config import ALPHA_TO_GROUP, FOUR_SET_COMPS, FOUR_SETS, BXYZ_LIKE, \
        SUPER_SCRIPTS, SUB_SCRIPTS, GROUP_TO_4SET
from .ar_types import Alpha, Pair, DelMultiVector


term = namedtuple('term', ['d', 'xi', 'alpha', 'sign', 'pair'])

# (α, ξ, ∂) -> component sign
CURL_SIGN = {
    ('x', 'z', 'y'): 1, ('x', 'y', 'z'): -1,
    ('y', 'x', 'z'): 1, ('y', 'z', 'x'): -1,
    ('z', 'y', 'x'): 1, ('z', 'x', 'y'): -1
}


def _filter_on_partials(terms, partials, n=0):
    '''n is the index for xi.partials'''
    filtered_terms = [
        term(t.xi.partials[n].index, t.xi.val, t.alpha.index, t.xi.sign, t)
        for t in terms if t.xi.partials[n].index in partials
    ]
    return sorted(filtered_terms, key=lambda t: ALPHA_TO_GROUP[t.xi])


def _present_4sets(pairs):
    '''
    Return all four sets with at least one partial derivative component
    in `pairs`
    '''
    return set([FOUR_SETS[p.xi.partials[0].index] for p in pairs])


def del_grouped(mvec):
    '''
    Group the components of a multivector into del notation
    '''
    output = []
    alpha_grouped = groupby(mvec, lambda p: ALPHA_TO_GROUP[p.alpha.index])
    for group, components in alpha_grouped:
        components = [c for c in components]
        replaced, components = replace_curl(components)
        output.extend(replaced)
        replaced, components = replace_grad(components)
        output.extend(replaced)
        replaced, components = replace_div(components)
        output.extend(replaced)
        replaced, components = replace_partials(components)
        output.extend(replaced)

        for component in components:
            component.alpha.index = ALPHA_TO_GROUP[component.alpha.index]
        output.extend(components)
    return DelMultiVector(output)


def replace_curl(pairs):
    '''Curl F = αx[dFz/dy-dFy/dz] + αy[dFx/dz-dFz/dx] + αz[dFy/dx-dFx/dy]'''
    replaced = []

    for fourset in _present_4sets(pairs):
        comps = [FOUR_SET_COMPS[fourset][k] for k in ['x', 'y', 'z']]
        candidates = []

        for t in _filter_on_partials(pairs, comps):
            aix, xix, dix = [BXYZ_LIKE[k] for k in [t.alpha, t.xi, t.d]]
            curl_sign = CURL_SIGN.get((aix, xix, dix))
            if curl_sign:
                candidates.append({'term': t, 'curl_sign': curl_sign})

        grouped_cands = groupby(candidates,
                                lambda t: ALPHA_TO_GROUP[t['term'].xi])

        for group, cands in grouped_cands:
            cands = [c for c in cands]
            if len(cands) != 6:
                continue

            if all([c['curl_sign'] == c['term'].sign for c in cands]):
                sign = 1
            elif all([c['curl_sign'] == -c['term'].sign for c in cands]):
                sign = -1
            else:
                continue

            _fourset = '' if fourset == 'A' else SUPER_SCRIPTS[fourset]
            alpha = ALPHA_TO_GROUP[cands[0]['term'].alpha]
            group = GROUP_TO_4SET[group]
            replaced.append(
                Pair(Alpha(alpha, sign), '∇{}x{}'.format(_fourset, group))
            )

            for candidate in [c['term'].pair for c in cands]:
                pairs.remove(candidate)

    return replaced, pairs


def replace_grad(pairs):
    '''Grad f = αx[df/dx] + αy[df/dy] + αz[df/dz]'''
    replaced = []
    for fourset in _present_4sets(pairs):
        comps = [FOUR_SET_COMPS[fourset][k] for k in ['x', 'y', 'z']]

        candidates = _filter_on_partials(pairs, comps)
        sorted_candidates = sorted(candidates, key=lambda t: t.xi)
        grouped = groupby(sorted_candidates, lambda t: t.xi)

        for xi, candidates in grouped:
            candidates = [c for c in candidates]
            if len(candidates) != 3:
                continue

            if set(c.d for c in candidates) == set(comps):
                if all(c.sign == 1 for c in candidates):
                    sign = 1
                elif all(c.sign == -1 for c in candidates):
                    sign = -1
                else:
                    continue

                xi = candidates[0].xi
                alpha = ALPHA_TO_GROUP[candidates[0].alpha]
                _fourset = '' if fourset == 'A' else SUPER_SCRIPTS[fourset]
                replaced.append(
                    Pair(Alpha(alpha, sign), '∇{}Ξ{}'.format(_fourset, xi))
                )

                for candidate in candidates:
                    pairs.remove(candidate.pair)
    return replaced, pairs


def replace_div(pairs):
    '''Div F = dFx/dx + dFy/dy + dFz/dz'''
    replaced = []
    for fourset in _present_4sets(pairs):
        comps = [FOUR_SET_COMPS[fourset][k] for k in ['x', 'y', 'z']]

        # partial from this 4set and x-partial with x-3vec component, etc
        candidates = [
            c for c in _filter_on_partials(pairs, comps)
            if BXYZ_LIKE[c.xi] == BXYZ_LIKE[c.d]
        ]

        if len(candidates) == 3:
            alpha = candidates[0].alpha
            xi = GROUP_TO_4SET[ALPHA_TO_GROUP[candidates[0].xi]]
            # The 'A' 3Vector calculus operators are the standard ones
            _fourset = '' if fourset == 'A' else SUPER_SCRIPTS[fourset]

            if all(c.sign == 1 for c in candidates):
                sign = 1
            elif all(c.sign == -1 for c in candidates):
                sign = -1
            else:
                continue

            replaced.append(
                Pair(Alpha(alpha, sign), '∇{}•{}'.format(_fourset, xi))
            )
            for candidate in candidates:
                pairs.remove(candidate.pair)
    return replaced, pairs


def replace_partials(pairs):
    ''' Partial F = d{comp}F'''
    def key(p):
        return ALPHA_TO_GROUP[p.alpha.index]

    replaced = []

    for fourset in _present_4sets(pairs):
        for _, grouped_pairs in groupby(sorted(pairs, key=key), key):
            grouped_pairs = [p for p in grouped_pairs]
            for blade in FOUR_SET_COMPS[fourset].values():
                candidates = [
                    p for p in grouped_pairs
                    if p.xi.partials[0].index == blade
                ]
                if len(candidates) == 3:
                    ix = candidates[0].xi.val
                    component_4set = FOUR_SETS[ix]
                    needed = {
                        FOUR_SET_COMPS[component_4set][k]
                        for k in ['x', 'y', 'z']
                    }
                    have = {c.xi.val for c in candidates}
                    if have != needed:
                        continue

                    if all([c.xi.sign == 1 for c in candidates]):
                        sign = 1
                    elif all([c.xi.sign == -1 for c in candidates]):
                        sign = -1
                    else:
                        continue

                    alpha = ALPHA_TO_GROUP[candidates[0].alpha.index]
                    blade = SUB_SCRIPTS[blade]

                    replaced.append(
                        Pair(Alpha(alpha, sign),
                             '∂{}{}'.format(blade, component_4set)))
                    for candidate in candidates:
                        pairs.remove(candidate)
    return replaced, pairs
