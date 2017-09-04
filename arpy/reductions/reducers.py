'''
arpy (Absolute Relativity in Python)
Copyright (C) 2016-2017 Innes D. Anderson-Morrison All rights reserved.

This module provides helper functions for writing reductions on multivectors.

::NOTE::
It looks like the issue with parsing Force style equations is that the current
Term implementation only looks at partials for the full term, not the
possibility that there may be a XiProduct with nested partials inside.

TODO:
    Second Derivatives
    Product and derivatives (F dF etc)
'''
from itertools import groupby

from ..algebra.ar_types import Xi, Pair, Alpha
from ..algebra.config import config as cfg
from ..utils.utils import SUPER_SCRIPTS, SUB_SCRIPTS


class FailedMatch(Exception):
    pass


class Term:
    def __init__(self, sign, alpha, partials, xis):
        '''
        sign :    + -
        alpha:    [bxyz][group]    (group is a capital letter)
        partials: list of strings
        xis:      list of strings
        '''
        self.sign = sign
        self.alpha_bxyz = alpha[0]
        self.alpha_group = '_' if len(alpha) == 1 else alpha[1]
        self.partials = None if partials is None else tuple(partials)
        self.xis = None if xis is None else tuple(xis)


class Replacement:
    def __init__(self, unallowed, required, termfunc):
        self.unallowed = unallowed
        self.required = required
        self.termfunc = termfunc


class Template:
    def __init__(self, terms=[], replacements=[], can_match=[]):
        self.terms = terms
        self.replacements = replacements
        self.can_match = can_match

    def replace(self, terms, cfg):
        '''Use a template to replace terms in a multivector'''
        self.match_map = {t: [] for t in self.terms}
        self.non_matching = []
        self.bind(terms)
        return self.validate_and_substitute(cfg)

    def match(self, value):
        '''
        Check to see if a given value matches on of the template terms.
        '''
        print('here')
        for term in self.terms:
            # Check that we correctly have either a Xi or a
            # XiProduct with the correct number of components
            if term.xis not in [tuple("_"), None]:
                pairs = ([], [])

                if isinstance(value.xi, Xi):
                    if len(term.xis) > 1:
                        continue
                    pairs = zip(term.xis, [value.xi])
                else:
                    if len(term.xis) != len(value.xi.components):
                        continue
                    pairs = zip(term.xis, value.xi.components)

                # Confirm that bxyz-ness is correct
                ok = (cfg.bxyz_like.get(have.val) == want[0]
                      for want, have in pairs if want[0] in 'bxyz')
                if not all(ok):
                    continue

            # Check that we have the correct kind of alpha
            if term.alpha_bxyz != "_" and term.alpha_bxyz in 'bxyz':
                if cfg.bxyz_like[value.alpha.index] != term.alpha_bxyz:
                    continue

            # Check partials (Xi partials are Alpha objects)
            if term.partials not in [tuple("_"), None]:
                if len(value.xi.partials) == len(term.partials):
                    ok = (cfg.bxyz_like.get(have.index) == want[0]
                          for want, have
                          in zip(term.partials, value.xi.partials))
                    if not all(ok):
                        continue
                else:
                    continue

            return term
        return None

    def bind(self, values):
        '''
        Bind input terms to the individual term templates provided.
        '''
        for v in values:
            try:
                matching_template = self.match(v)
                if matching_template is not None:
                    self.match_map[matching_template].append(v)
                else:
                    self.non_matching.append(v)
            except:
                # If we get any errors then keep the terms as provided
                self.non_matching.append(v)

    def validate_and_substitute(self, cfg):
        '''
        Check that we a consistent matches for the terms in the template. For
        any complete match that we are able to build, remove the matching terms
        and substitute the replacement.
        '''
        output = self.non_matching

        # If there is only one term in the pattern then we are done. Just
        # substitute any matches and return.
        if len(self.terms) == 1:
            for term, matches in self.match_map.items():
                for match in matches:
                    output.append(self.replace(match, term))
            return output

        # Otherwise, take candidates from the first term template and try
        # to build complete matches of the template.
        for t1_candidate in tuple(self.match_map.get(self.terms[0], [])):
            matching_terms = {self.terms[0]: t1_candidate}

            # Get the initial requirements for this match
            try:
                requirements = self.update_match_or_fail(
                    self.terms[0], t1_candidate, {'+_sign': None})
            except FailedMatch:
                continue

            # Try to match each term in the template as we build up the
            # known requirements.
            for term in self.terms[1:]:
                match = None
                for candidate in self.match_map.get(term, []):
                    try:
                        new_requirements = self.update_match_or_fail(
                                term, candidate, requirements)
                        # We need to do this so we don't accidentally mutate
                        # the requirements as part of a failed match.
                        requirements = new_requirements
                        match = candidate
                        break
                    except FailedMatch:
                        # That term broke the match so try again
                        pass

                # If we matched then build up the matching terms
                if match is not None:
                    matching_terms[term] = match
                else:
                    # If we failed to match one of the terms in the template
                    # then we know that the t1_cadidate is not part of a
                    # complete match. However, the remaining terms could still
                    # match another t1_candidate so don't discard them yet
                    output.append(t1_candidate)
                    # reset the requirements
                    # requirements = self.update_match_or_fail(
                    #         self.terms[0], t1_candidate, {'+_sign': None})
                    break
            else:
                # If there were no breaks then we should have a match for each
                # of the terms in the template and we can now substitute.
                replacement = self.generate_replacement(requirements, cfg)
                if replacement is not None:
                    # Some matches may result in the terms cancelling
                    output.append(replacement)

                # Remove any replaced terms from the term lists
                for k, v in matching_terms.items():
                    self.match_map[k].remove(v)

        # Append all remaining terms to the output
        for term in self.terms[1:]:
            output.extend(self.match_map.get(term, []))

        return output

    def update_match_or_fail(self, term, candidate, reqs):
        '''
        Given a term and the current requirements, check to see if we
        can continue. Returns the new requirements.

        NOTE: here, a candiate is a pair.
        '''
        def check_group(want, have, reqs):
            required_group = reqs.get(want)
            a_group = cfg.four_sets[have]

            if required_group is None:
                reqs[want] = a_group
            else:
                if required_group != a_group:
                    raise FailedMatch(
                        '{} != {} ({} {})'.format(
                            want, have, required_group, a_group))

            return reqs

        # It's ok for terms to be negated so long as _all_ terms are
        sign = reqs['+_sign']
        tsign = 1 if term.sign == '+' else -1

        if sign is None:
            reqs['+_sign'] = 1 if tsign == candidate.xi.sign else -1
        else:
            if (tsign * sign) != candidate.xi.sign:
                raise FailedMatch

        # At this stage we have now checked that the bxyz nature of each
        # element of the term is correct, that it has the right number
        # of partials and xi components and that the sign is correct.
        # All we need to do now is confirm that the match groups (capital
        # letters in the patterns) are consistent.
        if term.alpha_group is not '_':
            reqs = check_group(term.alpha_group, candidate.alpha.index, reqs)

        if term.alpha_bxyz is not '_':
            if term.alpha_bxyz not in 'bxyz':
                required = reqs.get(term.alpha_bxyz)
                if required is None:
                    reqs[term.alpha_bxyz] = candidate.alpha.index
                else:
                    if required != candidate.alpha.index:
                        raise FailedMatch

        if term.partials is not None:
            for want, have in zip(term.partials, candidate.xi.partials):
                reqs = check_group(want[1], have.index, reqs)

        if term.xis is not None:
            for want, have in zip(term.xis, candidate.xi.components):
                if len(want) > 1:
                    reqs = check_group(want[1], have.val, reqs)
                else:
                    required = reqs.get(want)
                    if required is None:
                        reqs[want] = have.val
                    else:
                        if required != have.val:
                            raise FailedMatch

        return reqs

    def generate_replacement(self, requirements, cfg):
        '''
        Given the assembled requirements, build the replacement term.
        '''
        for r in self.replacements:
            unallowed_bindings = [requirements[r] for r in r.unallowed]
            if len(set(unallowed_bindings)) != len(unallowed_bindings):
                continue

            required_bindings = [requirements[r] for r in r.required]
            if len(set(required_bindings)) > 1:
                continue

            return r.termfunc(requirements, cfg)
        return None


def cancel_like_terms(terms):
    '''
    For each alpha in the multivector, cancel terms that match their
    negative and return a new multivector of the remaining terms.
    '''
    filtered_pairs = []
    for g in groupby(terms, lambda p: p.alpha):
        alpha, pairs = g
        seen = {}

        for p in pairs:
            if -p.xi in seen:
                # We already have this term's negative so cancel
                if len(seen[-p.xi]) == 1:
                    del seen[-p.xi]
                else:
                    # Cancel one term but leave the rest
                    seen[-p.xi] = seen[-p.xi][1:]
            elif p.xi in seen:
                # We have more than one term the same so append it
                seen[p.xi].append(p)
            else:
                # Store this term for now and use it to check for
                # future cancelling terms
                seen[p.xi] = [p]

        # Add in all of the values we have left
        for v in seen.values():
            filtered_pairs.extend(v)

    return filtered_pairs


# Termfuncs need to take a requirements dict and a config, and return
# a new pair based on the requirements.
def grad_termfunc(reqs, cfg):
    xi = ''.join(SUB_SCRIPTS[x] for x in reqs['k'])
    tex_xi = reqs['k']
    alpha = cfg.alpha_to_group[cfg.four_set_comps[reqs['G']]['x']]

    tex_fourset = '' if reqs['H'] == 'A' else '^' + reqs['H']
    _fourset = '' if reqs['H'] == 'A' else SUPER_SCRIPTS[reqs['H']]
    return Pair(
        Alpha(alpha, reqs['+_sign'], cfg=cfg),
        Xi('∇{}Ξ{}'.format(_fourset, xi),
           tex='\\nabla' + tex_fourset + '\\Xi_{' + tex_xi + '}'),
        cfg=cfg)


def div_termfunc(reqs, cfg):
    proxy_alpha = cfg.four_set_comps[reqs['H']]['x']
    xi = cfg.group_to_4set[cfg.alpha_to_group[proxy_alpha]]
    alpha = cfg.alpha_to_group[cfg.four_set_comps[reqs['F']]['b']]

    tex_fourset = '' if reqs['G'] == 'A' else '^' + reqs['G']
    _fourset = '' if reqs['G'] == 'A' else SUPER_SCRIPTS[reqs['G']]
    return Pair(
        Alpha(alpha, reqs['+_sign'], cfg=cfg),
        Xi('∇{}•{}'.format(_fourset, xi),
           tex='\\nabla{}\\cdot {}'.format(tex_fourset, xi)),
        cfg=cfg)


def curl_termfunc(reqs, cfg):
    proxy_alpha = cfg.four_set_comps[reqs['H']]['x']
    xi = cfg.group_to_4set[cfg.alpha_to_group[proxy_alpha]]
    alpha = cfg.alpha_to_group[cfg.four_set_comps[reqs['F']]['x']]

    tex_fourset = '' if reqs['G'] == 'A' else '^' + reqs['G']
    _fourset = '' if reqs['G'] == 'A' else SUPER_SCRIPTS[reqs['G']]
    return Pair(
        Alpha(alpha, reqs['+_sign'], cfg=cfg),
        Xi('∇{}x{}'.format(_fourset, xi),
           tex='\\nabla{}\\times {}'.format(tex_fourset, xi)),
        cfg=cfg)


def partial_termfunc(reqs, cfg):
    proxy_alpha = cfg.four_set_comps[reqs['H']]['x']
    xi = cfg.group_to_4set[cfg.alpha_to_group[proxy_alpha]]
    alpha = cfg.alpha_to_group[cfg.four_set_comps[reqs['F']]['x']]

    partial = cfg.four_set_comps[reqs['G']]['b']
    _partial = ''.join(SUB_SCRIPTS[b] for b in partial)
    return Pair(
        Alpha(alpha, reqs['+_sign'], cfg=cfg),
        Xi('∂{}{}'.format(_partial, xi),
           tex='\\partial_{}{}'.format(partial, xi)),
        cfg=cfg)


def dot_termfunc(reqs, cfg):
    alpha = cfg.alpha_to_group[cfg.four_set_comps[reqs['F']]['b']]
    return Pair(
        Alpha(alpha, reqs['+_sign'], cfg=cfg),
        Xi('{}•{}'.format(reqs['G'], reqs['H']),
           tex='{}\\cdot {}'.format(reqs['G'], reqs['H'])),
        cfg=cfg)


def wedge_termfunc(reqs, cfg):
    alpha = cfg.alpha_to_group[cfg.four_set_comps[reqs['F']]['x']]
    return Pair(
        Alpha(alpha, reqs['+_sign'], cfg=cfg),
        Xi('{}Λ{}'.format(reqs['G'], reqs['H']),
           tex='{}\\Lambda {}'.format(reqs['G'], reqs['H'])),
        cfg=cfg)


def blade_3vec_termfunc(reqs, cfg):
    b_map = {frozenset('p'): 'p', frozenset('0123'): 'q',
             frozenset('0'): 't', frozenset('123'): 'h'}

    blade = b_map[frozenset(cfg.four_set_comps[reqs['G']]['b'])]
    alpha = cfg.alpha_to_group[cfg.four_set_comps[reqs['F']]['x']]
    return Pair(
        Alpha(alpha, reqs['+_sign'], cfg=cfg),
        Xi('{}{}'.format(blade, reqs['H']),
           tex='{}{}'.format(blade, reqs['H'])),
        cfg=cfg)


def whole_3vec_termfunc(reqs, cfg):
    return Pair(
        Alpha(reqs['k'], reqs['+_sign'], cfg=cfg),
        Xi('{}'.format(reqs['G']),
           tex='{}'.format(reqs['G'])),
        cfg=cfg)


def whole_3vec_squared_termfunc(reqs, cfg):
    return Pair(
        Alpha(reqs['k'], reqs['+_sign'], cfg=cfg),
        Xi('{}²'.format(reqs['G']),
           tex='{}^2'.format(reqs['G'])),
        cfg=cfg)


def dot_square_termfunc(reqs, cfg):
    alpha = cfg.alpha_to_group[cfg.four_set_comps[reqs['F']]['b']]
    return Pair(
        Alpha(alpha, reqs['+_sign'], cfg=cfg),
        Xi('{}²'.format(reqs['G']),
           tex='{}^2'.format(reqs['G'])),
        cfg=cfg)


# Terms are specified according to their component parts.
# Special characters are:
#   b,x,y,z     --> bxyz-like in general
#   bG,xG,yG,zG --> where G is a group (any capital letter)
#   lowercase   --> a Xi value
#   _           --> don't care
# NOTE :: the sets that form the first and second elements of the
#         replacements are the unallowed and required group matches.
grad_template = Template(
    terms=[
        Term('+', 'xG', ('xH',), ('k',)),
        Term('+', 'yG', ('yH',), ('k',)),
        Term('+', 'zG', ('zH',), ('k',))
    ],
    replacements=[Replacement(set(), set(), grad_termfunc)]
)

div_template = Template(
    terms=[
        Term('+', 'bF', ('xG',), ('xH',)),
        Term('+', 'bF', ('yG',), ('yH',)),
        Term('+', 'bF', ('zG',), ('zH',))
    ],
    replacements=[Replacement(set(), set(), div_termfunc)]
)

curl_template = Template(
    terms=[
        Term('+', 'xF', ('yG',), ('zH',)),
        Term('-', 'xF', ('zG',), ('yH',)),
        Term('+', 'yF', ('zG',), ('xH',)),
        Term('-', 'yF', ('xG',), ('zH',)),
        Term('+', 'zF', ('xG',), ('yH',)),
        Term('-', 'zF', ('yG',), ('xH',))
    ],
    replacements=[Replacement(set(), set(), curl_termfunc)]
)

partial_template = Template(
    terms=[
        Term('+', 'xF', ('bG',), ('xH',)),
        Term('+', 'yF', ('bG',), ('yH',)),
        Term('+', 'zF', ('bG',), ('zH',))
    ],
    replacements=[Replacement(set(), set(), partial_termfunc)]
)

dot_template = Template(
    terms=[
        Term('+', 'bF', tuple(), ('xG', 'xH')),
        Term('+', 'bF', tuple(), ('yG', 'yH')),
        Term('+', 'bF', tuple(), ('zG', 'zH'))
    ],
    replacements=[Replacement(set(), set({'G', 'H'}), dot_square_termfunc),
                  Replacement(set(), set(), dot_termfunc)]
)

wedge_template = Template(
    terms=[
        Term('+', 'xF', tuple(), ('yG', 'zH')),
        Term('-', 'xF', tuple(), ('zG', 'yH')),
        Term('+', 'yF', tuple(), ('zG', 'xH')),
        Term('-', 'yF', tuple(), ('xG', 'zH')),
        Term('+', 'zF', tuple(), ('xG', 'yH')),
        Term('-', 'zF', tuple(), ('yG', 'xH'))
    ],
    replacements=[Replacement(set(), set(), wedge_termfunc)]
)

blade_3vec_template = Template(
    terms=[
        Term('+', 'xF', tuple(), ('bG', 'xH')),
        Term('+', 'yF', tuple(), ('bG', 'yH')),
        Term('+', 'zF', tuple(), ('bG', 'zH'))
    ],
    replacements=[Replacement(set(), set(), blade_3vec_termfunc)]
)

# This is to allow for Sb as well as bS (which is caught by the
# template above)
blade_3vec_flipped_template = Template(
    terms=[
        Term('+', 'xF', tuple(), ('xH', 'bG')),
        Term('+', 'yF', tuple(), ('yH', 'bG')),
        Term('+', 'zF', tuple(), ('zH', 'bG'))
    ],
    replacements=[Replacement(set(), set(), blade_3vec_termfunc)]
)

whole_3vec_template = Template(
    terms=[
        Term('+', 'k', tuple(), ('xG',)),
        Term('+', 'k', tuple(), ('yG',)),
        Term('+', 'k', tuple(), ('zG',))
    ],
    replacements=[Replacement(set(), set(), whole_3vec_termfunc)]
)

whole_3vec_squared_template = Template(
    terms=[
        Term('+', 'k', tuple(), ('xG', 'xG')),
        Term('+', 'k', tuple(), ('yG', 'yG')),
        Term('+', 'k', tuple(), ('zG', 'zG'))
    ],
    replacements=[Replacement(set(), set(), whole_3vec_squared_termfunc)]
)


def replace_all(terms, cfg):
    '''Run all known conversions on a Multivector'''
    # Cancel first to cut down the search space
    terms = cancel_like_terms(terms)
    terms = partial_template.replace(terms, cfg)
    terms = grad_template.replace(terms, cfg)
    terms = div_template.replace(terms, cfg)
    terms = curl_template.replace(terms, cfg)
    terms = blade_3vec_template.replace(terms, cfg)
    terms = blade_3vec_flipped_template.replace(terms, cfg)
    terms = dot_template.replace(terms, cfg)
    terms = wedge_template.replace(terms, cfg)
    terms = whole_3vec_squared_template.replace(terms, cfg)
    terms = whole_3vec_template.replace(terms, cfg)
    return terms
