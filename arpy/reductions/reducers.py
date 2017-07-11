'''
arpy (Absolute Relativity in Python)
Copyright (C) 2016-2017 Innes D. Anderson-Morrison All rights reserved.

This module provides helpr functions for writing reductions on multivectors.


TODO:
    Simplified Div/Grad/Curl/Partials
    Second Derivatives
    Cross, Dot, Wedge products
'''
from itertools import groupby

from ..algebra.ar_types import Xi, Pair, Alpha
from ..algebra.config import config as cfg
from ..algebra.multivector import MultiVector
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
        for term in self.terms:
            # Check that we correctly have either a Xi or a
            # XiProduct with the correct number of components
            if term.xis not in [tuple("_"), None]:
                pairs = None

                if isinstance(value.xi, Xi):
                    if len(term.xis) > 1:
                        continue
                    pairs = zip(term.xis, [value.xi.val])
                else:
                    if len(term.xis) != len(value.xi.components):
                        continue
                    pairs = zip(term.xis, value.xi.components)

                # Confirm that bxyz-ness is correct
                if pairs:
                    for want, have in pairs:
                        if have[0] in 'bxyz':
                            if cfg.bxyz_like.get(have.val) != want:
                                continue

            # Check that we have the correct kind of alpha
            if term.alpha_bxyz != "_":
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
            return MultiVector(output)

        # Otherwise, take candidates from the first term template and try
        # to build complete matches of the template.
        for t1_candidate in tuple(self.match_map.get(self.terms[0], [])):
            requirements = {'+_sign': None}
            matching_terms = {self.terms[0]: t1_candidate}

            # Get the initial requirements for this match
            requirements = self.update_match_or_fail(
                    self.terms[0], t1_candidate, requirements)

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
                    requirements = self.update_match_or_fail(
                            self.terms[0], t1_candidate, {'+_sign': None})
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
        def check_group(want, have):
            required_group = reqs.get(want)
            a_group = cfg.four_sets[have]

            if required_group is None:
                reqs[want] = a_group
            else:
                if required_group != a_group:
                    raise FailedMatch

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
            reqs = check_group(term.alpha_group, candidate.alpha.index)

        if term.partials is not None:
            for want, have in zip(term.partials, candidate.xi.partials):
                reqs = check_group(want[1], have.index)

        if term.xis is not None:
            for want, have in zip(term.xis, candidate.xi.components):
                if len(want) > 1:
                    reqs = check_group(want[1], have.val)
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

    return MultiVector(filtered_pairs)


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


def replace_all(terms, cfg):
    '''Run all known conversions on a Multivector'''
    terms = partial_template.replace(terms, cfg)
    terms = grad_template.replace(terms, cfg)
    terms = div_template.replace(terms, cfg)
    terms = curl_template.replace(terms, cfg)
    return terms
