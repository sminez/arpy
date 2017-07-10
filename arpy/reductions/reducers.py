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

from ..algebra.ar_types import Xi
from ..algebra.config import config as cfg
from ..algebra.multivector import MultiVector
# from ..utils.utils import SUPER_SCRIPTS, SUB_SCRIPTS


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


class Template:
    def __init__(self, terms=[], replacements=[], can_match=[]):
        self.terms = terms
        self.replacements = replacements
        self.can_match = can_match

    def replace(self, values):
        '''Use a template to replace terms in a multivector'''
        self.match_map = {t: [] for t in self.terms}
        self.non_matching = []
        self.bind(values)
        return self.validate_and_substitute()

    def match(self, value):
        '''
        Check to see if a given value matches on of the template terms.
        '''
        for term in self.terms:
            # Check that we correctly have either a Xi or a XiProduct with
            # the correct number of components
            if term.xis != tuple("_"):
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
                        if term[0] in 'bxyz':
                            if cfg.bxyz_like.get(have.val) != want:
                                continue

            # Check that we have the correct kind of alpha
            if term.alpha != "_":
                if cfg.bxyz_like[value.alpha.index] != term.alpha[0]:
                    continue

            # Check partials (Xi partials are Alpha objects)
            if term.partials != tuple("_"):
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
            matching_template = self.match(v)
            if matching_template is not None:
                self.match_map[matching_template].append(v)
            else:
                self.non_matching.append(v)

    def validate_and_substitute(self):
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
        for t1_candidate in self.match_map.get(self.terms[0], []):
            requirements = {}
            matching_terms = {self.terms[0]: t1_candidate}

            # Get the initial requirements for this match
            requirements = self.update_match_or_fail(
                    t1_candidate, requirements)

            # Try to match each term in the template as we build up the
            # known requirements.
            for term in self.terms[1:]:
                match = None
                for candidate in self.match_map.get(term, []):
                    new_requirements = self.update_match_or_fail(
                            term, candidate, requirements)
                    if new_requirements is not None:
                        requirements = new_requirements
                        match = candidate
                        break

                # If we matched then build up the matching terms
                if match is not None:
                    matching_terms[term] = match
                else:
                    # If we failed to match one of the terms in the template
                    # then we know that the t1_cadidate is not part of a
                    # complete match. However, the remaining terms could still
                    # match another t1_candidate so don't discard them yet
                    output.append(t1_candidate)
                    break
            else:
                # If there were no breaks then we should have a match for each
                # of the terms in the template and we can now substitute.
                replacement = self.generate_replacement(requirements)
                if replacement is not None:
                    # Some matches may result in the terms cancelling
                    output.append(replacement)

                # Remove any replaced terms from the term lists
                for key, term in matching_terms:
                    self.terms[key].remove(term)

        # Append all remaining terms to the output
        for term in self.terms[1:]:
            output.extend(self.match_map.get(term, []))

        return MultiVector(output)

    def update_match_or_fail(self, term, candidate, reqs):
        '''
        Given a term and the current requirements, check to see if we
        can continue. Returns the new requirements.

        NOTE: here, a candiate is a pair.
        '''
        # It's ok for terms to be negated so long as _all_ terms are
        sign = reqs.get('+')
        if sign is None:
            csign = '+' if candidate.xi.sign == 1 else '-'
            reqs['+'] = 1 if term.sign == csign else -1
        else:
            if sign != candidate.xi.sign:
                return None

        # At this stage we have now checked that the bxyz nature of each
        # element of the term is correct, that it has the right number
        # of partials and xi components and that the sign is correct.
        # All we need to do now is confirm that the match groups (capital
        # letters in the patterns) are consistent.
        # TODO :: This is trickier than it sounds!

        return None

    def generate_replacement(self, requirements):
        '''
        Given the assembled requirements, build the replacement term.
        '''
        pass


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


# Terms are specified according to their component parts.
# Special characters are:
#   b,x,y,z     --> bxyz-like in general
#   bG,xG,yG,zG --> where G is a group (any capital letter)
#   lowercase   --> a Xi value
#   _           --> don't care
grad_template = Template(
    terms=[
        Term('+', 'xG', ('xH',), ('k',)),
        Term('+', 'yG', ('yH',), ('k',)),
        Term('+', 'zG', ('zH',), ('k',))
    ],
    replacements=[Term('+', 'G', None, ('∇^{H}Ξk',))],
    can_match=[{'G', 'H'}],
)

div_template = Template(
    terms=[
        Term('+', 'bF', ('xG',), ('xH',)),
        Term('+', 'bF', ('yG',), ('yH',)),
        Term('+', 'bF', ('zG',), ('zH',))
    ],
    replacements=[Term('+', 'G', None, ('∇^{G}•Ξk',))],
    can_match=[{'F', 'G', 'H'}],
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
    replacements=[Term('+', 'G', None, ('∇^{G}xΞH',))],
    can_match=[{'F', 'G'}, {'F', 'H'}],
)
