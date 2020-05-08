"""
This module implements a pattern matching algorithm for use with AR calculation
terms. This is intended to be a generic engine that is extensible as we identify
further patterns and groupings with higher order meaning.
"""
from dataclasses import dataclass
from typing import Callable, List

from ..algebra.data_types import Alpha, Term, Xi
from ..config import config as cfg
from ..consts import Orientation, Zet


class FailedMatch(Exception):
    pass


@dataclass
class Target:
    sign: int
    zet: Zet
    orientation: Orientation
    xis: List[Xi]
    partials: List[Alpha]


@dataclass
class Replacement:
    unallowed: List
    required: List
    term_func: Callable


class Template:
    def __init__(self, targets: List[Target], replacements: List[Replacement]) -> "Template":
        self.targets = targets
        self.replacements = replacements

        self.non_matching = []
        self.match_map = {t: [] for t in self.targets}

    def bind(self, terms: List[Term]):
        """Bind input terms to the individual term templates provided"""
        for t in terms:
            try:
                target = self.match(t)
                if target is not None:
                    self.match_map[target].append(t)
                else:
                    self.non_matching.append(t)

            except Exception:
                # If we get any errors then keep the terms as provided
                self.non_matching.append(t)

    # TODO: Update and rewrite
    def match(self, term):
        """Check to see if a given term matches on of our targets"""
        for target in self.targets:
            # Check that we correctly have either a Xi or a
            # XiProduct with the correct number of components
            if term.xis not in [tuple("_"), None]:
                pairs = ([], [])

                if isinstance(term.xi, Xi):
                    if len(term.xis) > 1:
                        continue
                    pairs = zip(term.xis, [term.xi])
                else:
                    if len(term.xis) != len(term.xi.components):
                        continue
                    pairs = zip(term.xis, term.xi.components)

                # Confirm that exyz-ness is correct
                ok = (Nat(have.val) == want[0] for want, have in pairs if want[0] in "exyz")
                if not all(ok):
                    continue

            # Check that we have the correct kind of alpha
            if term.alpha_exyz != "_" and term.alpha_exyz in "exyz":
                if Nat(value.alpha) != term.alpha_exyz:
                    continue

            # Check partials (Xi partials are Alpha objects)
            if term.partials not in [tuple("_"), None]:
                if len(value.xi.partials) == len(term.partials):
                    ok = (
                        Nat(have) == want[0] for want, have in zip(term.partials, value.xi.partials)
                    )
                    if not all(ok):
                        continue
                else:
                    continue

            return term
        return None

    def rewrite(self, terms, config):
        """Rewrite the given terms using the matchings captured by this template"""
        self.bind(terms)

        return self.validate_and_substitute(config)
