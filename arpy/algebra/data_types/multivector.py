from collections import defaultdict
from copy import deepcopy
from itertools import groupby
from typing import List, Union

from ..config import ARConfig
from ..config import config as cfg
from .alpha import Alpha
from .term import Term

# Custom type union to allow MultiVectors to be initialised using either
# A list of terms, a list of strings or a single string.
TermsOrStrings = Union[List[Term], List[str], str]


class MultiVector:
    """
    A MultiVector is an unordered collection of a Terms representing a particular
    composite quantity within the Algebra. In its simplest form, a MultiVector is
    a simple linear sum of Alphas, though it is possible for there to be significantly
    more structure.

    In practice, almost all arpy computations are done using MultiVectors as their
    primary data structure so there are a number of methods designed for aiding in
    simplifying such computations.

    NOTE: The __invert__ method on MultiVecors is defined in __init__.py
          as it requires the use of the full product function which results
          in a cyclic import.
    """

    def __init__(self, terms: TermsOrStrings = [], cfg: ARConfig = cfg):
        if isinstance(terms, str):
            terms = [Term(t) for t in terms.split()]

        _terms = []

        for t in terms:
            if isinstance(t, (Alpha, str)):
                t = Term(t, cfg=cfg)

            if not isinstance(t, Term):
                raise ValueError("Arguments must be Terms or strings")

            if t.index not in cfg.allowed:
                raise ValueError(f"Invalid alpha ({t.alpha}): allowed values are {cfg.allowed}")

            _terms.append(t)

        self._terms = _terms
        self.cfg = cfg
        self.__ensure_standard_form()

    def __eq__(self, other):
        if not isinstance(other, MultiVector):
            return False

        return sorted(self._terms) == sorted(other._terms)

    def __len__(self):
        return len(self._terms)

    def __add__(self, other):
        if isinstance(other, Term):
            terms = self._terms + [other]
        elif isinstance(other, MultiVector):
            terms = self._terms + other._terms
        else:
            raise TypeError()

        return MultiVector(terms, cfg=self.cfg)

    def __sub__(self, other):
        return self + -other

    def __neg__(self):
        terms = deepcopy(self._terms)
        for t in terms:
            t._sign *= -1

        return MultiVector(terms, cfg=self.cfg)

    def __contains__(self, other):
        if isinstance(other, Term):
            return other in self._terms
        elif isinstance(other, Alpha):
            return other in set(t._alpha for t in self._terms)

        return False

    def __getitem__(self, key):
        key = self.__ensure_key_is_alpha(key)
        terms = [t for t in self._terms if t._alpha == key]
        return MultiVector(terms, cfg=self.cfg)

    def __delitem__(self, key):
        key = self.__ensure_key_is_alpha(key)
        self._terms = list(filter(lambda t: t._alpha != key, self._terms))

    def __iter__(self):
        yield from self._terms

    def __repr__(self):
        rep = []
        for alpha, terms in groupby(self._terms, lambda t: t._alpha):
            xis = " ".join(t._repr_no_alpha() for t in terms)
            if xis.startswith("+ "):
                xis = xis[2:]
            rep.append(f"  {repr(alpha).ljust(5)}( {xis} )")

        return "\n".join(["{"] + rep + ["}"])

    def __ensure_key_is_alpha(self, key):
        """
        Helper to allow for shorthand strings to be used in place of Alphas.
        This is an instance method to allow us to also enforce using the correct
        config (allowed/metric) when constructing the Alpha version of the key.
        """
        if isinstance(key, str):
            key = Alpha(key, cfg=self.cfg)

        if not isinstance(key, Alpha):
            raise KeyError()

        return key

    def __ensure_standard_form(self):
        """
        Ensure that the ordering of the terms in this MultiVector are in standard
        form ordering and that all term cancellations have been carried out.
        """
        if len(self._terms) == 0:
            return

        seen = defaultdict(list)

        for term in sorted(self._terms):
            neg = -term
            if neg in seen:
                if len(seen[neg]) == 1:
                    del seen[neg]
                else:
                    seen[neg] = seen[neg][1:]
            else:
                seen[term].append(term)

        self._terms = sum(seen.values(), [])

    def iter_alphas(self):
        """
        Iterate over the contents of a MultiVector by Alpha yielding tuples
        of the Alpha and a list of Terms. The iteration order for the Alphas
        is defined to be the same as the order specified in the ARConfig used
        to create this MultiVector.
        """
        groups = {k: list(v) for k, v in groupby(self._terms, lambda t: t._alpha)}

        for alpha in self.cfg.allowed:
            key = Alpha(alpha, cfg=self.cfg)
            terms = groups.get(key)
            if terms:
                yield key, terms

    # =================================================== #
    # Alternative string representations for MultiVectors
    # =================================================== #

    @property
    def with_factored_terms(self):
        rep = []

        for alpha, terms in groupby(self._terms, lambda t: t._alpha):
            rep.append(f"  {(repr(alpha) + ':').ljust(5)}")
            for factor, others in groupby(terms, lambda t: t._components[0]):
                factored = f"      {repr(factor).ljust(2)}"
                xis = " ".join(t._repr_no_alpha(ix=1) for t in others)

                if xis.startswith("+ "):
                    xis = xis[2:]
                if xis:
                    factored += f".({xis})"

                rep.append(factored)
            rep.append("")

        return "\n".join(["{"] + rep + ["}"])
