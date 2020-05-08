"""
arpy (Absolute Relativity in Python)
Copyright (C) 2016-2018 Innes D. Anderson-Morrison All rights reserved.

NOTE:: To avoid cyclic imports, the __invert__ method on multivectors (which
       returns the Hermitian conjugate) is written in the arpy __init__ file!
       (I am fully aware of how horrible this is...)
"""
import collections.abc
import re
from collections import Counter, namedtuple
from copy import deepcopy
from itertools import groupby

from ..reductions.reducers import del_grouped, replace_all
from ..utils.utils import Nat
from .ar_types import Alpha, Pair, Xi, XiProduct
from .config import config as cfg

MvecLabel = namedtuple("MvecLabel", "label originals")
explicit_xi = r"[-p0123]*\[.*\]"


class MultiVector(collections.abc.Set):
    """
    A custom container type for working efficiently with multivectors.
    """

    # NOTE: The __invert__ method on MultiVecors is defined in __init__.py
    #       as it requires the use of the full product function which results
    #       in a cyclic import.

    def __init__(self, components=[], cfg=cfg):
        # Given a list of pairs, build the mulitvector by binding the ξ values
        self.cfg = cfg
        self.components = {Alpha(a, cfg=cfg): [] for a in cfg.allowed}
        self.replacements = []

        if isinstance(components, str):  # Allow for single string input
            components = components.split()

        for comp in components:
            if isinstance(comp, Alpha):
                comp = Pair(comp, cfg=cfg)
            elif isinstance(comp, str):
                if re.match(explicit_xi, comp):
                    # This is something of the form 012[Sin(kx-ωt)] so split
                    # it into `012` and `Sin(kx-ωt)`
                    alpha, xi = comp.split("[")
                    # Remove the trailing `]`
                    comp = Pair(alpha, xi[:-1], cfg=cfg)
                else:
                    comp = Pair(comp, cfg=cfg)

            if not isinstance(comp, Pair):
                raise ValueError("Arguments must be Alphas, Pairs or Strings")

            if comp.alpha.index in cfg.allowed:
                _comp = deepcopy(comp)
                try:
                    self.components[_comp.alpha].append(_comp.xi)
                except KeyError:
                    # Negative Alpha value
                    alpha, xi = _comp.alpha, _comp.xi
                    alpha.sign = 1
                    xi.sign *= -1
                    self.components[alpha].append(xi)
            else:
                raise ValueError("Invalid alpha: allowed values are {}".format(self.cfg.allowed))

    def __repr__(self):
        comps = [
            "  {}{}".format(repr(a).ljust(5), self._nice_xi(a))
            for a in (Alpha(x, cfg=self.cfg) for x in self.cfg.allowed)
            if self.components[a]
        ]
        return "{\n" + "\n".join(comps) + "\n}"

    def __tex__(self):
        comps = "\n".join(
            ("  \\alpha_{" + str(a) + "}").ljust(17)
            + self._nice_xi(Alpha(a, cfg=self.cfg), tex=True)
            + r"+ \nonumber\\"
            for a in self.cfg.allowed
            if self.components[Alpha(a, cfg=self.cfg)]
        )
        return r"\{ \nonumber\\" + "\n" + comps + "\n\\}" + r"\nonumber\\"

    def show(self, ordering):
        """Print the components of the MultiVector in a specified ordering"""
        if isinstance(ordering, str):
            ordering = ordering.split()

        if not all([o in self.cfg.allowed for o in ordering]):
            raise ValueError("Invalid index in ordering")

        comps = [
            "  α{}{}".format(str(a).ljust(5), self._nice_xi(Alpha(a, cfg=cfg)))
            for a in ordering
            if self.components[Alpha(a, cfg=cfg)]
        ]
        print("{\n" + "\n".join(comps) + "\n}")

    def __len__(self):
        # All allowed values are initialised with [] so we are
        # only counting componets that have a Xi value set.
        return sum([len(v) for v in self.components.values()])

    def __add__(self, other):
        # Allow for the addition of Multivectors and Pairs.
        # This will always return a new MultiVector.
        if not isinstance(other, (Pair, MultiVector)):
            raise TypeError()

        comps = [p for p in self]
        if isinstance(other, Pair):
            comps.append(other)
        elif isinstance(other, MultiVector):
            comps.extend(p for p in other)

        res = MultiVector(comps, cfg=self.cfg)
        return res

    def __sub__(self, other):
        res = self + -other
        return res

    def __eq__(self, other):
        if not isinstance(other, MultiVector):
            return False

        if self.cfg != other.cfg:
            # Multivectors from different algebras are never equal
            return False

        for alpha in self.components:
            if self.components[alpha] != other.components[alpha]:
                return False
        return True

    def __neg__(self):
        # -mvec to negate all Xis
        negated = deepcopy(self)
        for xis in negated.components.values():
            for xi in xis:
                xi.sign *= -1
        return negated

    def __contains__(self, other):
        if isinstance(other, Alpha):
            return self.components[other] != []
        elif isinstance(other, Pair):
            return other.xi in self.components[other.alpha]
        else:
            return False

    def __getitem__(self, key):
        if isinstance(key, str):
            # Allow retreval by string as well as Alpha
            key = Alpha(key, cfg=self.cfg)
        if not isinstance(key, Alpha):
            raise KeyError
        return [Pair(key, xi) for xi in self.components[key]]

    def __delitem__(self, key):
        if isinstance(key, str):
            key = Alpha(key, cfg=self.cfg)
        if not isinstance(key, Alpha):
            raise KeyError

        self.components[key] = []

    def __iter__(self):
        for alpha in self.cfg.allowed:
            key = Alpha(alpha, cfg=self.cfg)
            try:
                for xi in self.components[key]:
                    yield Pair(alpha, xi, cfg=self.cfg)
            except KeyError:
                pass

    def iter_alphas(self):
        """
        Iterate over the contents of a MultiVector by alpha. This method
        yields tuples of the Alpha and a list of Pairs.
        """
        for alpha in self.cfg.allowed:
            key = Alpha(alpha, cfg=self.cfg)
            try:
                pairs = [Pair(alpha, xi, cfg=self.cfg) for xi in self.components[key]]
                if pairs:
                    yield key, pairs
            except KeyError:
                pass

    def cancel_terms(self):
        """Remove terms that cancel following a calculation"""
        # XXX: optimise! This is quadratic in len(components.values())
        #      Maybe groupby would work?
        for xis in self.components.values():
            i = 0
            while True:
                try:
                    current = xis[i]
                    rest = xis[i + 1 :]
                    for r in rest:
                        if current == -r:
                            xis.remove(current)
                            xis.remove(r)
                            break
                    else:
                        # Triggered if we don't break
                        i += 1
                except IndexError:
                    break

    def _nice_xi(self, alpha, raise_key_error=False, for_print=True, tex=False):
        """Single element xi lists return their value raw"""

        def n_xi(comps):
            '''Convert [x, x] -> "2x"'''
            c = Counter([str(comp) for comp in comps])

            str_comps = []
            for k, v in c.items():
                if k.startswith("-"):
                    if v > 1:
                        str_comps.append("- %s%s" % (v, k[1:]))
                    else:
                        str_comps.append("- %s" % k[1:])
                else:
                    if v > 1:
                        str_comps.append("+ %s%s" % (v, k))
                    else:
                        str_comps.append("+ %s" % k)
            return " ".join(str_comps)

        try:
            xi = sorted(self.components[alpha], key=lambda x: x.partials)
        except KeyError:
            if raise_key_error:
                raise KeyError
        if len(xi) == 1:
            if tex:
                return "( " + xi[0].__tex__() + " )"
            else:
                return "( " + str(xi[0]) + " )"
        else:
            if for_print:
                if tex:
                    return "( " + " ".join(x.__tex__() for x in xi) + " )"
                else:
                    return "( " + n_xi(xi) + " )"
            else:
                return xi

    def copy(self):
        """Return a copy of this multivector"""
        return deepcopy(self)

    def zet_grouped(self):
        """
        Print an zet grouped representation of the MultiVector
        NOTE:: This deliberately does not return a new MultiVector as we
               should always be working with strict alpha values not grouped.
        """
        by_alpha = groupby(self, key=lambda x: self.cfg.alpha_to_group[x.alpha.index])
        zet = [(group, tuple(c.xi for c in components)) for (group, components) in by_alpha]
        print("{")
        for group, comps in zet:
            print("  α{}".format(group).ljust(7), comps)
        print("}")

    def del_notation(self):
        """
        Print a del grouped version of the multivector
        """
        return DelMultiVector(self, cfg=self.cfg)

    def simplified(self, ix=0, exyz=False, sign=False):
        """
        Display the multivector with simplified Xi values.
        """

        def key(x):
            return x.alpha.index

        ix2 = 1 if ix == 0 else 1
        alpha_grouped = groupby(sorted(self, key=key), key)
        seen = []

        try:
            print("{")
            for _, full in alpha_grouped:
                full = sorted(full, key=lambda x: x.xi.components[ix])
                xi_grouped = groupby(full, key=lambda x: x.xi.components[ix])
                for common_xi, comps in xi_grouped:
                    comps = [
                        Pair(
                            c.alpha,
                            Xi(c.xi.components[ix2].val, c.xi.components[ix2].partials, c.xi.sign),
                        )
                        for c in comps
                    ]
                    comps = sorted(comps, key=key)
                    for component in comps:
                        if component.alpha.sign == -1:
                            component.alpha.sign = 1
                            component.xi.sign *= -1

                    if comps[0].alpha not in seen:
                        print("  {}:".format(comps[0].alpha))
                        seen.append(comps[0].alpha)

                    if exyz:
                        x = Nat(common_xi.val).ljust(3)
                        if sign:
                            signs = [c.xi.exyz()[0] for c in comps]
                            blocks = ["■" if s == "-" else "□" for s in signs]
                            comp_str = " ".join(blocks)
                        else:
                            comp_str = ", ".join([c.xi.exyz() for c in comps])
                    else:
                        x = str(common_xi).ljust(6)
                        comp_str = ", ".join([str(c.xi) for c in comps])
                    print("    {}( {} )".format(x, comp_str))
            print("}")
        except:
            print("Unable to simplify MultiVector Xi values")

    def relabel(self, index, replacement):
        """
        Manually relabel the components of a multivector. This is intended
        for simplifying results following manual analysis and returns a new
        MultiVector when called.
        """
        new_mvec = deepcopy(self)

        if index.startswith("-"):
            index = index[1:]
            xi_sign = -1
        else:
            xi_sign = 1

        if replacement.startswith("-"):
            replacement = replacement[1:]
            xi_sign *= -1

        if index in self.cfg.allowed:
            new_xi = Xi(replacement, sign=xi_sign)
            originals = deepcopy(new_mvec.components[Alpha(index, cfg=cfg)])
            replacements = [MvecLabel(new_xi, originals)]
            new_mvec.components[Alpha(index, cfg=cfg)] = [new_xi]
        else:
            try:
                replacements = []
                indices = zip(["₁", "₂", "₃"], self.cfg.xi_groups[index])
                for comp, index in indices:
                    new_xi = Xi(replacement + comp, sign=xi_sign)
                    originals = deepcopy(new_mvec.components[Alpha(index, cfg=cfg)])
                    replacements.append(MvecLabel(new_xi, originals))
                    new_mvec.components[Alpha(index, cfg=cfg)] = [new_xi]
            except KeyError:
                raise ValueError("{} is not a valid index".format(index))

        new_mvec.replacements.extend(replacements)
        return new_mvec

    def relabel_many(self, pairs):
        """Relabel each case in pairs: (index, replacement)"""
        new_mvec = deepcopy(self)
        for index, replacement in pairs:
            new_mvec = new_mvec.relabel(index, replacement)
        return new_mvec

    def remove_labels(self):
        """Generate a new MultiVector without the current replacements"""

        def replace_or_keep(alpha, xi, return_pairs=True):
            for r in self.replacements:
                if xi.val == r.label.val:
                    vals = deepcopy(r.originals)
                    new_comps = []
                    for comp in vals:
                        comp.partials.extend(xi.partials)
                        comp.sign *= xi.sign
                        if return_pairs:
                            new_comps.append(Pair(alpha, comp, cfg=cfg))
                        else:
                            new_comps.append(comp)
                    return new_comps
            if return_pairs:
                return [Pair(alpha, xi, cfg=cfg)]
            else:
                return [xi]

        if self.replacements == []:
            return deepcopy(self)

        new_comps = []

        for pair in self:
            if isinstance(pair.xi, Xi):
                new_comps.extend(replace_or_keep(pair.alpha, pair.xi))

            elif isinstance(pair.xi, XiProduct):
                new_pair_components = []
                for component in pair.xi.components:
                    new_pair_components.extend(replace_or_keep(pair.alpha, component, False))
                pair.xi.components = tuple(new_pair_components)
                new_comps.append(pair)

        return MultiVector(new_comps, cfg=cfg)

    @property
    def v(self):
        return self.del_notation()

    @property
    def s(self, ix=0, exyz=False, sign=False):
        self.simplified(ix, exyz, sign)


class DelMultiVector(MultiVector):
    def __init__(self, components=[], cfg=cfg):
        # Given a list of pairs, build the mulitvector by binding the ξ values
        self.cfg = cfg
        self._raw_components = components

        self.components = {Alpha(a, cfg=cfg): [] for a in self.cfg.allowed_groups}

        for comp in del_grouped(self._raw_components, cfg=self.cfg):
            if isinstance(comp, (str, Alpha)):
                comp = Pair(comp, cfg=cfg)
            if not isinstance(comp, Pair):
                raise ValueError("Arguments must be Alphas, Pairs or Strings")
            if comp.alpha.index in self.cfg.allowed_groups:
                try:
                    self.components[comp.alpha].append(comp.xi)
                except KeyError:
                    # Negative Alpha value
                    alpha, xi = comp.alpha, comp.xi
                    alpha.sign = 1
                    xi.sign *= -1
                    self.components[alpha].append(xi)

    def __eq__(self, other):
        if isinstance(other, MultiVector):
            other = GroupedMultiVector(other)
        if not isinstance(other, DelMultiVector):
            return False
        for alpha in self.components:
            if self.components[alpha] != other.components[alpha]:
                return False
        return True

    def __repr__(self):
        comps = [
            "  {}{}".format(repr(a).ljust(5), self._nice_xi(a))
            for a in (Alpha(x, cfg=self.cfg) for x in self.cfg.allowed_groups)
            if self.components[a]
        ]
        return "{\n" + "\n".join(comps) + "\n}"

    def __tex__(self):
        comps = "\n".join(
            ("  \\alpha_{" + str(a) + "}").ljust(17)
            + self._nice_xi(Alpha(a, cfg=self.cfg), tex=True)
            + r"+ \nonumber\\"
            for a in self.cfg.allowed_groups
            if self.components[Alpha(a, cfg=self.cfg)]
        )

        return r"\{ \nonumber\\" + "\n" + comps + "\n\\}" + r"\nonumber\\"


class GroupedMultiVector(MultiVector):
    def __init__(self, components=[], cfg=cfg):
        # Given a list of pairs, build the mulitvector by binding the ξ values
        self.cfg = cfg
        if isinstance(components, (MultiVector, DelMultiVector, GroupedMultiVector)):
            cfg = components.cfg
            components = [p for p in components]

        self._raw_components = components
        self.alphas = cfg.allowed + cfg.allowed_groups[-4:]

        self.components = {Alpha(a, cfg=cfg): [] for a in self.alphas}

        for comp in replace_all(self._raw_components, cfg=self.cfg):
            if isinstance(comp, (str, Alpha)):
                comp = Pair(comp, cfg=cfg)
            if not isinstance(comp, Pair):
                raise ValueError("Arguments must be Alphas, Pairs or Strings")
            if comp.alpha.index in self.alphas:
                try:
                    self.components[comp.alpha].append(comp.xi)
                except KeyError:
                    # Negative Alpha value
                    alpha, xi = comp.alpha, comp.xi
                    alpha.sign = 1
                    xi.sign *= -1
                    self.components[alpha].append(xi)

    def __eq__(self, other):
        if not isinstance(other, GroupedMultiVector):
            return False
        for alpha in self.components:
            if self.components[alpha] != other.components[alpha]:
                return False
        return True

    def __repr__(self):
        comps = [
            "  {}{}".format(repr(a).ljust(5), self._nice_xi(a))
            for a in (Alpha(x, cfg=self.cfg) for x in self.alphas)
            if self.components[a]
        ]
        return "{\n" + "\n".join(comps) + "\n}"

    def __tex__(self):
        comps = "\n".join(
            ("  \\alpha_{" + str(a) + "}").ljust(17)
            + self._nice_xi(Alpha(a, cfg=self.cfg), tex=True)
            + r"+ \nonumber\\"
            for a in self.alphas
            if self.components[Alpha(a, cfg=self.cfg)]
        )

        return r"\{ \nonumber\\" + "\n" + comps + "\n\\}" + r"\nonumber\\"
