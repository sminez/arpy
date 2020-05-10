from collections import Counter
from copy import copy
from dataclasses import dataclass
from itertools import groupby
from typing import Callable, List, Optional, Set, Tuple

from ..algebra.data_types import Alpha, MultiVector, Term, Xi
from ..config import ARConfig
from ..consts import Orientation, Zet
from ..utils.utils import SUB_SCRIPTS
from .helpers import (
    PartialReplacement,
    alpha_to_group,
    by_alpha_group,
    by_xi,
    by_xi_group,
    filter_by_partials,
    first_partial_str,
    for_all_zets,
    partial_in,
    partial_is,
)


@dataclass
class TaggedCurlTerm:
    term: Term
    xi: str
    is_positive_curl: bool


def as_curl_term(term: Term) -> Optional[TaggedCurlTerm]:
    """(α, ξ, ∂) -> component sign"""

    alpha = Orientation.from_index(term.index).name
    xi = Orientation.from_index(term._components[0].val).name
    partial = Orientation.from_index(first_partial_str(term)).name
    key = f"{alpha}{xi}{partial}"

    if key in ["XYZ", "YZX", "ZXY"]:
        sign = -1
    elif key in ["XZY", "YXZ", "ZYX"]:
        sign = 1
    else:
        return None

    return TaggedCurlTerm(term, term._components[0].val, term.sign == sign)


def del_grouped_terms(mvec: MultiVector) -> List[Term]:
    """
    Group the components of a MultiVector into vector calculus del notation as
    a flat list of Terms using zet grouped Alpha values.
    """
    alpha_grouped = groupby(mvec, lambda t: alpha_to_group(t.index, mvec.cfg))
    output = []

    for group, iter_components in alpha_grouped:
        components = list(iter_components)
        rep_partials, components = for_all_zets(replace_partials, components)
        rep_grad, components = for_all_zets(replace_grad, components)
        rep_div, components = for_all_zets(replace_div, components)
        rep_curl, components = for_all_zets(replace_curl, components)

        output.extend(rep_partials + rep_div + rep_grad + rep_curl)

        for component in components:
            c = copy(component)
            c.alpha.index = alpha_to_group(c.index, mvec.cfg)
            output.append(c)

    return output


def del_grouped(mvec: MultiVector):
    """Print a MultiVector using vector calculus del notation"""
    print("{")

    for alpha, terms in groupby(del_grouped_terms(mvec), lambda t: t._alpha):
        xis = " ".join([term._repr_no_alpha(count=count) for term, count in Counter(terms).items()])
        if xis.startswith("+"):
            xis = xis[2:]

        print(f"  {repr(alpha).ljust(5)}( {xis} )")

    print("}")


def replace_partials(terms: List[Term], zet: Zet) -> PartialReplacement:
    """Partial F = d<comp>F"""
    if len(terms) == 0:
        return [], []

    cfg = terms[0].cfg
    elements = zet.elements(cfg).all
    replaced = []

    for alpha_group, group_terms in groupby(sorted(terms, key=by_alpha_group), by_alpha_group):
        for blade in elements:
            candidates = [t for t in group_terms if partial_is(t, blade)]

            if len(candidates) != 3 or len(set(c.sign for c in candidates)) > 1:
                continue

            have = sorted([c._components[0].val for c in candidates])
            comp_zet = Zet.from_index(have[0])
            needed = sorted(comp_zet.elements(cfg).all[1:])

            if have != needed:
                continue

            sign = candidates[0].sign
            blade = "".join(SUB_SCRIPTS[b] for b in blade)

            replaced.append(
                Term(
                    Alpha(alpha_group, sign, cfg=cfg),
                    [Xi(f"∂{blade}{comp_zet.name}", tex=f"\\partial_{blade}{comp_zet.name}")],
                    cfg=cfg,
                )
            )

            for candidate in candidates:
                terms.remove(candidate)

    return replaced, terms


def replace_div(terms: List[Term], zet: Zet) -> PartialReplacement:
    """Div F = dFx/dx + dFy/dy + dFz/dz"""
    if len(terms) == 0:
        return [], []

    cfg = terms[0].cfg
    elements = zet.elements(cfg).all[1:]
    replaced: List[Term] = []
    candidates = []

    for c in filter_by_partials(terms, elements):
        val_orientation = Orientation.from_index(c._components[0].val)
        partial_orientation = Orientation.from_index(first_partial_str(c))
        if val_orientation == partial_orientation:
            candidates.append(c)

    if len(candidates) != 3 or len(set(c.sign for c in candidates)) > 1:
        return replaced, terms

    sign = candidates[0].sign
    alpha_group = alpha_to_group(candidates[0].index, cfg)
    xi = Zet.from_index(candidates[0]._components[0].val).name

    # The 'A' 3Vector calculus operators are the standard ones
    _zet = zet.superscript
    tex_zet = zet.name if zet != Zet.A else ""

    replaced.append(
        Term(
            Alpha(alpha_group, sign, cfg=cfg),
            [Xi(f"∇{_zet}•{xi}", tex=f"\\nabla_{tex_zet}\\cdot {xi}")],
            cfg=cfg,
        )
    )

    for candidate in candidates:
        terms.remove(candidate)

    return replaced, terms


def replace_grad(terms: List[Term], zet: Zet) -> PartialReplacement:
    """Grad f = αx[df/dx] + αy[df/dy] + αz[df/dz]"""
    if len(terms) == 0:
        return [], []

    cfg = terms[0].cfg
    elements = zet.elements(cfg).all[1:]
    replaced = []

    for xi, group_terms in groupby(sorted(filter_by_partials(terms, elements), key=by_xi), by_xi):
        candidates = [t for t in group_terms]
        consistent_sign = len(set(c.sign for c in candidates)) == 1
        correct_partials = sorted(first_partial_str(c) for c in candidates) == sorted(elements)

        if not (len(candidates) == 3 and consistent_sign and correct_partials):
            continue

        sign = candidates[0].sign
        alpha_group = alpha_to_group(candidates[0].index, cfg)
        xi = Zet.from_index(candidates[0]._components[0].val).time_like

        tex_zet = "" if zet == Zet.A else f"_{{{zet.name}}}"
        _zet = zet.superscript

        replaced.append(
            Term(
                Alpha(alpha_group, sign, cfg=cfg),
                [Xi(f"∇Ξ{_zet}{xi}", tex=f"\\nabla{tex_zet}\\{xi}")],
                cfg=cfg,
            )
        )

        for candidate in candidates:
            terms.remove(candidate)

    return replaced, terms


def replace_curl(terms: List[Term], zet: Zet) -> PartialReplacement:
    """Curl F = αx[dFz/dy-dFy/dz] + αy[dFx/dz-dFz/dx] + αz[dFy/dx-dFx/dy]"""
    if len(terms) == 0:
        return [], []

    cfg = terms[0].cfg
    elements = zet.elements(cfg).all[1:]
    curl_like = [
        c for t in filter_by_partials(terms, elements) if (c := as_curl_term(t)) is not None
    ]
    replaced = []

    for group, group_terms in groupby(curl_like, lambda c: alpha_to_group(c.xi, cfg)):
        candidates = [t for t in group_terms]

        if len(candidates) != 6 or len(set(c.is_positive_curl for c in candidates)) > 1:
            continue

        sign = 1 if candidates[0].is_positive_curl else -1
        alpha_group = alpha_to_group(candidates[0].term.index, cfg)
        xi = Zet.from_index(candidates[0].term._components[0].val).name

        _zet = zet.superscript
        tex_zet = "" if zet == "A" else f"^{zet.name}"

        replaced.append(
            Term(
                Alpha(alpha_group, sign, cfg=cfg),
                [Xi(f"∇{_zet}x{xi}", tex=f"\\nabla{tex_zet}\\times {group}")],
                cfg=cfg,
            )
        )

        for term in [c.term for c in candidates]:
            terms.remove(term)

    return replaced, terms
