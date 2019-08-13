"""
arpy (Absolute Relativity in Python)
Copyright (C) 2016-2018 Innes D. Anderson-Morrison All rights reserved.
"""
from itertools import groupby
from collections import namedtuple

from .reducers import alpha_to_group
from ..algebra.ar_types import Alpha, Xi, Pair
from ..algebra.config import config as cfg
from ..utils.utils import SUPER_SCRIPTS, SUB_SCRIPTS, Zet, Nat


term = namedtuple("term", ["d", "xi", "alpha", "sign", "pair"])

# (α, ξ, ∂) -> component sign
CURL_SIGN = {
    ("x", "z", "y"): 1,
    ("x", "y", "z"): -1,
    ("y", "x", "z"): 1,
    ("y", "z", "x"): -1,
    ("z", "y", "x"): 1,
    ("z", "x", "y"): -1,
}


def _filter_on_partials(terms, partials, cfg, n=0):
    """n is the index for xi.partials"""
    filtered_terms = [
        term(t.xi.partials[n].index, t.xi.val, t.alpha.index, t.xi.sign, t)
        for t in terms
        if t.xi.partials and t.xi.partials[n].index in partials and not isinstance(t.xi.val, tuple)
    ]
    return sorted(filtered_terms, key=lambda t: alpha_to_group(t.xi))


def _present_zets(pairs, cfg):
    """
    Return all four sets with at least one partial derivative component
    in `pairs`
    """
    return set([Zet(p.xi.partials[0]) for p in pairs if p.xi.partials])


def del_grouped(mvec, cfg=cfg):
    """
    Group the components of a multivector into del notation
    """
    output = []
    alpha_grouped = groupby(mvec, lambda p: alpha_to_group(p.alpha.index))
    for group, components in alpha_grouped:
        components = [c for c in components]
        replaced, components = replace_partials(components, cfg)
        output.extend(replaced)
        replaced, components = replace_grad(components, cfg)
        output.extend(replaced)
        replaced, components = replace_div(components, cfg)
        output.extend(replaced)
        replaced, components = replace_curl(components, cfg)
        output.extend(replaced)

        for component in components:
            component.alpha.index = alpha_to_group(component.alpha.index)
        output.extend(components)
    return output


def replace_curl(pairs, cfg):
    """Curl F = αx[dFz/dy-dFy/dz] + αy[dFx/dz-dFz/dx] + αz[dFy/dx-dFx/dy]"""
    replaced = []

    for zet in _present_zets(pairs, cfg):
        comps = [cfg.zet_comps[zet][k] for k in ["x", "y", "z"]]
        candidates = []

        for t in _filter_on_partials(pairs, comps, cfg):
            aix, xix, dix = [Nat(k) for k in [t.alpha, t.xi, t.d]]
            curl_sign = CURL_SIGN.get((aix, xix, dix))
            if curl_sign:
                candidates.append({"term": t, "curl_sign": curl_sign})

        grouped_cands = groupby(candidates, lambda t: alpha_to_group(t["term"].xi))

        for group, cands in grouped_cands:
            cands = [c for c in cands]
            if len(cands) != 6:
                continue

            if all([c["curl_sign"] == c["term"].sign for c in cands]):
                sign = 1
            elif all([c["curl_sign"] == -c["term"].sign for c in cands]):
                sign = -1
            else:
                continue

            alpha = alpha_to_group(cands[0]["term"].alpha)
            group = cfg.group_to_zet[group]

            _zet = "" if zet == "A" else SUPER_SCRIPTS[zet]
            tex_zet = "" if zet == "A" else "^" + zet

            replaced.append(
                Pair(
                    Alpha(alpha, sign, cfg=cfg),
                    Xi(
                        "∇{}x{}".format(_zet, group),
                        tex="\\nabla{}\\times {}".format(tex_zet, group),
                    ),
                    cfg=cfg,
                )
            )

            for candidate in [c["term"].pair for c in cands]:
                pairs.remove(candidate)

    return replaced, pairs


def replace_grad(pairs, cfg):
    """Grad f = αx[df/dx] + αy[df/dy] + αz[df/dz]"""
    replaced = []
    for zet in _present_zets(pairs, cfg):
        comps = [cfg.zet_comps[zet][k] for k in ["x", "y", "z"]]

        candidates = _filter_on_partials(pairs, comps, cfg)
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

                xi = "".join(SUB_SCRIPTS[x] for x in candidates[0].xi)
                tex_xi = candidates[0].xi
                alpha = alpha_to_group(candidates[0].alpha)

                tex_zet = "" if zet == "A" else "^" + zet
                _zet = "" if zet == "A" else SUPER_SCRIPTS[zet]

                replaced.append(
                    Pair(
                        Alpha(alpha, sign, cfg=cfg),
                        Xi(
                            "∇{}Ξ{}".format(_zet, xi),
                            tex="\\nabla" + tex_zet + "\\Xi_{" + tex_xi + "}",
                        ),
                        cfg=cfg,
                    )
                )

                for candidate in candidates:
                    pairs.remove(candidate.pair)

    return replaced, pairs


def replace_div(pairs, cfg):
    """Div F = dFx/dx + dFy/dy + dFz/dz"""
    replaced = []
    for zet in _present_zets(pairs, cfg):
        comps = [cfg.zet_comps[zet][k] for k in ["x", "y", "z"]]

        # partial from this zet and x-partial with x-3vec component, etc
        candidates = [c for c in _filter_on_partials(pairs, comps, cfg) if Nat(c.xi) == Nat(c.d)]

        if len(candidates) == 3:
            alpha = candidates[0].alpha
            xi = Zet(candidates[0].xi)

            if all(c.sign == 1 for c in candidates):
                sign = 1
            elif all(c.sign == -1 for c in candidates):
                sign = -1
            else:
                continue

            # The 'A' 3Vector calculus operators are the standard ones
            _zet = "" if zet == "A" else SUPER_SCRIPTS[zet]
            tex_zet = "" if zet == "A" else SUPER_SCRIPTS[zet]

            replaced.append(
                Pair(
                    Alpha(alpha, sign, cfg=cfg),
                    Xi("∇{}•{}".format(_zet, xi), tex="\\nabla{}\\cdot {}".format(tex_zet, xi)),
                    cfg=cfg,
                )
            )

            for candidate in candidates:
                pairs.remove(candidate.pair)

    return replaced, pairs


def replace_partials(pairs, cfg):
    """ Partial F = d{comp}F"""

    def key(p):
        return alpha_to_group(p.alpha.index)

    replaced = []

    for zet in _present_zets(pairs, cfg):
        for _, grouped_pairs in groupby(sorted(pairs, key=key), key):
            grouped_pairs = [p for p in grouped_pairs]
            for blade in cfg.zet_comps[zet].values():
                candidates = [
                    p for p in grouped_pairs if p.xi.partials and p.xi.partials[0].index == blade
                ]
                if len(candidates) == 3:
                    ix = candidates[0].xi.val
                    component_zet = Zet(ix)
                    needed = {cfg.zet_comps[component_zet][k] for k in ["x", "y", "z"]}
                    have = {c.xi.val for c in candidates}
                    if have != needed:
                        continue

                    if all([c.xi.sign == 1 for c in candidates]):
                        sign = 1
                    elif all([c.xi.sign == -1 for c in candidates]):
                        sign = -1
                    else:
                        continue

                    alpha = alpha_to_group(candidates[0].alpha.index)

                    _blade = "".join(SUB_SCRIPTS[b] for b in blade)
                    replaced.append(
                        Pair(
                            Alpha(alpha, sign, cfg=cfg),
                            Xi(
                                "∂{}{}".format(_blade, component_zet),
                                tex="\\partial_{}{}".format(blade, component_zet),
                            ),
                            cfg=cfg,
                        )
                    )

                    for candidate in candidates:
                        pairs.remove(candidate)

    return replaced, pairs
