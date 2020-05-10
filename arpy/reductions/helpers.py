from typing import Callable, List, Set, Tuple

from ..algebra.data_types import Term
from ..config import ARConfig
from ..consts import Orientation, Zet

PartialReplacement = Tuple[List[Term], List[Term]]
ReplacementFunc = Callable[[List[Term], Zet], PartialReplacement]


def alpha_to_group(index: str, cfg: ARConfig) -> str:
    E_key = "0i" if cfg._E[0][0] == "0" else "i0"
    groups = {"B": "jk", "T": "0jk", "A": "i", "E": E_key}

    if Orientation.from_index(index) == Orientation.T:
        return index

    return groups[Zet.from_index(index).name]


def by_alpha_group(t: Term) -> str:
    return alpha_to_group(t.index, t.cfg)


def by_xi_group(t: Term) -> str:
    return alpha_to_group(t._components[0].val, t.cfg)


def by_xi(t: Term) -> str:
    return t._components[0].val


def present_zets(terms: List[Term]) -> Set[Zet]:
    return set(
        [
            Zet.from_index(t._component_partials[0].index)
            for t in terms
            if len(t._component_partials) > 0
        ]
    )


def first_partial_str(t: Term) -> str:
    return t._components[0].partials[0]._index


def partial_is(t: Term, partial: str) -> bool:
    p = t._components[0].partials
    if len(p) == 0:
        return False
    return first_partial_str(t) == partial


def partial_in(t: Term, partials: List[str]) -> bool:
    p = t._components[0].partials
    if len(p) == 0:
        return False
    return first_partial_str(t) in partials


def filter_by_partials(terms: List[Term], partials: List[str]) -> List[Term]:
    return sorted([t for t in terms if partial_in(t, partials)], key=by_xi_group,)


def for_all_zets(func: ReplacementFunc, terms: List[Term]) -> PartialReplacement:
    replaced = []
    for zet in [Zet.B, Zet.T, Zet.A, Zet.E]:
        _replaced, terms = func(terms, zet)
        replaced.extend(_replaced)

    return replaced, terms
