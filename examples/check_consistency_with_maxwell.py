"""
This program uses the arpy library to verify whether or not AR candidate algebras are
consistent with Maxwell's laws of Electromagnetism. Algebras are rejected if they do
not give a superset of the Maxwell Equations as the result of Dmu(G) and they are
accepted if they give Maxwell in either the conventional sign or entirely negated.

There are some Algebras that can be _made_ to work by mapping the αjk bivectors as the
negative of the magnetic field as opposed to mapping them to the magnetic field
directly. To include checking such candidates, pass the '--allow-neg-B' flag on the
command line.
"""
from argparse import ArgumentParser
from collections import namedtuple
from itertools import permutations

from arpy import ARContext
from arpy.utils.utils import SUB_SCRIPTS

Candidate = namedtuple("candidate", "allowed div metric for_maxwell")
Result = namedtuple("result", "sign negate_B metric division allowed pivot_terms")


def allowed_repr(allowed):
    def as_group(ix):
        remap = {"0": "0", "1": "i", "2": "j", "3": "k"}
        return "".join(remap[c] for c in allowed[ix])

    B = as_group(1)
    T = as_group(5)
    E = as_group(13)
    h = as_group(8)
    q = as_group(12)

    return f"p 0 {h} {q} | {B} {T} i {E}"


def all_candidates(include_redundant):
    def configurations(seed, n):
        return [
            ["".join(element[i] for i in p) for element in seed] for p in permutations(range(n))
        ]

    Bs = [["23", "31", "12"], ["32", "13", "21"]]
    Es = [["01", "02", "03"], ["10", "20", "30"]]

    if include_redundant:
        Ts = configurations(["023", "031", "012"], n=3)
        hs = configurations(["123"], n=3)
        qs = configurations(["0123"], n=4)
        metrics = ["----", "+---", "--++", "-+++", "++++"]
    else:
        Ts = [["023", "031", "012"], ["032", "013", "021"]]
        hs = [["123"], ["321"]]
        qs = [["0123"], ["1230"]]
        metrics = ["+---", "-+++"]

    return [
        Candidate(
            allowed=["p"] + B + ["0"] + T + h + ["1", "2", "3"] + q + E,
            div=division,
            metric=metric,
            for_maxwell=["0", "1", "2", "3"] + h + T,
        )
        for B in Bs
        for E in Es
        for T in Ts
        for h in hs
        for q in qs
        for metric in metrics
        for division in ["by", "into"]
    ]


def compact_term_repr(t):
    sign = "+" if t.sign == 1 else ""  # be explicit about +/-
    comps = ".".join(map(str, sorted(t._components)))
    partial_strs = [
        "∂{}".format("".join(SUB_SCRIPTS[i] for i in p._index))
        for p in sorted(t._component_partials)
    ]
    partials = "".join(partial_strs)

    return f"{sign}{t.alpha}{partials}{comps}"


def check_candidate(candidate, allow_neg_B):
    def _check(Dg, negated_B=False):
        maxwell_signs = [
            "".join(
                "+" if t.sign == 1 else "-"
                for t in sorted(Dg[blade], key=lambda t: t._component_partials)
                if t._components[0].val not in ["p", ar.allowed[-4]]
            )
            for blade in candidate.for_maxwell
        ]

        pivot_terms = [
            [
                compact_term_repr(t)
                for t in Dg[blade]
                if t._components[0].val in ["p", ar.allowed[-4]]
            ]
            for blade in candidate.for_maxwell
        ]

        if maxwell_signs not in [maxwell, negated_maxwell]:
            print(".", end="", flush=True)
            return None

        if maxwell_signs == maxwell:
            print("+", end="", flush=True)
            sign = "+"
        else:
            print("-", end="", flush=True)
            sign = "-"

        return Result(
            sign=sign,
            negate_B=negated_B,
            metric=candidate.metric,
            division=candidate.div,
            allowed=candidate.allowed,
            pivot_terms=" ".join(sum(pivot_terms, [])),
        )

    def _is_B(s):
        return len(s) == 2 and "0" not in s

    maxwell = ["---", "+-+", "-++", "+-+", "---", "+-+", "++-", "+-+"]
    negated_maxwell = ["+++", "-+-", "+--", "-+-", "+++", "-+-", "--+", "-+-"]

    context = ARContext(allowed=candidate.allowed, div=candidate.div, metric=candidate.metric,)

    with context as ar:
        XiG = "{%s}" % " ".join(ar.allowed)
        Dg = ar("<0 1 2 3> %s" % XiG)

        res = _check(Dg)
        if res is not None:
            return res

        if allow_neg_B:
            XiG = "{%s}" % " ".join(f"-{a}" if _is_B(a) else a for a in ar.allowed)
            Dg = ar("<0 1 2 3> %s" % XiG)
            return _check(Dg, negated_B=True)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "--show-pivot-signs",
        action="store_true",
        help="show the sign of the additional terms beyond Maxwell",
    )
    parser.add_argument(
        "--include-redundant",
        action="store_true",
        help="construct all possible algebras including those that are isomorphisms",
    )
    parser.add_argument(
        "--allow-neg-B",
        action="store_true",
        help="also check for candidates that produce Maxwell with negated jk bivectors",
    )
    args = parser.parse_args()

    candidates = all_candidates(args.include_redundant)

    print(f"Checking {len(candidates)} candidate Algebras for consistency with Maxwell")
    if args.allow_neg_B:
        print(f"Allowing negated jk bivectors: up to {len(candidates) * 2} cases will be checked")

    results = [check_candidate(c, args.allow_neg_B) for c in candidates]

    hits = sorted([r for r in results if r is not None])
    print(f"\n\nFound {len(hits)} candidate Algebras that support Maxwell")

    for h in hits:
        p = f"  [{h.pivot_terms}]" if args.show_pivot_signs else ""
        if args.allow_neg_B:
            neg_B = " -B " if h.negate_B else "  B "
        else:
            neg_B = ""

        print(f"[{h.sign}]{neg_B}{h.metric} {h.division.ljust(4)} {allowed_repr(h.allowed)} {p}")
