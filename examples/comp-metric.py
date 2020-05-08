"""
An attempt at categorising the differences between the +--- and -+++
metrics when used with the 0i ("logical") space-time bivector ordering.

At present, this is NOT complete and should be added to as further cases
and examples are found.
"""
from itertools import groupby

from arpy import ARContext, __version__, commutator, dagger, sign_cayley, sign_distribution

PRINT_ALL = False
ALLOWED = "p 23 31 12 0 023 031 012 123 1 2 3 0123 01 02 03".split()
EXYZ_ALLOWED = "p 0 123 0123 23 023 1 01 31 031 2 02 12 012 3 03".split()


print("Computed using arpy version %s" % __version__)


def compare_configs(allowed=ALLOWED, func=None):
    """
    Take a function and feed it both configs so we can compare the two.
    This is _not_ a normal decorator, it will call the functions immediately
    when the script is run.
    """
    # Hack to add arguments to a decorator
    if func is None:
        return lambda f: compare_configs(allowed=allowed, func=f)

    pmmm = ARContext(allowed, "+---", "into", print_all=PRINT_ALL)
    mppp = ARContext(allowed, "-+++", "by", print_all=PRINT_ALL)

    # Print the name of the function and its docstring
    qname = func.__qualname__
    print("\n.: %s :." % qname)
    print("=" * (len(qname) + 6))
    print(func.__doc__)

    # Iterate over the configs
    for ctx in [pmmm, mppp]:
        print("[%s]" % "".join("+" if m == 1 else "-" for m in ctx.metric))
        func(ctx)
        print("")


@compare_configs()
def squaring_to_neg_alpha_p(ctx):
    """
    One of the primary differences between metrics is in which elements
    square to -αp. This affects which elements support oscillations amongst
    other things.
    """
    pos = []
    neg = []

    with ctx as ar:
        for a in ar.allowed:
            res = ar("a%s a%s" % (a, a)).sign
            if res == 1:
                pos.append("α%s" % a)
            else:
                neg.append("α%s" % a)

        print(" αp:  %s" % ", ".join(pos))
        print("-αp:  %s" % ", ".join(neg))


@compare_configs()
def square_of_a_vector(ctx):
    """
    The square of a 4Vector should be an invariant scalar (== 0 for light).
    """
    with ctx as ar:
        v = ar("{0 1 2 3}")
        v_2 = ar("v v")
        print("v = %s" % v)
        print("\nv² = %s" % v_2)


@compare_configs()
def handedness_under_product(ctx):
    """
    The 3D handedness of each of the 3Vectors when forming products
    between the elements of the 3Vector. This is computed using:
        x y = z    --> Right Handed
        x y = -z   --> Left Handed
    """
    with ctx as ar:
        for vec in "BTAE":
            comps = ar.cfg.zet_comps[vec]
            sign = ar("a%s a%s" % (comps["x"], comps["y"])).sign
            handedness = "RH" if sign == 1 else "LH"
            print("%s:  %s" % (vec, handedness))


@compare_configs()
def handedness_under_rotation(ctx):
    """
    The 3D handedness of each of the 3Vectors when forming products
    between the elements of the 3Vector and a space-space bivector.
        α23 y = z     --> Right Handed
        α23 y = -z    --> Left Handed
    """
    with ctx as ar:
        for vec in "BTAE":
            comps = ar.cfg.zet_comps[vec]
            sign = ar("a%s a%s" % ("23", comps["y"])).sign
            handedness = "RH" if sign == 1 else "LH"
            print("%s:  %s" % (vec, handedness))


@compare_configs()
def full_handedness(ctx):
    """
    Look at the handedness of all 3Vector interactions.
        x y = z   --> Right Handed (□)
        x y = -z  --> Left Handed  (■)
    """
    with ctx as ar:
        print("   B T A E")
        for vec1 in "BTAE":
            print(vec1, end="  ")
            for vec2 in "BTAE":
                comps = ar.cfg.zet_comps
                v1, v2 = comps[vec1], comps[vec2]
                sign = ar("a%s a%s" % (v1["x"], v2["y"])).sign
                handedness = "□" if sign == 1 else "■"
                print(handedness, end=" ")
            print("")


@compare_configs()
def zet_sign_distribution(ctx):
    """
    The full sign Cayley table above can be broken down into five groups that
    always have a consistent sign:

      ∂e |■ □ □ □|  ∂Ξ |□ ■ ■ ■|   ∇ |□ □ □ □|  ∇• |□ □ □ □|  ∇x |□ □ □ □|
         |□ □ □ □|     |□ □ □ □|     |■ □ □ □|     |□ ■ □ □|     |□ □ ■ ■|
         |□ □ □ □|     |□ □ □ □|     |■ □ □ □|     |□ □ ■ □|     |□ ■ □ ■|
         |□ □ □ □|     |□ □ □ □|     |■ □ □ □|     |□ □ □ ■|     |□ ■ ■ □|

    (Note that, in the case of the curl (∇x), □ represents a Right Handed
     set and ■ a Left Handed set.)
    In the following diagrams, the rows and columns are the 3Vector(+extra)
    sets: `B T A E`.

    """
    with ctx as ar:
        sign_distribution(cfg=ar.cfg)


@compare_configs(allowed=EXYZ_ALLOWED)
def exyz_sign_distribution(ctx):
    """
    When the grouping / ordering of the components is changed to group
    by exyz first and then by 3Vector, the sign distribution may shed
    some light on the sign distribution of the algebra.

    This time, the rows/columns are `e x y z` and the meaning of the diagrams
    is no longer the traditional vector calculus indicated by the notation.
    (Strictly, I should write another function that uses different notation but
     I'm not sure what to use yet / if this is at all useful!)
    """
    with ctx as ar:
        sign_distribution(cfg=ar.cfg)


@compare_configs()
def maxwell_is_maintained(ctx):
    """
    Maxwell's Equations should have the correct or inverted signs. The sign
    of vector, trivector, pivot, hedgehog and quedgehog terms are unknown
    as of now.
    """
    # These are the signs of the relevant components, computed using division
    # into and i0 using the following code:
    # ------------------------------------------------------------------------
    # Dg = Dmu(G)
    # MAXWELL = [
    #     [
    #         k.sign
    #         for k in sorted(Dg[blade], key=lambda k: k._component_partials)
    #         if k._components[0].val not in ["p", "0123"]
    #     ]
    #     for blade in "123 1 2 3 0 023 031 012".split()
    # ]
    # NEG_MAXWELL = [
    #     [
    #         k.sign * -1
    #         for k in sorted(Dg[blade], key=lambda k: k._component_partials)
    #         if k._components[0].val not in ["p", "0123"]
    #     ]
    #     for blade in "123 1 2 3 0 023 031 012".split()
    # ]
    MAXWELL = [
        [-1, -1, -1],
        [1, -1, 1],
        [-1, 1, 1],
        [1, -1, 1],
        [-1, -1, -1],
        [1, -1, 1],
        [1, 1, -1],
        [1, -1, 1],
    ]
    NEG_MAXWELL = [
        [1, 1, 1],
        [-1, 1, -1],
        [1, -1, -1],
        [-1, 1, -1],
        [1, 1, 1],
        [-1, 1, -1],
        [-1, -1, 1],
        [-1, 1, -1],
    ]

    with ctx as ar:
        XiG = "{%s}" % " ".join(ar.allowed)
        Dg = ar("<0 1 2 3> %s" % XiG)

        # Ignoring the sign of pivot terms as they are not governed by maxwell
        maxwell_signs = [
            [
                k.sign
                for k in sorted(Dg[blade], key=lambda k: k._component_partials)
                if k._components[0].val not in ["p", ar.cfg.allowed[-4]]
            ]
            for blade in ["0", "1", "2", "3", "123", "023", "031", "012"]
        ]

        if maxwell_signs == MAXWELL:
            print("Gives conventional Maxwell")
        elif maxwell_signs == NEG_MAXWELL:
            print("Gives inverted Maxwell")
        else:
            print("Does not give Maxwell")

        print(Dg)


@compare_configs()
def FF_dagger_force_equation(ctx):
    """
    F F_dagger style products should exhibit cancelation and lead to force
    equations.
    """
    with ctx as ar:
        F = ar("{ 23 31 12 01 02 03 }")
        res = ar("F F!")
        print("F = %s\n" % F)
        print("F_dagger = %s\n" % dagger(F, cfg=ar.cfg))
        print("F F_dagger = %s" % res)


@compare_configs()
def GG_dagger_force_equation(ctx):
    """
    G G_dagger style products should exhibit cancelation and lead to force
    equations.
    """
    with ctx as ar:
        G = ar("{%s}" % " ".join(ar.cfg.allowed))
        res = ar("G G!")
        print("G = %s\n" % G)
        print("G_dagger = %s\n" % dagger(G, cfg=ar.cfg))
        print("G G_dagger = %s" % res)


@compare_configs()
def group_commutator_full(ctx):
    """
    Compute the group commutator [a, b] = (a b a^-1 b^-1) for each combination
    of a and b in the algebra.

    [a, b] =  αp  --> □
    [a, b] = -αp  --> ■
    """
    with ctx as ar:
        sign_cayley(commutator, cfg=ar.cfg)


@compare_configs()
def group_commutator_distribution(ctx):
    """
    This is the same idea as above, combined with the sign distribution format
    from before that looks at vector calculus groupings.

    [a, b] =  αp  --> □
    [a, b] = -αp  --> ■
    """
    with ctx as ar:
        sign_distribution(commutator, cfg=ar.cfg)


@compare_configs()
def FF_dagger_signs(ctx):
    """
    Sign distributions for FF! force equations
    """

    with ctx as ar:
        for alpha, terms in groupby(ar("F F!"), lambda t: t._alpha):
            signs = " ".join(["□" if t._sign == -1 else "■" for t in terms])
            print(f"{str(alpha).ljust(5)} {signs}")
