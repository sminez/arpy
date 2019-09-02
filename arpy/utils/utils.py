SUPER_SCRIPTS = {"B": "ᴮ", "A": "ᴬ", "T": "ᵀ", "E": "ᴱ"}
SUB_SCRIPTS = {"0": "₀", "1": "₁", "2": "₂", "3": "₃", "p": "ₚ", "i": "ᵢ", "j": "ⱼ", "k": "ₖ"}
ZET_MAP = {
    # B: p 23 31 12
    frozenset("p"): {"direction": "t", "zet": "B"},
    frozenset("23"): {"direction": "x", "zet": "B"},
    frozenset("31"): {"direction": "y", "zet": "B"},
    frozenset("12"): {"direction": "z", "zet": "B"},
    # T: 0 023 031 012
    frozenset("0"): {"direction": "t", "zet": "T"},
    frozenset("023"): {"direction": "x", "zet": "T"},
    frozenset("031"): {"direction": "y", "zet": "T"},
    frozenset("012"): {"direction": "z", "zet": "T"},
    # A: 123 1 2 3
    frozenset("123"): {"direction": "t", "zet": "A"},
    frozenset("1"): {"direction": "x", "zet": "A"},
    frozenset("2"): {"direction": "y", "zet": "A"},
    frozenset("3"): {"direction": "z", "zet": "A"},
    # E: 0123 01 02 03
    frozenset("0123"): {"direction": "t", "zet": "E"},
    frozenset("01"): {"direction": "x", "zet": "E"},
    frozenset("02"): {"direction": "y", "zet": "E"},
    frozenset("03"): {"direction": "z", "zet": "E"},
}


def Tex(obj):
    """
    Convert the string representation of an object to TeX and print.
    """
    print(obj.__tex__())


def Zet(alpha):
    """Return the Zet of a given alpha value or term"""
    if isinstance(alpha, str):
        ix = alpha
    else:
        ix = alpha.index

    return ZET_MAP[ix]["zet"]


def Dir(alpha):
    """Return the Space-Time 'direction' of a given alpha or term"""
    if isinstance(alpha, str):
        ix = alpha
    else:
        ix = alpha.index

    return ZET_MAP[ix]["direction"]


def Nat(alpha):
    """Return the Nature of a given Alpha."""
    # Element sets for each e,x,y,z nature
    nat_map = {
        frozenset("p"): "e",
        frozenset("123"): "e",
        frozenset("0"): "e",
        frozenset("0123"): "e",
        frozenset("1"): "x",
        frozenset("23"): "x",
        frozenset("023"): "x",
        frozenset("01"): "x",
        frozenset("2"): "y",
        frozenset("31"): "y",
        frozenset("031"): "y",
        frozenset("02"): "y",
        frozenset("3"): "z",
        frozenset("12"): "z",
        frozenset("012"): "z",
        frozenset("03"): "z",
    }

    # Allow for raw string indices to be passed
    if isinstance(alpha, str):
        ix = alpha
    else:
        ix = alpha.index

    return nat_map.get(frozenset(ix))


def reorder_allowed(allowed, order="pBtThAqE"):
    """
    Shuffle the ordering of allowed, keeping 3-Vectors together.
    NOTE: This assumes that the input is in pBtThAqE order to start.
    """
    p = ["p"]
    t = ["0"]
    h = [a for a in allowed if len(a) == 3 and "0" not in a]
    q = [a for a in allowed if len(a) == 4]
    B = [a for a in allowed if len(a) == 2 and "0" not in a]
    T = [a for a in allowed if len(a) == 3 and "0" in a]
    A = [a for a in allowed if len(a) == 1 and a not in ["p", "0"]]
    E = [a for a in allowed if len(a) == 2 and "0" in a]

    groups = {"p": p, "t": t, "h": h, "q": q, "B": B, "T": T, "A": A, "E": E}
    new = []
    for group in order:
        new += groups[group]
    return new
