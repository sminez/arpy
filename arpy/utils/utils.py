from collections import Counter
from typing import Any, Union
from ..consts import Zet, Orientation

SUPER_SCRIPTS = {"B": "ᴮ", "A": "ᴬ", "T": "ᵀ", "E": "ᴱ"}
SUB_SCRIPTS = {"0": "₀", "1": "₁", "2": "₂", "3": "₃", "p": "ₚ", "i": "ᵢ", "j": "ⱼ", "k": "ₖ"}


def tex(obj: Any) -> str:
    """
    Convert the string representation of an object to TeX and print.
    """
    print(obj.__tex__())


def zet(alpha: Union[str, "Alpha"]) -> Zet:
    ix = alpha.index if isinstance(alpha, "Alpha") else alpha
    if not isinstance(ix, str):
        raise ValueError("argument must be a string or Alpha")

    return Zet.from_index(ix)


def orientation(alpha: Union[str, "Alpha"]) -> Orientation:
    ix = alpha.index if isinstance(alpha, "Alpha") else alpha
    if not isinstance(ix, str):
        raise ValueError("argument must be a string or Alpha")

    return Orientation.from_index(ix)


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


def power_notation(lst):
    """
    Given a list of elements (typically representing values of some sort,
    relating to a string representation of the result of a computation)
    express repeated elements in a^b power notation.
    """
    result = []

    for item, power in Counter(lst).items():
        if power == 1:
            result.append(str(item))
        else:
            result.append(f"{item}^{power}")

    return result
