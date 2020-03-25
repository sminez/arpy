import pytest

from .. import Alpha, Term, Xi, config
from ..algebra.differential import term_partial
from ..reductions.del_grouping import replace_curl, replace_div, replace_grad, replace_partials


def test_new_component_partial():
    """
    component_partial returns a new object rather than modifying
    the original and passing it back
    """
    original = Term("012")
    wrt = Alpha("2")
    differentiated = term_partial(original, wrt, config, "by")
    assert differentiated is not original


# @pytest.mark.parametrize("sign", [1, -1])
# def test_replace_curl(sign):
#     """
#     A valid set of curl terms gets replaced correctly
#     """
#     curl_like = [
#         Term(Alpha("1"), [Xi("31", partials=[Alpha("3")])], sign=-sign),
#         Term(Alpha("1"), [Xi("12", partials=[Alpha("2")])], sign=sign),
#         Term(Alpha("2"), [Xi("12", partials=[Alpha("1")])], sign=-sign),
#         Term(Alpha("2"), [Xi("23", partials=[Alpha("3")])], sign=sign),
#         Term(Alpha("3"), [Xi("23", partials=[Alpha("2")])], sign=-sign),
#         Term(Alpha("3"), [Xi("31", partials=[Alpha("1")])], sign=sign),
#     ]
#     replaced, left_over = replace_curl(curl_like, config)
#     assert left_over == []
#     assert replaced == [Term(Alpha("i", sign), Xi("∇xB"))]


# @pytest.mark.parametrize("sign", [1, -1])
# def test_replace_grad(sign):
#     """
#     A valid set of grad terms gets replaced correctly
#     """
#     grad_like = [
#         Term(Alpha("1"), [Xi("p", partials=[Alpha("1")])], sign=sign),
#         Term(Alpha("2"), [Xi("p", partials=[Alpha("2")])], sign=sign),
#         Term(Alpha("3"), [Xi("p", partials=[Alpha("3")])], sign=sign),
#     ]
#     replaced, left_over = replace_grad(grad_like, config)
#     assert left_over == []
#     assert replaced == [Term(Alpha("i", sign), [Xi("∇Ξₚ")])]


# @pytest.mark.parametrize("sign", [1, -1])
# def test_replace_div(sign):
#     """
#     A valid set of div terms gets replaced correctly
#     """
#     div_like = [
#         Term(Alpha("0"), [Xi("01", partials=[Alpha("1")])], sign=sign),
#         Term(Alpha("0"), [Xi("02", partials=[Alpha("2")])], sign=sign),
#         Term(Alpha("0"), [Xi("03", partials=[Alpha("3")])], sign=sign),
#     ]
#     replaced, left_over = replace_div(div_like, config)
#     assert left_over == []
#     assert replaced == [Term(Alpha("0", sign), [Xi("∇•E")])]


# @pytest.mark.parametrize("sign", [1, -1])
# def test_replace_partial(sign):
#     """
#     A valid set of partial terms gets replaced correctly
#     """
#     partial_like = [
#         Term(Alpha("23"), [Xi("23", partials=[Alpha("p")])], sign=sign),
#         Term(Alpha("31"), [Xi("31", partials=[Alpha("p")])], sign=sign),
#         Term(Alpha("12"), [Xi("12", partials=[Alpha("p")])], sign=sign),
#     ]
#     replaced, left_over = replace_partials(partial_like, config)
#     assert left_over == []
#     assert replaced == [Term(Alpha("jk", sign), [Xi("∂ₚB")])]
