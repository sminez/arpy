import pytest
from arpy import Alpha, Xi, Pair
from arpy.algebra.differential import component_partial
from arpy.algebra.del_grouping import replace_div, replace_grad, \
    replace_partials


def test_new_component_partial():
    '''
    component_partial returns a new object rather than modifying
    the original and passing it back
    '''
    original = Pair('012')
    wrt = Alpha('2')
    differentiated = component_partial(original, wrt, 'by', (-1, 1, 1, 1))
    assert differentiated is not original


def test_replace_curl():
    '''
    A valid set of curl terms gets replaced correctly
    '''
    pass


@pytest.mark.parametrize('sign', [1, -1])
def test_replace_grad(sign):
    '''
    A valid set of grad terms gets replaced correctly
    '''
    grad_like = [
        Pair(Alpha(i), Xi('p', partials=[Alpha(i)], sign=sign))
        for i in ['1', '2', '3']
    ]
    replaced, left_over = replace_grad(grad_like)
    assert left_over == []
    assert replaced == [Pair(Alpha('i', sign), Xi('∇Ξp'))]


@pytest.mark.parametrize('sign', [1, -1])
def test_replace_div(sign):
    '''
    A valid set of div terms gets replaced correctly
    '''
    div_like = [
        Pair(Alpha('p'), Xi(i, partials=[Alpha(i)], sign=sign))
        for i in ['1', '2', '3']
    ]
    replaced, left_over = replace_div(div_like)
    assert left_over == []
    assert replaced == [Pair(Alpha('p', sign), Xi('∇•A'))]


@pytest.mark.parametrize('sign', [1, -1])
def test_replace_partial(sign):
    '''
    A valid set of partial terms gets replaced correctly
    '''
    partial_like = [
        Pair(Alpha(i), Xi(i, partials=[Alpha('p')], sign=sign))
        for i in ['1', '2', '3']
    ]
    replaced, left_over = replace_partials(partial_like)
    assert left_over == []
    assert replaced == [Pair(Alpha('i', sign), Xi('∂ₚA'))]
