'''
Experimental Symbolic manipulation of algebraic expressions.

To create new symbolic variables you need to do the following:
>>> x, y, z = symbols('x, y, z')
'''
from arpy import Alpha, config, full

# Bring in some standard sympy mathematical functions and operators
from sympy import symbols, Add, Mul, sin, sinc, cos, cosh, tan, tanh, log, exp


# Vars for sympy work
ALPHA_SYMS = symbols(' '.join('a{}'.format(a) for a in config.allowed))
ALPHA_MAP = {s: Alpha(str(s)[1:]) for s in ALPHA_SYMS}
SYM_MAP = dict(reversed(item) for item in ALPHA_MAP.items())


# Create our symbolic alphas for use in scripts
# NOTE: In case we are making anything else later, we MUST set commutative=0
#       for ordering to me maintained in sympy.args calls
_a_str = ', '.join('a{}'.format(a) for a in config.allowed)
exec('{a} = symbols("{a}", commutative=0)'.format(a=_a_str))


def sym_diff(*alphas):
    '''
    Symbolically differentiate a sympy function with respect to the
    given alpha values.
    '''
    def _diff(func):
        return arpify(sum([func.diff(symbols('a'+a.ix)) for a in alphas]))

    return _diff


def arpify(expr):
    '''
    Perform recursive arpy alpha reductions on a sympy expression.
    '''
    # Map over the sub expressions of a sum
    if isinstance(expr, Add):
        return sum([arpify(arg) for arg in expr.args])

    pending = None
    others = []

    for arg in expr.args:
        if arg in ALPHA_SYMS:
            if pending is None:
                pending = arg
            else:
                others = [arg]
                break

        # Compute for each branch of a sum
        elif isinstance(arg, Add):
            if all([a in ALPHA_SYMS for a in arg.args]):
                if pending is None:
                    pending = arg.args
                else:
                    others = arg.args
                    break

        # Simplify what we can from a product
        elif isinstance(arg, Mul):
            pass

    if others == []:
        return expr

    for other in others:
        res = full(ALPHA_MAP[pending], ALPHA_MAP[other])
        if res.sign == -1:
            res.sign *= -1
            res = SYM_MAP[res]
            expr = expr.expand().subs(pending * other, -1 * res).simplify()
        else:
            res = SYM_MAP[res]
            expr = expr.expand().subs(pending * other, res).simplify()

    return arpify(expr)
