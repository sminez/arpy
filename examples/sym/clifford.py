'''
A sympy based version of 4th order Clifford Algebra
'''
from sympy import sympify
from sympy import Add, Mul
from sympy.core.expr import Expr

from arpy import find_prod, Alpha


def _generic_mul(a, b):
    '''
    Run the arpy full product algorithm
    '''
    # Handle special cases first
    if isinstance(a, SymAlpha) and isinstance(b, SymMultiVector):
        mvec = {alpha[1:]: getattr(b, alpha) for alpha in b._alphas}

        prods = {}
        for ix, vals in mvec.items():
            res = find_prod(Alpha(a.ix), Alpha(ix))
            prods['a' + res.index] = vals * res.sign

        return SymMultiVector(**prods)

    elif isinstance(a, SymMultiVector) and isinstance(b, Alpha):
        mvec = {alpha[1:]: getattr(a, alpha) for alpha in a._alphas}

        prods = {}
        for ix, vals in mvec.items():
            res = find_prod(Alpha(ix), Alpha(b.ix))
            prods['a' + res.index] = vals * res.sign

        return SymMultiVector(**prods)

    elif all([isinstance(m, SymAlpha) for m in [a, b]]):
        res = find_prod(Alpha(a.ix), Alpha(b.ix))
        ix = res.index
        if res.sign == -1:
            ix = '-' + ix

        return Alpha(ix)

    elif all([isinstance(m, SymMultiVector) for m in [a, b]]):
        # Compute the full product
        prods = {}
        for i, ival in a._alphas.items():
            for j, jval in b._alphas.items():
                res = find_prod(Alpha(i[1:]), Alpha(j[1:]))
                key = 'a' + res.index
                current = prods.get(key, 0)
                current += ival * jval * res.sign
                prods[key] = current

        return SymMultiVector(**prods)

    return Mul(a, b)


class SymAlpha(Expr):
    '''
    A symbolic alpha value
    '''
    _op_priority = 11.0
    is_commutative = False
    is_number = False

    def __new__(cls, ix):
        obj = Expr.__new__(cls, sympify('a' + ix))
        obj._ix = ix
        return obj

    def __repr__(self):
        s = str(self._ix)
        if s[0] == '-':
            return '-α' + s[1:]
        else:
            return 'α' + s

    @property
    def ix(self):
        return self._ix

    # def __add__(self, other):
    #     return self.add(other)

    # def __radd__(self, other):
    #     return self.add(other)

    # def __sub__(self, other):
    #     return self.add(other*-1)

    def __mul__(self, other):
        return _generic_mul(self, other)

    def __rmul__(self, other):
        return _generic_mul(other, self)

    def __neg__(self):
        if self.ix[0] == '-':
            ix = self.ix[1:]
        else:
            ix = '-' + self.ix

        return SymAlpha(ix)


class SymMultiVector(Expr):
    '''
    A 4th order Clifford Algebra Multivector.

    Currently hard coded to use the default arpy config:
        Allowed Alphas: αp, α23, α31, α12, α0, α023, α031, α012,
                        α123, α1, α2, α3, α0123, α01, α02, α03
        Metric:         +---
        Division type:  into

    Ultimately, this needs to be like namedtuple and allow for custom
    properties at init with verification. Or alternatively, this
    class should be redefined when the config is changed (if that is
    possible).
    '''
    _op_priority = 11.0
    is_commutative = False

    def __new__(
            cls,
            ap=0, a0=0, a1=0, a2=0, a3=0, a23=0, a31=0, a12=0, a01=0,
            a02=0, a03=0, a023=0, a031=0, a012=0, a123=0, a0123=0,
            real_field=True
    ):
        alphas = {
            "ap": ap, "a0": a0, "a1": a1, "a2": a2, "a3": a3,
            "a23": a23, "a31": a31, "a12": a12, "a01": a01,
            "a02": a02, "a03": a03, "a023": a023, "a031": a031,
            "a012": a012, "a123": a123, "a0123": a0123,
        }
        alphas = {k: sympify(v) for k, v in alphas.items()}

        if any(a.is_commutative is False for a in alphas.values()):
            raise ValueError("arguments have to be commutative")
        else:
            obj = Expr.__new__(cls, *alphas.values())

            for k, v in alphas.items():
                setattr(obj, k, v)

            obj.real_field = real_field
            obj._alphas = alphas

            return obj

    def __repr__(self):
        comps = [
            '    α{}({})'.format(a[1:].ljust(5), getattr(self, a))
            for a in self._alphas
            if getattr(self, a) != 0
        ]
        return '{\n' + '\n'.join(comps) + '\n}'

    def __add__(self, other):
        return self.add(other)

    def __radd__(self, other):
        return self.add(other)

    def __sub__(self, other):
        return self.add(other*-1)

    def __mul__(self, other):
        return _generic_mul(self, other)

    def __rmul__(self, other):
        return _generic_mul(other, self)

    def __pow__(self, p):
        return self.pow(p)

    def __neg__(self):
        return SymMultiVector(
            -self.ap, -self.a0, -self.a1, -self.a2, -self.a3,
            -self.a23, -self.a31, -self.a12, -self.a01, -self.a02,
            -self.a03, -self.a023, -self.a031, -self.a012,
            -self.a123, -self.a0123,
        )

    def _eval_Integral(self, *args):
        return self.integrate(*args)

    def _eval_diff(self, *symbols, **kwargs):
        return self.diff(*symbols)

    def add(self, other):
        '''Add MultiVectors.'''
        m1 = self
        m2 = sympify(other)

        # If q2 is a number or a sympy expression instead of a quaternion
        if not isinstance(m2, SymMultiVector):
            return Add(m1, m2)

        args = [getattr(self, a) + getattr(other, a) for a in self._alphas]
        return SymMultiVector(*args)

    def mul(self, other):
        '''Multiply MultiVectors.'''
        return _generic_mul(self, other)
