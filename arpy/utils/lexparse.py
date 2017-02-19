'''
Lexing and Parsing of a more mathematical syntax for performing calculations
with the arpy Absolute Relativity library.

SLY is rather magical - docs are here:
    http://sly.readthedocs.io/en/latest/
    https://github.com/dabeaz/sly
'''
import sly
from sys import _getframe
from ..algebra.config import METRIC, DIVISION_TYPE
from ..algebra.ar_types import Alpha, Pair  # MultiVector
from ..algebra.operations import wedge, dot, full, div_by, div_into, project


class ArpyLexer(sly.Lexer):
    tokens = {'ALPHA', 'PAIR', 'INDEX', 'VAR'}
    ignore = ' \t'
    literals = {'+', '^', '.', '*', '/', '\\', '(', ')', '<', '>'}

    ALPHA = r'-a[0123]{1,4}|a[0123]{1,4}|-?ap'
    PAIR = r'-p[0123]{1,4}|p[0123]{1,4}'
    VAR = r'[a-zA-Z_][a-zA-Z_0-9]*'
    INDEX = r'[01234]'

    def __init__(self, _globals=None):
        self._globals = _globals

    def ALPHA(self, token):
        if token.value.startswith('-'):
            token.value = Alpha(token.value[2:], -1)
        else:
            token.value = Alpha(token.value[1:])
        return token

    def PAIR(self, token):
        if token.value.startswith('-'):
            token.value = Pair(Alpha(token.value[2:], -1))
        else:
            token.value = Pair(Alpha(token.value[1:]))
        return token

    def INDEX(self, token):
        token.value = int(token.value)
        return token

    def VAR(self, token):
        token.value = eval(token.value, self._globals)
        return token


class ArpyParser(sly.Parser):
    tokens = ArpyLexer.tokens
    metric = METRIC

    precedence = (
        ('left', '*', '^'),
        ('left', '.', '/'),
        ('left', '\\', 'ALPHA',),
    )

    @_('expr')
    def result(self, p):
        return p.expr

    @_('"<" expr ">" INDEX')
    def expr(self, p):
        return project(p.expr, p.INDEX)

    @_('expr "^" expr')
    def expr(self, p):
        return wedge(p.expr0, p.expr1, self.metric)

    @_('expr "." expr')
    def expr(self, p):
        return dot(p.expr0, p.expr1, self.metric)

    @_('expr "*" expr')
    def expr(self, p):
        return full(p.expr0, p.expr1, self.metric)

    @_('expr "/" expr')
    def expr(self, p):
        return div_by(p.expr0, p.expr1, self.metric)

    @_('expr "\\" expr')
    def expr(self, p):
        return div_into(p.expr0, p.expr1, self.metric)

    @_('ALPHA')
    def expr(self, p):
        return p.ALPHA

    @_('PAIR')
    def expr(self, p):
        return p.PAIR

    @_('VAR')
    def expr(self, p):
        return p.VAR

    @_('"(" expr ")"')
    def expr(self, p):
        return p.expr

##############################################################################


class ARContext:
    '''
    User interface class for working with the library.
    Create an instance and then use by calling ar as a function.
    i.e.
    >>> ar = ARContext()
    >>> ar("a12 ^ a23")
    >>> α31
    '''
    def __init__(self, metric=METRIC, division=DIVISION_TYPE):
        self._metric = METRIC
        self._division = DIVISION_TYPE
        self._lexer = ArpyLexer()
        self._parser = ArpyParser()

    @property
    def metric(self):
        return self._metric

    @metric.setter
    def metric(self, sign_str):
        signs = sign_str.split()
        if not all(sign in ["+", "-"] for sign in signs):
            raise TypeError("metric must be comprised of +/- only")
        if len(signs) != 4:
            raise ValueError(
                "metric should be a 4 element space deliminated string.\n"
                "i.e. 'ar.metric = \"+ - - -\"'"
            )
        self._metric = tuple(1 if s == "+" else -1 for s in signs)
        self._parser.metric = self._metric
        print("metric has been set to '{}'".format(sign_str))

    @property
    def division(self):
        return self._division

    @division.setter
    def division(self, div_type):
        if div_type not in ["by", "into"]:
            raise ValueError("Division type must be either 'by' or 'into'")
        self._division = div_type

    def __call__(self, user_input):
        # NOTE:: The following is a horrible hack that allows you to
        #        inject local variables into the parser.
        stack_frame = _getframe(1)
        self._lexer._globals = stack_frame.f_locals
        try:
            return self._parser.parse(self._lexer.tokenize(user_input))
        except sly.lex.LexError as e:
            raise SyntaxWarning(e)
