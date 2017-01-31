'''
Lexing and Parsing of a more mathematical syntax for performing calculations
with the arpy Absolute Relativity library.

SLY is rather magical - docs are here:
    http://sly.readthedocs.io/en/latest/
    https://github.com/dabeaz/sly
'''
import sly
from sys import _getframe
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

    def __init__(self, _globals):
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
        return wedge(p.expr0, p.expr1)

    @_('expr "." expr')
    def expr(self, p):
        return dot(p.expr0, p.expr1)

    @_('expr "*" expr')
    def expr(self, p):
        return full(p.expr0, p.expr1)

    @_('expr "/" expr')
    def expr(self, p):
        return div_by(p.expr0, p.expr1)

    @_('expr "\\" expr')
    def expr(self, p):
        return div_into(p.expr0, p.expr1)

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


def ar(user_input):
    '''
    Parse user input under a nicer API.
    Call ar_help() for details on syntax.
    '''
    # NOTE:: The following is a horrible hack that allows you to
    #        inject local variables into the parser.
    frame = _getframe(1)
    lexer = ArpyLexer(frame.f_locals)
    parser = ArpyParser()
    try:
        return parser.parse(lexer.tokenize(user_input))
    except sly.lex.LexError as e:
        raise SyntaxWarning(e)
