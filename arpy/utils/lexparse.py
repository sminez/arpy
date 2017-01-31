'''
Lexing and Parsing of a more mathematical syntax for performing calculations
with the arpy Absolute Relativity library.

SLY is rather magical - docs are here:
    http://sly.readthedocs.io/en/latest/
    https://github.com/dabeaz/sly
'''
from sly import Lexer, Parser
from ..algebra.ar_types import Alpha, Pair, Xi_vecs
from ..algebra.operations import wedge, div_by, div_into
from ..algebra.differential import Dmu


class ArpyLexer(Lexer):
    tokens = {'ALPHA', 'PAIR'}  # TODO:: add multivectors
    ignore = ' \t'
    literals = set(r'+ ^ . * / \ ( ) < > { }'.split())

    @_(r'-a[0123]{1,4}', r'a[0123]{1,4}', r'-?ap')
    def ALPHA(self, t):
        try:
            if t.value.startswith('-'):
                t.value = Alpha(t.value[2:], -1)
            else:
                t.value = Alpha(t.value[1:])
            return t
        except ValueError as e:
            print(e)

    @_(r'-p[0123]{1,4}', r'p[0123]{1,4}')
    def PAIR(self, t):
        try:
            if t.value.startswith('-'):
                t.value = Pair(Alpha(t.value[2:], -1))
            else:
                t.value = Pair(Alpha(t.value[1:]))
            return t
        except ValueError as e:
            print(e)

    def error(self, value):
        pass


class ArpyParser(Parser):
    tokens = ArpyLexer.tokens

    precedence = (
        ('left', '^', '/'),
        ('left', '\\', 'ALPHA'),
    )

    @_('expr')
    def result(self, p):
        return p.expr

    @_('"(" expr ")"')
    def expr(self, p):
        return p.expr

    @_('expr "^" expr')
    def expr(self, p):
        return wedge(p.expr0, p.expr1)

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

##############################################################################

ar_lexer = ArpyLexer()
ar_parser = ArpyParser()


def ar(user_input):
    '''
    Parse user input under a nicer API.
    Call ar_help() for details on syntax.
    '''
    return ar_parser.parse(ar_lexer.tokenize(user_input))


if __name__ == '__main__':
    while True:
        try:
            text = input('ξα > ')
        except (EOFError, KeyboardInterrupt):
            break
        if text:
            ar(text)
