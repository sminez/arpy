'''
Lexing and Parsing of a more mathematical syntax for performing calculations
with the arpy Absolute Relativity library.

SLY is rather magical - docs are here:
    http://sly.readthedocs.io/en/latest/
    https://github.com/dabeaz/sly

TODO::
    Differentials
    Help text
        This might be as simple as matching '? <word>' or simply '?'
    Keywords to swap out the metric and allowed values?
        'SET metric +---'
    Keyword to dump calculation to file?
        Worthwhile keeping a buffer of the computation so far?
    Keyword to load equations that have been saved
    Keyword to dump expr for use with further real valued calulations
'''
from sly import Lexer, Parser
from algebra.ar_types import Alpha, Pair, Xi_vecs
from algebra.operations import wedge, div_by, div_into
from algebra.differential import Dmu, show


class ArpyLexer(Lexer):
    tokens = {'NAME', 'ALPHA', 'PAIR', 'DMU', 'SHOW', 'VECSHOW'}
    ignore = ' \t'
    literals = set(r'= + ^ . * / \ ( ) < > { } [ ]'.split())

    # Token regex definitions
    DMU = r'Dmu'
    SHOW = r'show'
    VECSHOW = r'vshow'
    NAME = r'[A-Z][a-zA-Z0-9_]*'

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

    @_(r'\n+')
    def newline(self, t):
        self.lineno += t.value.count('\n')
        return ""

    def error(self, value):
        pass


class ArpyParser(Parser):
    tokens = ArpyLexer.tokens

    precedence = (
        ('left', '=', '^'),
        ('left', '/', '\\'),
        ('left', 'NAME', 'DMU'),
    )

    def __init__(self):
        self.names = Xi_vecs

    @_('NAME "=" expr')
    def statement(self, p):
        self.names[p.NAME] = p.expr
        print('{} -> {}'.format(p.NAME, p.expr))

    @_('DMU expr')
    def expr(self, p):
        try:
            return Dmu(p.expr)
        except:
            print("Invalid differential")

    @_('SHOW DMU expr')
    def statement(self, p):
        try:
            show(Dmu(p.expr), False)
        except:
            print("Invalid differential")

    @_('VECSHOW DMU expr')
    def statement(self, p):
        try:
            show(Dmu(p.expr), True)
        except:
            print("Invalid differential")

    @_('expr')
    def statement(self, p):
        print(p.expr)

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

    @_('NAME')
    def expr(self, p):
        try:
            return self.names[p.NAME]
        except LookupError:
            print("Undefined name '%s'" % p.NAME)
            return "ERROR"


if __name__ == '__main__':
    lexer = ArpyLexer()
    parser = ArpyParser()
    while True:
        try:
            text = input('ξα > ')
        except (EOFError, KeyboardInterrupt):
            break
        if text:
            parser.parse(lexer.tokenize(text))
