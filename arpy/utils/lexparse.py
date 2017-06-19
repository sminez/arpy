'''
arpy (Absolute Relativity in Python)
Copyright (C) 2016-2017 Innes D. Anderson-Morrison All rights reserved.

Lexing and Parsing of a more mathematical syntax for performing calculations
with the arpy Absolute Relativity library.
'''
import re
from operator import add
from sys import _getframe, stderr
from collections import namedtuple
from ..algebra.ar_types import Alpha, Pair
from ..algebra.config import config as cfg
from ..algebra.multivector import MultiVector
from ..algebra.differential import differential_operator
from ..algebra.operations import full, div_by, div_into, project, \
        dagger, commutator


tags = [
    ('ALPHA', r'-?a[0123]{1,4}|-?ap'),
    ('PAIR',  r'-?p[0123]{1,4}'),
    ('VAR',   r'[a-zA-Z_][a-zA-Z_0-9]*'),
    ('INDEX', r'[01234]')
]
literals = [
    ('PAREN_OPEN',  r'\('), ('PAREN_CLOSE',  r'\)'),
    ('ANGLE_OPEN',  r'\<'), ('ANGLE_CLOSE',  r'\>'),
    ('SQUARE_OPEN', r'\['), ('SQUARE_CLOSE', r'\]'),
    ('PLUS', r'\+'), ('COMMA', r','), ('FULL', r'\^'),
    ('BY', r'/'), ('INTO', r'\\'), ('DOT', r'\.'),
    ('DAG', r'!')
]

_tags = '|'.join('(?P<{}>{})'.format(t[0], t[1]) for t in tags + literals)
Token = namedtuple('token', ['tag', 'val'])


class AR_Error(Exception):
    pass


class ArpyLexer:
    tags = re.compile(_tags)
    literals = [tag_regex[0] for tag_regex in literals]

    def __init__(self, _globals=None):
        self._globals = _globals
        self.context_vars = {}

    def lex(self, string, context_vars=None):
        string = re.sub(' \t\n', '', string)  # remove whitespace
        matched_text = []

        if context_vars:
            self.context_vars = context_vars

        for match in re.finditer(self.tags, string):
            lex_tag = match.lastgroup
            group = [g for g in match.groups() if g is not None]
            text = group[1] if len(group) == 2 else match.group(lex_tag)
            matched_text.append(text)

            if lex_tag == 'ALPHA':
                if text.startswith('-'):
                    token = Token('EXPR', Alpha(text[2:], -1))
                else:
                    token = Token('EXPR', Alpha(text[1:]))
            elif lex_tag == 'PAIR':
                if text.startswith('-'):
                    token = Token('EXPR', Pair(Alpha(text[2:], -1)))
                else:
                    token = Token('EXPR', Pair(Alpha(text[1:])))
            elif lex_tag == 'INDEX':
                token = Token('INDEX', int(text))
            elif lex_tag == 'VAR':
                try:
                    # Use definitions from the context over the global values
                    token = Token('EXPR', eval(text, self.context_vars))
                except NameError:
                    try:
                        token = Token('EXPR', eval(text, self._globals))
                    except:
                        stderr.write(
                            '"{}" is not currently defined\n'.format(text)
                        )
                        raise AR_Error()
            elif lex_tag in self.literals:
                token = Token(lex_tag, text)
            else:
                message = (
                    'Input contains invalid syntax for the ar() function: {}\n'
                )
                stderr.write(message.format(text))
                raise AR_Error()

            yield token


class ArpyParser:
    unops = {'DAG': dagger}
    binops = {'FULL': full, 'BY': div_by, 'INTO': div_into, 'PLUS': add}

    def __init__(self, cfg=cfg):
        self.cfg = cfg
        self.context_vars = {}

    def sub_expr(self, tokens, deliminator_tag):
        '''Pull tokens until we hit the specified deliminator.'''
        token = next(tokens)
        sub_expression = []
        while token.tag != deliminator_tag:
            try:
                sub_expression.append(token)
                token = next(tokens)
            except StopIteration:
                message = 'Invalid subexpression: missing "{}"\n'
                stderr.write(message.format(deliminator_tag))
                raise AR_Error()
        return (s for s in sub_expression)

    def parse(self, tokens, raw_text, compound=[], context_vars=None):
        '''Naive recursive decent parsing of the input.'''
        previous_token = None

        try:
            while True:
                token = next(tokens)

                if token.tag == 'PAREN_OPEN':
                    sub_expression = self.sub_expr(tokens, 'PAREN_CLOSE')
                    sub_expr_token = self.parse(sub_expression, raw_text)
                    if previous_token:
                        val = full(
                            previous_token.val, sub_expr_token.val,
                            cfg=self.cfg
                        )
                        previous_token = Token('EXPR', val)
                    else:
                        previous_token = sub_expr_token

                elif token.tag == 'EXPR':
                    if previous_token:
                        # default to forming the full product
                        val = full(
                            previous_token.val, token.val,
                            cfg=self.cfg
                        )
                        previous_token = Token('EXPR', val)
                    else:
                        # Store the token and then check the next token to
                        # determine what we should do next.
                        previous_token = token

                elif token.tag in self.binops:
                    if previous_token is None:
                        msg = 'Missing left argument to "{}" in "{}"\n'
                        stderr.write(msg.format(token.val, raw_text))
                        raise AR_Error()
                    else:
                        LHS, previous_token = previous_token, None
                        op = self.binops.get(token.tag)
                        RHS = self.parse(tokens, raw_text)
                        if RHS is None:
                            err = 'Missing right argument to "{}" in "{}"\n'
                            stderr.write(err.format(token.val, raw_text))
                            raise AR_Error()
                        if RHS.tag != 'EXPR':
                            # BinaryOps take a single LHS and RHS expression
                            msg = 'Invalid argument for {}: {}\n'
                            stderr.write(msg.format(token.val, RHS.val))
                            raise AR_Error()
                        else:
                            if token.tag == 'PLUS':
                                val = op(LHS.val, RHS.val)
                            else:
                                val = op(LHS.val, RHS.val, cfg=self.cfg)
                            previous_token = Token('EXPR', val)

                elif token.tag == 'ANGLE_OPEN':
                    sub_expression = self.sub_expr(tokens, 'ANGLE_CLOSE')
                    arg = self.parse(sub_expression, raw_text)

                    index = next(tokens)
                    if index.tag != 'INDEX':
                        msg = 'Missing index for projection: {}\n'
                        stderr.write(msg.format(raw_text))
                        raise AR_Error()
                    else:
                        val = project(arg.val, index.val)
                        previous_token = Token('EXPR', val)

                elif token.tag == 'SQUARE_OPEN':
                    sub_expression = self.sub_expr(tokens, 'COMMA')
                    LHS = self.parse(sub_expression, raw_text)
                    sub_expression = self.sub_expr(tokens, 'SQUARE_CLOSE')
                    RHS = self.parse(sub_expression, raw_text)
                    val = commutator(LHS.val, RHS.val)
                    previous_token = Token('EXPR', val)

                else:
                    stderr.write('Invalid input: {}\n'.format(raw_text))
                    return None

        except StopIteration:
            if previous_token:
                return previous_token
            else:
                stderr.write('Unable to parse input: {}\n'.format(raw_text))
                return None
        except AR_Error:
            return None


class ARContext:
    '''
    User interface class for working with the library.
    Create an instance and then use by calling ar as a function.
    i.e.
    >>> ar = ARContext()
    >>> ar("a12 ^ a23")
    >>> α31
    '''
    def __init__(self, cfg=cfg):
        self.cfg = cfg
        self._lexer = ArpyLexer()
        self._parser = ArpyParser(cfg=cfg)
        self._initialise_vars()

    def _initialise_vars(self):
        '''Set all of the standard variables'''
        # Check that we have a (roughly) valid set of values
        _h = [a for a in self.cfg.allowed if len(a) == 3 and '0' not in a]
        assert len(_h) == 1, 'h is a single element: {}'.format(_h)
        _h = _h[0]
        _q = [a for a in self.cfg.allowed if len(a) == 4]
        assert len(_q) == 1, 'q is a single element: {}'.format(_q)
        _q = _q[0]
        _B = [a for a in self.cfg.allowed if len(a) == 2 and '0' not in a]
        assert len(_B) == 3, 'B is a 3-vector: {}'.format(_B)
        _T = [a for a in self.cfg.allowed if len(a) == 3 and '0' in a]
        assert len(_T) == 3, 'T is a 3-vector: {}'.format(_T)
        _A = [a for a in self.cfg.allowed if len(a) == 1 and a not in 'p0']
        assert len(_A) == 3, 'A is a 3-vector: {}'.format(_A)
        _E = [a for a in self.cfg.allowed if len(a) == 2 and '0' in a]
        assert len(_E) == 3, 'E is a 3-vector: {}'.format(_E)

        self._vars = {
            # Multivectors
            'h': MultiVector(_h, cfg=cfg),
            'q': MultiVector(_q, cfg=cfg),
            'B': MultiVector(_B, cfg=cfg),
            'E': MultiVector(_E, cfg=cfg),
            'F': MultiVector(_E + _B, cfg=cfg),
            'T': MultiVector(_T, cfg=cfg),
            'G': MultiVector(self.cfg.allowed, cfg=cfg),
            'B4': MultiVector(['p'] + _B, cfg=cfg),
            'T4': MultiVector(['0'] + _T, cfg=cfg),
            'A4': MultiVector([_h] + _A, cfg=cfg),
            'E4': MultiVector([_q] + _E, cfg=cfg),
            'Fp': MultiVector(['p'] + _B + _E, cfg=cfg),
            'Fpq': MultiVector(['p'] + _B + [_q] + _E, cfg=cfg),
            # Differentials
            'DG': differential_operator(self.cfg.allowed, cfg=cfg),
            'DF': differential_operator(_B + _E, cfg=cfg),
            'DB': differential_operator(['p'] + _B, cfg=cfg),
            'DT': differential_operator(['0'] + _T, cfg=cfg),
            'DA': differential_operator([_h] + _A, cfg=cfg),
            'DE': differential_operator([_q] + _E, cfg=cfg),
        }

    @property
    def metric(self):
        return self.cfg.metric

    @metric.setter
    def metric(self, signs):
        if all(sign in ["+", "-"] for sign in signs):
            if len(signs) != 4:
                raise ValueError(
                    "metric should be a 4 element string.\n"
                    "i.e. 'ar.metric = \"+---\"'"
                )
            metric = tuple(1 if s == "+" else -1 for s in signs)
        elif all(sign in [1, -1] for sign in signs):
            metric = signs
        else:
            raise TypeError("metric must be comprised of +/- only")

        self.cfg.metric = metric
        self._parser.metric = metric

    @property
    def allowed(self):
        return self.cfg.allowed

    @allowed.setter
    def allowed(self, allowed):
        if len(allowed) != 16:
            raise ValueError('Must provide all 16 elements for allowed')

        self.cfg.allowed = allowed
        self._parser.allowed = allowed
        self._initialise_vars()

    def __call__(self, text):
        # NOTE:: The following is a horrible hack that allows you to
        #        inject local variables into the parser.
        stack_frame = _getframe(1)
        self._lexer._globals = stack_frame.f_locals
        try:
            result = self._parser.parse(
                self._lexer.lex(text, context_vars=self._vars), text)
        except AR_Error:
            return None
        if result:
            # Result is an internal Token so pull of the value for returning
            # If there was an error we have printed the error and returned None
            return result.val
