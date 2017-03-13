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
from ..algebra.config import ALLOWED, METRIC
from ..algebra.ar_types import Alpha, Pair
from ..algebra.operations import full, div_by, div_into, project, commutator


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
    ('PLUS', r'\+'), ('COMMA', r','), ('WEDGE', r'\^'),
    ('BY', r'/'), ('INTO', r'\\'), ('DOT', r'\.')
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

    def lex(self, string):
        string = re.sub(' \t\n', '', string)  # remove whitespace
        matched_text = []

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
                    token = Token('EXPR', eval(text, self._globals))
                except NameError:
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
    metric = METRIC
    allowed = ALLOWED
    operations = {'WEDGE': full, 'BY': div_by, 'INTO': div_into, 'PLUS': add}

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

    def parse(self, tokens, raw_text, compound=[]):
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
                            metric=self.metric, allowed=self.allowed
                        )
                        previous_token = Token('EXPR', val)
                    else:
                        previous_token = sub_expr_token

                elif token.tag == 'EXPR':
                    if previous_token:
                        # default to forming the full product
                        val = full(
                            previous_token.val, token.val,
                            metric=self.metric, allowed=self.allowed
                        )
                        previous_token = Token('EXPR', val)
                    else:
                        # Store the token and then check the next token to
                        # determine what we should do next.
                        previous_token = token

                elif token.tag in self.operations:
                    if previous_token is None:
                        msg = 'Missing left argument to "{}" in "{}"\n'
                        stderr.write(msg.format(token.val, raw_text))
                        raise AR_Error()
                    else:
                        LHS, previous_token = previous_token, None
                        op = self.operations.get(token.tag)
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
                                val = op(
                                    LHS.val, RHS.val, metric=self.metric,
                                    allowed=self.allowed
                                )
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
    def __init__(self, metric=METRIC, allowed=ALLOWED):
        self._metric = METRIC
        self._allowed = ALLOWED
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

    def __call__(self, text):
        # NOTE:: The following is a horrible hack that allows you to
        #        inject local variables into the parser.
        stack_frame = _getframe(1)
        self._lexer._globals = stack_frame.f_locals
        try:
            result = self._parser.parse(self._lexer.lex(text), text)
        except AR_Error:
            return None
        if result:
            # Result is an internal Token so pull of the value for returning
            # If there was an error we have printed the error and returned None
            return result.val
