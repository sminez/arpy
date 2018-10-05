'''
arpy (Absolute Relativity in Python)
Copyright (C) 2016-2018 Innes D. Anderson-Morrison All rights reserved.

Lexing and Parsing of a more mathematical syntax for performing calculations
with the arpy Absolute Relativity library.
'''
import re
from operator import add
from sys import _getframe, stderr
from collections import namedtuple
from itertools import permutations

from ..algebra.ar_types import Alpha, Pair
from ..algebra.config import ARConfig
from ..algebra.config import config as cfg
from ..algebra.multivector import MultiVector
from ..algebra.differential import differential_operator
from ..algebra.operations import full, div_by, div_into, project, \
        dagger, commutator


tags = [
    ('MVEC',  r'\{(.*)\}$'),
    ('DIFF',  r'<([p0213, -]*)\>'),
    ('ALPHA', r'-?a[0123]{1,4}|-?ap'),
    ('PAIR',  r'-?p[0123]{1,4}'),
    ('VAR',   r'-?[a-zA-Z_][a-zA-Z_0-9]*'),
    ('INDEX', r'[01234]')
]
literals = [
    ('PAREN_OPEN',  r'\('), ('PAREN_CLOSE',  r'\)'),
    ('ANGLE_OPEN',  r'\<'), ('ANGLE_CLOSE',  r'\>'),
    ('SQUARE_OPEN', r'\['), ('SQUARE_CLOSE', r'\]'),
    ('CURLY_OPEN',  r'\{'), ('CURLY_CLOSE', r'\}'),
    ('PLUS', r'\+'), ('MINUS', r'-'), ('COMMA', r','),
    ('FULL', r'\^'), ('BY', r'/'), ('INTO', r'\\'),
    ('DOT', r'\.'), ('DAG', r'!')
]

_tags = '|'.join('(?P<{}>{})'.format(t[0], t[1]) for t in tags + literals)
Token = namedtuple('token', ['tag', 'val'])


class AR_Error(Exception):
    pass


class ArpyLexer:
    tags = re.compile(_tags)
    literals = [tag_regex[0] for tag_regex in literals]

    def __init__(self, cfg=cfg, _globals=None):
        self._globals = _globals
        self.cfg = cfg
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

            if lex_tag == 'MVEC':
                alphas = re.split(', |,| ', text.strip())
                token = Token('EXPR', MultiVector(alphas, cfg=self.cfg))

            elif lex_tag == 'DIFF':
                alphas = re.split(', |,| ', text.strip())
                token = Token(
                    'EXPR',
                    differential_operator(alphas, cfg=self.cfg))

            elif lex_tag == 'ALPHA':
                if text.startswith('-'):
                    token = Token('EXPR', Alpha(text[2:], -1, cfg=self.cfg))
                else:
                    token = Token('EXPR', Alpha(text[1:], cfg=self.cfg))

            elif lex_tag == 'PAIR':
                if text.startswith('-'):
                    token = Token(
                        'EXPR', Pair(Alpha(text[2:], -1, cfg=self.cfg),
                                     cfg=self.cfg))
                else:
                    token = Token(
                        'EXPR', Pair(Alpha(text[1:], cfg=self.cfg),
                                     cfg=self.cfg))

            elif lex_tag == 'INDEX':
                token = Token('INDEX', int(text))

            elif lex_tag == 'VAR':
                is_negated = False
                if text.startswith('-'):
                    is_negated = True
                    text = text[1:]
                try:
                    # Use definitions from the context over the global values
                    val = eval(text, self.context_vars)
                except NameError:
                    try:
                        val = eval(text, self._globals)
                    except Exception:
                        stderr.write(
                            '"{}" is not currently defined\n'.format(text)
                        )
                        raise AR_Error()

                token = Token('EXPR', -val if is_negated else val)

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
    unops_postfix = {'DAG': dagger}
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

                elif token.tag in self.unops_postfix:
                    if previous_token is None:
                        msg = 'Missing argument to "{}" in "{}"\n'
                        stderr.write(msg.format(token.val, raw_text))
                        raise AR_Error()
                    else:
                        LHS, previous_token = previous_token, None
                        op = self.unops_postfix.get(token.tag)
                        val = op(LHS.val, cfg=self.cfg)
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
    def __init__(self, allowed=None, metric=None, div=None,
                 cfg=None, print_all=False):
        self._print = print_all
        if cfg is None:
            cfg = ARConfig(allowed, metric, div)
        self.cfg = cfg
        self._lexer = ArpyLexer(cfg=cfg)
        self._parser = ArpyParser(cfg=cfg)
        self._initialise_vars()

    def __repr__(self):
        return str(self.cfg)

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
            'h': MultiVector(_h, cfg=self.cfg),
            'q': MultiVector(_q, cfg=self.cfg),
            'B': MultiVector(_B, cfg=self.cfg),
            'E': MultiVector(_E, cfg=self.cfg),
            'F': MultiVector(_E + _B, cfg=self.cfg),
            'T': MultiVector(_T, cfg=self.cfg),
            'G': MultiVector(self.cfg.allowed, cfg=self.cfg),
            'zet_B': MultiVector(['p'] + _B, cfg=self.cfg),
            'zet_T': MultiVector(['0'] + _T, cfg=self.cfg),
            'zet_A': MultiVector([_h] + _A, cfg=self.cfg),
            'zet_E': MultiVector([_q] + _E, cfg=self.cfg),
            'Fp': MultiVector(['p'] + _B + _E, cfg=self.cfg),
            'zet_F': MultiVector(['p'] + _B + [_q] + _E, cfg=self.cfg),
            # Differentials
            'Dmu': differential_operator(['0', '1', '2', '3'], cfg=self.cfg),
            'DG': differential_operator(self.cfg.allowed, cfg=self.cfg),
            'DF': differential_operator(_B + _E, cfg=self.cfg),
            'DB': differential_operator(['p'] + _B, cfg=self.cfg),
            'DT': differential_operator(['0'] + _T, cfg=self.cfg),
            'DA': differential_operator([_h] + _A, cfg=self.cfg),
            'DE': differential_operator([_q] + _E, cfg=self.cfg),
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

    @property
    def allowed(self):
        return self.cfg.allowed

    @allowed.setter
    def allowed(self, allowed):
        if len(allowed) != 16:
            raise ValueError('Must provide all 16 elements for allowed')

        self.cfg.allowed = allowed
        self._initialise_vars()

    def decompose(self):
        '''Decompose the algebra into Zets'''
        # Bring the component definitions into scope
        self.cfg.update_env()

        # Define the additional components required
        quedgehog = 'a{}'.format([p for p in self.cfg.q][0].alpha.index)
        bases = (('zet_{}', 'ζ'), ('(zet_{}!)', 'ζ†'),
                 # This is easier than hacking the '-' into the correct place!
                 ('-zet_{}', '-ζ'), ('(-zet_{}!)', '-ζ†'))
        zets = ['B', 'T', 'A', 'E']

        def _decompose_zet(base_zet, zet):
            # Pull out the target set of alphas (always +ve)
            target = {z[0] for z in self('zet_{}'.format(zet)).iter_alphas()}

            for tmp, str_rep in bases:
                base = tmp.format(base_zet)
                # Try all versions but try to float a0 to the front and the
                # quedgehog to the back if possible.
                all_comps = [
                    ['a0', base], [base, quedgehog], ['a0', base, quedgehog]
                ]
                for components in all_comps:
                    for permutation in permutations(components):
                        expr = ' ^ '.join(permutation)
                        res = self(expr)
                        signs = [p.xi.sign for p in res]
                        if all(map(lambda s: s == 1, signs)):
                            candidate = {z[0] for z in res.iter_alphas()}
                            if candidate == target:
                                return '{} = {}'.format(
                                    zet, expr.replace(base, str_rep))

            raise AR_Error('Unable to decompose in terms of {}'.format(
                base_zet))

        # Attempt to decompose everything in terms of everything else!
        print(self)
        print('-' * 20)

        for base_zet in zets:
            decompositions = []
            for zet in filter(lambda z: z != base_zet, zets):
                decompositions.append(_decompose_zet(base_zet, zet))

            print('zet_{} = ζ'.format(base_zet))
            for decomp in decompositions:
                print(decomp)
            print('-' * 20)

    def __call__(self, text, *, cancel_terms=False):
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
            result = result.val
            if cancel_terms and isinstance(result, MultiVector):
                result.cancel_terms()

            if self._print:
                print('"{}": {}'.format(text, result))

            return result

    # Allow ARContext to be used as a context manager
    def __enter__(self):
        self.cfg.update_env()
        return self

    def __exit__(self, *args):
        pass
