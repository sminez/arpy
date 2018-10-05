import re
from io import StringIO
from collections import namedtuple

from ..algebra.config import config
from .lexparse import ARContext
from arpy import *


mvec_pattern = r'([a-zA-Z_][a-zA-Z_0-9]*)\s?=\s?\{(.*)\}$'
operator_pattern = r'([a-zA-Z_][a-zA-Z_0-9]*)\s?=\s?\<([p0213, -]*)\>$'
modifier_map = {
    'DEL NOTATION': '.v',
    'SIMPLIFIED': '.simplified()',
    'TEX': '.__tex__()'
    }

raw = namedtuple('raw', 'lnum var')
comment = namedtuple('raw', 'lnum text')
step = namedtuple('step', 'lnum var args')
context_update = namedtuple('context_update', 'lnum param val')
mvec_def = namedtuple('mvec_def', 'lnum var alphas')
operator_def = namedtuple('op_def', 'lnum var alphas')


def parse_calculation_script(script, default_allowed=config.allowed,
                             default_metric=config.metric):
    '''
    Parse the contents of a calculation file and onfigure the ARContext for
    carrying out the calculation.

    `script` should be an iterable of valid calcualtion directive lines.
    '''
    def convert_metric(s):
        if not all([c in '+-' for c in s]):
            raise RuntimeError('Invalid metric: {}', s)

        metric = [1 if c == '+' else -1 for c in s]
        return tuple(metric)

    # Set paramaters to default to start
    metric, allowed = None, None
    lines = []
    modifiers = {}

    for lnum, line in enumerate(script):
        line = line.strip()
        lnum += 1

        if line == '':
            lines.append(comment(lnum, ''))

        # Check and set paramaters
        elif line.startswith('// METRIC:'):
            m = convert_metric(line.split('// METRIC: ')[1])
            if metric is None:
                metric = m
            lines.append(comment(lnum, line))
            lines.append(context_update(lnum, 'metric', m))

        elif line.startswith('// ALLOWED:'):
            a = line.split('// ALLOWED: ')[1].split()
            if allowed is None:
                allowed = a
            lines.append(comment(lnum, line))
            lines.append(context_update(lnum, 'allowed', a))

        elif line.startswith('// '):
            action = line[3:].strip()
            modifiers[lnum + 1] = modifier_map[action]

        # extract comments
        elif line.startswith('#'):
            lines.append(comment(lnum, line))

        # extract steps
        else:
            if '=' not in line:
                lines.append(raw(lnum, line))
            else:
                # Check for multivector assignent
                mvec_match = re.match(mvec_pattern, line)
                operator_match = re.match(operator_pattern, line)
                if mvec_match:
                    var, alphas = mvec_match.groups()
                    alphas = re.split(', |,| ', alphas.strip())
                    lines.append(mvec_def(lnum, var, alphas))
                elif operator_match:
                    var, alphas = operator_match.groups()
                    alphas = re.split(', |,| ', alphas.strip())
                    lines.append(operator_def(lnum, var, alphas))
                else:
                    # Try to parse an ar command
                    var, args = line.split(' = ')
                    lines.append(step(lnum, var, args))

    # Fall back to defaults if metric/allowed were not specified
    if metric is None:
        metric = default_metric
        m = ''.join('+' if x == 1 else '-' for x in metric)
        lines = [comment(0, '// METRIC: ' + m)] + lines

    if allowed is None:
        allowed = default_allowed
        lines = [comment(0, '// ALLOWED: ' + ' '.join(allowed))] + lines

    config.allowed = allowed
    config.metric = metric
    context = ARContext(cfg=config)
    return context, lines, modifiers


def run_calculation(script, modifier=''):
    '''
    Run a calculation from a list of calculation lines.

    NOTE: This makes extensive use of eval!
    '''
    context, lines, modifiers = parse_calculation_script(script)
    output = StringIO()

    for l in lines:
        step_modifier = modifiers.get(l.lnum)

        if isinstance(l, comment):
            print(l.text, file=output)

        elif isinstance(l, raw):
            mod = step_modifier if step_modifier else modifier
            eval('''print('{} = ', context("{}"){}, file=output)'''.format(
                l.var, l.var, mod))

        elif isinstance(l, context_update):
            # The update will have a matching comment line to show when
            # it occured in the calculation.
            if l.param == 'metric':
                context.metric = l.val
            elif l.param == 'allowed':
                context.allowed = l.val

        elif isinstance(l, mvec_def):
            # exec('{} = MultiVector({})'.format(l.var, l.alphas))
            exec('%s = context("{%s}")' % (l.var, ' '.join(l.alphas)))
            mod = step_modifier if step_modifier else modifier
            eval('''print('{} = ', {}{}, file=output)'''.format(
                l.var, l.var, mod))

        elif isinstance(l, operator_def):
            exec('{} = context("<{}>")'.format(l.var, ' '.join(l.alphas)))
            eval('''print('{} = ', {}, file=output)'''.format(l.var, l.var))

        elif isinstance(l, step):
            print('{} = context("{}")'.format(l.var, l.args))
            exec('{} = context("{}")'.format(l.var, l.args))
            print('{} = {}'.format(l.var, l.args), file=output)
            mod = step_modifier if step_modifier else modifier
            eval('print({}{}, file=output)'.format(l.var, mod))

    # Split back into lines and remove the trailing newline so that
    # the caller can decide how to format the result.
    return output.getvalue().split('\n')[:-1]
