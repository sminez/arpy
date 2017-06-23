'''
Attempt at making this work:

    $ python3 -m arpy <calculation_file> [--vector --simplify --latex]

The file is fed through and ARContext with some aditional pre-parsing in
order to deal with assignment and paramater setting. All Variables are printed
when at the head of the output and calculations are labelled.
'''
import re
import argparse
from collections import namedtuple
from .utils.lexparse import ARContext
from .algebra.config import config

from . import *  # Bring in all of arpy


description = '''\
.: arpy :: Absolute Relativity and the Algebra of Reality :.
------------------------------------------------------------
'''

epilog = '''\
------------------------------------------------------------------------------

This is the command line interface for arpy. Input should be provided in the
form of a *.arp calculation file. The format of an arp file is as follows:

```
// ALLOWED: p 23 31 12 0 023 031 012 123 1 2 3 0123 01 02 03
// METRIC: +---

# Define a multivector
odd = {1, 2, 3, 023, 031, 012}

# Compute FDF
DF = Dmu ^ F
FDF = F ^ DF

# Compute FFdagger
FFdag = F ^ F!
```

Lines begining "//" set paramaters for the calculation. At present the only
permitted values are "ALLOWED" and "METRIC" as in the example.

Lines begining "#" are calculation comments and will be printed as given.

Assigning a result to a name will compute a result and display it.
'''


raw = namedtuple('raw', 'lnum var')
comment = namedtuple('raw', 'lnum text')
step = namedtuple('step', 'lnum var args')
context_update = namedtuple('context_update', 'lnum param val')
mvec_def = namedtuple('mvec_def', 'lnum var alphas')

mvec_pattern = r'([a-zA-Z_][a-zA-Z_0-9]*)\s?=\s?\{([p0213, -]*)\}$'
modifier_map = {
    'DEL NOTATION': '.v',
    'SIMPLIFIED': '.simplified()',
    'TEX': '.__tex__()'
    }


def parse_calculation_file(fname, default_allowed=config.allowed,
                           default_metric=config.metric):
    '''Configure the ARContext for carrying out the calculation'''
    def convert_metric(s):
        s = s[:-1]  # remove trailing newline
        if not all([c in '+-' for c in s]):
            raise CalculationError('Invalid metric: {}', s)

        metric = [1 if c == '+' else -1 for c in s]
        return tuple(metric)

    # Set paramaters to default to start
    metric, allowed = None, None
    lines = []
    modifiers = {}

    with open(fname, 'r') as f:
        for lnum, line in enumerate(f):
            lnum += 1

            if line == '\n':
                lines.append(comment(lnum, ''))

            # Check and set paramaters
            elif line.startswith('// METRIC:'):
                m = convert_metric(line.split('// METRIC: ')[1])
                if metric is None:
                    metric = m
                lines.append(comment(lnum, line.strip()))
                lines.append(context_update(lnum, 'metric', m))

            elif line.startswith('// ALLOWED:'):
                a = line.split('// ALLOWED: ')[1].split()
                if allowed is None:
                    allowed = a
                lines.append(comment(lnum, line.strip()))
                lines.append(context_update(lnum, 'allowed', a))

            elif line.startswith('// '):
                action = line[3:].strip()
                modifiers[lnum + 1] = modifier_map[action]

            # extract comments
            elif line.startswith('#'):
                lines.append(comment(lnum, line.strip()))

            # extract steps
            else:
                if '=' not in line:
                    lines.append(raw(lnum, line.strip()))
                else:
                    # Check for multivector assignent
                    mvec_match = re.match(mvec_pattern, line)
                    if mvec_match:
                        var, alphas = mvec_match.groups()
                        alphas = re.split(', |,| ', alphas.strip())
                        lines.append(mvec_def(lnum, var, alphas))
                    else:
                        # Try to parse an ar command
                        var, args = line.split(' = ')
                        lines.append(step(lnum, var, args.strip()))

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


parser = argparse.ArgumentParser(
    description=description,
    epilog=epilog,
    formatter_class=argparse.RawDescriptionHelpFormatter
)

parser.add_argument(
    '-v',
    '--vector',
    action='store_true',
    help='display the result in vector calculus notation'
)
# parser.add_argument(
#     '-s',
#     '--simplify',
#     action='store_true',
#     help='simplify results before printing'
# )
parser.add_argument(
    '-l',
    '--latex',
    action='store_true',
    help='print results as LaTex instead of unicode'
)
parser.add_argument('script')
args = parser.parse_args()

modifier = ''
if args.vector:
    modifier = '.v'
if args.latex:
    modifier += '.__tex__()'

context, lines, modifiers = parse_calculation_file(args.script)

for l in lines:
    step_modifier = modifiers.get(l.lnum)

    if isinstance(l, comment):
        print(l.text)

    elif isinstance(l, raw):
        mod = step_modifier if step_modifier else modifier
        eval('''print('{} = ', {}{})'''.format(l.var, l.var, mod))

    elif isinstance(l, context_update):
        # The update will have a matching comment line to show when
        # it occured in the calculation.
        if l.param == 'metric':
            context.metric = l.val
        elif l.param == 'allowed':
            context.allowed = l.val

    elif isinstance(l, mvec_def):
        exec('{} = MultiVector({})'.format(l.var, l.alphas))
        mod = step_modifier if step_modifier else modifier
        eval('''print('{} = ', {}{})'''.format(l.var, l.var, mod))

    elif isinstance(l, step):
        exec('{} = context("{}")'.format(l.var, l.args))
        print('{} = {}'.format(l.var, l.args))
        mod = step_modifier if step_modifier else modifier
        eval('print({}{})'.format(l.var, mod))
