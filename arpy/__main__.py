'''
Attempt at making this work:

    $ python3 -m arpy <calculation_file> [--vector --simplify --latex]

The file is fed through and ARContext with some aditional pre-parsing in
order to deal with assignment and paramater setting. All Variables are printed
when at the head of the output and calculations are labelled.
'''
# import sys
import argparse
from collections import namedtuple
from .utils.lexparse import ARContext
from .algebra.config import ALLOWED, METRIC
from . import *


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


class CalculationError(Exception):
    pass


comment = namedtuple('comment', 'lnum text')
step = namedtuple('step', 'lnum var args')


def parse_calculation_file(fname):
    '''Configure the ARContext for carrying out the calculation'''
    def convert_metric(s):
        s = s[:-1]  # remove trailing newline
        if not all([c in '+-' for c in s]):
            raise CalculationError('Invalid metric: {}', s)

        metric = [1 if c == '+' else -1 for c in s]
        return tuple(metric)

    # Set paramaters to default to start
    metric, allowed = METRIC, ALLOWED
    lines = []

    with open(fname, 'r') as f:
        for lnum, line in enumerate(f):
            # Skip blank lines
            if line == '\n':
                continue
            # Check and set paramaters
            elif line.startswith('// METRIC:'):
                metric = convert_metric(line.split('// METRIC: ')[1])
                lines.append(comment(lnum, line[:-1]))
            elif line.startswith('// ALLOWED:'):
                allowed = line.split('// ALLOWED: ')[1].split()
                lines.append(comment(lnum, line[:-1]))
            # extract comments
            elif line.startswith('#'):
                if line.startswith('# '):
                    text = line.split('# ')[1]
                else:
                    text = line.split('#')[1]
                lines.append(comment(lnum, text[:-1]))
            # extract steps
            else:
                if '=' not in line:
                    tmp = 'Steps must assign to a variable:\n[{}] {}'
                    raise CalculationError(tmp.format(lnum, line))
                else:
                    var, args = line.split(' = ')
                    lines.append(step(lnum, var, args[:-1]))

    context = ARContext(metric=metric, allowed=allowed)
    return context, lines


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
parser.add_argument(
    '-s',
    '--simplify',
    action='store_true',
    help='simplify results before printing'
)
parser.add_argument(
    '-l',
    '--latex',
    action='store_true',
    help='print results as LaTex instead of unicode'
)
parser.add_argument('script')
args = parser.parse_args()


context, lines = parse_calculation_file(args.script)

for l in lines:
    if isinstance(l, comment):
        print(l.text)
    elif isinstance(l, step):
        exec('{} = context("{}")'.format(l.var, l.args))
        print('{} = {}'.format(l.var, l.args))
        if args.latex:
            eval('''print('{} = ', {}.__tex__())'''.format(l.var, l.var))
        else:
            eval('''print('{} = ', {})'''.format(l.var, l.var))
        # Split calculation steps with newlines
        print('')
