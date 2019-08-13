"""
Attempt at making this work:

    $ python3 -m arpy <calculation_file> [--vector --simplify --latex]

The file is fed through and ARContext with some aditional pre-parsing in
order to deal with assignment and paramater setting. All Variables are printed
when at the head of the output and calculations are labelled.
"""
import argparse
from .utils.calc_file import run_calculation


description = """\
.: arpy :: Absolute Relativity and the Algebra of Reality :.
------------------------------------------------------------
"""

epilog = """\
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
"""


parser = argparse.ArgumentParser(
    description=description, epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter
)

parser.add_argument(
    "-v", "--vector", action="store_true", help="display the result in vector calculus notation"
)
# parser.add_argument(
#     '-s',
#     '--simplify',
#     action='store_true',
#     help='simplify results before printing'
# )
parser.add_argument(
    "-l", "--latex", action="store_true", help="print results as LaTex instead of unicode"
)
parser.add_argument("script")
args = parser.parse_args()

# Read in the contents of the script
with open(args.script, "r") as f:
    script = [s.strip() for s in f.readlines() if s != "\n"]

modifier = ""

if args.vector:
    modifier = ".v"

if args.latex:
    modifier += ".__tex__()"


output = run_calculation(script, modifier)
print("\n".join(output))
