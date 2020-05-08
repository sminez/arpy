[![Build Status](https://travis-ci.org/sminez/arpy.svg?branch=master)](https://travis-ci.org/sminez/arpy) [![Coverage Status](https://coveralls.io/repos/github/sminez/arpy/badge.svg?branch=master)](https://coveralls.io/github/sminez/arpy?branch=master)

arpy (Absolute Relativity in Python)

Copyright (C) 2016-2020 Innes D. Anderson-Morrison All rights reserved.

![Cayley Table for AR](readme_icon.png)


## Overview
arpy is a module for performing calculations within the theory of Absolute Relativity
as devised by [Dr J.G.Williamson](https://quicycle.com).

This repository is under active development. If you have any questions or
suggestions for features / bug fixes please [open an issue](https://github.com/sminez/arpy/issues).

### Installing and updating arpy
To install the module on your system you will need Python version 3.6 or
greater. From this directory run `make install` and the script
will take care of the rest.

In order to update the the latest version of the code it is _strongly_ advised
to clone this repo using git rather than manually downloading a zip each time
- Though the second option _is_ possible if using git is not something you want
to do. To update using git run the following commands from the root directory of
this repo (where this README file is located):

```bash
$ git pull
$ make install
```

After install, the module can be imported into a Python repl session in the
usual way. While it is generally regarded as bad practice when developing
stand-alone programs, the `from arpy import *` command is the recommended way to
work with arpy when in an interactive terminal session.

### Installing the Jupyter qtconsole
The default Python command line interpreter leaves a lot to be desired. It is
strongly recommended that you also install the [Jupyter QT console](https://qtconsole.readthedocs.io/en/latest/)
for working with arpy as it provides syntax highlighting, interactive help and
access to your system shell. Install using the following commands:

```bash
$ python3 -m pip install python3-pyqt5
$ python3 -m pip install jupyter
```

Start the console using:
```bash
$ python3 -m qtconsole
# or if you want a larger font size: in this example 18pt
$ python3 -m qtconsole --JupyterWidget.font_size=18
```


### Reporting issues and requesting features
Github has a built in issue tracker that makes it much easier to keep on top of
changes and bug fixes. In either case, please [create a new issue](https://github.com/sminez/arpy/issues) using this link or by going to the
`Issues` panel at the top of the page in Github.

(If it is a bug report, please copy in and error messages you get when describing
the problem.)


### Running the test suite
If you would like to run the current test suite (you're working on a patch /
you're just curious!) then simply run the following command in the root
directory of the repo:
```bash
$ make test
```

If you are trying to write new tests for the module then please follow the
example of the tests written so far. I am using [pytest](http://doc.pytest.org/en/latest/)
as the test runner and wherever possible I am trying to paramaterize test cases
to cover all possible inputs to eliminate edge cases as far as possible.
All tests should have a clear docstring description of what it is the test is
attempting to prove about the code.
