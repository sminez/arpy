Working With The Dynamics Of A Relativistic Fluid
=================================================

![Cayley Table for the Williamson Algebra](readme_icon.png)

The code in this repository is a work to calculate with and simulate the dynamics
of the relativistic fluid theorised by [Dr J.G.Williamson](http://www.gla.ac.uk/schools/engineering/staff/johnwilliamson/).

### Installing and updating arpy
To install the module on your system you will need Python version 3.5 or
greater. From this directory run `sudo python3 setup.py install` and the script
will take care of the rest.

In order to update the the latest version of the code it is _strongly_ advised
to [clone this repo using git](https://confluence.atlassian.com/bitbucket/create-and-clone-a-repository-800695642.html#Createandclonearepository-CloningaGitrepository)
rather than manually downloading the zip each time - though the second option
_is_ possible if using git is not something you want to do. To update using git
run the following commands from the root directory of this repo (where this
README file is located):

```bash
$ git pull   # You will be prompted for your BitBucket login details
$ git checkout master
$ sudo python3 setup.py install
```

To update via .zip, go to the [downloads page](https://bitbucket.org/sminez/arpy/downloads/)
and click on `Download repository`. This will give you a zip file that you can
extract and run the `sudo python3 setup.py install` command in.

After either method, the module can be imported into a Python repl session in
the usual way. While it is generally regarded as bad practice when developing
stand-alone programs, the `from arpy import *` command is the recommended way to
work with arpy when in an interactive terminal session.

### Installing the Jupyter qtconsole
The default Python command line interpreter leaves a lot to be desired. It is
strongly recommended that you also install the [Jupyter QT console](https://qtconsole.readthedocs.io/en/latest/)
for working with arpy as it provides syntax highlighting, interactive help and
access to your system shell. Install using the following commands:

```bash
$ sudo python3 -m pip install python3-pyqt5
$ sudo python3 -m pip install qtconsole
```

Start the console using:
```bash
$ python3 -m qtconsole
# or if you want a larger font size: in this example 18pt
$ python3 -m qtconsole --JupyterWidget.font_size=18
```


### Reporting issues and requesting features
BitBucket has a built in issue tracker that makes it much easier to keep on top
of changes and bug fixes. In either case, please [create anew issue](https://bitbucket.org/sminez/arpy/issues/new)
using this link or by going to the `Issues` panel on the left of the page in
BitBucket. Please remember to tag Innes as the Assignee so that he gets notified
automatically!

(If it is a bug report, please copy in and error messages you get when describing
the problem.)

### Documentation
Please see the [docs](docs/) directory for markdown files detailing the use of the
module.
