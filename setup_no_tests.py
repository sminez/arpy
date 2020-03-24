import os

from setuptools import find_packages, setup

from arpy import __version__

# allow setup.py to be run from any path
os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

setup(
    name="arpy",
    version=__version__,
    description="A library for computing with Absolute Relativity",
    url="https://github.com/sminez/arpy",
    author="Innes Anderson-Morrison",
    author_email="innesdmorrison@gmail.com",
    install_requires=[],
    packages=find_packages(),
    package_dir={"arpy": "arpy"},
    classifiers=["Programming Language :: Python :: 3", "Development Status :: 4 - Beta"],
)
