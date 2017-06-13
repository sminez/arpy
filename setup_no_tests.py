from setuptools import setup, find_packages
import os

# allow setup.py to be run from any path
os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

setup(
    name='arpy',
    version="0.1.7",
    description="A library for computing with Absolute Relativity",
    url="https://github.com/sminez/arpy",
    author="Innes Anderson-Morrison",
    author_email='innesdmorrison@gmail.com',
    install_requires=[],
    packages=find_packages(),
    package_dir={'arpy': 'arpy'},
    classifiers=[
        'Programming Language :: Python :: 3',
        'Development Status :: 4 - Beta'
    ]
)