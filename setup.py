from setuptools import setup, find_packages
import os

# allow setup.py to be run from any path
os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

setup(
    name='arpy',
    version="0.0.1",
    description="A REPL for computing with Absolute Relativity",
    url="https://bitbucket.com/sminez/arpy",
    author="Innes Anderson-Morrison",
    author_email='innesdmorrison@gmail.com',
    install_requires=[
        'sly',
        'prompt_toolkit==1.0.0',
    ],
    tests_require=[
        'pytest',
        ],
    packages=find_packages(),
    package_dir={'arpy': 'arpy'},
    zip_safe=False,
    classifiers=[
        'Programming Language :: Python',
        'Development Status :: 4 - Beta'
    ],
    entry_points={
        'console_scripts': [
            'arpy = arpy.cli:run',
        ]
    },
)