#!/usr/bin/env python3

import os
import sys
import re
from setuptools import setup, find_packages


if sys.argv[-1] == 'publish':
    os.system('python3 setup.py sdist')
    os.system('twine upload dist/*')
    sys.exit()


def readme():
    with open('README.rst') as readme_file:
        return readme_file.read()


def find_version():
    with open('pdb_eda/__init__.py', 'r') as fd:
        version = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                            fd.read(), re.MULTILINE).group(1)
    if not version:
        raise RuntimeError('Cannot find version information')
    return version


REQUIRES = [
    "docopt",
    "jsonpickle",
    "numpy",
    "scipy",
    "pandas",
    "biopython"
]


setup(
    name='pdb_eda',
    version=find_version(),
    author='Sen Yao',
    author_email='yaosen1313@gmail.com',
    description='Methods for analyzing electron density maps in wwPDB',
    keywords='PDB, electron densiy map',
    license='Modified Clear BSD License',
    url='https://github.com/MoseleyBioinformaticsLab/pdb_eda',
    packages=find_packages(),
    package_data={'pdb_eda': ['conf/*.json']},
    platforms='any',
    long_description=readme(),
    install_requires=REQUIRES,
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ]
)
