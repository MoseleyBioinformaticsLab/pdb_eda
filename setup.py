#!/usr/bin/env python3

import os
import sys
import re
from setuptools import setup, find_packages, Extension
import Cython.Build
import numpy

if sys.argv[-1] == 'publish':
    os.system('python3 setup.py sdist')
    os.system('twine upload dist/*')
    sys.exit()

def readme():
    with open('README.rst') as readme_file:
        return readme_file.read()

def find_version():
    with open('pdb_eda/__init__.py', 'r') as fd:
        version = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',fd.read(), re.MULTILINE).group(1)
    if not version:
        raise RuntimeError('Cannot find version information')
    return version

REQUIRES = [
    "docopt >= 0.6.2",
    "jsonpickle >= 1.4.1",
    "numpy >= 1.16.2",
    "scipy >= 1.2.1",
    "biopython >= 1.77",
    "PyCifRW >= 4.4.2",
    "cython >= 0.29.21"
]

SETUP_REQUIRES = [
    "numpy >= 1.14.5",
    "cython >= 0.29.21"
]

EXTENSIONS = [
    Extension("pdb_eda.cutils",sources=["pdb_eda/cutils.pyx"],include_dirs=[numpy.get_include()])
]

setup(
    name='pdb_eda',
    version=find_version(),
    author='Sen Yao, Hunter N.B. Moseley',
    author_email='hunter.moseley@gmail.com',
    description='Methods for analyzing electron density maps in wwPDB',
    keywords='PDB, electron densiy map',
    license='Modified Clear BSD License',
    url='https://github.com/MoseleyBioinformaticsLab/pdb_eda',
    cmdclass={'build_ext': Cython.Build.build_ext},
    packages=find_packages(),
    package_data={'pdb_eda': ['conf/*.json*']},
    platforms='any',
    long_description=readme(),
    setup_requires=SETUP_REQUIRES,
    install_requires=REQUIRES,
    ext_modules=EXTENSIONS,
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
    ],
    entry_points={"console_scripts": ["pdb_eda = pdb_eda.__main__:main"]},
)
