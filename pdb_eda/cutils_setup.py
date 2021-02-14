# Compile using:
# $ python3 cutils_setup.py build_ext --inplace
# $ cp build/lib.linux-x86_64-3.6/pdb_eda/cutils.cpython-36m-x86_64-linux-gnu.so cutils.so

from setuptools import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("cutils.pyx")
)
