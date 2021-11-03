# Compile using:
# $ python3 cutils_setup.py build_ext --inplace

from setuptools import setup, Extension
import Cython.Build

setup(
  name = 'cutils',
  ext_modules=[
    Extension('cutils',
              sources=['cutils.pyx'],
              extra_compile_args=['-O3'],
              language='c++')
    ],
  cmdclass = {'build_ext': Cython.Build.build_ext}
)
