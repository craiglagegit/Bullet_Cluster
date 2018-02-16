# Run as:
#    python setup.py build_ext --inplace

import numpy
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
# from Cython.Build import cythonize

ext_modules = [Extension(
    name="fill",
    sources=["fill.pyx"],
    include_dirs = [numpy.get_include()],  # .../site-packages/numpy/core/include
    language="c",
    )]

setup(
    name = 'fill',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules, )
