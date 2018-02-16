# Run as:
#    python setup.py build_ext --inplace

NAME = "fill"
EXT_SOURCES = []
EXT_LIBRARIES=[]
EXT_LIBRARY_DIRS=[]
EXT_INCLUDE_DIRS=['/export/ursa1/csl336/Software/yt-x86_64/lib/python2.7/site-packages/numpy/core/include']
DEFINES = []

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension(NAME,
                         [NAME+".pyx"] + EXT_SOURCES,
                         libraries = EXT_LIBRARIES,
                         library_dirs = EXT_LIBRARY_DIRS,
                         include_dirs = EXT_INCLUDE_DIRS,
                         define_macrso = DEFINES)
]

setup(
  name = NAME,
  cmdclass = {'build_ext':build_ext},
  ext_modules = ext_modules
)
