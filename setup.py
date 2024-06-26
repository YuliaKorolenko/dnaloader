from distutils.core import setup
from Cython.Build import cythonize
import numpy
import torch

setup(ext_modules=cythonize(["dnaloader/sequences/toonehot.pyx", "dnaloader/characteristics/bwhelper.pyx"]),
      include_dirs=[numpy.get_include()]
      )
