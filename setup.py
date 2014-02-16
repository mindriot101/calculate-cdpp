from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from numpy import get_include

exts = Extension('calculate_cdpp', ['cython/calculate_cdpp.pyx',
    'src/calculate_cdpp.c'],
                libraries=['cfitsio'],
                include_dirs=[get_include(), 'include'],
                )

setup(
        name='calculate_cdpp',
        author='Simon Walker',
        author_email='s.r.walker101@googlemail.com',
        version='0.0.1',
        ext_modules = cythonize(exts),
        )
