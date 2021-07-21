import os
from setuptools import setup
from distutils.file_util import copy_file

from pybind11.setup_helpers import Pybind11Extension, build_ext

# Make a copy of README.rst ... required by setup.cfg
t_dir = os.path.dirname(__file__)
par_dir = os.path.normpath(os.path.join(t_dir, os.pardir))
readme_name = 'README.rst'
copy_file(
    os.path.join(par_dir, readme_name),
    os.path.join(t_dir, readme_name),
    update=True
)

# execute gsl-config --libs to see what is required
gsl_lib_dirs = ['/usr/lib/x86_64-linux-gnu']
gsl_libs = ['gsl', 'gslcblas', 'm']
gsl_include = ['']

ext_modules = [
    Pybind11Extension(
        "libtracy",
        sorted(['src/tracy_py.cc']),
        include_dirs=['../tracy/inc'] + gsl_include,
        library_dirs=['../tracy/src/.libs'] + gsl_lib_dirs,
        libraries=['tracy'] + gsl_libs
    ),
]

setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=ext_modules
)
