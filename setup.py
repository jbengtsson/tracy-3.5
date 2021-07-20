from glob import glob
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

# execute gsl-config --libs to see what is required
gsl_lib_dirs = ['/usr/lib/x86_64-linux-gnu']
gsl_libs = ['gsl', 'gslcblas', 'm']
ext_modules = [
    Pybind11Extension(
        "libtracy",
        sorted(["tracy/src/tracy_py.cc"]),
        include_dirs=['tracy/inc'],
        library_dirs=['tracy/src/.libs'] + gsl_lib_dirs,
        libraries=['tracy'] + gsl_libs
    ),
]

setup(
    cmdclass={"build_ext": build_ext},
    ext_modules=ext_modules
)
