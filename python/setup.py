import logging
import os
import subprocess
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
# inspired by https://github.com/mattpitkin/lintegrate
def gsl_config(*args, **kwargs):
    """Run gsl-config and return pre-formatted output
    """
    gsl_config = 'gsl-config'

    if os.name == 'nt':
        cmd = "{} {}".format([gsl_config, " ".join(args)])
        kwargs.setdefault("shell", True)

    else:
        cmd = [gsl_config] + list(args)

    output = subprocess.check_output(cmd, **kwargs)
    return output.decode("utf-8").strip()


# execute gsl-config --cflags to get include dir
gsl_include = gsl_config("--cflags")[2:]


# execute gsl-config --libs for dir and libs
# assuming posix compatabilitly
_gsl_lib = gsl_config("--libs").split(" ")
# assumption :
gsl_lib_dir = _gsl_lib[0][2:]
gsl_libs = [ll[2:] for ll in _gsl_lib if ll[:2] == '-l']
del _gsl_lib

logger = logging.getLogger('setup')
logger.info('Using GSL: include %s library_dir %s libraries %s',
            gsl_include, gsl_lib_dir, gsl_libs)

ext_modules = [
    Pybind11Extension(
        "libtracy",
        sorted(['src/tracy_py.cc']),
        include_dirs=['../tracy/inc'] + [gsl_include],
        library_dirs=['../tracy/src/.libs'] + [gsl_lib_dir],
        libraries=['tracy'] + gsl_libs
    ),
]

setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=ext_modules
)
