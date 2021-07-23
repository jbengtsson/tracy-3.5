'''GSL config: read export to python
'''
import os
import subprocess
import logging

logger = logging.getLogger('setup')


# execute gsl-config --libs to see what is required
# inspired by https://github.com/mattpitkin/lintegrate
def _gsl_config(*args, **kwargs):
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


def gsl_config():
    '''provide a dictionary with GSL dirs and libraries
    '''
    # execute gsl-config --cflags to get include dir
    gsl_include = _gsl_config("--cflags")[2:]

    # execute gsl-config --libs for dir and libs
    # assuming posix compatabilitly
    gsl_lib = _gsl_config("--libs").split(" ")

    # assumption :
    gsl_lib_dir = gsl_lib[0][2:]
    gsl_libs = [ll[2:] for ll in gsl_lib if ll[:2] == '-l']

    logger = logging.getLogger('setup')
    logger.info('Using GSL: include %s library_dir %s libraries %s',
                gsl_include, gsl_lib_dir, gsl_libs)

    return dict(gsl_include=gsl_include, gsl_lib_dir=gsl_lib_dir,
                gsl_libs=gsl_libs)


if __name__ == '__main__':
    print(gsl_config())
