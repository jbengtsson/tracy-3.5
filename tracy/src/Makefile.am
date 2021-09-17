
#lib_LIBRARIES          = libtracy.a libtracy_ptc.a
#lib_LIBRARIES          = libtracy.a
lib_LTLIBRARIES        = libtracy.la

#libtracy_a_SOURCES     = tracy_lib.cc
libtracy_la_SOURCES     = tracy_lib.cc

#libtracy_ptc_a_SOURCES = tracy_ptc_lib.cc

# AM_LDFLAGS flags:
#  {-static | -shared}
#  {-no-undefined | -undefined dynamic_lookup}
AM_LDFLAGS  = -shared

LIBS	    = -lstdc++ \
              -lgsl \
              -larmadillo 
              $(GSL_LIBS)

## LIBS	    = -L/usr/local/Cellar/gcc/11.1.0_1/lib/gcc/11 \
#              -lstdc++ \
#              -L/usr/local/Cellar/gsl/2.7/lib \
#              -lgsl \
#              -L/usr/local/Cellar/armadillo/10.6.1/lib \
#              -larmadillo 
##              $(GSL_LIBS)

# # C Preprocessor flags.
# AM_CPPFLAGS = -I../inc \
#               -I/usr/local/Cellar/armadillo/10.6.1/include \
#               -I/usr/local/Cellar/gsl/2.7/include \
#               -I/usr/local/Cellar/pybind11/2.6.2/include \
#               -I/usr/local/Cellar/python@3.9/3.9.6/Frameworks/Python.framework/Versions/Current/include/python3.9 \
#               -I/usr/local/carma/include \
#               -I/usr/local/Cellar/numpy/1.21.1/lib/python3.9/site-packages/numpy/core/include
# 
#              $(GSL_CFLAGS)

AM_CPPFLAGS = -I../inc \
              $(GSL_CFLAGS)

# C++ flags.
AM_CXXFLAGS = -g -O0 -Wall -Wno-error=all -std=gnu++14 -fPIC \
              -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC
#              $(GSL_CFLAGS)

# Fortan flags.
FFLAGS      = -g -O2 -Wall -fbounds-check -fdefault-integer-8 -mcmodel=large

# AM_CXXFLAGS flags:
#  {-g}                      debugging
#  {-pg}                     profiling
#  {-O0}                     optimization
#  {-Wall}                   warnings
#  {-Wno-error=all}          errors
#  {-fno-implicit-templates}
#  {-fPIC}
#  {-frepo}
