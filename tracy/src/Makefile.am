
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

LIBS	    = -L/usr/local/Cellar/gcc/11.1.0_1/lib/gcc/11 \
              -lstdc++ \
              -L/usr/local/Cellar/armadillo/10.5.2/lib \
              -larmadillo \
              $(GSL_LIBS)

# C Preprocessor flags.
AM_CPPFLAGS = -I../inc \
              -I/usr/local/Cellar/armadillo/10.5.2/include \
              $(GSL_CFLAGS)

# C++ flags.
AM_CXXFLAGS = -g -O0 -Wall -Wno-error=all -std=gnu++11 -fPIC \
              $(GSL_CFLAGS)

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
