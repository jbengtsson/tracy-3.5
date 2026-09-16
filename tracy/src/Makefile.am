AUTOMAKE_OPTIONS = subdir-objects

# lib_LIBRARIES = libtracy.a libtracy_ptc.a
# Make shared libraries - use glibtoolize to initialise.
lib_LTLIBRARIES = libtracy.la libtracy_ptc.la

libtracy_la_SOURCES = tracy_lib.cc

libtracy_ptc_la_SOURCES = tracy_ptc_lib.cc

# C Preprocessor flags.
# GSL_CFLAGS is resolved by configure via gsl-config.
AM_CPPFLAGS = \
	$(GSL_CFLAGS) \
	-I../inc \
	-I$(NUM_REC)/inc

# FORTRAN 77 flags.
#FFLAGS = \
#	-g -O2 -Wall -fbounds-check -fdefault-integer-8 -mcmodel=medium \
#	-w -fallow-argument-mismatch

# C++ flags.
#AM_CXXFLAGS = -g -O0 -Wall -Wno-error=all
AM_CXXFLAGS = \
	-g -O0 -pedantic -Wall -std=c++11 \
	-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -fPIC

AM_FFLAGS = -ggdb -O -fbacktrace -Wall -fbounds-check -fPIC \
            -std=gnu -fdefault-integer-8 $(MCMODEL_FFLAGS)

