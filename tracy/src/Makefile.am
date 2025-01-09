AUTOMAKE_OPTIONS = subdir-objects

# lib_LIBRARIES = libtracy.a libtracy_ptc.a

# Make shared libraries - use glibtoolize to initialise.
lib_LTLIBRARIES = libtracy.la libtracy_ptc.la

libtracy_la_SOURCES = tracy_lib.cc

libtracy_ptc_la_SOURCES = \
	tracy_ptc_lib.cc

GNU = /usr/local/opt

# C Preprocessor flags.
AM_CPPFLAGS = \
	-I$(GNU)/gcc/include \
	-I$(GNU)/gsl/include \
	-I../inc \
	-I$(NUM_REC)/inc

# Load flags.
LDFLAGS = $(THOR_LIB)/TPSA/TPSALib.o $(THOR_LIB)/TPSA/LieLib.o

# FORTRAN 77 flags.
FFLAGS = \
	-g -O2 -Wall -fbounds-check -fdefault-integer-8 -mcmodel=medium \
	-w -fallow-argument-mismatch

# C++ flags.
AM_CXXFLAGS = -g -O0 -Wall -Wno-error=all

AM_LDFLAGS = -F/System/Library/Frameworks -framework CoreServices
