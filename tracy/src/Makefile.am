AUTOMAKE_OPTIONS = subdir-objects

lib_LIBRARIES = libtracy.a libtracy_ptc.a

libtracy_a_SOURCES = tracy_lib.cc

libtracy_ptc_a_SOURCES = \
	tracy_ptc_lib.cc \
	../../TPSA/TPSALib.f ../../TPSA/LieLib.f

GNU = /usr/local/opt

# C Preprocessor flags.
AM_CPPFLAGS = \
	-I$(GNU)/gcc/include \
	-I$(GNU)/gsl/include \
	-I../inc \
	-I$(NUM_REC)/inc

# FORTRAN 77 flags.
FFLAGS = \
	-g -O2 -Wall -fbounds-check -fdefault-integer-8 -mcmodel=medium \
	-w -fallow-argument-mismatch

# C++ flags.
AM_CXXFLAGS = -g -O0 -Wall -Wno-error=all -fPIC
