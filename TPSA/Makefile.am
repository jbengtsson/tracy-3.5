AUTOMAKE_OPTIONS = subdir-objects

lib_LIBRARIES = libTPSALib.a libLieLib.a
# Make shared libraries - use glibtoolize to initialise.
# lib_LTLIBRARIES = libTPSALib.la libLieLib.la

libTPSALib_a_SOURCES = TPSALib.f

libLieLib_a_SOURCES  = LieLib.f

# FORTRAN 77 flags.
# Use -O2 to avoid Clang issues on MacBook.
#FFLAGS = \
#	-ggdb -O2 -Wall -fbounds-check -std=gnu -fdefault-integer-8 \
#	-mcmodel=medium

# FORTRAN 77 flags (The package maintainer’s variant of FFLAGS).
AM_FFLAGS = -ggdb -O -fbacktrace -Wall -fbounds-check -fPIC \
            -std=gnu -fdefault-integer-8 $(MCMODEL_FFLAGS)
