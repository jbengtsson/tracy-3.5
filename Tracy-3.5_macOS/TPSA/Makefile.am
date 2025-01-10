AUTOMAKE_OPTIONS = subdir-objects

lib_LIBRARIES = libTPSALib.a libLieLib.a

libTPSALib_a_SOURCES = TPSALib.f

libLieLib_a_SOURCES  = LieLib.f

# FORTRAN 77 flags.
# Use -O2 to avoid Clang issued on MacBook.
FFLAGS = \
	-ggdb -O2 -Wall -fbounds-check -std=gnu -fdefault-integer-8 \
	-mcmodel=medium -w -fallow-argument-mismatch
#FFLAGS	= -ggdb -O2 -Wall -fbounds-check -fdefault-integer-8 -mcmodel=large

# FORTRAN 77 flags (The package maintainer’s variant of FFLAGS).
AM_FFLAGS = -ggdb -O2 -Wall -fbounds-check -fPIC
