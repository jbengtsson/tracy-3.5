#lib_LIBRARIES          = libtracy.a libtracy_ptc.a
lib_LIBRARIES          = libtracy.a

libtracy_a_SOURCES     = tracy_lib.cc

libtracy_ptc_a_SOURCES = tracy_ptc_lib.cc

# C Preprocessor flags.
AM_CPPFLAGS = -I../inc \
              -I/usr/local/Cellar/armadillo/10.5.2/include \
              -I/usr/local/Cellar/gsl/2.7/include

# Fortan flags.
FFLAGS      = -g -O2 -Wall -fbounds-check -fdefault-integer-8 -mcmodel=large
# C++ flags.
AM_CXXFLAGS = -g -O0 -Wall -Wno-error=all -std=c++11 -fPIC
# AM_CXXFLAGS = -g -O0 -Wall -Wno-error=all -fno-implicit-templates -fPIC
#CXXFLAGS = -pg -O2 -Wall -fno-implicit-templates -fPIC
#CXXFLAGS = -g -O2 -Wall -frepo
