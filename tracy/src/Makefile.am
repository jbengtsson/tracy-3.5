lib_LIBRARIES          = libtracy.a libtracy_ptc.a

libtracy_a_SOURCES     = tracy_lib.cc

libtracy_ptc_a_SOURCES = tracy_ptc_lib.cc

#INCLUDES = -I../inc -I$(NUM_REC)/inc
# GSL_INCLUDE = -I/usr/local/Cellar/gsl/2.7.1/include
GSL_INCLUDE = -I/net/home/lli/OPAL/OPAL-2.0.1/include
AM_CPPFLAGS = $(GSL_INCLUDE) \
              -I../inc -I$(NUM_REC)/inc
#INCLUDES = -I../inc -I$(NUM_REC)/inc \
#	   -I/usr/include/c++/3.3.6

FFLAGS      = -g -O2 -Wall -fbounds-check -fdefault-integer-8 -mcmodel=large
AM_CXXFLAGS = -g -O0 -Wall -Wno-error=all -fPIC
# AM_CXXFLAGS = -g -O0 -Wall -Wno-error=all -fno-implicit-templates -fPIC
#CXXFLAGS = -pg -O2 -Wall -fno-implicit-templates -fPIC
#CXXFLAGS = -g -O2 -Wall -frepo
