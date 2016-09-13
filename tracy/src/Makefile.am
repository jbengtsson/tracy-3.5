lib_LIBRARIES          = libtracy.a libtracy_ptc.a

libtracy_a_SOURCES     = tracy_lib.cc

libtracy_ptc_a_SOURCES = tracy_ptc_lib.cc

#INCLUDES = -I../inc -I$(NUM_REC)/inc
AM_CPPFLAGS = -I../inc -I$(NUM_REC)/inc
#INCLUDES = -I../inc -I$(NUM_REC)/inc \
#	   -I/usr/include/c++/3.3.6

AM_CXXFLAGS    = -g -O0 -Wall -Wno-error=all -fno-implicit-templates -fPIC
#CXXFLAGS   = -pg -O2 -Wall -fno-implicit-templates -fPIC
#CXXFLAGS   = -g -O2 -Wall -frepo
