AUTOMAKE_OPTIONS = subdir-objects

check_PROGRAMS = tracking_element_test
TESTS = $(check_PROGRAMS)
EXTRA_DIST = README.md

tracking_element_test_SOURCES = \
	tracking_element_test.cc \
	tracking_support.cc \
	tracking_support.h \
	tracking_elements.cc
tracking_element_test_LDADD = \
	$(top_builddir)/tracy/src/libtracy.la \
	-L$(NUM_REC)/lib -lnum_rec \
	$(GSL_LIBS)

AM_CPPFLAGS = \
	$(GSL_CFLAGS) \
	-I$(top_srcdir)/tracy/inc \
	-I$(NUM_REC)/inc

AM_CXXFLAGS = \
	-g -O0 -pedantic -Wall -std=c++11 \
	-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC
