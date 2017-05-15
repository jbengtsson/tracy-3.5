bin_PROGRAMS = fmap get_dnu vlad tracy_2_3

fmap_SOURCES      = fmap.cc
get_dnu_SOURCES   = get_dnu.cc
vlad_SOURCES      = vlad.cc
tracy_2_3_SOURCES = tracy_2_3.cc

LIBS     = -L$(TRACY_LIB)/tracy/lib -ltracy_ptc -lTPSALib -lLieLib \
           -L$(NUM_REC)/lib -lnum_rec \
	   -lstdc++ \
           -lm \
           -lgfortran
#INCLUDES = -I$(TRACY_LIB)/tracy/inc -I$(NUM_REC)/inc
AM_CPPFLAGS = -I$(TRACY_LIB)/tracy/inc -I$(NUM_REC)/inc
AM_FFLAGS   = -g -O0 -Wall -fbounds-check
AM_CXXFLAGS = -g -O0 -Wall -Wno-non-template-friend -fPIC
