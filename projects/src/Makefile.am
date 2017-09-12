SUBDIRS      = ptc

bin_PROGRAMS = optics dynap leac dnu touschek ibs vlad tbt \
               max_iv_1 max_iv_2 ibs_for_les delta bare_da match \
               carbon_map dls 6BA

#sc_f_LINK    = gfortran -g -O0 -Wall -fbounds-check -o sc_f

#main_SOURCES      = main.cc
optics_SOURCES      = optics.cc
dynap_SOURCES       = dynap.cc
dnu_SOURCES         = dnu.cc
leac_SOURCES        = leac.cc
touschek_SOURCES    = touschek.cc
ibs_SOURCES         = ibs.cc
vlad_SOURCES        = vlad.cc
tbt_SOURCES         = tbt.cc
#m4_1_SOURCES      = m4_1.cc
max_iv_1_SOURCES    = max_iv_1.cc
max_iv_2_SOURCES    = max_iv_2.cc
ibs_for_les_SOURCES = ibs_for_les.cc
delta_SOURCES       = delta.cc
bare_da_SOURCES     = bare_da.cc
match_SOURCES       = match.cc
carbon_map_SOURCES  = carbon_map.cc
dls_SOURCES         = dls.cc
6BA_SOURCES         = 6BA.cc

#LIBS     = -L$(TRACY_LIB)/tracy -ltracy
# gfortran needs to be before stdc_++ for Fortran programs
LIBS	 = -L$(TRACY_LIB)/tracy/lib -ltracy \
	   -L$(NUM_REC)/lib -lnum_rec \
           -lstdc++ \
           -lm
#	   -lgfortran
#INCLUDES = -I$(TRACY_LIB)/tracy/inc -I$(NUM_REC)/inc
AM_CPPFLAGS = -I$(TRACY_LIB)/tracy/inc -I$(NUM_REC)/inc
AM_FFLAGS   = -g -O0 -Wall -fbounds-check
AM_CXXFLAGS = -g -O0 -Wall -Wno-non-template-friend
#CXXFLAGS = -pg -O2 -Wall -Wno-non-template-friend
