bin_PROGRAMS = fmap get_dnu vlad tracy_2_3 carbon_map drv_terms get_nl_optics \
               opt_des_1_ptc opt_des_2_ptc

fmap_SOURCES          = fmap.cc
get_dnu_SOURCES       = get_dnu.cc
vlad_SOURCES          = vlad.cc
tracy_2_3_SOURCES     = tracy_2_3.cc
carbon_map_SOURCES    = carbon_map.cc
drv_terms_SOURCES     = drv_terms.cc
get_nl_optics_SOURCES = get_nl_optics.cc
opt_des_1_ptc_SOURCES = opt_des_1_ptc.cc
opt_des_2_ptc_SOURCES = opt_des_2_ptc.cc

LIBS     = -L$(TRACY_LIB)/tracy/lib -ltracy_ptc -lTPSALib -lLieLib \
           -L$(NUM_REC)/lib -L/usr/local/Cellar/gcc/9.3.0/lib/gcc/9 \
           -lnum_rec \
	   -lstdc++ \
           -lm \
           -lgfortran
#INCLUDES = -I$(TRACY_LIB)/tracy/inc -I$(NUM_REC)/inc
AM_CPPFLAGS = -I$(TRACY_LIB)/tracy/inc -I$(NUM_REC)/inc
AM_FFLAGS   = -g -O0 -Wall -fbounds-check
AM_CXXFLAGS = -g -O0 -Wall -Wno-non-template-friend -fPIC
#AM_CXXFLAGS = -g -O0 -Wall -Wno-non-template-friend -fPIC -Wl,--no-relax
