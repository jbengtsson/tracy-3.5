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

LIBS     = -L$(TRACY_LIB)/tracy/lib
           -L$(TRACY_LIB)/tracy/lib64
           -ltracy_ptc -lTPSALib -lLieLib \
           -L$(NUM_REC)/lib -lnum_rec \
           -L/usr/local/Cellar/gsl/2.7.1/lib \
           -L/net/home/lli/OPAL/OPAL-2.0.1/include \
           -L/opt/OPAL-2.2.1/lib \
           -lstdc++ \
           -lm \
           -lgsl -lgslcblas \
           -lgfortran
AM_CPPFLAGS = -I/opt/OPAL-2.2.1/lib \
              -I$(TRACY_LIB)/tracy/inc -I$(NUM_REC)/inc
AM_FFLAGS   = -g -O0 -Wall -fbounds-check
AM_CXXFLAGS = -g -O0 -Wall -fPIC
#AM_CXXFLAGS = -g -O0 -Wall -Wno-non-template-friend -fPIC -Wl,--no-relax
