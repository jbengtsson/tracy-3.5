bin_PROGRAMS = fmap get_dnu vlad tracy_2_3 carbon_map drv_terms get_nl_optics \
               opt_des_1_ptc opt_des_2_ptc poincare_map wakefield map_norm \
               lielib map_norm_jb stat_mom moments tune_conf bessy_ii max_4u \
               max_4u_scan_dnu

fmap_SOURCES            = fmap.cc
get_dnu_SOURCES         = get_dnu.cc
vlad_SOURCES            = vlad.cc
tracy_2_3_SOURCES       = tracy_2_3.cc
carbon_map_SOURCES      = carbon_map.cc
drv_terms_SOURCES       = drv_terms.cc
get_nl_optics_SOURCES   = get_nl_optics.cc
opt_des_1_ptc_SOURCES   = opt_des_1_ptc.cc
opt_des_2_ptc_SOURCES   = opt_des_2_ptc.cc
poincare_map_SOURCES    = poincare_map.cc
wakefield_SOURCES       = wakefield.cc
map_norm_SOURCES        = map_norm.cc
lielib_SOURCES          = lielib.cc
map_norm_jb_SOURCES     = map_norm_jb.cc
stat_mom_SOURCES        = stat_mom.cc
moments_SOURCES         = moments.cc
tune_conf_SOURCES       = tune_conf.cc
bessy_ii_SOURCES        = bessy_ii.cc
max_4u_SOURCES          = max_4u.cc
max_4u_scan_dnu_SOURCES = max_4u_scan_dnu.cc

GSL_INCLUDE = /usr/local/Cellar/gsl/2.7.1
# GSL_INCLUDE = /net/home/lli/OPAL/OPAL-2.0.1

LIBS        = -L$(TRACY_LIB)/tracy/lib \
              -L$(TRACY_LIB)/tracy/lib64 \
              -ltracy_ptc -lTPSALib -lLieLib \
              -L$(NUM_REC)/lib -lnum_rec \
              -lstdc++ \
              -lm \
              -L$(GSL_INCLUDE)/lib -lgsl -lgslcblas \
              -L/usr/local/Cellar/gcc/14.1.0_1/lib/gcc/current \
              -lgfortran

AM_CPPFLAGS = -I/opt/OPAL-2.2.1/lib \
              -I$(TRACY_LIB)/tracy/inc -I$(NUM_REC)/inc \
              -I$(GSL_INCLUDE)/include \
              -std=c++11

AM_FFLAGS   = -g -O0 -Wall -fbounds-check

AM_CXXFLAGS = -g -O0 -Wall -fPIC
#AM_CXXFLAGS = -g -O0 -Wall -Wno-non-template-friend -fPIC -Wl,--no-relax
