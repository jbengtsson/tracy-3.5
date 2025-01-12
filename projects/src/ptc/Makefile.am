AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = \
  fmap get_dnu vlad tracy_2_3 carbon_map drv_terms get_nl_optics opt_des_1_ptc \
  opt_des_2_ptc poincare_map wakefield lielib map_norm_jb \
  stat_mom tune_conf bessy_ii max_4u lie_series_chk

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
# map_norm_SOURCES        = map_norm.cc
lielib_SOURCES          = lielib.cc
map_norm_jb_SOURCES     = map_norm_jb.cc
stat_mom_SOURCES        = stat_mom.cc
# moments_SOURCES         = moments.cc
tune_conf_SOURCES       = tune_conf.cc
bessy_ii_SOURCES        = bessy_ii.cc
max_4u_SOURCES          = max_4u.cc
# max_4u_scan_dnu_SOURCES = max_4u_scan_dnu.cc
lie_series_chk_SOURCES  = lie_series_chk.cc

LOCAL = /usr/local/opt

# Explicit reference to libc++ required for CC front-end vs. C++.
# Cope with macOS by using static libraries.
LIBS =	$(TRACY_LIB)/tracy/lib/libtracy_ptc.a \
	$(THOR_LIB)/thor/lib/libTPSALib.a $(THOR_LIB)/thor/lib/libLieLib.a \
	-L$(NUM_REC)/lib -lnum_rec \
	-L$(LOCAL)/gsl/lib -lgsl \
	-L$(LOCAL)/gfortran/lib/gcc/current -lgfortran

# Load flags.
# LDFLAGS = $(THOR_LIB)/TPSA/TPSALib.o $(THOR_LIB)/TPSA/LieLib.o

# C Preprocessor flags.
AM_CPPFLAGS = \
	-I$(TRACY_LIB)/tracy/inc \
	-I$(NUM_REC)/inc

# FORTRAN 77 flags (The package maintainer’s variant of FFLAGS).
AM_FFLAGS = -ggdb -O0 -Wall -fbounds-check

# C++ flags.
AM_CXXFLAGS = \
	-g -O0 -pedantic -Wall \
	-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -fPIC 
# AM_CXXFLAGS = -ggdb -O0 -Wall -std=c++11
#AM_CXXFLAGS = -ggdb -O0 -Wall -Wno-non-template-friend -fPIC -Wl,--no-relax

AM_LDFLAGS = -F/System/Library/Frameworks -framework CoreServices
