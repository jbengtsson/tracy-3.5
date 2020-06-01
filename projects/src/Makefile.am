SUBDIRS      = ptc

bin_PROGRAMS = optics dynap leac dnu touschek ibs vlad tbt \
               max_iv_1 max_iv_2 ibs_for_les delta bare_da match \
               carbon_map dls 6BA track_fft m4 m4_mod dls_ii track opt_des_1 \
	       opt_des_2 id_model sigma tech_note dx_dj tune_scan sxt \
               chk_drv_terms

#sc_f_LINK    = gfortran -g -O0 -Wall -fbounds-check -o sc_f

#main_SOURCES        = main.cc
optics_SOURCES        = optics.cc
dynap_SOURCES         = dynap.cc
dnu_SOURCES           = dnu.cc
leac_SOURCES          = leac.cc
touschek_SOURCES      = touschek.cc
ibs_SOURCES           = ibs.cc
vlad_SOURCES          = vlad.cc
tbt_SOURCES           = tbt.cc
max_iv_1_SOURCES      = max_iv_1.cc
max_iv_2_SOURCES      = max_iv_2.cc
ibs_for_les_SOURCES   = ibs_for_les.cc
delta_SOURCES         = delta.cc
bare_da_SOURCES       = bare_da.cc
match_SOURCES         = match.cc
carbon_map_SOURCES    = carbon_map.cc
dls_SOURCES           = dls.cc
6BA_SOURCES           = 6BA.cc
track_fft_SOURCES     = track_fft.cc
m4_SOURCES            = m4.cc
m4_mod_SOURCES        = m4_mod.cc
dls_ii_SOURCES        = dls_ii.cc
track_SOURCES         = track.cc
opt_des_1_SOURCES     = opt_des_1.cc
#opt_des_2_SOURCES     = opt_des_2.cc Powell/src/newuoa.f
opt_des_2_SOURCES     = opt_des_2.cc
id_model_SOURCES      = id_model.cc
sigma_SOURCES         = sigma.cc
tech_note_SOURCES     = tech_note.cc
dx_dj_SOURCES         = dx_dj.cc
tune_scan_SOURCES     = tune_scan.cc
sxt_SOURCES           = sxt.cc
chk_drv_terms_SOURCES = chk_drv_terms.cc

#opt_des_2_LDADD = -L/usr/lib/x86_64-linux-gnu \
#	  	  -lgfortran

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
AM_CXXFLAGS = -g -O0 -Wall
#CXXFLAGS = -pg -O2 -Wall -Wno-non-template-friend
FFLAGS      = -g -O2 -Wall -fbounds-check -fdefault-integer-8
