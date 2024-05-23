#SUBDIRS      = ptc

bin_PROGRAMS = optics dynap leac dnu touschek ibs vlad tbt tbt_ref \
               max_iv_1 max_iv_2 ibs_for_les delta bare_da match \
               carbon_map dls 6BA track_fft m4 m4_mod dls_ii track opt_des_1 \
	       opt_des_2 id_model sigma tech_note dx_dj tune_scan sxt \
               chk_drv_terms fit_nu opt_des_bessy_iii for_pierre tst_it \
               ctrl_cod wakefield compute_maps lin_opt map_diag_3D bessy_ii \
               rev_bend max_4u

#sc_f_LINK    = gfortran -g -O0 -Wall -fbounds-check -o sc_f

#main_SOURCES              = main.cc
optics_SOURCES            = optics.cc
dynap_SOURCES             = dynap.cc
dnu_SOURCES               = dnu.cc
leac_SOURCES              = leac.cc
touschek_SOURCES          = touschek.cc
ibs_SOURCES               = ibs.cc
vlad_SOURCES              = vlad.cc
tbt_SOURCES               = tbt.cc
tbt_ref_SOURCES           = tbt_ref.cc
max_iv_1_SOURCES          = max_iv_1.cc
max_iv_2_SOURCES          = max_iv_2.cc
ibs_for_les_SOURCES       = ibs_for_les.cc
delta_SOURCES             = delta.cc
bare_da_SOURCES           = bare_da.cc
match_SOURCES             = match.cc
carbon_map_SOURCES        = carbon_map.cc
dls_SOURCES               = dls.cc
6BA_SOURCES               = 6BA.cc
track_fft_SOURCES         = track_fft.cc
m4_SOURCES                = m4.cc
m4_mod_SOURCES            = m4_mod.cc
dls_ii_SOURCES            = dls_ii.cc
track_SOURCES             = track.cc
opt_des_1_SOURCES         = opt_des_1.cc
#opt_des_2_SOURCES         = opt_des_2.cc Powell/src/newuoa.f
opt_des_2_SOURCES         = opt_des_2.cc
id_model_SOURCES          = id_model.cc
sigma_SOURCES             = sigma.cc
tech_note_SOURCES         = tech_note.cc
dx_dj_SOURCES             = dx_dj.cc
tune_scan_SOURCES         = tune_scan.cc
sxt_SOURCES               = sxt.cc
chk_drv_terms_SOURCES     = chk_drv_terms.cc
fit_nu_SOURCES            = fit_nu.cc
opt_des_bessy_iii_SOURCES = opt_des_bessy_iii.cc
for_pierre_SOURCES        = for_pierre.cc
tst_it_SOURCES            = tst_it.cc
ctrl_cod_SOURCES          = ctrl_cod.cc
wakefield_SOURCES         = wakefield.cc
compute_maps_SOURCES      = compute_maps.cc
lin_opt_SOURCES           = lin_opt.cc
map_diag_3D_SOURCES       = map_diag_3D.cc
bessy_ii_SOURCES          = bessy_ii.cc
rev_bend_SOURCES          = rev_bend.cc
max_4u_SOURCES            = max_4u.cc

GSL_INCLUDE = /usr/local/Cellar/gsl/2.7.1
# GSL_INCLUDE = /net/home/lli/OPAL/OPAL-2.0.1

LIBS	    = -L$(TRACY_LIB)/tracy/lib \
              -L$(TRACY_LIB)/tracy/lib64 \
              -ltracy \
              -L$(NUM_REC)/lib -lnum_rec \
              -lstdc++ \
              -L$(GSL_INCLUDE)/lib -lgsl -lgslcblas

# C Preprocessor flags.
AM_CPPFLAGS = -I/usr/local/Cellar/gsl/2.7.1/include \
              -I$(TRACY_LIB)/tracy/inc \
              -I$(NUM_REC)/inc \
              -I$(GSL_INCLUDE)/include

#              -I/opt/OPAL-2.2.1/lib 

# C++ flags.
AM_CXXFLAGS = -g -O0 -pedantic -Wall -std=c++11 \
              -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC
#CXXFLAGS = -pg -O2 -Wall -Wno-non-template-friend
