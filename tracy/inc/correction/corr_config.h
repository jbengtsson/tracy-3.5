#ifndef CORRECTION_CORR_CONFIG_H
#define CORRECTION_CORR_CONFIG_H

// Correction configuration and reference state — first slice of the
// param_data_type config extraction. See correction_refactor.md.
//
// Two things live here and are deliberately kept apart:
//
//   orbit_cfg   — knobs read from param.dat. Input.
//   bare_optics — the design-lattice Twiss at the sextupoles, MEASURED from the
//                 error-free machine by capture(). Not input: it is computed
//                 from the lattice, and must be re-captured whenever the design
//                 optics change (e.g. after a tune/chroma fit). It is the
//                 reference every beta-beat / tune-shift number is quoted
//                 against.
//
// Lumping both into one "config" bag would be a category error — the refactor
// note's original sketch called for a single config_data; splitting it keeps the
// measured reference from masquerading as a parsed knob.

// Max no of quad corrector families listable on param.dat's ID_quads line — a
// config sizing limit, so it lives here (it also sizes id_corr::b2). Global
// scope, not namespaced, to match the sibling limits in id_corr.h.
const int N_Fam_max = 25;

namespace corr {

struct coupling_cfg;  // correction/loco/coupling_corr.h — included after this.

// Design-lattice Twiss at the sextupoles, captured before errors are applied.
struct bare_optics {
  int    n_sext = 0;           // sextupoles found
  int    sexts[max_elem];      // their Cell[] indices
  double betas0_[max_elem][2], // beta_x/y at each
         nus0_[max_elem][2];   // nu_x/y   at each

  // Scan the lattice and record the above. Call on the bare (error-free)
  // machine, after Ring_GetTwiss. Was param_data_type::get_bare.
  void capture(void);
};

// Closed-orbit-correction knobs from param.dat.
struct orbit_cfg {
  std::string              loc_Fam_name;      // beam-threading start family
  std::vector<std::string> bpm_Fam_names,     // BPM families
                           corr_Fam_names[2]; // h/v corrector families
  int                      n_thread = -1;     // beam-threading iterations
  int                      n_orbit  = 5;      // SVD orbit-correction iterations
};

// Everything param.dat can set. Every default below is the value the
// corresponding param_data_type static carried; get_param overwrites only the
// keywords the file actually names, so the defaults are load-bearing.
//
// get_param is NOT a pure parser — it also loads the lattice (Read_Lattice /
// rdmfile / rdmfile_at), seeds the RNG cut (setrancut), and resolves family
// names to indices in globval (gs/ge/bpm/hcorr/vcorr/qt). Those side effects on
// global state are kept exactly as they were; untangling them is not Phase 1.
struct config_data {
  // Input files.
  std::string in_dir, ae_file, fe_file, ap_file, lat_FileName;

  // Error model / seeding.
  int    n_stat = 1;        // no of seeds
  int    n_meth = 0;        // error model: 0 = standard, 1 = cormisal (girder)
  int    n_scale = 1;       // no of steps the rms errors are ramped over
  bool   bba = false;       // beam-based alignment (code-only, no keyword)

  // Orbit correction.
  std::string loc_Fam_name = "";
  int         n_cell = -1, n_thread = -1, n_orbit = 5;
  double      h_maxkick = 1.0e-3, v_maxkick = 1.0e-3;
  double      h_cut = 1.0e-4, v_cut = 1.0e-4;
  std::vector<std::string> bpm_Fam_names, corr_Fam_names[2];

  // Coupling / vertical-dispersion correction (the LOCO off-diagonal block).
  int    n_lin = 3;
  int    SQ_per_scell = 1, BPM_per_scell = 10;
  int    HCM_per_scell = 10, VCM_per_scell = 10;
  double kick = 0.01e-3;    // trim kick used to measure the response matrix
  double VDweight = 1e3, HVweight = 1e0, VHweight = 1e0;
  double disp_wave_y = 0e0, disp_wave_o = 0e0, qt_s_cut = 1e0;
  int    qt_from_file = 0;

  // ID (insertion-device) correction.
  int    N_calls = 0, N_steps = 0;
  int    N_Fam = 0, Q_Fam[N_Fam_max];
  double ID_s_cut = 1e1;

  // Target tunes / chromaticities (fitted in err_and_corr_init when set).
  double TuneX = 0e0, TuneY = 0e0;
  double ChromX = 1e6, ChromY = 1e6;

  // Dynamic aperture.
  int    n_track_DA = 512, n_aper_DA = 15, n_delta_DA = 12;
  double delta_DA = 3e-2;
  bool   DA_bare = false;

  // Frequency map. NTURN is set to 10000 (2*NTURN for diffusion) in naffutils.h.
  bool   freq_map = false;
  int    n_x = 50, n_y = 30, n_dp = 25, n_tr = 2064;
  double x_max_FMA = 20e-3, y_max_FMA = 6e-3, delta_FMA = 3e-2;

  // Parse param.dat (and load the lattice — see the note above).
  void get_param(const std::string &param_file);

  // Project the knobs each corrector needs. Same shape for both: the corrector
  // takes its config by value and never sees the rest of the bag.
  orbit_cfg    orbit_config(void) const;
  coupling_cfg coupling_config(void) const;
};

}  // namespace corr

#endif  // CORRECTION_CORR_CONFIG_H
