#ifndef CORRECTION_CORR_CONFIG_H
#define CORRECTION_CORR_CONFIG_H

// Per-corrector configuration slices and the bare-optics reference.
//
// The param.dat reader itself (param_data_type::get_param and all the knobs it
// sets) lives in the param module, NOT here — reading param.dat is "config for a
// full-machine study", not "correction". What stays here is only what a corrector
// under correction/ actually consumes:
//
//   orbit_cfg   — the closed-orbit knobs cod_corr needs. Input, projected from
//                 param_data_type by orbit_config().
//   bare_optics — the design-lattice Twiss at the sextupoles, MEASURED from the
//                 error-free machine by capture(). Not input: it is computed
//                 from the lattice, and must be re-captured whenever the design
//                 optics change (e.g. after a tune/chroma fit). It is the
//                 reference every beta-beat / tune-shift number is quoted
//                 against.
//
// (coupling_cfg is the matching slice for the skew corrector; it lives with that
// corrector in correction/loco/coupling_corr.h.)

// Max no of quad corrector families listable on param.dat's ID_quads line — a
// config sizing limit. It sizes both param_data_type::Q_Fam and id_corr::b2, so
// it lives here (a header both include) rather than in either owner. Global
// scope, not namespaced, to match the sibling limits in id_corr.h.
const int N_Fam_max = 25;

namespace corr {

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

}  // namespace corr

#endif  // CORRECTION_CORR_CONFIG_H
