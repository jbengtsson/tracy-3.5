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
