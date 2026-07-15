#ifndef CORRECTION_ORBIT_CORR_H
#define CORRECTION_ORBIT_CORR_H

// Closed-orbit-distortion (COD) correction driver.
//
// Extracted from param_data_type. This is a thin orchestration layer over the
// already-modular orb_corr_type (orb_corr.h) and the lsoc response-matrix
// primitives (gcmat/gtcmat): build the corrector->BPM response matrices, then
// per seed find/thread the closed orbit and SVD-correct it.

namespace corr {

// Build the horizontal/vertical orbit response matrices (and their transposes)
// from the BPM and corrector families. Self-contained (globals only).
void ini_COD_corr(const int n_bpm_Fam, const std::string bpm_names[],
                  const int n_hcorr_Fam, const std::string hcorr_names[],
                  const int n_vcorr_Fam, const std::string vcorr_names[],
                  const bool svd);

// Correct the closed orbit for the current machine state: clear trims, find the
// COD (threading the beam if none exists), then SVD-correct via orb_corr[].
// Takes only its config (cfg) and the bare-lattice reference it reports the
// residual beta-beat against (bare) — no god-class dependency.
bool cod_corr(const orbit_cfg &cfg, const bare_optics &bare, const int n_cell,
              const double scl, const double h_maxkick, const double v_maxkick,
              orb_corr_type orb_corr[]);

// Report orbit-at-sextupole and trim-strength statistics. Self-contained.
void Orb_and_Trim_Stat(orb_corr_type orb_corr[]);

}  // namespace corr

#endif  // CORRECTION_ORBIT_CORR_H
