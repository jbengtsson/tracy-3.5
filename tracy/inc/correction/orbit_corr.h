#ifndef CORRECTION_ORBIT_CORR_H
#define CORRECTION_ORBIT_CORR_H

// Closed-orbit-distortion (COD) correction driver.
//
// Extracted from param_data_type. This is a thin orchestration layer over the
// already-modular orb_corr_type (orb_corr.h) and the lsoc response-matrix
// primitives (gcmat/gtcmat): build the corrector->BPM response matrices, then
// per seed find/thread the closed orbit and SVD-correct it.
// See correction_refactor.md.

class param_data_type;  // transitional; see corr::cod_corr below.

namespace corr {

// Build the horizontal/vertical orbit response matrices (and their transposes)
// from the BPM and corrector families. Self-contained (globals only).
void ini_COD_corr(const int n_bpm_Fam, const std::string bpm_names[],
                  const int n_hcorr_Fam, const std::string hcorr_names[],
                  const int n_vcorr_Fam, const std::string vcorr_names[],
                  const bool svd);

// Correct the closed orbit for the current machine state: clear trims, find the
// COD (threading the beam if none exists), then SVD-correct via orb_corr[].
// Takes param_data_type& only to read orbit-control config members
// (loc_Fam_name/n_thread/n_orbit/bpm_/corr_Fam_names) and the bare-lattice
// reference (n_sext/sexts/betas0_/nus0_); it no longer calls back into any
// param method. This shrinks to a plain config struct once config is extracted.
bool cod_corr(param_data_type &p, const int n_cell, const double scl,
              const double h_maxkick, const double v_maxkick,
              orb_corr_type orb_corr[]);

// Report orbit-at-sextupole and trim-strength statistics. Self-contained.
void Orb_and_Trim_Stat(orb_corr_type orb_corr[]);

}  // namespace corr

#endif  // CORRECTION_ORBIT_CORR_H
