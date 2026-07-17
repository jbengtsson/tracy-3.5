#ifndef CORRECTION_ORBIT_CORR_H
#define CORRECTION_ORBIT_CORR_H

// Closed-orbit-distortion (COD) correction driver.
//
// A thin orchestration layer over orb_corr_type (orb_corr.h) and the lsoc
// response-matrix primitives (gcmat/gtcmat): build the corrector->BPM response
// matrices, then per seed find/thread the closed orbit and SVD-correct it.

namespace corr {

// Build the horizontal/vertical orbit response matrices (and their transposes)
// from the BPM and corrector families. Self-contained: it configures the lsoc
// globals, not the orbit_corr state below.
void ini_COD_corr(const int n_bpm_Fam, const std::string bpm_names[],
		  const int n_hcorr_Fam, const std::string hcorr_names[],
		  const int n_vcorr_Fam, const std::string vcorr_names[],
		  const bool svd);

// The horizontal/vertical orbit correctors, owning the orb_corr_type pair. It
// is held here rather than by each application so that orb_corr_type stays out
// of param.h and dynap.h.
class orbit_corr {
private:
  orb_corr_type orb_corr[2];

public:
  orbit_corr(void) {}
  orbit_corr(const orbit_corr &) = delete;
  orbit_corr &operator=(const orbit_corr &) = delete;

  // Build both planes' response matrices from the BPM/corrector families.
  void alloc(const std::vector<string> &bpm_Fam_names,
	     const std::vector<string> corr_Fam_names[]);
  void dealloc(void);

  // Correct the closed orbit for the current machine state: clear trims, find
  // the COD (threading the beam if none exists), then SVD-correct. Takes its
  // config (cfg) and the bare-lattice reference it reports the residual
  // beta-beat against (bare).
  bool cod_corr(const orbit_cfg &cfg, const bare_optics &bare, const int n_cell,
		const double scl, const double h_maxkick,
		const double v_maxkick);

  // Report orbit-at-sextupole and trim-strength statistics.
  void Orb_and_Trim_Stat(void);

  // Dump both planes' SVD matrices (diagnostic, trace only).
  void prt_svdmat(void);
};

}  // namespace corr

#endif  // CORRECTION_ORBIT_CORR_H
