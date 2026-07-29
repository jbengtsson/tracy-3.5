#ifndef CORRECTION_LOCO_COUPLING_CORR_H
#define CORRECTION_LOCO_COUPLING_CORR_H

// Coupling / vertical-dispersion correction — LOCO for the off-diagonal block:
//
//   1. build the MODEL response of each skew quad on
//        [vertical dispersion at the BPMs,
//         h-trim -> vertical BPM orbit,
//         v-trim -> horizontal BPM orbit]                   (find_model_matrix)
//   2. SVD it, cutting small singular values                (corr_linalg)
//   3. MEASURE the same three blocks on the real machine    (find_coup_vector,
//                                                            over loco/orm.h)
//   4. solve for the skew strengths and apply, iterating    (corr_eps_y)
//
// Index convention is Numerical-Recipes 1-based for the response matrix and the
// vectors; the element-location arrays (bpm_loc/h_corr/v_corr) are 0-based.

// Sizing caps for the BPM / corrector position arrays. Used only on the n_lin>0
// (coupling) path.
// TODO: BPM/corrector position arrays should become dynamic.
const int max_corr = 800, max_bpm = 800;

// Files written by the coupling corrector.
const char SkewMatFileName[]    = "skewmat.out";
const char skew_FileName[]      = "skew";
const char eta_y_FileName[]     = "eta_y";
const char deta_y_FileName[]    = "deta_y.out";

namespace corr {

// Knobs set by param.dat.
struct coupling_cfg {
  double VDweight;      // weight on the vertical-dispersion block
  double HVweight;      // weight on the h-trim -> vertical-BPM block
  double VHweight;      // weight on the v-trim -> horizontal-BPM block
  double qt_s_cut;      // SVD singular-value cut for the skew-quad fit
  double kick;          // trim kick used to measure the response [rad]
  int    n_lin;         // number of correction iterations
  int    SQ_per_scell;  // skew quads per super-cell (target eta_y wave)
  int    qt_from_file;  // if set, load skew strengths from qt_file.dat instead
};

struct coupling_corr {
  int    N_BPM, N_HCOR, N_VCOR, N_SKEW, N_COUPLE;
  int    h_corr[max_corr], v_corr[max_corr], bpm_loc[max_bpm];
  double **SkewRespMat, *VertCouple, *SkewStrengthCorr, *eta_y;
  double *b, *w, **V, **U;

  coupling_corr()
    : N_BPM(0), N_HCOR(0), N_VCOR(0), N_SKEW(0), N_COUPLE(0), SkewRespMat(0),
      VertCouple(0), SkewStrengthCorr(0), eta_y(0), b(0), w(0), V(0), U(0) {}
  ~coupling_corr();
  // Owns raw dvector/dmatrix allocations, freed in the dtor.
  coupling_corr(const coupling_corr &) = delete;
  coupling_corr &operator=(const coupling_corr &) = delete;
  // Read the target vertical dispersion from TolFileName into eta_y (mm -> m).
  void read_eta(const char *TolFileName);

  // Build the model response of each skew quad on the vertical dispersion at
  // the BPMs and on the two off-diagonal trim->orbit blocks. Also builds the
  // eta_y target wave from deta_y_max / deta_y_offset.
  void find_model_matrix(const coupling_cfg &cfg, const double deta_y_max,
			 const double deta_y_offset);

  // Locate the skew quads, BPMs and trims, allocate the response matrix, build
  // it, and SVD-decompose it with a singular-value cut. Call once before
  // corr_eps_y.
  void ini_skew_cor(const coupling_cfg &cfg, const double deta_y_max,
		    const double deta_y_offset);

  // Measure the same three blocks on the real machine: kick each trim and
  // finite-difference the resulting closed orbit at the BPMs. Fills VertCouple.
  void find_coup_vector(const coupling_cfg &cfg, double *VertCouple);

  // Report the residual coupling in VertCouple. cnt < 0 prints to stdout, else
  // writes skew_<cnt>.out.
  void skew_stat(const coupling_cfg &cfg, double VertCouple[], const int cnt);

  // Solve for the skew-quad strengths and apply them, re-measuring and
  // iterating cfg.n_lin times. Needs ini_skew_cor first.
  void corr_eps_y(const coupling_cfg &cfg, const int cnt);
};

}  // namespace corr

#endif  // CORRECTION_LOCO_COUPLING_CORR_H
