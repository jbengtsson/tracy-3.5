#ifndef CORRECTION_ID_CORR_H
#define CORRECTION_ID_CORR_H

// Insertion-device linear-optics correction (NOT LOCO).
//
// Extracted from param_data_type. As an ID is ramped on, its focusing perturbs
// the linear optics; this corrector zeroes the beta-beat and tune shift at the
// sextupoles by fitting thin quadrupole trims (b_2) via SVD of an analytic
// response matrix (Bet/Nus). See id_correction.md and correction_refactor.md.
//
// The struct owns the ID-correction working state (the response matrix A1, the
// distortion vector Xsext, the SVD scratch U1/w1/V1, and the per-sext/-quad
// Twiss samples). The quad-family configuration (N_Fam/Q_Fam) and the number
// of iterations (N_calls/N_steps) and the SVD cut (ID_s_cut) still live in the
// param_data_type façade (config, extracted later) and are passed in as
// explicit arguments — the file/config is an input, not hidden member state.

// Sizing limits and ID-correction weights (were in param.h; kept at global
// scope because sxt.cc/dnu_dJ.cc reference n_b3_max unqualified).
// N_Fam_max sizes a param.dat knob (config_data::Q_Fam) as well as b2 below, so
// it lives in correction/corr_config.h, which is included before this header.
const int n_b2_max  = 1500;   // max no of quad correctors
const int n_b3_max  = 1500;   // max no of sextupoles
const int max_ID_Fams = 25;   // max no of ID families

// Weights for ID correction.
const double scl_nu = 1e2, scl_dbeta = 1.0, scl_dnu = 0.1, ID_step = 0.5;

namespace corr {

// Beta response Bet and tune response Nus at phase nus to a kick at phase nuq,
// for a ring tune NuQ; bq is the beta at the kicked element.
double Bet(double bq, double nus, double nuq, double NuQ);
double Nus(double bq, double nus, double nuq, double NuQ);

// Owns the ID-correction working state (see file header).
struct id_corr {
  long int S_locs[n_b3_max];
  int      Nsext, Nquad, Nconstr, quad_prms[n_b2_max];
  int      n_ID_Fams, ID_Fams[max_ID_Fams];
  double   Ss[n_b3_max], Sq[n_b2_max], sb[2][n_b3_max], sNu[2][n_b3_max];
  double   qb[2][n_b2_max], qb0[2][n_b2_max], sNu0[2][n_b3_max];
  double   qNu0[2][n_b2_max], qNu[2][n_b2_max];
  double   Nu_X, Nu_Y, Nu_X0, Nu_Y0;
  double   **A1, *Xsext, *Xsext0, *b2Ls_, *w1, **U1, **V1;
  double   b2[N_Fam_max];  // per-family design b_2, captured by quad_config
  Vector2  dnu0, nu_0;

  id_corr() : A1(0), Xsext(0), Xsext0(0), b2Ls_(0), w1(0), U1(0), V1(0) {}
  ~id_corr();
  // Non-copyable: owns raw dvector/dmatrix allocations (freed in the dtor).
  id_corr(const id_corr &) = delete;
  id_corr &operator=(const id_corr &) = delete;

  void get_IDs(void);
  void set_IDs(const double scl);
  void reset_quads(const int N_Fam, const int Q_Fam[]);
  void SVD(const int m, const int n, double **M, double beta_nu[],
	   double b2Ls_[], const bool first, const double ID_s_cut);
  void quad_config(const int N_Fam, const int Q_Fam[]);
  bool get_SQ(void);
  void A_matrix(void);
  void X_vector(const bool first);
  void ini_ID_corr(const bool IDs, const int N_Fam, const int Q_Fam[]);
  void W_diag(void);
  bool ID_corr(const int N_calls, const int N_steps, const bool IDs,
	       const int cnt, const int N_Fam, const int Q_Fam[],
	       const double ID_s_cut);
};

}  // namespace corr

#endif  // CORRECTION_ID_CORR_H
