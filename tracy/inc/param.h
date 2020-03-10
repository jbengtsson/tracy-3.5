#ifndef PARAM_H
#define PARAM_H

const int N_Fam_max = 20, max_corr = 100, max_bpm = 200;

// Computation result files
const char beam_envelope_file[] = "beam_envelope.out";

// Lattice error and correction files
const char CodCorLatFileName[]  = "codcorlat.out";

const char SkewMatFileName[]    = "skewmat.out";
const char eta_y_FileName[]     = "eta_y.out";
const char deta_y_FileName[]    = "deta_y.out";

const int n_b2_max    = 1500;  // max no of quad corrector families
const int n_b3_max    = 1500; // max no of sextupoles
const int max_ID_Fams = 25;   // max no of ID families

// Weights for ID correction
const double scl_nu = 1e2, scl_dbeta = 1.0, scl_dnu = 0.1, ID_step = 0.5;

class param_data_type {
 private:

 public:
  string ae_file, fe_file, ap_file, in_dir, lat_FileName;

  static bool DA_bare,
              freq_map;
  static int  n_orbit,
              n_scale;

  static int n_lin,
             SQ_per_scell,
             BPM_per_scell,
             HCM_per_scell,
             VCM_per_scell;

  static double kick;   // 0.01 mrad kick for trims
  static int    n_stat; // number of statistics

  int h_corr[max_corr], v_corr[max_corr], bpm_loc[max_bpm];

  static double VDweight, // weight for vertical dispersion
                HVweight, // weight for coupling Htrim vertical BPM
                VHweight; // weight for coupling Vtrim horizontal BPM

  // Parameters for dynamic aperture
  static int    n_track_DA,
                n_aper_DA,
                n_delta_DA;
  static double delta_DA;

  // Parameters for frequency map
  // Note NTURN is set to 10000 (2*NTURN for diffusion)) in "naffutils.h".
   static int    n_x, n_y, n_dp, n_tr;
   static double x_max_FMA, y_max_FMA, delta_FMA;

  int                      N_BPM, N_HCOR, N_VCOR, N_SKEW, N_COUPLE;
  // Orbit control.
  static std::string       loc_Fam_name;
  static int               n_cell, n_thread;
  std::vector<std::string> bpm_Fam_names, corr_Fam_names[2];
  static bool              bba;

  // ID control.
  int                      N_calls, N_steps, N_Fam, Q_Fam[N_Fam_max];
  int                      n_sext, sexts[max_elem];
  double                   betas0_[max_elem][2], nus0_[max_elem][2], nu0_[2];
  double                   b2[N_Fam_max], ID_s_cut;
  double                   **SkewRespMat, *VertCouple, *SkewStrengthCorr;
  double                   *eta_y;
  double                   *b, *w, **V, **U;
  double                   disp_wave_y;

  // ID_corr global variables
  long int S_locs[n_b3_max];
  int      Nsext, Nquad, Nconstr, NconstrO, quad_prms[n_b2_max], id_loc;
  int      n_ID_Fams, ID_Fams[max_ID_Fams];
  double   Ss[n_b3_max], Sq[n_b2_max], sb[2][n_b3_max], sNu[2][n_b3_max];
  double   qb[2][n_b2_max], qb0[2][n_b2_max], sNu0[2][n_b3_max];
  double   qNu0[2][n_b2_max], qNu[2][n_b2_max], IDb[2], IDNu[2];
  double   Nu_X, Nu_Y, Nu_X0, Nu_Y0;
  double   **A1, *Xsext, *Xsext0, *b2Ls_, *w1, **U1, **V1;
  double   *Xoct, *b4s, **Aoct;
  Vector2  dnu0, nu_0;

  void get_param(const string &param_file);
  void get_bare(void);
  void get_dbeta_dnu(double m_dbeta[], double s_dbeta[], double m_dnu[],
		     double s_dnu[]);
  
// Control of vertical beam size.
  void FindSQ_SVDmat(double **SkewRespMat, double **U, double **V, double *w,
		     int N_COUPLE, int N_SKEW);
  void FindMatrix(double **SkewRespMat, const double deta_y_max);
  void ini_skew_cor(const double deta_y_max);
  void FindCoupVector(double *VertCouple);
  void SkewStat(double VertCouple[]);
  void corr_eps_y(void);

  // Control of IDs.
  void get_IDs(void);
  void set_IDs(const double scl);
  void reset_quads(void);
  void SVD(const int m, const int n, double **M, double beta_nu[],
	   double b2Ls_[], const bool first);
  void quad_config();
  bool get_SQ(void);
  double Bet(double bq, double nus, double nuq, double NuQ);
  double Nus(double bq, double nus, double nuq, double NuQ);
  void A_matrix(void);
  void X_vector(const bool first);
  void ini_ID_corr(const bool IDs);
  void W_diag(void);
  bool ID_corr(const bool IDs, orb_corr_type orb_corr[]);
  
  void LoadAlignTol(const bool Scale_it, const double Scale, const bool new_rnd,
		    const int seed) const;
  void LoadFieldErr(const bool Scale_it, const double Scale,
		    const bool new_rnd) const;
  void LoadApers(const double scl_x, const double scl_y) const;

  void Align_BPMs(const int n) const;
  bool CorrectCOD_N(const int n_orbit, const int k);
  void ini_COD_corr(const int n_bpm_Fam, const std::string bpm_names[],
		    const int n_hcorr_Fam, const std::string hcorr_names[],
		    const int n_vcorr_Fam, const std::string vcorr_names[],
		    const bool svd);

  bool cod_corr(const double scl, orb_corr_type orb_corr[]);

  void Orb_and_Trim_Stat(void);

  void prt_cod_corr_lat(void);

  void err_and_corr_init(const string &param_file, orb_corr_type orb_corr[]);

  void err_and_corr_exit(orb_corr_type orb_corr[]);
};

void get_bn2(const string file_name1, const string file_name2, int n,
	     const bool prt);

#endif
