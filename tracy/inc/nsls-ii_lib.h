// global params

const int     max_elem = Cell_nLocMax;

extern char   in_dir[];

extern int n_aper, n_track;

extern ss_vect<tps>  map;
extern MNF_struct    MNF;

extern double chi_m;


const int N_Fam_max = 15, max_corr = 100, max_bpm = 200;

// Computation result files
const char beam_envelope_file[] = "beam_envelope.out";

// Lattice error and correction files
const char CodCorLatFileName[]  = "codcorlat.out";

const char SkewMatFileName[]    = "skewmat.out";
const char eta_y_FileName[]     = "eta_y.out";
const char deta_y_FileName[]    = "deta_y.out";

const int n_b2_max    = 500;  // max no of quad corrector families
const int n_b3_max    = 1000; // max no of sextupoles
const int max_ID_Fams = 25;   // max no of ID families

// Weights for ID correction
const double scl_nu = 1e2, scl_dbeta = 1.0, scl_dnu = 0.1, ID_step = 0.5;

class param_data_type {
 public:
  string ae_file, fe_file, ap_file, in_dir, lat_FileName;

  bool   DA_bare       = false,
         freq_map      = false;
  int    n_orbit       = 3,
         n_scale       = 1;

  int    n_lin         =  3,
         SQ_per_scell  =  2,
         BPM_per_scell = 12,
         HCM_per_scell = 12,
         VCM_per_scell = 12;

  double kick          = 0.01e-3; // 0.01 mrad kick for trims
  int    n_stat        = 1;       // number of statistics

  int h_corr[max_corr], v_corr[max_corr], bpm_loc[max_bpm];

  double VDweight      = 1e3,     // weight for vertical dispersion
         HVweight      = 1e0,     // weight for coupling Htrim vertical BPM
         VHweight      = 1e0;     // weight for coupling Vtrim horizontal BPM

  // Parameters for dynamic aperture
  double  delta_DA_    = 5.0e-2;

  // Parameters for frequency map
  // Note NTURN is set to 10000 (2*NTURN for diffusion)) in "naffutils.h".
  int    n_x = 50, n_y = 30, n_dp = 25, n_tr = 2064;
  double x_max_FMA = 20e-3, y_max_FMA = 6e-3, delta_FMA = 3e-2;
  //double x_max_FMA = 20e-3, y_max_FMA = 3e-3, delta_FMA = 3e-2;

  int                      N_BPM, N_HCOR, N_VCOR, N_SKEW, N_COUPLE;
  // Orbit control.
  std::vector<std::string> bpm_Fam_names, corr_Fam_names[2];
  bool   bba     = false;

  // ID control.
  int                      N_calls, N_steps, N_Fam, Q_Fam[N_Fam_max];
  int                      n_sext, sexts[max_elem];
  double                   betas0_[max_elem][2], nus0_[max_elem][2], nu0_[2];
  double                   b2[N_Fam_max];
  double                   **SkewRespMat, *VertCouple, *SkewStrengthCorr;
  double                   *eta_y;
  double                   *b, *w, **V, **U;
  double                   disp_wave_y;

  void get_param(const string &param_file);
  void get_bare(void);
  void get_dbeta_dnu(double m_dbeta[], double s_dbeta[], double m_dnu[],
		     double s_dnu[]);
  
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
	   double b2Ls_[], const double s_cut, const bool first);
  void quad_config();
  bool get_SQ(void);
  double Bet(double bq, double nus, double nuq, double NuQ);
  double Nus(double bq, double nus, double nuq, double NuQ);
  void A_matrix(void);
  void X_vector(const bool first);
  void ini_ID_corr(const bool IDs);
  void W_diag(void);
  bool ID_corr(const int N_calls, const int N_steps, const bool IDs);
  
  void LoadAlignTol(const bool Scale_it, const double Scale, const bool new_rnd,
		    const int seed) const;
  void LoadFieldErr(const bool Scale_it, const double Scale,
		    const bool new_rnd) const;
  void LoadApers(const double scl_x, const double scl_y) const;

  void Align_BPMs(const int n) const;
  bool CorrectCOD_N(const int n_orbit,
		    const int n_scale, const int k);
  void ini_COD_corr(const int n_bpm_Fam, const std::string bpm_names[],
		    const int n_hcorr_Fam, const std::string hcorr_names[],
		    const int n_vcorr_Fam, const std::string vcorr_names[],
		    const bool svd);
  void Orb_and_Trim_Stat(void);
  void prt_codcor_lat(void);
};


void lwr_case(char str[]);

void upr_case(char str[]);

//void prt_trace (void);

void file_rd(std::ifstream &inf, const char file_name[]);

void file_wr(std::ofstream &outf, const char file_name[]);

FILE* file_read(const char file_name[]);

FILE* file_write(const char file_name[]);

void chk_cod(const bool cod, const char *proc_name);

void no_sxt(void);

void get_map(const bool cod);

tps get_h(void);

void get_m2(const ss_vect<tps> &ps, tps m2[]);

ss_vect<tps> get_S(const int n_DOF);

ss_vect<tps> tp_S(const int n_DOF, const ss_vect<tps> &A);

void get_dnu(const int n, const ss_vect<tps> &A, double dnu[]);

ss_vect<tps> get_A_CS(const int n, const ss_vect<tps> &A, double dnu[]);

void prt_lin_map(const int n_DOF, const ss_vect<tps> &map);

void get_twoJ(const int n_DOF, const ss_vect<double> &ps,
	      const ss_vect<tps> &A, double twoJ[]);

double get_curly_H(const double alpha_x, const double beta_x,
		   const double eta_x, const double etap_x);

double get_eps_x(void);

void GetEmittance(const int Fnum, const bool prt);

void prt_lat(const char *fname, const int Fnum, const bool all);

void Cell_Twiss(const long int i0, const long int i1);

void prt_lat(const char *fname, const int Fnum, const bool all, const int n);

void prt_chrom_lat(void);

void prt_cod(const char *file_name, const int Fnum, const bool all);

void prt_beampos(const char *file_name);

void CheckAlignTol(const char *OutputFile);

void misalign_rms_elem(const int Fnum, const int Knum,
		       const double dx_rms, const double dy_rms,
		       const double dr_rms, const bool new_rnd);

void misalign_sys_elem(const int Fnum, const int Knum,
		       const double dx_sys, const double dy_sys,
		       const double dr_sys);

void misalign_rms_fam(const int Fnum,
		      const double dx_rms, const double dy_rms,
		      const double dr_rms, const bool new_rnd);

void misalign_sys_fam(const int Fnum,
		      const double dx_sys, const double dy_sys,
		      const double dr_sys);

void misalign_rms_type(const int type,
		       const double dx_rms, const double dy_rms,
		       const double dr_rms, const bool new_rnd);

void misalign_sys_type(const int type,
		       const double dx_sys, const double dy_sys,
		       const double dr_sys);

void misalign_rms_girders(const int gs, const int ge,
			  const double dx_rms, const double dy_rms,
			  const double dr_rms, const bool new_rnd);

void misalign_sys_girders(const int gs, const int ge,
			  const double dx_sys, const double dy_sys,
			  const double dr_sys);

void LoadAlignTol(const char *AlignFile, const bool Scale_it,
		  const double Scale, const bool new_rnd, const int k);

void set_aper_elem(const int Fnum, const int Knum,
		   const double Dxmin, const double Dxmax,
		   const double Dymin, const double Dymax);

void set_aper_fam(const int Fnum,
		  const double Dxmin, const double Dxmax,
		  const double Dymin, const double Dymax);

void set_aper_type(const int type, const double Dxmin, const double Dxmax,
		   const double Dymin, const double Dymax);

void LoadApers(const char *AperFile, const double scl_x, const double scl_y);

double get_L(const int Fnum, const int Knum);

void set_L(const int Fnum, const int Knum, const double L);

void set_L(const int Fnum, const double L);

void set_dL(const int Fnum, const int Knum, const double dL);

void set_dL(const int Fnum, const double dL);

void get_bn_design_elem(const int Fnum, const int Knum,
			const int n, double &bn, double &an);

void get_bnL_design_elem(const int Fnum, const int Knum,
			 const int n, double &bnL, double &anL);

void set_bn_design_elem(const int Fnum, const int Knum,
			const int n, const double bn, const double an);

void set_dbn_design_elem(const int Fnum, const int Knum,
			 const int n, const double dbn, const double dan);

void set_bn_design_fam(const int Fnum,
		       const int n, const double bn, const double an);

void set_dbn_design_fam(const int Fnum,
			const int n, const double dbn, const double dan);

void set_bnL_design_elem(const int Fnum, const int Knum,
			 const int n, const double bnL, const double anL);

void set_dbnL_design_elem(const int Fnum, const int Knum,
			 const int n, const double dbnL, const double danL);

void set_bnL_design_fam(const int Fnum,
			const int n, const double bnL, const double anL);

void set_dbnL_design_fam(const int Fnum,
			 const int n, const double dbnL, const double danL);

void set_bnL_design_type(const int type,
			 const int n, const double bnL, const double anL);

void set_bnL_sys_elem(const int Fnum, const int Knum,
		      const int n, const double bnL, const double anL);

void set_bnL_sys_fam(const int Fnum,
		     const int n, const double bnL, const double anL);

void set_bnL_sys_type(const int type,
		      const int n, const double bnL, const double anL);

void set_bnL_rms_elem(const int Fnum, const int Knum,
		      const int n, const double bnL, const double anL,
		      const bool new_rnd);

void set_bnL_rms_fam(const int Fnum,
		     const int n, const double bnL, const double anL,
		     const bool new_rnd);

void set_bnL_rms_type(const int type,
		      const int n, const double bnL, const double anL,
		      const bool new_rnd);

void set_bnr_sys_elem(const int Fnum, const int Knum,
		      const int n, const double bnr, const double anr);

void set_bnr_sys_fam(const int Fnum,
		     const int n, const double bnr, const double anr);

void set_bnr_sys_type(const int type,
		      const int n, const double bnr, const double anr);

void set_bnr_rms_elem(const int Fnum, const int Knum,
		      const int n, const double bnr, const double anr,
		      const bool new_rnd);

void set_bnr_rms_fam(const int Fnum,
		     const int n, const double bnr, const double anr,
		     const bool new_rnd);

void set_bnr_rms_type(const int type,
		      const int n, const double bnr, const double anr,
		      const bool new_rnd);

double get_Wiggler_BoBrho(const int Fnum, const int Knum);

void set_Wiggler_BoBrho(const int Fnum, const int Knum, const double BoBrhoV);

void set_Wiggler_BoBrho(const int Fnum, const double BoBrhoV);

void set_ID_scl(const int Fnum, const int Knum, const double scl);

void SetFieldValues_fam(const int Fnum, const bool rms, const double r0,
			const int n, const double Bn, const double An,
			const bool new_rnd);

void SetFieldValues_type(const int N, const bool rms, const double r0,
			 const int n, const double Bn, const double An,
			 const bool new_rnd);

void SetFieldErrors(const char *name, const bool rms, const double r0,
		    const int n, const double Bn, const double An,
		    const bool new_rnd);

void LoadFieldErr(const char *FieldErrorFile, const bool Scale_it,
		  const double Scale, const bool new_rnd);

bool CorrectCOD(const int n_orbit, const double scl);

void Align_BPMs(const int n);

void get_bare();

void get_dbeta_dnu(double m_dbeta[], double s_dbeta[],
		   double m_dnu[], double s_dnu[]);

bool CorrectCOD_N(const char *ae_file, const int n_orbit,
		  const int n_scale, const int k);

void ini_skew_cor(const double deta_y_max);

void corr_eps_y(void);

void get_IDs(void);

void set_IDs(const double scl);

void reset_quads(void);

void ini_ID_corr(const bool IDs);

bool ID_corr(const int N_calls, const int N_steps, const bool IDs);

void ini_COD_corr(const int n_bpm_Fam, const std::string bpm_names[],
		  const int n_hcorr_Fam, const std::string hcorr_names[],
		  const int n_vcorr_Fam, const std::string vcorr_names[],
		  const bool svd);

void get_param(const char *param_file, param_data_type &params);

void error_and_correction(const char *param_file, bool rd_lat);

void prt_codcor_lat(void);

void prt_beamsizes();

double Touschek(const double Qb, const double delta_RF,
		const double eps_x, const double eps_y,
		const double sigma_delta, const double sigma_s);

double Touschek(const double Qb, const double delta_RF,const bool consistent,
		const double eps_x, const double eps_y,
		const double sigma_delta, double sigma_s,
		const int n_turn, const bool aper_on,
		double sum_delta[][2], double sum2_delta[][2]);

double f_IBS(const double chi_m);

double get_int_IBS(void);

void IBS(const double Qb, const double eps_SR[], double eps[],
	 const bool prt1, const bool prt2);

void IBS_BM(const double Qb, const double eps_SR[], double eps[],
	    const bool prt1, const bool prt2);

void rm_space(char *name);

void get_bn(const char file_name[], int n, const bool prt);

double get_dynap(const double delta, const bool cod);

double get_chi2(long int n, double x[], double y[], long int m, psVector b);

void pol_fit(int n, double x[], double y[], int order, psVector &b,
	     double &sigma, const bool prt);

void get_ksi2(const double d_delta);

bool find_nu(const int n, const double nus[], const double eps, double &nu);

bool get_nu(const double Ax, const double Ay, const double delta,
	    double &nu_x, double &nu_y);

void dnu_dA(const double Ax_max, const double Ay_max, const double delta, const int n_ampl);

bool orb_corr(const int n_orbit);

double get_code(CellType &Cell);

void get_alphac(void);

void get_alphac2(void);

void bend_cal_Fam(const int Fnum);

void bend_cal(void);

double h_ijklm(const tps &h, const int i, const int j, const int k,
	       const int l, const int m);

ss_vect<tps> get_A(const double alpha[], const double beta[],
		   const double eta[], const double etap[]);


void get_ab(const ss_vect<tps> &A,
	    double alpha[], double beta[], double nu[],
	    double eta[], double etap[]);

void set_tune(const char file_name1[], const char file_name2[], const int n);

void prt_H_long(const int n, const double phi_max, const double delta_max,
		const double U0);
