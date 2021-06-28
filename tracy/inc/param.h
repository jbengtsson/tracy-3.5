#ifndef PARAM_H
#define PARAM_H

const int N_Fam_max = 25, max_corr = 150, max_bpm = 150;

// Lattice error and correction files
const char CodCorLatFileName[]  = "codcorlat.out";

const char SkewMatFileName[]    = "skewmat.out";
const char skew_FileName[]      = "skew";
const char eta_y_FileName[]     = "eta_y";
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
  static double h_maxkick; // Default 1 mrad
  static double v_maxkick; // Default 1 mrad
  static double h_cut;  // weigthing factor cut (Default 1.0e-4)
  static double v_cut;  // weigthing factor cut (Default 1.0e-4)
  static int    n_stat; // number of statistics
  static int    n_meth; // machine errors (0=standard,1=cormisal)
  static int    n_bits; // PS resolution in amplitude in number of bits. 
  
  int h_corr[max_corr], v_corr[max_corr], bpm_loc[max_bpm];

  static double VDweight, // weight for vertical dispersion
                HVweight, // weight for coupling Htrim vertical BPM
                VHweight; // weight for coupling Vtrim horizontal BPM
  static double disp_wave_y, disp_wave_o, qt_s_cut;
  static int    qt_from_file;

  static double TuneX, // target tunes and chromaticities
                TuneY,
                ChromX,
                ChromY;

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
  double                   b2[N_Fam_max];
  static double            ID_s_cut;
  double                   **SkewRespMat, *VertCouple, *SkewStrengthCorr;
  double                   *eta_y;
  double                   *b, *w, **V, **U;

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

//-------------------------------------------------------------------
// types and variables used by GirderSetup and SetCorMis

#define reportflag      true
#define plotflag        true
#define igrmax          2000
#define ilatmax        10000
#define iseednrmax        20
 
  typedef struct girdertype {
    double gsp[2], gdx[2], gdy[2], gdt;
    long ilat[2], igir[2], gco[2], level;
   } girdertype;
  girdertype Girder[igrmax];

  long NGirderLevel [3];

  typedef struct latticetype {
    long igir;
    double smid;
  } latticetype;

  latticetype Lattice[ilatmax];

  void GirderSetup(LatticeType &lat);
  void SetCorMis(LatticeType &lat, double gxrms, double gyrms, double gtrms,
		 double jxrms, double jyrms, double exrms, double eyrms,
		 double etrms, double rancutx, double rancuty, double rancutt,
		 long iseed);
  void CorMis_in(double *gdxrms, double *gdzrms, double *gdarms,
		 double *jdxrms, double *jdzrms, double *edxrms,
		 double *edzrms, double *edarms, double *bdxrms,
		 double *bdzrms, double *bdarms, double *rancutx,
		 double *rancuty, double *rancutt, long *iseed, long *iseednr);
  
  void get_param(LatticeType &lat, const string &param_file);
  void get_bare(LatticeType &lat);
  void get_dbeta_dnu(LatticeType &lat, double m_dbeta[], double s_dbeta[],
		     double m_dnu[], double s_dnu[]);
  
// Control of vertical beam size.
  void FindSQ_SVDmat(double **SkewRespMat, double **U, double **V, double *w,
		     int N_COUPLE, int N_SKEW);
  void FindMatrix(LatticeType &lat, double **SkewRespMat,
		  const double deta_y_max, const double deta_y_offset);
  void ini_skew_cor(LatticeType &lat, const double deta_y_max,
		    const double deta_y_offset);
  void FindCoupVector(LatticeType &lat, double *VertCouple);
  void SkewStat(LatticeType &lat, double VertCouple[], const int cnt);
  void corr_eps_y(LatticeType &lat, const int cnt);
  void ReadEta(const char *TolFileName);

  // Control of IDs.
  void get_IDs(LatticeType &lat);
  void set_IDs(LatticeType &lat, const double scl);
  void reset_quads(LatticeType &lat);
  void SVD(const int m, const int n, double **M, double beta_nu[],
	   double b2Ls_[], const bool first);
  void quad_config(LatticeType &lat);
  bool get_SQ(LatticeType &lat);
  double Bet(double bq, double nus, double nuq, double NuQ);
  double Nus(double bq, double nus, double nuq, double NuQ);
  void A_matrix(void);
  void X_vector(LatticeType &lat, const bool first);
  void ini_ID_corr(LatticeType &lat, const bool IDs);
  void W_diag(LatticeType &lat);
  bool ID_corr(LatticeType &lat, const int N_calls, const int N_steps,
	       const bool IDs, const int cnt);
  void ReadCorMis(LatticeType &lat, const bool Scale_it, const double Scale)
    const;
  void LoadAlignTol(LatticeType &lat, const bool Scale_it, const double Scale,
		    const bool new_rnd, const int seed) const;
  void LoadFieldErr(LatticeType &lat, const bool Scale_it, const double Scale,
		    const bool new_rnd) const;
  void LoadApers(LatticeType &lat, const double scl_x, const double scl_y)
    const;

  void Align_BPMs(LatticeType &lat, const int n, const double bdxrms,
		  const double bdzrms, const double bdarms) const;
  bool CorrectCOD_N(LatticeType &lat, const int n_orbit, const int k);
  void ini_COD_corr(LatticeType &lat, const int n_bpm_Fam,
		    const std::string bpm_names[], const int n_hcorr_Fam,
		    const std::string hcorr_names[], const int n_vcorr_Fam,
		    const std::string vcorr_names[], const bool svd);

  bool cod_corr(LatticeType &lat, const int n_cell, const double scl,
		const double h_maxkick, const double v_maxkick,
		const long n_bits, orb_corr_type orb_corr[]);

  void Orb_and_Trim_Stat(LatticeType &lat, orb_corr_type orb_corr[]);

  void prt_cod_corr_lat(LatticeType &lat);

  void err_and_corr_init(LatticeType &lat, const string &param_file,
			 orb_corr_type orb_corr[]);

  void err_and_corr_exit(orb_corr_type orb_corr[]);
};

void get_bn2(LatticeType &lat, const string file_name1,
	     const string file_name2, int n, const bool prt);

#endif
