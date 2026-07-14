#ifndef PARAM_H
#define PARAM_H

// N_Fam_max lives in correction/id_corr.h; max_bpm/max_corr and the skew/eta_y
// output file names in correction/loco/coupling_corr.h. Both are included
// before this header.

// Computation result files
const char beam_envelope_file[] = "beam_envelope";

// Lattice error and correction files
const char CodCorLatFileName[]  = "codcorlat.out";

// N_Fam_max, n_b2_max, n_b3_max, max_ID_Fams and the ID-correction weights
// (scl_nu/scl_dbeta/scl_dnu/ID_step) moved to correction/id_corr.h (included
// before this header), where the corr::id_corr working state now lives.

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

  static double kick;      // 0.01 mrad kick for trims
  static double h_maxkick; // Default 1 mrad
  static double v_maxkick; // Default 1 mrad
  static double h_cut;     // weigthing factor cut (Default 1.0e-4)
  static double v_cut;     // weigthing factor cut (Default 1.0e-4)
  static int    n_stat;    // number of statistics
  static int    n_meth;    // machine errors (0=standard,1=cormisal)
  
  std::vector<double> bn_an[2*HOMmax+1];

  static double VDweight,  // weight for vertical dispersion
                HVweight,  // weight for coupling Htrim vertical BPM
                VHweight;  // weight for coupling Vtrim horizontal BPM
  static double disp_wave_y, disp_wave_o, qt_s_cut;
  static int    qt_from_file;

  static double TuneX,     // target tunes and chromaticities
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

  // Orbit control.
  static std::string       loc_Fam_name;
  static int               n_cell, n_thread;
  std::vector<std::string> bpm_Fam_names, corr_Fam_names[2];
  static bool              bba;

  // ID control.
  int                      N_calls, N_steps, N_Fam, Q_Fam[N_Fam_max];
  static double            ID_s_cut;

  // Bare-lattice reference optics at the sextupoles (n_sext/sexts/betas0_/nus0_
  // moved into correction/config as corr::bare_optics, which get_bare() now
  // fills). Dropped with the move: nu0_, which get_bare() wrote and nothing ever
  // read.
  corr::bare_optics bare;

  // ID (insertion-device) linear-optics correction. The working state (response
  // matrix, distortion vector, SVD scratch, per-sext/-quad Twiss, per-family b2)
  // moved to correction/id_corr. N_calls/N_steps/N_Fam/Q_Fam above and the
  // static ID_s_cut stay here as config (extracted later) and are passed in.
  corr::id_corr id;

//-------------------------------------------------------------------
// Cormisal (girder) error model (n_meth == 1). The girder-tree state and the
// GirderSetup/SetCorMis algorithms were extracted to correction/girder_model
// (types, Girder[]/Lattice[]/NGirderLevel, and the igrmax/ilatmax/iseednrmax
// limits now live there). The methods below remain as façade delegators.

  corr::girder_model girders;

  void GirderSetup();
  void SetCorMis(double gxrms, double gyrms, double gtrms, double jxrms,
		 double jyrms, double exrms, double eyrms, double etrms,
		 double rancutx, double rancuty, double rancutt, long iseed);
  void CorMis_in(double *gdxrms, double *gdzrms, double *gdarms,
		 double *jdxrms, double *jdzrms, double *edxrms,
		 double *edzrms, double *edarms, double *bdxrms,
		 double *bdzrms, double *bdarms, double *rancutx,
		 double *rancuty, double *rancutt, long *iseed, long *iseednr);
  
  void get_param(const string &param_file);
  void get_bare(void);
  void get_dbeta_dnu(double m_dbeta[], double s_dbeta[], double m_dnu[],
		     double s_dnu[]);

// Orbit-correction knobs (loc_Fam_name, n_thread, n_orbit + the BPM/corrector
// family names) stay here as config; orbit_config() packs them for the
// corrector, as coupling_config() does for the skew corrector.
  corr::orbit_cfg orbit_config(void) const;
  
// Control of vertical beam size. The knobs (n_lin, the three weights, qt_s_cut,
// kick, SQ_per_scell, qt_from_file) stay here as config; coupling_config()
// packs them for the corrector.
  corr::coupling_corr skew;

  corr::coupling_cfg coupling_config(void) const;
  void ini_skew_cor(const double deta_y_max, const double deta_y_offset);
  void corr_eps_y(const int cnt);

  // Control of IDs — thin façades delegating to the corr::id_corr member
  // above (still called by dynap/leac/touschek and err_and_corr_init). The
  // remaining ID methods (get_IDs/set_IDs/SVD/quad_config/get_SQ/Bet/Nus/
  // A_matrix/X_vector/W_diag) had no external callers and now live only on the
  // corr::id_corr struct.
  void reset_quads(void);
  void ini_ID_corr(const bool IDs);
  bool ID_corr(const int N_calls, const int N_steps, const bool IDs,
	       const int cnt);
  void ReadCorMis(const bool Scale_it, const double Scale) const;
  void LoadAlignTol(const bool Scale_it, const double Scale,
		    const bool new_rnd,
		    const int seed) const;
  void LoadFieldErr(const bool Scale_it, const double Scale,
		    const bool new_rnd) const;
  void LoadApers(const double scl_x, const double scl_y) const;
  void zero_mult(void);
  void restore_mult(void);
  void Align_BPMs(const int n, const double bdxrms, const double bdzrms,
		  const double bdarms) const;
  bool CorrectCOD_N(const int n_orbit, const int k);
  void ini_COD_corr(const int n_bpm_Fam, const std::string bpm_names[],
		    const int n_hcorr_Fam, const std::string hcorr_names[],
		    const int n_vcorr_Fam, const std::string vcorr_names[],
		    const bool svd);

  bool cod_corr(const int n_cell, const double scl, const double h_maxkick,
		const double v_maxkick, orb_corr_type orb_corr[]);

  void Orb_and_Trim_Stat(orb_corr_type orb_corr[]);

  void prt_cod_corr_lat(void);

  void err_and_corr_init(const string &param_file, orb_corr_type orb_corr[]);

  void err_and_corr_exit(orb_corr_type orb_corr[]);
};

void get_bn2(const string file_name1, const string file_name2, int n,
	     const bool prt);

#endif
