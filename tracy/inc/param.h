#ifndef PARAM_H
#define PARAM_H

// N_Fam_max lives in correction/corr_config.h; n_b2_max/n_b3_max/max_ID_Fams and
// the ID-correction weights in correction/id_corr.h; max_bpm/max_corr and the
// skew/eta_y output file names in correction/loco/coupling_corr.h. All are
// included before this header.

// Computation result files
const char beam_envelope_file[] = "beam_envelope";

// Lattice error and correction files
const char CodCorLatFileName[] = "codcorlat.out";

// The param.dat reader (and, transitionally, a correction façade).
//
// This class owns get_param and every param.dat knob — reading the full-machine
// study config is a param-module job, not a correction one, so it lives here and
// not under correction/. The apps and tracy/src/dynap.cc read the knobs as
// params.<knob> (params.n_cell, params.h_maxkick, ...) in ~90 places. Each
// corrector under correction/ takes only the thin slice it needs, projected by
// orbit_config()/coupling_config(); nothing under correction/ sees this whole bag.
//
// Still transitional: the corrector objects below plus a wall of delegators.
// A later step repoints callers at the correction classes directly and strips
// the delegators, leaving this class as the param reader (NOT deleted).

class param_data_type
{
private:
public:
  // ------------------------------------------------------------------
  // param.dat knobs. Every default below is the value the corresponding
  // param_data_type static once carried; get_param overwrites only the keywords
  // the file actually names, so the defaults are load-bearing.

  // Input files.
  std::string in_dir, ae_file, fe_file, ap_file, lat_FileName;

  // Error model / seeding.
  int n_stat = 1;   // no of seeds
  int n_meth = 0;   // error model: 0 = standard, 1 = cormisal (girder)
  int n_scale = 1;  // no of steps the rms errors are ramped over
  bool bba = false; // beam-based alignment (code-only, no keyword)

  // Orbit correction.
  std::string loc_Fam_name = "";
  int n_cell = -1, n_thread = -1, n_orbit = 5;
  double h_maxkick = 1.0e-3, v_maxkick = 1.0e-3;
  double h_cut = 1.0e-4, v_cut = 1.0e-4;
  std::vector<std::string> bpm_Fam_names, corr_Fam_names[2];

  // Coupling / vertical-dispersion correction (the LOCO off-diagonal block).
  int n_lin = 3;
  int SQ_per_scell = 1, BPM_per_scell = 10;
  int HCM_per_scell = 10, VCM_per_scell = 10;
  double kick = 0.01e-3; // trim kick used to measure the response matrix
  double VDweight = 1e3, HVweight = 1e0, VHweight = 1e0;
  double disp_wave_y = 0e0, disp_wave_o = 0e0, qt_s_cut = 1e0;
  int qt_from_file = 0;

  // ID (insertion-device) correction.
  int N_calls = 0, N_steps = 0;
  int N_Fam = 0, Q_Fam[N_Fam_max];
  double ID_s_cut = 1e1;

  // Target tunes / chromaticities (fitted in err_and_corr_init when set).
  double TuneX = 0e0, TuneY = 0e0;
  double ChromX = 1e6, ChromY = 1e6;

  // Fit-knob families (matched by element-name prefix). Default to the historical
  // SLS-2 tune quads / chroma sextupoles so pre-keyword param files are a drop-in
  // (they behave exactly as before); newer lattices (e.g. m4U: q1_n1/q2_n1,
  // s2_n1/s4_n1) name their families via the tune_fams / chrom_fams keywords.
  std::string tune_fam[2]  = {"qax", "qay"};
  std::string chrom_fam[2] = {"sf", "sd"};

  // Dynamic aperture.
  int n_track_DA = 512, n_aper_DA = 15, n_delta_DA = 12;
  double delta_DA = 3e-2;
  bool DA_bare = false;

  // Frequency map. NTURN is set to 10000 (2*NTURN for diffusion) in naffutils.h.
  bool freq_map = false;
  int n_x = 50, n_y = 30, n_dp = 25, n_tr = 2064;
  double x_max_FMA = 20e-3, y_max_FMA = 6e-3, delta_FMA = 3e-2;

  // Parse param.dat. NOT a pure parser — it also loads the lattice
  // (Read_Lattice / rdmfile / rdmfile_at), seeds the RNG cut (setrancut), and
  // resolves family names to indices in globval (gs/ge/bpm/hcorr/vcorr/qt).
  // Those side effects on global state are kept exactly as they were; untangling
  // them is a backlog item, not Phase 1.
  void get_param(const std::string &param_file);

  // Project the knobs each corrector needs. Same shape for both: the corrector
  // takes its config by value and never sees the rest of the bag.
  corr::orbit_cfg orbit_config(void) const;
  corr::coupling_cfg coupling_config(void) const;

  // ------------------------------------------------------------------

  // Sextupole b_3 save buffer, shared by the paired zero_mult/restore_mult
  // façades (ctrl_cod.cc calls them as a pair on the same instance).
  std::vector<double> bn_an[2 * HOMmax + 1];

  // Bare-lattice reference optics at the sextupoles — measured, not a parsed
  // knob, so it is a corr:: slice rather than a get_param field. Filled by
  // get_bare().
  corr::bare_optics bare;

  // The correctors. Each owns its own working state.
  corr::id_corr id;           // ID (insertion-device) optics correction
  corr::girder_model girders; // cormisal girder error model (n_meth == 1)
  corr::coupling_corr skew;   // coupling / vertical dispersion (LOCO off-diag)

  //-------------------------------------------------------------------
  // Delegators. Each forwards to the correction/ module named in its body; they
  // exist only so the callers above do not have to move yet.

  // Girder error model.
  void GirderSetup();
  void SetCorMis(double gxrms, double gyrms, double gtrms, double jxrms,
                 double jyrms, double exrms, double eyrms, double etrms,
                 double rancutx, double rancuty, double rancutt, long iseed);
  void CorMis_in(double *gdxrms, double *gdzrms, double *gdarms,
                 double *jdxrms, double *jdzrms, double *edxrms,
                 double *edzrms, double *edarms, double *bdxrms,
                 double *bdzrms, double *bdarms, double *rancutx,
                 double *rancuty, double *rancutt, long *iseed, long *iseednr);

  // Bare-lattice reference.
  void get_bare(void);
  void get_dbeta_dnu(double m_dbeta[], double s_dbeta[], double m_dnu[],
                     double s_dnu[]);

  // Coupling / vertical beam size.
  void ini_skew_cor(const double deta_y_max, const double deta_y_offset);
  void corr_eps_y(const int cnt);

  // ID correction.
  void reset_quads(void);
  void ini_ID_corr(const bool IDs);
  bool ID_corr(const int N_calls, const int N_steps, const bool IDs,
               const int cnt);

  // Error model.
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

  // Orbit correction.
  bool CorrectCOD_N(const int n_orbit, const int k);
  void ini_COD_corr(const int n_bpm_Fam, const std::string bpm_names[],
                    const int n_hcorr_Fam, const std::string hcorr_names[],
                    const int n_vcorr_Fam, const std::string vcorr_names[],
                    const bool svd);

  bool cod_corr(const int n_cell, const double scl, const double h_maxkick,
                const double v_maxkick, orb_corr_type orb_corr[]);

  void Orb_and_Trim_Stat(orb_corr_type orb_corr[]);

  void prt_cod_corr_lat(void);

  // Driver.
  void err_and_corr_init(const string &param_file, orb_corr_type orb_corr[]);

  void err_and_corr_exit(orb_corr_type orb_corr[]);
};

void get_bn2(const string file_name1, const string file_name2, int n,
             const bool prt);

#endif
