#ifndef PARAM_H
#define PARAM_H

// N_Fam_max lives in correction/corr_config.h; n_b2_max/n_b3_max/max_ID_Fams and
// the ID-correction weights in correction/id_corr.h; max_bpm/max_corr and the
// skew/eta_y output file names in correction/loco/coupling_corr.h. All are
// included before this header.

// Computation result files
const char beam_envelope_file[] = "beam_envelope";

// Lattice error and correction files
const char CodCorLatFileName[]  = "codcorlat.out";

// The transitional façade.
//
// Every param.dat knob (and get_param itself) now lives in corr::config_data,
// which this class DERIVES FROM rather than holds: the apps and tracy/src/dynap.cc
// read the knobs through the instance (params.n_cell, params.h_maxkick, ...) in
// ~90 places, and inheriting keeps every one of them compiling untouched while
// the state itself has already moved. That is scaffolding, not design — step 4
// of the refactor repoints those callers at corr::config_data directly (a
// param_data_type& already binds to a const config_data&), after which step 5
// deletes this class. See correction_refactor.md.
//
// What is left here that is NOT inherited config: the four corrector objects,
// the bare-lattice reference, one save buffer, and a wall of delegators.

class param_data_type : public corr::config_data {
 private:

 public:
  // Sextupole b_3 save buffer, shared by the paired zero_mult/restore_mult
  // façades (ctrl_cod.cc calls them as a pair on the same instance).
  std::vector<double> bn_an[2*HOMmax+1];

  // Bare-lattice reference optics at the sextupoles — measured, not config, so
  // it is NOT part of corr::config_data. Filled by get_bare().
  corr::bare_optics bare;

  // The correctors. Each owns its own working state.
  corr::id_corr       id;      // ID (insertion-device) optics correction
  corr::girder_model  girders; // cormisal girder error model (n_meth == 1)
  corr::coupling_corr skew;    // coupling / vertical dispersion (LOCO off-diag)

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
