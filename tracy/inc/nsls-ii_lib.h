// global params

#ifndef NSLS_II_LIB_H
#define NSLS_II_LIB_H

const int           max_elem = Cell_nLocMax;

extern ss_vect<tps> map;
extern MNF_struct   MNF;

extern double       chi_m;


void lwr_case(char str[]);

void upr_case(char str[]);

//void prt_trace (void);

void file_rd(std::ifstream &inf, const string &file_name);

void file_wr(std::ofstream &outf, const string &file_name);

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

void get_I(double I[], const bool prt);

void get_eps_x(double &eps_x, double &sigma_delta, double &U_0, double J[]);

void GetEmittance(const int Fnum, const bool prt);

void prt_lat(const int loc1, const int loc2, const char *fname, const int Fnum,
	     const bool all);

void prt_lat(const char *fname, const int Fnum, const bool all);

void Cell_Twiss(const long int i0, const long int i1);

void prt_lat(const int loc1, const int loc2, const char *fname, const int Fnum,
	     const bool all, const int n);

void prt_lat(const char *fname, const int Fnum, const bool all, const int n);

void prt_chrom_lat(void);

void prt_cod(const char *file_name, const int Fnum, const bool all);

void printcod(const char *file_name);

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

void set_aper_elem(const int Fnum, const int Knum,
		   const double Dxmin, const double Dxmax,
		   const double Dymin, const double Dymax);

void set_aper_fam(const int Fnum,
		  const double Dxmin, const double Dxmax,
		  const double Dymin, const double Dymax);

void set_aper_type(const int type, const double Dxmin, const double Dxmax,
		   const double Dymin, const double Dymax);

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

bool CorrectCOD(const int n_orbit, const double scl);

void prt_beamsizes(const int cnt);

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

double get_dynap(const double delta, const int n_aper, const int n_track,
		 const bool cod);

double get_chi2(long int n, double x[], double y[], long int m, psVector b);

void pol_fit(int n, double x[], double y[], int order, psVector &b,
	     double &sigma, const bool prt);

void get_ksi2(const double d_delta);

bool find_nu(const int n, const double nus[], const double eps, double &nu);

bool get_nu(const double Ax, const double Ay, const double delta,
	    double &nu_x, double &nu_y);

void dnu_dA(const double Ax_max, const double Ay_max, const double delta,
	    const int n_ampl);

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

void get_map_twiss(const ss_vect<tps> &M,
		   double beta0[], double beta1[], double nu[], bool stable[]);

void set_map(MapType *Map);

void set_map(const int Fnum, const double dnu[]);

void set_map_per(MapType *Map,
		 const double alpha0[], const double beta0[],
		 const double eta0[], const double etap0[]);

void set_map_per(const int Fnum,
		 const double alpha0[], const double beta0[],
		 const double eta0[], const double etap0[]);

void set_map_reversal(CellType &Cell);

void set_map_reversal(const long int Fnum);

void setmp(long ilat, long m, long n, double rr, double bnoff, double cmn);

void setmpall (double rref);

#endif
