// global params

#ifndef NSLS_II_LIB_H
#define NSLS_II_LIB_H


void lwr_case(char str[]);

void upr_case(char str[]);

//void prt_trace (void);

void chk_cod(const bool cod, const char *proc_name);

void no_sxt(LatticeType &lat);

void get_map(LatticeType &lat, const bool cod);

tps get_h(void);

void get_m2(const ss_vect<tps> &ps, tps m2[]);

double get_curly_H(const double alpha_x, const double beta_x,
		   const double eta_x, const double etap_x);

void printcod(LatticeType &lat, const char *file_name);


//------------------------------------------------------------------------------

double get_Wiggler_BoBrho(LatticeType &lat, const int Fnum, const int Knum);

void set_Wiggler_BoBrho(LatticeType &lat, const int Fnum, const int Knum,
			const double BoBrhoV);
void set_Wiggler_BoBrho(LatticeType &lat, const int Fnum, const double BoBrhoV);

void set_ID_scl(LatticeType &lat, const int Fnum, const int Knum,
		const double scl);

void SetFieldValues_fam(LatticeType &lat, const int Fnum, const bool rms,
			const double r0, const int n, const double Bn,
			const double An, const bool new_rnd);

void SetFieldValues_type(LatticeType &lat, const int N, const bool rms,
			 const double r0, const int n, const double Bn,
			 const double An, const bool new_rnd);

void SetFieldErrors(LatticeType &lat, const char *name, const bool rms,
		    const double r0, const int n, const double Bn,
		    const double An, const bool new_rnd);

bool CorrectCOD(LatticeType &lat, const int n_orbit, const double scl);

double Touschek(LatticeType &lat, const double Qb, const double delta_RF,
		const double eps_x, const double eps_y,
		const double sigma_delta, const double sigma_s);

double Touschek(LatticeType &lat, const double Qb, const double delta_RF,
		const bool consistent, const double eps_x, const double eps_y,
		const double sigma_delta, double sigma_s, const int n_turn,
		const bool aper_on, double sum_delta[][2],
		double sum2_delta[][2]);

double f_IBS(const double chi_m);

double get_int_IBS(void);

void IBS(LatticeType &lat, const double Qb, const double eps_SR[], double eps[],
	 const bool prt1, const bool prt2);

void IBS_BM(LatticeType &lat, const double Qb, const double eps_SR[],
	    double eps[], const bool prt1, const bool prt2);

void rm_space(char *name);

void get_bn(LatticeType &lat, const char file_name[], int n, const bool prt);

double get_chi2(long int n, double x[], double y[], long int m, psVector b);

void pol_fit(int n, double x[], double y[], int order, psVector &b,
	     double &sigma, const bool prt);

void get_ksi2(LatticeType &lat, const double d_delta);

bool find_nu(const int n, const double nus[], const double eps, double &nu);

bool get_nu(LatticeType &lat, const double Ax, const double Ay,
	    const double delta, double &nu_x, double &nu_y);

void dnu_dA(LatticeType &lat, const double Ax_max, const double Ay_max,
	    const double delta, const int n_ampl);

bool orb_corr(LatticeType &lat, const int n_orbit);

void get_alphac(LatticeType &lat);

void get_alphac2(LatticeType &lat);

/* void bend_cal_Fam(const int Fnum); */

/* void bend_cal(void); */

double h_ijklm(const tps &h, const int i, const int j, const int k,
	       const int l, const int m);

void set_tune(LatticeType &lat, const char file_name1[],
	      const char file_name2[], const int n);

void prt_H_long(const int n, const double phi_max, const double delta_max,
		const double U0);

void get_map_twiss(const ss_vect<tps> &M,
		   double beta0[], double beta1[], double nu[], bool stable[]);

void set_map(MapType *Map);

void set_map(LatticeType &lat, const int Fnum, const double dnu[]);

void set_map_per(MapType *Map,
		 const double alpha0[], const double beta0[],
		 const double eta0[], const double etap0[]);

void set_map_per(LatticeType &lat, const int Fnum, const double alpha0[],
		 const double beta0[], const double eta0[],
		 const double etap0[]);

void set_map_reversal(LatticeType &lat, CellType &Cell);

void set_map_reversal(LatticeType &lat, const long int Fnum);

void setmp(long ilat, long m, long n, double rr, double bnoff, double cmn);

void setmpall(LatticeType &lat, double rref);

#endif
