/* Tracy-3

   Johan Bengtsson.

   NO   1   link to tlinear TPSA (nv_tps = 1)
        >1  link to Fortran-77 TPSA                                           */


// C standard library.
#include <stdio.h>
#include <stddef.h>
#include <setjmp.h>
#include <time.h>
#include <memory.h>
#include <malloc.h>
//#include <execinfo.h>


// C++ standard library.
#include <cstdlib>
#include <cfloat>
#include <cctype>
#include <cmath>
#include <complex>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>

//#undef max
#include <vector>

using namespace std;

// Tracy-3
#include "field.h"
#include "mathlib.h"

#if NO == 1
  // linear TPSA.
  #include "tpsa_lin.h"
  #include "tpsa_lin_pm.h"
#else
  // interface to Fortran-77 TPSA.
  #include "tpsa_for.h"
  #include "tpsa_for_pm.h"
#endif

#include "tracy.h"
#include "tracy_global.h"
#include "ety.h"
#include "eigenv.h"

#include "num_rec.h"

#include "radia2tracy.h"
#include "pascalio.h"

#include "t2elem.h"
#include "t2cell.h"
#include "t2lat.h"
#include "t2ring.h"

#include "fft.h"

//#include "physlib.h"
//#include "nsls-ii_lib.h"
#include "orb_corr.h"
#include "param.h"
#include "dynap.h"

#include "lsoc.h"

#include "modnaff.h"

#include "naffutils.h"
#include "complexeheader_naff.h"

//#include "soleillib.h"


extern Lattice_Type Lattice;

// Truncated Power Series Algebra (TPSA).
extern const int  nv_tps, nd_tps, iref_tps;
extern int        no_tps, ndpt_tps;
extern double     eps_tps;

extern ElemFamType ElemFam[];

extern CellType Cell[];

extern globvalrec globval;


// From nsls-ii_lib.h.

// Global parameters.

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

double f_IBS(const double chi_m);

double get_int_IBS(void);

void rm_space(char *name);

void get_bn(const char file_name[], int n, const bool prt);

double get_chi2(long int n, double x[], double y[], long int m, psVector b);

void pol_fit(int n, double x[], double y[], int order, psVector &b,
	     double &sigma, const bool prt);

bool find_nu(const int n, const double nus[], const double eps, double &nu);

bool get_nu(const double Ax, const double Ay, const double delta,
	    double &nu_x, double &nu_y);

double get_code(CellType &Cell);

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


// from physlib.h.

// For tune fitting.
#define nueps    1e-6 // Precision.
#define nudkL    0.01 // Step.
#define nuimax   10   // Maximum number of iterations.

// For chromaticity fitting.
#define ksieps   1e-5
#define ksidkpL  0.01
#define ksiimax  10

//* For dispersion fitting.
#define dispeps  1e-10
#define dispdkL  0.001
#define dispimax 10
#define npeakmax 10

// Dynamical aperture.
#define px_0      0.0
#define py_0      0.0

// 80% sigma coupling.
#define sigma_eps sqrt((25.0/16.0-1.0)/(25.0/16.0+1.0))

#define writetrack      true   /*protocol from tracking*/

// getfloq.
#define nfloq     4

// inibump.
#define dnux      0.02
#define dnuy      0.01

// TraceABN.
#define ntrace    4

typedef long ipeakbuf[npeakmax];
typedef double peakbuf[npeakmax];

void rm_mean(long int n, double x[]);

void printglob(void);

void PrintMat(long n, psVector *A);

void PrintVec(long n, double *X);

void recalc_S();

double Circumference(void);

void GetMean(long n, double *x);

void prt_sigma(void);

struct LOC_getdynap {
  double phi, delta;
  long   nturn;
  bool   floqs, lost;
} ;

double get_aper(int n, double x[], double y[]);

void GetTrack(const char *file_name,
	      long *n, double *x, double *px, double *y, double *py);

void Getj(long n, double *x, double *px, double *y, double *py);

double Fract(double x);

double GetArg(double x, double px, double nu);

void GetPhi(long n, double *x, double *px, double *y, double *py);

void Sinfft(int n, double *xr);

void sin_FFT(int n, double xr[]);

void sin_FFT(int n, double xr[], double xi[]);

void GetInd(int n, int k, int *ind1, int *ind3);

void GetInd1(int n, int k, int *ind1, int *ind3);

void GetPeak(int n, double *x, int *k);

void GetPeak1(int n, double *x, int *k);

double Int2snu(int n, double *x, int k);

double Sinc(double omega);

double intsampl(int n, double *x, double nu, int k);

double linint(int n, int k, double nu, double *x);


extern double  FindRes_eps;

struct LOC_findres {
  int n;
  double nux, nuy, f;
  int *nx, *ny;
  double eps;
  bool found;
} ;

void FndRes(struct LOC_findres *LINK);

void FindRes(int n_, double nux_, double nuy_, double f_, int *nx_, int *ny_);

void GetPeaks(int n, double *x, int nf, double *nu, double *A);

void GetPeaks1(int n, double *x, int nf, double *nu, double *A);

void SetTol(int Fnum, double dxrms, double dyrms, double drrms);

void Scale_Tol(int Fnum, double dxrms, double dyrms, double drrms);

void SetaTol(int Fnum, int Knum, double dx, double dy, double dr);

void ini_aper(const double Dxmin, const double Dxmax, 
              const double Dymin, const double Dymax);

void set_aper(const int Fnum, const double Dxmin, const double Dxmax,
		     const double Dymin, const double Dymax);

void LoadApertures(const char *ChamberFileName);

void LoadTolerances(const char *TolFileName);

void ScaleTolerances(const char *TolFileName, const double scl);

void SetKpar(int Fnum, int Knum, int Order, double k);

void SetdKpar(int Fnum, int Knum, int Order, double k);

void SetL(int Fnum, int Knum, double L);

void SetL(int Fnum, double L);

void SetKLpar(int Fnum, int Knum, int Order, double kL);

void SetdKLpar(int Fnum, int Knum, int Order, double dkL);

void SetdKrpar(int Fnum, int Knum, int Order, double dkrel);

void Setbn(int Fnum, int order, double bn);

void SetbnL(int Fnum, int order, double bnL);

void Setdbn(int Fnum, int order, double dbn);

void SetdbnL(int Fnum, int order, double dbnL);

void Setbnr(int Fnum, int order, double bnr);

void SetbnL_sys(int Fnum, int Order, double bnL_sys);

void set_dbn_rel(const int type, const int n, const double dbn_rel);

double GetKpar(int Fnum, int Knum, int Order);

double GetL(int Fnum, int Knum);

double GetKLpar(int Fnum, int Knum, int Order);

void SetdKLsys(int Fnum, int Order, double dkLsys);

void SetdKLrms(int Fnum, int Order, double dkLrms);

void Setdkrrms(int Fnum, int Order, double dkrrms);

void SetKL(int Fnum, int Order);

void set_dx(const int type, const double sigma_x, const double sigma_y);

void SetBpmdS(int Fnum, double dxrms, double dyrms);

void codstat(double *mean, double *sigma, double *xmax, long lastpos,
		    bool all);

void CodStatBpm(double *mean, double *sigma, double *xmax, long lastpos,
                long bpmdis[]);
                
double Sgn (double x);

double digitize(double x, double maxkick, double maxsamp);

//svdarray xmemo[2];

double digitize2(long plane, long inum, double x, double maxkick,
			double maxsamp);

void Dis_In(long *bpmdis, long *vcorrdis, long *hcorrdis,
                   long *wvdis, long *whdis);


/* high level functions for reading lattice file */
long get_bpm_number(void);
long get_hcorr_number(void);
long get_vcorr_number(void);
long get_qt_number(void);


/* tracking */
void GetChromTrac(long Nb, long Nbtour, double emax, double *xix, double *xiz);

void GetTuneTrac(long Nbtour, double emax, double *nux, double *nuz);

void Trac(double x, double px, double y, double py, double dp, double ctau,
                 long nmax, long pos, long &lastn, long &lastpos, FILE *outf1);

/* close orbit */
// simple precision

void findcodS(double dP);

void computeFandJS(double *x, int n, double **fjac, double *fvect);

void Newton_RaphsonS(int ntrial, double x[], int n, double tolx);

// double precision

void findcod(double dP);

void computeFandJ(int n, double *x, psVector *fjac, double *fvect);

int Newton_Raphson(int n, psVector &x, int ntrial, double tolx);

/* Vacuum chamber */

void PrintCh(void);

void ChamberOff(void);


// From soleillib.h.

/**** Protypes ****/

void SetErr(void);

void InducedAmplitude(long spos);

void Hfonction(long pos, double dP, Vector2 H);

void Hcofonction(long pos, double dP, Vector2 H);

void Get_Disp_dp(void);

void read_corrh(void);

void set_vectorcod(psVector codvector[], double dP);

void SetDecapole(void);

/* Tracking */
void Phase(double x, double xp, double y, double yp, double energy, double ctau,
	   long Nbtour);

void Phase2(long pos, double x, double xp, double y, double yp, double energy,
	    double ctau, long Nbtour);

void PhasePoly(long pos, double x0, double px0, double z0, double pz0,
	       double delta0, double ctau0, long Nbtour);

void Check_Trac(double x, double px, double y, double py, double dp);

void PhasePortrait(double x0, double px0, double z0, double pz0,
		   double delta0, double ctau, double end, long Nb,
		   long Nbtour, int num);

void PhasePortrait2(long pos, double x0, double px0, double z0, double pz0,
		    double delta0, double ctau, double end, long Nb,
		    long Nbtour, int num);

void Multipole(void);

void MomentumAcceptance(long deb, long fin, double ep_min, double ep_max,
			long nstepp, double em_min, double em_max, long nstepm);

void Trac_Tab(double x, double px, double y, double py, double dp, long nmax,
	      long pos, long *lastn, long *lastpos, FILE *outf1,
	      double Tx[][NTURN]);

void SetSkewQuad(void);

void TracCO(double x, double px, double y, double py, double dp, double ctau,
	    long nmax, long pos, long *lastn, long *lastpos, FILE *outf1);

void Dyna(long Nbx, long Nbz, long Nbtour, double xmax, double zmax,
	  double energy, bool diffusion);

/* Frequency map analysis */

void NuDp(long Nb, long Nbtour, double emax);

void fmap(long Nbx, long Nbz, long Nbtour, double xmax, double zmax,
	  double energy, bool diffusion, bool matlab);

void fmapdp(long Nbx, long Nbe, long Nbtour, double xmax, double emax, double z,
	    bool diffusion, bool matlab);

void Nu_Naff(void);

void NuDx(long Nbx, long Nbz, long Nbtour, double xmax, double ymax,
	  double energy);

/* Vacuum chamber */

void DefineCh(void);

void Enveloppe(double x, double px, double y, double py, double dp,
	       double nturn);

void ChamberOn(void);

/* Longitudinal Hamiltonian*/

void PhaseLongitudinalHamiltonien(void);

void PassA(double *phi, double delta0, double step);

void PassB(double phi0, double *delta, double step);

double Hsynchrotron(double phi, double delta);

/* Miscelleneous */ 

void Enveloppe2(double x, double px, double y, double py, double dp,
		double nturn);

void Phase3(long pos, double x, double px, double y, double py, double energy,
            double ctau, long Nbtour);

double EnergySmall(double *X, double irho);

double EnergyDrift(double *X);

void getA4antidamping();

void fmapfull(long Nbx, long Nbz, long Nbtour, double xmax, double zmax,
              double energy, bool diffusion);

void spectrum(long Nbx, long Nbz, long Nbtour, double xmax, double zmax,
              double energy, bool diffusion);
