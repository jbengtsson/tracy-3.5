/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

*/

/* For tune fitting */
#define nueps           1e-6      //precision
#define nudkL           0.01   //step
#define nuimax          10        // maximum number of iterations

/* For chromaticity fitting */
#define ksieps          1e-5
#define ksidkpL         0.01
#define ksiimax         10

/* For dispersion fitting */
#define dispeps         1e-10
#define dispdkL         0.001
#define dispimax        10
#define npeakmax        10

// Dynamical aperture
#define px_0      0.0
#define py_0      0.0

/* 80% sigma coupling */
#define sigma_eps        sqrt((25.0/16.0-1.0)/(25.0/16.0+1.0))

#define writetrack      true   /*protocol from tracking*/

// getfloq
#define nfloq     4

// inibump
#define dnux      0.02
#define dnuy      0.01

// TraceABN
#define ntrace    4

typedef long ipeakbuf[npeakmax];
typedef double peakbuf[npeakmax];

double int_curly_H(long int n);

void rm_mean(long int n, double x[]);

void printglob(void);

void printlatt(void);

void PrintMat(long n, psVector *A);

void PrintVec(long n, double *X);

void recalc_S();

double Circumference(void);

void GetMean(long n, double *x);

bool getcod(double dP, long &lastpos);

void get_twiss3(long int loc,
                       Vector2 alpha[], Vector2 beta[], Vector2 nu[],
                       Vector2 eta[], Vector2 etap[]);

void getabn(double *alpha, double *beta, double *nu);

void TraceABN(long i0, long i1, const Vector2 &alpha, const Vector2 &beta,
	      const Vector2 &eta, const Vector2 &etap, const double dP);

void ttwiss(const Vector2 &alpha, const Vector2 &beta,
                   const Vector2 &eta, const Vector2 &etap, const double dP);

void prt_sigma(void);

/* 2 parameter fitting routines */
void FitTune(long qf, long qd, double nux, double nuy);

void FitChrom(long sf, long sd, double ksix, double ksiy);

void FitDisp(long q, long  pos, double eta);

void inibump(long coh, long cov);

void getfloqs(psVector &x);

void track(const char* file_name,
	   double ic1, double ic2, double ic3, double ic4, double dp,
	   long int nmax, long int &lastn, long int &lastpos, int floqs,
	   double f_rf);

struct LOC_getdynap {
  double phi, delta;
  long nturn;
  bool floqs, lost;
} ;

void track_(double r, struct LOC_getdynap *LINK);

void getdynap(double &r, double phi, double delta, double eps,
	      int nturn, bool floqs);

void getcsAscr(void);

void dynap(FILE *fp, double r, const double delta,
	   const double eps, const int npoint, const int nturn,
       double x[], double y[], const bool floqs, const bool cod, const bool print);

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

void FindRes(int n_, double nux_, double nuy_, double f_,
		    int *nx_, int *ny_);

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


/* high level functions for reading lattice file*/  
void Read_Lattice(const char *fic);
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

