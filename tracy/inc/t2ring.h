/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -                                         */

#ifndef T2RING_H
#define T2RING_H

// Computation result files
const char beam_envelope_file[] = "beam_envelope";


void GetNu(Vector2 &nu, arma::mat &M);

bool Cell_GetABGN(arma::mat &M, Vector2 &alpha, Vector2 &beta, Vector2 &gamma,
		  Vector2 &nu);

void Cell_Geteta(long i0, long i1, bool ring, double dP);

void Cell_Twiss(long i0, long i1, ss_vect<tps> &Ascr, bool chroma, bool ring,
		double dP);

void Ring_Getchrom(double dP);

void Ring_GetTwiss(bool chroma, double dP);

/* void Ring_Fittune(Vector2 &nu, double eps, iVector2 &nq, long qf[], long qd[], */
/*                   double dkL, long imax); */

/* void Ring_Fitchrom(Vector2 &ksi, double eps, iVector2 &ns, */
/* 		   long sf[], long sd[], double dkpL, long imax); */

/* void Ring_FitDisp(long pos, double eta, double eps, long nq, long q[], */
/*                   double dkL, long imax); */

void get_dI_eta_5(const int k, std::vector<ElemType*> Elem);

double get_code(const ConfigType &conf, const ElemType &Cell);

void Cell_Twiss(const long int i0, const long int i1);

void findcod(LatticeType &lat, double dP);

void prt_lin_map(const int n_DOF, const ss_vect<tps> &map);

ss_vect<tps> get_A(const std::vector<double> &alpha,
		   const std::vector<double> &beta,
		   const std::vector<double> &eta,
		   const std::vector<double> &etap);

ss_vect<tps> get_A_CS(const int n, const ss_vect<tps> &A,
		      std::vector<double> &dnu);

void get_ab(const ss_vect<tps> &A, std::vector<double> &alpha,
	    std::vector<double> &beta, std::vector<double> &dnu,
	    std::vector<double> &eta, std::vector<double> &etap);

void get_twoJ(const int n_DOF, const ss_vect<double> &ps,
	      const ss_vect<tps> &A, double twoJ[]);

void Trac(LatticeType &lat, double x, double px, double y, double py, double dp,
	  double ctau, long nmax, long pos, long &lastn, long &lastpos,
	  FILE *outf1);

void SetKLpar(LatticeType &lat, int Fnum, int Knum, int Order, double kL);

double GetKpar(LatticeType &lat, int Fnum, int Knum, int Order);

double get_dynap(LatticeType &lat, const double delta, const int n_aper,
		 const int n_track, const bool cod);

// Same C asctime \n.
char *asctime2(const struct tm *timeptr);
struct tm* GetTime();

void computeFandJ(LatticeType &lat, int n, double *x, ss_vect<double> *fjac,
		  double *fvect);
int Newton_Raphson(LatticeType &lat, int n, ss_vect<double> &x, int ntrial,
		   double tolx);

void get_dnu(const int n, const ss_vect<tps> &A, std::vector<double> &dnu);

ss_vect<tps> tp_S(const int n_DOF, const ss_vect<tps> &A);

void dynap(FILE *fp, LatticeType &lat, double r, const double delta,
	   const double eps, const int npoint, const int nturn,double x[],
	   double y[], const bool floqs, const bool cod, const bool print);

double get_aper(int n, double x[], double y[]);

ss_vect<tps> get_S(const int n_DOF);

void getdynap(LatticeType &lat, double &r, double phi, double delta, double eps,
	      int nturn, bool floqs);

#endif
