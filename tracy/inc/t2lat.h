/* Tracy-2:

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version,
                 SLS, PSI      1995 - 1997,
   M. Boege      SLS, PSI      1998          Pascal to C translation,
   L. Nadolski   SOLEIL        2002          Link to NAFF and Radia field maps,
   J. Bengtsson  NSLS-II, BNL  2004 - 2015,
   J. Bengtsson                2017          C to C++ re-factoring.           */


// maximum number of LEGO blocks (Cell_nLoc)
#define Cell_nLocMax 20000

// maximum number of families for Elem_NFam
#define Elem_nFamMax 3000


class Lattice_Type {
 private:
 public:
  ElemFamType ElemFam[Elem_nFamMax];
  CellType    Cell[Cell_nLocMax+1];

  bool Lattice_Read(FILE **fi_, FILE **fo_);
  void Read_Lattice(const char *fic);

  void rdmfile(const char *mfile_dat);
  void prtmfile(const char mfile_dat[]);

  long Elem_Index(const std::string &name1);

  // From t2ring.cc.

  void GetNu(Vector2 &nu, Matrix &M);

  void Cell_Twiss(long i0, long i1, ss_vect<tps> &Ascr, bool chroma, bool ring,
		  double dP);

  void Cell_GetABGN(Matrix &M,
		    Vector2 &alpha, Vector2 &beta, Vector2 &gamma, Vector2 &nu);

  void Cell_Geteta(long i0, long i1, bool ring, double dP);

  void Ring_Getchrom(double dP);

  void Ring_Twiss(bool chroma, double dP);

  void Ring_GetTwiss(bool chroma, double dP);

  void Ring_Fittune(Vector2 &nu, double eps, iVector2 &nq, long qf[], long qd[],
		    double dkL, long imax);

  void Ring_Fitchrom(Vector2 &ksi, double eps, iVector2 &ns,
		     long sf[], long sd[], double dkpL, long imax);

  void Ring_FitDisp(long pos, double eta, double eps, long nq, long q[],
		    double dkL, long imax);

  // From physlib.cc.

  bool getcod(double dP, long &lastpos);

  void get_twiss3(long int loc,
		  Vector2 alpha[], Vector2 beta[], Vector2 nu[],
		  Vector2 eta[], Vector2 etap[]);

  double int_curly_H1(long int n);

  void getabn(Vector2 &alpha, Vector2 &beta, Vector2 &nu);

  void TraceABN(long i0, long i1, const Vector2 &alpha, const Vector2 &beta,
		const Vector2 &eta, const Vector2 &etap, const double dP);

  void FitTune(long qf, long qd, double nux, double nuy);

  void FitChrom(long sf, long sd, double ksix, double ksiy);

  void FitDisp(long q, long pos, double eta);

  void getfloqs(psVector &x);

  void track(const char* file_name,
	     double ic1, double ic2, double ic3, double ic4, double dp,
	     long int nmax, long int &lastn, long int &lastpos, int floqs,
	     double f_rf);

  void track_(double r, struct LOC_getdynap *LINK);

  void getdynap(double &r, double phi, double delta, double eps,
		int nturn, bool floqs);

  void getcsAscr(void);

  void dynap(FILE *fp, double r, const double delta,
	     const double eps, const int npoint, const int nturn,
	     double x[], double y[], const bool floqs, const bool cod,
	     const bool print);

  void ttwiss(const Vector2 &alpha, const Vector2 &beta,
	      const Vector2 &eta, const Vector2 &etap, const double dP);
};

