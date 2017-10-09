

class Lattice_Type {
 private:
 public:
  Lattice_Param param;
  ElemFamType   ElemFam[Elem_nFamMax];
  CellType      Cell[Cell_nLocMax+1];

  bool Lattice_Read(FILE **fi_, FILE **fo_);
  void Read_Lattice(const char *fic);

  void rdmfile(const char *mfile_dat);
  void prtmfile(const char mfile_dat[]);

  long Elem_Index(const std::string &name1);

  int GetnKid(const int Fnum1);

  long Elem_GetPos(const int Fnum1, const int Knum1);

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

  // From nsls-ii_lib.cc.

  double get_eps_x(void);

  void GetEmittance(const int Fnum, const bool prt);

  void prt_lat(const int loc1, const int loc2, const char *fname,
	       const int Fnum, const bool all);

  void prt_lat(const char *fname, const int Fnum, const bool all);

  void Cell_Twiss(const long int i0, const long int i1);

  void prt_lat(const int loc1, const int loc2, const char *fname,
	       const int Fnum, const bool all, const int n);

  void prt_lat(const char *fname, const int Fnum, const bool all, const int n);

  void prt_chrom_lat(void);

  void prt_cod(const char *file_name, const int Fnum, const bool all);

  void prt_beampos(const char *file_name);

  void CheckAlignTol(const char *OutputFile);

  bool CorrectCOD(const int n_orbit, const double scl);

  void prt_beamsizes();

  double Touschek(const double Qb, const double delta_RF,
		  const double eps_x, const double eps_y,
		  const double sigma_delta, const double sigma_s);

  double Touschek(const double Qb, const double delta_RF,const bool consistent,
		  const double eps_x, const double eps_y,
		  const double sigma_delta, double sigma_s,
		  const int n_turn, const bool aper_on,
		  double sum_delta[][2], double sum2_delta[][2]);

  void IBS(const double Qb, const double eps_SR[], double eps[],
	   const bool prt1, const bool prt2);

  void IBS_BM(const double Qb, const double eps_SR[], double eps[],
	      const bool prt1, const bool prt2);

  double get_dynap(const double delta, const int n_aper, const int n_track,
		   const bool cod);

  void get_ksi2(const double d_delta);

  void dnu_dA(const double Ax_max, const double Ay_max, const double delta,
	      const int n_ampl);

  bool orb_corr(const int n_orbit);

  void get_alphac(void);

  void get_alphac2(void);
};

