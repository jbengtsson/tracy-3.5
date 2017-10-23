/* Tracy-2:

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version,
   SLS, PSI      1995 - 1997,
   M. Boege      SLS, PSI      1998          Pascal to C translation,
   L. Nadolski   SOLEIL        2002          Link to NAFF and Radia field maps,
   J. Bengtsson  NSLS-II, BNL  2004 - 2015,
   J. Bengtsson                2017          Re-factored from C to C++.       */

// LEGO structure for lattice.
//
// LEGO:
//   Leg Godt (play well); invented by Master Cabinet Maker Ole Kirk
//   Christiansen, 1934.

#define Cell_nLocMax 20000 // Max number of LEGO blocks (Cell_nLoc).
#define Elem_nFamMax 3000  // Max number of families for Elem_NFam.
#define nKidMax      5000  // Max number of kids.

// ID Laurent.
#define IDXMAX       200
#define IDZMAX       100

/* definitions for old soleilcommon.h */
#define NTURN        10000 // 2*NTURN for diffusion
#define DIM          6

const int max_elem     = Cell_nLocMax;
const int n_harm_max   = 10;
const int Spreader_max = 10;


#define DOF          (ss_dim/2)
#define nv_          6
#define maxtransfmat 2000
#define maxkicks     5

#define M_PI         3.14159265358979323846    // pi.

#define debug true

const double c0    = 2.99792458e8;             // speed of light in vacuum.
const double q_e   = 1.602e-19;                // electron charge.
const double m_e   = 0.51099906e6;             // electron rest mass [eV/c^2].
const double mu_0  = 4.0*M_PI*1e-7;            // permittivity of free space.
const double eps_0 = 1.0/(sqr(c0)*mu_0);       // permeability of free space.
const double r_e   = q_e/(4.0*M_PI*eps_0*m_e); // classical electron radius.
const double h_bar = 6.58211899e-16;           // reduced Planck constant [eVs].

const double max_ampl = 10.0; // [m].


enum pthicktype {
  thick = 0, thin = 1 };
enum PartsKind {
  drift = 0, Wigl = 1, Mpole = 2, Cavity = 3, marker = 4, undef = 5,
  Insertion = 6, FieldMap = 7, Spreader = 8, Recombiner = 9, Solenoid = 10 };
enum { All = 0, Dip = 1, Quad = 2, Sext = 3, Oct = 4, Dec = 5, Dodec = 6 };
enum { Horizontal = 1, Vertical = 2 };
enum {
  Meth_Linear = 0, Meth_First = 1, Meth_Second = 2, Meth_Fourth = 4,
  Meth_genfun = 5 };

#define fitvectmax   200
#define NameLength   150  // maximum length of identifiers (e.g. file names)
#define DBNameLen    39
#define HOMmax       21
#define SymbolLength 15   // maximum length of element name

typedef char   alfa_[NameLength];

typedef long   iVector2[2];
typedef double Vector2[2];
typedef double Vector3[3];

typedef long   fitvect[fitvectmax];
typedef char   partsName[NameLength];
typedef char   DBNameType[DBNameLen];
typedef double mpolArray[HOMmax+HOMmax+1];


typedef struct statusrec{
  bool tuneflag, chromflag, codflag, mapflag, passflag, overflag, chambre;
  int lossplane; /* lost in: horizontal    1
		             vertical      2
			     longitudinal  3 */
} statusrec;


// Virtual base class for LEGO blocks.
class ElemType {
 private:
 public:
  partsName Name;    // Element name.
  double    L;       // Length [m].
  bool      Reverse; // Reversed element.
  PartsKind Kind;    // Enumeration for magnet types.

  //explicit ElemType(std::string &Name1, double L1, bool Reverse1)
  //  : Name(Name1), L(L1), Reverse(Reverse1) {}

  //virtual void f() = 0;

  // Templates can't be virtual functions. 
  //template<typename T>
  //  virtual void propagate(ss_vect<T> &ps) const; 

  std::ostream& Show(std::ostream &str);
  template<typename T>
    void Marker_Pass(ss_vect<T> &X);
};


class CellType : public ElemType {
 private:
 public:
  int      Fnum;       // Element Family #.
  int      Knum;       // Element Kid #.
  double   S;          // Position in the ring.
  CellType *next_ptr;  // pointer to next cell (for tracking).
  Vector2  dS,         // Transverse displacement.
           dT;         // dT = (cos(dT), sin(dT)).
  Vector2  Nu,         // Phase advances.
           Alpha,      // Alpha functions (redundant).
           Beta,       // beta fonctions (redundant).
           Eta, Etap;  // dispersion and its derivative (redundant).
  psVector BeamPos;    // Last position of the beam this cell.
  Matrix   A,          // Floquet space to phase space transformation.
           sigma;      // sigma matrix (redundant).
  Vector2  maxampl[2]; /* Horizontal and vertical physical apertures:
				maxampl[X_][0] < x < maxampl[X_][1]
				maxampl[Y_][0] < y < maxampl[Y_][1]. */

  void Cell_Init(void);
  std::ostream& Show(std::ostream &str);

  template<typename T>
    T get_p_s(const ss_vect<T> &);

  template<typename T>
    void GtoL(ss_vect<T> &X, const Vector2 &S, const Vector2 &R,
	      const double c0, const double c1, const double s1);

  template<typename T>
    void LtoG(ss_vect<T> &X, const Vector2 &S, const Vector2 &R,
	      const double c0, const double c1, const double s1);

  template<typename T>
    void Pass(ss_vect<T> &x) { this->Pass(x); }

  void Elem_Print(FILE *f, int Fnum);
};


class MarkerType : public CellType {
 private:
 public:

  MarkerType(void);
  void Init(int Fnum);
  template<typename T>
    void Pass(ss_vect<T> &x);

  void Marker_Print(FILE *f, int Fnum);
};


class DriftType : public CellType {
 private:
 public:

  DriftType(void);
  void Init(int Fnum);
  template<typename T>
    void Pass(ss_vect<T> &x);

  void Drift_Print(FILE *f, int Fnum);
};


class MpoleType : public CellType {
 private:
 public:
  int        method;     // Integration Method.
  int        N;          // Number of integration steps.
  // Displacement Errors
  Vector2    dSsys;      // systematic [m].
  Vector2    dSrms;      // rms [m].
  Vector2    dSrnd;      // random number.
  // Roll angle
  double     dTpar;      // design [deg].
  double     dTsys;      // systematic [deg].
  double     dTrms;      // rms [deg].
  double     dTrnd;      // random number.
  // Multipole strengths
  mpolArray  Bpar;       // design.
  mpolArray  Bsys;       // systematic.
  mpolArray  Brms;       // rms.
  mpolArray  Brnd;       // random number.
  mpolArray  B;          // total.
  int        order;      // The highest order in PB.
  int        n_design;   // multipole order (design).
  pthicktype thick;
  // Bending Angles
  double     Tx1;        // horizontal entrance angle [deg].
  double     Tx2;        // horizontal exit angle [deg].
  double     gap;        // total magnet gap [m].
  double     irho;       // 1/rho [1/m].
  double     c0, c1, s1; // corrections for roll error of bend.

  MpoleType(void);
  void Init(int Fnum);
  std::ostream& Show(std::ostream &str);
  template<typename T>
    void Pass(ss_vect<T> &x);

  void SI_init(void);

  void Mpole_Print(FILE *f, int Fnum);
};


class WigglerType : public CellType {
 private:
 public:
  int       method;              // Integration Method.
  int       N;                   // number of integration steps.
  // Displacement Error
  Vector2   dSsys;               // systematic [m].
  Vector2   dSrms;               // rms [m].
  Vector2   dSrnd;               // random number.
  // Roll angle
  double    dTpar;               // design [deg].
  double    dTsys;               // systematic [deg].
  double    dTrms;               // rms [deg].
  double    dTrnd;               // random number.
  double    lambda;              // lambda.
  int       n_harm;              // no of harmonics.
  int       harm[n_harm_max];    // harmonic number.
  double    BoBrhoV[n_harm_max]; // B/Brho vertical.
  double    BoBrhoH[n_harm_max]; // B/Brho horizontal.
  double    kxV[n_harm_max];     // kx.
  double    kxH[n_harm_max];     // kx.
  double    phi[n_harm_max];     // phi.
  mpolArray BW;
  int       order;               // The highest order in PB.

  WigglerType(void);
  void Init(int Fnum);
  void Wiggler_Print(FILE *f, int Fnum);
  template<typename T>
    void Pass(ss_vect<T> &X);
};


class InsertionType : public CellType {
 private:
 public:
  int    method;       // Integration Method.
  int    N;            // number of integration steps.
  char   fname1[100];  // Filename for insertion description: first ordre.
  char   fname2[100];  // Filename for insertion description: second ordre.
  int    nx;           // Horizontal point number.
  int    nz;           // Vertical point number.
  double scaling;      // static scaling factor as in BETA ESRF.
  bool   linear;       // if true linear interpolation else spline.
  bool   firstorder;   // true if first order kick map loaded.
  bool   secondorder;  // true if second order kick map loaded.
  double phi;          // Bend angle.
  double tabx[IDXMAX]; // spacing in H-plane.
  double tabz[IDZMAX]; // spacing in V-plane.
  double thetax[IDZMAX][IDXMAX], thetax1[IDZMAX][IDXMAX]; // 1 for first order.
  double thetaz[IDZMAX][IDXMAX], thetaz1[IDZMAX][IDXMAX];
  bool   long_comp;          // flag for longitudinal comp.
  double B2[IDZMAX][IDXMAX]; // B^2_perp.
  double **tx, **tz, **f2x, **f2z;
  double **tx1, **tz1, **f2x1, **f2z1; // a voir.
  double *tab1, *tab2;       // tab of x and z meshes from Radia code.

  // Displacement Error
  Vector2 dSsys;   // systematic [m]
  Vector2 dSrms;   // rms [m]
  Vector2 dSrnd;   // random number
  // Roll angle
  double  dTpar;   // design [deg]
  double  dTsys;   // systematic [deg]
  double  dTrms;   // rms [deg]
  double  dTrnd;   // random number
  // Strength
  //  double lperiod;  // Length Period [m]
  //  int nperiod;     // Number of periods
  //  double BoBrho;   // B/Brho
  //  double Kx;       // kx
  //  mpolArray BW;
  int order;        // The highest order in PB

  InsertionType(void);
  void Init(int Fnum);
  void Insertion_Print(FILE *f, int Fnum);
  template<typename T>
    void Pass(ss_vect<T> &x);
};


class FieldMapType : public CellType {
 private:
 public:
  int    n_step;                       // number of integration steps.
  int    n[3];                         // no of steps.
  int    cut;                          // cut in z direction.
  double scl, phi, x0, Lr, Lm, Ld, L1;
  double dx[3], *x[3];                 // [dx, dy, dz], [x, y, z].
  double ***BoBrho[3], ***BoBrho2[3];  // [B_x, B_y, B_z].
  double ***AoBrho[2], ***AoBrho2[2];  /* [Ax(x, y, z), Ay(x, y, z)],
					  spline info. */
  FieldMapType(void);
  void get_B(const char *file_name, FieldMapType *FM);
  void Init(int Fnum);
  template<typename T>
    void Pass(ss_vect<T> &x);
};


class SolenoidType : public CellType {
 private:
 public:
  int     N;      // Number of integration steps
  // Displacement Errors
  Vector2 dSsys; // systematic [m]
  Vector2 dSrms; // rms [m]
  Vector2 dSrnd; // random number
  // Roll angle
  double  dTpar;  // design [deg]
  double  dTsys;  // systematic [deg]
  double  dTrms;  // rms [deg]
  double  dTrnd;  // random number
  double  BoBrho; // normalized field strength

  SolenoidType(void);
  void Init(int Fnum);
  template<typename T>
    void Pass(ss_vect<T> &ps);
};


class CavityType : public CellType {
 private:
 public:
  int    N;           // Number of integration steps.
  double volt;        // Vrf [V].
  double freq;        // Vrf [Hz].
  double phi;         // RF phase.
  int    h;           // Harmonic number.
  bool   entry_focus; // Edge focusing at entry.
  bool   exit_focus;  // Edge focusing at exit.

  CavityType(void);
  void Init(int Fnum);
  template<typename T>
    void Pass(ss_vect<T> &X);

};


class SpreaderType : public CellType {
 private:
 public:
  double   E_max[Spreader_max];      // energy levels in increasing order
  CellType *Cell_ptrs[Spreader_max];

  SpreaderType(void);
  void Init(int Fnum);
};


class RecombinerType : public CellType {
 private:
 public:
  double E_min;
  double E_max;

  RecombinerType(void);
  void Init(int Fnum);
};


class ElemFamType : public ElemType {
 private:
 public:
  CellType   CellF;
  int        nKid;             // Kid number.
  int        KidList[nKidMax];
  int        NoDBN;
  DBNameType DBNlist[nKidMax];

  void Fam_Init(void);
};


struct LatticeParam {
  double   dPcommon,       // dp for numerical differentiation.
           dPparticle;     // energy deviation.
  double   delta_RF;       // RF acceptance.
  Vector2  TotalTune;      // transverse tunes.
  double   Omega,
           U0,             // energy lost per turn in keV.
           Alphac;         // alphap.
  Vector2  Chrom;          // chromaticities.
  double   Energy;         // ring energy.
  long     Cell_nLoc,      // number of elements.
           Elem_nFam,      // number of families.
           CODimax;        // Max number of cod search before failing.
  double   CODeps;         // precision for closed orbit finder.
  psVector CODvect;        // closed orbit.
  int      bpm;            // bpm number.
  int      hcorr;          // horizontal corrector number.
  int      vcorr;          // vertical corrector number.
  int      qt;             // vertical corrector number.
  int      gs;             // girder start marker.
  int      ge;             // girder end marker.
  Matrix   OneTurnMat,     // oneturn matrix.
           Ascr,
           Ascrinv,
           Vr,             // real part of the eigenvectors.
           Vi;             // imaginal par of the eigenvectors.

  bool     Cavity_on,      // if true, cavity turned on.
           radiation,      // if true, radiation turned on.
           emittance,
           dip_fringe,     // dipole hard-edge fringe field.
           quad_fringe,    // quadrupole hard-edge fringe field.
           H_exact,        // "small ring" Hamiltonian.
           pathlength,     // absolute path length.
           stable,
           Aperture_on,
           EPU,
           wake_on;

  double   dE,              // energy loss.
           alpha_rad[DOF],  // damping coeffs.
           D_rad[DOF],      // diffusion coeffs (Floquet space).
           J[DOF],          // partition numbers.
           tau[DOF];        // damping times.
  bool     IBS;             // intrabeam scattering.
  double   Qb,              // bunch charge.
           D_IBS[DOF];      // diffusion matrix (Floquet space).
  psVector wr, wi;          // real and imaginary part of eigenvalues.
  double   eps[DOF],        // 3 motion invariants.
           epsp[DOF],       /* transverse and longitudinal projected
			      emittances. */
           alpha_z, beta_z, // longitudinal alpha and beta.
           beta0, gamma0;   // Relativistic factors.
  int      RingType;        // 1 if a ring (0 if transfer line).
};


class LatticeType {
 private:
 public:
  LatticeParam param;
  ElemFamType  ElemFam[Elem_nFamMax];
  CellType     Cell[Cell_nLocMax+1];

  bool Lattice_Read(FILE **fi_, FILE **fo_);
  void Read_Lattice(const char *fic);

  void rdmfile(const char *mfile_dat);
  void prtmfile(const char mfile_dat[]);

  std::ostream& Show_ElemFam(std::ostream &str);
  std::ostream& Show(std::ostream &str);

  // From t2lat.cc.

  long Elem_Index(const std::string &name1);

  // From t2elem.cc.

  int GetnKid(const int Fnum);

  long Elem_GetPos(const int Fnum, const int Knum);


  void getelem(long i, CellType *cellrec);
  void putelem(long i, CellType *cellrec);
 
  double Elem_GetKval(int Fnum, int Knum, int Order);

  void Mpole_SetB(int Fnum, int Knum, int Order);
  double Mpole_GetB(int Fnum, int Knum, int Order);

  void Mpole_DefBpar(int Fnum, int Knum, int Order, double Bpar);
  void Mpole_DefBsys(int Fnum, int Knum, int Order, double Bsys);

  void Mpole_SetdS(int Fnum, int Knum);
  void Mpole_SetdT(int Fnum, int Knum);

  double Mpole_GetdT(int Fnum, int Knum);

  void Mpole_DefdTpar(int Fnum, int Knum, double PdTpar);
  void Mpole_DefdTsys(int Fnum, int Knum, double PdTsys);

  void Wiggler_SetB(int Fnum, int Knum, int Order);

  void Wiggler_SetdS(int Fnum, int Knum);
  void Wiggler_SetdT(int Fnum, int Knum);


  // From t2ring.cc.

  template<typename T>
    void Elem_Pass(const long i, ss_vect<T> &x);

  template<typename T>
    void Cell_Pass(const long i0, const long i1, ss_vect<T> &x, long &lastpos);

  void Cell_Pass(const long i0, const long i1, tps &sigma, long &lastpos);

  void Cell_Init(void);

  bool Cell_getCOD(const long imax, const double eps, const double dP,
		   long &lastpos);
  bool GetCOD(long imax, double eps, double dP, long &lastpos);


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

  void recalc_S();

  double Circumference(void);

  void prt_sigma(void);

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

  /* high level functions for reading lattice file */
  long get_bpm_number(void);
  long get_hcorr_number(void);
  long get_vcorr_number(void);
  long get_qt_number(void);

  /* tracking */
  void GetChromTrac(long Nb, long Nbtour, double emax,
		    double *xix, double *xiz);

  void GetTuneTrac(long Nbtour, double emax, double *nux, double *nuz);

  void Trac(double x, double px, double y, double py, double dp, double ctau,
	    long nmax, long pos, long &lastn, long &lastpos, FILE *outf1);

  /* close orbit */
  // simple precision

  void findcodS(double dP);

  void computeFandJS(double *x, int n, double **fjac, double *fvect);

  // double precision

  void findcod(double dP);

  void computeFandJ(int n, double *x, psVector *fjac, double *fvect);

  /* Vacuum chamber */

  void PrintCh(void);

  void ChamberOff(void);

  // From eigenv.cc.
  void geigen(int n, Matrix &fm, Matrix &Vre, Matrix &Vim,
	      psVector &wr, psVector &wi);


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

  void bend_cal_Fam(const int Fnum);

  void bend_cal(void);

  void set_tune(const char file_name1[], const char file_name2[], const int n);
};


extern double    Fdrift1, Fkick1, Fdrift2, Fkick2, crad, cfluc;
extern tps       sigma_;
extern statusrec status;
extern bool      reverse_elem;
extern double    FindRes_eps;

struct LOC_findres {
  int n;
  double nux, nuy, f;
  int *nx, *ny;
  double eps;
  bool found;
};


// Truncated Power Series Algebra (TPSA).
extern const int nv_tps, nd_tps, iref_tps;
extern int       no_tps, ndpt_tps;
extern double    eps_tps;

extern ss_vect<tps> map;
extern MNF_struct   MNF;

extern double chi_m;


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
#define px_0     0.0
#define py_0     0.0

// 80% sigma coupling.
#define sigma_eps sqrt((25.0/16.0-1.0)/(25.0/16.0+1.0))

#define writetrack true   /*protocol from tracking*/

// getfloq.
#define nfloq    4

// inibump.
#define dnux     0.02
#define dnuy     0.01

// TraceABN.
#define ntrace   4

typedef long   ipeakbuf[npeakmax];
typedef double peakbuf[npeakmax];

struct LOC_getdynap {
  double phi, delta;
  long   nturn;
  bool   floqs, lost;
};


void t2init(void);

template<typename T>
void p_rot(double phi, ss_vect<T> &x);

int Updateorder(CellType &Cell);
