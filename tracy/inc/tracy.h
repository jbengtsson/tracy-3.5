
#ifndef TRACY_H
#define TRACY_H

#define Cell_nLocMax    20000 // maximum number of LEGO blocks (Cell_nLoc).

#ifndef LONG_MAX
# define LONG_MAX       ((long)(((unsigned long) -1) >> 1))
# define LONG_MIN       (~LONG_MAX)
#endif

#define S_SIZE          200  // max size for file name of a lattice file

#define NameLength      150  // maximum length of identifiers (e.g. file names)
#define SymbolLength    15   // maximum length of element name

#define blankname       "               "

#define maxincl         5
#define maxfil          10
#define bigvectmax      4096

// Dynamic aperture (chk_if_lost).
#define px_0            0.0
#define py_0            0.0


typedef char str80[80];


const int max_elem = Cell_nLocMax;

// extern ss_vect<tps> map;
extern MNF_struct   MNF;

extern double       chi_m;


// Macros.

#define degtorad(x) ((x)*M_PI/180.0)
#define radtodeg(x) ((x)*180.0/M_PI)

#define sqr(x)   ((x)*(x))
#define cube(x)  ((x)*(x)*(x))

#define fract(x) ((x)-(int)(x))
#define nint(x) ((x) < 0 ? ((long)(x-0.5)) : ((long)(x+0.5))) 

#define sgn(n) ((n > 0) ? 1 : ((n < 0) ? -1 : 0)) 


// Inline.
inline arma::vec pstoarma(ss_vect<double> &ps)
{
  arma::vec ps_vec = {ps[x_], ps[px_], ps[y_], ps[py_], ps[ct_], ps[delta_]};
  return ps_vec;
}

inline ss_vect<double> armatops(arma::vec &ps_vec)
{
  ss_vect<double> ps(ps_vec[x_], ps_vec[px_], ps_vec[y_], ps_vec[py_],
		     ps_vec[ct_], ps_vec[delta_]);
  return ps;
}


#define HOMmax   21     // [a_n, b_n] <=> [-HOMmax..HOMmax].

#define IDXMAX  200
#define IDZMAX  100

#define DOF     (ss_dim/2)
#define nv_     6

#define debug   false

// Physics constants
const double  c0    = 2.99792458e8;             // speed of light in vacuum
const double  q_e   = 1.602e-19;                // electron charge
const double  m_e   = 0.51099906e6;             // electron rest mass [eV/c^2]
const double  mu_0  = 4.0*M_PI*1e-7;            // permittivity of free space
const double  eps_0 = 1.0/(sqr(c0)*mu_0);       // permeability of free space
const double  r_e   = q_e/(4.0*M_PI*eps_0*m_e); // classical electron radius
const double  h_bar = 6.58211899e-16;           /* reduced Planck constant
						   [eV s] */
typedef char alfa_[NameLength];

typedef long   iVector2[2];
typedef double Vector2[2];
typedef double Vector3[3];

#define fitvectmax      200
typedef long   fitvect[fitvectmax];

extern double Fdrift1, Fkick1, Fdrift2, Fkick2, crad, cfluc;

extern str80 finame,           // input  data file
             foname,           // output data file
             fname;            // temp file name

extern FILE  *fi,              // lattice input  file
             *fo,              // lattice output file
             *psin[],          // program input file
             *psout,           // program output file
             *prr[];           // prr[1] : input, prr[2] : output

extern int P_eof(FILE *f);
extern int P_eoln(FILE *f);

extern void NormEigenVec(arma::mat &Vr, arma::mat &Vi, double *wr, double *wi,
			 arma::mat &t6a);

extern void t2init(void);

extern void prt_gcmat(int bpm, int corr, int plane);

extern void gcmat(int bpm, int corr, int plane);

extern void lsoc(int niter, int bpm, int corr, int plane);


typedef char                partsName[NameLength];
typedef std::vector<double> MpoleArray;

enum PartsKind
  { drift = 0, Wigl = 1, Mpole = 2, Cavity = 3, marker = 4, undef = 5,
    Insertion = 6, FieldMap = 7, Spreader = 8, Recombiner = 9, Solenoid = 10,
    Map = 11 };
enum pthicktype
  { thick = 0, thin = 1 };
enum
  { All = 0, Dip = 1, Quad = 2, Sext = 3, Oct = 4, Dec = 5, Dodec = 6 };
enum
  { Horizontal = 1, Vertical = 2 };
enum
  { Meth_Linear = 0, Meth_First = 1, Meth_Second = 2, Meth_Fourth = 4 };
enum set_mpole { b_n_ = 0, db_n_ = 1, b_nL_ = 2, db_nL_ = 3 };


const int
  n_harm_max   = 10,
  Spreader_max = 10;

const double
  max_ampl = 10.0; // [m]


class ConfigType {
 public:
  bool
    trace,
    reverse_elem,
    stable,
    ErrFlag,
    Cavity_on,                 // if true, cavity turned on
    radiation,                 // if true, radiation turned on
    emittance,
    quad_fringe,               // quadrupole hard-edge fringe field.
    H_exact,                   // "Small Ring" Hamiltonian.
    Cart_Bend,
    dip_edge_fudge,            // Dipole Edge fudge.
    pathlength,                // Absolute Path Length.
    Aperture_on,
    EPU,
    mat_meth,                  // Matrix method.
    IBS,                       // Intrabeam Scattering.
    tuneflag,
    chromflag,
    codflag,
    mapflag,
    passflag,
    overflag,
    chambre;
  long int
    Cell_nLoc,                 // Number of Elements.
    Elem_nFam,                 // Number of Families.
    CODimax;                   // closed Orbit Finder: max number of iterations,
  int
    bpm,                       // BPM Number.
    hcorr,                     // Corrector: Horizontal number,
    vcorr,                     //            Vertical number.
    qt,                        // Vertical corrector number.
    gs,                        // Girder: start marker,
    ge,                        //         end marker.
    RingType,                  // 1 if a ring (0 if transfer line).
    lossplane;                 /* lost in: horizontal    1
		                           vertical      2
			                   longitudinal  3 */
  double
    dPcommon,                  // dp for numerical differentiation.
    dPparticle,                // Energy deviation.
    delta_RF,                  // RF Acceptance.
    TotalTune[DOF],            // Transverse Tunes.
    Omega,                     // Synchrotron Frequency.
    U0,                        // Energy Loss per turn [keV].
    Alphac,                    // Linear Momentum Compaction.
    Energy,                    // Beam Energy.
    dE,                        // Energy Loss.
    alpha_rad[DOF],            // Damping Coeffs.
    D_rad[DOF],                // Diffusion Coeffs (Floquet Space).
    J[DOF],                    // Partition Numbers.
    tau[DOF],                  // Damping Times.
    CODeps,                    //                      precision.
    Qb,                        // Bunch Charge.
    D_IBS[DOF],                // Diffusion Matrix (Floquet Ipace).
    eps[DOF],                  // Eigenemittances.
    epsp[DOF],                 // Trans. & Long. projected Emittances.
    alpha_z,                   // Long. alpha and beta.
    beta_z,
    beta0,                     // Relativistic factors.
    gamma0,
    Chrom[2];                  // Linear Chromaticities.
  ss_vect<double>
    CODvect,                   // Closed Orbit.
    wr,
    wi;                        // Eigenvalues: Real and Imaginary part.
  arma::mat
    OneTurnMat = arma::mat(ss_dim, ss_dim), // Linear Poincare Map.
    Ascr = arma::mat(ss_dim, ss_dim),
    Ascrinv = arma::mat(ss_dim, ss_dim),
    Vr = arma::mat(ss_dim, ss_dim),         // Eigenvectors: Real part, 
    Vi = arma::mat(ss_dim, ss_dim);         //               Imaginary part.
};


// Beam line class.

// LEGO lattice structure.
class CellType {
 public:
  int
    Fnum,                      // Element Family #.
    Knum;                      // Element Kid #.
  double
    S,                         // Position in the ring.
    curly_dH_x,                // Contribution to curly_H_x.
    dI[6],                     // Contribution to I[1..5].
    dS[2],                     // Transverse displacement.
    dT[2],                     // dT = (cos(dT), sin(dT)).
    Nu[2],                     // Phase advances.
    Alpha[2],                  // Alpha functions (redundant).
    Beta[2],                   // beta fonctions (redundant).
    Eta[2],                    // dispersion and its derivative (redundant).
    Etap[2],
    maxampl[2][2];             /* Horizontal and vertical physical apertures:
				  maxampl[X_][0] < x < maxampl[X_][1]
				  maxampl[Y_][0] < y < maxampl[Y_][1].        */
  ss_vect<double>
    BeamPos;                   // Last position of the beam this cell.
  arma::mat
    A = arma::mat(ss_dim, ss_dim),     // Floquet space to phase space
                                       // transformation.
    sigma = arma::mat(ss_dim, ss_dim); // sigma matrix (redundant).
  CellType
    *next_ptr;                         // pointer to next cell (for tracking).
};

// Element base class.
class ElemType : public CellType {
 public:
  partsName
    PName;                     // Element name.
  bool
    Reverse;                   // Reverse element.
  double
    PL;                        // Length[m].
  PartsKind
    Pkind;                     // Enumeration for magnet types.

  virtual ElemType* Elem_Init(const ConfigType &conf, const bool reverse)
  { return NULL; };
  virtual void print(void) {};

  virtual void SetdS(void) {};
  virtual void SetdT(void) {};
  virtual void SetPB(const int n) {};
  virtual double GetdT(void) { return 0e0; };
  virtual double GetPB(const int n) { return 0e0; };

  // C++ templates not supported for virtual functions.
  virtual void Elem_Pass(ConfigType &conf, ss_vect<double> &ps) {};
  virtual void Elem_Pass(ConfigType &conf, ss_vect<tps> &ps) {};

  template<typename T>
  bool CheckAmpl(ConfigType &conf, const ss_vect<T> &x);
  template<typename T>
  void Cell_Pass(ConfigType &conf, ss_vect<T> &ps);
};

// Index for lattice families & elements.
class ElemFamType {
 public:
  ElemType
    *ElemF;
  int
    nKid,                      // No of kids.
    NoDBN;
  std::vector<int>
    KidList;
  std::vector<string>
    DBNlist;                   // For control system.
};

class LatticeType {
 public:
  std::vector<ElemFamType> elemf;
  std::vector<ElemType*>   elems;
  ConfigType               conf;

  void Lat_Init(void);
  void SI_init();

  void prtName(FILE *fp, const int i, const int type, const int method,
	       const int N, const bool reverse);
  void prt_fam(void);
  void prt_elem(void);

  friend long int ElemIndex(const std::string &name);

  void SetdS(const int Fnum, const int Knum);
  void SetdT(const int Fnum, const int Knum);
  void SetPB(const int Fnum, const int Knum, const int Order);
  double GetdT(const int Fnum, const int Knum);
  double GetPB(const int Fnum, const int Knum, const int Order);

  double Elem_GetKval(const int Fnum, const int Knum, const int Order);


  void Mpole_DefPBpar(const int Fnum, const int Knum, const int Order,
		      const double PBpar);
  void Mpole_DefPBsys(const int Fnum, const int Knum, const int Order,
		      const double PBsys);
  void Mpole_DefdTpar(const int Fnum, const int Knum, const double PdTpar);
  void Mpole_DefdTsys(const int Fnum, const int Knum, const double PdTsys);

  // nsls-ii_lib.
  void set_b_n(const set_mpole set_, const int Fnum, const int Knum,
	       const int n, const double b_n);

  bool Lattice_Read(FILE *inf, FILE *outf);
  void prtmfile(const char mfile_dat[]);
  void rdmfile(const char *mfile_dat);

  // t2elem.
  // Obsolete.
  void getelem(long i, ElemType *cellrec);
  void putelem(long i, ElemType *cellrec);

  int GetnKid(const int Fnum);
  long Elem_GetPos(const int Fnum, const int Knum);

  void get_mats(const double delta);

  template<typename T>
  friend T get_p_s(const ConfigType &conf, const ss_vect<T> &ps);
  template<typename T>
  friend void radiate(const ConfigType &conf, ss_vect<T> &ps, const double L,
		      const double h_ref, const T B[]);
  template<typename T>
  void radiate_ID(const ConfigType &conf, ss_vect<T> &ps, const double L,
		  const T &B2_perp);
  friend void emittance(ConfigType &conf, const tps &B2_perp, const tps &ds,
			const tps &p_s0, const ss_vect<tps> &A);

  // t2cell.
  template<typename T>
  void Cell_Pass(const long i0, const long i1, ss_vect<T> &ps, long &lastpos);
  void Cell_Pass(const long i0, const long i1, tps &sigma, long &lastpos);
  bool Cell_getCOD(long imax, double eps, double dp, long &lastpos);
  bool GetCOD(long imax, double eps, double dp, long &lastpos);
  bool getcod(double dp, long &lastpos);

  // t2ring.
  void GDiag(int n_, double C, arma::mat &A, arma::mat &Ainv_, arma::mat &R,
	     arma::mat &M, double &Omega, double &alphac);

  void Cell_Geteta(long i0, long i1, bool ring, double dp);
  void Cell_Twiss(long i0, long i1, ss_vect<tps> &Ascr, bool chroma, bool ring,
		  double dp);
  void Cell_Twiss(const long int i0, const long int i1);
  void TraceABN(long i0, long i1, const Vector2 &alpha, const Vector2 &beta,
		const Vector2 &eta, const Vector2 &etap, const double dp);
  void ttwiss(const Vector2 &alpha, const Vector2 &beta, const Vector2 &eta,
	      const Vector2 &etap, const double dp);
  void Ring_Twiss(bool chroma, double dp);
  void Ring_GetTwiss(bool chroma, double dp);
  void Ring_Getchrom(double dp);

  void Ring_Fittune(Vector2 &nu, double eps, iVector2 &nq, long qf[], long qd[],
		    double dkL, long imax);
  void Ring_Fitchrom(Vector2 &ksi, double eps, iVector2 &ns, long sf[],
		     long sd[], double dkpL, long imax);
  void Ring_FitDisp(long pos, double eta, double eps, long nq, long q[],
		    double dkL, long imax);

  void get_I(double I[], const bool prt);
  template<typename T>
  void Elem_Pass_Lin(ss_vect<T> ps);
  void get_eps_x(double &eps_x, double &sigma_delta, double &U_0,
		 double J[], double tau[], double I[], const bool prt);

  friend void get_dI_eta_5(const int k, ElemType *Elem[]);

  void prt_lat(const int loc1, const int loc2, const char *fname,
	       const bool all);
  void prt_lat(const char *fname, const bool all);
  void prt_lat(const int loc1, const int loc2, const char *fname,
	       const bool all, const int n);
  void prt_lat(const char *fname, const bool all, const int n);
  void prt_chrom_lat(void);
  void prt_cod(const char *file_name, const bool all);
  void prt_beampos(const char *file_name);
  void prt_beamsizes(const int cnt);

  void checkifstable_(struct LOC_Ring_FitDisp *LINK);
  void checkifstable(struct LOC_Ring_Fittune *LINK);

  // Lattice.
  void shiftk(long Elnum, double dk, struct LOC_Ring_Fittune *LINK);
  void shiftk_(long Elnum, double dk, struct LOC_Ring_FitDisp *LINK);
  void shiftkp(long Elnum, double dkp);

  // Vacuum chamber.
  void ChamberOff(void);
  void PrintCh(void);
};

class DriftType : public ElemType {
 public:
  friend DriftType* Drift_Alloc(void);
  ElemType* Elem_Init(const ConfigType &conf, const bool reverse);
  void print(void);

  void SetdS(void) {};
  void SetdT(void) {};
  void SetPB(const int n) {};
  double GetdT(void) { return 0e0; };
  double GetPB(const int n) { return 0e0; };

  template<typename T>
  void Drift_Pass(ConfigType &conf, ss_vect<T> &ps);
  void Elem_Pass(ConfigType &conf, ss_vect<double> &ps)
  { Drift_Pass(conf, ps); };
  void Elem_Pass(ConfigType &conf, ss_vect<tps> &ps)
  { Drift_Pass(conf, ps); };
};

class MpoleType : public ElemType {
 public:
  int
    Pmethod,                   // Integration Method.
    PN,                        // Number of integration steps.
    Porder,                    // The highest order in PB.
    n_design;                  // multipole order (design).
  double
    // Roll angle.
    PdTpar,                    // design [deg].
    PdTsys,                    // systematic [deg].
    PdTrms,                    // rms [deg].
    PdTrnd,                    // random number.
    // Bending Angles.
    PTx1,                      // horizontal entrance angle [deg].
    PTx2,                      // horizontal exit angle [deg].
    Pgap,                      // total magnet gap [m].
    Pirho,                     // 1/rho [1/m].
    Pc0,                       // corrections for roll error of bend.
    Pc1,
    Ps1,
    // Displacement Errors.
    PdSsys[2],                 // systematic [m].
    PdSrms[2],                 // rms [m].
    PdSrnd[2];                 // random number.
  MpoleArray
    // Multipole strengths.
    PBpar,                     // design.
    PBsys,                     // systematic.
    PBrms,                     // rms.
    PBrnd,                     // random number.
    PB;                        // total.
  pthicktype
    Pthick;
  arma::mat
    M_lin = arma::mat(ss_dim, ss_dim); // Linear Map for Element.

  friend MpoleType* Mpole_Alloc(void);
  ElemType* Elem_Init(const ConfigType &conf, const bool reverse);
  void print(void);

  void SetdS(void);
  void SetdT(void);
  void SetPB(const int n);
  double GetdT(void);
  double GetPB(const int n);

  template<typename T>
  void Mpole_Pass(ConfigType &conf, ss_vect<T> &ps);
  void Elem_Pass(ConfigType &conf, ss_vect<double> &ps)
  { Mpole_Pass(conf, ps); };
  void Elem_Pass(ConfigType &conf, ss_vect<tps> &ps)
  { Mpole_Pass(conf, ps); };
};

class CavityType : public ElemType {
 public:
  bool
    entry_focus,               // Edge focusing at entry.
    exit_focus;                // Edge focusing at exit.
  int
    PN,                        // Number of integration steps.
    Ph;                        // Harmonic number.
  double
    Pvolt,                     // Vrf [V].
    Pfreq,                     // Vrf [Hz].
    phi;                       // RF phase.

  friend CavityType* Cavity_Alloc(void);
  ElemType* Elem_Init(const ConfigType &conf, const bool reverse);
  void print(void);

  void SetdS(void) {};
  void SetdT(void) {};
  void SetPB(const int n) {};
  double GetdT(void) { return 0e0; };
  double GetPB(const int n) { return 0e0; };

  template<typename T>
  void Cavity_Pass(ConfigType &conf, ss_vect<T> &ps);
  void Elem_Pass(ConfigType &conf, ss_vect<double> &ps)
  { Cavity_Pass(conf, ps); };
  void Elem_Pass(ConfigType &conf, ss_vect<tps> &ps)
  { Cavity_Pass(conf, ps); };
};

class MarkerType : public ElemType {
 public:
  friend MarkerType* Marker_Alloc(void);
  ElemType* Elem_Init(const ConfigType &conf, const bool reverse);
  void print(void);

  void SetdS(void) {};
  void SetdT(void) {};
  void SetPB(const int n) {};
  double GetdT(void) { return 0e0; };
  double GetPB(const int n) { return 0e0; };

  template<typename T>
  void Marker_Pass(ConfigType &conf, ss_vect<T> &ps);
  void Elem_Pass(ConfigType &conf, ss_vect<double> &ps)
  { Marker_Pass(conf, ps); };
  void Elem_Pass(ConfigType &conf, ss_vect<tps> &ps)
  { Marker_Pass(conf, ps); };
};

class WigglerType : public ElemType {
 public:
  int
    Pmethod,                   // Integration Method.
    PN,                        // number of integration steps.
    n_harm,                    // No of harmonics.
    harm[n_harm_max],          // Harmonic number.
    Porder;                    // The highest order in PB.
    // Roll angle
  double
    PdTpar,                    // Design [deg].
    PdTsys,                    // Systematic [deg].
    PdTrms,                    // RMS [deg].
    PdTrnd,                    // Random number.
    Lambda,                    // lambda.
    // Displacement Error.
    PdSsys[2],                 // Systematic [m].
    PdSrms[2],                 // RMS [m].
    PdSrnd[2],                 // Random number.
    BoBrhoV[n_harm_max],       // B/Brho vertical.
    BoBrhoH[n_harm_max],       // B/Brho horizontal.
    kxV[n_harm_max],           // kx.
    kxH[n_harm_max],           // kx.
    phi[n_harm_max];           // phi.
  MpoleArray
    PBW;

  friend WigglerType* Wiggler_Alloc(void);
  ElemType* Elem_Init(const ConfigType &conf, const bool reverse);
  void print(void);

  void SetdS(void);
  void SetdT(void);
  void SetPB(const int n) {};
  double GetdT(void);
  double GetPB(const int n);

  template<typename T>
  void Wiggler_Pass(ConfigType &conf, ss_vect<T> &ps);
  void Elem_Pass(ConfigType &conf, ss_vect<double> &ps)
  { Wiggler_Pass(conf, ps); };
  void Elem_Pass(ConfigType &conf, ss_vect<tps> &ps)
  { Wiggler_Pass(conf, ps); };
};

class InsertionType : public ElemType {
 public:
  char
    fname1[100],               // Filename for insertion description: 1st order.
    fname2[100];               // Filename for insertion description: 2nd order.
  bool
    linear,                    // if true linear interpolation else spline.
    firstorder,                // true if first order kick map loaded.
    secondorder,               // true if second order kick map loaded.
    long_comp;                 // flag for longitudinal comp.
  int
    Pmethod,                   // Integration Method.
    PN,                        // number of integration steps.
    nx,                        // Horizontal point number.
    nz,                        // Vertical point number.
    Porder;                    // The highest order in PB.
  double
    scaling,                   // static scaling factor as in BETA ESRF.
    phi,                       // Bend angle.
    tabx[IDXMAX],              // spacing in H-plane.
    tabz[IDZMAX],              // spacing in V-plane.
    thetax[IDZMAX][IDXMAX],
    thetax1[IDZMAX][IDXMAX],   // 1 for first order.
    thetaz[IDZMAX][IDXMAX],
    thetaz1[IDZMAX][IDXMAX],
    B2[IDZMAX][IDXMAX],        // B^2_perp.
    **tx,
    **tz,
    **f2x,
    **f2z,
    **tx1,
    **tz1,
    **f2x1,
    **f2z1,                    // a voir.
    *tab1,
    *tab2,                     // tab of x and z meshes from Radia code.
    // Roll angle.
    PdTpar,                    // design [deg].
    PdTsys,                    // systematic [deg].
    PdTrms,                    // rms [deg].
    PdTrnd,                    // random number.
//  Strength
//  double Plperiod;           // Length Period [m].
//  int Pnperiod;              // Number of periods.
//  double PBoBrho;            // B/Brho.
//  double PKx;                // kx.
//  mpolArray PBW;
  // Displacement Error.
    PdSsys[2],                 // systematic [m].
    PdSrms[2],                 // rms [m].
    PdSrnd[2];                 // random number.

  friend InsertionType* Insertion_Alloc(void);
  ElemType* Elem_Init(const ConfigType &conf, const bool reverse);
  void print(void);

  void SetdS(void) {};
  void SetdT(void) {};
  void SetPB(const int n) {};
  double GetdT(void) { return 0e0; };
  double GetPB(const int n) { return 0e0; };

  template<typename T>
  void Insertion_Pass(ConfigType &conf, ss_vect<T> &ps);
  void Elem_Pass(ConfigType &conf, ss_vect<double> &ps)
  { Insertion_Pass(conf, ps); };
  void Elem_Pass(ConfigType &conf, ss_vect<tps> &ps)
  { Insertion_Pass(conf, ps); };
};

class FieldMapType : public ElemType {
 public:
  int
    n_step,                    // number of integration steps.
    n[3],                      // no of steps.
    cut;                       // cut in z direction.
  double
    scl,
    phi,
    x0,
    Lr,
    Lm,
    Ld,
    L1,
    dx[3],
    *x[3],                     // [dx, dy, dz], [x, y, z].
    ***BoBrho[3],
    ***BoBrho2[3],             // [B_x, B_y, B_z].
    ***AoBrho[2],
    ***AoBrho2[2];             // [Ax(x, y, z), Ay(x, y, z)], spline info.

  friend FieldMapType* FieldMap_Alloc(void);
  ElemType* Elem_Init(const ConfigType &conf, const bool reverse);
  void print(void);

  void SetdS(void) {};
  void SetdT(void) {};
  void SetPB(const int n) {};
  double GetdT(void) { return 0e0; };
  double GetPB(const int n) { return 0e0; };

  template<typename T>
  void FieldMap_Pass(ConfigType &conf, ss_vect<T> &ps);
  void Elem_Pass(ConfigType &conf, ss_vect<double> &ps)
  { FieldMap_Pass(conf, ps); };
  void Elem_Pass(ConfigType &conf, ss_vect<tps> &ps)
  { FieldMap_Pass(conf, ps); };
};

class SpreaderType : public ElemType {
 public:
  double
    E_max[Spreader_max];       // energy levels in increasing order.
  CellType
    *Cell_ptrs[Spreader_max];

  friend SpreaderType* Spreader_Alloc(void);
  ElemType* Elem_Init(const ConfigType &conf, const bool reverse);
  void print(void);

  void SetdS(void) {};
  void SetdT(void) {};
  void SetPB(const int n) {};
  double GetdT(void) { return 0e0; };
  double GetPB(const int n) { return 0e0; };

  template<typename T>
  void Spreader_Pass(ConfigType &conf, ss_vect<T> &ps);
  void Elem_Pass(ConfigType &conf, ss_vect<double> &ps)
  { Spreader_Pass(conf, ps); };
  void Elem_Pass(ConfigType &conf, ss_vect<tps> &ps)
  { Spreader_Pass(conf, ps); };
};

class RecombinerType : public ElemType {
 public:
  double
    E_min,
    E_max;

  friend RecombinerType* Recombiner_Alloc(void);
  ElemType* Elem_Init(const ConfigType &conf, const bool reverse);
  void print(void);

  void SetdS(void) {};
  void SetdT(void) {};
  void SetPB(const int n) {};
  double GetdT(void) { return 0e0; };
  double GetPB(const int n) { return 0e0; };

  template<typename T>
  void Recombiner_Pass(ConfigType &conf, ss_vect<T> &ps);
  void Elem_Pass(ConfigType &conf, ss_vect<double> &ps)
  { Recombiner_Pass(conf, ps); };
  void Elem_Pass(ConfigType &conf, ss_vect<tps> &ps)
  { Recombiner_Pass(conf, ps); };
};

class SolenoidType : public ElemType {
 public:
  int
    N;                         // Number of integration steps.
  double
    BoBrho,                    // normalized field strength.
    // Roll angle.
    dTpar,                     // design [deg].
    dTsys,                     // systematic [deg].
    dTrms,                     // rms [deg].
    dTrnd,                     // random number.
    // Displacement Errors.
    PdSsys[2],                 // systematic [m].
    PdSrms[2],                 // rms [m].
    PdSrnd[2];                 // random number.

  friend SolenoidType* Solenoid_Alloc(void);
  ElemType* Elem_Init(const ConfigType &conf, const bool reverse);
  void print(void);

  void SetdS(void) {};
  void SetdT(void) {};
  void SetPB(const int n) {};
  double GetdT(void) { return 0e0; };
  double GetPB(const int n) { return 0e0; };

  template<typename T>
  void Solenoid_Pass(ConfigType &conf, ss_vect<T> &ps);
  void Elem_Pass(ConfigType &conf, ss_vect<double> &ps)
  { Solenoid_Pass(conf, ps); };
  void Elem_Pass(ConfigType &conf, ss_vect<tps> &ps)
  { Solenoid_Pass(conf, ps); };
};

class MapType : public ElemType {
 public:
  double
    dnu[2],
    alpha[2],
    beta[2],
    eta_x,
    etap_x;
  ss_vect<tps>
    M;

  friend MapType* Map_Alloc(void);
  ElemType* Elem_Init(const ConfigType &conf, const bool reverse);
  void print(void);

  void SetdS(void) {};
  void SetdT(void) {};
  void SetPB(const int n) {};
  double GetdT(void) { return 0e0; };
  double GetPB(const int n) { return 0e0; };

  template<typename T>
  void Map_Pass(ConfigType &conf, ss_vect<T> &ps);
  void Elem_Pass(ConfigType &conf, ss_vect<double> &ps)
  { Map_Pass(conf, ps); };
  void Elem_Pass(ConfigType &conf, ss_vect<tps> &ps)
  { Map_Pass(conf, ps); };
};


void file_rd(std::ifstream &inf, const string &file_name);

void file_wr(std::ofstream &outf, const string &file_name);

void file_rd(std::ifstream &inf, const char file_name[]);

void file_wr(std::ofstream &outf, const char file_name[]);

FILE* file_read(const char file_name[]);

FILE* file_write(const char file_name[]);


void t2init(void);

void exit_(int exit_code);

double xabs(long n, ss_vect<double> &x);

#endif
