
#ifndef TRACY_H
#define TRACY_H

#define Cell_nLocMax    20000 // maximum number of LEGO blocks (Cell_nLoc).

#ifndef LONG_MAX
# define LONG_MAX       ((long)(((unsigned long) -1) >> 1))
# define LONG_MIN       (~LONG_MAX)
#endif

#define maxincl         5
#define maxfil          10

// Dynamic aperture (chk_if_lost).
#define px_0            0e0
#define py_0            0e0


// Macros.

#define degtorad(x) ((x)*M_PI/180.0)
#define radtodeg(x) ((x)*180.0/M_PI)

#define sqr(x)      ((x)*(x))
#define cube(x)     ((x)*(x)*(x))

#define fract(x)    ((x)-(int)(x))
#define nint(x)     ((x) < 0 ? ((long)(x-0.5)) : ((long)(x+0.5))) 

#define sgn(n)      ((n > 0) ? 1 : ((n < 0) ? -1 : 0)) 


// Inline functions.

inline arma::vec pstovec(const ss_vect<double> &vec)
{ return {vec[x_], vec[px_], vec[y_], vec[py_], vec[delta_], vec[ct_], 1e0}; }

inline ss_vect<double> vectops(const arma::vec vec)
{ return {vec[x_], vec[px_], vec[y_], vec[py_], vec[delta_], vec[ct_]}; }

inline std::vector<double> vectostlvec(const arma::vec &vec)
{ return {vec[x_], vec[px_], vec[y_], vec[py_], vec[delta_], vec[ct_], 1e0}; }

inline std::vector<double> stlvectovec(const arma::vec &vec)
{ return {vec[x_], vec[px_], vec[y_], vec[py_], vec[delta_], vec[ct_], 1e0}; }

inline std::vector<double> pstostlvec(const ss_vect<double> &vec)
{ return {vec[x_], vec[px_], vec[y_], vec[py_], vec[delta_], vec[ct_], 1e0}; }

inline ss_vect<double> stlvectops(const std::vector<double> &vec)
{ return {vec[x_], vec[px_], vec[y_], vec[py_], vec[delta_], vec[ct_]}; }


#define HOMmax   21     // [a_n, b_n] <=> [-HOMmax..HOMmax].

#define IDXMAX  200
#define IDZMAX  100

// Physics constants
const double  c0    = 2.99792458e8;             // speed of light in vacuum
const double  q_e   = 1.602e-19;                // electron charge
const double  m_e   = 0.51099906e6;             // electron rest mass [eV/c^2]
const double  mu_0  = 4.0*M_PI*1e-7;            // permittivity of free space
const double  eps_0 = 1.0/(sqr(c0)*mu_0);       // permeability of free space
const double  r_e   = q_e/(4.0*M_PI*eps_0*m_e); // classical electron radius
const double  h_bar = 6.58211899e-16;           /* reduced Planck constant
						   [eV s] */
#define fitvectmax 200
typedef long   fitvect[fitvectmax];

extern double Fdrift1, Fkick1, Fdrift2, Fkick2, crad, cfluc;

extern int P_eof(FILE *f);
extern int P_eoln(FILE *f);

extern void t2init(void);

extern void prt_gcmat(int bpm, int corr, int plane);

extern void gcmat(int bpm, int corr, int plane);

extern void lsoc(int niter, int bpm, int corr, int plane);


typedef std::vector<double> MpoleArray;


enum PartsKind
  { drift      = 0,
    Wigl       = 1,
    Mpole      = 2,
    Cavity     = 3,
    marker     = 4,
    undef      = 5,
    Insertion  = 6,
    FieldMap   = 7,
    Spreader   = 8,
    Recombiner = 9,
    Solenoid   = 10,
    Map        = 11 };

enum pthicktype
  { thick = 0,
    thin  = 1 };

enum MpoleKind
  { All   = 0,
    Dip   = 1,
    Quad  = 2,
    Sext  = 3,
    Oct   = 4,
    Dec   = 5,
    Dodec = 6 };

enum PlaneKind
  { Horizontal = 1,
    Vertical   = 2 };

enum IntMethKind
  { Meth_Linear = 0,
    Meth_First  = 1,
    Meth_Second = 2,
    Meth_Fourth = 4 };

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
    Cavity_on,                    // if true, cavity turned on
    radiation,                    // if true, radiation turned on
    emittance,
    quad_fringe,                  // quadrupole hard-edge fringe field.
    H_exact,                      // "Small Ring" Hamiltonian.
    Cart_Bend,
    dip_edge_fudge,               // Dipole Edge fudge.
    pathlength,                   // Absolute Path Length.
    Aperture_on,
    EPU,
    mat_meth,                     // Matrix method.
    IBS,                          // Intrabeam Scattering.
    tuneflag,
    chromflag,
    codflag,
    mapflag,
    passflag,
    overflag,
    chambre;
  long int
    Cell_nLoc,                    // Number of Elements.
    Elem_nFam,                    // Number of Families.
    CODimax;                      /* closed Orbit Finder: max number of
                                     iterations. */
  int
    bpm,                          // BPM Number.
    hcorr,                        // Corrector: Horizontal number,
    vcorr,                        //            Vertical number.
    qt,                           // Vertical corrector number.
    gs,                           // Girder: start marker,
    ge,                           //         end marker.
    RingType,                     // 1 if a ring (0 if transfer line).
    lossplane;                    /* lost in: horizontal    1
		                              vertical      2
			                      longitudinal  3 */
  double
    dPcommon,                     // dp for numerical differentiation.
    dPparticle,                   // Energy deviation.
    delta_RF,                     // RF Acceptance.
    Omega,                        // Synchrotron Frequency.
    U0,                           // Energy Loss per turn [keV].
    Alphac,                       // Linear Momentum Compaction.
    Energy,                       // Beam Energy.
    dE,                           // Energy Loss.
    CODeps,                       // Closed Orbit precision.
    Qb,                           // Bunch Charge.
    alpha_z,                      // Long. alpha and beta.
    beta_z,
    beta0,                        // Relativistic factors.
    gamma0;
  std::vector<double>
    TotalTune{0e0, 0e0, 0e0},     // Transverse tunes.
    Chrom{0e0, 0e0},              // Linear chromaticities.
    alpha_rad{0e0, 0e0, 0e0},     // Damping coeffs.
    D_rad{0e0, 0e0, 0e0},         // Diffusion coeffs (Floquet space).
    J{0e0, 0e0, 0e0},             // Partition numbers.
    tau{0e0, 0e0, 0e0},           // Damping times.
    D_IBS{0e0, 0e0, 0e0},         // Diffusion matrix (Floquet space).
    eps{0e0, 0e0, 0e0},           // Eigenemittances.
    epsp{0e0, 0e0, 0e0},          // Trans. & Long. projected emittances.
    CODvect                       // Closed orbit.
      {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
    wr                            // Eigenvalues Re & Im part.
      {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
    wi
      {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0};
  std::vector< std::vector<double> >
    OneTurnMat
      {{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0}},
    Ascr
      {{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0}},
    Ascrinv
      {{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0}},
    Vr                            // Eigenvectors: Real part,
      {{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0}},
  
    Vi                            //               Imaginary part.
      {{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0}};
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
    curly_dH_x;                // Contribution to curly_H_x.
  std::vector<double>
    dI                         // Contribution to I[1..5].
    {0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
    dS{0e0, 0e0},              // Transverse displacement.
    dT{0e0, 0e0},              // dT = (cos(dT), sin(dT)).
    Eta{0e0, 0e0},             // Eta & eta' (redundant).
    Etap{0e0, 0e0},
    Alpha{0e0, 0e0},           // Twiss parameters (redundant).
    Beta{0e0, 0e0},
    Nu{0e0, 0e0},
    BeamPos                    // Last position of the beam this cell.
      {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0};
  std::vector< std::vector<double> >
    maxampl                    // Hor & ver physical aperture.
      {{0e0, 0e0},
       {0e0, 0e0}},
    A                          // Floquet space to phase space transformation.
      {{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0}},
    sigma                      // sigma matrix (redundant).
      {{0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0}};
    CellType
      *next_ptr;               // pointer to next cell (for tracking).
};

// Element virtual base class.
class ElemType : public CellType {
 public:
  std::string
    Name;                      // Element name.
  bool
    Reverse;                   // Reverse element.
  double
    PL;                        // Length[m].
  PartsKind
    Pkind;                     // Enumeration for magnet types.

  void prt_elem(void);

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
  std::vector<std::string>
    DBNlist;                   // For control system.
};

class LatticeType {
 public:
  std::vector<ElemFamType> elemf;
  std::vector<ElemType*>   elems;
  ConfigType               conf;

  void Lat_Init(void);
  void SI_init();

  void prt_fams(void);
  void prt_elems(void);

  int GetnKid(const int Fnum);
  long int ElemIndex(const std::string &name);
  long int Elem_GetPos(const int Fnum, const int Knum);

  void get_transp_mats(const double delta);

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

  bool Lat_Read(const std::string &filnam);
  bool Lat_Read_py(const std::string &filnam);
  void prtmfile(const std::string &mfile_dat);
  void rdmfile(const std::string &mfile_dat);

  // t2elem.

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
  void TraceABN(long i0, long i1, const std::vector<double> &alpha,
		const std::vector<double> &beta, const std::vector<double> &eta,
		const std::vector<double> &etap, const double dp);
  void ttwiss(const std::vector<double> &alpha, const std::vector<double> &beta,
	      const std::vector<double> &eta, const std::vector<double> &etap,
	      const double dp);
  void Ring_Twiss(bool chroma, double dp);
  void Ring_GetTwiss(bool chroma, double dp);
  void Ring_Getchrom(double dp);

  void Ring_Fittune(std::vector<double> &nu, double eps, std::vector<int> &nq,
		    long qf[], long qd[], double dkL, long imax);
  void Ring_Fitchrom(std::vector<double> &ksi, double eps, std::vector<int> &ns,
		     long sf[], long sd[], double dkpL, long imax);
  void Ring_FitDisp(long pos, double eta, double eps, long nq, long q[],
		    double dkL, long imax);

  void get_I(std::vector<double> &I, const bool prt);
  template<typename T>
  void Elem_Pass_Lin(ss_vect<T> ps);
  void get_eps_x(double &eps_x, double &sigma_delta, double &U_0,
		 std::vector<double> &J, std::vector<double> &tau,
		 std::vector<double> &I, const bool prt);

  friend void get_dI_eta_5(const int k, ElemType *Elem[]);

  void prt_lat1(const int loc1, const int loc2, const std::string &fname,
		const bool all);
  void prt_lat2(const std::string &fname, const bool all);
  void prt_lat3(const int loc1, const int loc2, const std::string &fname,
		const bool all, const int n);
  void prt_lat4(const std::string &fname, const bool all, const int n);
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

  void print(void);

  void GetEmittance(const int Fnum, const bool prt);
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
    PdTpar,                    // Roll angle [deg]: design
    PdTsys,                    //                    systematic
    PdTrms,                    //                    rms
    PdTrnd,                    //                    random number.
    PTx1,                      // Bend angle [deg]: hor. entrance angle
    PTx2,                      //                   hor. exit angle.
    Pgap,                      // Total magnet gap [m].
    Pirho,                     // 1/rho [1/m].
    Pc0,                       // Corrections for roll error of bend.
    Pc1,
    Ps1;
  std::vector<double>
    PdSsys{0e0, 0e0},          // Displacement errors [m]: systematic
    PdSrms{0e0, 0e0},          //                          rms
    PdSrnd{0e0, 0e0};          //                          random number.
  MpoleArray
    PBpar,                     // Multipole strengths: design
    PBsys,                     //                      systematic
    PBrms,                     //                      rms
    PBrnd,                     //                      random number
    PB;                        //                      total.
  pthicktype
    Pthick;
  std::vector< std::vector<double> >
    M_elem                     // Transport matrix & orbit.
      {{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0}};


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
  double
    PdTpar,                    // Roll angle [deg]: design
    PdTsys,                    //                   systematic
    PdTrms,                    //                   RMS
    PdTrnd,                    //                   random number.
    Lambda,                    // lambda.
    BoBrhoV[n_harm_max],       // B/Brho: ver.
    BoBrhoH[n_harm_max],       //         hor.
    kxV[n_harm_max],           // kx.
    kxH[n_harm_max],           // kx.
    phi[n_harm_max];           // phi.
  std::vector<double>
    PdSsys{0e0, 0e0},          // Displacement error [m]: systematic
    PdSrms{0e0, 0e0},          //                         RMS
    PdSrnd{0e0, 0e0};          //                         random number.
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
    PdTrnd;                    // random number.
//  Strength
//  double Plperiod;           // Length Period [m].
//  int Pnperiod;              // Number of periods.
//  double PBoBrho;            // B/Brho.
//  double PKx;                // kx.
//  mpolArray PBW;
  std::vector<double>
    PdSsys{0e0, 0e0},          // Displacement error [m]: systematic
    PdSrms{0e0, 0e0},          //                         rms
    PdSrnd{0e0, 0e0};          //                         random number.

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
    dTrnd;                     // random number.
  std::vector<double>
    PdSsys{0e0, 0e0},          // Displacement errors [m]: systematic
    PdSrms{0e0, 0e0},          //                          rms
    PdSrnd{0e0, 0e0};          //                          random number.

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
    eta_x,
    etap_x;
  std::vector<double>
    dnu{0e0, 0e0},
    alpha{0e0, 0e0},
    beta{0e0, 0e0};
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


void file_rd(std::ifstream &inf, const std::string &file_name);

void file_wr(std::ofstream &outf, const std::string &file_name);

void file_rd(std::ifstream &inf, const char file_name[]);

void file_wr(std::ofstream &outf, const char file_name[]);

FILE* file_read(const char file_name[]);

FILE* file_write(const char file_name[]);


void t2init(void);

void exit_(int exit_code);

double xabs(long n, ss_vect<double> &x);

#endif
