
#ifndef TRACY_GLOBAL_H
#define TRACY_GLOBAL_H

#define HOMmax          21     // [a_n, b_n] <=> [-HOMmax..HOMmax].
#define nKidMax       1000     // maximum number of kids.

#define IDXMAX         200
#define IDZMAX         100

#define DOF     (ss_dim/2)
#define nv_              6

#define debug       false

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

const int
  n_harm_max   = 10,
  Spreader_max = 10;


typedef class globvalrec {
 public:
  bool
    Cavity_on,                 // if true, cavity turned on
    radiation,                 // if true, radiation turned on
    emittance,
    quad_fringe,               // quadrupole hard-edge fringe field.
    H_exact,                   // "Small Ring" Hamiltonian.
    Cart_Bend,
    dip_edge_fudge,            // Dipole Edge fudge.
    pathlength,                // Absolute Path Length.
    stable,
    Aperture_on,
    EPU,
    mat_meth,                  // Matrix method.
    IBS;                       // Intrabeam Scattering.
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
    RingType;                  // 1 if a ring (0 if transfer line).
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
  psVector
    CODvect,                   // Closed Orbit.
    wr,
    wi;                        // Eigenvalues: Real and Imaginary part.
  Matrix
    OneTurnMat,                // Linear Poincare Map.
    Ascr,
    Ascrinv,
    Vr,                        // Eigenvectors: Real part, 
    Vi;                        //               Imaginary part.
} globvalrec;



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
  psVector
    BeamPos;                   // Last position of the beam this cell.
  Matrix
    A,                         // Floquet space to phase space transformation.
    sigma;                     // sigma matrix (redundant).
  CellType
    *next_ptr;                 // pointer to next cell (for tracking).
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

  // Wrapper functions; becuase C++ does not support templates for virtual
  // functions.
  virtual void Elem_Pass(ss_vect<double> &ps) {};
  virtual void Elem_Pass(ss_vect<tps> &ps) {};

  virtual void print(void) {};

  template<typename T>
  bool CheckAmpl(const ss_vect<T> &x);
  template<typename T>
  void Cell_Pass(ss_vect<T> &ps);
};

// Index for lattice elements.
class ElemFamType {
 public:
  ElemType
    *ElemF;
  int
    nKid,                      // No of kids.
    KidList[nKidMax],
    NoDBN;
  std::vector<string>
    DBNlist;                   // For control system.
};

class LatticeType {
 public:
  std::vector<ElemFamType> elemf;
  std::vector<ElemType*>   elems;

  void Drift_Init(const int Fnum);
  void Mpole_Init(const int Fnum);
  void Cavity_Init(const int Fnum);
  void Marker_Init(const int Fnum);
  void Wiggler_Init(const int Fnum);
  void Insertion_Init(const int Fnum);
  void FieldMap_Init(const int Fnum);
  void Spreader_Init(const int Fnum);
  void Recombiner_Init(const int Fnum);
  void Solenoid_Init(const int Fnum);
  void Map_Init(const int Fnum);
  void Lat_Init(void);

  bool Lattice_Read(FILE *inf, FILE *outf);
  friend long int ElemIndex(const std::string &name1);
  void prt_fam(void);
  void prt_elem(void);

  // t2elem.
  void getelem(long i, ElemType *cellrec);
  void putelem(long i, ElemType *cellrec);
  int GetnKid(const int Fnum1);
  long Elem_GetPos(const int Fnum1, const int Knum1);
  double Elem_GetKval(int Fnum1, int Knum1, int Order);
  void get_lin_maps(const double delta);
  void Mpole_SetPB(int Fnum1, int Knum1, int Order);
  double Mpole_GetPB(int Fnum1, int Knum1, int Order);
  void Mpole_DefPBpar(int Fnum1, int Knum1, int Order, double PBpar);
  void Mpole_DefPBsys(int Fnum1, int Knum1, int Order, double PBsys);
  void Mpole_SetdS(int Fnum1, int Knum1);
  void Mpole_SetdT(int Fnum1, int Knum1);
  double Mpole_GetdT(int Fnum1, int Knum1);
  void Mpole_DefdTpar(int Fnum1, int Knum1, double PdTpar);
  void Mpole_DefdTsys(int Fnum1, int Knum1, double PdTsys);
  void Wiggler_SetPB(int Fnum1, int Knum1, int Order);
  void Wiggler_SetdS(int Fnum1, int Knum1);
  void Wiggler_SetdT(int Fnum1, int Knum1);

  // t2cell.
  template<typename T>
  void Cell_Pass(const long i0, const long i1, ss_vect<T> &ps, long &lastpos);
  void Cell_Pass(const long i0, const long i1, tps &sigma, long &lastpos);
  bool Cell_getCOD(long imax, double eps, double dp, long &lastpos);
  bool GetCOD(long imax, double eps, double dp, long &lastpos);
  bool getcod(double dp, long &lastpos);

  // t2ring.
  void shiftk(long Elnum, double dk, struct LOC_Ring_Fittune *LINK);
  void shiftkp(long Elnum, double dkp);
  void shiftk_(long Elnum, double dk, struct LOC_Ring_FitDisp *LINK);

  void Cell_Geteta(long i0, long i1, bool ring, double dp);
  void Cell_Twiss(long i0, long i1, ss_vect<tps> &Ascr, bool chroma, bool ring,
		  double dp);
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
  void get_eps_x(double &eps_x, double &sigma_delta, double &U_0,
		 double J[], double tau[], double I[], const bool prt);

  friend void get_dI_eta_5(const int k, ElemType *Elem[]);
};

class DriftType : public ElemType {
 public:
  friend DriftType* Drift_Alloc(void);

  template<typename T>
  void Drift_Pass(ss_vect<T> &ps);

  void Elem_Pass(ss_vect<double> &ps) { Drift_Pass(ps); };
  void Elem_Pass(ss_vect<tps> &ps) { Drift_Pass(ps); };

  void print(void);
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
  ss_vect<tps>
    M_lin;                     // Linear Map for Element.

  friend MpoleType* Mpole_Alloc(void);

  template<typename T>
  void Mpole_Pass(ss_vect<T> &ps);

  void Elem_Pass(ss_vect<double> &ps) { Mpole_Pass(ps); };
  void Elem_Pass(ss_vect<tps> &ps) { Mpole_Pass(ps); };

  void print(void);
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

  template<typename T>
  void Cavity_Pass(ss_vect<T> &ps);

  void Elem_Pass(ss_vect<double> &ps) { Cavity_Pass(ps); };
  void Elem_Pass(ss_vect<tps> &ps) { Cavity_Pass(ps); };

  void print(void);
};

class MarkerType : public ElemType {
 public:
  friend MarkerType* Marker_Alloc(void);

  template<typename T>
  void Marker_Pass(ss_vect<T> &ps);

  void Elem_Pass(ss_vect<double> &ps) { Marker_Pass(ps); };
  void Elem_Pass(ss_vect<tps> &ps) { Marker_Pass(ps); };

  void print(void);
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

  template<typename T>
  void Wiggler_Pass(ss_vect<T> &ps);

  void Elem_Pass(ss_vect<double> &ps) { Wiggler_Pass(ps); };
  void Elem_Pass(ss_vect<tps> &ps) { Wiggler_Pass(ps); };

  void print(void);
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

  template<typename T>
  void Insertion_Pass(ss_vect<T> &ps);

  void Elem_Pass(ss_vect<double> &ps) { Insertion_Pass(ps); };
  void Elem_Pass(ss_vect<tps> &ps) { Insertion_Pass(ps); };

  void print(void);
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

  template<typename T>
  void FieldMap_Pass(ss_vect<T> &ps);

  void Elem_Pass(ss_vect<double> &ps) { FieldMap_Pass(ps); };
  void Elem_Pass(ss_vect<tps> &ps) { FieldMap_Pass(ps); };

  void print(void);
};

class SpreaderType : public ElemType {
 public:
  double
    E_max[Spreader_max];       // energy levels in increasing order.
  CellType
    *Cell_ptrs[Spreader_max];

  friend SpreaderType* Spreader_Alloc(void);

  template<typename T>
  void Spreader_Pass(ss_vect<T> &ps);

  void Elem_Pass(ss_vect<double> &ps) { Spreader_Pass(ps); };
  void Elem_Pass(ss_vect<tps> &ps) { Spreader_Pass(ps); };

  void print(void);
};

class RecombinerType : public ElemType {
 public:
  double
    E_min,
    E_max;

  friend RecombinerType* Recombiner_Alloc(void);

  template<typename T>
  void Recombiner_Pass(ss_vect<T> &ps);

  void Elem_Pass(ss_vect<double> &ps) { Recombiner_Pass(ps); };
  void Elem_Pass(ss_vect<tps> &ps) { Recombiner_Pass(ps); };

  void print(void);
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

  template<typename T>
  void Solenoid_Pass(ss_vect<T> &ps);

  void Elem_Pass(ss_vect<double> &ps) { Solenoid_Pass(ps); };
  void Elem_Pass(ss_vect<tps> &ps) { Solenoid_Pass(ps); };

  void print(void);
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

  void Elem_Pass(ss_vect<double> &ps);
  void Elem_Pass(ss_vect<tps> &ps);
  void print(void);
};

#endif
