/* Tracy-2:

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version,
                 SLS, PSI      1995 - 1997,
   M. Boege      SLS, PSI      1998          Pascal to C translation,
   L. Nadolski   SOLEIL        2002          Link to NAFF and Radia field maps,
   J. Bengtsson  NSLS-II, BNL  2004 - 2015,
   J. Bengtsson                2017          C to C++ re-factoring.           */


#define PLANES 2


// Max number of LEGO blocks (Cell_nLoc).
#define Cell_nLocMax 20000

// Max number of families for Elem_NFam.
#define Elem_nFamMax 3000

// Max number of kids.
#define nKidMax 5000

// ID Laurent.
#define IDXMAX 200
#define IDZMAX 100

const int max_elem     = Cell_nLocMax;
const int n_harm_max   = 10;
const int Spreader_max = 10;


// LEGO Element pieces for Lattice.

class DriftType;
class MpoleType;
class WigglerType;
class InsertionType;
class FieldMapType;
class CavityType;
class SpreaderType;
class RecombinerType;
class SolenoidType;

class elemtype {
 private:
 public:
  partsName Name;        // Element name.
  double    L;           // Length[m].
  bool      Reverse;     // Reverse element.
  PartsKind Kind;        // Enumeration for magnet types.
  union
  {
    DriftType      *D;   // Drift.
    MpoleType      *M;   // Multipole.
    WigglerType    *W;   // Wiggler.
    InsertionType  *ID;  // Kick Map.
    FieldMapType   *FM;  // Field Map.
    SolenoidType   *Sol; // Solenoid.
    CavityType     *C;   // Cavity.
    SpreaderType   *Spr; // Spreader.
    RecombinerType *Rec; // Recombiner.
  };
};


class CellType {
 private:
 public:
  int       Fnum;            // Element Family #.
  int       Knum;            // Element Kid #.
  double    S;               // Position in the ring.
  CellType* next_ptr;        // pointer to next cell (for tracking).
  Vector2   dS,              // Transverse displacement.
            dT;              // dT = (cos(dT), sin(dT)).
  elemtype  Elem;            // Structure (name, type).
  Vector2   Nu,              // Phase advances.
            Alpha,           // Alpha functions (redundant).
            Beta,            // beta fonctions (redundant).
            Eta, Etap;       // dispersion and its derivative (redundant).
  psVector  BeamPos;         // Last position of the beam this cell.
  Matrix    A,               // Floquet space to phase space transformation.
            sigma;           // sigma matrix (redundant).
  Vector2   maxampl[PLANES]; /* Horizontal and vertical physical apertures:
				  maxampl[X_][0] < x < maxampl[X_][1]
				  maxampl[Y_][0] < y < maxampl[Y_][1]. */

  void Cell_Init(void);

  void Elem_Print(FILE *f, int Fnum1);

  template<typename T>
  void GtoL(ss_vect<T> &X, const Vector2 &S, const Vector2 &R,
	    const double c0, const double c1, const double s1);

  template<typename T>
  void LtoG(ss_vect<T> &X, const Vector2 &S, const Vector2 &R,
	    const double c0, const double c1, const double s1);
    
  template<typename T>
  void Marker_Pass(CellType &Cell, ss_vect<T> &X);

  template<typename T>
  void Elem_Pass(const long i, ss_vect<T> &x);

  template<typename T>
  void Cell_Pass(const long i0, const long i1, ss_vect<T> &x, long &lastpos);

  void Cell_Pass(const long i0, const long i1, tps &sigma, long &lastpos);

  bool Cell_getCOD(long imax, double eps, double dP, long &lastpos);

  bool GetCOD(long imax, double eps, double dP, long &lastpos);
};


class DriftType {
 private:
 public:
  void Drift_Init(int Fnum1);

  void Drift_Print(FILE *f, int Fnum1);

  template<typename T>
  void Drift_Pass(CellType &Cell, ss_vect<T> &x);
};


class MpoleType {
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

  void Mpole_Init(int Fnum1);

  void Mpole_Print(FILE *f, int Fnum1);

  template<typename T>
  void Mpole_Pass(CellType &Cell, ss_vect<T> &x);
};


class WigglerType {
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

  void Wiggler_Init(int Fnum1);

  void Wiggler_Print(FILE *f, int Fnum1);

  template<typename T>
  void Wiggler_Pass(CellType &Cell, ss_vect<T> &X);
};


class InsertionType {
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

  void Insertion_Init(int Fnum1);
  
  void Insertion_Print(FILE *f, int Fnum1);

  template<typename T>
  void Insertion_Pass(CellType &Cell, ss_vect<T> &x);
};


class FieldMapType {
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

  void FieldMap_Init(int Fnum1);

  template<typename T>
  void FieldMap_Pass(CellType &Cell, ss_vect<T> &ps);
};


class SolenoidType {
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

  void Solenoid_Init(int Fnum1);

  template<typename T>
  void Solenoid_Pass(CellType &Cell, ss_vect<T> &ps);
};


class CavityType {
 private:
 public:
  int    N;           // Number of integration steps.
  double volt;        // Vrf [V].
  double freq;        // Vrf [Hz].
  double phi;         // RF phase.
  int    h;           // Harmonic number.
  bool   entry_focus; // Edge focusing at entry.
  bool   exit_focus;  // Edge focusing at exit.

  void Cav_Init(int Fnum1);

  template<typename T>
  void Cav_Pass(CellType &Cell, ss_vect<T> &X);
};


class SpreaderType {
 private:
 public:
  double   E_max[Spreader_max];      // energy levels in increasing order
  CellType *Cell_ptrs[Spreader_max];

  void Spreader_Init(int Fnum1);
};


class RecombinerType {
 private:
 public:
  double E_min;
  double E_max;

  void Recombiner_Init(int Fnum1);
};


class ElemFamType {
 private:
 public:
  elemtype   ElemF;            // Structure (name, type).
  int        nKid;             // Kid number.
  int        KidList[nKidMax];
  int        NoDBN;
  DBNameType DBNlist[nKidMax];
};


struct Lattice_Param {
  double   dPcommon,        // dp for numerical differentiation.
           dPparticle;      // energy deviation.
  double   delta_RF;        // RF acceptance.
  Vector2  TotalTune;       // transverse tunes.
  double   Omega,
           U0,              // energy lost per turn in keV.
           Alphac;          // alphap.
  Vector2  Chrom;           // chromaticities.
  double   Energy;          // ring energy.
  long     Cell_nLoc,       // number of elements.
           Elem_nFam,       // number of families.
           CODimax;         // Max number of cod search before failing.
  double   CODeps;          // precision for closed orbit finder.
  psVector CODvect;         // closed orbit.
  int      bpm;             // bpm number.
  int      hcorr;           // horizontal corrector number.
  int      vcorr;           // vertical corrector number.
  int      qt;              // vertical corrector number.
  int      gs;              // girder start marker.
  int      ge;              // girder end marker.
  Matrix   OneTurnMat,      // oneturn matrix.
           Ascr,
           Ascrinv,
           Vr,              // real part of the eigenvectors.
           Vi;              // imaginal par of the eigenvectors.

  bool     Cavity_on,       // if true, cavity turned on.
           radiation,       // if true, radiation turned on.
           emittance,
           dip_fringe,      // dipole hard-edge fringe field.
           quad_fringe,     // quadrupole hard-edge fringe field.
           H_exact,         // "small ring" Hamiltonian.
           pathlength,      // absolute path length.
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
