#define PLANES 2

typedef struct globvalrec {
  double        dPcommon,        // dp for numerical differentiation
                dPparticle;      // energy deviation
  double        delta_RF;        // RF acceptance
  Vector2       TotalTune;       // transverse tunes
  double        Omega,
                U0,              // energy lost per turn in keV
                Alphac;          // alphap
  Vector2       Chrom;           // chromaticities
  double        Energy;          // ring energy
  long          Cell_nLoc,       // number of elements
                Elem_nFam,       // number of families
                CODimax;         /* maximum number of cod search before
				    failing */
  double        CODeps;          // precision for closed orbit finder
  psVector        CODvect;         // closed orbit
  int           bpm;             // bpm number
  int           hcorr;           // horizontal corrector number
  int           vcorr;           // vertical corrector number
  int           qt;              // vertical corrector number
  int           gs;              // girder start marker
  int           ge;              // girder end marker
  Matrix        OneTurnMat,      // oneturn matrix
                Ascr,
                Ascrinv,
                Vr,              // real part of the eigenvectors
                Vi;              // imaginal par of the eigenvectors

  bool          Cavity_on,       // if true, cavity turned on
                radiation,       // if true, radiation turned on
                emittance,
                quad_fringe,     // quadrupole hard-edge fringe field.
                H_exact,         // "small ring" Hamiltonian.
                Cart_Bend,
                dip_edge_fudge,  // Dipole edge fudge.
                pathlength,      // absolute path length
                stable,
                Aperture_on,
                EPU;

  double        dE,              // energy loss
                alpha_rad[DOF],  // damping coeffs.
                D_rad[DOF],      // diffusion coeffs (Floquet space)
                J[DOF],          // partition numbers
                tau[DOF];        // damping times
  bool          IBS;             // intrabeam scattering
  double        Qb,              // bunch charge
                D_IBS[DOF];      // diffusion matrix (Floquet space)
  psVector        wr, wi;          // real and imaginary part of eigenvalues
  double        eps[DOF],        // 3 motion invariants
		epsp[DOF],       /* transverse and longitudinal projected
				   emittances */
                alpha_z, beta_z, // longitudinal alpha and beta
                beta0, gamma0;   // Relativistic factors.
  int           RingType;        // 1 if a ring (0 if transfer line)
} globvalrec;


struct DriftType { };


struct MpoleType {
  int         Pmethod;   // Integration Method
  int         PN;        // Number of integration steps
  // Displacement Errors
  Vector2     PdSsys;    // systematic [m]
  Vector2     PdSrms;    // rms [m]
  Vector2     PdSrnd;    // random number
  // Roll angle
  double      PdTpar;    // design [deg]
  double      PdTsys;    // systematic [deg]
  double      PdTrms;    // rms [deg]
  double      PdTrnd;    // random number
  // Multipole strengths
  mpolArray   PBpar;     // design
  mpolArray   PBsys;     // systematic
  mpolArray   PBrms;     // rms
  mpolArray   PBrnd;     // random number
  mpolArray   PB;        // total
  int         Porder;    // The highest order in PB
  int         n_design;  // multipole order (design)
  pthicktype  Pthick;
  // Bending Angles
  double PTx1;           // horizontal entrance angle [deg]
  double PTx2;           // horizontal exit angle [deg]
  double Pgap;           // total magnet gap [m]
  double Pirho;          // 1/rho [1/m]
  double Pc0, Pc1, Ps1;  // corrections for roll error of bend
};

const int  n_harm_max = 10;

struct WigglerType {
  int Pmethod;                // Integration Method
  int PN;                     // number of integration steps
  // Displacement Error
  Vector2 PdSsys;             // systematic [m]
  Vector2 PdSrms;             // rms [m]
  Vector2 PdSrnd;             // random number
  // Roll angle
  double PdTpar;              // design [deg]
  double PdTsys;              // systematic [deg]
  double PdTrms;              // rms [deg]
  double PdTrnd;              // random number
  double lambda;              // lambda
  int    n_harm;              // no of harmonics
  int    harm[n_harm_max];    // harmonic number
  double BoBrhoV[n_harm_max]; // B/Brho vertical
  double BoBrhoH[n_harm_max]; // B/Brho horizontal
  double kxV[n_harm_max];     // kx
  double kxH[n_harm_max];     // kx
  double phi[n_harm_max];     // phi
  mpolArray PBW;
  int Porder;                 // The highest order in PB
};


struct FieldMapType {
  int     n_step;                       // number of integration steps
  int     n[3];                         // no of steps
  int     cut;                          // cut in z direction
  double  scl, phi, x0, Lr, Lm, Ld, L1;
  double  dx[3], *x[3];                 // [dx, dy, dz], [x, y, z]
  double  ***BoBrho[3], ***BoBrho2[3];  // [B_x, B_y, B_z]
  double  ***AoBrho[2], ***AoBrho2[2];  /* [Ax(x, y, z), Ay(x, y, z)],
					   spline info */
};


// ID Laurent
#define IDXMAX 200
#define IDZMAX 100

struct InsertionType {
  int    Pmethod;      // Integration Method
  int    PN;           // number of integration steps
  char   fname1[100];  // Filename for insertion description: first ordre
  char   fname2[100];  // Filename for insertion description: second ordre
  int    nx;           // Horizontal point number
  int    nz;           // Vertical point number
  double scaling;      // static scaling factor as in BETA ESRF
  bool   linear;       // if true linear interpolation else spline
  bool   firstorder;   // true if first order kick map loaded
  bool   secondorder;  // true if second order kick map loaded
  double phi;          // Bend angle.
  double tabx[IDXMAX]; // spacing in H-plane
  double tabz[IDZMAX]; // spacing in V-plane
  double thetax[IDZMAX][IDXMAX], thetax1[IDZMAX][IDXMAX]; // 1 for first order
  double thetaz[IDZMAX][IDXMAX], thetaz1[IDZMAX][IDXMAX];
  bool   long_comp;    // flag for longitudinal comp
  double B2[IDZMAX][IDXMAX]; // B^2_perp
  double **tx, **tz, **f2x, **f2z;
  double **tx1, **tz1, **f2x1, **f2z1; // a voir
  double *tab1, *tab2; // tab of x and z meshes from Radia code

  // Displacement Error
  Vector2 PdSsys;   // systematic [m]
  Vector2 PdSrms;   // rms [m]
  Vector2 PdSrnd;   // random number
  // Roll angle
  double PdTpar;    // design [deg]
  double PdTsys;    // systematic [deg]
  double PdTrms;    // rms [deg]
  double PdTrnd;    // random number
  // Strength
//  double Plperiod;  // Length Period [m]
//  int Pnperiod;    // Number of periods
//  double PBoBrho;   // B/Brho
//  double PKx;       // kx
//  mpolArray PBW;
  int Porder;        // The highest order in PB
};

struct CavityType {
  int    PN;           // Number of integration steps
  double Pvolt;        // Vrf [V]
  double Pfreq;        // Vrf [Hz]
  double phi;          // RF phase
  int    Ph;           // Harmonic number
  bool   entry_focus;  // Edge focusing at entry.
  bool   exit_focus;   // Edge focusing at exit.
};

struct CellType;

const int  Spreader_max = 10;

struct SpreaderType {
  double    E_max[Spreader_max];      // energy levels in increasing order
  CellType  *Cell_ptrs[Spreader_max];
};

struct RecombinerType {
  double    E_min;
  double    E_max;
};

struct SolenoidType {
  int         N;         // Number of integration steps
  // Displacement Errors
  Vector2     PdSsys;    // systematic [m]
  Vector2     PdSrms;    // rms [m]
  Vector2     PdSrnd;    // random number
  // Roll angle
  double      dTpar;     // design [deg]
  double      dTsys;     // systematic [deg]
  double      dTrms;     // rms [deg]
  double      dTrnd;     // random number
  double      BoBrho;    // normalized field strength
};

struct MapType {
  double       dnu[2], alpha[2], beta[2], eta_x, etap_x;
  ss_vect<tps> M;
};

struct elemtype {
  partsName PName;       // Element name.
  double    PL;          // Length[m].
  bool      Reverse;     // Reverse element.
  PartsKind Pkind;       // Enumeration for magnet types.
  union
  {
    DriftType      *D;   // Drift.
    MpoleType      *M;   // Multipole.
    WigglerType    *W;   // Wiggler.
    FieldMapType   *FM;  // Field Map.
    InsertionType  *ID;  // Insertion.
    CavityType     *C;   // Cavity.
    SpreaderType   *Spr; // Spreader.
    RecombinerType *Rec; // Recombiner.
    SolenoidType   *Sol; // Solenoid.
    MapType        *Map; // Map.
  };
};

struct ElemFamType {
  elemtype    ElemF;                // Structure (name, type)
  int         nKid;                 // Kid number
  int         KidList[nKidMax];
  int         NoDBN;
  DBNameType  DBNlist[nKidMax];
};

// LEGO block structure for each element of the lattice
struct CellType {
  int              Fnum;        // Element Family #
  int              Knum;        // Element Kid #
  double           S;           // Position in the ring
  CellType*        next_ptr;    // pointer to next cell (for tracking)
  Vector2          dS,          // Transverse displacement
                   dT;          // dT = (cos(dT), sin(dT))
  elemtype         Elem;        // Structure (name, type)
  Vector2          Nu,          // Phase advances
                   Alpha,       // Alpha functions (redundant)
                   Beta,        // beta fonctions (redundant)
                   Eta, Etap;   // dispersion and its derivative (redundant)
  psVector         BeamPos;     // Last position of the beam this cell
  Matrix           A,           // Floquet space to phase space transformation
                   sigma;       // sigma matrix (redundant)
  Vector2          maxampl[PLANES]; /* Horizontal and vertical physical
				       apertures:
				         maxampl[X_][0] < x < maxampl[X_][1]
					 maxampl[Y_][0] < y < maxampl[Y_][1] */
};
