#define PLANES 2

typedef struct globvalrec {
  double   dPcommon,        // dp for numerical differentiation.
           dPparticle,      // Energy deviation.
           delta_RF,        // RF Acceptance.
           TotalTune[3];    // Transverse Tunes.
  double   Omega,           // Synchrotron Frequency.
           U0,              // Energy Loss per turn [keV].
           Alphac;          // Linear Momentum Compaction.
  Vector2  Chrom;           // Linear Chromaticities.
  double   Energy;          // Beam Energy.
  long     Cell_nLoc,       // Number of Elements.
           Elem_nFam,       // Number of Families.
           CODimax;         // closed Orbit Finder: max number of iterations,
  double   CODeps;          //                      precision.
  psVector CODvect;         // Closed Orbit.
  int      bpm,             // BPM Number.
           hcorr,           // Corrector: Horizontal number,
           vcorr,           //            Vertical number.
           qt,              // Vertical corrector number.
           gs,              // Girder: start marker,
           ge;              //         end marker.
  Matrix   OneTurnMat,      // Linear Poincare Map.
           Ascr,
           Ascrinv,
           Vr,              // Eigenvectors: Real part, 
           Vi;              //               Imaginary part.

  bool     Cavity_on,       // if true, cavity turned on
           radiation,       // if true, radiation turned on
           emittance,
           quad_fringe,     // quadrupole hard-edge fringe field.
           H_exact,         // "Small Ring" Hamiltonian.
           Cart_Bend,
           dip_edge_fudge,  // Dipole Edge fudge.
           pathlength,      // Absolute Path Length.
           stable,
           Aperture_on,
           EPU;

  double   dE,              // Energy Loss.
           alpha_rad[DOF],  // Damping Coeffs.
           D_rad[DOF],      // Diffusion Coeffs (Floquet Space).
           J[DOF],          // Partition Numbers.
           tau[DOF];        // Damping Times.
  bool     IBS;             // Intrabeam Scattering.
  double   Qb,              // Bunch Charge.
           D_IBS[DOF];      // Diffusion Matrix (Floquet Ipace).
  psVector wr, wi;          // Eigenvalues: Real and Imaginary part.
  double   eps[DOF],        // Eigenemittances.
           epsp[DOF],       // Trans. & Long. projected Emittances.
           alpha_z, beta_z, // Long. alpha and beta.
           beta0, gamma0;   // Relativistic factors.
  int      RingType;        // 1 if a ring (0 if transfer line)
} globvalrec;


struct DriftType { };


struct MpoleType {
  int         Pmethod,   // Integration Method
              PN;        // Number of integration steps
  // Displacement Errors
  Vector2     PdSsys,    // systematic [m]
              PdSrms,    // rms [m]
              PdSrnd;    // random number
  // Roll angle
  double      PdTpar,    // design [deg]
              PdTsys,    // systematic [deg]
              PdTrms,    // rms [deg]
              PdTrnd;    // random number
  // Multipole strengths
  mpolArray   PBpar,     // design
              PBsys,     // systematic
              PBrms,     // rms
              PBrnd,     // random number
              PB;        // total
  int         Porder,    // The highest order in PB
              n_design;  // multipole order (design)
  pthicktype  Pthick;
  // Bending Angles
  double      PTx1,           // horizontal entrance angle [deg]
              PTx2,           // horizontal exit angle [deg]
              Pgap,           // total magnet gap [m]
              Pirho,          // 1/rho [1/m]
              Pc0, Pc1, Ps1;  // corrections for roll error of bend
};

const int n_harm_max = 10;

struct WigglerType {
  int       Pmethod,             // Integration Method
            PN;                  // number of integration steps
  // Displacement Error
  Vector2   PdSsys,              // Systematic [m].
            PdSrms,              // RMS [m].
            PdSrnd;              // Random number.
  // Roll angle
  double    PdTpar,              // Design [deg].
            PdTsys,              // Systematic [deg].
            PdTrms,              // RMS [deg].
            PdTrnd,              // Random number.
            lambda;              // lambda.
  int       n_harm,              // No of harmonics.
            harm[n_harm_max];    // Harmonic number.
  double    BoBrhoV[n_harm_max], // B/Brho vertical.
            BoBrhoH[n_harm_max], // B/Brho horizontal.
            kxV[n_harm_max],     // kx.
            kxH[n_harm_max],     // kx.
            phi[n_harm_max];     // phi.
  mpolArray PBW;
  int       Porder;              // The highest order in PB.
};


struct FieldMapType {
  int    n_step,                       // number of integration steps
         n[3],                         // no of steps
         cut;                          // cut in z direction
  double scl, phi, x0, Lr, Lm, Ld, L1,
         dx[3], *x[3],                 // [dx, dy, dz], [x, y, z]
         ***BoBrho[3], ***BoBrho2[3],  // [B_x, B_y, B_z]
         ***AoBrho[2], ***AoBrho2[2];  /* [Ax(x, y, z), Ay(x, y, z)],
					  spline info */
};


// ID Laurent
#define IDXMAX 200
#define IDZMAX 100

struct InsertionType {
  int    Pmethod,      // Integration Method
         PN;           // number of integration steps
  char   fname1[100],  // Filename for insertion description: first ordre
         fname2[100];  // Filename for insertion description: second ordre
  int    nx,           // Horizontal point number
         nz;           // Vertical point number
  double scaling;      // static scaling factor as in BETA ESRF
  bool   linear,       // if true linear interpolation else spline
         firstorder,   // true if first order kick map loaded
         secondorder;  // true if second order kick map loaded
  double phi,          // Bend angle.
         tabx[IDXMAX], // spacing in H-plane
         tabz[IDZMAX], // spacing in V-plane
         thetax[IDZMAX][IDXMAX], thetax1[IDZMAX][IDXMAX], // 1 for first order
    thetaz[IDZMAX][IDXMAX], thetaz1[IDZMAX][IDXMAX];
  bool   long_comp;    // flag for longitudinal comp
  double B2[IDZMAX][IDXMAX], // B^2_perp
         **tx, **tz, **f2x, **f2z,
         **tx1, **tz1, **f2x1, **f2z1, // a voir
         *tab1, *tab2; // tab of x and z meshes from Radia code

  // Displacement Error
  Vector2 PdSsys,   // systematic [m]
          PdSrms,   // rms [m]
          PdSrnd;   // random number
  // Roll angle
  double  PdTpar,    // design [deg]
          PdTsys,    // systematic [deg]
          PdTrms,    // rms [deg]
          PdTrnd;    // random number
  // Strength
//  double Plperiod;  // Length Period [m]
//  int Pnperiod;    // Number of periods
//  double PBoBrho;   // B/Brho
//  double PKx;       // kx
//  mpolArray PBW;
  int     Porder;        // The highest order in PB
};

struct CavityType {
  int    PN;           // Number of integration steps
  double Pvolt,        // Vrf [V]
         Pfreq,        // Vrf [Hz]
         phi;          // RF phase
  int    Ph;           // Harmonic number
  bool   entry_focus,  // Edge focusing at entry.
         exit_focus;   // Edge focusing at exit.
};

struct CellType;

const int Spreader_max = 10;

struct SpreaderType {
  double   E_max[Spreader_max];      // energy levels in increasing order
  CellType *Cell_ptrs[Spreader_max];
};

struct RecombinerType {
  double E_min,
         E_max;
};

struct SolenoidType {
  int     N;         // Number of integration steps
  // Displacement Errors
  Vector2 PdSsys,    // systematic [m]
          PdSrms,    // rms [m]
          PdSrnd;    // random number
  // Roll angle
  double  dTpar,     // design [deg]
          dTsys,     // systematic [deg]
          dTrms,     // rms [deg]
          dTrnd,     // random number
          BoBrho;    // normalized field strength
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
  elemtype   ElemF;                // Structure (name, type)
  int        nKid,                 // Kid number
             KidList[nKidMax],
             NoDBN;
  DBNameType DBNlist[nKidMax];
};

// LEGO Block Structure for each Element of the Lattice.
struct CellType {
  int       Fnum,        // Element Family #
            Knum;        // Element Kid #
  double    S;           // Position in the ring
  CellType* next_ptr;    // pointer to next cell (for tracking)
  Vector2   dS,          // Transverse displacement
            dT;          // dT = (cos(dT), sin(dT))
  elemtype  Elem;        // Structure (name, type)
  Vector2   Nu,          // Phase advances
            Alpha,       // Alpha functions (redundant)
            Beta,        // beta fonctions (redundant)
            Eta, Etap;   // dispersion and its derivative (redundant)
  psVector  BeamPos;     // Last position of the beam this cell
  Matrix    A,           // Floquet space to phase space transformation
            sigma;       // sigma matrix (redundant)
  Vector2   maxampl[PLANES]; /* Horizontal and vertical physical
				apertures:
				maxampl[X_][0] < x < maxampl[X_][1]
				maxampl[Y_][0] < y < maxampl[Y_][1] */
  double    curly_dH_x,  // Contribution to curly_H_x.
            dI[6];       // Contribution to I[2..5].
};
