#ifndef TRACY_GLOBAL_H
#define TRACY_GLOBAL_H

#define PLANES 2

typedef struct globvalrec {
  bool
    Cavity_on,      // If true, cavity turned on.
    radiation,      // If true, radiation turned on.
    emittance,
    quad_fringe,    // Quadrupole hard-edge fringe field.
    H_exact,        // "Small Ring" Hamiltonian.
    Cart_Bend,
    dip_edge_fudge, // Dipole edge fudge.
    pathlength,     // Absolute path length.
    stable,
    Aperture_on,
    EPU,
    mat_meth,       // Matrix method.
    IBS,            // Intrabeam scattering.
    rad_D;          // Compute radiation diffusion matrix.
  long
    Cell_nLoc,      // Number of elements.
    Elem_nFam,      // Number of families.
    CODimax;        // Closed orbit finder: max number of iterations.
  int
    bpm,            // BPM number.
    hcorr, vcorr,   // Corrector: horizontal number, vertical number.
    qt,             // Vertical corrector number.
    gs, ge,         // Girder: [start, end] marker.
    RingType;       // Ring - 1, transfer line - 0.
  double
    dPcommon,       // Dp for numerical differentiation.
    dPparticle,     // Momentum deviation.
    delta_RF,       // RF Acceptance.
    TotalTune[3],   // Tunes.
    Omega,          // Synchrotron frequency.
    U0,             // Energy loss per turn [keV].
    Alphac,         // Linear momentum compaction.
    Energy,         // Beam energy.
    CODeps,         // Precision for fixed point finder.
    dE,             // Energy loss.
    alpha_rad[DOF], // Damping coeffs.
    D_rad[DOF],     // Diffusion coeffs, Floquet Space.
    J[DOF],         // Partition numbers.
    tau[DOF],       // Damping times.
    eps[DOF],       // Eigen emittances.
    epsp[DOF],      // Projected emittances.
    alpha_z,        // Longitudinal alpha and beta.
    beta_z,
    beta0,          // Relativistic factors.
    gamma0,
    Qb,             // Bunch charge.
    D_IBS[DOF];     // IBS diffusion matrix, Floquet space.
  Vector2
    Chrom;          // Linear chromaticity.
  psVector
    CODvect,        // Closed orbit.
    wr, wi;         // Eigenvalues: [real, imaginary].
  Matrix
    OneTurnMat,     // Linear Poincare map.
    Ascr,
    Ascrinv,
    Vr, Vi;         // Eigenvectors: [real, imaginary].
  ss_vect<tps>
    Diff_mat;       // Diffusion matrix, phase-space.
} globvalrec;


struct DriftType { };


struct MpoleType {
  int
    Pmethod,          // Integration method.
    PN;               // Number of integration steps.
  // Displacement Errors.
  Vector2
    PdSsys,           // Systematic [m].
    PdSrms,           // RMS [m].
    PdSrnd;           // Random number.
  // Roll angle.
  double
    PdTpar,           // Design [deg].
    PdTsys,           // Systematic [deg].
    PdTrms,           // RMS [deg].
    PdTrnd;           // Random number.
  // Multipole strengths.
  mpolArray
    PBpar,            // Design.
    PBsys,            // Systematic.
    PBrms,            // RMS.
    PBrnd,            // Random number.
    PB;               // Total.
  int
    Porder,           // Max multipole order.
    n_design;         // Multipole order, design.
  pthicktype
    Pthick;
  // Bending Angles.
  double
    PTx1,             // Horizontal entrance angle [deg].
    PTx2,             // Horizontal exit angle [deg].
    Pgap,             // Total magnet gap [m].
    Pirho,            // 1/rho [1/m].
    Pc0, Pc1, Ps1;    // Bend roll error corrections.
  ss_vect<tps>
    M_lin,            // Linear map.
    dD_mat;           // Diffusion matrix increment, phase-space.
};

const int n_harm_max = 10;

struct WigglerType {
  int
    Pmethod,             // Integration Method.
    PN;                  // number of integration steps.
  // Displacement Error.
  Vector2
    PdSsys,              // Systematic [m].
    PdSrms,              // RMS [m].
    PdSrnd;              // Random number.
  // Roll angle
  double
    PdTpar,              // Design [deg].
    PdTsys,              // Systematic [deg].
    PdTrms,              // RMS [deg].
    PdTrnd,              // Random number.
    Lambda;              // lambda.
  int
    n_harm,              // Number of harmonics.
    harm[n_harm_max];    // Harmonic number.
  double
    BoBrhoV[n_harm_max], // B/Brho vertical.
    BoBrhoH[n_harm_max], // B/Brho horizontal.
    kxV[n_harm_max],     // Kx.
    kxH[n_harm_max],     // Kx.
    phi[n_harm_max];     // Phi.
  mpolArray
    PBW;
  int
    Porder;              // Max multipole order.
};


struct FieldMapType {
  int
    n_step,              // number of integration steps.
    n[3],                // no of steps.
    cut;                 // cut in z direction.
  double
    scl,
    phi,
    x0,
    Lr,
    Lm,
    Ld,
    L1,
    dx[3],
    *x[3],               // [dx, dy, dz], [x, y, z].
    ***BoBrho[3],
    ***BoBrho2[3],       // [B_x, B_y, B_z].
    ***AoBrho[2],
    ***AoBrho2[2];       // [Ax(x, y, z), Ay(x, y, z)], spline info.
};


// ID Laurent
#define IDXMAX 200
#define IDZMAX 100

struct InsertionType {
  int
    Pmethod,                 // Integration Method.
    PN;                      // number of integration steps.
  char
    fname1[100],             // Filename for insertion description: 1st order.
    fname2[100];             // Filename for insertion description: 2nd order.
  int
    nx,                      // Horizontal point number.
    nz;                      // Vertical point number.
  double
    scaling;                 // static scaling factor as in BETA, ESRF.
  bool
    linear,                  // If true linear interpolation else spline.
    firstorder,              // True if first order kick map loaded.
    secondorder;             // True if second order kick map loaded.
  double
    phi,                     // Bend angle.
    tabx[IDXMAX],            // Spacing in H-plane.
    tabz[IDZMAX],            // Spacing in V-plane.
    thetax[IDZMAX][IDXMAX],
    thetax1[IDZMAX][IDXMAX], // First order - 1.
    thetaz[IDZMAX][IDXMAX],
    thetaz1[IDZMAX][IDXMAX];
  bool
    long_comp;               // flag for longitudinal comp.
  double
    B2[IDZMAX][IDXMAX],      // B^2_perp.
    **tx,
    **tz,
    **f2x,
    **f2z,
    **tx1,
    **tz1,
    **f2x1,
    **f2z1,                  // a voir.
    *tab1,
    *tab2;                   // tab of x and z meshes from Radia code.
  // Displacement Error.
  Vector2
    PdSsys,                  // systematic [m].
    PdSrms,                  // rms [m].
    PdSrnd;                  // random number.
  // Roll angle.
  double
    PdTpar,                  // design [deg].
    PdTsys,                  // systematic [deg].
    PdTrms,                  // rms [deg].
    PdTrnd;                  // random number.
//  Strength
//  double Plperiod;         // Length Period [m].
//  int Pnperiod;            // Number of periods.
//  double PBoBrho;          // B/Brho.
//  double PKx;              // kx.
//  mpolArray PBW;
  int
    Porder;                  // The highest order in PB.
};

struct CavityType {
  bool
    entry_focus,       // Edge focusing at entry.
    exit_focus;        // Edge focusing at exit.
  int
    N_step,            // Number of integration steps.
    harm_num;          // RF harmonic number.
  double
  // Fundamental.
    f_RF,              // RF frequency [Hz].
    V_RF,              // RF voltage [V].
    phi_RF,            // RF phase [deg].
  // For beam loading compensation for the fundamental mode.
    R_sh_fund,         // Shunt impedance.
    Q_RF,              // Quality factor.
    beta_RF,           /* Power coupling factor, ring definition:
			    P = V^2/(2*R_s).                                  */
    coord_scl[3]
    = {1.0, 1.0, 1.0}; // Scale factor for coordinates [x, y, z] for HOMs.

  // Higher order modes.
  std::vector<double>
  // Transverse.
    HOM_f_trans[2],    // Frequency
    HOM_R_sh_trans[2], // Shunt impedance.
    HOM_Q_trans[2],    // Quality factor.
    beta_trans[2],     /* Power coupling factor, ring definition:
			    P = V^2/(2*R_s).                                  */
  // Longitudinal.
    HOM_f_long,        // Frequency
    HOM_R_sh_long,     // Shunt impedance.
    HOM_Q_long,        // Quality factor.
    beta_long;         /* Power coupling factor, ring definition:
			    P = V^2/(2*R_s).                                  */
  std::vector< std::complex<double> >
  // Transverse.
    HOM_V_trans[2],    // Amplitude & phase.
  // Longitudinal.
    HOM_V_long;        // Amplitude & phase.
};

struct CellType;

const int Spreader_max = 10;

struct SpreaderType {
  double
    E_max[Spreader_max];      // energy levels in increasing order.
  CellType
    *Cell_ptrs[Spreader_max];
};

struct RecombinerType {
  double
    E_min,
    E_max;
};

struct SolenoidType {
  int
    N;         // Number of integration steps.
  // Displacement Errors.
  Vector2
    PdSsys,    // systematic [m].
    PdSrms,    // rms [m].
    PdSrnd;    // random number.
  // Roll angle.
  double
    dTpar,     // design [deg].
    dTsys,     // systematic [deg].
    dTrms,     // rms [deg].
    dTrnd,     // random number.
    BoBrho;    // normalized field strength.
};

struct MapType {
  double
    dnu[2], alpha[2], beta[2], eta_x, etap_x;
  ss_vect<tps>
    M;
};

struct elemtype {
  partsName
    PName;               // Element name.
  double
    PL;                  // Length[m].
  bool
    Reverse;             // Reverse element.
  PartsKind
    Pkind;               // Enumeration for magnet types.
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
  elemtype
    ElemF;            // Structure (name, type).
  int
    nKid,             // Kid number.
    KidList[nKidMax],
    NoDBN;
  DBNameType
    DBNlist[nKidMax];
};

// LEGO Block Structure for each Element of the Lattice.
struct CellType {
  int
    Fnum,            // Element Family #.
    Knum;            // Element Kid #.
  double
    S,               // Position in the ring.
    curly_dH_x,      // Contribution to curly_H_x.
    dI[6],           // Contribution to I[2..5].
    dD[3];           /* Increamental change of the linear actions by the
		        radiation reaction                                    */
  CellType*
    next_ptr;        // pointer to next cell (for tracking).
  Vector2
    dS,              // Transverse displacement.
    dT,              // dT = (cos(dT), sin(dT)).
    Nu,              // Phase advances.
    Alpha,           // Alpha functions (redundant).
    Beta,            // beta fonctions (redundant).
    Eta,             // dispersion and its derivative (redundant).
    Etap,
    maxampl[PLANES]; /* Horizontal and vertical physical apertures:
			maxampl[X_][0] < x < maxampl[X_][1]
			maxampl[Y_][0] < y < maxampl[Y_][1]. */
  elemtype
    Elem;            // Structure (name, type).
  psVector
    BeamPos;         // Last position of the beam this cell.
  Matrix
    A,               // Floquet space to phase space transformation.
    sigma;           // sigma matrix (redundant).
};

#endif
