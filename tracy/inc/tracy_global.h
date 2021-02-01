#ifndef TRACY_GLOBAL_H
#define TRACY_GLOBAL_H

// ID Laurent
#define IDXMAX 200
#define IDZMAX 100

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
    dI[6],                     // Contribution to I[2..5].
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
class elemtype : public CellType {
 public:
  partsName
    PName;                     // Element name.
  bool
    Reverse;                   // Reverse element.
  double
    PL;                        // Length[m].
  PartsKind
    Pkind;                     // Enumeration for magnet types.

  virtual void Elem_Pass(ss_vect<double> &x);
  virtual void Elem_Pass(ss_vect<tps> &x);
};

// Index for lattice elements.
class ElemFamType {
 public:
  int
    nKid,                      // Kid number.
    KidList[nKidMax],
    NoDBN;
  elemtype
    *ElemF;                    // Structure (name, type).
  DBNameType
    DBNlist[nKidMax];
};

class DriftType : public elemtype {
  template<typename T>
  void Elem_Pass(ss_vect<T> &x);
};

class MpoleType : public elemtype {
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
  mpolArray
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

  template<typename T>
  void Elem_Pass(ss_vect<T> &x);
};

class CavityType : public elemtype {
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

  template<typename T>
  void Elem_Pass(ss_vect<T> &x);
};

class MarkerType : public elemtype {
  template<typename T>
  void Elem_Pass(ss_vect<T> &x);
};

class WigglerType : public elemtype {
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
  mpolArray
  PBW;

  template<typename T>
  void Elem_Pass(ss_vect<T> &x);
};

class InsertionType : public elemtype {
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

  template<typename T>
  void Elem_Pass(ss_vect<T> &x);
};

class FieldMapType : public elemtype {
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

  template<typename T>
  void Elem_Pass(ss_vect<T> &x);
};

class CellType;

class SpreaderType : public elemtype {
 public:
  double
    E_max[Spreader_max];       // energy levels in increasing order.
  CellType
    *Cell_ptrs[Spreader_max];

  template<typename T>
  void Elem_Pass(ss_vect<T> &x);
};

class RecombinerType : public elemtype {
 public:
  double
    E_min,
    E_max;

  template<typename T>
  void Elem_Pass(ss_vect<T> &x);
};

class SolenoidType : public elemtype {
 public:
  int
    N;                         // Number of integration steps.
  // Roll angle.
  double
    dTpar,                     // design [deg].
    dTsys,                     // systematic [deg].
    dTrms,                     // rms [deg].
    dTrnd,                     // random number.
    BoBrho,                    // normalized field strength.
  // Displacement Errors.
    PdSsys[2],                 // systematic [m].
    PdSrms[2],                 // rms [m].
    PdSrnd[2];                 // random number.

  template<typename T>
  void Elem_Pass(ss_vect<T> &x);
};

class MapType : public elemtype {
 public:
  double
    dnu[2],
    alpha[2],
    beta[2],
    eta_x,
    etap_x;
  ss_vect<tps>
    M;

  void Elem_Pass(ss_vect<double> &x);
  void Elem_Pass(ss_vect<tps> &x);
};

#endif
