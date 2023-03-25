#define NO 2

#include <random>
#include <complex>

#include "tracy_lib.h"


int
  no_tps   = NO,
  ndpt_tps = 5;


class MomentType {
private:
  double       q;     // Charge.
public:
  ss_vect<tps> Sigma; // Statistical moments for charge distribution.

  void propagate_cav(void);
  void propagate_lat(void);
  void propagate(void);
};

  
void  MomentType::propagate_cav(void)
{
  long int lastpos;

  const int loc = globval.Cell_nLoc;

  printf("\n  From: %10s\n", Cell[loc].Elem.PName);
  Cell_Pass(loc, loc, Sigma, lastpos);
}


void  MomentType::propagate_lat(void)
{
  long int lastpos;

  printf("\n  From: %10s\n  To:   %10s\n", Cell[0].Elem.PName,
	 Cell[globval.Cell_nLoc-1].Elem.PName);
  Cell_Pass(0, globval.Cell_nLoc-1, Sigma, lastpos);
}


void  MomentType::propagate(void)
{
  // Assumes that the RF cavity is at the end of the lattice.

  this->propagate_lat();
  this->propagate_cav();
}


void prt_vec(const int n_dof, const string &str, const ss_vect<double> vec)
{
  const int n_dec = 6;

  printf("%s\n", str.c_str());
  for (int k = 0; k < 2*n_dof; k++)
    printf("%*.*e", n_dec+8, n_dec, vec[k]);
  printf("\n");
}


void prt_map(const int n_dof, const string &str, const ss_vect<tps> map)
{
  const int n_dec = 6;

  printf("%s\n", str.c_str());
  for (int i = 0; i < 2*n_dof; i++) {
    for (int j = 0; j < 2*n_dof; j++)
      printf("%*.*e", n_dec+8, n_dec, map[i][j]);
    printf("\n");
  }
}


ss_vect<tps> tp_map(const int n_dof, const ss_vect<tps> &A)
{
  // Matrix transpose.
  Matrix A_mat;

  getlinmat(2*n_dof, A, A_mat);
  TpMat(2*n_dof, A_mat);
  return putlinmat(2*n_dof, A_mat);
}


ss_vect<tps> compute_M(const ss_vect<double> fix_point)
{
  // Compute the Poincaré map for the fix point.
  long int lastpos;
  ss_vect<tps> M;

  M.identity();
  M += fix_point;
  Cell_Pass(0, globval.Cell_nLoc, M, lastpos);
  M -= fix_point;
  return M;
}


ss_vect<tps> compute_A
(const ss_vect<tps> &M, double alpha_rad[], ss_vect<tps> &R)
{
  // Poincaré map Diagonalization.
  // The damping coeffients alpha[] are obtained from the eigen values.
  double       nu_s, alpha_c;
  Matrix       M_mat, A_mat, A_inv_mat, R_mat;
  ss_vect<tps> A;

  const int    n_dof = 3;
  const double C     = Cell[globval.Cell_nLoc].S;

  globval.Cavity_on = globval.radiation = true;

  getlinmat(2*n_dof, M, M_mat);
  GDiag(2*n_dof, C, A_mat, A_inv_mat, R_mat, M_mat, nu_s, alpha_c);
  for (int k = 0; k < n_dof; k++)
    alpha_rad[k] = globval.alpha_rad[k];

  A = putlinmat(2*n_dof, A_mat);
  if (n_dof < 3) {
    // Coasting longitudinal plane.
    A[ct_] = tps(0e0, ct_+1);
    A[delta_] = tps(0e0, delta_+1);
  }
  R = putlinmat(2*n_dof, R_mat);
  return A;
}


void compute_tau(const double alpha_rad[], double tau[])
{
  // Compute the damping times.
  
  const int    n_dof = 3;
  const double C     = Cell[globval.Cell_nLoc].S;

  for (int k = 0; k < n_dof; k++)
    tau[k] = -C/(c0*alpha_rad[k]);
}


void compute_D
(const ss_vect<double> &fix_point, const ss_vect<tps> &A, double D[])
{
  // Compute the diffusion coefficients.

  int long     lastpos;
  ss_vect<tps> As, D_diag;

  const int n_dof = 3;

  globval.Cavity_on = globval.radiation = globval.emittance = true;

  As = A + fix_point;
  Cell_Pass(0, globval.Cell_nLoc, As, lastpos);

  globval.emittance = false;

  for (int k = 0; k < n_dof; k++)
    D[k] = globval.D_rad[k];
}

void compute_eps(const double tau[], const double D[], double eps[])
{
  // Compute the eigen emittances.

  const int    n_dof = 3;
  const double C     = Cell[globval.Cell_nLoc].S;

  for (int k = 0; k < n_dof; k++)
    eps[k] = D[k]*tau[k]*c0/(2e0*C);
}


void compute_bunch_size
(const ss_vect<tps> &A, const double eps[], double &sigma_s,
 double &sigma_delta)
{
  // Compute the bunch size: sigma_s & sigma_delta.
  double alpha_s, beta_s, gamma_s;

  alpha_s = -A[ct_][ct_]*A[delta_][ct_] - A[ct_][delta_]*A[delta_][delta_];
  beta_s = sqr(A[ct_][ct_]) + sqr(A[ct_][delta_]);
  gamma_s = (1e0+sqr(alpha_s))/beta_s;

  sigma_s = sqrt(beta_s*eps[Z_]);
  sigma_delta = sqrt(gamma_s*eps[Z_]);
}


ss_vect<tps> compute_M_diff(const ss_vect<tps> &A, const double D[])
{
  // Compute the diffusion matrix.
  ss_vect<tps> D_diag;

  const int n_dof = 3;

  D_diag.zero();
  for (int k = 0; k < 2*n_dof; k++)
    D_diag[k] = D[k/2]*tps(0e0, k+1);
  return A*D_diag*tp_map(n_dof, A);
}


ss_vect<tps> compute_M_Chol_t(const ss_vect<tps> &M_diff)
{
  // Compute the Cholesky decomposition: D = L^T L.
  int          j, k, j1, k1;
  double       *diag, **d1, **d2;
  ss_vect<tps> M_Chol_t;

  const int
    n_dof = 3,
    n     = 2*n_dof;

  diag = dvector(1, n);
  d1 = dmatrix(1, n, 1, n);
  d2 = dmatrix(1, n, 1, n);

  // Remove vertical plane.
  for (j = 1; j <= 4; j++) {
    j1 = (j < 3)? j : j+2;
    for (k = 1; k <= 4; k++) {
      k1 = (k < 3)? k : k+2;
      d1[j][k] = M_diff[j1-1][k1-1];
    }
  }

  dcholdc(d1, 4, diag);

  for (j = 1; j <= 4; j++)
    for (k = 1; k <= j; k++)
      d2[j][k] = (j == k)? diag[j] : d1[j][k];

  M_Chol_t.zero();
  for (j = 1; j <= 4; j++) {
    j1 = (j < 3)? j : j+2;
    for (k = 1; k <= 4; k++) {
      k1 = (k < 3)? k : k+2;
      M_Chol_t[j1-1] += d2[j][k]*tps(0e0, k1);
    }
  }

  free_dvector(diag, 1, n);
  free_dmatrix(d1, 1, n, 1, n);
  free_dmatrix(d2, 1, n, 1, n);

  return M_Chol_t;
}


void compute_maps
(ss_vect<tps> &M, ss_vect<tps> &M_Chol, ss_vect<tps> &M_Chol_t)
{
  int k;

  const int n_dof = 3;

  long int
    lastpos;
  double
    dnu[n_dof], alpha_rad[n_dof], tau[n_dof], D[n_dof], eps[n_dof],
    sigma_s, sigma_delta;
  ss_vect<double>
    fix_point;
  ss_vect<tps>
    M_t, R, A, A_t, A_inv, A_t_inv, M_diff;

  globval.radiation = true;
  globval.Cavity_on = true;

  getcod(0e0, lastpos);
  fix_point = globval.CODvect;
  printf("\nFix point:");
  for (k = 0; k < 6; k++)
    if (k != ct_)
      printf(" %10.3e", globval.CODvect[k]);
    else
      printf(" %10.3e", globval.CODvect[k]/c0);
  printf("\n");

  M = compute_M(fix_point);
  M_t = tp_map(n_dof, M);
  prt_map(n_dof, "\nM:", M);

  A = compute_A(M, alpha_rad, R);
  compute_tau(alpha_rad, tau);
  A = get_A_CS(n_dof, A, dnu);
  A_t = tp_map(n_dof, A);
  A_inv = Inv(A);
  A_t_inv = Inv(A_t);

  prt_map(n_dof, "\nA:", A);
  printf("\ntau [msec]          = [%5.3f, %5.3f, %5.3f]\n",
	 1e3*tau[X_], 1e3*tau[Y_], 1e3*tau[Z_]);

  compute_D(fix_point, A, D);
  compute_eps(tau, D, eps);
  compute_bunch_size(A, eps, sigma_s, sigma_delta);
  M_diff = compute_M_diff(A, D);

  printf("D                   = [%9.3e, %9.3e, %9.3e]\n",
	 D[X_], D[Y_], D[Z_]);
  printf("eps                 = [%9.3e, %9.3e, %9.3e]\n",
	 eps[X_], eps[Y_], eps[Z_]);
  printf("sigma_x [micro m]   = %5.3f\n", 1e6*sqrt(eps[X_]*Cell[0].Beta[X_]));
  printf("sigma_x'[micro rad] = %5.3f\n",
	 1e6*sqrt(eps[X_]*(1e0+sqr(Cell[0].Alpha[X_]))/Cell[0].Beta[X_]));
  printf("sigma_t [ps]        = %5.3f\n", 1e12*sigma_s/c0);
  printf("sigma_delta         = %9.3e\n", sigma_delta);
  prt_map(n_dof, "\nDiffusion Matrix:", M_diff);

  M_Chol_t = compute_M_Chol_t(M_diff);
  M_Chol = tp_map(n_dof, M_Chol_t);
  prt_map(n_dof, "\nCholesky Decomposition:", M_Chol);
}


void compute_stochastic_part
(ss_vect<double> &X, ss_vect<tps> &X_map)
{
  // Compute stochastic part of map.
  static std::default_random_engine rand;
  std::normal_distribution<double>  norm_ranf(0e0, 1e0);

  const int n_dof = 3;

  for (int k = 0; k < 2*n_dof; k++) {
    X[k] = norm_ranf(rand);
    X_map[k] = X[k]*tps(0e0, k+1);
  }
}


void propagate_rad
(tps &sigma, const ss_vect<tps> &M_Chol, const ss_vect<tps> &M_Chol_t)
{
  ss_vect<double> X;
  ss_vect<tps>    X_map;

  std::default_random_engine       rand;
  std::normal_distribution<double> norm_ranf(0e0, 1e0);

  compute_stochastic_part(X, X_map);
  // sigma += (M_Chol_t*X).cst() + M_Chol_t*sqr(X_map)*M_Chol;
}


void set_HOM_long
(const string &name, const double beta, const double f, const double R_sh,
 const double Q)
{
  const long int loc = Elem_GetPos(ElemIndex(name.c_str()), 1);

  Cell[loc].Elem.C->beta_long.push_back(beta);
  Cell[loc].Elem.C->HOM_f_long.push_back(f);
  Cell[loc].Elem.C->HOM_R_sh_long.push_back(R_sh);
  Cell[loc].Elem.C->HOM_Q_long.push_back(Q);
  Cell[loc].Elem.C->HOM_V_long.push_back(0e0);
}


void prt_HOM(const int n, const string &name)
{
  const long int
    loc = Elem_GetPos(ElemIndex(name.c_str()), 1);

  printf("%3d", n);
  printf("  %22.15e %22.15e %22.15e\n",
	 abs(Cell[loc].Elem.C->HOM_V_long[0]),
	 arg(Cell[loc].Elem.C->HOM_V_long[0]),
	 real(Cell[loc].Elem.C->HOM_V_long[0]));
}


void prt_sigma(ofstream &outf, const int n, const tps &sigma)
{
  const long int
    x[]           = {1, 0, 0, 0, 0, 0, 0},
    p_x[]         = {0, 1, 0, 0, 0, 0, 0},
    y[]           = {0, 0, 1, 0, 0, 0, 0},
    p_y[]         = {0, 0, 0, 1, 0, 0, 0},
    ct[]          = {0, 0, 0, 0, 0, 1, 0},
    delta[]       = {0, 0, 0, 0, 1, 0, 0},

    x_x[]         = {2, 0, 0, 0, 0, 0, 0},
    p_x_p_x[]     = {0, 2, 0, 0, 0, 0, 0},
    y_y[]         = {0, 0, 2, 0, 0, 0, 0},
    p_y_p_y[]     = {0, 0, 0, 2, 0, 0, 0},
    ct_ct[]       = {0, 0, 0, 0, 0, 2, 0},
    delta_delta[] = {0, 0, 0, 0, 2, 0, 0};

  outf << scientific << setprecision(3)
       << setw(3) << n
       << setw(11) << sigma[x] << setw(11) << sigma[p_x]
       << setw(11) << sigma[y] << setw(11) << sigma[p_y]
       << setw(11) << sigma[delta] << setw(11) << sigma[ct]/c0 << " |"
       << setw(10) << sqrt(sigma[p_x_p_x]) << setw(10) << sqrt(sigma[x_x])
       << setw(10) << sqrt(sigma[p_y_p_y]) << setw(10) << sqrt(sigma[y_y])
       << setw(10) << sqrt(sigma[ct_ct])
       << setw(10) << sqrt(sigma[delta_delta])/c0 << "\n";
}


void propagate_mag_lat(tps &sigma, const ss_vect<tps> &M_inv)
{
  sigma = sigma*M_inv;
}


double get_ct(tps &sigma)
{
  const long int jj[] = {0, 0, 0, 0, 0, 1, 0};

  return sigma[jj];
}


void propagate_cav_HOM_transv
(const int n, const double Q_b, tps &sigma)
{

  const string
    cav_name = "cav";
  const long int
    loc = Elem_GetPos(ElemIndex(cav_name.c_str()), 1);
  const double
    Circ        = Cell[globval.Cell_nLoc].S,
    beta_RF     = Cell[loc].Elem.C->beta_RF,
    f           = Cell[loc].Elem.C->HOM_f_long[0],
    R_sh        = Cell[loc].Elem.C->HOM_R_sh_long[0],
    Q           = Cell[loc].Elem.C->HOM_Q_long[0],
    k_loss      = 2e0*M_PI*f*R_sh/Q,
    R_sh_loaded = R_sh/(1e0+beta_RF),
    Q_loaded    = Q/(1e0+beta_RF);
  const complex<double>
    I = complex<double>(0e0, 1e0);

  // Update RF cavity HOM phasor.
  Cell[loc].Elem.C->HOM_V_long[0] *=
    exp(2e0*M_PI*f*(Circ+get_ct(sigma))/c0*(-1e0/(2e0*Q_loaded)+I));

  // First half increment of HOM phasor.
  Cell[loc].Elem.C->HOM_V_long[0] += Q_b*k_loss/2e0;

  // Propagate through wake field.
  sigma -= Q_b*real(Cell[loc].Elem.C->HOM_V_long[0])/2e0;

  // Second half increment of HOM phasor.
  Cell[loc].Elem.C->HOM_V_long[0] += Q_b*k_loss/2e0;

  if (false)
    prt_HOM(n, "cav");
}


void propagate_cav_HOM_long(const int n, const double Q_b, tps &sigma)
{
  ss_vect<tps> M;

  const string
    cav_name = "cav";
  const long int
    loc = Elem_GetPos(ElemIndex(cav_name.c_str()), 1);
  const double
    Circ       = Cell[globval.Cell_nLoc].S,
    beta_RF    = Cell[loc].Elem.C->beta_RF,
    f          = Cell[loc].Elem.C->HOM_f_long[0],
    R_sh       = Cell[loc].Elem.C->HOM_R_sh_long[0],
    Q          = Cell[loc].Elem.C->HOM_Q_long[0],
    k_loss     = 2e0*M_PI*f*R_sh/Q,
    R_shloaded = R_sh/(1e0+beta_RF),
    Q_loaded   = Q/(1e0+beta_RF);
  const complex<double>
    I = complex<double>(0e0, 1e0);

  // Update RF cavity HOM phasor.
  Cell[loc].Elem.C->HOM_V_long[0] *=
    exp(2e0*M_PI*f*(Circ+get_ct(sigma))/c0*(-1e0/(2e0*Q_loaded)+I));

  // First half increment of HOM phasor.
  Cell[loc].Elem.C->HOM_V_long[0] += Q_b*k_loss/2e0;

  // Propagate through wake field.
  sigma -= Q_b*real(Cell[loc].Elem.C->HOM_V_long[0])/2e0*tps(0e0, delta_+1);

  // Second half increment of HOM phasor.
  Cell[loc].Elem.C->HOM_V_long[0] += Q_b*k_loss/2e0;

  if (!false)
    prt_HOM(n, "cav");
}


void propagate_lat
(const int n, const double Q_b, tps &sigma, const ss_vect<tps> &M_inv,
 const ss_vect<tps> &M_Chol, const ss_vect<tps> &M_Chol_t)
{
  if (true)
    propagate_cav_HOM_long(n, Q_b, sigma);

  // Propagate through magnetic lattice.
  propagate_mag_lat(sigma, M_inv);

  if (false)
    // Radiate.
    propagate_rad(sigma, M_Chol, M_Chol_t);
}


tps compute_sigma
(const double eps[], const double sigma_s, const double sigma_delta,
 const ss_vect<tps> &A)
{
  int          k;
  tps          sigma;
  ss_vect<tps> Id;

  Id.identity();
  sigma = 0e0;
  for (k = 0; k < 2; k++)
    sigma += eps[k]*(sqr(Id[2*k])+sqr(Id[2*k+1]));
  sigma = sigma*Inv(A);
  sigma += sqr(sigma_delta*Id[ct_]) + sqr(sigma_s*Id[delta_]);
  return sigma;
}


void test_case_2(const string &name)
{
  // Single bunch, 1 longitudinal HOM.
  long int        lastpos;
  int             k;
  tps             sigma;
  ss_vect<tps>    M, M_inv, A, M_Chol, M_Chol_t;
  ofstream        outf;

  const bool
    rad         = false;

  const int
    n_turn      = 9;

  const double
    eps[]       = {161.7e-12, 8e-12},
    sigma_s     = 3.739e-3,
    sigma_delta = 9.353e-04,

    Q_b         = -0.6e-9,
    beta_HOM    = 1e0,
    f           = 1e9,
    R_sh        = 1e3,
    Q           = 1e8;

  const string
    file_name = "moments.out";

  file_wr(outf, file_name.c_str());

  globval.Cavity_on = globval.radiation = rad;
  globval.pathlength = false;

  if (!rad) {
    danot_(1);

    M.identity();
    Cell_Pass(0, globval.Cell_nLoc, M, lastpos);
    M_inv = Inv(M);
    MNF = MapNorm(M, 1);

    danot_(NO);
  } else {
    compute_maps(M, M_Chol, M_Chol_t);
  }

  set_HOM_long("cav", beta_HOM, f, R_sh, Q);

  sigma = compute_sigma(eps, sigma_s, sigma_delta, MNF.A1);

  prt_sigma(outf, 0, sigma);
  printf("\n");
  for (k = 1; k <= n_turn; k++) {
    propagate_lat(k, Q_b, sigma, M_inv, M_Chol, M_Chol_t);
    prt_sigma(outf, k, sigma);
  }
}


void set_state(void)
{
  globval.H_exact        = false;
  globval.quad_fringe    = false;
  globval.Cavity_on      = false;
  globval.radiation      = false;
  globval.emittance      = false;
  globval.IBS            = false;
  globval.pathlength     = false;
  globval.Aperture_on    = false;
  globval.Cart_Bend      = false;
  globval.dip_edge_fudge = true;
}


int main(int argc, char *argv[])
{

  reverse_elem     = true;
  globval.mat_meth = false;

  // disable from TPSALib- and LieLib log messages
  idprset_(-1);

  if (false)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  set_state();

  if (false) {
    danot_(1);

    Ring_GetTwiss(true, 0e0);
    printglob();

    if (!false)
      GetEmittance(ElemIndex("cav"), true);

    prt_lat("linlat1.out", globval.bpm, true);
    prt_lat("linlat.out", globval.bpm, true, 10);
  }

  if (!false)
    test_case_2("cav");
}
