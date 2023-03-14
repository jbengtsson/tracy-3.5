#define NO 1

#include <random>
#include <complex>

#include "tracy_lib.h"


int no_tps = NO;


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


void compute_C_S_long(double &alpha_z, double &beta_z)
{
  alpha_z =
    -globval.Ascr[ct_][ct_]*globval.Ascr[delta_][ct_]
    - globval.Ascr[ct_][delta_]*globval.Ascr[delta_][delta_];
  beta_z = sqr(globval.Ascr[ct_][ct_]) + sqr(globval.Ascr[ct_][delta_]);
}


ss_vect<tps> compute_A_A_t(void)
{
  int          k;
  double       alpha_z, beta_z;
  ss_vect<tps> Id, A_A_t;

  Id.identity();

  if (nd_tps > 2) compute_C_S_long(alpha_z, beta_z);

  A_A_t.zero();
  for (k = 0; k < nd_tps; k++) {
    if (k < 2) {
      A_A_t[2*k] = Cell[0].Beta[k]*Id[2*k] - Cell[0].Alpha[k]*Id[2*k+1];
      A_A_t[2*k+1] =
	-Cell[0].Alpha[k]*Id[2*k]
	+ (1e0+sqr(Cell[0].Alpha[k]))/Cell[0].Beta[k]*Id[2*k+1];
    } else {
      A_A_t[ct_] = beta_z*Id[ct_] - alpha_z*Id[delta_];
      A_A_t[delta_] = -alpha_z*Id[ct_]	+ (1e0+sqr(alpha_z))/beta_z*Id[delta_];
    }
  }

  return A_A_t;
}


tps compute_twoJ(const double eps[], const ss_vect<tps> &A_A_t)
{
  int          j, k;
  tps          twoJ;
  ss_vect<tps> Id, quad_form;

  const ss_vect<tps> omega = get_S(nd_tps);

  Id.identity();

  quad_form = tp_S(nd_tps, omega)*A_A_t*omega;
  twoJ = 0e0;
  for (j = 0; j < 2*nd_tps; j++)
    for (k = 0; k < 2*nd_tps; k++)
      twoJ += Id[j]*sqrt(eps[j/2])*quad_form[j][k]*sqrt(eps[k/2])*Id[k];
  return twoJ;
}


void compute_Twiss(const tps &twoJ, double alpha[], double beta[])
{
  long int jj[ss_dim];
  int      k;

  for (k = 0; k < ss_dim; k++)
    jj[k] = 0;
  for (k = 0; k < 2; k++) {
    jj[2*k]   = 1;
    jj[2*k+1] = 1;
    alpha[k] = twoJ[jj]/2e0;
    jj[2*k]   = 0;
    jj[2*k+1] = 0;
    jj[2*k+1] = 2;
    beta[k] = twoJ[jj];
    jj[2*k+1] = 0;
  }
}


void tst_moment(void)
{
  long int     lastpos;
  ss_vect<tps> M, A_A_t;

  const int
    nd_tps = 3,
#if 0
    loc   = globval.Cell_nLoc;
#else
    loc   = 10;
#endif

  globval.Cavity_on = !false;
  globval.radiation = !false;

  Ring_GetTwiss(true, 0e0);
  printglob();

  M.identity();
  Cell_Pass(0, loc, M, lastpos);
  printf("\nM:");
  prt_lin_map(nd_tps, M);

  A_A_t = compute_A_A_t();
  printf("\nInitial A_A_t beta = [%5.3f, %5.3f]:",
	 A_A_t[x_][x_], A_A_t[y_][y_]);
  prt_lin_map(nd_tps, A_A_t);

  A_A_t = M*tp_S(2, M*A_A_t);

  printf("\nFinal A_A_t beta = [%5.3f, %5.3f]:",
	 A_A_t[x_][x_], A_A_t[y_][y_]);
  prt_lin_map(nd_tps, A_A_t);

  A_A_t = compute_A_A_t();
  printf("\nInitial A_A_t beta = [%5.3f, %5.3f]:",
	 A_A_t[x_][x_], A_A_t[y_][y_]);
  prt_lin_map(nd_tps, A_A_t);

  Cell_Pass(0, loc, A_A_t, lastpos);
  A_A_t = tp_S(2, A_A_t);

  Cell_Pass(0, loc, A_A_t, lastpos);
  printf("\nFinal A_A_t beta = [%5.3f, %5.3f]:",
	 A_A_t[x_][x_], A_A_t[y_][y_]);
  prt_lin_map(nd_tps, A_A_t);
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


ss_vect<tps> compute_M(const ss_vect<double> fixed_point)
{
  // Compute the Poincaré map for the fixed point.
  long int lastpos;
  ss_vect<tps> M;

  M.identity();
  M += fixed_point;
  Cell_Pass(0, globval.Cell_nLoc, M, lastpos);
  M -= fixed_point;
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
(const ss_vect<double> &fixed_point, const ss_vect<tps> &A, double D[])
{
  // Compute the diffusion coefficients.

  int long     lastpos;
  ss_vect<tps> As, D_diag;

  const int n_dof = 3;

  globval.Cavity_on = globval.radiation = globval.emittance = true;

  As = A + fixed_point;
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
  const int n_dof = 3;

  long int
    lastpos;
  double
    dnu[n_dof], alpha_rad[n_dof], tau[n_dof], D[n_dof], eps[n_dof],
    sigma_s, sigma_delta;
  ss_vect<double>
    fixed_point;
  ss_vect<tps>
    M_t, R, A, A_t, A_inv, A_t_inv, M_diff;

  globval.radiation = true;
  globval.Cavity_on = true;

  getcod(0e0, lastpos);
  fixed_point = globval.CODvect;
  prt_vec(n_dof, "\nFixed Point:", fixed_point);

  M = compute_M(fixed_point);
  M_t = tp_map(n_dof, M);
  prt_map(n_dof, "\nM:", M);

  A = compute_A(M, alpha_rad, R);
  compute_tau(alpha_rad, tau);
  A = get_A_CS(n_dof, A, dnu);
  A_t = tp_map(n_dof, A);
  A_inv = Inv(A);
  A_t_inv = Inv(A_t);

  prt_map(n_dof, "\nA:", A);
  printf("\ntau [msec]   = [%5.3f, %5.3f, %5.3f]\n",
	 1e3*tau[X_], 1e3*tau[Y_], 1e3*tau[Z_]);

  compute_D(fixed_point, A, D);
  compute_eps(tau, D, eps);
  compute_bunch_size(A, eps, sigma_s, sigma_delta);
  M_diff = compute_M_diff(A, D);

  printf("D            = [%9.3e, %9.3e, %9.3e]\n", D[X_], D[Y_], D[Z_]);
  printf("eps          = [%9.3e, %9.3e, %9.3e]\n", eps[X_], eps[Y_], eps[Z_]);
  printf("sigma_s [mm] = %5.3f\n", 1e3*sigma_s);
  printf("sigma_delta  = %9.3e\n", sigma_delta);
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


void propagate_rad(ss_vect<double> &ps, ss_vect<tps> &M_Chol_t)
{
  ss_vect<double> X;
  ss_vect<tps>    X_map;

  std::default_random_engine       rand;
  std::normal_distribution<double> norm_ranf(0e0, 1e0);

  compute_stochastic_part(X, X_map);
  ps += (M_Chol_t*X).cst();
}


void set_HOM
(const string &name, const double f, const double Q, const double R_s)
{
  const long int loc = Elem_GetPos(ElemIndex(name.c_str()), 1);

  Cell[loc].Elem.C->HOM_f_long.push_back(f);
  Cell[loc].Elem.C->HOM_Q_long.push_back(Q);
  Cell[loc].Elem.C->HOM_Z_long.push_back(R_s);
  Cell[loc].Elem.C->HOM_V_long.push_back(0e0);
}


void prt_HOM(const string &name)
{
  const long int
    loc = Elem_GetPos(ElemIndex(name.c_str()), 1);

  printf("  U = %22.15e %22.15e %22.15e\n",
	 abs(Cell[loc].Elem.C->HOM_V_long[0]),
	 arg(Cell[loc].Elem.C->HOM_V_long[0]),
	 real(Cell[loc].Elem.C->HOM_V_long[0]));
}


void prt_ps(FILE *outf, const int n, const ss_vect<double> &ps)
{
  int k;

  fprintf(outf, "%3d", n);
  for (k = 0; k < 6; k++)
    if (k < 5)
      fprintf(outf, " %22.15e", ps[k]);
    else
      fprintf(outf, " %22.15e", ps[k]/c0);
  fprintf(outf, "\n");
}


void propagate_mag_lat(ss_vect<double> &ps, ss_vect<tps> &M)
{
  ps = (M*ps).cst();
}


void propagate_lat
(FILE *outf, const int n, const double Q_b, ss_vect<double> &ps,
 ss_vect<tps> &M, ss_vect<tps> &M_Chol_t, const bool map_track)
{
  long int lastpos;

  const string cav_name = "cav";

  const long int
    loc = Elem_GetPos(ElemIndex(cav_name.c_str()), 1);

  const double
    Circ   = Cell[globval.Cell_nLoc].S,
    Brho   = globval.Energy*1e9/c0,
    f      = Cell[loc].Elem.C->HOM_f_long[0],
    Q      = Cell[loc].Elem.C->HOM_Q_long[0],
    R_s    = Cell[loc].Elem.C->HOM_Z_long[0],
    k_loss = 2e0*M_PI*f*R_s/Q;

  const complex<double>
    I = complex<double>(0e0, 1e0);

  if (! map_track)
    Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
  else {
    // Propagate through magnetic lattice.
    propagate_mag_lat(ps, M);

    if (!true)
      // Radiate.
      propagate_rad(ps, M_Chol_t);

    if (!true) {
      // Update RF cavity HOMs voltage.
      Cell[loc].Elem.C->HOM_V_long[0] *=
	exp(2e0*M_PI*f*(Circ+ps[ct_])/c0*(-1e0/Q+I));
      Cell[loc].Elem.C->HOM_V_long[0] += Q_b*k_loss;

      prt_HOM("cav");
    }

    if (!true)
      // Propagate throug wake field.
      ps[delta_] -= Q_b*c0*real(Cell[loc].Elem.C->HOM_V_long[0])/(2e0*Brho);
  }
  prt_ps(outf, n, ps);
}


void test_case_2(const string &name)
{
  // Single bunch, 1 longitudinal HOM.
  long int        lastpos;
  int             k;
  ss_vect<double> ps;
  ss_vect<tps>    M, M_Chol, M_Chol_t;
  FILE            *outf;

  const int
    n_turn = 20000;

  const double
    Q_b    = -0.6e-9,
    f      = 1e9,
    Q      = 1e8,
    R_s    = 1e3;

  const string
    file_name = "track.out";

  outf = file_write(file_name.c_str());

  globval.Cavity_on = globval.radiation = true;
  globval.pathlength = false;

  const complex<double>
    u1 = polar(3.769911184307752e-05, -3.141592643589793e+00),
    u2 = polar(6.490144782493358e-05,  2.607734888987444e+00);

  printf("\n du = %22.15e + i %22.15e\n", abs(u2-u1), arg(u2-u1));

  compute_maps(M, M_Chol, M_Chol_t);

  set_HOM("cav", f, Q, R_s);

  globval.Cavity_on = globval.radiation = true;
  getcod(0e0, lastpos);
  printf("\nFixed point:");
  for (k = 0; k < 6; k++)
    if (k != ct_)
      printf(" %10.3e", globval.CODvect[k]);
    else
      printf(" %10.3e", globval.CODvect[k]/c0);
  printf("\n");

  ps.zero();
  ps += globval.CODvect;
  prt_ps(outf, 0, ps);
  printf("\n");
  for (k = 1; k <= n_turn; k++)
    propagate_lat(outf, k, Q_b, ps, M, M_Chol_t, false);
}


int main(int argc, char *argv[])
{
  double     alpha[2], beta[2], dnu[2], eta[2], etap[2];
  MomentType m;

  globval.H_exact    = false; globval.quad_fringe    = false;
  globval.Cavity_on  = false; globval.radiation      = false;
  globval.emittance  = false; globval.IBS            = false;
  globval.pathlength = false; globval.bpm            = 0;
  globval.Cart_Bend  = false; globval.dip_edge_fudge = true;
  globval.mat_meth   = false;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  if (!false) {
    Ring_GetTwiss(true, 0e0);
    printglob();
    if (false)
      GetEmittance(ElemIndex("cav"), true);
  }

  if (!false)
    test_case_2("cav");

  if (false)
    tst_moment();

  if (false) {
    m.Sigma = putlinmat(6, globval.Ascr);
    prt_lin_map(3, get_A_CS(2, m.Sigma, dnu));
    get_ab(m.Sigma, alpha, beta, dnu, eta, etap);
    printf("\n  alpha = [%5.3f, %5.3f] beta = [%5.3f, %5.3f]\n",
	   alpha[X_], alpha[Y_], beta[X_], beta[Y_]);

    m.propagate();

    prt_lin_map(3, get_A_CS(2, m.Sigma, dnu));
    get_ab(m.Sigma, alpha, beta, dnu, eta, etap);
    printf("\n  alpha = [%5.3f, %5.3f] beta = [%5.3f, %5.3f]\n",
	   alpha[X_], alpha[Y_], beta[X_], beta[Y_]);
  }
}
