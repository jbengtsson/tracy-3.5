
#include <random>


#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


class SigmaType {
private:
  int
    n_dof;
  double
    C,                 // Circumference.
    tau[3],            // Damping times.
    D[3],              // Diffusion coefficients.
    eps[3],            // Emittance.
    sigma_s,           // Bunch length.
    sigma_delta;       // Momentum spread.
  ss_vect<double>
    fixed_point,       // Closed orbit.
    mean;              // 1st moment.
  ss_vect<tps>
    sigma,             // 2nd moment; beam envelope.
    M,                 // Poincaré map.
    M_t,               // M^T.
    A,                 // M = A R A^-1.
    A_t,               // A^T.
    A_inv,             // A^-1.
    A_t_inv,           // (A^T)^-1.
    R,                 // Phase-space rotation.
    M_diff,            // Diffusion matrix.
    M_Chol,            // Cholesky decomposition.
    M_Chol_t;          // M_Chol^T.
  FILE
    *outf;             // Output file.
public:
  void set_params(const int n_dof, const double C, const string &file_name);
  ss_vect<double> get_eps(void) const;
  void prt_eps(const int n) const;
  void prt_sigma(void);
  void get_M(void);
  void get_A_and_tau(void);
  void get_D_and_eps(void);
  void get_sigma_s_and_delta(void);
  void get_M_diff(void);
  void get_M_Chol_t(void);
  void get_maps(void);
  void init_sigma(const double eps_x, const double eps_y, const double sigma_s,
		  const double sigma_delta);
  void propagate(const int n);
};


void prt_vec(const int n_dof, const string &str, ss_vect<double> vec)
{
  int k;

  const int n_dec = 6;

  printf("%s\n", str.c_str());
  for (k = 0; k < 2*n_dof; k++)
    printf("%*.*e", n_dec+8, n_dec, vec[k]);
  printf("\n");
}


void prt_map(const int n_dof, const string &str, ss_vect<tps> map)
{
  int i, j;

  const int n_dec = 6+6;

  printf("%s\n", str.c_str());
  for (i = 0; i < 2*n_dof; i++) {
    for (j = 0; j < 2*n_dof; j++)
      printf("%*.*e", n_dec+8, n_dec, map[i][j]);
    printf("\n");
  }
}


void prt_lin_map(const int n, const string &str, const ss_vect<tps> A)
{
  printf("%s", str.c_str());
  prt_lin_map(n, A);
}


double det_map(const int n_dof, const ss_vect<tps> &A)
{
  // Compute the determinant of a matrix.
  Matrix A_mat;

  getlinmat(2*n_dof, A, A_mat);
  return DetMat(2*n_dof, A_mat);
}


ss_vect<tps> tp_map(const int n_dof, const ss_vect<tps> &A)
{
  // Compute the transpose of a matrix.
  Matrix A_mat;

  getlinmat(2*n_dof, A, A_mat);
  TpMat(2*n_dof, A_mat);
  return putlinmat(2*n_dof, A_mat);
}


void SigmaType::set_params(const int n_dof, const double C,
			   const string &file_name)
{
  this->n_dof = n_dof;
  this->C = C;
  outf = file_write(file_name.c_str());
}


ss_vect<double> SigmaType::get_eps(void) const
{
  // Compute the emittances.
  int             k;
  ss_vect<double> eps;
  ss_vect<tps>    diag;

  diag = A_inv*sigma*A_t_inv;
  for (k = 0; k < 2*n_dof; k++)
    eps[k] = diag[k][k];
  return eps;
}


void SigmaType::prt_eps(const int n) const
{
  int             k;
  ss_vect<double> eps;

  eps = get_eps();
  fprintf(outf, "%7d", n);
  for (k = 0; k < 2*n_dof; k++)
    fprintf(outf, "%13.5e", eps[k]);
  fprintf(outf, "%13.5e %13.5e\n",
	  sqrt(sigma[delta_][delta_]), sqrt(sigma[ct_][ct_]));
}


void SigmaType::prt_sigma(void)
{
  int             k;
  ss_vect<double> eps;

  eps = get_eps();
  prt_vec(n_dof, "\neps:", eps);
  prt_map(n_dof, "\nsigma:", sigma);
  printf("\nsigma_kk:\n ");
  for (k = 0; k < 2*n_dof; k++)
    printf(" %11.5e", sqrt(sigma[k][k]));
  printf("\n");
}


void SigmaType::get_M(void)
{
  // Compute the Poincaré map for the fixed point.
  long int lastpos;

  M.identity();
  M += fixed_point;
  Cell_Pass(0, globval.Cell_nLoc, M, lastpos);
  M -= fixed_point;
}


void SigmaType::get_A_and_tau(void)
{
  // Compute the Floquet to phase space transformation and damping times.
  // The damping coeffients are obtained from the eigen values.
  int    k;
  double nu_s, alpha_c;
  Matrix M_mat, A_mat, A_inv_mat, R_mat;

  globval.Cavity_on = true;
  globval.radiation = true;

  getlinmat(2*n_dof, M, M_mat);
  GDiag(2*n_dof, C, A_mat, A_inv_mat, R_mat, M_mat, nu_s, alpha_c);

  for (k = 0; k < n_dof; k++)
    tau[k] = -C/(c0*globval.alpha_rad[k]);

  A = putlinmat(2*n_dof, A_mat);
  if (n_dof < 3) {
    // Coasting longitudinal plane.
    A[ct_] = tps(0e0, ct_+1);
    A[delta_] = tps(0e0, delta_+1);
  }
  R = putlinmat(2*n_dof, R_mat);
}


void SigmaType::get_D_and_eps(void)
{
  // Compute the diffusion coefficients and emittances.
  int long     lastpos;
  int          k;
  ss_vect<tps> As, D_diag;

  globval.Cavity_on = true;
  globval.radiation = true;
  globval.emittance = true;

  As = A + fixed_point;
  Cell_Pass(0, globval.Cell_nLoc, As, lastpos);

  globval.emittance = false;

  for (k = 0; k < n_dof; k++) {
    D[k] = globval.D_rad[k];
    eps[k] = D[k]*tau[k]*c0/(2e0*C);
  }
}

void SigmaType::get_sigma_s_and_delta(void)
{
  // Compute the bunch size: sigma_s & sigma_delta.
  double alpha_s, beta_s, gamma_s;

  alpha_s = -A[ct_][ct_]*A[delta_][ct_] - A[ct_][delta_]*A[delta_][delta_];
  beta_s = sqr(A[ct_][ct_]) + sqr(A[ct_][delta_]);
  gamma_s = (1e0+sqr(alpha_s))/beta_s;

  // Bunch size.
  sigma_s = sqrt(beta_s*eps[Z_]);
  sigma_delta = sqrt(gamma_s*eps[Z_]);
}


void SigmaType::get_M_diff(void)
{
  // Compute the diffusion matrix.
  int          k;
  ss_vect<tps> As, D_diag;

  D_diag.zero();
  for (k = 0; k < 2*n_dof; k++)
    D_diag[k] = D[k/2]*tps(0e0, k+1);
  M_diff = A*D_diag*tp_map(n_dof, A);
}


void SigmaType::get_M_Chol_t(void)
{
  // Compute the Cholesky decomposition: D = L^T L.
  int          j, k, j1, k1;
  double       *diag, **d1, **d2;

  const int n = 2*n_dof;

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
}


void SigmaType::get_maps(void)
{
  long int lastpos;
  double   dnu[3];

  globval.radiation = true;
  globval.Cavity_on = true;

  getcod(0e0, lastpos);
  fixed_point = globval.CODvect;
  prt_vec(n_dof, "\nFixed Point:", fixed_point);

  get_M();
  M_t = tp_map(n_dof, M);
  prt_map(n_dof, "\nM:", M);

  get_A_and_tau();
  A = get_A_CS(n_dof, A, dnu);
  A_t = tp_map(n_dof, A);
  A_inv = Inv(A);
  A_t_inv = Inv(A_t);

  prt_map(n_dof, "\nA:", A);
  printf("\ntau [msec]   = [%5.3f, %5.3f, %5.3f]\n",
	 1e3*tau[X_], 1e3*tau[Y_], 1e3*tau[Z_]);

  get_D_and_eps();
  get_sigma_s_and_delta();
  get_M_diff();

  printf("D            = [%9.3e, %9.3e, %9.3e]\n", D[X_], D[Y_], D[Z_]);
  printf("eps          = [%9.3e, %9.3e, %9.3e]\n", eps[X_], eps[Y_], eps[Z_]);
  printf("sigma_s [mm] = %5.3f\n", 1e3*sigma_s);
  printf("sigma_delta  = %9.3e\n", sigma_delta);
  prt_map(n_dof, "\nDiffusion Matrix:", M_diff);

  get_M_Chol_t();
  M_Chol = tp_map(n_dof, M_Chol_t);
  prt_map(n_dof, "\nCholesky Decomposition:", M_Chol);
}


void SigmaType::init_sigma(const double eps_x, const double eps_y,
			   const double sigma_s, const double sigma_delta)
{
  int          k;
  ss_vect<tps> Id;

  const double eps0[] = {eps_x, eps_y};

  Id.identity();
  sigma.zero();
  for (k = 0; k < 4; k++)
    sigma[k] += eps0[k/2]*Id[k];
  sigma = A*sigma*A_t;
  sigma[ct_] += sqr(sigma_s)*Id[ct_];
  sigma[delta_] += sqr(sigma_delta)*Id[delta_];
}


void SigmaType::propagate(const int n)
{
  int             j, k, n_rnd;
  double          rnd;
  ss_vect<double> sum, sum2, m, s, X;
  ss_vect<tps>    Id, X_map;

  std::default_random_engine       rand;
  std::normal_distribution<double> norm_ranf(0e0, 1e0);

  Id.identity();
  n_rnd = 0;
  for (j = 1; j <= n; j++) {
    n_rnd++;
    for (k = 0; k < 2*n_dof; k++) {
      rnd = norm_ranf(rand);
      sum[k] += rnd;
      sum2[k] += sqr(rnd);
      X[k] = rnd;
      X_map[k] = X[k]*Id[k];
    }
    // Average of stochastic term is zero.
    mean = (M*mean).cst() + 0*(M_Chol_t*X).cst();
    sigma = M*tp_map(n_dof, M*sigma) + M_Chol_t*sqr(X_map)*M_Chol;

    prt_eps(j);
  }
}


void track(void)
{
  SigmaType s;

  const string
    file_name = "wake_field.out";
  const int
    n_dof       = 3,
    n           = 30000;
  const double
    eps0[]      = {0e-9, 0.2e-9, 0e-3},
    sigma_s     = 0e0,
    sigma_delta = 1e-3;

  s.set_params(n_dof, Cell[globval.Cell_nLoc].S, file_name);
  s.get_maps();
  s.init_sigma(eps0[X_], eps0[Y_], sigma_s, sigma_delta);
  s.propagate(n);
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

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  set_state();

  Ring_GetTwiss(true, 0e0);
  printglob();

  track();
}
