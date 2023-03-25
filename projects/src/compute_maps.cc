#define NO 1

#include "tracy_lib.h"

int no_tps = NO;

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


ss_vect<tps> compute_M(void)
{
  // Compute the Poincaré map.
  long int     lastpos;
  ss_vect<tps> M;

  M.identity();
  Cell_Pass(0, globval.Cell_nLoc, M, lastpos);
  return M;
}


ss_vect<tps> compute_M(const ss_vect<double> fix_point)
{
  // Compute the Poincaré map for the fix point.
  long int     lastpos;
  ss_vect<tps> M;

  M.identity();
  M += fix_point;
  Cell_Pass(0, globval.Cell_nLoc, M, lastpos);
  M -= fix_point;
  return M;
}


ss_vect<tps> compute_A
(const ss_vect<tps> &M, double alpha_rad[], ss_vect<tps> &R, const bool rad)
{
  // Poincaré map Diagonalization.
  // The damping coeffients alpha[] are obtained from the eigen values.
  double       nu_s, alpha_c;
  Matrix       M_mat, A_mat, A_inv_mat, R_mat;
  ss_vect<tps> A;

  const int    n_dof = (rad)? 3 : 2;
  const double C     = Cell[globval.Cell_nLoc].S;

  globval.radiation = globval.Cavity_on = rad;

  getlinmat(2*n_dof, M, M_mat);
  GDiag(2*n_dof, C, A_mat, A_inv_mat, R_mat, M_mat, nu_s, alpha_c);
  for (int k = 0; k < n_dof; k++)
    alpha_rad[k] = globval.alpha_rad[k];

  A = putlinmat(2*n_dof, A_mat);
  if (n_dof != 3) {
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
  
  const double C = Cell[globval.Cell_nLoc].S;

  for (int k = 0; k < nd_tps; k++)
    tau[k] = -C/(c0*alpha_rad[k]);
}


void compute_D
(const ss_vect<double> &fix_point, const ss_vect<tps> &A, double D[])
{
  // Compute the diffusion coefficients.

  int long     lastpos;
  ss_vect<tps> As, D_diag;

  globval.Cavity_on = globval.radiation = globval.emittance = true;

  As = A + fix_point;
  Cell_Pass(0, globval.Cell_nLoc, As, lastpos);

  globval.emittance = false;

  for (int k = 0; k < nd_tps; k++)
    D[k] = globval.D_rad[k];
}

void compute_eps(const double tau[], const double D[], double eps[])
{
  // Compute the eigen emittances.

  const double C = Cell[globval.Cell_nLoc].S;

  for (int k = 0; k < nd_tps; k++)
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


ss_vect<tps> compute_D_mat(const ss_vect<tps> &A, const double D[])
{
  // Compute the diffusion matrix.
  ss_vect<tps> D_diag;

  D_diag.zero();
  for (int k = 0; k < 2*nd_tps; k++)
    D_diag[k] = D[k/2]*tps(0e0, k+1);
  return A*D_diag*tp_map(nd_tps, A);
}


ss_vect<tps> compute_M_Chol(const ss_vect<tps> &M_diff)
{
  // Compute the Cholesky decomposition: D = L^T L.
  int          j, k, j1, k1;
  double       *diag, **d1, **d2;
  ss_vect<tps> M_Chol_t;

  const int n = 2*nd_tps;

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

  return tp_map(3, M_Chol_t);
}


void prt_map(const string file_name, const ss_vect<tps> map)
{
  ofstream outf;

  file_wr(outf, file_name.c_str());
  outf << scientific << setprecision(3) << setw(11) << map;
  outf.close();
}


void compute_maps(void)
{
  int k;

  const string file_name = "compute_maps";

  long int
    lastpos;
  double
    dnu[3], alpha_rad[3], tau[3], D[3], eps[3],
    sigma_s, sigma_delta;
  ss_vect<double>
    fix_point;
  ss_vect<tps>
    M, R, A, M_rad, A_rad, D_mat, M_Chol;

  globval.radiation = globval.Cavity_on = false;

  M = compute_M();
  prt_map(nd_tps, "\nM:", M);
  prt_map(file_name+"_M.dat", M);

  A = compute_A(M, alpha_rad, R, false);
  A = get_A_CS(nd_tps, A, dnu);
  prt_map(nd_tps, "\nA:", A);
  prt_map(file_name+"_A.dat", A);

  globval.radiation = globval.Cavity_on = true;

  getcod(0e0, lastpos);
  fix_point = globval.CODvect;
  printf("\nFix point:");
  for (k = 0; k < 6; k++)
    if (k != ct_)
      printf(" %10.3e", globval.CODvect[k]);
    else
      printf(" %10.3e", globval.CODvect[k]/c0);
  printf("\n");

  M_rad = compute_M(fix_point);
  prt_map(nd_tps, "\nM_rad:", M_rad);
  prt_map(file_name+"_M_rad.dat", M_rad);

  A_rad = compute_A(M_rad, alpha_rad, R, true);
  compute_tau(alpha_rad, tau);
  A_rad = get_A_CS(nd_tps, A_rad, dnu);

  prt_map(nd_tps, "\nA_rad:", A_rad);
  prt_map(file_name+"_A_rad.dat", A_rad);
  printf("\ntau [msec]          = [%5.3f, %5.3f, %5.3f]\n",
	 1e3*tau[X_], 1e3*tau[Y_], 1e3*tau[Z_]);

  compute_D(fix_point, A_rad, D);
  compute_eps(tau, D, eps);
  compute_bunch_size(A_rad, eps, sigma_s, sigma_delta);
  D_mat = compute_D_mat(A_rad, D);

  printf("D                   = [%9.3e, %9.3e, %9.3e]\n",
	 D[X_], D[Y_], D[Z_]);
  printf("eps                 = [%9.3e, %9.3e, %9.3e]\n",
	 eps[X_], eps[Y_], eps[Z_]);
  printf("sigma_x [micro m]   = %5.3f\n", 1e6*sqrt(eps[X_]*Cell[0].Beta[X_]));
  printf("sigma_x'[micro rad] = %5.3f\n",
	 1e6*sqrt(eps[X_]*(1e0+sqr(Cell[0].Alpha[X_]))/Cell[0].Beta[X_]));
  printf("sigma_t [ps]        = %5.3f\n", 1e12*sigma_s/c0);
  printf("sigma_delta         = %9.3e\n", sigma_delta);
  prt_map(nd_tps, "\nDiffusion Matrix:", D_mat);

  M_Chol = compute_M_Chol(D_mat);
  prt_map(nd_tps, "\nCholesky Decomposition:", M_Chol);
  prt_map(file_name+"_M_Chol.dat", M_Chol);
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

  daeps_(1e-15);
  compute_maps();
}
