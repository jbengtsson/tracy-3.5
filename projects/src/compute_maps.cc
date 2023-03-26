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


ss_vect<double> compute_eta(const ss_vect<tps> &M)
{
  long int        jj[ss_dim];
  int             k;
  ss_vect<double> eta;
  ss_vect<tps>    Id, M_x;

  Id.identity();
  for (k = 0; k < ss_dim; k++)
    jj[k] = 0;
  // Only use 2x2 part of hor. matrix.
  M_x.zero();
  for (k = 0; k < 2; k++) {
    M_x[k] = M[k][x_]*Id[x_] + M[k][px_]*Id[px_];
    jj[k] = 1;
  }
  eta.zero();
  for (k = 0; k < 2; k++)
    eta[k] = M[k][delta_];
  eta = (PInv(Id-M_x, jj)*eta).cst();
  return eta;
}


ss_vect<tps> compute_A0(const ss_vect<double> &eta)
{
  // Canonical Transfomation to delta dependent Fixed Point.
  int          k;
  ss_vect<tps> Id, A0;

  Id.identity(); A0.identity();
  for (k = 0; k < 2; k++)
    A0[k] += eta[k]*Id[delta_];
  // Symplectic flow.
  A0[ct_] += eta[px_]*Id[x_] - eta[x_]*Id[px_];
  return A0;
}


ss_vect<tps> compute_A
(const double C, const ss_vect<tps> &M, double alpha_rad[], const bool rad)
{
  // Poincaré map Diagonalization.
  // The damping coeffients alpha[] are obtained from the eigen values.
  double       nu_s, alpha_c;
  Matrix       M_mat, A_mat, A_inv_mat, R_mat;
  ss_vect<tps> A;

  const int n_dof = (rad)? 3 : 2;

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
  return A;
}


ss_vect<tps> compute_M_cav(const string &cav_name, double &phi0)
{
  ss_vect<tps> M_cav;

  const int
    loc = Elem_GetPos(ElemIndex(cav_name.c_str()), 1);
  const double
    V_RF      = Cell[loc].Elem.C->V_RF,
    f_RF      = Cell[loc].Elem.C->f_RF,
    E0        = 1e9*globval.Energy,
    U0        = globval.dE*E0,
    delta_cav = U0/E0;

  phi0 = asin(delta_cav*E0/V_RF);
  M_cav.identity();
  M_cav[delta_] -=
    V_RF*2e0*M_PI*f_RF*cos(phi0)/(1e9*globval.Energy*c0)*tps(0e0, ct_+1);
  return M_cav;
}


void compute_tau
(const double C, const double alpha_rad[], double tau[])
{
  // Compute the damping times.
  
  for (int k = 0; k < nd_tps; k++)
    tau[k] = -C/(c0*alpha_rad[k]);
}


ss_vect<tps> compute_M_delta(const double alpha_rad[])
{
  int          k;
  double       delta_rad;
  ss_vect<tps> M_delta;

  M_delta.identity();
  for (k = 0; k < nd_tps; k++) {
    delta_rad = exp(alpha_rad[Z_]) - 1e0;
    M_delta[2*k] *= 1e0 + delta_rad;
    M_delta[2*k+1] /= 1e0 + delta_rad;
  }
  return M_delta;
}


ss_vect<tps> compute_M_tau(const double alpha[])
{
  int          k;
  ss_vect<tps> M_tau;
 
  M_tau.zero();
  for (k = 0; k < nd_tps; k++) {
    M_tau[2*k] = exp(alpha[k])*tps(0e0, 2*k+1);
    M_tau[2*k+1] = exp(alpha[k])*tps(0e0, 2*k+2);
  }
  return M_tau;
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

void compute_eps
(const double C, const double tau[], const double D[], double eps[])
{
  // Compute the eigen emittances.

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
  int    k;
  double U0;

  const string
    file_name = "compute_maps";
  const double
    C         = Cell[globval.Cell_nLoc].S,
    E0        = 1e9*globval.Energy;

  long int
    lastpos;
  double
    dnu[3], alpha_rad[3], tau[3], phi0, D[3], eps[3],
    sigma_s, sigma_delta;
  ss_vect<double>
    fix_point, eta;
  ss_vect<tps>
    M, A, M_delta, M_tau, M_cav, D_mat, M_Chol, A0;

  globval.radiation = globval.Cavity_on = true;

  getcod(0e0, lastpos);
  fix_point = globval.CODvect;
  U0 = globval.dE*E0;

  M = compute_M(fix_point);

  A = compute_A(C, M, alpha_rad, true);
  A = get_A_CS(nd_tps, A, dnu);

  compute_tau(C, alpha_rad, tau);
  M_delta = compute_M_delta(alpha_rad);

  M_tau = compute_M_tau(alpha_rad);
  M_tau = A*M_tau*Inv(A);

  M_cav = compute_M_cav("cav", phi0);

  compute_D(fix_point, A, D);
  compute_eps(C, tau, D, eps);
  compute_bunch_size(A, eps, sigma_s, sigma_delta);
  D_mat = compute_D_mat(A, D);

  M_Chol = compute_M_Chol(D_mat);

  prt_map(nd_tps, "\nM:", M);
  prt_map(file_name+"_M.dat", M);
  prt_map(nd_tps, "\nA:", A);
  prt_map(file_name+"_A.dat", A);
  prt_map(nd_tps, "\nM_delta:", M_delta);
  prt_map(file_name+"_M_delta.dat", M_delta);
  prt_map(nd_tps, "\nM_tau:", M_tau);
  prt_map(file_name+"_M_tau.dat", M_tau);
  prt_map(nd_tps, "\nDiffusion Matrix:", D_mat);
  prt_map(nd_tps, "\nCholesky Decomposition:", M_Chol);
  prt_map(file_name+"_M_Chol.dat", M_Chol);

  M = Inv(M_delta)*Inv(M_cav)*Inv(M_tau)*M;
  prt_map(nd_tps, "\nM_delta^-1.M_cav^-1.M_tau^-1.M:", M);
  eta = compute_eta(M);
  A0 = compute_A0(eta);
  prt_map(nd_tps, "\nA0:", A0);
  M = Inv(A0)*M*A0;
  prt_map(nd_tps, "\nA0^-1.M.A0:", M);

  prt_map(nd_tps, "\nM^t.Omega.M-Omega:", tp_map(3, M)*get_S(3)*M-get_S(3));

  printf("\nFix point           =");
  for (k = 0; k < 6; k++)
    if (k != ct_)
      printf(" %10.3e", globval.CODvect[k]);
    else
      printf(" %10.3e", globval.CODvect[k]/c0);
  printf("\n");
  cout << scientific << setprecision(3)
       << "eta                 ="  << setw(11) << eta << "\n";
  printf("U0 [keV]            = %3.1f\n", 1e-3*U0);
  printf("phi0 [deg]          = 180 - %4.2f\n", fabs(phi0*180e0/M_PI));
  printf("tau [msec]          = [%5.3f, %5.3f, %5.3f]\n",
	 1e3*tau[X_], 1e3*tau[Y_], 1e3*tau[Z_]);
  printf("D                   = [%9.3e, %9.3e, %9.3e]\n",
	 D[X_], D[Y_], D[Z_]);
  printf("eps                 = [%9.3e, %9.3e, %9.3e]\n",
	 eps[X_], eps[Y_], eps[Z_]);
  printf("sigma_x [micro m]   = %5.3f\n", 1e6*sqrt(eps[X_]*Cell[0].Beta[X_]));
  printf("sigma_x'[micro rad] = %5.3f\n",
	 1e6*sqrt(eps[X_]*(1e0+sqr(Cell[0].Alpha[X_]))/Cell[0].Beta[X_]));
  printf("sigma_t [ps]        = %5.3f\n", 1e12*sigma_s/c0);
  printf("sigma_delta         = %9.3e\n", sigma_delta);
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
