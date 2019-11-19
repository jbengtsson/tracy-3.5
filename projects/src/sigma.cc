#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


const int
  mat_dim     = 6,
  mat_vec_dim = mat_dim*(mat_dim-1)/2 + mat_dim;


ss_vect<tps> get_sympl_form(const int dof)
{
  int          k;
  ss_vect<tps> Id, omega;

  Id.identity(); omega.zero();
  for (k = 0; k < dof; k++) {
    omega[2*k] = Id[2*k+1]; omega[2*k+1] = -Id[2*k];
  }
  return omega;
}


ss_vect<tps> lin_map_tp(const int n_dim, const ss_vect<tps> &M)
{
  int          i, j;
  ss_vect<tps> Id, M_tp;

  Id.identity();
  M_tp.zero();
  for (i = 0; i < 2*n_dim; i++)
    for (j = 0; j < 2*n_dim; j++) {
      M_tp[i] += M[j][i]*Id[j];
    }
  return M_tp;
}


void prt_vec(const int n, double *a)
{
  int i;

  printf("\nvector:\n");
  for (i = 1; i <= n; i++)
    printf("%11.3e", a[i]);
  printf("\n");
}


void prt_mat(const int n, double **A)
{
  int i, j;

  printf("\nmatrix:\n");
  for (i = 1; i <= n; i++) {
    for (j = 1; j <= n; j++)
      printf("%11.3e", A[i][j]);
    printf("\n");
  }
}


void prt_mat(const std::vector< std::vector<double> > &A)
{
  int i, j, m, n;

  m = (int)A.size(); n = (int)A[0].size();
  printf("\nmatrix (%d, %d):\n", m, n);
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++)
      printf("%11.3e", A[i][j]);
    printf("\n");
  }
}


void prt_mat_vec(const std::vector<double> &A)
{
  int i, n;

  n = (int)A.size();
  printf("\nvector (%d):\n", n);
  for (i = 0; i < n; i++)
    printf("%10.3e", A[i]);
  printf("\n");
}


void get_lin_mat(const int n_dim, const ss_vect<tps> &map,
		 std::vector< std::vector<double> > &mat)
{
  int                 i, j;
  std::vector<double> row;

  for (i = 0; i < n_dim; i++) {
    row.clear();
    for (j = 0; j < n_dim; j++)
      row.push_back(map[i][j]);
    mat.push_back(row);
  }
}


void Kronecker_prod(const std::vector< std::vector<double> > &A,
		    const std::vector< std::vector<double> > &B,
		    std::vector< std::vector<double> > &C)
{
  int                 i, j, k, l, m, n;
  std::vector<double> row;

  m = (int)B.size(); n = (int)B[0].size();

  row.clear();
  for (j = 0; j < n*(int)A[0].size(); j++)
    row.push_back(0e0);
  for (i = 0; i < m*(int)A.size(); i++)
    C.push_back(row);

  for (i = 0; i < (int)A.size(); i++)
    for (j = 0; j < (int)A[0].size(); j++)
      for (k = 0; k < (int)B.size(); k++)
	for (l = 0; l < (int)B[0].size(); l++)
	  C[i*m+k][j*n+l] = A[i][j]*B[k][l];
}


void get_mat(std::vector< std::vector<double> > &A, double **B)
{
  int i, j;

  for (i = 0; i < (int)A.size(); i++)
    for (j = 0; j < (int)A[0].size(); j++)
      B[i+1][j+1] = A[i][j];
}


void vech(const ss_vect<tps> &M, double *M_vec)
{
  // Symmetric matrix vectorization. 
  int i, j, k;

  k = 0;
  for (i = 0; i < mat_dim; i++)
    for (j = i; j < mat_dim; j++) {
      k++;
      M_vec[k] = M[i][j];
    }
}


void inv_vech(const double *M_vec, ss_vect<tps> &M)
{
  // Inverse of symmetric matrix vectorization.
  int          i, j, k;
  ss_vect<tps> Id;

  Id.identity();
  M.zero();
  k = 0;
  for (i = 0; i < mat_dim; i++)
    for (j = i; j < mat_dim; j++) {
      k++;
      M[i] += M_vec[k]*Id[j];
    }
  for (i = 0; i < mat_dim; i++)
    for (j = 0; j < i; j++)
      M[i] += M[j][i]*Id[j];
}


void bubble_sort(std::vector<int> &order)
{
  bool   swapped;
  int    k, ind;

  do {
    swapped = false;
    for (k = 0; k < (int)order.size()-1; k++) {
      if (order[k] < order[k+1]) {
	ind = order[k]; order[k] = order[k+1]; order[k+1] = ind;
	swapped = true;
      }
    }
  } while (swapped);
}


void red_sym_mat_1(const int n_dim, const int i, const int j,
		   std::vector< std::vector<double> > &M)
{
  // Matrix dim reduction from symmetry.
  int k;

  for (k = 0; k < (int)M.size(); k++) {
    M[k][i-1] += M[k][j-1];
    M[k][j-1] = 0e0;
  }
}


void red_sym_mat_2(const std::vector<int> ind,
		   std::vector< std::vector<double> > &M)
{
  // Matrix dim reduction from symmetry.
  int j, k;

  for (j = 0; j < (int)ind.size(); j++) {
    M.erase(M.begin()+ind[j]-1);
    for (k = 0; k < (int)M.size(); k++)
      M[k].erase(M[k].begin()+ind[j]-1);
  }
}


void get_M_M_tp(const ss_vect<tps> &M)
{
  double **M_M_tp;

  M_M_tp = dmatrix(1, mat_vec_dim, 1, mat_vec_dim);

  M_M_tp[1][1] = sqr(M[x_][x_]);
  M_M_tp[1][2] = 2e0*M[x_][x_]*M[x_][px_];
  M_M_tp[1][3] = sqr(M[x_][px_]);

  M_M_tp[2][1] = M[x_][x_]*M[px_][x_];
  M_M_tp[2][2] = M[x_][px_]*M[px_][x_] + M[x_][x_]*M[px_][px_];
  M_M_tp[2][3] = M[x_][px_]*M[px_][px_];

  M_M_tp[3][1] = sqr(M[px_][x_]);
  M_M_tp[3][2] = 2e0*M[px_][x_]*M[px_][px_];
  M_M_tp[3][3] = sqr(M[px_][px_]);

  prt_mat(3, M_M_tp);

  free_dmatrix(M_M_tp, 1, mat_vec_dim, 1, mat_vec_dim);
}


void get_M_M_tp(const ss_vect<tps> &M, double **M_M_tp)
{
  // Roth's Relationship (W. Roth "On Direct Product Matrices" Bul. Amer. Math.
  // Soc. 40, 461-468 (1934):
  //   vec{A B C} = (C^T x_circ A) vec(B)
  int                                i, j, i1, j1;
  std::vector<int>                   ind;
  std::vector< std::vector<double> > M_mat, M_prod;

  get_lin_mat(mat_dim, M, M_mat);
  Kronecker_prod(M_mat, M_mat, M_prod);

  for (i = 1; i <= mat_dim-1; i++) {
    for (j = 1; j <= mat_dim-i; j++) {
      i1 = (i-1)*(mat_dim+1) + j + 1; j1 = i1 + j*(mat_dim-1);
      ind.push_back(j1);
      red_sym_mat_1(mat_dim, i1, j1, M_prod);
    }
  }
  bubble_sort(ind);
  red_sym_mat_2(ind, M_prod);
  get_mat(M_prod, M_M_tp);
}


void get_emit(double **M_M_tp, double *D_vec)
{
  int    i, j;
  double       *sigma_vec, **Id, **MmI, **MmI_inv;
  ss_vect<tps> sigma;

  sigma_vec = dvector(1, mat_vec_dim);
  Id = dmatrix(1, mat_vec_dim, 1, mat_vec_dim);
  MmI = dmatrix(1, mat_vec_dim, 1, mat_vec_dim);
  MmI_inv = dmatrix(1, mat_vec_dim, 1, mat_vec_dim);

  for (i = 1; i <= mat_vec_dim; i++)
    for (j = 1; j <= mat_vec_dim; j++)
      Id[i][j] = (i == j)? 1e0 : 0e0;
 
  dmsub(Id, mat_vec_dim, mat_vec_dim, M_M_tp, MmI);
  dinverse(MmI, mat_vec_dim, MmI_inv);  
  dmvmult(MmI_inv, mat_vec_dim, mat_vec_dim, D_vec, mat_vec_dim, sigma_vec);
  inv_vech(sigma_vec, sigma);

  printf("\nget_emit:\n");
  prt_lin_map(3, sigma);

  free_dvector(sigma_vec, 1, mat_vec_dim);
  free_dmatrix(Id, 1, mat_vec_dim, 1, mat_vec_dim);
  free_dmatrix(MmI, 1, mat_vec_dim, 1, mat_vec_dim);
  free_dmatrix(MmI_inv, 1, mat_vec_dim, 1, mat_vec_dim);
}


void prt_rad(const double E_0, const double rho)
{
  int    k;
  double C, gamma, U_0, P_gamma, P_gamma_avg, u_c, sigma_delta, eps_x, D[3];
  double gamma_z, J[3], tau[3], N, N_u2_avg, Q_p_t, T_0, I[6];

  C = Cell[globval.Cell_nLoc].S;
  T_0 = C/c0;

  gamma = 1e9*E_0/m_e;
  u_c = 3e0*h_bar*c0*cube(gamma)/(2e0*rho);
  P_gamma = 1e9*C_gamma*c0*pow(E_0, 4)/(2e0*M_PI*sqr(rho));
  N = 15e0*sqrt(3e0)*P_gamma/(8e0*u_c);

  printf("\n  m_e                  = %11.5e\n", m_e);
  printf("  r_e                  = %11.5e\n", r_e);
  printf("  E [GeV]              = %11.5e gamma = %11.5e rho = %11.5e\n",
	 E_0, gamma, rho);
  printf("  C_gamma [m/GeV^3]    = %11.5e\n", C_gamma);
  printf("  C_u []               = %11.5e\n", C_u);
  printf("  C_q [m]              = %11.5e\n", C_q);
  printf("  u_c [keV]            = %11.5e\n", 1e-3*u_c);
  printf("  P_gamma [GeV]        = %11.5e\n", 1e-9*P_gamma);
  printf("  N [sec^-1]           = %11.5e\n", N);
  printf("  C_q []               = %11.5e\n", 1e18*C_q/sqr(m_e));

  get_I(I, false);

  get_eps_x();

  P_gamma_avg = 1e9*C_gamma*c0*pow(E_0, 4)*I[2]/(2e0*M_PI*C);
  U_0 = 1e9*C_gamma*pow(E_0, 4)*I[2]/(2e0*M_PI);
  eps_x = C_q*sqr(gamma)*I[5]/(I[2]-I[4]);
  sigma_delta = sqrt(C_q*I[3]/(2e0*I[2]+I[4]))*gamma;
  N_u2_avg =
    1e36*3e0*C_u*C_gamma*h_bar*sqr(c0)*pow(E_0, 7)/(4e0*M_PI*cube(m_e*rho));
  Q_p_t = 3e0*C_u*h_bar*c0*cube(gamma)*P_gamma_avg*I[3]/(2e0*I[2]);

  J[X_] = 1e0 - I[4]/I[2]; J[Z_] = 2e0 + I[4]/I[2]; J[Y_] = 4e0 - J[X_] - J[Z_];
  for (k = 0; k < 3; k++)
    tau[k] = 4e0*M_PI*T_0/(C_gamma*cube(E_0)*J[k]*I[2]);
  D[X_] = C_q*C_gamma*sqr(gamma)*cube(E_0)*I[5]/(2e0*M_PI);
  D[Y_] = 0e0;
  D[Z_] = C_q*C_gamma*sqr(gamma)*cube(E_0)*I[3]/(2e0*M_PI);

  printf("\n  C [m]                 = %7.3f\n", C);
  printf("  E_0 [GeV]             = %3.1f\n", E_0);
  printf("  U_0 [keV]             = %5.1f\n", 1e-3*U_0);
  printf("  <P_gamma> [GeV/sec]   = %9.3e\n", 1e-9*c0*U_0/C);
  printf("  <P_gamma> [GeV/sec]   = %9.3e\n", 1e-9*P_gamma_avg);
  printf("  eps_x [m.rad]         = %9.3e\n", eps_x);
  printf("  sigma_delta           = %9.3e\n", sigma_delta);
  printf("  J                     = [%5.3f, %5.3f, %5.3f]\n",
	 J[X_], J[Y_], J[Z_]);
  printf("  tau [msec]            = [%5.3f, %5.3f, %5.3f]\n",
	 1e3*tau[X_], 1e3*tau[Y_], 1e3*tau[Z_]);
  printf("  N<u^2> [eV^2 sec^-1]  = %9.3e\n", C_u*u_c*P_gamma);
  printf("  N<u^2> [eV^2 sec^-1]  = %9.3e\n", N_u2_avg);
  printf("  Q_p_t [eV^2 sec^-1]   = %9.3e\n",
	 4e0*sqr(1e9*E_0*sigma_delta)/tau[Z_]);
  printf("  Q_p_t [eV^2 sec^-1]   = %9.3e\n", Q_p_t);
  printf("  [D_J_x, D_J_y, D_p_t] = [%9.3e, %9.3e, %9.3e]\n",
	 D[X_], D[Y_], D[Z_]);

  globval.Cavity_on = true; globval.radiation = true;
  Ring_GetTwiss(true, 0e0);
  globval.alpha_z =
    -globval.Ascr[ct_][ct_]*globval.Ascr[delta_][ct_]
    - globval.Ascr[ct_][delta_]*globval.Ascr[delta_][delta_];
  globval.beta_z = sqr(globval.Ascr[ct_][ct_]) + sqr(globval.Ascr[ct_][delta_]);
  gamma_z = (1.0+sqr(globval.alpha_z))/globval.beta_z;

  printf("  [D_J_x, D_J_y, D_J_s] = [%9.3e, %9.3e, %9.3e]\n",
	 D[X_], D[Y_], D[Z_]/gamma_z);
}


void get_lin_map()
{
  int          k;
  double       C, mu[3], alpha[3], beta[3], gamma[3], tau[3];
  ss_vect<tps> Id, A, M, M_rad;

  const bool rad = true;

  Id.identity();

  C = Cell[globval.Cell_nLoc].S;

  globval.Cavity_on = true; globval.radiation = rad;
  Ring_GetTwiss(true, 0e0); printglob();
  putlinmat(6, globval.Ascr, A);

  for (k = 0; k < 2; k++) {
    mu[k] = 2e0*M_PI*globval.TotalTune[k];
    alpha[k] = Cell[0].Alpha[k]; beta[k] = Cell[0].Beta[k];
    gamma[k] = (1e0+sqr(alpha[k]))/beta[k];
  }

  mu[Z_] = -2e0*M_PI*globval.Omega;
  alpha[Z_] =
    -globval.Ascr[ct_][ct_]*globval.Ascr[delta_][ct_]
    - globval.Ascr[ct_][delta_]*globval.Ascr[delta_][delta_];
  beta[Z_] = sqr(globval.Ascr[ct_][ct_]) + sqr(globval.Ascr[ct_][delta_]);
  gamma[Z_] = (1e0+sqr(alpha[Z_]))/beta[Z_];

  M.zero();
  for (k = 0; k < 3; k++) {
    if (k < 2) {
      M[2*k] =
	(cos(mu[k])+alpha[k]*sin(mu[k]))*Id[2*k] + beta[k]*sin(mu[k])*Id[2*k+1];
      M[2*k+1] =
	-gamma[k]*sin(mu[k])*Id[2*k]
	+ (cos(mu[k])-alpha[k]*sin(mu[k]))*Id[2*k+1];
    } else {
      M[2*k+1] =
	(cos(mu[k])+alpha[k]*sin(mu[k]))*Id[2*k+1] + beta[k]*sin(mu[k])*Id[2*k];
      M[2*k] =
	-gamma[k]*sin(mu[k])*Id[2*k+1]
	+ (cos(mu[k])-alpha[k]*sin(mu[k]))*Id[2*k];
    }
  }
  M[x_] += globval.OneTurnMat[x_][delta_]* Id[delta_];
  M[px_] += globval.OneTurnMat[px_][delta_]*Id[delta_];
  M[ct_] +=
    (M[x_][x_]*M[px_][delta_]-M[px_][x_]*M[x_][delta_])*Id[x_]
    +(M[x_][px_]*M[px_][delta_]-M[px_][px_]*M[x_][delta_])*Id[px_]
    - M[x_][delta_]*M[px_][delta_]*Id[delta_];
  M[delta_] +=
    -sqr(mu[Z_])*M[ct_][x_]/M[ct_][delta_]*Id[x_]
    -sqr(mu[Z_])*M[ct_][px_]/M[ct_][delta_]*Id[px_];

  if (rad) {
    M_rad.zero();
    for (k = 0; k < 3; k++) {
      tau[k] = -C/(c0*globval.alpha_rad[k]);
      M_rad[2*k] = exp(-C/(c0*tau[k]))*Id[2*k];
      M_rad[2*k+1] = exp(-C/(c0*tau[k]))*Id[2*k+1];
    }
    M = M_rad*M;

    printf("\n tau [msec]: %6.3f %6.3f %6.3f\n",
	   1e3*tau[X_], 1e3*tau[Y_], 1e3*tau[Z_]);
  }

  prt_lin_map(3, M);
}


void get_sigma(void)
{
  int          k;
  double       C, *sigma_vec, *sigma_vec1, *D_vec, **M_M_tp;
  ss_vect<tps> Id, J, A, D, sigma, sigma1, M, M_tp, A_A_tp;

  const int    n_turn = 500000;
  const double
    L   = 2.62,
    phi = 6.0,
    rho = L/(phi*M_PI/180.0);

  printf("\nget_sigma: n = %d\n", mat_vec_dim);

  sigma_vec = dvector(1, mat_vec_dim);
  sigma_vec1 = dvector(1, mat_vec_dim);
  D_vec = dvector(1, mat_vec_dim);
  M_M_tp = dmatrix(1, mat_vec_dim, 1, mat_vec_dim);

  Id.identity();

  J.zero();
  for (k = 0; k < 3; k++) {
    J[2*k]   =  Id[2*k+1];
    J[2*k+1] = -Id[2*k];
  }

  C = Cell[globval.Cell_nLoc].S;

  if (!false) prt_rad(globval.Energy, rho);

  GetEmittance(ElemIndex("cav"), true);

  prtmat(6, Cell[0].sigma);

  globval.Cavity_on = true; globval.radiation = true;
  Ring_GetTwiss(true, 0e0);
  putlinmat(ss_dim, globval.OneTurnMat, M);
  putlinmat(6, globval.Ascr, A);
  A_A_tp = A*lin_map_tp(3, A);

  printf("\n  tau   = [%10.3e, %10.3e, %10.3e]\n",
	 globval.tau[X_], globval.tau[Y_], globval.tau[Z_]);
  printf("  eps   = [%10.3e, %10.3e, %10.3e]\n",
	 globval.eps[X_], globval.eps[Y_], globval.eps[Z_]);
  printf("  D     = [%10.3e, %10.3e, %10.3e]\n",
	 globval.D_rad[X_], globval.D_rad[Y_], globval.D_rad[Z_]);

  if (!globval.radiation)
    for (k = 0; k < ss_dim; k++)
      M[k] *= exp(globval.alpha_rad[k/2]);

  prt_lin_map(3, lin_map_tp(3, M)*J*M);

  D.zero();
  // Floquet space.
  for (k = 0; k < ss_dim; k++)
    D[k] = globval.D_rad[k/2]*Id[k];
  prt_lin_map(3, D);
  D = D*A_A_tp;
  prt_lin_map(3, D);

  get_M_M_tp(M, M_M_tp);
  vech(sigma, sigma_vec);
  vech(D, D_vec);

  get_emit(M_M_tp, D_vec);

  exit(0);

  A.identity();
  putlinmat(4, globval.Ascr, A);
  A[ct_] = sqrt(globval.beta_z)*Id[ct_];
  A[delta_] =
    -globval.alpha_z/sqrt(globval.beta_z)*Id[ct_]
    + 1e0/sqrt(globval.beta_z)*Id[delta_];

  sigma = A*lin_map_tp(3, A);
  for (k = 0; k < ss_dim; k++)
    sigma[k] *= globval.eps[k/2];

  prt_lin_map(3, sigma);
  for (k = 0; k < 1; k++)
    // sigma = M*sigma*M_tp + D;
    sigma1 = M*sigma*M_tp;
  prt_lin_map(3, sigma1);

  dmvmult(M_M_tp, mat_vec_dim, mat_vec_dim, sigma_vec, mat_vec_dim, sigma_vec1);
  dvadd(sigma_vec1, mat_vec_dim, D_vec, sigma_vec1);

  inv_vech(sigma_vec1, sigma1);
  prt_lin_map(3, sigma1);

  free_dvector(sigma_vec, 1, mat_vec_dim);
  free_dvector(sigma_vec1, 1, mat_vec_dim);
  free_dvector(D_vec, 1, mat_vec_dim);
  free_dmatrix(M_M_tp, 1, mat_vec_dim, 1, mat_vec_dim);
}


int main(int argc, char *argv[])
{
  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  globval.H_exact        = false;  globval.quad_fringe = false;
  globval.Cavity_on      = false;  globval.radiation   = false;
  globval.emittance      = false;  globval.IBS         = false;
  globval.pathlength     = false;  globval.bpm         = 0;
  globval.dip_edge_fudge = true;

  Ring_GetTwiss(true, 0e0); printglob();
  prt_lat("linlat1.out", globval.bpm, true);

  get_sigma();
  // get_lin_map();
}
