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
  printf("\n  %9.3e %9.3e\n", sqrt(sigma[x_][x_]), sqrt(sigma[px_][px_]));
  printf("  %9.3e %9.3e\n", sqrt(sigma[ct_][ct_]), sqrt(sigma[delta_][delta_]));

  free_dvector(sigma_vec, 1, mat_vec_dim);
  free_dmatrix(Id, 1, mat_vec_dim, 1, mat_vec_dim);
  free_dmatrix(MmI, 1, mat_vec_dim, 1, mat_vec_dim);
  free_dmatrix(MmI_inv, 1, mat_vec_dim, 1, mat_vec_dim);
}


void get_sigma(void)
{
  int          k;
  double       C_u, C_gamma, C_q, gamma, P_gamma, eps_c, sigma_delta, D_delta;
  double       C, gamma_s, *sigma_vec, *sigma_vec1, *D_vec, **M_M_tp;
  ss_vect<tps> Id, J, A, D, sigma, sigma1, M, M_tp;

  const int    n_turn = 500000;
  const double
    h_bar = 6.582119e-16,
    alpha = 1e0/137.035999084,
    r_e   = 2.817940e-15,
    rho   = 2.62/(6.0*M_PI/180.0);

  printf("\nget_sigma: n = %d\n", mat_vec_dim);

  sigma_vec = dvector(1, mat_vec_dim);
  sigma_vec1 = dvector(1, mat_vec_dim);
  D_vec = dvector(1, mat_vec_dim);
  M_M_tp = dmatrix(1, mat_vec_dim, 1, mat_vec_dim);

  Id.identity();

  J.zero();
  for (k = 0; k < 3; k++) {
    J[2*k] = Id[2*k+1];
    J[2*k+1] = -Id[2*k];
  }

  C = Cell[globval.Cell_nLoc].S;

  GetEmittance(ElemIndex("cav"), true);

  prtmat(6, Cell[0].sigma);

  globval.Cavity_on = true; globval.radiation = true;
  Ring_GetTwiss(true, 0e0);
  putlinmat(ss_dim, globval.OneTurnMat, M);
  A.identity();
  putlinmat(4, globval.Ascr, A);

  printf("\n  alpha = [%10.3e, %10.3e, %10.3e]\n",
	 globval.alpha_rad[X_], globval.alpha_rad[Y_], globval.alpha_rad[Z_]);
  printf("  tau   = [%10.3e, %10.3e, %10.3e]\n",
	 -C/(c0*globval.alpha_rad[X_]), -C/(c0*globval.alpha_rad[Y_]),
	 -C/(c0*globval.alpha_rad[Z_]));
  printf("  eps   = [%10.3e, %10.3e, %10.3e]\n",
	 -globval.D_rad[X_]/(2.0*globval.alpha_rad[X_]),
	 -globval.D_rad[Y_]/(2.0*globval.alpha_rad[Y_]),
	 -globval.D_rad[Z_]/(2.0*globval.alpha_rad[Z_]));
  printf("  D     = [%10.3e, %10.3e, %10.3e]\n",
	 globval.D_rad[X_], globval.D_rad[Y_], globval.D_rad[Z_]);

  gamma_s = (1e0+sqr(globval.alpha_z))/globval.beta_z;
  printf("\nD_delta = %9.3e\n", gamma_s*globval.D_rad[Z_]);

  gamma = 1e9*globval.Energy/m_e;
  C_u = 55e0/(24e0*sqrt(3e0));
  C_gamma = 4e0*M_PI*r_e/(3e0*cube(m_e));
  C_q = 3e0*C_u*h_bar*c0/(4e0*m_e);
  eps_c = 3e0*h_bar*c0*cube(gamma)/(2e0*rho);
  P_gamma = C_gamma*c0*pow(1e9*globval.Energy, 4)/(2e0*M_PI*sqr(rho));
#if 0
  // Correct.
  sigma_delta = sqrt(C_q/(globval.J[Z_]*rho))*gamma;
#else
  // Correct.
  sigma_delta = sqrt(C_u*eps_c/(2e0*globval.J[Z_]*1e9*globval.Energy));
#endif
#if 0
  D_delta = C_u*P_gamma*eps_c/sqr(1e9*globval.Energy);
#else
  D_delta =
    2e0*C_u*eps_c/(2e0*globval.J[Z_]*1e9*globval.Energy*globval.tau[Z_])*C/c0;
#endif

  printf("\nm_e               = %9.3e q_e = %9.3e\n", m_e, q_e);
  printf("E [GeV]           = %9.3e gamma = %9.3e rho = %9.3e\n",
	 globval.Energy, gamma, rho);
  printf("C_gamma [m/GeV^3] = %9.3e\n", 1e27*C_gamma);
  printf("C_q [m]           = %9.3e\n", C_q);
  printf("eps_c [keV]       = %9.3e\n", 1e-3*eps_c);
  printf("P_gamma           = %9.3e\n", P_gamma);
  printf("sigma_delta       = %9.3e\n", sigma_delta);
  printf("D_delta           = %9.3e\n", D_delta);

#if 0
  for (k = 0; k < ss_dim; k++)
    M[k] *= exp(globval.alpha_rad[k/2]);
#endif
  prt_lin_map(3, lin_map_tp(3, M)*J*M);

  D.zero();
  // Floquet space.
  D[x_] = globval.D_rad[X_]*Id[x_];
  D[px_] = globval.D_rad[X_]*Id[px_];
  D[delta_] = 2e0*M_PI*globval.D_rad[Z_]*Id[delta_];
  prt_lin_map(3, D);
  D = A*D;
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
}
