#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


const int mat_vec_dim = 3;


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


void get_lin_mat(const int n_dim, const ss_vect<tps> &map, double **mat)
{
  int i, j;

  for (i = 0; i < n_dim; i++)
    for (j = 0; j < n_dim; j++)
      mat[i+1][j+1] = map[i][j];
}


void Kronecker_prod(const int n_dim, double **A, const int B_dim,
		    double **B, double **C)
{
  int i, j, k, l;

  for (i = 1; i <= n_dim; i++)
    for (j = 1; j <= n_dim; j++)
      for (k = 1; k <= n_dim; k++)
	for (l = 1; l <= n_dim; l++)
	  C[(i-1)*n_dim+k][(j-1)*n_dim+l] = A[i][j]*B[k][l];
}


void red_sym_mat(double **A, double **B)
{
}


void vech(const ss_vect<tps> &M, double *M_vec)
{
  // Symmetric matrix vectorization. 
  int i, j, k;

  k = 0;
  for (i = 1; i <= 2; i++)
    for (j = i; j <= 2; j++) {
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
  i = j = k = 0;
  for (i = 1; i <= mat_vec_dim; i++) {
    k++; j++;
    if (j == ss_dim) {
      i++; j = 0;
    }
    M[i] += M_vec[k]*Id[i];
    if (i != j) M[j] += M_vec[k]*Id[j];
  }
}


void get_M_M_tp(const ss_vect<tps> &M, double **M_M_tp)
{
  // Roth's Relationship (W. Roth "On Direct Product Matrices" Bul. Amer. Math.
  // Soc. 40, 461-468 (1934):
  //   vec{A B C} = (C^T x_circ A) vec(B)
  int    i, j;
  double **M_mat, **M_prod, *M_tp_vec;

  int
    mat_dim      = 2,
    prod_mat_dim = sqr(mat_dim);

  M_mat = dmatrix(1, mat_dim, 1, mat_dim);
  M_prod = dmatrix(1, prod_mat_dim, 1, prod_mat_dim);

  get_lin_mat(mat_dim, M, M_mat);
  dmdump(stdout, "\nmatrix:\n", M_mat, mat_dim, mat_dim, "%11.3e");
  Kronecker_prod(mat_dim, M_mat, mat_dim, M_mat, M_prod);
  dmdump(stdout, "\nmatrix:\n", M_prod, prod_mat_dim, prod_mat_dim, "%11.3e");
  free_dmatrix(M_mat, 1, mat_dim, 1, mat_dim);
  free_dmatrix(M_prod, 1, prod_mat_dim, 1, prod_mat_dim);
}


void get_sigma(void)
{
  int          k;
  double       C, **M_M_tp;
  ss_vect<tps> Id, J, A, D, sigma, M, M_tp;

  const int n_turn = 500000;

  M_M_tp = dmatrix(1, mat_vec_dim, 1, mat_vec_dim);

  Id.identity();

  J.zero();
  for (k = 0; k < 3; k++) {
    J[2*k] = Id[2*k+1];
    J[2*k+1] = -Id[2*k];
  }

  C = Cell[globval.Cell_nLoc].S;

  GetEmittance(ElemIndex("cav"), true);

  globval.Cavity_on = true; globval.radiation = true;
  Ring_GetTwiss(true, 0e0); printglob();

  printf("\n  alpha = [%10.3e, %10.3e, %10.3e]\n",
	 globval.alpha_rad[X_], globval.alpha_rad[Y_], globval.alpha_rad[Z_]);
  printf("  tau   = [%10.3e, %10.3e, %10.3e]\n",
	 -C/(c0*globval.alpha_rad[X_]), -C/(c0*globval.alpha_rad[Y_]),
	 -C/(c0*globval.alpha_rad[Z_]));
  printf("  D     = [%10.3e, %10.3e, %10.3e]\n",
	 globval.D_rad[X_], globval.D_rad[Y_], globval.D_rad[Z_]);

  putlinmat(ss_dim, globval.OneTurnMat, M);
  for (k = 0; k < ss_dim; k++)
    M[k] *= exp(globval.alpha_rad[k/2]);
  M_tp = lin_map_tp(3, M);

  D.zero();
  D[px_] = 4.015e0/7e0*globval.D_rad[X_]*Id[px_];
  D[delta_] = 4e0*M_PI*globval.D_rad[Z_]*Id[delta_];
  prt_lin_map(3, D);

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

  get_M_M_tp(M, M_M_tp);

  // for (k = 0; k < n_turn; k++)
  //   sigma = M*sigma*M_tp + D;
  // prt_lin_map(3, sigma);

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
