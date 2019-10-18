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


void sym_mat_to_vec(const ss_vect<tps> &M, double *M_vec,
		    stringstream str[], const bool tp)
{
  // Symmetric matrix vectorization.
  int i, j, k;

  k = 0;
  for (i = 1; i <= 2; i++)
    for (j = i; j <= 2; j++) {
      k++;
      if (!tp) {
	M_vec[k] = M[i][j];
	str[k-1] << "(" << i << j << ")";
      } else {
	M_vec[k] = M[j][i];
	str[k-1] << "(" << j << i << ")";
      }
    }

  cout << "\n";
  for (i = 0; i < mat_vec_dim; i++)
    cout << " " << str[i].str();
  cout << "\n";
}


void vec_sym_mat(const double *M_vec, ss_vect<tps> &M)
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
  int    i, j;
  double *M_vec, *M_tp_vec;
  stringstream
    str1[mat_vec_dim], str2[mat_vec_dim], str3[mat_vec_dim][mat_vec_dim];


  M_vec = dvector(1, mat_vec_dim);
  M_tp_vec = dvector(1, mat_vec_dim);

  sym_mat_to_vec(M, M_vec, str1, false);
  sym_mat_to_vec(lin_map_tp(3, M), M_tp_vec, str2, false);
  for (i = 1; i <= mat_vec_dim; i++)
    for (j = 1; j <= mat_vec_dim; j++) {
      M_M_tp[i][j] = M_vec[i]*M_tp_vec[j];
      str3[i-1][j-1] << str1[i-1].str() << str2[j-1].str();
    }

  cout << "\n";
  for (i = 0; i < mat_vec_dim; i++) {
    for (j = 0; j < mat_vec_dim; j++)
      cout << " " << str3[i][j].str();
    cout << "\n\n";
  }

  free_dvector(M_vec, 1, mat_vec_dim);
  free_dvector(M_tp_vec, 1, mat_vec_dim);
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
