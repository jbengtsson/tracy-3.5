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


void prt_mat(const int m, const int n, double **A)
{
  int i, j;
  printf("\n");

  printf("\nmatrix:\n");
  for (i = 1; i <= m; i++) {
    for (j = 1; j <= n; j++)
      printf("%12.3e", A[i][j]);
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
      printf("%12.3e", A[i][j]);
    printf("\n");
  }
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


void get_mat(const int n_dim, double **A, std::vector< std::vector<double> > &B)
{
  int                 i, j;
  std::vector<double> row;

  for (i = 0; i < n_dim; i++) {
    row.clear();
    for (j = 0; j < n_dim; j++)
      row.push_back(A[i+1][j+1]);
    B.push_back(row);
  }
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


void prt_mat(const std::vector< std::vector<string> > &A)
{
  int i, j;

  cout << "\nmatrix:\n";
  for (i = 0; i < (int)A.size(); i++) {
    for (j = 0; j < (int)A[0].size(); j++)
      cout << " " << A[i][j];
    cout << endl;
  }
}


void Kronecker_prod(const int n_dim,
		    const std::vector< std::vector<string> > &A,
		    const std::vector< std::vector<string> > &B,
		    std::vector< std::vector<string> > &C)
{
  int                 i, j, k, l;
  std::vector<string> row;
  stringstream        str;

  for (i = 0; i < sqr(n_dim); i++) {
    row.clear();
    for (j = 0; j < sqr(n_dim); j++)
      row.push_back("0");
    C.push_back(row);
  }

  for (i = 0; i < n_dim; i++)
    for (j = 0; j < n_dim; j++)
      for (k = 0; k < n_dim; k++)
	for (l = 0; l < n_dim; l++) {
	  str.clear(); str.str("");
	  str << A[i][j] << "*" << B[k][l];
	  C[i*n_dim+k][j*n_dim+l] = str.str();
	}
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
		   std::vector< std::vector<string> > &M)
{
  // Matrix dim reduction from symmetry.
  int k;

  for (k = 0; k < (int)M.size(); k++) {
    M[k][i-1] += "+" + M[k][j-1];
    M[k][j-1] = "-";
  }
}


void red_sym_mat_2(const std::vector<int> ind,
		   std::vector< std::vector<string> > &M)
{
  // Matrix dim reduction from symmetry.
  int j, k;

  for (j = 0; j < (int)ind.size(); j++) {
    M.erase(M.begin()+ind[j]-1);
    for (k = 0; k < (int)M.size(); k++)
      M[k].erase(M[k].begin()+ind[j]-1);
  }
}


void get_M_M_tp(const int n_dim)
{
  int                                i, j, i1, j1;
  std::vector<string>                row;
  std::vector< std::vector<string> > M, Mp;
  std::vector<int>                   ind;
  stringstream                       str;

  for (i = 0; i < n_dim; i++) {
    row.clear();
    for (j = 0; j < n_dim; j++) {
      str.clear(); str.str("");
      str << "(" << i+1 << ", " << j+1 << ")";
      row.push_back(str.str());
    }
    M.push_back(row);
  }
  for (i = 0; i < n_dim; i++)
    for (j = 0; j < i; j++)
      M[i][j] = M[j][i];

  prt_mat(M);
  Kronecker_prod(n_dim, M, M, Mp);
  // prt_mat(Mp);

  printf("\n");
  for (i = 1; i <= n_dim-1; i++) {
    for (j = 1; j <= n_dim-i; j++) {
      i1 = (i-1)*(n_dim+1) + j + 1; j1 = i1 + j*(n_dim-1);
      ind.push_back(j1);
      printf(" (%2d, %2d)", i1, j1);
      red_sym_mat_1(n_dim, i1, j1, Mp);
    }
    printf("\n");
  }
  prt_mat(Mp);
  printf("\n");
  for (i = 0; i < (int)ind.size(); i++)
    printf("  %2d", ind[i]);
  bubble_sort(ind);
  printf("\n");
  for (i = 0; i < (int)ind.size(); i++)
    printf("  %2d", ind[i]);
  printf("\n");
  red_sym_mat_2(ind, Mp);
  prt_mat(Mp);
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


void get_M_M_tp(const ss_vect<tps> &M, double **M_M_tp)
{
  // Roth's Relationship (W. Roth "On Direct Product Matrices" Bul. Amer. Math.
  // Soc. 40, 461-468 (1934):
  //   vec{A B C} = (C^T x_circ A) vec(B)
  int                                i, j, i1, j1;
  std::vector<int>                   ind;
  std::vector< std::vector<double> > M_mat, M_prod;

  const int
    mat_dim          = 2,
    prod_mat_dim     = sqr(mat_dim),
    red_prod_mat_dim = mat_dim*(mat_dim-1)/2+mat_dim;

  // get_M_M_tp(mat_dim);
  // exit(0);

  get_lin_mat(mat_dim, M, M_mat);
  prt_mat(M_mat);
  Kronecker_prod(M_mat, M_mat, M_prod);
  prt_mat(M_prod);

  printf("\n");
  for (i = 1; i <= mat_dim-1; i++) {
    for (j = 1; j <= mat_dim-i; j++) {
      i1 = (i-1)*(mat_dim+1) + j + 1; j1 = i1 + j*(mat_dim-1);
      ind.push_back(j1);
      printf(" (%2d, %2d)", i1, j1);
      red_sym_mat_1(mat_dim, i1, j1, M_prod);
    }
    printf("\n");
  }
  prt_mat(M_prod);
  printf("\n");
  for (i = 0; i < (int)ind.size(); i++)
    printf("  %2d", ind[i]);
  bubble_sort(ind);
  printf("\n");
  for (i = 0; i < (int)ind.size(); i++)
    printf("  %2d", ind[i]);
  printf("\n");
  red_sym_mat_2(ind, M_prod);
  prt_mat(M_prod);
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
