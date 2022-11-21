#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


class SigmaType {
private:
  int
    n_dof;
  ss_vect<double>
    mean;              // 1st moment.
  ss_vect<tps>
    sigma,             // 2nd moment; beam envelope.
    M,                 // Poincaré map.
    M_t,
    A,
    A_t,
    A_inv, 
    A_t_inv;
  FILE
    *outf;             // Output file.
public:
  void set_n_dof(const int n_dof) { this->n_dof = n_dof; }
  ss_vect<double> get_eps(void) const;
  void set_file_name(const string &file_name);
  void prt_eps(const int n) const;
  void prt_sigma(void);
  void propagate(const int n);
  void init_sigma(const double eps[]);
  void get_M(void);
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
  Matrix A_mat;

  getlinmat(2*n_dof, A, A_mat);
  return DetMat(2*n_dof, A_mat);
}


ss_vect<tps> tp_map(const int n_dof, const ss_vect<tps> &A)
{
  Matrix A_mat;

  getlinmat(2*n_dof, A, A_mat);
  TpMat(2*n_dof, A_mat);
  return putlinmat(2*n_dof, A_mat);
}


void SigmaType::set_file_name(const string &file_name)
{
  outf = file_write(file_name.c_str());
}


ss_vect<double> SigmaType::get_eps(void) const
{
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


void SigmaType::propagate(const int n)
{
  long int lastpos;
  int      k;

  for (k = 0; k < n; k++) {
    if (true) {
      sigma = M*sigma*M_t;
    } else {
      Cell_Pass(0, globval.Cell_nLoc, sigma, lastpos);
      sigma = tp_map(n_dof, sigma);
      Cell_Pass(0, globval.Cell_nLoc, sigma, lastpos);
    }
    this->prt_eps(k);
  }
}


void SigmaType::init_sigma(const double eps[])
{
  int          k;
  ss_vect<tps> Id;

  Id.identity();
  sigma.zero();
  for (k = 0; k < 2*n_dof; k++)
    sigma[k] += eps[k/2]*Id[k];
  sigma = A*sigma*A_t;
}


void SigmaType::get_M(void)
{
  long int lastpos;
  double   dnu[3];

  const bool prt = false;

  globval.radiation = true;
  globval.Cavity_on = true;

  Ring_GetTwiss(true, 0e0);
  printglob();

  globval.U0 = globval.dE*1e9*globval.Energy;
  printf("\nU0 [keV] = %3.1f\n", 1e-3*globval.U0);
  printf("alpha    = [%10.3e, %10.3e, %10.3e]\n",
	 globval.alpha_rad[X_], globval.alpha_rad[Y_], globval.alpha_rad[Z_]);
  
  A = get_A_CS(n_dof, putlinmat(6, globval.Ascr), dnu);
  A_t = tp_map(n_dof, A);
  A_inv = Inv(A);
  A_t_inv = Inv(A_t);

  M.identity();
  M += globval.CODvect;
  Cell_Pass(0, globval.Cell_nLoc, M, lastpos);
  M -= globval.CODvect;
  M_t = tp_map(n_dof, M);

  printf("Det{M}-1 = %10.3e\n", det_map(n_dof, M)-1e0);

  if (prt) {
    prt_map(n_dof, "\nA:", A);
    prt_map(n_dof, "\nA^T:", A_t);
    prt_map(n_dof, "\nA^-1:", A_inv);
    prt_map(n_dof, "\n(A^T)^-1:", A_t_inv);
    prt_map(n_dof, "\nM:", M);
    prt_map(n_dof, "\nM^T:", M_t);
  }
}


void track(void)
{
  SigmaType s;

  const string
    file_name = "wake_field.out";
  const int
    n_dof     = 3,
    n         = 30000;
  const double
    eps0[]    = {1e-9, 0.2e-9, 1e-3};

  s.set_n_dof(n_dof);
  s.get_M();
  s.init_sigma(eps0);
  s.set_file_name(file_name);
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
