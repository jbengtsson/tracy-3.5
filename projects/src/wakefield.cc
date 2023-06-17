
#include <random>
#include <complex>

#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


class BeamType;


class PoincareMapType {
private:
  int
    n_dof,             // Degrees of freedom.
    n_harm,            // RF cavity harmonic number.
    cav_fnum;          // RF cavity family number;
  CavityType
    *cav_ptr;          // Pointer to RF cavity.
  double
    C,                 // Circumference.
    alpha_rad[3],      // Damping coefficients.
    tau[3],            // Damping times: -C/(c0*alpha_rad).
    D[3],              // Diffusion coefficients.
    eps[3],            // Eigen emittances: eps = D*tau*c0/(2e0*C).
    sigma_s,           // Bunch length.
    sigma_delta;       // Momentum spread.
  ss_vect<double>
    fixed_point;       // Closed orbit.
  ss_vect<tps>
    M,                 // Poincaré map: deterministic part.
    M_diff,            // Diffusion matrix: stochastic part.

  // Cashe for the maps utilised by the computations.
  // Poincaré map.
    M_t,               // M^T.
    A,                 // Diagonalized: M = A R A^-1.
    A_t,               // A^T.
    A_inv,             // A^-1.
    A_t_inv,           // (A^T)^-1.
    R,                 // Phase-space rotation by 2*pi*nu.
  // Diffusion matrix.
    M_Chol,            // Cholesky decomposition of the diffusion matrix:
    M_Chol_t;          //   D = L^T L.
public:
  friend class BeamType;

  void set_params(const int n_dof, const double C);

  void compute_M(void);
  void compute_A(void);
  void compute_tau(void);
  void compute_D(void);
  void compute_eps(void);
  void compute_bunch_size(void);
  void compute_M_diff(void);
  void compute_M_Chol_t(void);
  void compute_maps(void);

  void compute_stochastic_part(ss_vect<double> &X, ss_vect<tps> &X_map);
  void propagate(const int n, BeamType &beam);

  void set_rf_cav_hom(const string &fam_name, const double f, const double Z,
		      const double Q);
  void propagate_rf_cav_hom(BeamType &beam);
};


class BeamType {
private:
  double
    Q_b;               // Bunch charge.
  ss_vect<double>
    mean;              // 1st moment: barycentre.
  ss_vect<tps>
    sigma;             // 2nd moment: beam envelope.
  FILE
    *outf;             // Output file.
public:
  friend class PoincareMapType;

  void set_file_name(const string &file_name);
  ss_vect<double> compute_eps(const PoincareMapType &map) const;
  void prt_eps(const int n, const PoincareMapType &map) const;
  void prt_beam(const int n, const PoincareMapType &map) const;
  void prt_sigma(const PoincareMapType &map);
  void init_sigma(const double eps_x, const double eps_y, const double sigma_s,
		  const double sigma_delta, const PoincareMapType &map);
};


void prt_vec(const int n_dof, const string &str, ss_vect<double> vec)
{
  const int n_dec = 6;

  printf("%s\n", str.c_str());
  for (int k = 0; k < 2*n_dof; k++)
    printf("%*.*e", n_dec+8, n_dec, vec[k]);
  printf("\n");
}


void prt_map(const int n_dof, const string &str, ss_vect<tps> map)
{
  const int n_dec = 6;

  printf("%s\n", str.c_str());
  for (int i = 0; i < 2*n_dof; i++) {
    for (int j = 0; j < 2*n_dof; j++)
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
  // Matrix determinant.
  Matrix A_mat;

  getlinmat(2*n_dof, A, A_mat);
  return DetMat(2*n_dof, A_mat);
}


ss_vect<tps> tp_map(const int n_dof, const ss_vect<tps> &A)
{
  // Matrix transpose.
  Matrix A_mat;

  getlinmat(2*n_dof, A, A_mat);
  TpMat(2*n_dof, A_mat);
  return putlinmat(2*n_dof, A_mat);
}


void PoincareMapType::set_params(const int n_dof, const double C)
{
  this->n_dof = n_dof;
  this->C = C;
}


void PoincareMapType::compute_M(void)
{
  // Compute the Poincaré map for the fixed point.
  long int lastpos;

  M.identity();
  M += fixed_point;
  Cell_Pass(0, globval.Cell_nLoc, M, lastpos);
  M -= fixed_point;
}


void PoincareMapType::compute_A(void)
{
  // Poincaré map Diagonalization.
  // The damping coeffients alpha[] are obtained from the eigen values.
  double nu_s, alpha_c;
  Matrix M_mat, A_mat, A_inv_mat, R_mat;

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
}


void PoincareMapType::compute_tau(void)
{
  // Compute the damping times.
  for (int k = 0; k < n_dof; k++)
    tau[k] = -C/(c0*alpha_rad[k]);
}


void PoincareMapType::compute_D(void)
{
  // Compute the diffusion coefficients.
  int long     lastpos;
  ss_vect<tps> As, D_diag;

  globval.Cavity_on = globval.radiation = globval.emittance = true;

  As = A + fixed_point;
  Cell_Pass(0, globval.Cell_nLoc, As, lastpos);

  globval.emittance = false;

  for (int k = 0; k < n_dof; k++)
    D[k] = globval.D_rad[k];
}

void PoincareMapType::compute_eps(void)
{
  // Compute the eigen emittances.
  for (int k = 0; k < n_dof; k++)
    eps[k] = D[k]*tau[k]*c0/(2e0*C);
}


void PoincareMapType::compute_bunch_size(void)
{
  // Compute the bunch size: sigma_s & sigma_delta.
  double alpha_s, beta_s, gamma_s;

  alpha_s = -A[ct_][ct_]*A[delta_][ct_] - A[ct_][delta_]*A[delta_][delta_];
  beta_s = sqr(A[ct_][ct_]) + sqr(A[ct_][delta_]);
  gamma_s = (1e0+sqr(alpha_s))/beta_s;

  sigma_s = sqrt(beta_s*eps[Z_]);
  sigma_delta = sqrt(gamma_s*eps[Z_]);
}


void PoincareMapType::compute_M_diff(void)
{
  // Compute the diffusion matrix.
  ss_vect<tps> As, D_diag;

  D_diag.zero();
  for (int k = 0; k < 2*n_dof; k++)
    D_diag[k] = D[k/2]*tps(0e0, k+1);
  M_diff = A*D_diag*tp_map(n_dof, A);
}


void PoincareMapType::compute_M_Chol_t(void)
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


void PoincareMapType::compute_maps(void)
{
  long int lastpos;
  double   dnu[3];

  globval.radiation = true;
  globval.Cavity_on = true;

  getcod(0e0, lastpos);
  fixed_point = globval.CODvect;
  prt_vec(n_dof, "\nFixed Point:", fixed_point);

  compute_M();
  M_t = tp_map(n_dof, M);
  prt_map(n_dof, "\nM:", M);

  compute_A();
  compute_tau();
  A = get_A_CS(n_dof, A, dnu);
  A_t = tp_map(n_dof, A);
  A_inv = Inv(A);
  A_t_inv = Inv(A_t);

  prt_map(n_dof, "\nA:", A);
  printf("\ntau [msec]   = [%5.3f, %5.3f, %5.3f]\n",
	 1e3*tau[X_], 1e3*tau[Y_], 1e3*tau[Z_]);

  compute_D();
  compute_eps();
  compute_bunch_size();
  compute_M_diff();

  printf("D            = [%9.3e, %9.3e, %9.3e]\n", D[X_], D[Y_], D[Z_]);
  printf("eps          = [%9.3e, %9.3e, %9.3e]\n", eps[X_], eps[Y_], eps[Z_]);
  printf("sigma_s [mm] = %5.3f\n", 1e3*sigma_s);
  printf("sigma_delta  = %9.3e\n", sigma_delta);
  prt_map(n_dof, "\nDiffusion Matrix:", M_diff);

  compute_M_Chol_t();
  M_Chol = tp_map(n_dof, M_Chol_t);
  prt_map(n_dof, "\nCholesky Decomposition:", M_Chol);
}


void PoincareMapType::compute_stochastic_part
(ss_vect<double> &X, ss_vect<tps> &X_map)
{
  // Compute stochastic part of map.
  static std::default_random_engine rand;
  std::normal_distribution<double>  norm_ranf(0e0, 1e0);

  for (int k = 0; k < 2*n_dof; k++) {
    X[k] = norm_ranf(rand);
    X_map[k] = X[k]*tps(0e0, k+1);
  }
}


void PoincareMapType::propagate(const int n, BeamType &beam)
{
  int             k;
  ss_vect<double> X;
  ss_vect<tps>    X_map;

  std::default_random_engine       rand;
  std::normal_distribution<double> norm_ranf(0e0, 1e0);

  for (k = 1; k <= n; k++) {
    compute_stochastic_part(X, X_map);
    // Average of stochastic term is zero.
    beam.mean = (M*beam.mean).cst() + (M_Chol_t*X).cst();
    beam.sigma = M*tp_map(n_dof, M*beam.sigma) + M_Chol_t*sqr(X_map)*M_Chol;

    propagate_rf_cav_hom(beam);

    beam.prt_beam(k, *this);
  }
}


void PoincareMapType::set_rf_cav_hom
(const string &fam_name, const double f, const double Z, const double Q)
{
  long int cav_loc;

  cav_fnum = ElemIndex(fam_name.c_str());
  cav_loc = Elem_GetPos(cav_fnum, 1);
  cav_ptr = Cell[cav_loc].Elem.C;

  cav_ptr->HOM_f_trans[0].push_back(f);
  cav_ptr->HOM_R_sh_trans[0].push_back(Z);
  cav_ptr->HOM_Q_trans[0].push_back(Q);
  cav_ptr->HOM_V_trans[0].push_back(0e0);

  printf("\nset_rf_cav_hom: %.6s h = %d\n",
	 Cell[cav_loc].Elem.PName, cav_ptr->harm_num);
}


void PoincareMapType::propagate_rf_cav_hom(BeamType &beam)
{
  if (false)
    cout << "propagate_rf_cav_hom: V = " << cav_ptr->HOM_V_trans[X_][0] << "\n";
}


void BeamType::set_file_name(const string &file_name)
{
  outf = file_write(file_name.c_str());
}


ss_vect<double> BeamType::compute_eps(const PoincareMapType &map) const
{
  // Compute the eigen emittances.
  ss_vect<double> eps;
  ss_vect<tps>    diag;

  diag = map.A_inv*sigma*map.A_t_inv;
  for (int k = 0; k < 2*map.n_dof; k++)
    eps[k] = diag[k][k];
  return eps;
}


void BeamType::prt_eps(const int n, const PoincareMapType &map) const
{
  ss_vect<double> eps;

  eps = compute_eps(map);
  fprintf(outf, "%7d", n);
  for (int k = 0; k < 2*map.n_dof; k++)
    fprintf(outf, "%13.5e", eps[k]);
  fprintf(outf, "%13.5e %13.5e\n",
	  sqrt(sigma[delta_][delta_]), sqrt(sigma[ct_][ct_]));
}


void BeamType::prt_beam(const int n, const PoincareMapType &map) const
{
  fprintf(outf, "%7d", n);
  for (int k = 0; k < 2*map.n_dof; k++)
    fprintf(outf, "%13.5e", mean[k]);
  for (int k = 0; k < 2*map.n_dof; k++)
    fprintf(outf, "%13.5e", sqrt(sigma[k][k]));
  fprintf(outf, "\n");
}


void BeamType::prt_sigma(const PoincareMapType &map)
{
  ss_vect<double> eps;

  eps = compute_eps(map);
  prt_vec(map.n_dof, "\neps:", eps);
  prt_map(map.n_dof, "\nsigma:", sigma);
  printf("\nsigma_kk:\n ");
  for (int k = 0; k < 2*map.n_dof; k++)
    printf(" %11.5e", sqrt(sigma[k][k]));
  printf("\n");
}


void BeamType::init_sigma
(const double eps_x, const double eps_y, const double sigma_s,
 const double sigma_delta, const PoincareMapType &map)
{
  const double eps0[] = {eps_x, eps_y};

  mean.zero();

  sigma.zero();
  for (int k = 0; k < 4; k++)
    sigma[k] += eps0[k/2]*tps(0e0, k+1);
  sigma = map.A*sigma*map.A_t;
  sigma[ct_] += sqr(sigma_s)*tps(0e0, ct_+1);
  sigma[delta_] += sqr(sigma_delta)*tps(0e0, delta_+1);
}


void track(const int n, const double eps[], const double sigma_s,
	   const double sigma_delta)
{
  PoincareMapType map;
  BeamType        beam;

  const string
    file_name = "wakefield.out";
  const int
    n_dof     = 3;

  map.set_params(n_dof, Cell[globval.Cell_nLoc].S);
  map.compute_maps();
  map.set_rf_cav_hom("cav", 800e6, 1e6, 4.8e4);

  beam.set_file_name(file_name);

  beam.init_sigma(eps[X_], eps[Y_], sigma_s, sigma_delta, map);

  map.propagate(n, beam);

  beam.prt_sigma(map);
}


void set_lat_state(void)
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
  double gamma_z, sigma_s, sigma_delta;

  reverse_elem     = true;
  globval.mat_meth = false;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  set_lat_state();

  Ring_GetTwiss(true, 0e0);
  printglob();

  GetEmittance(ElemIndex("cav"), true);
  gamma_z     = (1e0+sqr(globval.alpha_z))/globval.beta_z;
  sigma_s     = sqrt(globval.beta_z*globval.eps[Z_]);
  sigma_delta = sqrt(gamma_z*globval.eps[Z_]);
  printf("\nsigma_x [microns] = %5.3f\n",
	 1e6*sqrt(Cell[0].Beta[X_]*globval.eps[X_]));
  printf("sigma_y [microns] = %5.3f\n",
	 1e6*sqrt(Cell[0].Beta[Y_]*globval.eps[Y_]));
  printf("sigma_s [ps]      = %5.3f\n", 1e12*sigma_s/c0);
  printf("sigma_delta       = %9.3e\n", sigma_delta);

  track(50000, globval.eps, sigma_s, sigma_delta);
}
