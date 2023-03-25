#define NO 2

#include <random>
#include <complex>

#include "tracy_lib.h"


int
  no_tps   = NO,
  ndpt_tps = 5;


class MomentType {
private:
  double       q;     // Charge.
public:
  ss_vect<tps> Sigma; // Statistical moments for charge distribution.

  void propagate_cav(void);
  void propagate_lat(void);
  void propagate(void);
};

  
void  MomentType::propagate_cav(void)
{
  long int lastpos;

  const int loc = globval.Cell_nLoc;

  printf("\n  From: %10s\n", Cell[loc].Elem.PName);
  Cell_Pass(loc, loc, Sigma, lastpos);
}


void  MomentType::propagate_lat(void)
{
  long int lastpos;

  printf("\n  From: %10s\n  To:   %10s\n", Cell[0].Elem.PName,
	 Cell[globval.Cell_nLoc-1].Elem.PName);
  Cell_Pass(0, globval.Cell_nLoc-1, Sigma, lastpos);
}


void  MomentType::propagate(void)
{
  // Assumes that the RF cavity is at the end of the lattice.

  this->propagate_lat();
  this->propagate_cav();
}


void compute_stochastic_part
(ss_vect<double> &X, ss_vect<tps> &X_map)
{
  // Compute stochastic part of map.
  static std::default_random_engine rand;
  std::normal_distribution<double>  norm_ranf(0e0, 1e0);

  const int n_dof = 3;

  for (int k = 0; k < 2*n_dof; k++) {
    X[k] = norm_ranf(rand);
    X_map[k] = X[k]*tps(0e0, k+1);
  }
}


void propagate_rad
(tps &sigma, const ss_vect<tps> &M_Chol, const ss_vect<tps> &M_Chol_t)
{
  int             k;
  ss_vect<double> X;
  ss_vect<tps>    X_map;

  std::default_random_engine       rand;
  std::normal_distribution<double> norm_ranf(0e0, 1e0);

  compute_stochastic_part(X, X_map);
  for (k = 0; k < 6; k++)
    sigma += (M_Chol_t*X).cst()[k]*tps(0e0, k+1);
  // sigma += (M_Chol_t*X).cst() + M_Chol_t*sqr(X_map)*M_Chol;
}


void set_HOM_long
(const string &name, const double beta, const double f, const double R_sh,
 const double Q)
{
  const long int loc = Elem_GetPos(ElemIndex(name.c_str()), 1);

  Cell[loc].Elem.C->beta_long.push_back(beta);
  Cell[loc].Elem.C->HOM_f_long.push_back(f);
  Cell[loc].Elem.C->HOM_R_sh_long.push_back(R_sh);
  Cell[loc].Elem.C->HOM_Q_long.push_back(Q);
  Cell[loc].Elem.C->HOM_V_long.push_back(0e0);
}


void prt_HOM(const int n, const string &name)
{
  const long int
    loc = Elem_GetPos(ElemIndex(name.c_str()), 1);

  printf("%3d", n);
  printf("  %22.15e %22.15e %22.15e\n",
	 abs(Cell[loc].Elem.C->HOM_V_long[0]),
	 arg(Cell[loc].Elem.C->HOM_V_long[0]),
	 real(Cell[loc].Elem.C->HOM_V_long[0]));
}


void prt_sigma(ofstream &outf, const int n, const tps &sigma)
{
  const long int
    x[]           = {1, 0, 0, 0, 0, 0, 0},
    p_x[]         = {0, 1, 0, 0, 0, 0, 0},
    y[]           = {0, 0, 1, 0, 0, 0, 0},
    p_y[]         = {0, 0, 0, 1, 0, 0, 0},
    ct[]          = {0, 0, 0, 0, 0, 1, 0},
    delta[]       = {0, 0, 0, 0, 1, 0, 0},

    x_x[]         = {2, 0, 0, 0, 0, 0, 0},
    p_x_p_x[]     = {0, 2, 0, 0, 0, 0, 0},
    y_y[]         = {0, 0, 2, 0, 0, 0, 0},
    p_y_p_y[]     = {0, 0, 0, 2, 0, 0, 0},
    ct_ct[]       = {0, 0, 0, 0, 0, 2, 0},
    delta_delta[] = {0, 0, 0, 0, 2, 0, 0};

  outf << scientific << setprecision(3)
       << setw(3) << n
       << setw(11) << sigma[x] << setw(11) << sigma[p_x]
       << setw(11) << sigma[y] << setw(11) << sigma[p_y]
       << setw(11) << sigma[delta] << setw(11) << sigma[ct]/c0 << " |"
       << setw(10) << sqrt(sigma[p_x_p_x]) << setw(10) << sqrt(sigma[x_x])
       << setw(10) << sqrt(sigma[p_y_p_y]) << setw(10) << sqrt(sigma[y_y])
       << setw(10) << sqrt(sigma[ct_ct])
       << setw(10) << sqrt(sigma[delta_delta])/c0 << "\n";
}


void propagate_mag_lat(tps &sigma, const ss_vect<tps> &M_inv)
{
  sigma = sigma*M_inv;
}


double get_ct(tps &sigma)
{
  const long int jj[] = {0, 0, 0, 0, 0, 1, 0};

  return sigma[jj];
}


void propagate_cav_HOM_transv
(const int n, const double Q_b, tps &sigma)
{

  const string
    cav_name = "cav";
  const long int
    loc = Elem_GetPos(ElemIndex(cav_name.c_str()), 1);
  const double
    Circ        = Cell[globval.Cell_nLoc].S,
    beta_RF     = Cell[loc].Elem.C->beta_RF,
    f           = Cell[loc].Elem.C->HOM_f_long[0],
    R_sh        = Cell[loc].Elem.C->HOM_R_sh_long[0],
    Q           = Cell[loc].Elem.C->HOM_Q_long[0],
    k_loss      = 2e0*M_PI*f*R_sh/Q,
    R_sh_loaded = R_sh/(1e0+beta_RF),
    Q_loaded    = Q/(1e0+beta_RF);
  const complex<double>
    I = complex<double>(0e0, 1e0);

  // Update RF cavity HOM phasor.
  Cell[loc].Elem.C->HOM_V_long[0] *=
    exp(2e0*M_PI*f*(Circ+get_ct(sigma))/c0*(-1e0/(2e0*Q_loaded)+I));

  // First half increment of HOM phasor.
  Cell[loc].Elem.C->HOM_V_long[0] += Q_b*k_loss/2e0;

  // Propagate through wake field.
  sigma -= Q_b*real(Cell[loc].Elem.C->HOM_V_long[0])/2e0;

  // Second half increment of HOM phasor.
  Cell[loc].Elem.C->HOM_V_long[0] += Q_b*k_loss/2e0;

  if (false)
    prt_HOM(n, "cav");
}


void propagate_cav_HOM_long(const int n, const double Q_b, tps &sigma)
{
  ss_vect<tps> M;

  const string
    cav_name = "cav";
  const long int
    loc = Elem_GetPos(ElemIndex(cav_name.c_str()), 1);
  const double
    Circ        = Cell[globval.Cell_nLoc].S,
    beta_RF     = Cell[loc].Elem.C->beta_RF,
    f           = Cell[loc].Elem.C->HOM_f_long[0],
    R_sh        = Cell[loc].Elem.C->HOM_R_sh_long[0],
    Q           = Cell[loc].Elem.C->HOM_Q_long[0],
    k_loss      = 2e0*M_PI*f*R_sh/Q,
    R_sh_loaded = R_sh/(1e0+beta_RF),
    Q_loaded    = Q/(1e0+beta_RF);
  const complex<double>
    I = complex<double>(0e0, 1e0);

  // Update RF cavity HOM phasor.
  Cell[loc].Elem.C->HOM_V_long[0] *=
    exp(2e0*M_PI*f*(Circ+get_ct(sigma))/c0*(-1e0/(2e0*Q_loaded)+I));

  // First half increment of HOM phasor.
  Cell[loc].Elem.C->HOM_V_long[0] += Q_b*k_loss/2e0;

  // Propagate through wake field.
  sigma -= Q_b*real(Cell[loc].Elem.C->HOM_V_long[0])/2e0*tps(0e0, delta_+1);

  // Second half increment of HOM phasor.
  Cell[loc].Elem.C->HOM_V_long[0] += Q_b*k_loss/2e0;

  if (!false)
    prt_HOM(n, "cav");
}


void propagate_lat
(const int n, const double Q_b, tps &sigma, const ss_vect<tps> &M_inv,
 const ss_vect<tps> &M_Chol, const ss_vect<tps> &M_Chol_t)
{
  if (!true)
    propagate_cav_HOM_long(n, Q_b, sigma);

  // Propagate through magnetic lattice.
  propagate_mag_lat(sigma, M_inv);

  if (!false)
    // Radiate.
    propagate_rad(sigma, M_Chol, M_Chol_t);
}


tps compute_sigma
(const double eps[], const double sigma_s, const double sigma_delta,
 const ss_vect<tps> &A)
{
  int          k;
  tps          sigma;
  ss_vect<tps> Id;

  Id.identity();
  sigma = 0e0;
  for (k = 0; k < 2; k++)
    sigma += eps[k]*(sqr(Id[2*k])+sqr(Id[2*k+1]));
  sigma = sigma*Inv(A);
  sigma += sqr(sigma_delta*Id[ct_]) + sqr(sigma_s*Id[delta_]);
  return sigma;
}


istream& operator>>(std::istream &is, ss_vect<tps> &a)
{
  int k;

  const int ps_dim = 2*nd_tps;

  a.identity();
  for (k = 0; k < 2*nd_tps; k++)
    is >> a[k];
  a[2*nd_tps] = 0e0;
  return is;
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


void rd_maps(ss_vect<tps> &M, ss_vect<tps> &A, ss_vect<tps> &M_Chol)
{
  ifstream inf;

  const string
    file_name = "../compute_maps";

  file_rd(inf, (file_name+"_M.dat").c_str());
  inf >> M;
  inf.close();
  prt_map(3, "\nM:", M);

  file_rd(inf, (file_name+"_A.dat").c_str());
  inf >> A;
  inf.close();
  prt_map(3, "\nA:", A);

  file_rd(inf, (file_name+"_M_Chol.dat").c_str());
  inf >> M_Chol;
  inf.close();
  prt_map(3, "\nM_Chol:", M_Chol);
}


ss_vect<tps> tp_map(const int n_dof, const ss_vect<tps> &A)
{
  // Matrix transpose.
  Matrix A_mat;

  getlinmat(2*n_dof, A, A_mat);
  TpMat(2*n_dof, A_mat);
  return putlinmat(2*n_dof, A_mat);
}


void test_case_2(const string &name)
{
  // Single bunch, 1 longitudinal HOM.
  long int        lastpos;
  int             k;
  tps             sigma;
  ss_vect<tps>    M, M_inv, A, M_Chol, M_Chol_t;
  ofstream        outf;

  const bool
    rad         = !false;

  const int
    n_turn      = 9;

  const double
    eps[]       = {161.7e-12, 8e-12},
    sigma_s     = 3.739e-3,
    sigma_delta = 9.353e-04,

    Q_b         = -0.6e-9,
    beta_HOM    = 1e0,
    f           = 1e9,
    R_sh        = 1e3,
    Q           = 1e8;

  const string
    file_name = "moments.out";

  file_wr(outf, file_name.c_str());

  globval.Cavity_on = globval.radiation = rad;
  globval.pathlength = false;

  danot_(1);

  if (!rad) {
    M.identity();
    Cell_Pass(0, globval.Cell_nLoc, M, lastpos);
    outf << M << "\n";
    outf.close();
    exit(0);
  } else {
    rd_maps(M, A, M_Chol);
    M_Chol_t = tp_map(3, M_Chol);
  }
  M_inv = Inv(M);

  danot_(NO);

  set_HOM_long("cav", beta_HOM, f, R_sh, Q);

  sigma = compute_sigma(eps, sigma_s, sigma_delta, A);

  prt_sigma(outf, 0, sigma);
  printf("\n");
  for (k = 1; k <= n_turn; k++) {
    propagate_lat(k, Q_b, sigma, M_inv, M_Chol, M_Chol_t);
    prt_sigma(outf, k, sigma);
  }
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

  // disable from TPSALib- and LieLib log messages
  idprset_(-1);

  if (false)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  set_state();

  if (false) {
    danot_(1);

    Ring_GetTwiss(true, 0e0);
    printglob();

    if (!false)
      GetEmittance(ElemIndex("cav"), true);

    prt_lat("linlat1.out", globval.bpm, true);
    prt_lat("linlat.out", globval.bpm, true, 10);
  }

  if (!false)
    test_case_2("cav");
}
