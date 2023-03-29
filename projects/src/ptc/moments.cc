#define NO 2

#include <random>
#include <complex>

#include "tracy_lib.h"


int
  no_tps   = NO,
  ndpt_tps = 5;

// Moment indices.
const long int
  x[]           = {0, 1, 0, 0, 0, 0, 0},
  p_x[]         = {1, 0, 0, 0, 0, 0, 0},
  y[]           = {0, 0, 0, 1, 0, 0, 0},
  p_y[]         = {0, 0, 1, 0, 0, 0, 0},
  ct[]          = {0, 0, 0, 0, 1, 0, 0},
  delta[]       = {0, 0, 0, 0, 0, 1, 0},

  x_x[]         = {0, 2, 0, 0, 0, 0, 0},
  x_p_x[]       = {1, 1, 0, 0, 0, 0, 0},
  p_x_p_x[]     = {2, 0, 0, 0, 0, 0, 0},
  y_y[]         = {0, 0, 0, 2, 0, 0, 0},
  y_p_y[]       = {0, 0, 1, 1, 0, 0, 0},
  p_y_p_y[]     = {0, 0, 2, 0, 0, 0, 0},
  ct_ct[]       = {0, 0, 0, 0, 2, 0, 0},
  delta_delta[] = {0, 0, 0, 0, 0, 2, 0};
const int
  ps_index[]    = {px_, x_, py_, y_, ct_, delta_},
  ps_sign[]     = {-1, 1, -1, 1, -1, 1};


class MomentType {
  const string
    file_name = "moments.out"; // Output file name.

private:
  int
    n;
  double
    Q_b,       // Bunch charge.
    Circ,      // Circumference.
    phi_RF;    // RF phase [rad].
  ss_vect<tps>
    A,         // Transformation from Floquet to phase space without radiation.
    A_inv,
    M_delta,
    M_tau,
    M_Chol,
    M_Chol_t;
  ofstream
    outf;
public:
  long int
    cav_loc;
  ss_vect<tps>
    M,         // Poincaré map for lattice without radiation.
    M_inv;
  tps
    sigma;     // Statistical moments for charge distribution.

  void rd_maps(void);
  void init(const double Q_b, const double phi_RF, const string &cav_name);
  void set_HOM_long(const double beta, const double f, const double R_sh,
		    const double Q);
  void compute_sigma(const double eps[], const double sigma_s,
		     const double sigma_delta);
  void print_sigma(const int n);
  ss_vect<tps> compute_cav_HOM_long_M(const tps ct);
  void propagate_cav_HOM_long(const int n);
  void propagate_cav_HOM_transv(const int n);
  ss_vect<tps> compute_M_cav(void);
  void propagate_cav(void);
  void propagate_mag_lat(void);
  void propagate_delta(void);
  void propagate_tau(void);
  void propagate_qfluct(void);
  void propagate_lat(const int n);
};


void print_map(const int n_dof, const string &str, const ss_vect<tps> map)
{
  const int n_dec = 6;

  printf("%s\n", str.c_str());
  for (int i = 0; i < 2*n_dof; i++) {
    for (int j = 0; j < 2*n_dof; j++)
      printf("%*.*e", n_dec+8, n_dec, map[i][j]);
    printf("\n");
  }
}


istream& operator>>(std::istream &is, ss_vect<tps> &a)
{
  int k;

  a.zero();
  for (k = 0; k < 2*nd_tps; k++)
    is >> a[k];
  a[2*nd_tps] = 0e0;
  return is;
}


ss_vect<tps> rd_map(const string &file_name, const string &str)
{
  ss_vect<tps> map;
  ifstream     inf;

  const bool prt = !false;

  file_rd(inf, file_name.c_str());
  inf >> map;
  inf.close();
  if (prt)
    print_map(3, str.c_str(), map);
  return map;
}


ss_vect<tps> tp_map(const int n_dof, const ss_vect<tps> &A)
{
  // Matrix transpose.
  Matrix A_mat;

  getlinmat(2*n_dof, A, A_mat);
  TpMat(2*n_dof, A_mat);
  return putlinmat(2*n_dof, A_mat);
}


void MomentType::rd_maps(void)
{

  const string   file_name = "../compute_maps";

  M        = rd_map(file_name+"_M.dat", "\nM:\n");
  A        = rd_map(file_name+"_A.dat", "\nA:\n");
  M_delta  = rd_map(file_name+"_M_delta.dat", "\nM_delta:\n");
  M_tau    = rd_map(file_name+"_M_tau.dat", "\nM_tau:\n");
  M_Chol   = rd_map(file_name+"_M_Chol.dat", "\nM_Chol:\n");

  M_inv    = Inv(M);
  A_inv    = Inv(A);
  M_Chol_t = tp_map(3, M_Chol);
}


void MomentType::init
(const double Q_b, const double phi_RF, const string &cav_name)
{
  CavityType *C;

  this->cav_loc = Elem_GetPos(ElemIndex(cav_name.c_str()), 1);
  this->Q_b     = Q_b;
  this->phi_RF  = phi_RF*M_PI/180e0;
  this->Circ    = Cell[globval.Cell_nLoc].S;

  file_wr(outf, file_name.c_str());
  rd_maps();

  C = Cell[cav_loc].Elem.C;
  C->phi_RF = this->phi_RF;
}


void MomentType::set_HOM_long
(const double beta, const double f, const double R_sh, const double Q)
{
  Cell[cav_loc].Elem.C->beta_long.push_back(beta);
  Cell[cav_loc].Elem.C->HOM_f_long.push_back(f);
  Cell[cav_loc].Elem.C->HOM_R_sh_long.push_back(R_sh);
  Cell[cav_loc].Elem.C->HOM_Q_long.push_back(Q);
  Cell[cav_loc].Elem.C->HOM_V_long.push_back(0e0);
}


void MomentType::compute_sigma
(const double eps[], const double sigma_s, const double sigma_delta)
{
  int          k;
  ss_vect<tps> Id;

  Id.identity();
  sigma = 0e0;
  for (k = 0; k < 2; k++)
    sigma += eps[k]*(sqr(Id[2*k])+sqr(Id[2*k+1]));
  sigma = sigma*A_inv;
  sigma += sqr(sigma_delta*Id[ct_]) + sqr(sigma_s*Id[delta_]);
}


void MomentType::print_sigma(const int n)
{
  this->n = n;

  outf << scientific << setprecision(3)
       << setw(5) << n
       << setw(11) << -sigma[x] << setw(11) << sigma[p_x]
       << setw(11) << -sigma[y] << setw(11) << sigma[p_y]
       << setw(11) << -sigma[delta] << setw(11) << sigma[ct]/c0 << " |"
       << setw(10) << sqrt(sigma[x_x]) << setw(10) << sqrt(sigma[p_x_p_x])
       << setw(10) << sqrt(sigma[y_y]) << setw(10) << sqrt(sigma[p_y_p_y])
       << setw(10) << sqrt(sigma[delta_delta])
       << setw(10) << sqrt(sigma[ct_ct])/c0
       << setw(11) << abs(Cell[cav_loc].Elem.C->HOM_V_long[0])
       << setw(11) << arg(Cell[cav_loc].Elem.C->HOM_V_long[0]) << "\n";
}


void MomentType::propagate_cav_HOM_transv(const int n)
{
  complex<double> V0;

  const double
    beta_RF = Cell[cav_loc].Elem.C->beta_long[0],
    f       = Cell[cav_loc].Elem.C->HOM_f_long[0],
    R_sh    = Cell[cav_loc].Elem.C->HOM_R_sh_long[0],
    Q       = Cell[cav_loc].Elem.C->HOM_Q_long[0],
    k_loss  = 2e0*M_PI*f*R_sh/Q,                      /* Ring convention
							  P = V^2/(2Rs).  */
    Q_loaded = Q/(1e0+beta_RF);
  const complex<double>
    I = complex<double>(0e0, 1e0);

  V0 = Cell[cav_loc].Elem.C->HOM_V_long[0];

  // Update RF cavity HOM phasor.
  Cell[cav_loc].Elem.C->HOM_V_long[0] *=
    exp(2e0*M_PI*f*(Circ+sigma[ct])/c0*(-1e0/(2e0*Q_loaded)+I));

  // First half increment of HOM phasor.
  Cell[cav_loc].Elem.C->HOM_V_long[0] += Q_b*k_loss/2e0;

  // Propagate through wake field.
  sigma -= Q_b*real(Cell[cav_loc].Elem.C->HOM_V_long[0]);

  // Second half increment of HOM phasor.
  Cell[cav_loc].Elem.C->HOM_V_long[0] += Q_b*k_loss/2e0;
}


ss_vect<tps> MomentType::compute_cav_HOM_long_M(const tps ct)
{
  double          delta;
  complex<double> V0;
  ss_vect<tps>    M_cav;

  const double
    beta_RF  = Cell[cav_loc].Elem.C->beta_long[0],
    f        = Cell[cav_loc].Elem.C->HOM_f_long[0],
    R_sh     = Cell[cav_loc].Elem.C->HOM_R_sh_long[0],
    Q        = Cell[cav_loc].Elem.C->HOM_Q_long[0],
    k_loss   = 2e0*M_PI*f*R_sh/(2e0*Q),                 /* Ring convention
							   P = V^2/(2Rs).  */
    Q_loaded = Q/(1e0+beta_RF);
  const complex<double>
    I = complex<double>(0e0, 1e0);

  M_cav.identity();

  V0 = Cell[cav_loc].Elem.C->HOM_V_long[0];

  // Update RF cavity HOM phasor.
  Cell[cav_loc].Elem.C->HOM_V_long[0] *=
    exp(2e0*M_PI*f*(Circ+ct.cst())/c0*(-1e0/(2e0*Q_loaded)+I));

  // First half increment of HOM phasor: Q_b is negative.
  Cell[cav_loc].Elem.C->HOM_V_long[0] += Q_b*k_loss;

  // Propagate through wake field.
  delta = Q_b*real(Cell[cav_loc].Elem.C->HOM_V_long[0]);
  if (!false)
    delta *= 9.14e5;

  M_cav[delta_] += delta;

  // Second half increment of HOM phasor: Q_b is negative.
  Cell[cav_loc].Elem.C->HOM_V_long[0] += Q_b*k_loss;

  return M_cav;
}


void MomentType::propagate_cav_HOM_long(const int n)
{
  ss_vect<tps> M_cav;
   
  M_cav = compute_cav_HOM_long_M(tps(sigma[ct], ct_+1));
  sigma += M_cav[delta_].cst()*ps_sign[delta_]*tps(0e0, ps_index[delta_]+1);
}


ss_vect<tps> MomentType::compute_M_cav(void)
{
  tps          ct0, delta;
  ss_vect<tps> M_cav;

  const CavityType* C = Cell[cav_loc].Elem.C;

  danot_(1);
  ct0 = tps(sigma[ct], ct_+1);
  M_cav.identity();
#if 1
  M_cav[ct_] += ct0.cst();
  Cav_Pass(Cell[cav_loc], M_cav);
  M_cav[delta_] -= M_cav[delta_].cst();
#else
  delta =
    -C->V_RF/(1e9*globval.Energy)*sin(2e0*M_PI*C->f_RF*ct0/c0+this->phi_RF);
  M_cav[delta_] += delta - delta.cst();
#endif
  danot_(NO);
  return M_cav;
}


void MomentType::propagate_cav(void)
{
  ss_vect<tps> M_cav;
   
  M_cav = compute_M_cav();
  sigma = sigma*Inv(M_cav);
}


void MomentType::propagate_delta(void)
{
  sigma = sigma*M_delta;
}


void MomentType::propagate_tau(void)
{
  sigma = sigma*M_tau;
}


void compute_stochastic_part
(ss_vect<double> &X, ss_vect<tps> &X_map)
{
  // Compute stochastic part of map.
  static std::default_random_engine rand;
  std::normal_distribution<double>  norm_ranf(0e0, 1e0);

  for (int k = 0; k < 2*nd_tps; k++) {
    X[k] = norm_ranf(rand);
    X_map[k] = X[k]*tps(0e0, k+1);
  }
}


void MomentType::propagate_qfluct(void)
{
  int             j, k;
  ss_vect<double> X, qfluct;
  ss_vect<tps>    Id, X_map, qfluct2;

  std::default_random_engine       rand;
  std::normal_distribution<double> norm_ranf(0e0, 1e0);

  Id.identity();
  compute_stochastic_part(X, X_map);
  qfluct = (M_Chol_t*X).cst();
  qfluct2 = M_Chol_t*sqr(X_map)*M_Chol;
  for (j = 0; j < 2*nd_tps; j++) {
    sigma += ps_sign[j]*qfluct[j]*Id[ps_index[j]];
    for (k = 0; k < 2*nd_tps; k++)
      sigma +=
	ps_sign[j]*ps_sign[k]*qfluct2[j][k]*Id[ps_index[j]]*Id[ps_index[k]];
  }
}


void MomentType::propagate_mag_lat(void)
{
  sigma = sigma*M_inv;
}


void MomentType::propagate_lat(const int n)
{
  if (globval.Cavity_on)
    propagate_cav();

  // if (globval.radiation)
  //   propagate_delta();

  if (false)
    propagate_cav_HOM_long(n);

  propagate_mag_lat();

  if (globval.radiation)
    propagate_tau();

  if (globval.radiation)
    propagate_qfluct();
}


void track(MomentType &m, const int n, ss_vect<double> &ps)
{
  long int lastpos;
  int      k;

  ps[ct_] /= c0;
  cout << scientific << setprecision(3)
       << "\n" << setw(5) << 0 << setw(11) << ps << "\n";
  ps[ct_] *= c0;
  for (k = 1; k <= n; k++) {
    if (true)
      Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
    else {
      if (globval.Cavity_on)
	Cav_Pass(Cell[m.cav_loc], ps);
      ps = (m.M*ps).cst();
    }
    ps[ct_] /= c0;
    cout << scientific << setprecision(3)
	 << setw(5) << k << setw(11) << ps << "\n";
    ps[ct_] *= c0;
  }
}


void test_case(const string &cav_name)
{
  // Single bunch, 1 longitudinal HOM.
  int             k;
  ss_vect<double> ps;
  MomentType      m;

  const int
    n_turn      = 20000;
  const double
    eps[]       = {161.7e-12, 8e-12},
    sigma_s     = 3.739e-3,
    sigma_delta = 9.353e-04,

    Q_b         = -0.6e-9,
    phi_RF      = -30.63,
    beta_HOM    = 1e0,
    f           = 1e9,
    R_sh        = 1e3,
#if 0
    Q           = 1e8;
#else
    Q           = 5e3;
#endif  

  globval.pathlength = false;

  danot_(1);

  m.init(Q_b, phi_RF, cav_name);

  danot_(NO);

  m.set_HOM_long(beta_HOM, f, R_sh, Q);

  globval.radiation = true;
  globval.Cavity_on = true;

  ps.zero();
  if (!false) {
    ps[x_]     =  0e-6;
    ps[px_]    =  0e-6;
    ps[y_]     =  0e-6;
    ps[py_]    =  0e-6;
    ps[ct_]    =  0e-6;
    ps[delta_] =  0e-6;
  }

  m.sigma = 0e0;
  if (!false)
    m.compute_sigma(eps, sigma_s, sigma_delta);

  for (k = 0; k < 2*nd_tps; k++)
    m.sigma += ps_sign[k]*ps[k]*tps(0e0, ps_index[k]+1);

  m.print_sigma(0);
  printf("\n");
  for (k = 1; k <= n_turn; k++) {
    m.propagate_lat(k);
    m.print_sigma(k);
  }

  if (false)
    track(m, n_turn, ps);
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

  const string cav_name = "cav";

  reverse_elem     = true;
  globval.mat_meth = false;

  // disable from TPSALib- and LieLib log messages
  idprset_(-1);

  if (false)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  set_state();

  test_case(cav_name);
}
