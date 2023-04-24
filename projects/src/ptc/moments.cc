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
    t_q,       // Time stamp.
    Circ,      // Circumference.
    phi_RF,    // RF phase [rad].
    delta_RF;  // Energy change.
  ss_vect<tps>
    R,
    M_tau,
    A_CS,      // Courant-Snyder.
    A_sb,      // Synchro-betatron.
    M_Chol,    // Cholesky decomposition.
    M_Chol_t;
  ofstream
    outf;
public:
  bool
    lat_disp,
    RF_cav_linear;
  long int
    cav_loc;
  tps
    sigma;     // Statistical moments for charge distribution.
  ss_vect<double>
    fixed_point;
  ss_vect<tps>
    M,
    M_inv,
    M_cav;

  void rd_maps(void);
  void init(const double Q_b, const double phi_RF, const string &cav_name);
  void set_HOM_long(const double beta, const double f, const double R_sh,
		    const double Q);
  void compute_sigma(const double eps[], const double sigma_s,
		     const double sigma_delta);
  void print_sigma(const int n);
  ss_vect<tps> compute_cav_HOM_long_M(const tps &ct);
  ss_vect<tps> compute_cav_HOM_trans_M(const tps &ct);
  void propagate_cav_HOM_long(const int n);
  void propagate_cav_HOM_trans(const int n);
  void compute_M_cav(void);
  void compute_M_and_M_inv(const double fp[]);
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


ss_vect<tps> compute_transp(const int n_dof, const ss_vect<tps> &A)
{
  // Matrix transpose.
  Matrix A_mat;

  getlinmat(2*n_dof, A, A_mat);
  TpMat(2*n_dof, A_mat);
  return putlinmat(2*n_dof, A_mat);
}


double compute_det(const int n_dof, const ss_vect<tps> &A)
{
  // Matrix determinant.
  Matrix A_mat;

  getlinmat(2*n_dof, A, A_mat);
  return DetMat(2*n_dof, A_mat);
}


void MomentType::rd_maps(void)
{

  const string   file_name = "../compute_maps";

  R       = rd_map(file_name+"_R.dat", "\nR:");
  A_CS    = rd_map(file_name+"_A_CS.dat", "\nA_CS:");
  A_sb    = rd_map(file_name+"_A_sb.dat", "\nA_sb:");
  M_tau   = rd_map(file_name+"_M_tau.dat", "\nM_tau:");
  M_Chol  = rd_map(file_name+"_M_Chol.dat", "\nM_Chol:");

  M_Chol_t = compute_transp(3, M_Chol);
  print_map(3, "\nM_Chol_t:", M_Chol_t);
}


void MomentType::init
(const double Q_b, const double phi_RF, const string &cav_name)
{
  CavityType *C;

  this->cav_loc = Elem_GetPos(ElemIndex(cav_name.c_str()), 1);
  this->Q_b     = Q_b;
  this->t_q     = 0e0;
  this->phi_RF  = phi_RF*M_PI/180e0;
  this->Circ    = Cell[globval.Cell_nLoc].S;

  file_wr(outf, file_name.c_str());
  rd_maps();

  C = Cell[cav_loc].Elem.C;
  C->phi_RF = this->phi_RF;
  printf("\ninit: RF phase = %5.3f: \n", C->phi_RF*180e0/M_PI);
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
  for (k = 0; k < 2; k++)
    sigma += eps[k]*(sqr(Id[2*k])+sqr(Id[2*k+1]));
  sigma = sigma*Inv(A_CS);
  if (lat_disp)
    sigma = sigma*Inv(A_sb);
  sigma +=
    sqr(sigma_s*Id[ct_]) + sqr(sigma_delta*Id[delta_]);
}


void MomentType::print_sigma(const int n)
{
  this->n = n;

  outf << scientific << setprecision(3)
       << setw(5) << n
       << setw(11) << ps_sign[x_]*sigma[x]
       << setw(11) << ps_sign[px_]*sigma[p_x]
       << setw(11) << ps_sign[y_]*sigma[y]
       << setw(11) << ps_sign[py_]*sigma[p_y]
       << setw(11) << ps_sign[delta_]*sigma[delta]
       << setw(11) << ps_sign[ct_]*sigma[ct]
       << " |"
       << setw(10) << sqrt(sigma[x_x]) << setw(10) << sqrt(sigma[p_x_p_x])
       << setw(10) << sqrt(sigma[y_y]) << setw(10) << sqrt(sigma[p_y_p_y])
       << setw(10) << sqrt(sigma[delta_delta])
       << setw(10) << sqrt(sigma[ct_ct])
       << setw(11) << abs(Cell[cav_loc].Elem.C->HOM_V_long[0])
       << setw(11) << arg(Cell[cav_loc].Elem.C->HOM_V_long[0]) << "\n";
}


ss_vect<tps> MomentType::compute_cav_HOM_long_M(const tps &ct)
{
  double       delta;
  ss_vect<tps> M_cav;

  const double
    beta_RF  = Cell[cav_loc].Elem.C->beta_long[0],
    f        = Cell[cav_loc].Elem.C->HOM_f_long[0],
    R_sh     = Cell[cav_loc].Elem.C->HOM_R_sh_long[0],
    Q        = Cell[cav_loc].Elem.C->HOM_Q_long[0],
    k_loss   = 2e0*M_PI*f*R_sh/(2e0*Q),                 /* Ring convention
							   P = V^2/(2Rs).  */
    Q_loaded = Q/(1e0+beta_RF),

    E0       = 1e9*globval.Energy;

  const complex<double>
    I = complex<double>(0e0, 1e0);

  M_cav.identity();

  // Update RF cavity HOM phasor.
  Cell[cav_loc].Elem.C->HOM_V_long[0] *=
    exp(2e0*M_PI*f*(Circ+ct.cst()-t_q)/c0*(-1e0/(2e0*Q_loaded)+I));

  // First half increment of HOM phasor: Q_b is negative.
  Cell[cav_loc].Elem.C->HOM_V_long[0] += Q_b*k_loss;

  // Propagate through wake field.
  delta = real(Cell[cav_loc].Elem.C->HOM_V_long[0])/E0;

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


ss_vect<tps> MomentType::compute_cav_HOM_trans_M(const tps &ct)
{
  double       delta;
  ss_vect<tps> M_cav;

  const double
    beta_RF = Cell[cav_loc].Elem.C->beta_trans[0],
    f       = Cell[cav_loc].Elem.C->HOM_f_trans[0],
    R_sh    = Cell[cav_loc].Elem.C->HOM_R_sh_trans[0],
    Q       = Cell[cav_loc].Elem.C->HOM_Q_trans[0],
    k_loss  = 2e0*M_PI*f*R_sh/(2e0*Q),                 /* Ring convention
							  P = V^2/(2Rs).  */
    Q_loaded = Q/(1e0+beta_RF),

    E0       = 1e9*globval.Energy;

  const complex<double>
    I = complex<double>(0e0, 1e0);

  M_cav.identity();

  // Update RF cavity HOM phasor.
  Cell[cav_loc].Elem.C->HOM_V_trans[0] *=
    exp(2e0*M_PI*f*(Circ+ct.cst()-t_q)/c0*(-1e0/(2e0*Q_loaded)+I));

  // First half increment of HOM phasor.
  // Cell[cav_loc].Elem.C->HOM_V_trans[0] += [x, y]*Q_b*k_loss;

  // Propagate through wake field.
  delta = real(Cell[cav_loc].Elem.C->HOM_V_long[0])/E0;

  M_cav[delta_] += delta;

  // Second half increment of HOM phasor.
  // Cell[cav_loc].Elem.C->HOM_V_trans[0] += [x, y]*Q_b*k_loss;

  return M_cav;
}


void MomentType::propagate_cav_HOM_trans(const int n)
{
  ss_vect<tps> M_cav;
   
  M_cav = compute_cav_HOM_trans_M(tps(sigma[ct], ct_+1));
  sigma += M_cav[delta_].cst()*ps_sign[delta_]*tps(0e0, ps_index[delta_]+1);
}


template<typename T>
void RF_cav_pass(long int cav_loc, ss_vect<T> &ps)
{
  bool Cav_state;

  const CavityType* C = Cell[cav_loc].Elem.C;

#if 1
  Cav_state = globval.Cavity_on;
  globval.Cavity_on = true;
  Cav_Pass(Cell[cav_loc], ps);
  globval.Cavity_on = Cav_state;
#else
  ps[delta_] +=
    -C->V_RF/(1e9*globval.Energy)*sin(2e0*M_PI*C->f_RF*ct0/c0+this->phi_RF);
#endif
}


void MomentType::compute_M_cav(void)
{
  const double dct = ps_sign[ct_]*sigma[ct];

  danot_(1);

  M_cav.identity();
  M_cav[ct_] += dct;
  RF_cav_pass(cav_loc, M_cav);
  M_cav[ct_] -= dct;

  danot_(NO);
}


void MomentType::compute_M_and_M_inv(const double fp[])
{
  // Remark: Pseudo M^-1 for moment tracking.
  int k;
  
  danot_(1);

  for (k = 0; k < 2*nd_tps; k++)
    fixed_point[k] = fp[k];

  M = R;
  if (globval.radiation)
    M = M_tau*M;
  M = A_CS*M*Inv(A_CS);
  if (lat_disp)
    M = A_sb*M*Inv(A_sb);
  // Remove the RF cavity.
  M_cav.identity();
  RF_cav_pass(cav_loc, M_cav);
  cout << scientific << setprecision(3) << "\nM_cav constant part :\n"
       << setw(11) << M_cav.cst() << "\n";
  print_map(3, "\nM_cav:", M_cav);
  cout << scientific << setprecision(3) << "\nM_cav^-1 constant part :\n"
       << setw(11) << (Inv(M_cav)).cst() << "\n";
  print_map(3, "\nM_cav^-1:", Inv(M_cav));
  M = M*Inv(M_cav-M_cav.cst()) ;
    
  M_inv = Inv(R);
  if (globval.radiation)
    M_inv = M_tau*M_inv;
  M_inv = A_CS*M_inv*Inv(A_CS);
  if (lat_disp)
    M_inv = A_sb*M_inv*Inv(A_sb);
  M_inv = (M_cav-M_cav.cst())*M_inv;

  cout << scientific << setprecision(3) << "\nM constant part :\n"
       << setw(11) << M.cst() << "\n";
  print_map(3, "\nM:", M);
  printf("\ndet(M) - 1 = %11.5e\n", compute_det(nd_tps, M)-1e0);
  cout << scientific << setprecision(3)
       << "\nM^-1 constant part :\n" << setw(11) << M_inv.cst() << "\n";
  print_map(3, "\nM_inv:", M_inv);
  printf("\ndet(M_inv) - 1 = %11.5e\n", compute_det(nd_tps, M_inv)-1e0);

  danot_(NO);
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


void MomentType::propagate_lat(const int n)
{
  // Remark: The RF cavity is included - or not - in M_inv;
  double dct;

  if (globval.Cavity_on) {
    if (!RF_cav_linear) {
      compute_M_cav();
      dct = ps_sign[ct_]*sigma[ct];
      sigma -= ps_sign[ct_]*dct*tps(0e0, ps_index[ct_]+1);
      sigma = sigma*Inv(M_cav-M_cav.cst());
      sigma += ps_sign[ct_]*dct*tps(0e0, ps_index[ct_]+1);
    } else
      sigma = sigma*Inv(M_cav-M_cav.cst());
    // Include contribution from constant term.
    sigma += ps_sign[delta_]*M_cav[delta_].cst()*tps(0e0, ps_index[delta_]+1);
  }

  propagate_cav_HOM_long(n);
  // propagate_cav_HOM_trans(n);
  t_q = sigma[ct];

  sigma = sigma*M_inv;

  if (globval.radiation)
    propagate_qfluct();
}


void track(MomentType &m, const int n, ss_vect<double> &ps)
{
  long int lastpos;
  int      k;
  double   dct;
  ofstream outf;

  const string file_name = "track.out";

  file_wr(outf, file_name.c_str());

  ps[ct_] /= c0;
  outf << scientific << setprecision(5)
       << "\n" << setw(5) << 0 << setw(13) << ps << "\n";
  ps[ct_] *= c0;
  for (k = 1; k <= n; k++) {
#if 1
    if (globval.Cavity_on) {
      if (m.RF_cav_linear)
	ps = (m.M_cav*ps).cst();
      else {
	dct = ps[ct_];

	danot_(1);
	m.M_cav.identity();
	m.M_cav[ct_] += dct;
	RF_cav_pass(m.cav_loc, m.M_cav);
	m.M_cav[ct_] -= dct;
	danot_(NO);

	ps[ct_] -= dct;
	ps = (m.M_cav*ps).cst();
	ps[ct_] += dct;
      }
    }
    ps = (m.M*ps).cst();
#else
    // Absolute coordinates.
    Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
#endif
    ps[ct_] /= c0;
    outf << scientific << setprecision(5)
	 << setw(5) << k << setw(13) << ps << "\n";
    ps[ct_] *= c0;
  }

  outf.close();
}


ss_vect<double> get_barycentre(tps &sigma)
{
  int             k;
  ss_vect<double> bc;

  const double barycentre[] =
    {ps_sign[x_]*sigma[x],
     ps_sign[px_]*sigma[p_x],
     ps_sign[y_]*sigma[y],
     ps_sign[py_]*sigma[p_y],
     ps_sign[delta_]*sigma[delta],
     ps_sign[ct_]*sigma[ct]};

  for (k = 0; k < 2*nd_tps; k++)
    bc[k] = barycentre[k];

  return bc;
}


void test_case(const string &cav_name)
{
  // Single bunch, 1 longitudinal HOM.
  int             k;
  ss_vect<double> bc, ps;
  MomentType      m;

  const int
    n_turn        = 40000;
  const double
    fp[] =
    {-4.473e-08, 2.824e-08, 0.000e+00, 0.000e+00, -1.033e-04, -1.703e-10},
    eps[]         = {161.7e-12, 8e-12},
    sigma_s       = 3.739e-3,
    sigma_delta   = 9.353e-04,

    Q_b           = -0*0.6e-9,
    phi_RF        = -30.63,
    delta_RF      = 2.067e-04,
    beta_HOM      = 1e0,
    f             = 1e9,
    R_sh          = 1e3,
#if 0
    Q             = 1e8;
#else
    Q             = 5e3;
#endif  

  m.init(Q_b, phi_RF, cav_name);

  m.set_HOM_long(beta_HOM, f, R_sh, Q);

  globval.radiation = true;
  globval.Cavity_on = true;
  m.lat_disp        = true;
  m.RF_cav_linear   = !true;

  m.compute_M_and_M_inv(fp);

  m.sigma = 0e0;
  if (!false)
    m.compute_sigma(eps, sigma_s, sigma_delta);
  if (false)
    for (k = 0; k < 2*nd_tps; k++)
      m.sigma += ps_sign[k]*m.fixed_point[k]*tps(0e0, ps_index[k]+1);

  bc = get_barycentre(m.sigma);
  bc[ct_] /= c0;
  printf("\nBarycentre = [");
  for (k = 0; k < 2*nd_tps; k++)
    printf("%9.3e%s", bc[k], (k < 2*nd_tps-1)? ", " : "]\n");
  bc[ct_] *= c0;

  m.print_sigma(0);
  printf("\n");
  for (k = 1; k <= n_turn; k++) {
    m.propagate_lat(k);
    m.print_sigma(k);
  }

  if (false) {
    ps.zero();
    if (!false)
      ps = bc;
    track(m, n_turn, ps);
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

  danot_(1);
  test_case(cav_name);
}
