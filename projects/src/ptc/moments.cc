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
  bool
    quant_fluct, incl_cst;
  int
    n;
  double
    Q_b,         // Bunch charge.
    t_q,         // Time stamp.
    Circ,        // Circumference.
    phi_RF;      // RF phase [rad].
  ss_vect<tps>
    R_no_rad,
    A_CS_no_rad, // Courant-Snyder.
    A_sb_no_rad, // Synchro-betatron.
    R_rad,
    A_CS_rad,    // Courant-Snyder.
    A_sb_rad,    // Synchro-betatron.
    M_tau,
    M_Chol,      // Cholesky decomposition.
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
    sigma;       // Statistical moments for charge distribution.
  ss_vect<double>
    fixed_point_rad,
    M_fixed_point;
  ss_vect<tps>
    M_lat,
    M_lat_inv,
    M_cav,
    M_cav_inv;
  CavityType
    *C;

  void rd_maps(void);
  void init
  (const double Q_b, const double fixed_point[], const double phi_RF,
   const string &cav_name, const bool quant_fluct, const bool incl_cst);
  void set_HOM_long(const double beta, const double f, const double R_sh,
		    const double Q);
  void set_HOM_trans(const double beta, const double f, const double R_sh,
		     const double Q, const int plane);
  void compute_sigma(const double eps[], const double sigma_s,
		     const double sigma_delta);
  void print_sigma(const int n);
  ss_vect<tps> compute_cav_HOM_long_M(const tps &ct);
  ss_vect<tps> compute_cav_HOM_trans_M(const tps &ct);
  void propagate_cav_HOM_long(const int n);
  void propagate_cav_HOM_trans(const int n);
  void compute_M(const bool inverse);
  ss_vect<tps> compute_M_cav(void);
  ss_vect<tps> compute_M_cav_inv(void);
  void compute_M_and_M_inv(void);
  void propagate_qfluct(void);
  void propagate_lat(const int n);
};


double compute_det(const int n_dof, const ss_vect<tps> &A)
{
  // Matrix determinant.
  Matrix A_mat;

  getlinmat(2*n_dof, A, A_mat);
  return DetMat(2*n_dof, A_mat);
}


std::vector<double> get_nu(const ss_vect<tps> &M)
{
  double              nu_z, alpha_c;
  Matrix              A_mat, A_inv_mat, R_mat, M_mat;
  std::vector<double> nu;

  const int
    n_dof = (globval.Cavity_on)? 3 : 2;
  const double
    Circ  = Cell[globval.Cell_nLoc].S;

  getlinmat(2*n_dof, M, M_mat);
  GDiag(2*n_dof, Circ, A_mat, A_inv_mat, R_mat, M_mat, nu_z, alpha_c);
  for (int k =  0; k < n_dof; k++)
    nu.push_back(atan2(globval.wi[2*k], globval.wr[2*k])/(2e0*M_PI));
  return nu;
}


void prt_map(ostream &outf, const string &str, const ss_vect<tps> map)
{
  const int n_dec = 6;

  outf << str << "\n" << "Constant:\n";
  outf << scientific << setprecision(n_dec) << setw(8+n_dec)
       << map.cst() << "\n";
  outf << "Map:\n";
  for (int j = 0; j < 2*nd_tps; j++) {
    for (int k = 0; k < 2*nd_tps; k++)
      outf << scientific << setprecision(n_dec) << setw(8+n_dec) << map[j][k];
    outf << "\n";
  }
}


void prt_map(const string &str, const ss_vect<tps> map)
{
  prt_map(cout, str, map);
}


void prt_map(const string &file_name, const string &str, const ss_vect<tps> map)
{
  ofstream outf;

  const int n_dec = 6;

  file_wr(outf, file_name);
  outf << scientific << setprecision(n_dec) << setw(8+n_dec) << map;
  outf.close();
}


void prt_map_full(const string &str, const long int cav_loc,
		  const ss_vect<tps> &M)
{
  int                 k;
  std::vector<double> nu;

  printf("\n-----------------------------------\n");

  nu = get_nu(M);

  printf("%s", str.c_str());
  printf("Cavity_on    = %s\n", (globval.Cavity_on == true)?"true":"false");
  printf("radiation    = %s\n", (globval.radiation == true)?"true":"false");
  printf("phi_RF [deg] = %4.2f\n",
	 Cell[cav_loc].Elem.C->phi_RF*180e0/M_PI);
  prt_map("", M);
  printf("\ndet(M) - 1:\n %10.3e\n", compute_det(nd_tps, M)-1e0);
  printf("nu:\n");
  for (k = 0; k < nd_tps; k++)
    printf(" %7.5f", nu[k]);
  printf("\n");
  printf("-----------------------------------\n");
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

  const bool prt = false;

  file_rd(inf, file_name.c_str());
  inf >> map;
  inf.close();
  if (prt)
    prt_map(str.c_str(), map);
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


void MomentType::rd_maps(void)
{
  const string   file_name = "../compute_maps";

  R_no_rad    = rd_map(file_name+"_R_no_rad.dat", "\nR_no_rad:");
  A_CS_no_rad = rd_map(file_name+"_A_CS_no_rad.dat", "\nA_CS_no_rad:");
  A_sb_no_rad = rd_map(file_name+"_A_sb_no_rad.dat", "\nA_sb_no_rad:");
  R_rad       = rd_map(file_name+"_R_rad.dat", "\nR_rad:");
  A_CS_rad    = rd_map(file_name+"_A_CS_rad.dat", "\nA_CS_rad:");
  A_sb_rad    = rd_map(file_name+"_A_sb_rad.dat", "\nA_sb_rad:");
  M_tau       = rd_map(file_name+"_M_tau.dat", "\nM_tau:");
  M_Chol      = rd_map(file_name+"_M_Chol.dat", "\nM_Chol:");

  M_Chol_t = compute_transp(3, M_Chol);

  // prt_map("\nM_Chol_t:", M_Chol_t);
}


void MomentType::init
(const double Q_b, const double fixed_point[], const double phi_RF,
 const string &cav_name, const bool quant_fluct, const bool incl_cst)
{
  this->cav_loc = Elem_GetPos(ElemIndex(cav_name.c_str()), 1);
  this->C = Cell[this->cav_loc].Elem.C;

  this->phi_RF = phi_RF*M_PI/180e0;
  C->phi_RF = this->phi_RF;
  this->Q_b = Q_b;
  this->t_q = 0e0;
  this->Circ = Cell[globval.Cell_nLoc].S;

  for (int k =  0; k < 2*nd_tps; k++)
    fixed_point_rad[k] = fixed_point[k];

  this->quant_fluct = quant_fluct;
  this->incl_cst = incl_cst;

  printf("\ninit: RF phase = %5.3f: \n", C->phi_RF*180e0/M_PI);

  file_wr(outf, file_name.c_str());
  rd_maps();
}


void MomentType::set_HOM_long
(const double beta, const double f, const double R_sh, const double Q)
{
  C->beta_long.push_back(beta);
  C->HOM_f_long.push_back(f);
  C->HOM_R_sh_long.push_back(R_sh);
  C->HOM_Q_long.push_back(Q);
  C->HOM_V_long.push_back(0e0);
}


void MomentType::set_HOM_trans
(const double beta, const double f, const double R_sh, const double Q,
 const int plane)
{
  C->beta_trans[plane].push_back(beta);
  C->HOM_f_trans[plane].push_back(f);
  C->HOM_R_sh_trans[plane].push_back(R_sh);
  C->HOM_Q_trans[plane].push_back(Q);
  C->HOM_V_trans[plane].push_back(0e0);
}


void MomentType::compute_sigma
(const double eps[], const double sigma_s, const double sigma_delta)
{
  int          k;
  ss_vect<tps> Id;

  Id.identity();
  for (k = 0; k < 2; k++)
    sigma += eps[k]*(sqr(Id[2*k])+sqr(Id[2*k+1]));
  sigma = sigma*Inv(A_CS_no_rad);
  if (lat_disp)
    sigma = sigma*Inv(A_sb_no_rad);
  sigma +=
    sqr(sigma_s*Id[ps_index[ct_]]) + sqr(sigma_delta*Id[ps_index[delta_]]);
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
       << setw(11) << abs(C->HOM_V_long[0])
       << setw(11) << arg(C->HOM_V_long[0])
       << setw(11) << abs(C->HOM_V_trans[X_][0])
       << setw(11) << arg(C->HOM_V_trans[X_][0]) << "\n";
}


ss_vect<tps> MomentType::compute_cav_HOM_long_M(const tps &ct)
{
  double       delta;
  ss_vect<tps> M_HOM;

  const double
    beta_RF  = C->beta_long[0],
    f        = C->HOM_f_long[0],
    R_sh     = C->HOM_R_sh_long[0],
    Q        = C->HOM_Q_long[0],
    k_loss   = 2e0*M_PI*f*R_sh/(2e0*Q),                 /* Ring convention
							   P = V^2/(2Rs).  */
    Q_loaded = Q/(1e0+beta_RF),

    E0       = 1e9*globval.Energy;

  const complex<double>
    I = complex<double>(0e0, 1e0);

  M_HOM.identity();

  // Update RF cavity HOM phasor.
  C->HOM_V_long[0] *=
    exp(2e0*M_PI*f*(Circ+ct.cst()-t_q)/c0*(-1e0/(2e0*Q_loaded)+I));

  // First half increment of HOM phasor: Q_b is negative.
  C->HOM_V_long[0] += Q_b*k_loss;

  // Propagate through wake field.
  delta = real(C->HOM_V_long[0])/E0;

  M_HOM[delta_] += delta;

  // Second half increment of HOM phasor: Q_b is negative.
  C->HOM_V_long[0] += Q_b*k_loss;

  return M_HOM;
}


void MomentType::propagate_cav_HOM_long(const int n)
{
  ss_vect<tps> M_HOM;
   
  M_HOM = compute_cav_HOM_long_M(tps(sigma[ct], ct_+1));
  sigma += M_HOM[delta_].cst()*ps_sign[delta_]*tps(0e0, ps_index[delta_]+1);
}


ss_vect<tps> MomentType::compute_cav_HOM_trans_M(const tps &ct)
{
  double       delta;
  ss_vect<tps> M_HOM;

  const double
    beta_RF  = C->beta_trans[X_][0],
    f        = C->HOM_f_trans[X_][0],
    R_sh     = C->HOM_R_sh_trans[X_][0],
    Q        = C->HOM_Q_trans[X_][0],
    k_loss   = 2e0*M_PI*f*R_sh/(2e0*Q), /* Ring convention
						    P = V^2/(2Rs).            */
    Q_loaded = Q/(1e0+beta_RF),

    E0       = 1e9*globval.Energy;

  const complex<double>
    I = complex<double>(0e0, 1e0);

  M_HOM.identity();

  // Update RF cavity HOM phasor.
  C->HOM_V_trans[X_][0] *=
    exp(2e0*M_PI*f*(Circ+ct.cst()-t_q)/c0*(-1e0/(2e0*Q_loaded)+I));

  // First half increment of HOM phasor.
  C->HOM_V_trans[X_][0] += ps_sign[x_]*sigma[x]*Q_b*k_loss;

  // Propagate through wake field.
  delta = real(C->HOM_V_long[0])/E0;

  M_HOM[delta_] += delta;

  // Second half increment of HOM phasor.
  C->HOM_V_trans[X_][0] += ps_sign[x_]*sigma[x]*Q_b*k_loss;

  return M_HOM;
}


void MomentType::propagate_cav_HOM_trans(const int n)
{
  ss_vect<tps> M_HOM;
   
  M_HOM = compute_cav_HOM_trans_M(tps(sigma[x], x_+1));
  sigma += M_HOM[px_].cst()*ps_sign[px_]*tps(0e0, ps_index[px_]+1);
}


template<typename T>
void RF_cav_pass(long int cav_loc, const double phi_RF, ss_vect<T> &ps)
{
  const bool Cav_state = globval.Cavity_on;

  globval.Cavity_on = true;

  Cav_Pass(Cell[cav_loc], ps);

  globval.Cavity_on = Cav_state;
}


ss_vect<tps> MomentType::compute_M_cav(void)
{
  const int no = getno_();

  danot_(1);
  M_cav.identity();
  RF_cav_pass(cav_loc, phi_RF, M_cav);
  if (!incl_cst)
    M_cav -= M_cav.cst();
  danot_(no);
  return M_cav;
}


ss_vect<tps> MomentType::compute_M_cav_inv(void)
{
  const int no = getno_();

  danot_(1);
  M_cav_inv = Inv(M_cav-M_cav.cst());
  if (incl_cst)
    M_cav_inv -= M_cav.cst();
  danot_(no);
  return M_cav_inv;
}


void MomentType::compute_M_and_M_inv(void)
{
  long int        lastpos;
  ss_vect<double> fp1;
  ss_vect<tps>    R, A, map;

  const bool validation = false;

  if (!globval.radiation) {
    C->phi_RF = 0e0;
    M_fixed_point.zero();

    M_lat = A_CS_no_rad*R_no_rad*Inv(A_CS_no_rad);
    if (lat_disp)
      M_lat = A_sb_no_rad*M_lat*Inv(A_sb_no_rad);

    M_lat_inv = Inv(M_lat);
  } else {
    C->phi_RF = phi_RF;
    M_fixed_point = fixed_point_rad;

    M_lat = A_CS_rad*M_tau*R_rad*Inv(A_CS_rad);
    if (lat_disp)
      M_lat = A_sb_rad*M_lat*Inv(A_sb_rad);
 
    // Psuedo inverse for transformation theorem.
    M_lat_inv = A_CS_rad*M_tau*Inv(R_rad)*Inv(A_CS_rad);
    if (lat_disp)
      M_lat_inv = A_sb_rad*M_lat_inv*Inv(A_sb_rad);
  }

  prt_map_full("", cav_loc, M_lat);
  prt_map("\nM_lat:", M_lat);
  prt_map("\nM_lat^-1:", M_lat_inv);

  if (validation) {
    map.identity();
    map += M_fixed_point;
    Cell_Pass(0, globval.Cell_nLoc, map, lastpos);
    map -= M_fixed_point;
  
    prt_map("\nValidation:\n", map-M_lat);
  }

  // Remove the RF cavity.
  compute_M_cav();
  compute_M_cav_inv();

  M_lat = M_lat*M_cav_inv;
  M_lat_inv = M_cav*M_lat_inv;

  prt_map("\nM_cav:", M_cav);
  prt_map("\nM_cav^-1:", M_cav_inv);
  prt_map("\nM_lat.M_cav^-1:", M_lat);
  prt_map("\nM_cav.M^-1:", M_lat_inv);

  if (validation) {
    cout << scientific << setprecision(3) << "\nFixed point = "
	 << setw(11) << M_fixed_point << "\n";
    
    fp1 = (M_cav*M_fixed_point).cst();

    cout << scientific << setprecision(3) << "Fixed point = "
	 << setw(11) << fp1 << "\n";
    
    map.identity();
    map += fp1;
    Cell_Pass(2, globval.Cell_nLoc, map, lastpos);
    map -= map.cst();

    prt_map("\nValidation:\n", map-M_lat);
  }
}


void compute_stochastic_part(ss_vect<tps> &X_map)
{
  // Compute stochastic part of map.
  ss_vect<double>                   X;
  static std::default_random_engine rand;
  std::normal_distribution<double>  norm_ranf(0e0, 1e0);

  for (int k = 0; k < 2*nd_tps; k++) {
    X[k] = norm_ranf(rand);
    X_map[k] = X[k]*tps(0e0, k+1);
  }
}


void MomentType::propagate_qfluct(void)
{
  ss_vect<tps> Id, X_map, qfluct;

  Id.identity();
  compute_stochastic_part(X_map);
  qfluct = M_Chol_t*sqr(X_map)*M_Chol;
  for (int j = 0; j < 2*nd_tps; j++) {
    for (int k = 0; k < 2*nd_tps; k++)
      sigma +=
	ps_sign[j]*ps_sign[k]*qfluct[j][k]*Id[ps_index[j]]*Id[ps_index[k]];
  }
}


void MomentType::propagate_lat(const int n)
{
  // Moments are propageted by the inverse map.
  // Hence, the maps for the lattice & RF cavity appears in reversed order.

  double dct;

  sigma = sigma*M_lat_inv;

  if (globval.Cavity_on)
    sigma = sigma*M_cav_inv;

  propagate_cav_HOM_long(n);
  propagate_cav_HOM_trans(n);
  t_q = sigma[ct];

  if (globval.radiation & quant_fluct)
    propagate_qfluct();
}


void track(MomentType &m, const int n, ss_vect<double> &ps)
{
  long int        lastpos;
  double          dct;
  ss_vect<double> ps1;
  ofstream        outf;

  const bool   rel_to_fp = true;
  const string file_name = "track.out";

  file_wr(outf, file_name.c_str());

  ps1 = ps;
  outf << scientific << setprecision(5)
       << "\n" << setw(5) << 0 << setw(13) << ps1 << "\n";

  for (int k = 1; k <= n; k++) {
    if (rel_to_fp)
      ps1 += globval.CODvect;
    Cell_Pass(0, globval.Cell_nLoc, ps1, lastpos);
    if (rel_to_fp)
      ps1 -= globval.CODvect;
    if (lastpos != globval.Cell_nLoc)
      printf("\n particle lost: lastpos = %ld\n", lastpos);

    outf << scientific << setprecision(5)
	 << setw(5) << k << setw(13) << ps1 << "\n";
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
  bool            quant_fluct;
  long int        lastpos;
  int             k;
  double          phi_RF;
  ss_vect<double> bc, ps;
  MomentType      m;

  const int
    n_turn           = 3000;
  const double
    // Observation piont: long straight.
    fp_rad[] =
    {2.448e-09, 1.372e-07, 0.000e+00, 0.000e+00, -1.036e-04, -3.211e-16},
    eps[]            = {161.7e-12, 8e-12},
    sigma_s          = 3.739e-3,
    sigma_delta      = 9.353e-04,

    Q_b              = -0.6e-9,
    phi_RF_rad       = 30.62567,

    beta_HOM_long    = 1e0,
    f_long           = 1e9,
    R_sh_long        = (false)? 1e3 : 0e0,
#if 1
    Q_long           = 1e8,
#else
    Q_long           = 5e3,
#endif  

    beta_HOM_trans_x = 1e0,
    f_trans_x        = 1e9,
    R_sh_trans_x     = (false)? 1e3 : 0e0,
    Q_trans_x        = 1e8;

  danot_(1);

  globval.radiation = !true;
  quant_fluct       = true;
  globval.Cavity_on = true;
  m.lat_disp        = true;

  phi_RF = (!globval.radiation)? 0e0 : phi_RF_rad;

  m.init(Q_b, fp_rad, phi_RF, cav_name, quant_fluct, true);

  m.set_HOM_long(beta_HOM_long, f_long, R_sh_long, Q_long);
  m.set_HOM_trans(beta_HOM_trans_x, f_trans_x, R_sh_trans_x, Q_trans_x, X_);

  m.compute_M_and_M_inv();

  danot_(NO);

  m.sigma = 0e0;
  if (false) {
    m.sigma += ps_sign[x_]*10e-6*tps(0e0, ps_index[x_]+1);
    // m.sigma += ps_sign[y_]*2e-6*tps(0e0, ps_index[y_]+1);
    // m.sigma += ps_sign[delta_]*3e-6*tps(0e0, ps_index[delta_]+1);
  }
  if (false)
    m.compute_sigma(eps, sigma_s, sigma_delta);

  bc = get_barycentre(m.sigma);
  printf("\nBarycentre = [");
  for (k = 0; k < 2*nd_tps; k++)
    printf("%9.3e%s", bc[k], (k < 2*nd_tps-1)? ", " : "]\n");

  if (!false) {
    m.print_sigma(0);
    for (k = 1; k <= n_turn; k++) {
      m.propagate_lat(k);
      m.print_sigma(k);
    }
  }

  if (false) {
    printf("\ntrack:\n");
    printf("Cavity_on = %s\n", (globval.Cavity_on == true)?"true":"false");
    printf("radiation = %s\n", (globval.radiation == true)?"true":"false");

    ps = bc;

    getcod(0e0, lastpos);
    printf("\nFixed Point = [");
    for (k = 0; k < 2*nd_tps; k++)
      printf("%9.3e%s", globval.CODvect[k], (k < 2*nd_tps-1)? ", " : "]\n");
    printf("\n");

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

  test_case(cav_name);
}
