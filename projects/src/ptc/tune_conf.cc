#define NO 6

#include "tracy_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;


const bool
  mpole_zero = !false;

const double
#if 0
  beta_inj[] = {2.7, 3.5},
#else
  beta_inj[] = {2.8, 2.8},
#endif

  A_max[]    = {3e-3, 1.5e-3},
  delta_max  = 4e-2,
  twoJ[]     = {sqr(A_max[X_])/beta_inj[X_], sqr(A_max[Y_])/beta_inj[Y_]},

  scl_ksi[]  = {0e0, 1e2, 1e0, 1e0},
  scl_a[]    = {1e0};


class Lie_object {
private:
  std::vector<std::string>
    prm_name,
    K_label;
public:
  int
    prm_n;
  std::vector<int>
    prm_ord,
    prm_Fnum;
  std::vector<double>
    K_val;
  ss_vect<tps>
    Id_scl;

  void add_prm(const string &name, const int n);
  void prt_Lie_object(void);
  void get_term(const tps &h, const int i, const int j, const int k,
		const int l, const int m);
  void set_prms(const double *b_n);
  void print_term(const int k);
  void print_terms(void);
  void print_prm(const int k);
  void print_prms(void);
};


void Lie_object::add_prm(const string &name, const int n)
{
  prm_name.push_back(name);
  prm_n = prm_name.size();
  prm_Fnum.push_back(ElemIndex(name.c_str()));
  prm_ord.push_back(n);
}


void Lie_object::get_term
(const tps &h, const int i, const int j, const int k, const int l, const int m)
{
  std::ostringstream label;

  label << "k_" << i << j << k << l << m;
  K_label.push_back(label.str());
  K_val.push_back(h_ijklm(h, i, j, k, l, m));
}


void Lie_object::set_prms(const double *b_n)
{
  int k;

  for (k = 0; k < prm_n; k++)
    set_bn_design_fam(prm_Fnum[k], prm_ord[k], b_n[k+1], 0e0);
}


void Lie_object::print_term(const int k)
{
  printf("  %-7s %10.3e\n", K_label[k].c_str(), K_val[k]);
}


void Lie_object::print_terms(void)
{
  int k;

  printf("\nLinear chromaticity:\n");
  for (k = 0; k < 2; k++)
    print_term(k);

  printf("\nAnharmonic terms:\n");
  for (k = 2; k < 5; k++)
    print_term(k);

  printf("\nAnharmonic terms:\n");
  for (k = 5; k < 9; k++)
    print_term(k);

  printf("\n2nd order chromaticity:\n");
  for (k = 9; k < 11; k++)
    print_term(k);

  printf("\n3rd order chromaticity:\n");
  for (k = 11; k < 13; k++)
    print_term(k);

  printf("\n4th order chromaticity:\n");
  for (k = 13; k < 15; k++)
    print_term(k);

  printf("\n5th order chromaticity:\n");
  for (k = 15; k < 17; k++)
    print_term(k);
}


void Lie_object::print_prms(void)
{
  int    k;
  double b_n, a_n;

  printf("\nParameters:\n");
  for (k = 0; k < prm_n; k++) {
    get_bn_design_elem(prm_Fnum[k], 1, prm_ord[k], b_n, a_n);
    printf("  %-5s %1d %10.3e\n",
	   prm_name[k].c_str(), prm_ord[k], b_n);
  }
}


void no_mpoles(const int n)
{
  int k;

  printf("\nno_mpoles: zeroing %d\n", n);
  for (k = 0; k < globval.Cell_nLoc; k++)
    if (Cell[k].Elem.Pkind == Mpole)
      set_bn_design_elem(Cell[k].Fnum, Cell[k].Knum, n, 0e0, 0e0);
}


// Global object, to pass parameters to minimisation function.
Lie_object Lie;


void compute_Lie_terms(Lie_object &Lie, tps &K_re)
{
  Lie.K_val.clear();

  Lie.get_term(K_re, 1, 1, 0, 0, 1);
  Lie.get_term(K_re, 0, 0, 1, 1, 1);

  Lie.get_term(K_re, 2, 2, 0, 0, 0);
  Lie.get_term(K_re, 1, 1, 1, 1, 0);
  Lie.get_term(K_re, 0, 0, 2, 2, 0);

  Lie.get_term(K_re, 3, 3, 0, 0, 0);
  Lie.get_term(K_re, 2, 2, 1, 1, 0);
  Lie.get_term(K_re, 1, 1, 2, 2, 0);
  Lie.get_term(K_re, 0, 0, 3, 3, 0);

  Lie.get_term(K_re, 1, 1, 0, 0, 2);
  Lie.get_term(K_re, 0, 0, 1, 1, 2);

  Lie.get_term(K_re, 1, 1, 0, 0, 3);
  Lie.get_term(K_re, 0, 0, 1, 1, 3);

  Lie.get_term(K_re, 1, 1, 0, 0, 4);
  Lie.get_term(K_re, 0, 0, 1, 1, 4);

  Lie.get_term(K_re, 1, 1, 0, 0, 5);
  Lie.get_term(K_re, 0, 0, 1, 1, 5);
}


double compute_chi2(Lie_object &Lie)
{
  double chi2 = 0e0;

  chi2 += scl_ksi[1]*(sqr(Lie.K_val[0])+sqr(Lie.K_val[1]));

  chi2 +=
    scl_a[0]*
    (sqr(Lie.K_val[2]+Lie.K_val[5])
     +sqr(Lie.K_val[3]+Lie.K_val[6]+Lie.K_val[7])
     +sqr(Lie.K_val[4]+Lie.K_val[8]));

  chi2 +=
    scl_ksi[2]*
    (sqr(Lie.K_val[9]+Lie.K_val[13])+sqr(Lie.K_val[10]+Lie.K_val[14]));

  chi2 +=
    scl_ksi[3]*
    (sqr(Lie.K_val[11]+Lie.K_val[15])+sqr(Lie.K_val[12]+Lie.K_val[16]));

  return chi2;
}


double f_dnu(double *b_n)
{
  static double chi2_ref = 1e30;
  double        chi2;
  tps           g_re, g_im, K, K_re, K_im;
  ss_vect<tps>  nus;

  Lie.set_prms(b_n);

  danot_(no_tps-1);
  get_map(false);
  danot_(no_tps);
  MNF = MapNorm(map, 1);

  CtoR(MNF.g*Lie.Id_scl, g_re, g_im);
  CtoR(MNF.K*Lie.Id_scl, K_re, K_im);

  compute_Lie_terms(Lie, K_re);

  chi2 = compute_chi2(Lie);

  if (chi2 < chi2_ref) {
    printf("\nchi2 = %10.3e\n", chi2);

    chi2_ref = chi2;
    Lie.print_prms();
    Lie.print_terms();

    prtmfile("flat_file.fit");
  }

  return chi2;
}


void get_prms(Lie_object &Lie)
{
  if (!false) {
    Lie.add_prm("sf",  Sext);
    Lie.add_prm("sd",  Sext);
    // Lie.add_prm("sd2", Sext);
  }

  if (!false) {
    Lie.add_prm("uq1",  Oct);
    Lie.add_prm("uq2",  Oct);
    Lie.add_prm("uq3", Oct);
    Lie.add_prm("uq4", Oct);
  }

  if (!false) {
    Lie.add_prm("sf",  Oct);
    Lie.add_prm("sd",  Oct);
    // Lie.add_prm("sd2", Oct);
  }

  if (!false) {
    Lie.add_prm("sf",  Dec);
    Lie.add_prm("sd",  Dec);
    // Lie.add_prm("sd2", Dec);
  }
}


void fit_powell(void)
{
  int          j, k, iter;
  double       *b_n, **xi, fret, b, a;

  const double eps = 1e-1, ftol = 1e-6;

  // Set parameters.
  get_prms(Lie);

  b_n = dvector(1, Lie.prm_n);
  xi = dmatrix(1, Lie.prm_n, 1, Lie.prm_n);

  if (mpole_zero) {
    no_mpoles(Sext);
    no_mpoles(Oct);
  }

  // Initialise parameters.
  for (k = 0; k < Lie.prm_n; k++) {
    get_bn_design_elem(Lie.prm_Fnum[k], 1, Lie.prm_ord[k], b, a);
    b_n[k+1] = b;
  }

  Lie.Id_scl.identity();
  for (k = 0; k < 4; k++)
    Lie.Id_scl[k] *= sqrt(twoJ[k/2]);
  Lie.Id_scl[delta_] *= delta_max;

  f_dnu(b_n);

  // Set initial directions (unit vectors).
  for (j = 1; j <= Lie.prm_n; j++)
    for (k = 1; k <= Lie.prm_n; k++)
      xi[j][k] = (j == k)? eps : 0e0;

  dpowell(b_n, xi, Lie.prm_n, ftol, &iter, &fret, f_dnu);

  printf("\n  iter = %d fret = %12.5e\n", iter, fret);
  printf("bns:\n");

  free_dvector(b_n, 1, Lie.prm_n);
  free_dmatrix(xi, 1, Lie.prm_n, 1, Lie.prm_n);
}


void set_lat_state()
{
  globval.H_exact        = false;
  globval.quad_fringe    = false;
  globval.Cavity_on      = false;
  globval.radiation      = false;
  globval.emittance      = false;
  globval.IBS            = false;
  globval.pathlength     = false;
  globval.bpm            = 0;
  globval.Cart_Bend      = false;
  globval.dip_edge_fudge = true;
  globval.mat_meth       = false;
}


int main(int argc, char *argv[])
{

  set_lat_state();

  // disable from TPSALib and LieLib log messages
  idprset(-1);

  std::string home_dir = "";

  daeps_(eps_tps);

  if (!true)
    Read_Lattice((home_dir+argv[1]).c_str());
  else
    rdmfile(argv[1]);

  fit_powell();
}
