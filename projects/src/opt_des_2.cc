#define NO 1

#include "tracy_lib.h"

#include "Powell/src/newuoa.h"

int no_tps = NO;


const bool
  ps_rot        = !false, // Note, needs to be zeroed; after use.
  phi_spec_case = false,
  qf6_rb        = false,
  sp_short      = !true,
  sp_std        = true,
  relaxed       = true,
  pert_dip_cell = !false;
 

/* opt_case:
     opt_mI_std  1,
     match_ls    2,
     opt_mi_sp   3.                                                           */
const int opt_case = 3;

// From Center of Mid Straight: alpha, beta, eta, eta'.
const int    n_ic        = 4;
const double ic[n_ic][2] =
    {{-0.0000000000, -0.0000000000}, {2.7030084584, 1.9407337499},
     {0.0162153257, 0.0000000000}, {0.0, 0.0}};

#define LAT_CASE 6

const double
#if LAT_CASE == 1
  eps0_x             = 0.079,
  high_ord_achr_nu[] = {21.0/8.0, 7.0/8.0},
  twoJ[]             = {sqr(3e-3)/10.0, sqr(2e-3)/4.0},
  delta              = 2e-2,
#elif LAT_CASE == 2
  eps0_x             = 0.079,
  high_ord_achr_nu[] = {11.0/4.0+0.01, 3.0/4.0-0.01},
  twoJ[]             = {sqr(3e-3)/10.0, sqr(2e-3)/4.0},
  delta              = 2e-2,
#elif LAT_CASE == 3
  eps0_x             = 0.079,
  high_ord_achr_nu[] = {21.0/8.0+0.01, 3.0/4.0-0.01},
  twoJ[]             = {sqr(3e-3)/10.0, sqr(2e-3)/4.0},
  delta              = 2e-2,
#elif LAT_CASE == 4
  eps0_x             = 0.099,
  high_ord_achr_nu[] = {19.0/8.0+0.01, 6.0/8.0-0.01},
  twoJ[]             = {sqr(3e-3)/10.0, sqr(2e-3)/4.0},
  delta              = 2e-2,
#elif LAT_CASE == 5
  eps0_x             = 0.099,
  high_ord_achr_nu[] = {21.0/8.0+0.01, 6.0/8.0-0.01},
  twoJ[]             = {sqr(3e-3)/10.0, sqr(2e-3)/4.0},
  delta              = 2e-2,
#elif LAT_CASE == 6
  eps0_x             = 0.099,
  high_ord_achr_nu[] = {21.0/8.0+0.01, 7.0/8.0-0.01},
  twoJ[]             = {sqr(3e-3)/10.0, sqr(2e-3)/4.0},
  delta              = 2e-2,
#elif LAT_CASE == 7
  eps0_x             = 0.149,
  high_ord_achr_nu[] = {9.0/4.0+0.01, 3.0/4.0-0.01},
  twoJ[]             = {sqr(7e-3)/10.0, sqr(4e-3)/4.0},
  delta              = 2.5e-2,
#elif LAT_CASE == 8
  eps0_x             = 0.149,
  high_ord_achr_nu[] = {19.0/8.0+0.01, 3.0/4.0-0.01},
  twoJ[]             = {sqr(7e-3)/10.0, sqr(4e-3)/4.0},
  delta              = 2.5e-2,
#elif LAT_CASE == 9
  eps0_x             = 0.149,
  high_ord_achr_nu[] = {9.0/4.0+0.01, 5.0/4.0-0.01},
  twoJ[]             = {sqr(7e-3)/10.0, sqr(4e-3)/4.0},
  delta              = 2.5e-2,
#endif
  mI_dnu[]           = {0.0, 0.0},
  mI_nu_ref[]        = {1.5-mI_dnu[X_], 0.5-mI_dnu[Y_]};

double rad2deg(const double a) { return a*180e0/M_PI; }

double deg2rad(const double a) { return a*M_PI/180e0; }


double eps_x;


double bn_internal(const double bn_bounded,
		   const double bn_min, const double bn_max);
double bn_bounded(const double bn_internal,
		  const double bn_min, const double bn_max);
void set_ds(const int Fnum, const double ds);
void set_phi(const int Fnum, const double phi);
void set_dphi(const int Fnum, const double dphi);
void get_S(void);


struct param_type {
  // Parameter Type:
  //   Bend Angle  -3,
  //   Length      -2,
  //   Position    -1,
  //   Quadrupole   2.
private:

public:
  int
    n_prm;
  double
    bn_tol,
    step;
  std::vector<double>
    bn_min,
    bn_max,
    bn_scl,
    l0;
  std::vector<int>
    Fnum,
    n;
  std::vector< std::vector<int> >
    grad_dip_Fnum;
  std::vector< std::vector<double> >
    grad_dip_scl;

  void add_prm(const std::string Fname, const int n,
	       const double bn_min, const double bn_max,
	       const double bn_scl);
  void add_prm(const std::vector<int> &grad_dip_Fnum,
	       const std::vector<double> &grad_dip_scl,
	       const int n, const double bn_min, const double bn_max,
	       const double bn_scl);
  void ini_prm(double *bn);
  void set_prm(double *bn);
  void prt_prm(double *bn) const;
};


struct constr_type {
  // Constraint Type:
  //   alpha_x, alpha_y, beta_x, beta_y, eta_x, eta'_x.
private:

public:
  bool
    ring;
  int
    n_loc, n_constr, n_b3, n_b1,
    n_iter;
  double
    ic[4][2],
    chi2, chi2_prt,
    eps_x_scl,
    eps0_x,              // Hor. emittance [nm.rad].
    drv_terms_simple_scl,
    drv_terms_simple[2],
    nu[2],
    ksi1[2],
    ksi1_ctrl_scl[3],
    ksi1_svd_scl,
    ksi1_svd[2],
    phi_scl,
    phi_tot, phi0,       // Cell bend angle.
    high_ord_achr_scl,
    high_ord_achr_nu[2], // Higher-Order-Achromat Phase Advance.
    mI_scl[2],
    mI0[2],              // -I Transformer.
    alpha_c_scl,         // alpha_c.
    L_scl,
    L0;                  // Cell length.
  std::vector<double> 
    ksi1_ctrl;
  std::vector< std::vector<double> >
    value,
    value_scl,
    mI,
    high_ord_achr_dnu;
  double **Jacobian;
  std::vector<int>
    Fnum,
    Fnum_b3,
    Fnum_b1,
    high_ord_achr_Fnum,
    loc,
    type;
  std::vector< std::vector<int> >
    grad_dip_Fnum_b1,
    mI_loc;
  

  constr_type(void) {
    n_iter = 0; chi2 = 1e30; chi2_prt = 1e30;
    eps_x_scl = 0e0; phi_scl = 0e0;
    ksi1_ctrl_scl[0] = ksi1_ctrl_scl[1] = ksi1_ctrl_scl[2] = 0e0;
    ksi1_svd_scl = 0e0;
    drv_terms_simple_scl = 0e0;
    high_ord_achr_scl = 0e0;
    mI_scl[X_] = mI_scl[Y_] = 0e0;
    L_scl = 0e0;
  }

  void add_constr(const int loc,
		  const double scl1, const double scl2, const double scl3,
		  const double scl4, const double scl5, const double scl6,
		  const double v1, const double v2, const double v3,
		  const double v4, const double v5, const double v6);
  void ini_constr(const bool ring);
  void get_Jacobian(param_type &lat_prms);
  double get_chi2(void) const;
  void get_dchi2(double *df) const;
  void prt_Jacobian(const int n) const;
  void prt_constr(const double chi2);
};


param_type  lat_prms;
constr_type lat_constr;


void prt_name(FILE *outf, const char *name, const string &str, const int len)
{
  int j, k;

  j = 0;
  do {
    fprintf(outf, "%c", name[j]);
    j++;
  } while (name[j] != ' ');
  fprintf(outf, str.c_str());
  for (k = j; k < len; k++)
    fprintf(outf, " ");
}


void param_type::add_prm(const std::string Fname, const int n,
			 const double bn_min, const double bn_max,
			 const double bn_scl)
{
  const bool prt = false;

  this->Fnum.push_back(ElemIndex(Fname.c_str()));
  this->n.push_back(n);
  this->bn_min.push_back(bn_min);
  this->bn_max.push_back(bn_max);
  this->bn_scl.push_back(bn_scl);
  this->l0.push_back(0e0);
  this->n_prm = this->Fnum.size();

  if (prt) printf("add_prm: %2d\n", this->n_prm);
}


void param_type::add_prm(const std::vector<int> &grad_dip_Fnum,
			 const std::vector<double> &grad_dip_scl,
			 const int n, const double bn_min, const double bn_max,
			 const double bn_scl)
{
  const bool prt = false;

  // Flag Gradient Dipole by Negative Family Number.
  this->Fnum.push_back(-grad_dip_Fnum[0]);
  this->grad_dip_Fnum.push_back(grad_dip_Fnum);
  this->grad_dip_scl.push_back(grad_dip_scl);
  this->n.push_back(n);
  this->bn_min.push_back(bn_min);
  this->bn_max.push_back(bn_max);
  this->bn_scl.push_back(bn_scl);
  this->l0.push_back(0e0);
  this->n_prm = this->Fnum.size();

  if (prt) printf("add_prm: %2d\n", this->n_prm);
}


double get_grad_dip_phi(std::vector<int> &Fnum)
{
  int    k, loc;
  double phi;

  phi = 0e0;
  for (k = 0; k < (int)Fnum.size(); k++) {
    loc = Elem_GetPos(Fnum[k], 1);
    phi += Cell[loc].Elem.PL*Cell[loc].Elem.M->Pirho;
  }
  return rad2deg(phi);
}


void param_type::ini_prm(double *bn)
{
  int    i, loc;
  double an;

  bool prt = false;

  for (i = 1; i <= n_prm; i++) {
    if (prt) printf("ini_prm: %2d %3d %2d\n", i, Fnum[i-1], n[i-1]);
    if (n[i-1] > 0)
      // Multipole.
      get_bn_design_elem(abs(Fnum[i-1]), 1, n[i-1], bn[i], an);
    else if (n[i-1] == -1) {
      // Position.
      loc = Elem_GetPos(Fnum[i-1], 1);
      if (Cell[loc+1].Elem.Pkind != drift) {
	printf("\nini_prm: upstream element of %s not a drift %s\n",
	       Cell[loc].Elem.PName, Cell[loc+1].Elem.PName);
	exit(1);
      }
      if (Cell[loc-1].Elem.Pkind != drift) {
	printf("\nini_prm: downstream element of %s not a drift %s\n",
	       Cell[loc].Elem.PName, Cell[loc-1].Elem.PName);
	exit(1);
      }
      bn[i] = 0e0;
      bn_min[i-1] = -(Cell[loc-1].Elem.PL-bn_min[i-1]);
      bn_max[i-1] = Cell[loc+1].Elem.PL - bn_max[i-1];
    } else if (n[i-1] == -2)
      // Length.
      bn[i] = get_L(Fnum[i-1], 1);
    else if (n[i-1] == -3) {
      // Bend angle; L is fixed.
      if (Fnum[i-1] > 0) {
	loc = Elem_GetPos(Fnum[i-1], 1);
	bn[i] = rad2deg(Cell[loc].Elem.PL*Cell[loc].Elem.M->Pirho);
      } else
	bn[i] = get_grad_dip_phi(grad_dip_Fnum[i-1]);
    }
    // Bounded.
    if ((bn_min[i-1] <= bn[i]) && (bn[i] <= bn_max[i-1]))
      bn[i] = bn_internal(bn[i], bn_min[i-1], bn_max[i-1]);
    else {
      loc = Elem_GetPos(Fnum[i-1], 1);
      printf("\nini_prm:\n outside range ");
      prt_name(stdout, Cell[loc].Elem.PName, ":", 8);
      printf(" %10.3e [%10.3e, %10.3e]\n", bn[i], bn_min[i-1], bn_max[i-1]);
      exit(1);
    }
  }
}


void set_phi(const int Fnum, const double phi,
	     const double phi0, const double phi1)
{
  int    loc, k;
  double i_rho;

  loc = Elem_GetPos(Fnum, 1);
  i_rho = deg2rad(phi)/Cell[loc].Elem.PL;
  for (k = 1; k <= GetnKid(Fnum); k++) {
    loc = Elem_GetPos(Fnum, k);
    Cell[loc].Elem.M->Pirho = i_rho;
    Cell[loc].Elem.M->PTx1 = phi0;
    Cell[loc].Elem.M->PTx2 = phi1;
  }
}


void set_grad_dip_phi(const std::vector<int> &Fnum,
		      const std::vector<double> & scl, const double phi_tot)
{
  int    k;
  double phi, phi0, phi1;

  phi1 = 0e0;
  for (k = 0; k < (int)Fnum.size(); k++) {
    phi0 = (k == 0)? 0e0 : -phi;
    phi = phi_tot*scl[k];
    phi1 += phi;
    set_phi(Fnum[k], phi, phi0, phi1);
  }
}


void set_grad_dip_b2(const std::vector<int> &Fnum, const double b2)
{
  int k;

  for (k = 0; k < (int)Fnum.size(); k++)
    set_bn_design_fam(Fnum[k], Quad, b2, 0e0);
}


double get_phi_tot(const int Fnum)
{
  int loc;

  loc = Elem_GetPos(Fnum, 1);
  return GetnKid(Fnum)*rad2deg(Cell[loc].Elem.PL*Cell[loc].Elem.M->Pirho);
}


double get_phi(constr_type &constr)
{
  int    j, k;
  double phi;

  const bool prt = false;

  if (prt) printf("\nget_phi:\n");
  phi = 0e0;
  for (j = 0; j < constr.n_b1; j++) {
    if (constr.Fnum_b1[j] > 0) {
      phi += get_phi_tot(constr.Fnum_b1[j]);
      if (phi_spec_case && (constr.Fnum_b1[j] == ElemIndex("qf4"))) {
	phi += get_phi_tot(ElemIndex("qf4l"));
	phi += get_phi_tot(ElemIndex("qf4_c1"));
      }
      if (prt) {
	printf("  ");
	prt_name(stdout, Cell[Elem_GetPos(constr.Fnum_b1[j], 1)].Elem.PName,
		 "", 8);
	printf("%6.3f\n", phi);
      }
    } else
      for (k = 0; k < (int)constr.grad_dip_Fnum_b1[j].size(); k++) {
	phi += get_phi_tot(constr.grad_dip_Fnum_b1[j][k]);
	if (prt) {
	  printf("  ");
	  prt_name(stdout,
		   Cell[Elem_GetPos(constr.grad_dip_Fnum_b1[j][k], 1)]
		   .Elem.PName, "", 8);
	  printf("%6.3f\n", phi);
	}
      }
  }

  return phi;
}


void phi_corr(constr_type &constr)
{
  // Correct total bend angle.
  int    loc;
  double phi1, phi;

  const bool prt = false;

  loc = Elem_GetPos(constr.Fnum_b1[constr.n_b1-1], 1);
  phi = get_phi(constr);
  phi1 =
    rad2deg(Cell[loc].Elem.PL*Cell[loc].Elem.M->Pirho)
    - (phi-constr.phi0)
    /GetnKid(constr.Fnum_b1[constr.n_b1-1]);
  set_phi(constr.Fnum_b1[constr.n_b1-1], phi1);
 
  constr.phi_tot = get_phi(constr);

  if (prt) {
    printf("\nphi_corr:\n  ");
    prt_name(stdout, Cell[loc].Elem.PName, "", 8);
    printf("\n  %5.3f -> %5.3f (%5.3f)\n", phi, constr.phi_tot, constr.phi0);
  }
}


void set_phi_spec(const char *name,
		  const double phi, const double phi0, double &phi1)
{
  int Fnum, loc;

  Fnum = ElemIndex(name);
  loc = Elem_GetPos(Fnum, 1);
  printf("  %10s %2d %3d\n", Cell[loc].Elem.PName, Fnum, loc);
  phi1 += phi;
  set_phi(Fnum, phi, phi0, phi1);
}


void param_type::set_prm(double *bn)
{
  int    i;
  double bn_ext;

  const bool prt = false;

  const int  n_prt = 6;

  if (prt) printf("\nset_prm:\n  ");
  for (i = 1; i <= n_prm; i++) {
    // Bounded.
    bn_ext = bn_bounded(bn[i], bn_min[i-1], bn_max[i-1]);
    if (n[i-1] > 0)
      // Multipole strength.
      if (Fnum[i-1] > 0)
	set_bn_design_fam(Fnum[i-1], n[i-1], bn_ext, 0e0);
      else
	set_grad_dip_b2(grad_dip_Fnum[i-1], bn_ext);
    else if (n[i-1] == -1) {
      // Position.
      set_ds(Fnum[i-1], bn_ext-l0[i-1]);
      l0[i-1] = bn_ext;
    } else if (n[i-1] == -2) {
      // Length.
      set_L(Fnum[i-1], bn_ext); get_S();
    } else if (n[i-1] == -3) {
      // Bend angle; L is fixed.
      if (Fnum[i-1] > 0) {
	set_phi(Fnum[i-1], bn_ext);
	if (phi_spec_case && (Fnum[i-1] == ElemIndex("qf4"))) {
	  set_phi(ElemIndex("qf4l"), bn_ext);
	  set_phi(ElemIndex("qf4_c1"), bn_ext);
	}
      } else
	set_grad_dip_phi(grad_dip_Fnum[i-1], grad_dip_scl[i-1], bn_ext);
    }
    if (prt) {
      printf(" %12.5e", bn_ext);
      if (i % n_prt == 0) printf("\n  ");
    }
  }
  if (prt && (n_prm % n_prt != 0)) printf("\n");
}


void param_type::prt_prm(double *bn) const
{
  int    i, loc;
  double bn_ext;

  const string labels[] = {"phi", "L  ", "s  ", "   ", "   ", "b_2"};
  const int    labels_offset = 3;

  for (i = 1; i <= n_prm; i++) {
    // Bounded.
    bn_ext = bn_bounded(bn[i], bn_min[i-1], bn_max[i-1]);
    if (Fnum[i-1] > 0) {
      loc = Elem_GetPos(Fnum[i-1], 1);
      printf("    ");
      prt_name(stdout, Cell[loc].Elem.PName, ":", 16);
    } else {
      loc = Elem_GetPos(grad_dip_Fnum[i-1][0], 1);
      printf("    ");
      prt_name(stdout, Cell[loc].Elem.PName, "", -1);
      loc = Elem_GetPos(grad_dip_Fnum[i-1][grad_dip_Fnum[i-1].size()-1], 1);
      printf(" - ");
      prt_name(stdout, Cell[loc].Elem.PName, ":", 7);
    }
    printf("   %3s %9.5f [%9.5f, %9.5f]",
	   labels[n[i-1]+labels_offset].c_str(),
	   bn_ext, bn_min[i-1], bn_max[i-1]);
    if (n[i-1] == -1) {
      printf(" ");
      prt_name(stdout, Cell[loc-1].Elem.PName, ":", 8);
      printf(" %7.5f, ", Cell[loc-1].Elem.PL);
      prt_name(stdout, Cell[loc+1].Elem.PName, ":", 8);
      printf(" %7.5f", Cell[loc+1].Elem.PL);
    }
    printf("\n");
  }
}


double bn_internal(const double bn_bounded,
		   const double bn_min, const double bn_max)
{
  return asin((2e0*(bn_bounded-bn_min))/(bn_max-bn_min)-1e0);
}


double bn_bounded(const double bn_internal,
		  const double bn_min, const double bn_max)
{
  return bn_min + (sin(bn_internal)+1e0)*(bn_max-bn_min)/2e0;
}


void set_ds(const int Fnum, const int Knum, const double ds)
{
  int loc;

  const bool prt = false;

  loc = Elem_GetPos(Fnum, Knum);

  if ((Cell[loc+1].Elem.Pkind == drift) && (Cell[loc-1].Elem.Pkind == drift)) {
    if (Knum % 2 == 1) {
      set_dL(Cell[loc+1].Fnum, Knum, -ds);
      set_dL(Cell[loc-1].Fnum, Knum, ds);
    } else {
      set_dL(Cell[loc+1].Fnum, Knum, ds);
      set_dL(Cell[loc-1].Fnum, Knum, -ds);
    }
  } else {
    printf("\nset_ds: non drift element up or downstream of magnet %d %d\n",
	   Cell[loc+1].Elem.Pkind, Cell[loc-1].Elem.Pkind);
    exit(1);
  }

  if (prt)
    printf("\n  %5s %8.5f\n  %5s %8.5f\n  %5s %8.5f\n",
	   Cell[loc].Elem.PName, ds,
	   Cell[loc+1].Elem.PName, Cell[loc+1].Elem.PL,
	   Cell[loc-1].Elem.PName, Cell[loc-1].Elem.PL);
}


void set_ds(const int Fnum, const double ds)
{
  int k;

  for (k = 1; k <= GetnKid(Fnum); k++)
    set_ds(Fnum, k, ds);
}


void set_phi(const int Fnum, const double phi)
{
  int    loc, k;
  double i_rho;

  loc = Elem_GetPos(Fnum, 1);
  i_rho = deg2rad(phi)/Cell[loc].Elem.PL;
  for (k = 1; k <= GetnKid(Fnum); k++) {
    loc = Elem_GetPos(Fnum, k);
    Cell[loc].Elem.M->Pirho = i_rho;
    Cell[loc].Elem.M->PTx1 = phi/2e0;
    Cell[loc].Elem.M->PTx2 = phi/2e0;
  }
}


void set_dphi(const int Fnum, const double dphi)
{
  int    loc, k;
  double phi, i_rho;

  loc = Elem_GetPos(Fnum, 1);
  phi = rad2deg(Cell[loc].Elem.PL*Cell[loc].Elem.M->Pirho);
  i_rho = deg2rad(phi+dphi)/Cell[loc].Elem.PL;
  for (k = 1; k <= GetnKid(Fnum); k++) {
    loc = Elem_GetPos(Fnum, k);
    Cell[loc].Elem.M->Pirho = i_rho;
  }
}


void constr_type::add_constr(const int loc,
			     const double scl1, const double scl2,
			     const double scl3, const double scl4,
			     const double scl5, const double scl6,
			     const double v1, const double v2, const double v3,
			     const double v4, const double v5, const double v6)
{
  // Parameters are: alpha_x,y, beta_x,y, eta_x, eta'_x.
  const int                 n     = 6;
  const double              scl[] = {scl1, scl2, scl3, scl4, scl5, scl6};
  const double              val[] = {v1,   v2,   v3,   v4,   v5,   v6};
  const std::vector<double> vec1(scl, scl+n);
  const std::vector<double> vec2(val, val+n);

  this->loc.push_back(loc);
  this->value_scl.push_back(vec1);
  this->value.push_back(vec2);
  this->n_loc = this->loc.size();
}


double get_eps_x1(const bool track)
{
  // eps_x [nm.rad].
  long int     lastpos;
  double       I[6], eps_x;
  ss_vect<tps> A;

  const bool prt = false;

  if (track) {
    globval.emittance = true;
    // A = get_A(ic[0], ic[1], ic[2], ic[3]);
    putlinmat(6, globval.Ascr, A);
    Cell_Pass(0, globval.Cell_nLoc, A, lastpos);
    // Cell_Twiss(0, globval.Cell_nLoc, A, false, false, 0e0);
    globval.emittance = false;
  }

  get_I(I, false);

  eps_x = 1470e0*sqr(globval.Energy)*I[5]/(I[2]-I[4]);

  if (prt) {
    printf("\neps_x = %5.3f nm.rad\n", eps_x);
    printf("J_x   = %5.3f, J_z = %5.3f\n", 1.0-I[4]/I[2], 2.0+I[4]/I[2]);
  }

  return eps_x;
}


double get_lin_opt(constr_type &constr)
{
  double       eps_x;
  ss_vect<tps> A;

  if (lat_constr.ring) {
    Ring_GetTwiss(true, 0e0);
    eps_x = get_eps_x1(true);
  } else {
    globval.emittance = true;
    A = get_A(constr.ic[0], constr.ic[1], constr.ic[2], constr.ic[3]);
    Cell_Twiss(0, globval.Cell_nLoc, A, false, false, 0e0);
    eps_x = get_eps_x1(false);
    globval.emittance = false;
  }

  return eps_x;
}


void constr_type::ini_constr(const bool ring)
{
  this->ring = ring;
  n_b3 = Fnum_b3.size();
  n_b1 = Fnum_b1.size();
}


double get_constr(const int loc, const int k)
{
  double constr = 0e0;

  const bool prt = false;

  switch (k) {
  case 0:
    constr = Cell[loc].Alpha[X_];
    break;
  case 1:
    constr = Cell[loc].Alpha[Y_];
    break;
  case 2:
    constr = Cell[loc].Beta[X_];
    break;
  case 3:
    constr = Cell[loc].Beta[Y_];
    break;
  case 4:
    constr = Cell[loc].Eta[X_];
    break;
  case 5:
    constr = Cell[loc].Etap[X_];
    break;
  default:
    printf("\nget_constr: unknown constr %d\n", k);
    exit(1);
    break;
  }

  if (prt) printf("get_constr: %1d %10.3e\n", k, constr);

  return constr;
}


double get_eps(const int n)
{
  double eps = 0e0;

  const bool   prt       = false;
  const double prm_eps[] = { 1e-3, 1e-4, 1e-4, 1e-4 };

  switch (n) {
  case -3:
    // Bend angle.
    eps = prm_eps[0];
    break;
  case -2:
    // Length.
    eps = prm_eps[1];
    break;
  case -1:
    // Position.
    eps = prm_eps[2];
    break;
  case 2:
    // Quadrupole strength.
    eps = prm_eps[3];
    break;
  default:
    printf("\nget_eps: unknown parameter type %d\n", n);
    exit(1);
    break;
  }

  if (prt) printf("get_eps: %2d %10.3e\n", n, eps);

  return eps;
}


void constr_dparam(const int Fnum, const int n, const double eps)
{

  switch (n) {
  case -3:
    // Bend angle.
    set_dphi(Fnum, eps);
    phi_corr(lat_constr);
    break;
  case -2:
    // Length.
    set_dL(Fnum, eps);
    break;
  case -1:
    // Position.
    set_ds(Fnum, eps);
    break;
  case 2:
    // Quadrupole strength.
    set_dbn_design_fam(Fnum, n, eps, 0e0);
    break;
  default:
    printf("\nconstr_dparam: unknown parameter type %d\n", n);
    break;
  }
}


void constr_type::get_Jacobian(param_type &lat_prms)
{
  int                 i, j, k, ind = 0;
  double              eps;
  std::vector<double> row;

  const int n_type = 5;

  n_constr = 0;
  for (j = 0; j < n_loc; j++)
    for (k = 0; k < n_type; k++)
      if (value_scl[j][k] != 0e0) n_constr++;
    
  Jacobian = dmatrix(1, n_constr, 1, lat_prms.n_prm);

  printf("\nget_Jacobian: %d %d\n", n_constr, lat_prms.n_prm);

  for (i = 0; i < lat_prms.n_prm; i++) {
    eps = get_eps(lat_prms.n[i]);

    constr_dparam(lat_prms.Fnum[i], lat_prms.n[i], eps);
    get_lin_opt(lat_constr);
    for (j = 0; j < n_loc; j++) {
      ind = 1;
      for (k = 0; k < n_type; k++)
	if (value_scl[j][k] != 0e0) {
	  Jacobian[ind][i+1] = get_constr(loc[j], k);
	  ind++;
	}
    }
 
    constr_dparam(lat_prms.Fnum[i], lat_prms.n[i], -2e0*eps);
    get_lin_opt(lat_constr);
    for (j = 0; j < n_loc; j++) {
      ind = 1;
      for (k = 0; k < n_type; k++)
    	if (value_scl[j][k] != 0e0) {
    	  Jacobian[ind][i+1] -= get_constr(loc[j], k);
    	  Jacobian[ind][i+1] /= 2e0*eps;
    	  ind++;
    	}
    }

    constr_dparam(lat_prms.Fnum[i], lat_prms.n[i], eps);
  }

}


void constr_type::prt_Jacobian(const int n) const
{
  int i, j;

  printf("\n");
  for (i = 1; i <= n_constr; i++) {
    for (j = 1; j <= n; j++)
      printf(" %12.5e", Jacobian[i][j]);
    printf("\n");
  }
}


double constr_type::get_chi2(void) const
{
  int    j, k;
  double chi2, mean, geom_mean;

  const bool prt = false;

  chi2 = 0e0; 
  chi2 += eps_x_scl*sqr(eps_x-eps0_x);

  for (k = 0; k < n_loc; k++) {
    chi2 += value_scl[k][0]*sqr(Cell[loc[k]].Alpha[X_]-value[k][0]);
    chi2 += value_scl[k][1]*sqr(Cell[loc[k]].Alpha[Y_]-value[k][1]);
    chi2 += value_scl[k][2]*sqr(Cell[loc[k]].Beta[X_]-value[k][2]);
    chi2 += value_scl[k][3]*sqr(Cell[loc[k]].Beta[Y_]-value[k][3]);
    chi2 += value_scl[k][4]*sqr(Cell[loc[k]].Eta[X_]-value[k][4]);
    chi2 += value_scl[k][5]*sqr(Cell[loc[k]].Etap[X_]-value[k][5]);
  }

  chi2 += L_scl*sqr(Cell[globval.Cell_nLoc].S-L0);
  for (j = 0; j < (int)mI.size(); j++)
    for (k = 0; k < 2; k++)
      chi2 += mI_scl[k]*sqr(mI[j][k]-mI0[k]);

  for (k = 0; k < 2; k++)
    chi2 += drv_terms_simple_scl*drv_terms_simple[k];

  for (j = 0; j < (int)ksi1_ctrl.size(); j++)
    if (ksi1_ctrl_scl[j] != 0e0) {
      if (j < 2)
	chi2 += ksi1_ctrl_scl[j]/sqr(ksi1_ctrl[j]);
      else
	chi2 += ksi1_ctrl_scl[j]*sqr(ksi1_ctrl[j]);
    }
  
  if (ksi1_svd_scl != 0) {
    mean = (ksi1_svd[X_]+ksi1_svd[Y_])/2e0;
    geom_mean = sqrt(ksi1_svd[X_]*ksi1_svd[Y_]);
    chi2 += ksi1_svd_scl/sqr(geom_mean/mean);
  }

  for (j = 0; j < (int)high_ord_achr_dnu.size(); j++)
    for (k = 0; k < 2; k++)
      chi2 +=
	high_ord_achr_scl*sqr(high_ord_achr_dnu[j][k]-high_ord_achr_nu[k]);

  if (alpha_c_scl != 0)
    chi2 += alpha_c_scl*sqr(1e0/globval.Alphac);

  if (prt) printf("\nget_chi2: %11.5e (%11.5e)\n", chi2, this->chi2);

  return chi2;
}


void constr_type::get_dchi2(double *df) const
{
  int    k, loc;
  double eps;

  const bool prt = false;

  for (k = 0; k < lat_prms.n_prm; k++) {
    eps = get_eps(lat_prms.n[k]);

    if (prt) {
      printf("\nget_dchi2: ");
      loc = Elem_GetPos(lat_prms.Fnum[k], 1);
      prt_name(stdout, Cell[loc].Elem.PName, ":", 6);
      printf(" %2d %10.3e", lat_prms.n[k], eps);
    }

    constr_dparam(lat_prms.Fnum[k], lat_prms.n[k], eps);
    eps_x = get_lin_opt(lat_constr);
    df[k+1] = lat_constr.get_chi2();

    constr_dparam(lat_prms.Fnum[k], lat_prms.n[k], -2e0*eps);
    eps_x = get_lin_opt(lat_constr);
    df[k+1] -= lat_constr.get_chi2();
    df[k+1] /= 2e0*eps;

    constr_dparam(lat_prms.Fnum[k], lat_prms.n[k], eps);
  }

  // Avoid: "warning: deprecated conversion from string constant to ‘char*’".
  // dvdump(stdout,
  // 	 (char *)"\nget_dchi2:", df, lat_prms.n_prm, (char *)" %12.5e");
}


int not_zero(const double a)
{
  return (a != 0e0)? 1 : 0; 
}


void prt_val(const constr_type &constr)
{
  int    k, loc;

  printf("\n    Loc.        alpha_x   alpha_y  beta_x   beta_y"
	 "    eta_x    eta'_x\n");
  for (k = 0; k < constr.n_loc; k++) {
    loc = constr.loc[k];
    printf("    ");
    prt_name(stdout, Cell[loc].Elem.PName, ":", 8);
    printf("  %8.5f  %8.5f %8.5f %8.5f %8.5f %8.5f\n",
	   not_zero(constr.value_scl[k][0])*Cell[loc].Alpha[X_],
	   not_zero(constr.value_scl[k][1])*Cell[loc].Alpha[Y_],
	   not_zero(constr.value_scl[k][2])*Cell[loc].Beta[X_],
	   not_zero(constr.value_scl[k][3])*Cell[loc].Beta[Y_],
	   not_zero(constr.value_scl[k][4])*Cell[loc].Eta[X_],
	   not_zero(constr.value_scl[k][5])*Cell[loc].Etap[X_]);
  }
}


void prt_h(const constr_type &constr)
{
  int    k;
  double b3L, a3L, b3, a3;

  printf("    drv. terms   = [%9.3e, %9.3e]\n",
	 sqrt(constr.drv_terms_simple[X_]), sqrt(constr.drv_terms_simple[Y_]));

  printf("    b_3L         = [");
  for (k = 0; k < constr.n_b3; k++) {
    get_bnL_design_elem(constr.Fnum_b3[k], 1, Sext, b3L, a3L);
    printf("%10.3e", b3L);
    if (k != constr.n_b3-1) printf(", ");
  }
  printf("]\n");
  printf("    b_3          = [");
  for (k = 0; k < constr.n_b3; k++) {
    get_bn_design_elem(constr.Fnum_b3[k], 1, Sext, b3, a3);
    printf("%10.3e", b3);
    if (k != constr.n_b3-1) printf(", ");
  }
  printf("]\n");
}


void prt_high_ord_achr(const constr_type &constr)
{
  int j;

  printf("\n  Higher-Order-Achromat:\n    [%7.5f, %7.5f]\n\n",
	 constr.high_ord_achr_nu[X_], constr.high_ord_achr_nu[Y_]);
  for (j = 0; j < (int)constr.high_ord_achr_dnu.size(); j++)
    printf("    [%7.5f, %7.5f]\n",
	   constr.high_ord_achr_dnu[j][X_], constr.high_ord_achr_dnu[j][Y_]);
}


void constr_type::prt_constr(const double chi2)
{
  int    loc, k;
  double phi, mean, geom_mean;

  mean = (ksi1_svd[X_]+ksi1_svd[Y_])/2e0;
  geom_mean = sqrt(ksi1_svd[X_]*ksi1_svd[Y_]);
  printf("\n%3d chi2: %11.5e -> %11.5e\n", n_iter, this->chi2_prt, chi2);
  this->chi2_prt = chi2;
  printf("\n  Linear Optics:\n");
  printf("    eps_x        = %5.3f (%5.3f)\n"
	 "    nu           = [%5.3f, %5.3f]\n"
	 "    ksi1         = [%5.3f, %5.3f]\n",
	 eps_x, eps0_x, nu[X_], nu[Y_], ksi1[X_], ksi1[Y_]);
  if (lat_constr.ksi1_ctrl.size() != 0)
    printf("    ksi1_ctrl    = [%5.3f, %5.3f, %5.3f]\n",
	   lat_constr.ksi1_ctrl[0], lat_constr.ksi1_ctrl[1],
	   lat_constr.ksi1_ctrl[2]);
  printf("    svd          = [%9.3e, %9.3e] %6.4f\n",
	 lat_constr.ksi1_svd[X_], lat_constr.ksi1_svd[Y_], geom_mean/mean);
  if (phi_scl != 0e0) {
    loc = Elem_GetPos(Fnum_b1[n_b1-1], 1);
    phi = rad2deg(Cell[loc].Elem.PL*Cell[loc].Elem.M->Pirho);
    printf("    phi          = %7.5f (%7.5f)\n    ", phi_tot, phi0);
    prt_name(stdout, Cell[loc].Elem.PName, "_phi:", 7);
    printf(" = %7.5f\n", phi);
  }
  if (alpha_c_scl != 0e0)
    printf("    alpha_c      = %9.3e\n", globval.Alphac);
  if (L_scl != 0e0)
    printf("    L            = %7.5f (%7.5f)\n", Cell[globval.Cell_nLoc].S, L0);

  prt_h(*this);

  if (high_ord_achr_scl != 0e0) prt_high_ord_achr(*this);

  if ((mI_scl[X_] != 0e0) || (mI_scl[Y_] != 0e0)) {
    printf("\n  -I Transf.:\n");
    for (k = 0; k < (int)mI.size(); k++)
      printf("    [%7.5f, %7.5f] [%3d, %3d]\n",
	     mI[k][X_], mI[k][Y_], mI_loc[k][0], mI_loc[k][1]);
  }

  prt_val(*this);
}


void get_S(void)
{
  int k;

  Cell[0].S = 0e0;
  for (k = 1; k <= globval.Cell_nLoc; k++)
    Cell[k].S = Cell[k-1].S + Cell[k].Elem.PL;
}


void get_nu(ss_vect<tps> &A, const double delta, double nu[])
{
  long int     lastpos;
  ss_vect<tps> A_delta;

  const bool prt = !false;

  A_delta = get_A_CS(2, A, nu); A_delta[delta_] += delta;
  Cell_Pass(0, globval.Cell_nLoc, A_delta, lastpos);
  get_A_CS(2, A_delta, nu);
  if (prt) printf("\n[%18.16f, %18.16f]\n", nu[X_], nu[Y_]);
}


void get_ksi1(ss_vect<tps> &A, double ksi1[])
{
  int    k;
  double nu0[2], nu1[2];

  const double delta = 1e-8;

  get_nu(A, -delta, nu0); get_nu(A, delta, nu1);
  for (k = 0; k < 2; k++)
    ksi1[k] = (nu1[k]-nu0[k])/(2e0*delta);
}


void fit_ksi1(const std::vector<int> &Fnum_b3,
	      const double ksi_x, const double ksi_y, const double db3L,
	      double svd[])
{
  int                 n_b3, j, k, n_svd;
  double              **A, **U, **V, *w, *b, *x, b3L, a3L;

  const bool   prt = false;
  const int    m   = 2;
  const double
    ksi0[]  = {ksi_x, ksi_y},
    svd_cut = 1e-10;

  n_b3 = Fnum_b3.size();

  A = dmatrix(1, m, 1, n_b3); U = dmatrix(1, m, 1, n_b3);
  V = dmatrix(1, n_b3, 1, n_b3);
  w = dvector(1, n_b3); b = dvector(1, m); x = dvector(1, n_b3);

  // Zero sextupoles to track linear chromaticity.
  no_sxt();

  for (k = 1; k <= n_b3; k++) {
    set_dbnL_design_fam(Fnum_b3[k-1], Sext, db3L, 0e0);
    Ring_Getchrom(0e0);
    if (prt)
      printf("\nfit_ksi1: ksi1+ = [%9.5f, %9.5f]\n",
	     globval.Chrom[X_], globval.Chrom[Y_]);

    for (j = 1; j <= m; j++)
      A[j][k] = globval.Chrom[j-1];
    set_dbnL_design_fam(Fnum_b3[k-1], Sext, -2e0*db3L, 0e0);
    Ring_Getchrom(0e0);
    if (prt)
      printf("fit_ksi1: ksi1- = [%9.5f, %9.5f]\n",
	 globval.Chrom[X_], globval.Chrom[Y_]);
    for (j = 1; j <= 2; j++) {
      A[j][k] -= globval.Chrom[j-1]; A[j][k] /= 2e0*db3L;
    }

    set_dbnL_design_fam(Fnum_b3[k-1], Sext, db3L, 0e0);
  }

  Ring_Getchrom(0e0);
  if (prt)
    printf("\nfit_ksi1: ksi1  = [%9.5f, %9.5f]\n",
	   globval.Chrom[X_], globval.Chrom[Y_]);
  for (j = 1; j <= 2; j++) {
    b[j] = -(globval.Chrom[j-1]-ksi0[j-1]);
    lat_constr.ksi1[j-1] = globval.Chrom[j-1];
  }

  dmcopy(A, m, n_b3, U); dsvdcmp(U, m, n_b3, w, V);

  if (prt) printf("\nfit_ksi1:\n  singular values:\n  ");
  n_svd = 0;
  for (j = 1; j <= n_b3; j++) {
    if (prt) printf("%11.3e", w[j]);
    if (w[j] < svd_cut) {
      w[j] = 0e0;
      if (prt) printf(" (zeroed)");
    } else {
      if (n_svd > 2) {
	if (prt) printf("fit_ksi1: more than 2 non-zero singular values");
	exit(1);
      }
      svd[n_svd] = w[j];
      n_svd++;
    }
  }
  if (prt) printf("\n");

  dsvbksb(U, w, V, m, n_b3, b, x);

  for (k = 1; k <= n_b3; k++)
    set_dbnL_design_fam(Fnum_b3[k-1], Sext, x[k], 0e0);

  if (prt) {
    printf("  b3:\n  ");
    for (k = 0; k < n_b3; k++) {
      get_bn_design_elem(Fnum_b3[k], 1, Sext, b3L, a3L);
      printf(" %9.5f", b3L);
    }
    printf("\n");
  }

  free_dmatrix(A, 1, m, 1, n_b3); free_dmatrix(U, 1, m, 1, n_b3);
  free_dmatrix(V, 1, n_b3, 1, n_b3);
  free_dvector(w, 1, n_b3); free_dvector(b, 1, m); free_dvector(x, 1, n_b3);
}


void prt_drift(FILE *outf, CellType &Cell)
{
  fprintf(outf, " drift, l = %11.8f;\n", Cell.Elem.PL);
}


void prt_dip(FILE *outf, CellType &Cell)
{
  double phi;

  phi = rad2deg(Cell.Elem.PL*Cell.Elem.M->Pirho);
  fprintf(outf, " bending, l = %11.8f,"
	  "\n    t = %11.8f, t1 = %11.8f, t2 = %11.8f, k = %11.8f,"
	  "\n    n = nbend, method = 4;\n",
	  Cell.Elem.PL, phi, Cell.Elem.M->PTx1,
	  Cell.Elem.M->PTx2, Cell.Elem.M->PBpar[Quad+HOMmax]);
}


void prt_grad_dip(FILE *outf, const std::vector<int> &Fnum)
{
  int k, loc;

  const bool prt = false;

  if (prt) printf("  prt_grad_dip: %d\n", (int)Fnum.size());
  for (k = 0; k < (int)Fnum.size(); k++) {
    loc = Elem_GetPos(Fnum[k], 1);
    prt_name(outf, Cell[loc].Elem.PName, ":", 8);
    prt_dip(outf, Cell[loc]);
  }
}


void prt_elem(FILE *outf, const param_type &lat_prms, const int n)
{
  int loc;

  const bool prt = false;

  loc = Elem_GetPos(abs(lat_prms.Fnum[n-1]), 1);
  if (prt) {
    printf("prt_elem:\n  ");
    prt_name(stdout, Cell[loc].Elem.PName, "\n", 0);
  }
  if (Cell[loc].Elem.Pkind == drift) {
    prt_name(outf, Cell[loc].Elem.PName, ":", 8);
    prt_drift(outf, Cell[loc]);
  } else if (Cell[loc].Elem.Pkind == Mpole) {
    if (prt) printf("  n_design:     %1d\n", Cell[loc].Elem.M->n_design);
    if (Cell[loc].Elem.M->n_design == Dip) {
      if (lat_prms.Fnum[n-1] > 0) {
	prt_name(outf, Cell[loc].Elem.PName, ":", 8);
	prt_dip(outf, Cell[loc]);
      } else
      	prt_grad_dip(outf, lat_prms.grad_dip_Fnum[n-1]);
    } else if (Cell[loc].Elem.M->n_design == Quad) {
      prt_name(outf, Cell[loc].Elem.PName, ":", 8);
      fprintf(outf, " quadrupole, l = %11.8f, k = %11.8f, n = nquad"
	      ", method = 4;\n",
	      Cell[loc].Elem.PL, Cell[loc].Elem.M->PBpar[Quad+HOMmax]);
    }
  } else {
    printf("\nprt_elem: %s %d\n", Cell[loc].Elem.PName, Cell[loc].Elem.Pkind);
    exit(1);
  }
}


bool Fam_printed(const int Fnum, const std::vector<int> &Fam_prt)
{
  int k;

  for (k = 0; k < (int)Fam_prt.size(); k++)
    if (Fnum == Fam_prt[k]) return true;

  return false;
}


void prt_b2(const param_type &lat_prms)
{
  long int         loc;
  int              k;
  std::vector<int> Fam_prt;
  FILE             *outf;

  std::string file_name = "opt_des_1_b2.out";

  outf = file_write(file_name.c_str());

  for (k = 0; k < lat_prms.n_prm; k++) {
    if (!Fam_printed(abs(lat_prms.Fnum[k]), Fam_prt)) {
      loc = Elem_GetPos(abs(lat_prms.Fnum[k]), 1);
      if (lat_prms.n[k] != -1)
	prt_elem(outf, lat_prms, k+1);
      else {
	// Position.
	fprintf(outf, "\n  ");
	prt_name(outf, Cell[loc+1].Elem.PName, ":", 8);
	prt_drift(outf, Cell[loc+1]);
	fprintf(outf, "  ");
	prt_name(outf, Cell[loc-1].Elem.PName, ":", 8);
	prt_drift(outf, Cell[loc-1]);
      }
    }
    Fam_prt.push_back(abs(lat_prms.Fnum[k]));
  }

  fclose(outf);
}


void prt_b3(const std::vector<int> &Fnum)
{
  long int loc;
  int      k;
  FILE     *outf;

  std::string file_name = "opt_des_1_b3.out";

  outf = file_write(file_name.c_str());

  for (k = 0; k < (int)Fnum.size(); k++) {
    loc = Elem_GetPos(Fnum[k], 1);
    prt_name(outf, Cell[loc].Elem.PName, ":", 8);
    fprintf(outf, " sextupole, l = %11.8f, k = %13.8f, n = nsext"
	    ", method = 4;\n",
	    Cell[loc].Elem.PL, Cell[loc].Elem.M->PBpar[Sext+HOMmax]);
  }

  fclose(outf);
}


bool get_nu(double nu[])
{
  long int     lastpos;
  int          k;
  double       cosmu;
  ss_vect<tps> ps;

  bool prt = false;

  ps.identity();
  Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
  for (k = 0; k < 2; k++) {
    cosmu = (ps[2*k][2*k]+ps[2*k+1][2*k+1])/2e0;
    if (fabs(cosmu) < 1e0) {
      nu[k] = acos(cosmu)/(2e0*M_PI);
      if (ps[2*k][2*k+1] < 0e0) nu[k] = 1e0 - nu[k];
    } else {
      nu[X_] = NAN; nu[Y_] = NAN;
      printf("\nget_nu: unstable in plane %d %10.3e\n", k, cosmu);
      return false;
    }
  }
  if (prt) printf("\nget_nu: nu = [%7.5f, %7.5f]\n", nu[X_], nu[Y_]);
  return true;
}


double h_abs_ijklm(const tps &h_re, const tps &h_im, const int i, const int j,
		   const int k, const int l, const int m)
{
  int      i1;
  long int jj[ss_dim];

  for (i1 = 0; i1 < nv_tps; i1++)
    jj[i1] = 0;
  jj[x_] = i; jj[px_] = j; jj[y_] = k; jj[py_] = l; jj[delta_] = m;
  return sqrt(sqr(h_re[jj])+sqr(h_im[jj]));
}


void get_drv_terms(constr_type &constr)
{
  int    j, k, n_kid;
  double b3L, a3L;

  for (k = 0; k < 2; k++)
    constr.drv_terms_simple[k] = 0e0;
  for (k = 0; k < constr.n_b3; k++) {
    get_bnL_design_elem(constr.Fnum_b3[k], 1, Sext, b3L, a3L);
    n_kid = GetnKid(constr.Fnum_b3[k]);
    for (j = 0; j < 2; j++) {
      constr.drv_terms_simple[j] +=
	n_kid*sqr(b3L*Cell[Elem_GetPos(constr.Fnum_b3[k], 1)].Beta[j]);
    }
  }
}


void fit_powell(param_type &lat_prms, const double eps, double (*f)(double *))
{
  const int
    n_prm   = lat_prms.n_prm,
    n_pt    = n_prm + 2,
    n_w     = (n_pt+13)*(n_pt+n_prm)+3*n_prm*(n_prm+3)/2,
    n_prt   = 2,
    max_fun = 5000;

  const double
    // rho_beg = lat_prms.dbn[0],
    rho_beg = 1e-1,
    rho_end = 1e-6;

  int          n_b2, i, j, iter;
  double       *b2, **xi, fret, eps_x, w[n_w], x[n_prm];
  ss_vect<tps> A;

  const double ftol = 1e-8;

  n_b2 = lat_prms.n_prm;

  b2 = dvector(1, n_b2); xi = dmatrix(1, n_b2, 1, n_b2);

  lat_prms.ini_prm(b2);
  f(b2);

  if (false) {
    lat_constr.get_Jacobian(lat_prms);
    lat_constr.prt_Jacobian(n_b2);
    exit(0);
  }

  if (true) {
    // Set initial directions (unit vectors).
    for (i = 1; i <= n_b2; i++)
      for (j = 1; j <= n_b2; j++)
	xi[i][j] = (i == j)? eps : 0e0;

    dpowell(b2, xi, n_b2, ftol, &iter, &fret, f);
  } else {
    // for (i = 1; i <= n_prm; i++)
    //   x[i-1] = bn[i];
    // newuoa_(n_prm, n_pt, x, rho_beg, rho_end, n_prt, max_fun, w);
  }

  printf("\n  iter = %d fret = %12.5e\n", iter, fret);
  printf("b2s:\n");
  lat_prms.prt_prm(b2);
  // lat_prms.set_prm(b2);
  // eps_x = get_lin_opt(false);
  // f_prt(b2);

  free_dvector(b2, 1, n_b2); free_dmatrix(xi, 1, n_b2, 1, n_b2);
}


void f_der(double *b2, double *df)
{
  lat_prms.set_prm(b2);
  if (lat_constr.phi_scl != 0e0) phi_corr(lat_constr);

  lat_constr.get_dchi2(df);
}


void fit_conj_grad(param_type &lat_prms, double (*f)(double *))
{
  int          n_b2, iter;
  double       *b2, fret, eps_x;
  ss_vect<tps> A;

  const double ftol = 1e-8;

  n_b2 = lat_prms.n_prm;

  b2 = dvector(1, n_b2);

  lat_prms.ini_prm(b2);
  f(b2);

  dfrprmn(b2, n_b2, ftol, &iter, &fret, f, f_der);

  printf("\n  iter = %d fret = %12.5e\n", iter, fret);
  printf("b2s:\n");
  lat_prms.prt_prm(b2);
  // lat_prms.set_prm(b2);
  // eps_x = get_lin_opt(false);
  // f_prt(b2);

  free_dvector(b2, 1, n_b2);
}


void prt_f(double *b2, const double chi2, constr_type &lat_constr,
	   param_type &lat_prms, const bool chrom)
{
  lat_constr.prt_constr(chi2);

  printf("\n  Parameters:\n");
  lat_prms.prt_prm(b2);
  prtmfile("flat_file.fit");
  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);
  if (chrom) prt_chrom_lat();

  prt_b2(lat_prms);
  prt_b3(lat_constr.Fnum_b3);
}


void get_high_ord_achr(constr_type &constr)
{
  int j, k;

  for (j = 1; j < (int)constr.high_ord_achr_Fnum.size(); j++)
    for (k = 0; k < 2; k++)
      constr.high_ord_achr_dnu[j-1][k] =
	Cell[constr.high_ord_achr_Fnum[j]].Nu[k]
	- Cell[constr.high_ord_achr_Fnum[j-1]].Nu[k];
}


void get_mI(constr_type &constr)
{
  int                 j, k;
  std::vector<double> dnu;

  constr.mI.clear();
  for (j = 0; j < (int)constr.mI_loc.size(); j++) {
    for (k = 0; k < 2; k++)
      dnu.push_back(Cell[lat_constr.mI_loc[j][1]].Nu[k]
		    -Cell[lat_constr.mI_loc[j][0]].Nu[k]);
    constr.mI.push_back(dnu);
    dnu.clear();
  }
}


void get_ksi1_ctrl(constr_type &constr)
{
  int k, loc[3];

  for (k = 0; k < 3; k++)
    loc[k] = Elem_GetPos(constr.Fnum_b3[k], 1);
  constr.ksi1_ctrl.clear();
  constr.ksi1_ctrl.push_back((Cell[loc[0]].Beta[Y_]-Cell[loc[0]].Beta[X_])
			     *Cell[loc[0]].Eta[X_]);
  constr.ksi1_ctrl.push_back((Cell[loc[1]].Beta[X_]-Cell[loc[1]].Beta[Y_])
			     *Cell[loc[1]].Eta[X_]);
  constr.ksi1_ctrl.push_back((Cell[loc[0]].Beta[Y_]-Cell[loc[2]].Beta[Y_]));
}


double f_match(double *b2)
{
  double      chi2;
  static bool first = true;

  const int n_prt = 5;

  lat_prms.set_prm(b2);
  if (lat_constr.phi_scl != 0e0) phi_corr(lat_constr);

  eps_x = get_lin_opt(lat_constr);

  if (lat_constr.high_ord_achr_scl != 0e0) get_high_ord_achr(lat_constr);

  chi2 = lat_constr.get_chi2();

  if (first || (chi2 < lat_constr.chi2)) {
    first = false;
    if (lat_constr.n_iter % n_prt == 0)
      prt_f(b2, chi2, lat_constr, lat_prms, false);
    lat_constr.n_iter++;
    lat_constr.chi2 = chi2;
  }

  return chi2;
}


double f_achrom(double *b2)
{
  bool        stable;
  double      chi2;
  static bool first = true;

  const int n_prt = 5;

  lat_prms.set_prm(b2);
  if (lat_constr.phi_scl != 0e0) phi_corr(lat_constr);

  if (lat_constr.ring) stable = get_nu(lat_constr.nu);

  if ((lat_constr.ring && stable) || !lat_constr.ring) {
    eps_x = get_lin_opt(lat_constr);

    fit_ksi1(lat_constr.Fnum_b3, 0e0, 0e0, 1e1, lat_constr.ksi1_svd);

    if (lat_constr.drv_terms_simple_scl != 0e0) get_drv_terms(lat_constr);

    if (lat_constr.high_ord_achr_scl != 0e0) get_high_ord_achr(lat_constr);

    if ((lat_constr.mI_scl[X_] != 0e0) || (lat_constr.mI_scl[Y_] != 0e0))
      get_mI(lat_constr);

    if ((lat_constr.ksi1_ctrl_scl[0] != 0e0)
	|| (lat_constr.ksi1_ctrl_scl[1] != 0e0)
	|| (lat_constr.ksi1_ctrl_scl[2] != 0e0))
      get_ksi1_ctrl(lat_constr);

    chi2 = lat_constr.get_chi2();

    if (first || (chi2 < lat_constr.chi2)) {
      first = false;
      if (lat_constr.n_iter % n_prt == 0)
	prt_f(b2, chi2, lat_constr, lat_prms, true);

      lat_constr.n_iter++;
      lat_constr.chi2 = chi2;
    }
  } else
    chi2 = 1e30;

  return chi2;
}


void prt_prms(constr_type &constr)
{
  printf("\n  eps_x_scl            = %9.3e\n"
	 "  ksi1_svd_scl         = %9.3e\n"
	 "  drv_terms_simple_scl = %9.3e\n"
	 "  mI_nu_ref            = [%7.5f, %7.5f]\n"
	 "  mI_scl               = [%9.3e, %9.3e]\n"
	 "  high_ord_achr_scl    = %9.3e\n"
	 "  phi_scl              = %9.3e\n",
	 constr.eps_x_scl, constr.ksi1_svd_scl,
	 constr.drv_terms_simple_scl,
	 mI_nu_ref[X_], mI_nu_ref[Y_],
	 constr.mI_scl[X_], constr.mI_scl[Y_],
	 constr.high_ord_achr_scl,
	 constr.phi_scl);
}


void opt_mI_std(param_type &prms, constr_type &constr)
{
  // Parameter Type:
  //   Bend Angle  -3,
  //   Length      -2,
  //   Position    -1,
  //   Quadrupole   2.
  int                 k, n;
  std::vector<int>    grad_dip_Fnum, mI_loc;
  std::vector<double> grad_dip_scl;

  const bool
    dphi          = !false,
    long_grad_dip = !false;

  printf("\n opt_mI_std:\n relaxed: %d\n", relaxed);

  // Standard Cell.
  grad_dip_scl.push_back(0.129665);
  grad_dip_scl.push_back(0.149256);
  grad_dip_scl.push_back(0.174737);
  grad_dip_scl.push_back(0.216027);
  grad_dip_scl.push_back(0.330314);

  grad_dip_Fnum.push_back(ElemIndex("dl1a_1"));
  grad_dip_Fnum.push_back(ElemIndex("dl1a_2"));
  grad_dip_Fnum.push_back(ElemIndex("dl1a_3"));
  grad_dip_Fnum.push_back(ElemIndex("dl1a_4"));
  grad_dip_Fnum.push_back(ElemIndex("dl1a_5"));
  if (dphi)
    prms.add_prm(grad_dip_Fnum, grad_dip_scl, -3, -20.0, 20.0, 1.0);
  if (long_grad_dip)
    prms.add_prm(grad_dip_Fnum, grad_dip_scl,  2, -20.0, 20.0, 1.0);

  lat_constr.Fnum_b1.push_back(-ElemIndex("dl1a_1"));
  lat_constr.grad_dip_Fnum_b1.push_back(grad_dip_Fnum);

  grad_dip_Fnum.clear();
  grad_dip_Fnum.push_back(ElemIndex("dl2a_1"));
  grad_dip_Fnum.push_back(ElemIndex("dl2a_2"));
  grad_dip_Fnum.push_back(ElemIndex("dl2a_3"));
  grad_dip_Fnum.push_back(ElemIndex("dl2a_4"));
  grad_dip_Fnum.push_back(ElemIndex("dl2a_5"));
  if (dphi)
    prms.add_prm(grad_dip_Fnum, grad_dip_scl, -3, -20.0, 20.0, 1.0);
  if (long_grad_dip)
    prms.add_prm(grad_dip_Fnum, grad_dip_scl,  2, -20.0, 20.0, 1.0);

  lat_constr.Fnum_b1.push_back(-ElemIndex("dl2a_1"));
  lat_constr.grad_dip_Fnum_b1.push_back(grad_dip_Fnum);

  if (dphi) {
    // Dipole Cell.
    prms.add_prm("qf4", -3, -20.0, 20.0, 1.0);
    lat_constr.Fnum_b1.push_back(ElemIndex("qf4"));
    if (qf6_rb) {
      prms.add_prm("qf6", -3, -20.0, 20.0, 1.0);
      lat_constr.Fnum_b1.push_back(ElemIndex("qf6"));
    }
    prms.add_prm("qf8", -3, -20.0, 20.0, 1.0);
    lat_constr.Fnum_b1.push_back(ElemIndex("qf8"));

    // Commented out must be defined last.
    // prms.add_prm("dq1", -3, -20.0,   20.0,  1.0);
    lat_constr.Fnum_b1.push_back(ElemIndex("dq1"));
  }

  prms.add_prm("dq1",  2, -20.0, 20.0, 1.0);

  // Mid Straight.
  prms.add_prm("qf1", 2,   0.0, 20.0, 1.0);
  prms.add_prm("qd2", 2, -20.0,  0.0, 1.0);

  // Dipole Cell.
  prms.add_prm("qd3", 2, -20.0,  0.0, 1.0);
  prms.add_prm("qf4", 2,   0.0, 20.0, 1.0);
  prms.add_prm("qd5", 2, -20.0,  0.0, 1.0);

  // Standard Straight.
  // prms.add_prm("qf6", -1, -0.3,  0.01, 1.0);
  prms.add_prm("qf6",  2,  0.0, 20.0, 1.0);
  prms.add_prm("qf8",  2,  0.0, 20.0, 1.0);

 // Parameters are initialized in optimizer.

  // Lattice constraints are: alpha_x,y, beta_x,y, eta_x, eta'_x.
  if (relaxed) {
    // constr.add_constr(Elem_GetPos(ElemIndex("dl1a_5"), 1)-1,
    // 		      0e0, 0e0, 0e0, 0e0, 1e7, 1e7,
    // 		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    // constr.add_constr(Elem_GetPos(ElemIndex("dl1a_5"), 2),
    // 		      0e0, 0e0, 0e0, 0e0, 1e7, 1e7,
    // 		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("qf1"), 1)-1,
		      0e0, 0e0, 0e0, 0e0, 1e7, 1e7,
		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("qf1"), 2),
		      0e0, 0e0, 0e0, 0e0, 1e7, 1e7,
		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    // Default.
#if 0
    constr.add_constr(Elem_GetPos(ElemIndex("ms"), 1),
    		      1e7, 1e7, 1e-2, 1e-2, 1e7,   0e0,
    		      0.0, 0.0, 3.0,  1.5,  0.024, 0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("ss"), 1),
		      1e7, 1e7, 1e-2, 1e-2, 1e7, 1e7,
		      0.0, 0.0, 4.0,  2.5,  0.0, 0.0);
#else
    constr.add_constr(Elem_GetPos(ElemIndex("ms"), 1),
    		      1e7, 1e7, 5e2, 5e2, 1e7,   0e0,
    		      0.0, 0.0, 3.0, 1.5, 0.010, 0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("ss"), 1),
		      1e7, 1e7, 5e2, 5e2, 1e7, 1e7,
		      0.0, 0.0, 4.0, 2.5, 0.0, 0.0);
#endif
  } else {
    constr.add_constr(Elem_GetPos(ElemIndex("dl1a_5"), 1)-1,
		      0e0, 0e0, 0e0, 0e0, 1e7, 1e7,
		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("dl1a_5"), 2),
		      0e0, 0e0, 0e0, 0e0, 1e7, 1e7,
		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("ms"), 1),
		      1e7, 1e7, 1e0, 1e0, 1e7,   0e0,
		      0.0, 0.0, 3.0, 1.5, 0.024, 0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("ss"), 1),
		      1e7, 1e7, 1e0, 1e0, 1e7, 1e7,
		      0.0, 0.0, 4.0, 2.5, 0.0, 0.0);
  }

  if (false) {
    // Increase beta_x.
    constr.add_constr(Elem_GetPos(ElemIndex("sf1"), 1),
		      0e5, 0e5, 1e3,  1e1, 0e7, 0e7,
		      0.0, 0.0, 12.0, 1.0, 0.0, 0.0);
  }
  if (false) {
    // Increase beta_y.
    constr.add_constr(Elem_GetPos(ElemIndex("sd2"), 1),
		      0e5, 0e5, 1e1, 1e3, 0e7, 0e7,
		      0.0, 0.0, 1.0, 8.0, 0.0, 0.0);
  }

  lat_prms.bn_tol = 1e-5; lat_prms.step = 1.0;

  lat_constr.Fnum_b3.push_back(ElemIndex("sd1"));
  lat_constr.Fnum_b3.push_back(ElemIndex("sf1"));
  lat_constr.Fnum_b3.push_back(ElemIndex("sd2"));
  // lat_constr.Fnum_b3.push_back(ElemIndex("sh2"));

  lat_constr.eps0_x = eps0_x;

  mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 1));
  mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 2));
  lat_constr.mI_loc.push_back(mI_loc);

  for (k = 0; k < 2; k++)
    lat_constr.high_ord_achr_nu[k] = high_ord_achr_nu[k];

  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ss"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ss"), 2));

  n = lat_constr.high_ord_achr_Fnum.size() - 1;
  lat_constr.high_ord_achr_dnu.resize(n);
  for (k = 0; k < n; k++)
    lat_constr.high_ord_achr_dnu[k].resize(2, 0e0);

  for (k = 0; k < 2; k++)
    lat_constr.mI0[k] = mI_nu_ref[k];

  if (relaxed) {
    lat_constr.eps_x_scl              = 5e6;
    lat_constr.ksi1_ctrl_scl[0]       = 1e0;
#if 1
    // Default.
    // lat_constr.ksi1_ctrl_scl[1]       = 1e1;
    // lat_constr.ksi1_ctrl_scl[2]       = 1e1;
    lat_constr.ksi1_ctrl_scl[1]       = 1e0;
    lat_constr.ksi1_ctrl_scl[2]       = 1e0;
#else
    lat_constr.ksi1_ctrl_scl[1]       = 1e2;
    lat_constr.ksi1_ctrl_scl[2]       = 1e3;
#endif
#if 0
    // Default.
    lat_constr.drv_terms_simple_scl   = 1e-4;
#else
    // lat_constr.drv_terms_simple_scl   = 1e-2;
    lat_constr.drv_terms_simple_scl   = 5e-2;
    // lat_constr.drv_terms_simple_scl   = 1e-1;
#endif
    // Not useful.
    lat_constr.ksi1_svd_scl           = 0e3;
    // Default.
#if 0
    lat_constr.mI_scl[X_]             = 1e7;
    lat_constr.mI_scl[Y_]             = 1e7;
#else
    lat_constr.mI_scl[X_]             = 1e5;
    lat_constr.mI_scl[Y_]             = 1e5;
#endif
#if 0
    lat_constr.high_ord_achr_scl      = 1e7;
#else
    lat_constr.high_ord_achr_scl      = 1e4;
#endif
#if 1
    // lat_constr.alpha_c_scl            = 5e-7;
    lat_constr.alpha_c_scl            = 5e-6;
#else
    lat_constr.alpha_c_scl            = 1e-7;
#endif
  } else {
    lat_constr.eps_x_scl              = 5e6;
    lat_constr.ksi1_ctrl_scl[0]       = 0e-1;
    lat_constr.ksi1_ctrl_scl[1]       = 0e-2;
    lat_constr.ksi1_ctrl_scl[2]       = 0e-1;
    lat_constr.drv_terms_simple_scl   = 1e-3;
    // lat_constr.drv_terms_simple_scl = 1e-4;
    // Not useful.
    lat_constr.ksi1_svd_scl           = 0e3;
    lat_constr.mI_scl[X_]             = 1e5;
    lat_constr.mI_scl[Y_]             = 1e5;
    lat_constr.high_ord_achr_scl      = 1e5;
    lat_constr.alpha_c_scl            = 5e-7;
  }
  lat_constr.phi_scl                  = (dphi)? 1e0 : 0e0;

  prt_prms(lat_constr);

  // Super Period.
  lat_constr.phi0 = 15.0;
  lat_constr.L_scl = 0e-10; lat_constr.L0 = 10.0;

  lat_constr.ini_constr(true);
}


void opt_mI_sp(param_type &prms, constr_type &constr)
{
  // Parameter Type:
  //   Bend Angle  -3,
  //   Length      -2,
  //   Position    -1,
  //   Quadrupole   2.
  int                 k, n;
  std::vector<int>    grad_dip_Fnum, mI_loc;
  std::vector<double> grad_dip_scl;

  const bool
    dphi          = !false,
    long_grad_dip = !false,
    dip_cell      = true;

  printf("\n opt_mI_sp:\n relaxed: %d\n", relaxed);

  // Standard Cell.
  grad_dip_scl.push_back(0.129665);
  grad_dip_scl.push_back(0.149256);
  grad_dip_scl.push_back(0.174737);
  grad_dip_scl.push_back(0.216027);
  grad_dip_scl.push_back(0.330314);

  grad_dip_Fnum.push_back(ElemIndex("dl1a_1"));
  grad_dip_Fnum.push_back(ElemIndex("dl1a_2"));
  grad_dip_Fnum.push_back(ElemIndex("dl1a_3"));
  grad_dip_Fnum.push_back(ElemIndex("dl1a_4"));
  grad_dip_Fnum.push_back(ElemIndex("dl1a_5"));
  if (dphi)
    prms.add_prm(grad_dip_Fnum, grad_dip_scl, -3, -20.0, 20.0, 1.0);
  if (long_grad_dip)
    prms.add_prm(grad_dip_Fnum, grad_dip_scl,  2, -20.0, 20.0, 1.0);

  lat_constr.Fnum_b1.push_back(-ElemIndex("dl1a_1"));
  lat_constr.grad_dip_Fnum_b1.push_back(grad_dip_Fnum);

  if (true) {
    grad_dip_Fnum.clear();
    grad_dip_Fnum.push_back(ElemIndex("dl2a_1"));
    grad_dip_Fnum.push_back(ElemIndex("dl2a_2"));
    grad_dip_Fnum.push_back(ElemIndex("dl2a_3"));
    grad_dip_Fnum.push_back(ElemIndex("dl2a_4"));
    grad_dip_Fnum.push_back(ElemIndex("dl2a_5"));
    if (dphi)
      prms.add_prm(grad_dip_Fnum, grad_dip_scl, -3, -20.0, 20.0, 1.0);
    if (long_grad_dip)
      prms.add_prm(grad_dip_Fnum, grad_dip_scl,  2, -20.0, 20.0, 1.0);

    lat_constr.Fnum_b1.push_back(-ElemIndex("dl2a_1"));
    lat_constr.grad_dip_Fnum_b1.push_back(grad_dip_Fnum);
  }

  if (dphi) {
    // Dipole Cell.
    prms.add_prm("qf4", -3, -20.0, 20.0, 1.0);
    lat_constr.Fnum_b1.push_back(ElemIndex("qf4"));
    if (qf6_rb) {
      prms.add_prm("qf6", -3, -20.0, 20.0, 1.0);
      lat_constr.Fnum_b1.push_back(ElemIndex("qf6"));
    }
    prms.add_prm("qf8", -3, -20.0, 20.0, 1.0);
    lat_constr.Fnum_b1.push_back(ElemIndex("qf8"));

    // Commented out must be defined last.
    // prms.add_prm("dq1", -3, -20.0,   20.0,  1.0);
    lat_constr.Fnum_b1.push_back(ElemIndex("dq1"));
  }

  if (dip_cell) {
    prms.add_prm("dq1", 2, -20.0,   20.0,  1.0);

    // Mid Straight.
    prms.add_prm("qf1", 2, -20.0, 20.0, 1.0);
    prms.add_prm("qd2", 2, -20.0, 20.0, 1.0);

    // Dipole Cell.
    prms.add_prm("qd3", 2, -20.0, 20.0, 1.0);
    prms.add_prm("qf4", 2, -20.0, 20.0, 1.0);
    prms.add_prm("qd5", 2, -20.0, 20.0, 1.0);

    // Standard Straight.
    prms.add_prm("qf6", 2, -20.0, 20.0, 1.0);
    prms.add_prm("qf8", 2, -20.0, 20.0, 1.0);
  }

  // Long Straight.
  // Perturbation of Dipole Cell.
  prms.add_prm("qd3_c1",   2, -20.0, 20.0, 1.0);

  prms.add_prm("qf1_c1",   2, -20.0, 20.0, 1.0);
  prms.add_prm("qd2_c1",   2, -20.0, 20.0, 1.0);
  prms.add_prm("quad_add", 2, -20.0, 20.0, 1.0);

 // Parameters are initialized in optimizer.

  // Lattice constraints are: alpha_x,y, beta_x,y, eta_x, eta'_x.
  if (relaxed) {
#if 0
    constr.add_constr(Elem_GetPos(ElemIndex("dl1a_5"), 1)-1,
    		      0e0, 0e0, 0e0, 0e0, 1e7, 1e7,
    		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("dl1a_5"), 2),
    		      0e0, 0e0, 0e0, 0e0, 1e7, 1e7,
    		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
#else
    // Include quadrupoles before Std & Long Straights.
    // constr.add_constr(Elem_GetPos(ElemIndex("qf1_c1"), 1)-1,
    // 		      0e0, 0e0, 0e0, 0e0, 1e7, 1e7,
    // 		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("quad_add"), 1)-1,
		      0e0, 0e0, 0e0, 0e0, 1e7, 1e7,
		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("qf1"), 1),
		      0e0, 0e0, 0e0, 0e0, 1e7, 1e7,
		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("qf1"), 2)-1,
		      0e0, 0e0, 0e0, 0e0, 1e7, 1e7,
		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("qf1"), 3),
		      0e0, 0e0, 0e0, 0e0, 1e7, 1e7,
		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("qf1"), 4)-1,
		      0e0, 0e0, 0e0, 0e0, 1e7, 1e7,
		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    // constr.add_constr(Elem_GetPos(ElemIndex("qf1_c1"), 2),
    // 		      0e0, 0e0, 0e0, 0e0, 1e7, 1e7,
    // 		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("quad_add"), 2),
		      0e0, 0e0, 0e0, 0e0, 1e7, 1e7,
		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
#endif

#if 0
    // Include constraint on alpha; in case of using ps_rot.
    constr.add_constr(Elem_GetPos(ElemIndex("ms"), 1),
		      1e6, 1e6, 1e-2, 1e-2, 1e6,   1e7,
		      0.0, 0.0, 3.0,  1.5,  0.024, 0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("ms"), 2),
		      1e6, 1e6, 1e-2, 1e-2, 1e6,   1e7,
		      0.0, 0.0, 3.0,  1.5,  0.024, 0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("ms"), 3),
		      1e6, 1e6, 1e-2, 1e-2, 1e6,   1e7,
		      0.0, 0.0, 3.0,  1.5,  0.024, 0.0);
    // Both SS constraints are needed.
    constr.add_constr(Elem_GetPos(ElemIndex("ss"), 1),
		      1e6, 1e6, 1e-2, 1e-2, 1e7, 1e7,
		      0.0, 0.0, 4.0,  2.5,  0.0, 0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("ss"), 2),
		      1e6, 1e6, 1e-2, 1e-2, 1e7, 1e7,
		      0.0, 0.0, 4.0,  2.5,  0.0, 0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("ls"), 1),
		      1e6, 1e6, 1e-2, 1e-2, 1e7, 1e7,
		      0.0, 0.0, 10.0, 4.0,  0.0, 0.0);
#else
    // Include constraint on alpha; in case of using ps_rot.
    constr.add_constr(Elem_GetPos(ElemIndex("ms"), 1),
		      1e6, 1e6, 1e2, 1e2, 5e6,   1e7,
		      0.0, 0.0, 3.0, 1.5, 0.015, 0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("ms"), 2),
		      1e6, 1e6, 1e2, 1e2, 5e6,   1e7,
		      0.0, 0.0, 3.0, 1.5, 0.015, 0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("ms"), 3),
		      1e6, 1e6, 1e2, 1e2, 5e6,   1e7,
		      0.0, 0.0, 3.0, 1.5, 0.015, 0.0);
    // Both SS constraints are needed.
    constr.add_constr(Elem_GetPos(ElemIndex("ss"), 1),
		      1e6, 1e6, 1e2,  1e2, 1e7, 1e7,
		      0.0, 0.0, 4.0,  2.5, 0.0, 0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("ss"), 2),
		      1e6, 1e6, 1e2,  1e2, 1e7, 1e7,
		      0.0, 0.0, 4.0,  2.5, 0.0, 0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("ls"), 1),
		      1e6, 1e6, 1e2,  1e2, 1e7, 1e7,
		      0.0, 0.0, 10.0, 4.0, 0.0, 0.0);
#endif
  } else {
    constr.add_constr(Elem_GetPos(ElemIndex("dl1a_5"), 1)-1,
		      0e0, 0e0, 0e0, 0e0, 1e8, 1e8,
		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("dl1a_5"), 2),
		      0e0, 0e0, 0e0, 0e0, 1e8, 1e8,
		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    // Include constraint on alpha; in case of using ps_rot.
    constr.add_constr(Elem_GetPos(ElemIndex("ms"), 1),
		      1e6, 1e6, 1e-1, 1e-1, 1e6,   1e8,
		      0.0, 0.0, 3.0, 1.5, 0.024, 0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("ms"), 2),
		      1e6, 1e6, 1e-1, 1e-1, 1e6,   1e8,
		      0.0, 0.0, 3.0, 1.5, 0.024, 0.0);
    // Both SS constraints are needed.
    constr.add_constr(Elem_GetPos(ElemIndex("ss"), 1),
		      1e6, 1e6, 1e-1, 1e-1, 1e8, 1e8,
		      0.0, 0.0, 4.0, 2.5, 0.0, 0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("ss"), 2),
		      1e6, 1e6, 1e-1, 1e-1, 1e8, 1e8,
		      0.0, 0.0, 4.0, 2.5, 0.0, 0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("ls"), 1),
		      1e6, 1e6, 1e-1,  1e-1, 1e8, 1e8,
		      0.0, 0.0, 10.0, 4.0, 0.0, 0.0);
  }

  if (false) {
    // Increase beta_x.
    constr.add_constr(Elem_GetPos(ElemIndex("sf1"), 1),
		      0e5, 0e5, 1e3,  1e2, 0e7, 0e7,
		      0.0, 0.0, 12.0, 1.0, 0.0, 0.0);
  }
  if (false) {
    // Increase beta_y.
    constr.add_constr(Elem_GetPos(ElemIndex("sd2"), 1),
		      0e5, 0e5, 1e1, 1e3, 0e7, 0e7,
		      0.0, 0.0, 1.0, 8.0, 0.0, 0.0);
  }

  lat_prms.bn_tol = 1e-5; lat_prms.step = 1.0;

  lat_constr.Fnum_b3.push_back(ElemIndex("sd1"));
  lat_constr.Fnum_b3.push_back(ElemIndex("sf1"));
  lat_constr.Fnum_b3.push_back(ElemIndex("sd2"));
  // lat_constr.Fnum_b3.push_back(ElemIndex("sh2"));

  lat_constr.eps0_x = eps0_x;

  mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 1));
  mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 2));
  lat_constr.mI_loc.push_back(mI_loc);
  mI_loc.clear();
  mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 3));
  mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 4));
  lat_constr.mI_loc.push_back(mI_loc);

  for (k = 0; k < 2; k++)
    lat_constr.high_ord_achr_nu[k] = high_ord_achr_nu[k];

  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ls"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ss"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ss"), 2));

  n = lat_constr.high_ord_achr_Fnum.size() - 1;
  lat_constr.high_ord_achr_dnu.resize(n);
  for (k = 0; k < n; k++)
    lat_constr.high_ord_achr_dnu[k].resize(2, 0e0);

  for (k = 0; k < 2; k++)
    lat_constr.mI0[k] = mI_nu_ref[k];

  if (relaxed) {
    // lat_constr.eps_x_scl            = 1e7;
    lat_constr.eps_x_scl            = 5e6;
    // lat_constr.eps_x_scl            = 1e6;
    lat_constr.ksi1_ctrl_scl[0]     = 0e-1;
    lat_constr.ksi1_ctrl_scl[1]     = 0e0;
    lat_constr.ksi1_ctrl_scl[2]     = 0e-1;
#if 0
    // Default.
    lat_constr.drv_terms_simple_scl   = 1e-4;
#else
    lat_constr.drv_terms_simple_scl   = 1e-2;
    // lat_constr.drv_terms_simple_scl   = 5e-2;
    // lat_constr.drv_terms_simple_scl   = 1e-1;
#endif
    // Not useful.
    lat_constr.ksi1_svd_scl         = 0e3;
#if 0
    lat_constr.mI_scl[X_]           = 1e7;
    lat_constr.mI_scl[Y_]           = 1e7;
#else
    lat_constr.mI_scl[X_]           = 1e5;
    lat_constr.mI_scl[Y_]           = 1e5;
#endif
#if 0
    lat_constr.high_ord_achr_scl    = 1e7;
#else
    lat_constr.high_ord_achr_scl    = 1e4;
    // lat_constr.high_ord_achr_scl    = 1e5;
#endif
#if 1
    // 1e-7 is too small.
    // lat_constr.alpha_c_scl            = 5e-8;
    // lat_constr.alpha_c_scl            = 1e-7;
    lat_constr.alpha_c_scl            = 5e-7;
    // lat_constr.alpha_c_scl            = 5e-6;
#else
    lat_constr.alpha_c_scl            = 1e-7;
#endif
  } else {
    lat_constr.eps_x_scl            = 1e6;
    lat_constr.ksi1_ctrl_scl[0]     = 0e-1;
    lat_constr.ksi1_ctrl_scl[1]     = 0e0;
    lat_constr.ksi1_ctrl_scl[2]     = 0e-1;
#if 0
    // Default.
    lat_constr.drv_terms_simple_scl   = 1e-4;
#else
    // lat_constr.drv_terms_simple_scl   = 1e-2;
    lat_constr.drv_terms_simple_scl   = 5e-2;
    // lat_constr.drv_terms_simple_scl   = 1e-1;
#endif
    // Not useful.
    lat_constr.ksi1_svd_scl         = 0e3;
    lat_constr.mI_scl[X_]           = 1e4;
    lat_constr.mI_scl[Y_]           = 1e4;
    lat_constr.high_ord_achr_scl    = 1e4;
    lat_constr.alpha_c_scl          = 1e-9;
  }
  lat_constr.phi_scl                = (dphi)? 1e0 : 0e0;

  prt_prms(lat_constr);

  if (!sp_short)
    // Super Period.
    lat_constr.phi0 = 60.0;
  else
    // Short Super Period.
    lat_constr.phi0 = 45.0;

  lat_constr.L_scl = 0e-10; lat_constr.L0 = 10.0;

  lat_constr.ini_constr(true);
}


void match_ls(param_type &prms, constr_type &constr)
{
  // Parameter Type:
  //   Bend Angle  -3,
  //   Length      -2,
  //   Position    -1,
  //   Quadrupole   2.
  int                 j, k, n;
  std::vector<int>    grad_dip_Fnum;
  std::vector<double> grad_dip_scl;

  if (pert_dip_cell)
    // Perturbed symmetry at end of Dipole Cell: 
    //   1. Initialize with: [Qf1, Qd2, Qd3].
    //   2. Exclude for 1st pass.
    //   3. Include for fine tuning.
    prms.add_prm("qd3_c1",   2, -20.0, 20.0, 1.0);
  // Long Straight.
  prms.add_prm("qd2_c1",   2, -20.0, 20.0, 1.0);
  prms.add_prm("qf1_c1",   2, -20.0, 20.0, 1.0);
  prms.add_prm("quad_add", 2, -20.0, 20.0, 1.0);

  // Parameters are initialized in optimizer.

  if (!pert_dip_cell)
    // Eta_x & eta_x' are zero after Dipole Cell.
    constr.add_constr(Elem_GetPos(ElemIndex("ls"), 1),
		      1e2, 1e2, 1e-2, 1e-3, 1e-10, 1e-10,
		      0.0, 0.0, 10.0, 4.0,  0.0,   0.0);
  else {
    constr.add_constr(Elem_GetPos(ElemIndex("ls"), 1),
		      1e2, 1e2, 1e-3, 1e-3, 1e3, 1e3,
		      0.0, 0.0, 10.0, 4.0,  0.0, 0.0);
    // constr.add_constr(Elem_GetPos(ElemIndex("qf1_c1"), 1),
    // 		      0e0, 0e0, 0e0, 0e0, 1e3, 1e3,
    // 		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("quad_add"), 1),
		      0e0, 0e0, 0e0, 0e0, 1e3, 1e3,
		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  }

  lat_prms.bn_tol = 1e-5; lat_prms.step = 1.0;

  for (k = 0; k < 2; k++)
    lat_constr.high_ord_achr_nu[k] = high_ord_achr_nu[k]/2e0;

  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ms"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ls"), 1));

  n = lat_constr.high_ord_achr_Fnum.size() - 1;
  lat_constr.high_ord_achr_dnu.resize(n);
  for (k = 0; k < n; k++)
    lat_constr.high_ord_achr_dnu[k].resize(2, 0e0);

  if (!pert_dip_cell)
    lat_constr.high_ord_achr_scl = 1e-2;
  else
    lat_constr.high_ord_achr_scl = 1e-2;

  lat_constr.ini_constr(false);

  for (j = 0; j < n_ic; j++)
    for (k = 0; k < 2; k++)
      lat_constr.ic[j][k] = ic[j][k];
}


void match_als_u(param_type &prms, constr_type &constr)
{
  // Parameter Type:
  //   Bend Angle  -3,
  //   Length      -2,
  //   Position    -1,
  //   Quadrupole   2.
  int                 j, k, n;
  std::vector<int>    grad_dip_Fnum;
  std::vector<double> grad_dip_scl;

  // From Center of Mid Straight: alpha, beta, eta, eta'.
  const int    n_ic        = 4;
  const double ic[n_ic][2] =
    {{0.0, 0.0}, {0.1153496843, 3.2620137331}, {0.0002943186, 0.0}, {0.0, 0.0}};
 
  prms.add_prm("q1",   2, -20.0, 20.0, 1.0);
  prms.add_prm("q2",   2, -20.0, 20.0, 1.0);
  prms.add_prm("q3",   2, -20.0, 20.0, 1.0);
  prms.add_prm("q4",   2, -20.0, 20.0, 1.0);

  prms.add_prm("mqfh", 2, -20.0, 20.0, 1.0);

  // Parameters are initialized in optimizer.

  constr.add_constr(Elem_GetPos(ElemIndex("q1"), 1),
  		    0e0, 0e0, 0e0, 0e0, 1e3, 1e3,
  		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  constr.add_constr(Elem_GetPos(ElemIndex("ssh"), 1),
  		    1e0, 1e0, 1e-5, 1e-5, 1e5, 1e5,
  		    0.0, 0.0, 2.1,  2.3,  0.0, 0.0);

  lat_prms.bn_tol = 1e-5; lat_prms.step = 1.0;

  for (k = 0; k < 2; k++)
    lat_constr.high_ord_achr_nu[k] = high_ord_achr_nu[k];

  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ssh"), 1));

  n = lat_constr.high_ord_achr_Fnum.size() - 1;
  lat_constr.high_ord_achr_dnu.resize(n);
  for (k = 0; k < n; k++)
    lat_constr.high_ord_achr_dnu[k].resize(2, 0e0);

  lat_constr.ksi1_svd_scl      = 0e0;
  lat_constr.high_ord_achr_scl = 0e-1;

  lat_constr.ini_constr(false);

  for (j = 0; j < n_ic; j++)
    for (k = 0; k < 2; k++)
      lat_constr.ic[j][k] = ic[j][k];
}


void opt_tba(param_type &prms, constr_type &constr)
{
  // Parameter Type:
  //   Bend Angle  -3,
  //   Length      -2,
  //   Position    -1,
  //   Quadrupole   2.
  int k, n;

  // TBA Cell.
  prms.add_prm("b1",       -3, -20.0,   20.0,  1.0);
  prms.add_prm("b1",       -2,   0.1,    1.2,  1.0);
  // prms.add_prm("b1",       -1,   0.075,  0.07, 1.0);
  prms.add_prm("b1",        2, -20.0,   20.0,  1.0);
  prms.add_prm("b2",       -2,   0.1,    0.7,  1.0);
  prms.add_prm("b2",        2, -20.0,   20.0,  1.0);
  prms.add_prm("qf1a_tba",  2, -20.0,   20.0,  1.0);
  prms.add_prm("qf1a_tba", -2,   0.1,    0.25, 1.0);
  prms.add_prm("qf1b_tba",  2, -20.0,   20.0,  1.0);
  prms.add_prm("qf1b_tba", -2,   0.1,    0.25, 1.0);

  // Mid-Straight.
  prms.add_prm("d8",       -2,   0.075,  0.2,  1.0);
  prms.add_prm("qf1_ms",    2, -20.0,   20.0,  1.0);
  prms.add_prm("qf1_ms",   -2,   0.1,    0.25, 1.0);
  prms.add_prm("d9",       -2,   0.075,  0.2,  1.0);
  prms.add_prm("qd2_ms",    2, -20.0,   20.0,  1.0);
  prms.add_prm("qd2_ms",   -2,   0.1,    0.2,  1.0);

 // Parameters are initialized in optimizer.

  // Lattice constraints are: alpha_x,y, beta_x,y, eta_x, eta'_x.
  constr.add_constr(Elem_GetPos(ElemIndex("b1"), 2),
  		    0e0, 0e0, 0e0, 0e0, 1e5, 1e5,
  		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

  // Include constraint on alpha; in case of using ps_rot.
  constr.add_constr(Elem_GetPos(ElemIndex("ms"), 1),
  		    1e4, 1e4, 1e-10, 1e-10, 0e0, 0e0,
  		    0.0, 0.0, 3.0,   1.5,   0.0, 0.0);

  lat_prms.bn_tol = 1e-6; lat_prms.step = 1.0;

  lat_constr.Fnum_b3.push_back(ElemIndex("sd1a"));
  lat_constr.Fnum_b3.push_back(ElemIndex("sfh"));
  lat_constr.Fnum_b3.push_back(ElemIndex("sd1b"));

  lat_constr.Fnum_b1.push_back(ElemIndex("b1"));
  lat_constr.Fnum_b1.push_back(ElemIndex("b2"));

  lat_constr.high_ord_achr_scl = 1e3;
  for (k = 0; k < 2; k++)
    lat_constr.high_ord_achr_nu[k] = high_ord_achr_nu[k];

  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ms"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ms"), 2));

  n = lat_constr.high_ord_achr_Fnum.size() - 1;
  lat_constr.high_ord_achr_dnu.resize(n);
  for (k = 0; k < n; k++)
    lat_constr.high_ord_achr_dnu[k].resize(2, 0e0);

  lat_constr.eps_x_scl = 1e2; lat_constr.eps0_x = 0.190;

  lat_constr.mI_scl[X_] = 1e-10; lat_constr.mI_scl[Y_] = 1e-10;
  for (k = 0; k < 2; k++)
    lat_constr.mI0[k] = mI_nu_ref[k];

  // 2 TBA & 1 MS: phi = 7.5.
  lat_constr.phi_scl = 1e0; lat_constr.phi0 = 7.5;
  lat_constr.L_scl = 1e-10; lat_constr.L0 = 10.0;

  lat_constr.ini_constr(true);
}


void match_ms(param_type &prms, constr_type &constr)
{
  // Parameter Type:
  //   Bend Angle  -3,
  //   Length      -2,
  //   Position    -1,
  //   Quadrupole   2.
  int                 j, k;
  std::vector<int>    grad_dip_Fnum;
  std::vector<double> grad_dip_scl;

  // From Center of Mid Straight: alpha, beta, eta, eta'.
  const int    n_ic        = 4;
  const double ic[n_ic][2] =
    {{0.0, 0.0}, {7.08276, 3.03915}, {0.0, 0.0}, {0.0, 0.0}};
 
  // Standard Straight.
  grad_dip_scl.push_back(0.129665);
  grad_dip_scl.push_back(0.149256);
  grad_dip_scl.push_back(0.174737);
  grad_dip_scl.push_back(0.216027);
  grad_dip_scl.push_back(0.330314);

  grad_dip_Fnum.push_back(ElemIndex("bl1_1"));
  grad_dip_Fnum.push_back(ElemIndex("bl1_2"));
  grad_dip_Fnum.push_back(ElemIndex("bl1_3"));
  grad_dip_Fnum.push_back(ElemIndex("bl1_4"));
  grad_dip_Fnum.push_back(ElemIndex("bl1_5"));
  prms.add_prm(grad_dip_Fnum, grad_dip_scl, 2, -20.0, 20.0, 1.0);

  grad_dip_Fnum.clear();
  grad_dip_Fnum.push_back(ElemIndex("bl2_1"));
  grad_dip_Fnum.push_back(ElemIndex("bl2_2"));
  grad_dip_Fnum.push_back(ElemIndex("bl2_3"));
  grad_dip_Fnum.push_back(ElemIndex("bl2_4"));
  grad_dip_Fnum.push_back(ElemIndex("bl2_5"));
  prms.add_prm(grad_dip_Fnum, grad_dip_scl, 2, -20.0, 20.0, 1.0);

  // Parameters are initialized in optimizer.

  constr.add_constr(Elem_GetPos(ElemIndex("ms"), 1),
  		    1e0, 1e0, 1e-5, 1e-5, 1e0, 1e0,
  		    0.0, 0.0, 3.0,  1.5,  0.0, 0.0);

  lat_prms.bn_tol = 1e-6; lat_prms.step = 1.0;

  lat_constr.high_ord_achr_scl = 0e0;
  for (k = 0; k < 2; k++)
    lat_constr.high_ord_achr_nu[k] = high_ord_achr_nu[k]/2.0;

  lat_constr.ini_constr(false);

  for (j = 0; j < n_ic; j++)
    for (k = 0; k < 2; k++)
      lat_constr.ic[j][k] = ic[j][k];
}


void opt_std_cell(param_type &prms, constr_type &constr)
{
  // Parameter Type:
  //   Bend Angle  -3,
  //   Length      -2,
  //   Position    -1,
  //   Quadrupole   2.
  int k, n;

  // TBA Cell.
  prms.add_prm("b1",       -3, -20.0,   20.0,  1.0);
  prms.add_prm("b1",       -2,   0.1,    1.2,  1.0);
  // prms.add_prm("b1",       -1,   0.075,  0.07, 1.0);
  prms.add_prm("b1",        2, -20.0,   20.0,  1.0);
  prms.add_prm("b2",       -2,   0.1,    0.7,  1.0);
  prms.add_prm("b2",        2, -20.0,   20.0,  1.0);
  prms.add_prm("qf1a_tba",  2, -20.0,   20.0,  1.0);
  prms.add_prm("qf1a_tba", -2,   0.1,    0.25, 1.0);
  prms.add_prm("qf1b_tba",  2, -20.0,   20.0,  1.0);
  prms.add_prm("qf1b_tba", -2,   0.1,    0.25, 1.0);

  // Mid-Straight.
  prms.add_prm("d8",       -2,   0.075,  0.2,  1.0);
  prms.add_prm("qf1_ms",    2, -20.0,   20.0,  1.0);
  prms.add_prm("qf1_ms",   -2,   0.1,    0.25, 1.0);
  prms.add_prm("d9",       -2,   0.075,  0.2,  1.0);
  prms.add_prm("qd2_ms",    2, -20.0,   20.0,  1.0);
  prms.add_prm("qd2_ms",   -2,   0.1,    0.2,  1.0);

  // Standard Straight.
  prms.add_prm("d10",    -2,   0.075,  0.4, 1.0);
  prms.add_prm("qd3_ss",  2, -20.0,   20.0, 1.0);
  prms.add_prm("qd3_ss", -2,   0.1,    0.4, 1.0);
  prms.add_prm("d11",    -2,   0.075,  0.4, 1.0);
  prms.add_prm("qf2_ss",  2, -20.0,   20.0, 1.0);
  prms.add_prm("qf2_ss", -2,   0.15,   0.4, 1.0);
  prms.add_prm("d12",    -2,   0.075,  0.4, 1.0);
  prms.add_prm("qd1_ss",  2, -20.0,   20.0, 1.0);
  prms.add_prm("qd1_ss", -2,   0.1,    0.4, 1.0);

  // Parameters are initialized in optimizer.

  // Lattice constraints are: alpha_x,y, beta_x,y, eta_x, eta'_x.
  constr.add_constr(Elem_GetPos(ElemIndex("b2"), 1),
  		    1e3, 1e3, 0e0, 0e-1, 0e0, 1e6,
  		    0.0, 0.0, 0.0, 10.0, 0.0, 0.0);
  constr.add_constr(Elem_GetPos(ElemIndex("b1"), 1)-1,
  		    0e0, 0e0, 0e0, 0e0, 1e6, 1e6,
  		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  constr.add_constr(Elem_GetPos(ElemIndex("b1"), 2),
  		    0e0, 0e0, 0e0, 0e0, 1e6, 1e6,
  		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

  constr.add_constr(Elem_GetPos(ElemIndex("ms"), 1),
  		    1e3, 1e3, 1e-2, 1e-2, 0e0, 0e0,
  		    0.0, 0.0, 3.0,  1.5,  0.0, 0.0);

  constr.add_constr(Elem_GetPos(ElemIndex("ss"), 1),
  		    1e3, 1e3, 1e-2, 1e-2, 0e0, 0e0,
  		    0.0, 0.0, 4.0,  2.5,  0.0, 0.0);

 lat_prms.bn_tol = 1e-6; lat_prms.step = 1.0;

  lat_constr.Fnum_b3.push_back(ElemIndex("sd1a"));
  lat_constr.Fnum_b3.push_back(ElemIndex("sfh"));
  lat_constr.Fnum_b3.push_back(ElemIndex("sd1b"));

  lat_constr.Fnum_b1.push_back(ElemIndex("b1"));
  lat_constr.Fnum_b1.push_back(ElemIndex("b2"));

  lat_constr.high_ord_achr_scl = 1e3;
  for (k = 0; k < 2; k++)
    lat_constr.high_ord_achr_nu[k] = high_ord_achr_nu[k]/2.0;

  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ms"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("b2"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ss"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("b2"), 3));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ms"), 2));

  n = lat_constr.high_ord_achr_Fnum.size() - 1;
  lat_constr.high_ord_achr_dnu.resize(n);
  for (k = 0; k < n; k++)
    lat_constr.high_ord_achr_dnu[k].resize(2, 0e0);

  lat_constr.mI_scl[X_] = 1e-10; lat_constr.mI_scl[Y_] = 1e-10;
  for (k = 0; k < 2; k++)
    lat_constr.mI0[k] = mI_nu_ref[k];

  lat_constr.eps_x_scl = 1e3; lat_constr.eps0_x = 0.190;

  // 2 TBA: phi = 15.
  lat_constr.phi_scl = 1e0; lat_constr.phi0 = 15.0;

  lat_constr.L_scl = 1e-10; lat_constr.L0 = 10.0;

  lat_constr.ini_constr(true);
}


void opt_long_cell(param_type &prms, constr_type &constr)
{
  // Parameter Type:
  //   Bend Angle  -3,
  //   Length      -2,
  //   Position    -1,
  //   Quadrupole   2.
  int k, n;

  // Long Straight.
  prms.add_prm("d4",     -2,   0.075,  0.4, 1.0);
  prms.add_prm("qd4_ls",  2, -20.0,    0.0, 1.0);
  prms.add_prm("qd4_ls", -2,   0.1,    0.4, 1.0);
  prms.add_prm("d3",     -2,   0.075,  0.4, 1.0);
  prms.add_prm("qf3_ls",  2,   0.0,   20.0, 1.0);
  prms.add_prm("qf3_ls", -2,   0.1,    0.4, 1.0);
  prms.add_prm("d2",     -2,   0.075,  0.4, 1.0);
  prms.add_prm("qd2_ls",  2, -20.0,    0.0, 1.0);
  prms.add_prm("qd2_ls", -2,   0.1,    0.4, 1.0);
  prms.add_prm("d1",     -2,   0.075,  0.2, 1.0);
  prms.add_prm("qf1_ls",  2,   0.0,   20.0, 1.0);
  prms.add_prm("qf1_ls", -2,   0.1,    0.4, 1.0);

  // Parameters are initialized in optimizer.

  // Lattice constraints are: alpha_x,y, beta_x,y, eta_x, eta'_x.
  constr.add_constr(Elem_GetPos(ElemIndex("ls"), 1),
  		    1e0, 1e0, 1e-4, 1e-4, 0e0, 0e0,
  		    0.0, 0.0, 15.0, 4.0,  0.0, 0.0);

  lat_prms.bn_tol = 1e-6; lat_prms.step = 1.0;

  lat_prms.bn_tol = 1e-6; lat_prms.step = 1.0;

  lat_constr.Fnum_b3.push_back(ElemIndex("sd1a"));
  lat_constr.Fnum_b3.push_back(ElemIndex("sfh"));
  lat_constr.Fnum_b3.push_back(ElemIndex("sd1b"));

  lat_constr.Fnum_b1.push_back(ElemIndex("b1"));
  lat_constr.Fnum_b1.push_back(ElemIndex("b2"));

  lat_constr.high_ord_achr_scl = 1e2;
  for (k = 0; k < 2; k++)
    lat_constr.high_ord_achr_nu[k] = high_ord_achr_nu[k]/2.0;

  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ls"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("b2"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ms"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("b2"), 3));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ls"), 2));

  n = lat_constr.high_ord_achr_Fnum.size() - 1;
  lat_constr.high_ord_achr_dnu.resize(n);
  for (k = 0; k < n; k++)
    lat_constr.high_ord_achr_dnu[k].resize(2, 0e0);

  lat_constr.mI_scl[X_] = 1e-10; lat_constr.mI_scl[Y_] = 1e-10;
  for (k = 0; k < 2; k++)
    lat_constr.mI0[k] = mI_nu_ref[k];

  lat_constr.eps_x_scl = 1e3; lat_constr.eps0_x = 0.190;

  // 2 TBA: phi = 15.
  lat_constr.phi_scl = 1e0; lat_constr.phi0 = 15.0;

  lat_constr.L_scl = 1e-10; lat_constr.L0 = 10.0;

  lat_constr.ini_constr(true);
}


void opt_short_cell(param_type &prms, constr_type &constr)
{
  // Parameter Type:
  //   Bend Angle  -3,
  //   Length      -2,
  //   Position    -1,
  //   Quadrupole   2.
  int k, n;

  // TBA Cell.
  prms.add_prm("b1",       -3, -20.0,   20.0,  1.0);
  prms.add_prm("b1",       -2,   0.1,    1.2,  1.0);
  // prms.add_prm("b1",       -1,   0.075,  0.07, 1.0);
  prms.add_prm("b1",        2, -20.0,   20.0,  1.0);
  prms.add_prm("b2",       -2,   0.1,    0.7,  1.0);
  prms.add_prm("b2",        2, -20.0,   20.0,  1.0);
  prms.add_prm("qf1a_tba",  2, -20.0,   20.0,  1.0);
  prms.add_prm("qf1a_tba", -2,   0.1,    0.25, 1.0);
  prms.add_prm("qf1b_tba",  2, -20.0,   20.0,  1.0);
  prms.add_prm("qf1b_tba", -2,   0.1,    0.25, 1.0);

  // Mid-Straight.
  prms.add_prm("d8",       -2,   0.075,  0.2,  1.0);
  prms.add_prm("qf1_ms",    2, -20.0,   20.0,  1.0);
  prms.add_prm("qf1_ms",   -2,   0.1,    0.25, 1.0);
  prms.add_prm("d9",       -2,   0.075,  0.2,  1.0);
  prms.add_prm("qd2_ms",    2, -20.0,   20.0,  1.0);
  prms.add_prm("qd2_ms",   -2,   0.1,    0.2,  1.0);

  prms.add_prm("d_ms",     -2,   1.2,     1.7, 1.0);

  // Standard Straight.
  prms.add_prm("d10",    -2,   0.075,  0.4, 1.0);
  prms.add_prm("qd3_ss",  2, -20.0,   20.0, 1.0);
  prms.add_prm("qd3_ss", -2,   0.1,    0.4, 1.0);
  prms.add_prm("d11",    -2,   0.075,  0.4, 1.0);
  prms.add_prm("qf2_ss",  2, -20.0,   20.0, 1.0);
  prms.add_prm("qf2_ss", -2,   0.15,   0.4, 1.0);
  prms.add_prm("d12",    -2,   0.075,  0.4, 1.0);
  prms.add_prm("qd1_ss",  2, -20.0,   20.0, 1.0);
  prms.add_prm("qd1_ss", -2,   0.1,    0.4, 1.0);

  prms.add_prm("d_ss",    -2,   2.1,     2.6, 1.0);

  // Long Straight.
  prms.add_prm("d4",     -2,   0.075,  0.4, 1.0);
  prms.add_prm("qd4_ls",  2,   0.0,   20.0, 1.0);
  prms.add_prm("qd4_ls", -2,   0.1,    0.4, 1.0);
  prms.add_prm("d3",     -2,   0.075,  0.4, 1.0);
  prms.add_prm("qf3_ls",  2, -20.0,    0.0, 1.0);
  prms.add_prm("qf3_ls", -2,   0.1,    0.4, 1.0);
  prms.add_prm("d2",     -2,   0.075,  0.4, 1.0);
  prms.add_prm("qd2_ls",  2,   0.0,   20.0, 1.0);
  prms.add_prm("qd2_ls", -2,   0.1,    0.4, 1.0);
  prms.add_prm("d1",     -2,   0.075,  0.2, 1.0);
  prms.add_prm("qf1_ls",  2, -20.0,    0.0, 1.0);
  prms.add_prm("qf1_ls", -2,   0.1,    0.4, 1.0);

  prms.add_prm("d_ls",   -2,   3.7,    4.1,  1.0);

  // Parameters are initialized in optimizer.

  // Lattice constraints are: alpha_x,y, beta_x,y, eta_x, eta'_x.
  constr.add_constr(Elem_GetPos(ElemIndex("b2"), 1),
  		    1e4, 1e4, 0e0, 0e-1, 0e0, 1e5,
  		    0.0, 0.0, 0.0, 10.0, 0.0, 0.0);
  constr.add_constr(Elem_GetPos(ElemIndex("b1"), 1)-1,
  		    0e0, 0e0, 0e0, 0e0, 1e5, 1e5,
  		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  constr.add_constr(Elem_GetPos(ElemIndex("b1"), 2),
  		    0e0, 0e0, 0e0, 0e0, 1e5, 1e5,
  		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

  constr.add_constr(Elem_GetPos(ElemIndex("ms"), 1),
  		    1e4, 1e4, 1e-2, 1e-2, 0e0, 0e0,
  		    0.0, 0.0, 3.0,  1.5,  0.0, 0.0);

  constr.add_constr(Elem_GetPos(ElemIndex("ss"), 1),
  		    1e4, 1e4, 1e-2, 1e-2, 0e0, 0e0,
  		    0.0, 0.0, 4.0,  2.5,  0.0, 0.0);

  constr.add_constr(Elem_GetPos(ElemIndex("ls"), 1),
  		    1e4, 1e4, 1e-10, 1e-10, 0e0, 0e0,
  		    0.0, 0.0, 15.0,  4.0,   0.0, 0.0);

 lat_prms.bn_tol = 1e-6; lat_prms.step = 1.0;

  lat_constr.Fnum_b3.push_back(ElemIndex("sd1a"));
  lat_constr.Fnum_b3.push_back(ElemIndex("sfh"));
  lat_constr.Fnum_b3.push_back(ElemIndex("sd1b"));

  lat_constr.Fnum_b1.push_back(ElemIndex("b1"));
  lat_constr.Fnum_b1.push_back(ElemIndex("b2"));

  lat_constr.high_ord_achr_scl = 1e2;
  for (k = 0; k < 2; k++)
    lat_constr.high_ord_achr_nu[k] = high_ord_achr_nu[k]/2.0;

  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ls"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("b2"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ms"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("b2"), 3));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ss"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("b2"), 5));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ms"), 2));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("b2"), 7));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ls"), 2));

  n = lat_constr.high_ord_achr_Fnum.size() - 1;
  lat_constr.high_ord_achr_dnu.resize(n);
  for (k = 0; k < n; k++)
    lat_constr.high_ord_achr_dnu[k].resize(2, 0e0);

  lat_constr.mI_scl[X_] = 1e-10; lat_constr.mI_scl[Y_] = 1e-10;
  for (k = 0; k < 2; k++)
    lat_constr.mI0[k] = mI_nu_ref[k];

  lat_constr.eps_x_scl = 1e2; lat_constr.eps0_x = 0.190;

  // 2 TBA: phi = 15.
  lat_constr.phi_scl = 1e0; lat_constr.phi0 = 30.0;

  lat_constr.L_scl = 1e-10; lat_constr.L0 = 10.0;

  lat_constr.ini_constr(true);
}


void fit_ksi1(const double ksi_x, const double ksi_y)
{
  double           svd[2];
  std::vector<int> Fnum;

  Fnum.push_back(ElemIndex("sf1"));
  Fnum.push_back(ElemIndex("sd1"));
  Fnum.push_back(ElemIndex("sd2"));

  fit_ksi1(Fnum, 0e0, 0e0, 1e1, svd);

  prt_b3(Fnum);
}


int main(int argc, char *argv[])
{
  char   buffer[BUFSIZ];
  double eps_x, dnu[2];

  reverse_elem = !false;

  trace = false;

  // Unbuffered output.
  setvbuf(stdout, buffer, _IONBF, BUFSIZ);

  if (!true)
    Read_Lattice(argv[1]);
  else {
    rdmfile(argv[1]);
    globval.dPcommon = 1e-10;
  }

  globval.H_exact        = false;  globval.quad_fringe = false;
  globval.Cavity_on      = false;  globval.radiation   = false;
  globval.emittance      = false;  globval.IBS         = false;
  globval.pathlength     = false;  globval.bpm         = 0;
  globval.dip_edge_fudge = true;

  if (ps_rot) {
    Ring_GetTwiss(true, 0e0); printglob();
    dnu[X_] = 0.0; dnu[Y_] = 0.0;
    set_map(ElemIndex("ps_rot"), dnu);
    dnu[X_] = 0.1; dnu[Y_] = 0.2;
    set_map(ElemIndex("ps_rot"), dnu);
    Ring_GetTwiss(true, 0e0); printglob();
  }

  if (false) {
    if (true) {
      Ring_GetTwiss(true, 0e0); printglob();
    } else
      ttwiss(lat_constr.ic[0], lat_constr.ic[1],
	     lat_constr.ic[2], lat_constr.ic[3], 0e0);

    eps_x = get_eps_x1(true);
    printf("\neps_x = %5.3f nm.rad\n\n", eps_x);

    prt_lat("linlat1.out", globval.bpm, true);
    prt_lat("linlat.out", globval.bpm, true, 10);
    prt_chrom_lat();
  }

  if (false) fit_ksi1(0e0, 0e0);

  switch (opt_case) {
  case 1:
      // Optimize Standard Straight: mI & 2*nu.
      opt_mI_std(lat_prms, lat_constr);
    no_sxt();
    fit_powell(lat_prms, 1e-3, f_achrom);
    break;
  case 2:
    // Match Long Straight: mI.
    match_ls(lat_prms, lat_constr);
    fit_powell(lat_prms, 1e-3, f_match);
    break;
  case 3:
    opt_mI_sp(lat_prms, lat_constr);
    no_sxt();
    fit_powell(lat_prms, 1e-3, f_achrom);
    break;
  default:
    printf("\nmain: unknown opt_case %d\n", opt_case);
    exit(1);
    break;
  }

  if (false) {
    // Match ALS-U.
    match_als_u(lat_prms, lat_constr);
    fit_powell(lat_prms, 1e-3, f_match);
  }

  if (false) {
    // Optimize TBA & Mid Straight: Higher-Order-Achromat.
    opt_tba(lat_prms, lat_constr);
    no_sxt();
    fit_powell(lat_prms, 1e-3, f_achrom);
  }

  if (false) {
    // Match Mid Straight.
    match_ms(lat_prms, lat_constr);
    no_sxt();
    fit_powell(lat_prms, 1e-3, f_match);
  }

  if (false) {
    // Optimize Standard Cell: Higher-Order-Achromat.
    opt_std_cell(lat_prms, lat_constr);
    no_sxt();
    if (true)
      fit_powell(lat_prms, 1e-3, f_achrom);
    else
      fit_conj_grad(lat_prms, f_achrom);
  }

  if (false) {
    // Optimize Long Cell: Higher-Order-Achromat.
    opt_long_cell(lat_prms, lat_constr);
    no_sxt();
    if (true)
      fit_powell(lat_prms, 1e-3, f_achrom);
    else
      fit_conj_grad(lat_prms, f_achrom);
  }

  if (false) {
    // Optimize Short Cell: Higher-Order-Achromat.
    opt_short_cell(lat_prms, lat_constr);
    no_sxt();
    if (true)
      fit_powell(lat_prms, 1e-3, f_achrom);
    else
      fit_conj_grad(lat_prms, f_achrom);
  }
}
