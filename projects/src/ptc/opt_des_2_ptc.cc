#define NO 2

#include "tracy_lib.h"

//#include "Powell/src/newuoa.h"

int no_tps   = NO,
    ndpt_tps = 5;


const bool
  ps_rot          = !false, // Note, needs to be zeroed; after use.
  qf6_rb          = false,
  sp_short        = false,
  sp_std          = true,
  pert_dip_cell   = false,
  dphi            = !false,
  long_grad_dip[] = {!false, !false};

/* opt_case:
     opt_mI_std  1,
     match_ls    2,
     opt_mi_sp   3,
     match_ss    4.                                                           */
const int opt_case = 3;

// From Center of Mid Straight: alpha, beta, eta, eta'.
const int    n_ic = 4;
const double
  ic[n_ic][2] =
    {{-0.0033444062, 0.0175010857}, {3.3128077224, 1.0532911743},
     {0.0179820391, 0.0000000000}, {0.0, 0.0}},
  beta_ms[] = { 2.0, 1.5},
  beta_ss[] = { 2.0, 1.5},
  beta_ls[] = {10.0, 4.0},
  eta_ms_x  = 15e-3;

#define LAT_CASE 1

const double
#if LAT_CASE == 1
  eps0_x             = 0.097,
  dnu[]              = {0.0/6.0, 0.0/6.0},
  high_ord_achr_nu[] = {11.0/4.0, 7.0/8.0},
#elif LAT_CASE == 2
  eps0_x             = 0.147,
  high_ord_achr_nu[] = {9.0/4.0, 7.0/8.0},
#endif
  mI_dnu[]           = {0.0, 0.0},
  mI_nu_ref[]        = {1.5, 0.5};

double rad2deg(const double a) { return a*180e0/M_PI; }

double deg2rad(const double a) { return a*M_PI/180e0; }


double eps_x, sigma_E;


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
    ksi1_scl,
    ksi1_ctrl_scl[3],
    ksi1_svd_scl,
    ksi1_svd[2],
    phi_scl,
    phi_tot, phi0,       // Cell bend angle.
    high_ord_achr_scl[2],
    high_ord_achr_nu[2], // Higher-Order-Achromat Phase Advance.
    mI_scl[2],
    mI0[2],              // -I Transformer.
    alpha_c_scl,         // alpha_c.
    L_scl,
    ksi2[2],
    ksi2_scl,
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
    ksi1_scl = 0e0;
    ksi1_ctrl_scl[0] = ksi1_ctrl_scl[1] = ksi1_ctrl_scl[2] = 0e0;
    ksi1_svd_scl = 0e0;
    drv_terms_simple_scl = 0e0;
    high_ord_achr_scl[X_] = high_ord_achr_scl[Y_] = 0e0;
    mI_scl[X_] = mI_scl[Y_] = 0e0;
    L_scl = 0e0;
    ksi2_scl = 0e0;
  }

  void add_constr(const int loc,
		  const double scl1, const double scl2, const double scl3,
		  const double scl4, const double scl5, const double scl6,
		  const double v1, const double v2, const double v3,
		  const double v4, const double v5, const double v6);
  void ini_constr(const bool ring);
  void get_Jacobian(param_type &lat_prms);
  double get_chi2(const param_type &prms, double *bn, const bool prt);
  void get_dchi2(double *bn, double *df) const;
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


void set_lin_map(const int Fnum)
{
  int loc, k;

  loc = Elem_GetPos(Fnum, 1);
  if (globval.mat_meth && (Cell[loc].Elem.Pkind == Mpole)
      && (Cell[loc].Elem.M->Pthick == thick))
    for (k = 1; k <= GetnKid(Fnum); k++) {
      loc = Elem_GetPos(Fnum, k);
      Cell[loc].Elem.M->M_lin = get_lin_map(Cell[loc].Elem, 0e0);
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
  set_lin_map(Fnum);
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

  for (k = 0; k < (int)Fnum.size(); k++) {
    set_bn_design_fam(Fnum[k], Quad, b2, 0e0);
    set_lin_map(Fnum[k]);
  }
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
      if (Fnum[i-1] > 0) {
	set_bn_design_fam(Fnum[i-1], n[i-1], bn_ext, 0e0);
	if (n[i-1] == Quad) set_lin_map(Fnum[i-1]);
      } else
	set_grad_dip_b2(grad_dip_Fnum[i-1], bn_ext);
    else if (n[i-1] == -1) {
      // Position.
      set_ds(Fnum[i-1], bn_ext-l0[i-1]);
      l0[i-1] = bn_ext;
    } else if (n[i-1] == -2) {
      // Length.
      set_L(Fnum[i-1], bn_ext); get_S();
      set_lin_map(Fnum[i-1]);
    } else if (n[i-1] == -3) {
      // Bend angle; L is fixed.
      if (Fnum[i-1] > 0)
	set_phi(Fnum[i-1], bn_ext);
      else
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
      printf("    %2d ", i);
      prt_name(stdout, Cell[loc].Elem.PName, ":", 16);
    } else {
      loc = Elem_GetPos(grad_dip_Fnum[i-1][0], 1);
      printf("    %2d ", i);
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
      if (prt) printf("\nset_ds Knum = %d (odd):", Knum);
      set_dL(Cell[loc+1].Fnum, Knum, -ds);
      set_dL(Cell[loc-1].Fnum, Knum, ds);
    } else {
      if (prt) printf("\nset_ds Knum = %d (even):", Knum);
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


double get_eps_x1(const bool track, double &eps_x, double &sigma_E)
{
  // eps_x [nm.rad].
  long int     lastpos;
  int          k;
  double       I[6];
  ss_vect<tps> A;

  const bool   prt     = false;
  const double C_q_scl = 1e18*C_q/sqr(m_e);

  if (track) {
    globval.emittance = true;
    // A = get_A(ic[0], ic[1], ic[2], ic[3]);
    A = putlinmat(6, globval.Ascr);
    Cell_Pass(0, globval.Cell_nLoc, A, lastpos);
    // Cell_Twiss(0, globval.Cell_nLoc, A, false, false, 0e0);
    globval.emittance = false;
  }

  get_I(I, false);

  eps_x = 1e9*C_q_scl*sqr(globval.Energy)*I[5]/(I[2]-I[4]);
  sigma_E = sqrt(C_q_scl*sqr(globval.Energy)*I[3]/(2e0*I[2]+I[4]));

  if (prt) {
    printf("\nget_eps_x1:\n");
    printf("  I[2..5]:");
    for (k = 2; k <= 5; k++)
      printf(" %10.3e", I[k]);
    printf("\n");
    printf("  eps_x   = %9.3f nm.rad\n", eps_x);
    printf("  sigma_E = %9.3e nm.rad\n", sigma_E);
    printf("  J_x     = %9.3f, J_z = %5.3f\n", 1.0-I[4]/I[2], 2.0+I[4]/I[2]);
  }

  return eps_x;
}


double get_lin_opt(constr_type &constr)
{
  double       eps_x;
  ss_vect<tps> A;

  if (lat_constr.ring) {
    Ring_GetTwiss(true, 0e0);
    eps_x = get_eps_x1(true, eps_x, sigma_E);
  } else {
    globval.emittance = true;
    A = get_A(constr.ic[0], constr.ic[1], constr.ic[2], constr.ic[3]);
    Cell_Twiss(0, globval.Cell_nLoc, A, false, false, 0e0);
    eps_x = get_eps_x1(false, eps_x, sigma_E);
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


void get_ksi1_ctrl(constr_type &constr)
{
  int k, loc;

  constr.ksi1_ctrl.clear();
  for (k = 0; k < 3; k++) {
    loc = Elem_GetPos(constr.Fnum_b3[k], 1);
    constr.ksi1_ctrl.push_back((Cell[loc].Beta[Y_]-Cell[loc].Beta[X_])
			       *Cell[loc].Eta[X_]);
  }
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


double arccos_(double x)
{
  double y;

  if (x == 0e0)
    y = 1e0;
  else if (x >= 0e0)
    y = atan(sqrt(1e0/(x*x)-1e0));
  else
    y = (pi-atan(sqrt(1e0/(sqr(x))-1e0)));

  return y;
}


void get_nu(Matrix &M, double nu[])
{
  int    k;
  double tr2;

  const bool prt = false;

  for (k = 0; k < 2; k++) {
    tr2 = globval.OneTurnMat[2*k][2*k] + globval.OneTurnMat[2*k+1][2*k+1];
    nu[k] = arccos_(tr2/2e0)/(2e0*M_PI);
    if (globval.OneTurnMat[2*k][2*k+1] < 0e0) nu[k] = 1e0 - nu[k];
  }

  if (prt) printf("\nget_nu: [%10.8f, %10.8f]\n", nu[X_], nu[Y_]);
}


void get_ksi2_(const double delta, double ksi2[])
{
  long int lastpos;
  int      k;
  double   nu[3][2];

  const bool prt = false;

  getcod(-delta, lastpos);
  get_nu(globval.OneTurnMat, nu[0]);
  getcod(0e0, lastpos);
  get_nu(globval.OneTurnMat, nu[1]);
  getcod(delta, lastpos);
  get_nu(globval.OneTurnMat, nu[2]);
  for (k = 0; k < 2; k++)
    ksi2[k] = (nu[2][k]-2e0*nu[1][k]+nu[0][k])/(2e0*sqr(delta));

  if (prt) printf("\nget_ksi2_: [%10.8f, %10.8f]\n", ksi2[X_], ksi2[Y_]);
}


double constr_type::get_chi2(const param_type &prms, double *bn,
			     const bool prt)
{
  int    j, k;
  double chi2, dchi2[3], mean, geom_mean, bn_ext;

  const bool   extra = false;
  const double delta = 1e-6, scl_extra = 1e8;

  if (prt) printf("\nget_chi2:\n");

  chi2 = 0e0;

  if (extra) {
    if (prt) printf("\n  extra:\n");
    if (!true) {
      // Reduce QF1 #11.
      k = 11;
      bn_ext = bn_bounded(bn[k], prms.bn_min[k-1], prms.bn_max[k-1]);
    } else {
      // Tweak D_09.
      bn_ext = get_L(ElemIndex("d_10_1"), 1); 
      dchi2[0] = scl_extra*sqr(bn_ext+0.14);
      chi2 += dchi2[0];
      if (prt)
	printf("  D_10_1:            %10.3e (%10.3e)\n", dchi2[0], bn_ext);
      bn_ext = get_L(ElemIndex("d_10_2"), 1); 
      dchi2[0] = scl_extra*sqr(bn_ext+0.14);
      chi2 += dchi2[0];
      if (prt)
	printf("  D_10_2:            %10.3e (%10.3e)\n", dchi2[0], bn_ext);
    }
  }

  if (eps_x_scl != 0e0) {
    dchi2[X_] = eps_x_scl*sqr(eps_x-eps0_x);
    chi2 += dchi2[X_];
    if (prt) printf("  eps_x:             %10.3e\n", dchi2[X_]);
  }

  if (prt) printf("\n"); 
  for (j = 0; j < n_loc; j++) {
    if (prt) {
      printf("  ");
      prt_name(stdout, Cell[loc[j]].Elem.PName, ":", 8);
      printf("\n");
    }
    if ((value_scl[j][0] != 0e0) || (value_scl[j][1] != 0e0)) {
      for (k = 0; k < 2; k++) {
	dchi2[k] = value_scl[j][k]*sqr(Cell[loc[j]].Alpha[k]-value[j][k]);
	chi2 += dchi2[k];
      }
      if (prt) printf("    alpha:          [%10.3e, %10.3e]\n",
		      dchi2[X_], dchi2[Y_]);
    }
    if ((value_scl[j][2] != 0e0) || (value_scl[j][3] != 0e0)) {
      for (k = 0; k < 2; k++) {
	dchi2[k] = value_scl[j][k+2]*sqr(Cell[loc[j]].Beta[k]-value[j][k+2]);
	chi2 += dchi2[k];
      }
      if (prt) printf("    beta:           [%10.3e, %10.3e]\n",
		      dchi2[X_], dchi2[Y_]);
    }
    if ((value_scl[j][4] != 0e0) || (value_scl[j][5] != 0e0)) {
	dchi2[0] = value_scl[j][4]*sqr(Cell[loc[j]].Eta[X_]-value[j][4]);
	dchi2[1] = value_scl[j][5]*sqr(Cell[loc[j]].Etap[X_]-value[j][5]);
	chi2 += dchi2[0]; chi2 += dchi2[1];
	if (prt) printf("    eta_x eta'_x:   [%10.3e, %10.3e]\n",
			dchi2[0], dchi2[1]);
    }
  }

  if (L_scl != 0e0) {
    dchi2[0] = L_scl*sqr(Cell[globval.Cell_nLoc].S-L0);
    chi2 += dchi2[0];
    printf("\n  S-L0:            %10.3e\n", dchi2[0]);
  }

  if (prt) printf("\n");
  if (ksi1_scl != 0e0) {
    for (k = 0; k < 2; k++) {
      dchi2[k] = ksi1_scl*sqr(ksi1[k]);
      chi2 += dchi2[k];
    }
    if (prt) printf("  ksi1:             [%10.3e, %10.3e]\n",
		    dchi2[X_], dchi2[Y_]);
  }

  if (drv_terms_simple_scl != 0e0) {
    get_drv_terms(lat_constr);
    for (k = 0; k < 2; k++) {
      dchi2[k] = drv_terms_simple_scl*drv_terms_simple[k];
      chi2 += dchi2[k];
    }
    if (prt) printf("  drv_terms_simple: [%10.3e, %10.3e]\n",
		    dchi2[X_], dchi2[Y_]);
  }

  if ((ksi1_ctrl_scl[0] != 0e0) || (ksi1_ctrl_scl[1] != 0e0)
      || (ksi1_ctrl_scl[2] != 0e0)) {
    get_ksi1_ctrl(lat_constr);
    for (k = 0; k < 3; k++) {
      dchi2[k] = ksi1_ctrl_scl[k]/sqr(ksi1_ctrl[k]);
      chi2 += dchi2[k];
    }
    if (prt) printf("  ksi1_ctrl:        [%10.3e, %10.3e, %10.3e]\n",
		    dchi2[0], dchi2[1], dchi2[2]);
  }
  
  if (ksi1_svd_scl != 0) {
    mean = (ksi1_svd[X_]+ksi1_svd[Y_])/2e0;
    geom_mean = sqrt(ksi1_svd[X_]*ksi1_svd[Y_]);
    dchi2[0] = ksi1_svd_scl/sqr(geom_mean/mean);
    chi2 += dchi2[0];
    if (prt) printf("  ksi1_svd:         %10.3e\n", dchi2[0]);
  }

  if (ksi2_scl != 0e0) {
    get_ksi2_(delta, ksi2);
    for (k = 0; k < 2; k++) {
      dchi2[k] = ksi2_scl*sqr(ksi2[k]);
      chi2 += dchi2[k];
    }
    if (prt) printf("  ksi2:             [%10.3e, %10.3e]\n",
		    ksi2[X_], ksi2[Y_]);
  }

  if ((mI_scl[X_] != 0e0) || (mI_scl[Y_] != 0e0)) {
    get_mI(lat_constr);
    for (j = 0; j < (int)mI.size(); j++) {
      for (k = 0; k < 2; k++) {
	dchi2[k] = mI_scl[k]*sqr(mI[j][k]-mI0[k]);
	chi2 += dchi2[k];
      }
      if (prt) printf("  mI:               [%10.3e, %10.3e]\n",
		      dchi2[X_], dchi2[Y_]);
    }
  }

  if ((high_ord_achr_scl[X_] != 0e0) || (high_ord_achr_scl[Y_] != 0e0)) {
    get_high_ord_achr(lat_constr);
    for (j = 0; j < (int)high_ord_achr_dnu.size(); j++) {
      for (k = 0; k < 2; k++) {
	dchi2[k] =
	  high_ord_achr_scl[k]*sqr(high_ord_achr_dnu[j][k]-high_ord_achr_nu[k]);
	chi2 += dchi2[k];
      }
      if (prt) printf("  high_ord_achr:    [%10.3e, %10.3e]\n",
		      dchi2[X_], dchi2[Y_]);
    }
  }

  return chi2;
}


void constr_type::get_dchi2(double *bn, double *df) const
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
    df[k+1] = lat_constr.get_chi2(lat_prms, bn, false);

    constr_dparam(lat_prms.Fnum[k], lat_prms.n[k], -2e0*eps);
    eps_x = get_lin_opt(lat_constr);
    df[k+1] -= lat_constr.get_chi2(lat_prms, bn, false);
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
  printf("    ksi2         = [%10.3e, %10.3e]\n\n",
	 constr.ksi2[X_], constr.ksi2[Y_]);
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

  if ((high_ord_achr_scl[X_] != 0e0) || (high_ord_achr_scl[Y_] != 0e0))
    prt_high_ord_achr(*this);

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

  const bool prt = false;

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


void fit_sim_anneal(param_type &lat_prms, const double eps,
		    double (*f)(double *))
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

  lat_constr.get_dchi2(b2, df);
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
  if (chrom) prt_chrom_lat("chromlat.out");

  prt_b2(lat_prms);
  prt_b3(lat_constr.Fnum_b3);
}


double f_match(double *b2)
{
  double      chi2;
  static bool first = true;

  const int n_prt = 5;

  lat_prms.set_prm(b2);
  if (lat_constr.phi_scl != 0e0) phi_corr(lat_constr);

  eps_x = get_lin_opt(lat_constr);

  if ((lat_constr.high_ord_achr_scl[X_] != 0e0)
      || (lat_constr.high_ord_achr_scl[Y_] != 0e0))
    get_high_ord_achr(lat_constr);

  chi2 = lat_constr.get_chi2(lat_prms, b2, false);

  if (first || (chi2 < lat_constr.chi2)) {
    first = false;
    if (lat_constr.n_iter % n_prt == 0) {
      // Print dchi2.
      lat_constr.get_chi2(lat_prms, b2, true);
      prt_f(b2, chi2, lat_constr, lat_prms, false);
    }
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

    chi2 = lat_constr.get_chi2(lat_prms, b2, false);

    if (first || (chi2 < lat_constr.chi2)) {
      first = false;
      if (lat_constr.n_iter % n_prt == 0) {
	// Print dchi2.
	lat_constr.get_chi2(lat_prms, b2, true);
	prt_f(b2, chi2, lat_constr, lat_prms, true);
      }

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
	 "  mI_scl               = %9.3e %9.3e\n"
	 "  high_ord_achr_scl    = [%9.3e, %9.3e]\n"
	 "  phi_scl              = %9.3e\n",
	 constr.eps_x_scl, constr.ksi1_svd_scl,
	 constr.drv_terms_simple_scl,
	 mI_nu_ref[X_], mI_nu_ref[Y_], constr.mI_scl[X_], constr.mI_scl[Y_],
	 constr.high_ord_achr_scl[X_], constr.high_ord_achr_scl[Y_],
	 constr.phi_scl);
}


void set_dip_cell_std(param_type &prms, constr_type &constr)
{
  std::vector<int>    grad_dip_Fnum, mI_loc;
  std::vector<double> grad_dip_scl;

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
    prms.add_prm(grad_dip_Fnum, grad_dip_scl, -3, -15.0, 15.0, 1.0);
  if (long_grad_dip[0])
    prms.add_prm(grad_dip_Fnum, grad_dip_scl,  2, -15.0, 15.0, 1.0);

  lat_constr.Fnum_b1.push_back(-ElemIndex("dl1a_1"));
  lat_constr.grad_dip_Fnum_b1.push_back(grad_dip_Fnum);

  grad_dip_Fnum.clear();
  grad_dip_Fnum.push_back(ElemIndex("dl2a_1"));
  grad_dip_Fnum.push_back(ElemIndex("dl2a_2"));
  grad_dip_Fnum.push_back(ElemIndex("dl2a_3"));
  grad_dip_Fnum.push_back(ElemIndex("dl2a_4"));
  grad_dip_Fnum.push_back(ElemIndex("dl2a_5"));
  if (dphi)
    prms.add_prm(grad_dip_Fnum, grad_dip_scl, -3, -15.0, 15.0, 1.0);
  if (long_grad_dip[1])
    prms.add_prm(grad_dip_Fnum, grad_dip_scl,  2, -15.0, 15.0, 1.0);

  lat_constr.Fnum_b1.push_back(-ElemIndex("dl2a_1"));
  lat_constr.grad_dip_Fnum_b1.push_back(grad_dip_Fnum);

  if (dphi) {
    // Dipole Cell.
    prms.add_prm("qf4", -3, -15.0, 15.0, 1.0);
    lat_constr.Fnum_b1.push_back(ElemIndex("qf4"));
    if (qf6_rb) {
      prms.add_prm("qf6", -3, -15.0, 15.0, 1.0);
      lat_constr.Fnum_b1.push_back(ElemIndex("qf6"));
    }
    prms.add_prm("qf8", -3, -15.0, 15.0, 1.0);
    lat_constr.Fnum_b1.push_back(ElemIndex("qf8"));

    // Commented out must be defined last.
    // prms.add_prm("dq1", -3, -15.0,   15.0,  1.0);
    lat_constr.Fnum_b1.push_back(ElemIndex("dq1"));
  }

  prms.add_prm("dq1", 2, -15.0, 15.0, 1.0);

  prms.add_prm("qd3", 2, -15.0,  0.0, 1.0);
  prms.add_prm("qf4", 2,   0.0, 15.0, 1.0);
  prms.add_prm("qd5", 2, -15.0,  0.0, 1.0);

  // prms.add_prm("qd3", -1,  -0.4, -0.4, 1.0);
  // prms.add_prm("qd5", -1,  -0.4, -0.4, 1.0);

  lat_constr.phi_scl = (dphi)? 1e0 : 0e0;
}


void set_b2_ms_std(param_type &prms)
{
  // Mid-Straight.
  if (false)
    prms.add_prm("qp_q",  2, -15.0, 15.0, 1.0);
  prms.add_prm("qf1",  2,   0.0, 15.0, 1.0);
  prms.add_prm("qd2",  2, -15.0,  0.0, 1.0);

  // Standard-Straight.
  // prms.add_prm("qf6", -1, -0.3,  0.01, 1.0);
  prms.add_prm("qf6",  2,  0.0, 15.0, 1.0);
  prms.add_prm("qf8",  2,  0.0, 15.0, 1.0);
}


void set_constr_std(constr_type &constr)
{
  // constr.add_constr(Elem_GetPos(ElemIndex("qf1"), 1)-1,
  // 		    0e0, 0e0, 0e0, 0e0, 1e8, 1e7,
  // 		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  // constr.add_constr(Elem_GetPos(ElemIndex("qf1"), 2),
  // 		    0e0, 0e0, 0e0, 0e0, 1e8, 1e7,
  // 		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  constr.add_constr(Elem_GetPos(ElemIndex("dl1a_5"), 1)-1,
		    0e0, 0e0, 0e0, 0e0, 1e8, 1e7,
		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  constr.add_constr(Elem_GetPos(ElemIndex("dl1a_5"), 2),
		    0e0, 0e0, 0e0, 0e0, 1e8, 1e7,
		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  constr.add_constr(Elem_GetPos(ElemIndex("ms"), 1),
		    1e7, 1e7, 1e2,         1e2,         1e7,   0e0,
		    0.0, 0.0, beta_ms[X_], beta_ms[Y_], 0.018, 0.0);
  constr.add_constr(Elem_GetPos(ElemIndex("ss"), 1),
		    1e5, 1e5, 1e2,         1e2,         1e7, 1e7,
		    0.0, 0.0, beta_ss[X_], beta_ss[Y_], 0.0, 0.0);
  // Symmetry.
  constr.add_constr(Elem_GetPos(ElemIndex("sf1"), 1),
		    5e2, 1e3, 0e0, 0e0, 0e0, 0e0,
		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

  if (false)
    // Increase beta_x.
    constr.add_constr(Elem_GetPos(ElemIndex("sf1"), 1),
		      0e5, 0e5, 1e1,  1e2, 0e7, 0e7,
		      0.0, 0.0, 10.0, 1.0, 0.0, 0.0);
  if (false)
    // Increase beta_y.
    constr.add_constr(Elem_GetPos(ElemIndex("sd1"), 1),
		      0e5, 0e5, 5e2, 5e3,  0e7, 0e7,
		      0.0, 0.0, 5.0, 15.0, 0.0, 0.0);
  if (false)
    // Increase beta_y.
    constr.add_constr(Elem_GetPos(ElemIndex("sd2"), 1),
		      0e5, 0e5, 5e2, 5e3,  0e7, 0e7,
		      0.0, 0.0, 5.0, 10.0, 0.0, 0.0);
}


void set_b3_constr_std(constr_type &constr)
{
  int              k, n;
  std::vector<int> mI_loc;

  lat_constr.Fnum_b3.push_back(ElemIndex("sf1"));
  lat_constr.Fnum_b3.push_back(ElemIndex("sd1"));
  lat_constr.Fnum_b3.push_back(ElemIndex("sd2"));
  // lat_constr.Fnum_b3.push_back(ElemIndex("sh2"));

  mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 1));
  mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 3));
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
}


void set_weights_std(constr_type &constr)
{
  lat_constr.eps_x_scl             = 1e6;

  lat_constr.ksi1_scl              = 1e1;
  lat_constr.drv_terms_simple_scl  = 1e-2;
  lat_constr.ksi1_ctrl_scl[0]      = 1e-1;
  lat_constr.ksi1_ctrl_scl[1]      = 1e0*1e0;
  lat_constr.ksi1_ctrl_scl[2]      = 1e0*1e0;
  // Not useful.
  lat_constr.ksi1_svd_scl          = 0e3;
  lat_constr.mI_scl[X_]            = 5e6;
  lat_constr.mI_scl[Y_]            = 5e6;
  lat_constr.high_ord_achr_scl[X_] = 5e6;
  lat_constr.high_ord_achr_scl[Y_] = 0e6;

  lat_constr.alpha_c_scl           = 1e-6;

  // Super Period.
  lat_constr.phi0 = 15.0;
  lat_constr.L_scl = 0e-10; lat_constr.L0 = 10.0;
}


void opt_mI_std(param_type &prms, constr_type &constr)
{
  // Parameter Type:
  //   Bend Angle  -3,
  //   Length      -2,
  //   Position    -1,
  //   Quadrupole   2.

  // Set parameters; initialized by optimizer.
  lat_prms.bn_tol = 1e-5; lat_prms.step = 1e0;

  set_dip_cell_std(prms, constr);
  set_b2_ms_std(prms);

  lat_constr.eps0_x = eps0_x;
  set_constr_std(constr);
  set_b3_constr_std(constr);
  set_weights_std(constr);

  lat_constr.ini_constr(true);

  prt_prms(lat_constr);
}


void set_dip_cell_sp(param_type &prms, constr_type &constr)
{
  std::vector<int>    grad_dip_Fnum;
  std::vector<double> grad_dip_scl;

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
    prms.add_prm(grad_dip_Fnum, grad_dip_scl, -3, -15.0, 15.0, 1.0);
  if (long_grad_dip[0])
    prms.add_prm(grad_dip_Fnum, grad_dip_scl,  2, -15.0, 15.0, 1.0);

  lat_constr.Fnum_b1.push_back(-ElemIndex("dl1a_1"));
  lat_constr.grad_dip_Fnum_b1.push_back(grad_dip_Fnum);

  if (true) {
    // Different longitudinal gradient dipoles.
    grad_dip_Fnum.clear();
    grad_dip_Fnum.push_back(ElemIndex("dl2a_1"));
    grad_dip_Fnum.push_back(ElemIndex("dl2a_2"));
    grad_dip_Fnum.push_back(ElemIndex("dl2a_3"));
    grad_dip_Fnum.push_back(ElemIndex("dl2a_4"));
    grad_dip_Fnum.push_back(ElemIndex("dl2a_5"));
    if (dphi)
      prms.add_prm(grad_dip_Fnum, grad_dip_scl, -3, -15.0, 15.0, 1.0);
    if (long_grad_dip[1])
      prms.add_prm(grad_dip_Fnum, grad_dip_scl,  2, -15.0, 15.0, 1.0);

    lat_constr.Fnum_b1.push_back(-ElemIndex("dl2a_1"));
    lat_constr.grad_dip_Fnum_b1.push_back(grad_dip_Fnum);
  }

  if (dphi) {
    // Dipole Cell.
    prms.add_prm("qf4", -3, -15.0, 15.0, 1.0);
    lat_constr.Fnum_b1.push_back(ElemIndex("qf4"));
    if (qf6_rb) {
      prms.add_prm("qf6", -3, -15.0, 15.0, 1.0);
      lat_constr.Fnum_b1.push_back(ElemIndex("qf6"));
    }
    prms.add_prm("qf8", -3, -15.0, 15.0, 1.0);
    lat_constr.Fnum_b1.push_back(ElemIndex("qf8"));

    // Commented out must be defined last.
    // prms.add_prm("dq1", -3, -15.0,   15.0,  1.0);
    lat_constr.Fnum_b1.push_back(ElemIndex("dq1"));
  }

  prms.add_prm("dq1", 2, -15.0, 15.0, 1.0);

  prms.add_prm("qd3", 2, -15.0,  0.0, 1.0);
  prms.add_prm("qf4", 2,   0.0, 15.0, 1.0);
  prms.add_prm("qd5", 2, -15.0,  0.0, 1.0);

  lat_constr.phi_scl = (dphi)? 1e0 : 0e0;
}


void set_b2_ms_std_ls_sp(param_type &prms)
{
  // Mid Straight.
  if (false)
    prms.add_prm("qp_q", 2,  -1.0,  1.0, 1.0);
  prms.add_prm("qd2",  2, -15.0, 15.0, 1.0);
  prms.add_prm("qf1",  2, -15.0, 15.0, 1.0);

  // Standard Straight.
  prms.add_prm("qf6", 2, -15.0, 15.0, 1.0);
  prms.add_prm("qf8", 2, -15.0, 15.0, 1.0);

  // Long Straight.
  if (pert_dip_cell)
    prms.add_prm("qd3_c1",   2, -15.0, 15.0, 1.0);
  prms.add_prm("qd2_c1",   2, -15.0,  0.0, 1.0);
  prms.add_prm("qf1_c1",   2, -15.0,  0.0, 1.0);
  prms.add_prm("quad_add", 2,   0.0, 15.0, 1.0);
}


void set_constr_sp(constr_type &constr)
{
  int k;

  // Lattice constraints are: alpha_x,y, beta_x,y, eta_x, eta'_x.
  // Include Matching Sections for Std & Long Straights.
  constr.add_constr(Elem_GetPos(ElemIndex("quad_add"), 1)-1,
		    0e0, 0e0, 0e0, 0e0, 1e8, 1e8,
		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  constr.add_constr(Elem_GetPos(ElemIndex("qf1"), 1),
		    0e0, 0e0, 0e0, 0e0, 1e8, 1e8,
		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  // Include constraint on alpha; in case of using ps_rot.
  constr.add_constr(Elem_GetPos(ElemIndex("ms"), 1),
		    1e6, 1e6, 1e2,         1e3,         1e8,      1e8,
		    0.0, 0.0, beta_ms[X_], beta_ms[Y_], eta_ms_x, 0.0);
  for (k = 1; k <= 2; k++)
    constr.add_constr(Elem_GetPos(ElemIndex("ss"), k),
		      1e6, 1e6, 1e2,         1e3,         1e8, 1e8,
		      0.0, 0.0, beta_ss[X_], beta_ss[Y_], 0.0, 0.0);
  constr.add_constr(Elem_GetPos(ElemIndex("ls"), 1),
		    1e6, 1e6, 1e2,         1e0,         1e8, 1e8,
		    0.0, 0.0, beta_ls[X_], beta_ls[Y_], 0.0, 0.0);
  // Symmetry.
  constr.add_constr(Elem_GetPos(ElemIndex("sf1"), 1),
		    5e2, 1e3, 0e0, 0e0, 0e0, 0e0,
		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

  if (false)
    // Increase beta_x.
    constr.add_constr(Elem_GetPos(ElemIndex("sf1"), 1),
		      0e5, 0e5, 1e3,  5e2, 0e7, 0e7,
		      0.0, 0.0, 10.0, 1.0, 0.0, 0.0);
  if (false)
    // Increase beta_y.
    constr.add_constr(Elem_GetPos(ElemIndex("sd2"), 1),
		      0e5, 0e5, 1e1, 1e3, 0e7, 0e7,
		      0.0, 0.0, 1.0, 8.0, 0.0, 0.0);
  if (false)
    // Increase beta_y.
    constr.add_constr(Elem_GetPos(ElemIndex("qd5"), 1),
		      0e5, 0e5, 0e1, 1e3,  0e7, 0e7,
		      0.0, 0.0, 0.0, 10.0, 0.0, 0.0);
}


void set_b3_constr_sp(constr_type &constr)
{
  int              k, n;
  std::vector<int> mI_loc;

  lat_constr.Fnum_b3.push_back(ElemIndex("sf1"));
  lat_constr.Fnum_b3.push_back(ElemIndex("sd1"));
  lat_constr.Fnum_b3.push_back(ElemIndex("sd2"));
  // lat_constr.Fnum_b3.push_back(ElemIndex("sh2"));

  lat_constr.eps0_x = eps0_x;

  mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 1));
  mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 3));
  lat_constr.mI_loc.push_back(mI_loc);
  mI_loc.clear();
  mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 5));
  mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 7));
  lat_constr.mI_loc.push_back(mI_loc);
  mI_loc.clear();
  mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 9));
  mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 11));
  lat_constr.mI_loc.push_back(mI_loc);
  mI_loc.clear();
  mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 13));
  mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 15));
  lat_constr.mI_loc.push_back(mI_loc);

  for (k = 0; k < 2; k++)
    lat_constr.high_ord_achr_nu[k] = high_ord_achr_nu[k];

  // lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ls"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ss"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ss"), 2));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ss"), 3));

  n = lat_constr.high_ord_achr_Fnum.size() - 1;
  lat_constr.high_ord_achr_dnu.resize(n);
  for (k = 0; k < n; k++)
    lat_constr.high_ord_achr_dnu[k].resize(2, 0e0);

  for (k = 0; k < 2; k++)
    lat_constr.mI0[k] = mI_nu_ref[k];
}


void set_weights_sp(constr_type &constr)
{
  lat_constr.eps_x_scl             = 1e7;

  lat_constr.ksi1_scl              = 1e0;
  lat_constr.drv_terms_simple_scl  = 5e-4;
  lat_constr.ksi1_ctrl_scl[0]      = 0e-1;
  lat_constr.ksi1_ctrl_scl[1]      = 0e0;
  lat_constr.ksi1_ctrl_scl[2]      = 0e-1;
  // Not useful.
  lat_constr.ksi1_svd_scl          = 0e3;
  lat_constr.mI_scl[X_]            = 1e6;
  lat_constr.mI_scl[Y_]            = 1e6;
  lat_constr.high_ord_achr_scl[X_] = 1e-2*1e6;
  lat_constr.high_ord_achr_scl[Y_] = 0e6;
  lat_constr.ksi2_scl              = 5e0;

  lat_constr.alpha_c_scl           = 5e-7;

  if (!sp_short)
    // Super Period.
    lat_constr.phi0 = 60.0;
  else
    // Super Period with 1 Std Straight only.
    lat_constr.phi0 = 45.0;
  lat_constr.L_scl = 0e-10; lat_constr.L0 = 10.0;
}


void opt_mI_sp(param_type &prms, constr_type &constr)
{
  // Parameter Type:
  //   Bend Angle  -3,
  //   Length      -2,
  //   Position    -1,
  //   Quadrupole   2.

  // Set parameters; initialized by optimizer.
  lat_prms.bn_tol = 1e-5; lat_prms.step = 1.0;

  set_dip_cell_sp(prms, constr);
  set_b2_ms_std_ls_sp(prms);

  lat_constr.eps0_x = eps0_x;
  set_constr_sp(constr);
  set_b3_constr_sp(constr);
  set_weights_sp(constr);

  lat_constr.ini_constr(true);

  prt_prms(lat_constr);
}


void set_b2_ls(param_type &prms)
{
  if (pert_dip_cell)
    prms.add_prm("qd3_c1", 2, -15.0, 15.0, 1.0);
  // Long Straight.
  prms.add_prm("qd2_c1",   2, -15.0, 15.0, 1.0);
  prms.add_prm("qf1_c1",   2, -15.0, 15.0, 1.0);
  prms.add_prm("quad_add", 2, -15.0, 15.0, 1.0);

  // Parameters are initialized in optimizer.
}


void set_constr_ls(constr_type &constr)
{
  if (!pert_dip_cell)
    // Eta_x & eta_x' are zero after Dipole Cell.
    constr.add_constr(Elem_GetPos(ElemIndex("ls"), 1),
		      1e4, 1e4, 1e0,         1e-1,        1e-10, 1e-10,
		      0.0, 0.0, beta_ls[X_], beta_ls[Y_], 0.0,   0.0);
  else {
    constr.add_constr(Elem_GetPos(ElemIndex("ls"), 1),
		      1e2, 1e2, 1e0,         1e-1,        1e-10, 1e-10,
		      0.0, 0.0, beta_ls[X_], beta_ls[Y_], 0.0,   0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("quad_add"), 1),
		      0e0, 0e0, 0e0, 0e0, 1e5, 1e5,
		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  }
}


void set_b3_constr_ls(constr_type &constr)
{
  int k, n;

  for (k = 0; k < 2; k++)
    lat_constr.high_ord_achr_nu[k] = high_ord_achr_nu[k]/2e0;

  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ms"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ls"), 1));

  n = lat_constr.high_ord_achr_Fnum.size() - 1;
  lat_constr.high_ord_achr_dnu.resize(n);
  for (k = 0; k < n; k++)
    lat_constr.high_ord_achr_dnu[k].resize(2, 0e0);
}


void set_weights_ls(constr_type &constr)
{
  lat_constr.high_ord_achr_scl[X_] = 0e0;
  lat_constr.high_ord_achr_scl[Y_] = 0e0;
}


void match_ls(param_type &prms, constr_type &constr)
{
  // Parameter Type:
  //   Bend Angle  -3,
  //   Length      -2,
  //   Position    -1,
  //   Quadrupole   2.
  int j, k;

  // Perturbed symmetry at end of Dipole Cell: 
  //   1. Initialize with: [Qf1, Qd2, Qd3].
  //   2. Exclude for 1st pass: pert_dip_cell = false.
  //   3. Include for fine tuning.
  lat_prms.bn_tol = 1e-5; lat_prms.step = 1.0;
  set_b2_ls(prms);

  set_constr_ls(constr);
  set_b3_constr_ls(constr);
  set_weights_ls(constr);

  lat_constr.ini_constr(false);

  for (j = 0; j < n_ic; j++)
    for (k = 0; k < 2; k++)
      lat_constr.ic[j][k] = ic[j][k];
}


void set_b2_ss(param_type &prms)
{
  // Std Straight.
  if (false)
    prms.add_prm("qp_q", 2, -15.0, 15.0, 1.0);
  prms.add_prm("qd2",  2, -15.0, 15.0, 1.0);
  prms.add_prm("qf1",  2, -15.0, 15.0, 1.0);
  // prms.add_prm("qf1", -1,   0.0, -0.5, 1.0);

  // Parameters are initialized in optimizer.
}


void set_constr_ss(constr_type &constr)
{
  constr.add_constr(Elem_GetPos(ElemIndex("ss"), 1),
		    1e2, 1e2, 1e-1,        1e-1,         1e-10, 1e-10,
		    0.0, 0.0, beta_ss[X_], beta_ss[Y_],  0.0,   0.0);
  constr.add_constr(Elem_GetPos(ElemIndex("qf1"), 1),
		    0e0, 0e0, 0e0, 0e0, 1e6, 1e6,
		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
}


void set_b3_constr_ss(constr_type &constr)
{
  int k, n;

  for (k = 0; k < 2; k++)
    lat_constr.high_ord_achr_nu[k] = high_ord_achr_nu[k]/2e0;

  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ms"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ss"), 1));

  n = lat_constr.high_ord_achr_Fnum.size() - 1;
  lat_constr.high_ord_achr_dnu.resize(n);
  for (k = 0; k < n; k++)
    lat_constr.high_ord_achr_dnu[k].resize(2, 0e0);
}


void set_weights_ss(constr_type &constr)
{
  lat_constr.high_ord_achr_scl[X_] = 0e0;
  lat_constr.high_ord_achr_scl[Y_] = 0e0;
}


void match_ss(param_type &prms, constr_type &constr)
{
  // Parameter Type:
  //   Bend Angle  -3,
  //   Length      -2,
  //   Position    -1,
  //   Quadrupole   2.
  int j, k;

  // Perturbed symmetry at end of Dipole Cell: 
  //   1. Initialize with: [Qf1, Qd2, Qd3].
  //   2. Exclude for 1st pass: pert_dip_cell = false.
  //   3. Include for fine tuning.
  lat_prms.bn_tol = 1e-5; lat_prms.step = 1.0;
  set_b2_ss(prms);

  set_constr_ss(constr);
  set_b3_constr_ss(constr);
  set_weights_ss(constr);

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

  lat_constr.ksi1_svd_scl          = 0e0;
  lat_constr.high_ord_achr_scl[X_] = 0e0;
  lat_constr.high_ord_achr_scl[Y_] = 0e0;

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

  lat_constr.high_ord_achr_scl[X_] = 1e3;
  lat_constr.high_ord_achr_scl[Y_] = 1e3;
  for (k = 0; k < 2; k++)
    lat_constr.high_ord_achr_nu[k] = high_ord_achr_nu[k];

  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ms"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ms"), 2));

  n = lat_constr.high_ord_achr_Fnum.size() - 1;
  lat_constr.high_ord_achr_dnu.resize(n);
  for (k = 0; k < n; k++)
    lat_constr.high_ord_achr_dnu[k].resize(2, 0e0);

  lat_constr.eps_x_scl = 1e2; lat_constr.eps0_x = 0.190;

  lat_constr.mI_scl[X_] = lat_constr.mI_scl[Y_] = 1e-10;
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

  lat_constr.high_ord_achr_scl[X_] = 0e0;
  lat_constr.high_ord_achr_scl[Y_] = 0e0;
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

  lat_constr.high_ord_achr_scl[X_] = 1e3;
  lat_constr.high_ord_achr_scl[Y_] = 1e3;
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

  lat_constr.mI_scl[X_] = lat_constr.mI_scl[Y_] = 1e-10;
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

  lat_constr.high_ord_achr_scl[X_] = 1e2;
  lat_constr.high_ord_achr_scl[Y_] = 1e2;
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

  lat_constr.mI_scl[X_] = lat_constr.mI_scl[Y_] = 1e-10;
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

  lat_constr.high_ord_achr_scl[X_] = 1e2;
  lat_constr.high_ord_achr_scl[Y_] = 1e2;
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

  lat_constr.mI_scl[X_] = lat_constr.mI_scl[Y_] = 1e-10;
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

  globval.mat_meth = false;

  // disable from TPSALib- and LieLib log messages
  idprset_(-1);

  if (!true)
    Read_Lattice(argv[1]);
  else {
    rdmfile(argv[1]);
    globval.dPcommon = 1e-6;
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
    dnu[X_] = 0.0; dnu[Y_] = -0.05;
    set_map(ElemIndex("ps_rot"), dnu);
    Ring_GetTwiss(true, 0e0); printglob();
  }

  if (false) {
    if (true) {
      Ring_GetTwiss(true, 0e0); printglob();
    } else
      ttwiss(lat_constr.ic[0], lat_constr.ic[1],
	     lat_constr.ic[2], lat_constr.ic[3], 0e0);

    eps_x = get_eps_x1(true, eps_x, sigma_E);
    printf("\neps_x = %5.3f nm.rad\n\n", eps_x);

    prt_lat("linlat1.out", globval.bpm, true);
    prt_lat("linlat.out", globval.bpm, true, 10);
    prt_chrom_lat("chromlat.out");
  }

  if (false) fit_ksi1(0e0, 0e0);

  switch (opt_case) {
  case 1:
    // Optimize Standard Straight: mI & 2*nu.
    opt_mI_std(lat_prms, lat_constr);
    no_sxt();
    if (true)
      fit_powell(lat_prms, 1e-3, f_achrom);
    else
      fit_sim_anneal(lat_prms, 1e-3, f_achrom);
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
    // fit_sim_anneal(lat_prms, 1e-3, f_achrom);
    break;
  case 4:
    // Match Short Straight: mI.
    match_ss(lat_prms, lat_constr);
    fit_powell(lat_prms, 1e-3, f_match);
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
