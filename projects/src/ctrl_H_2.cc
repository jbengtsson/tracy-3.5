
double bn_internal(const double bn_bounded,
		   const double bn_min, const double bn_max);
double bn_bounded(const double bn_internal,
		  const double bn_min, const double bn_max);
void set_ds(const int Fnum, const double ds);
void set_phi(const int Fnum, const double phi);
void set_dphi(const int Fnum, const double dphi);
void get_S(void);

double rad2deg(const double a) { return a*180e0/M_PI; }
double deg2rad(const double a) { return a*M_PI/180e0; }

void f_der(double *bn, double *df);
double f_match(double *bn);
double f_achrom(double *bn);
double f_mult(double *bn);


double eps_x, sigma_E;


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
  n_loc, n_constr, n_b1, n_b3, n_b4,
    n_iter;
  double
    ic[4][2],
    chi2, chi2_prt,
    eps_x_scl,
    eps0_x,               // Hor. emittance [nm.rad].
    nu[2],                // Cell tune.
    nu_ref_scl,
    nu_ref[2],            // Desired cell tune.
    ksi1_ctrl_scl[3],
    phi_scl,
    phi_tot, phi0,        // Cell bend angle.
    high_ord_achr_scl[2],
    mI_scl[2],
    mI0[2],               // -I Transformer.
    alpha_c_scl,          // alpha_c.
    L_scl,
    L0;                   // Cell length.
  std::vector<double> 
    ksi1_ctrl;
  std::vector< std::vector<double> >
    value,
    value_scl,
    mI,
    high_ord_achr_nu,  // Higher-Order-Achromat Phase Advance.
    high_ord_achr_dnu;
  double **Jacobian;
  std::vector<int>
    Fnum,
    Fnum_b1,
    Fnum_b3,
    Fnum_b4,
    high_ord_achr_Fnum,
    loc,
    type;
  std::vector< std::vector<int> >
    grad_dip_Fnum_b1,
    mI_loc;
  drv_terms_type
    drv_terms;

  constr_type(void) {
    n_iter = 0;
    chi2 = chi2_prt = 1e30;
    eps_x_scl = nu_ref_scl = phi_scl = 0e0;
    ksi1_ctrl_scl[0] = ksi1_ctrl_scl[1] = ksi1_ctrl_scl[2] = 0e0;
    high_ord_achr_scl[X_] = high_ord_achr_scl[Y_] = 0e0;
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
  double get_chi2(const double twoJ[], const double delta,
		  const double twoJ_delta[], const param_type &prms, double *bn,
		  const bool comp_drv_terms, const bool prt);
  void get_dchi2(const double twoJ[], const double delta,
		 const double twoJ_delta[], double *bn, double *df) const;
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
    if (!Cell[loc].Elem.Reverse) {
      Cell[loc].Elem.M->PTx1 = phi0;
      Cell[loc].Elem.M->PTx2 = phi1;
    } else {
      Cell[loc].Elem.M->PTx1 = phi1;
      Cell[loc].Elem.M->PTx2 = phi0;
    }
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

  const string
    labels[] = {"phi", "L  ", "s  ", "   ", "   ", "b_2", "b_3", "b_4", "b_5"};
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
    printf("   %3s %10.3e [%10.3e, %10.3e]",
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
  n_b1 = Fnum_b1.size();
  n_b3 = Fnum_b3.size();
  n_b4 = Fnum_b4.size();
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


void get_high_ord_achr(constr_type &constr)
{
  int j, k;

  for (j = 0; j < (int)constr.high_ord_achr_Fnum.size()/2; j++)
    for (k = 0; k < 2; k++)
      constr.high_ord_achr_dnu[j][k] =
	Cell[constr.high_ord_achr_Fnum[2*j+1]].Nu[k]
	- Cell[constr.high_ord_achr_Fnum[2*j]].Nu[k];
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


double constr_type::get_chi2(const double twoJ[], const double delta,
			     const double twoJ_delta[], const param_type &prms,
			     double *bn, const bool comp_drv_terms,
			     const bool prt)
{
  int    j, k;
  double chi2, dchi2[3], dnu[2], bn_ext;

  const bool
    extra     = false;
  const double
    delta_eps = 1e-6,
    scl_extra = 1e8;

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
    if (prt) printf("  eps_x          =  %10.3e\n", dchi2[X_]);
  }

  if (alpha_c_scl != 0e0) {
    dchi2[X_] = alpha_c_scl*sqr(1e0/globval.Alphac);
    chi2 += dchi2[X_];
    if (prt) printf("  alpha_c        =  %10.3e\n", dchi2[X_]);
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
      if (prt) printf("    alpha        = [%10.3e, %10.3e]\n",
		      dchi2[X_], dchi2[Y_]);
    }
    if ((value_scl[j][2] != 0e0) || (value_scl[j][3] != 0e0)) {
      for (k = 0; k < 2; k++) {
	dchi2[k] = value_scl[j][k+2]*sqr(Cell[loc[j]].Beta[k]-value[j][k+2]);
	chi2 += dchi2[k];
      }
      if (prt) printf("    beta         = [%10.3e, %10.3e]\n",
		      dchi2[X_], dchi2[Y_]);
    }
    if ((value_scl[j][4] != 0e0) || (value_scl[j][5] != 0e0)) {
	dchi2[0] = value_scl[j][4]*sqr(Cell[loc[j]].Eta[X_]-value[j][4]);
	dchi2[1] = value_scl[j][5]*sqr(Cell[loc[j]].Etap[X_]-value[j][5]);
	chi2 += dchi2[0]; chi2 += dchi2[1];
	if (prt) printf("    eta_x eta'_x = [%10.3e, %10.3e]\n",
			dchi2[0], dchi2[1]);
    }
  }

  if (L_scl != 0e0) {
    dchi2[0] = L_scl*sqr(Cell[globval.Cell_nLoc].S-L0);
    chi2 += dchi2[0];
    printf("\n  S-L0           = %10.3e\n", dchi2[0]);
  }

  if (comp_drv_terms) {
    lat_constr.drv_terms.get_h(delta_eps, twoJ, delta, twoJ_delta);

    if (prt) printf("\n");
    for (k = 0; k < (int)drv_terms.h.size(); k++) {
      dchi2[k] = lat_constr.drv_terms.h_scl[k]*lat_constr.drv_terms.h[k];
      chi2 += dchi2[k];
      if (prt && (lat_constr.drv_terms.h_scl[k] != 0e0)) {
	if ((k == 1) || (k == 4) || (k == 9) || (k == 20) || (k == 22)
	    || (k == 24)) printf("\n");
	printf("  h[%2d] = %10.3e\n", k, dchi2[k]);
      }
    }
    if (prt) lat_constr.drv_terms.print();
  }

  if (nu_ref_scl != 0e0) {
    if (prt) printf("\n");
    for (k = 0; k < 2; k++) {
      dnu[k] = nu_ref_scl*sqr(nu[k]-nu_ref[k]);
      chi2 += dnu[k];
    }
    if (prt) printf("  nu_ref          =  [%10.3e, %10.3e]\n",
		    dnu[X_], dnu[Y_]);
  }

  if ((ksi1_ctrl_scl[0] != 0e0) || (ksi1_ctrl_scl[1] != 0e0)
      || (ksi1_ctrl_scl[2] != 0e0)) {
    get_ksi1_ctrl(lat_constr);
    for (k = 0; k < 3; k++) {
      dchi2[k] = ksi1_ctrl_scl[k]/sqr(ksi1_ctrl[k]);
      chi2 += dchi2[k];
    }
    if (prt) printf("  ksi1_ctrl =       [%10.3e, %10.3e, %10.3e]\n",
		    dchi2[0], dchi2[1], dchi2[2]);
  }
  
  if ((mI_scl[X_] != 0e0) || (mI_scl[Y_] != 0e0)) {
    get_mI(lat_constr);
    if (prt) printf("\n");
    for (j = 0; j < (int)mI.size(); j++) {
      for (k = 0; k < 2; k++) {
	dchi2[k] = mI_scl[k]*sqr(mI[j][k]-mI0[k]);
	chi2 += dchi2[k];
      }
      if (prt) printf("  mI              =  [%10.3e, %10.3e]\n",
		      dchi2[X_], dchi2[Y_]);
    }
  }

  if ((high_ord_achr_scl[X_] != 0e0) || (high_ord_achr_scl[Y_] != 0e0)) {
    get_high_ord_achr(lat_constr);
    if (prt) printf("\n");
    for (j = 0; j < (int)high_ord_achr_dnu.size(); j++) {
      for (k = 0; k < 2; k++) {
	dchi2[k] =
	  high_ord_achr_scl[k]
	  *sqr(high_ord_achr_dnu[j][k]-high_ord_achr_nu[j][k]);
	chi2 += dchi2[k];
      }
      if (prt) printf("  high_ord_achr   =  [%10.3e, %10.3e]\n",
		      dchi2[X_], dchi2[Y_]);
    }
  }

  return chi2;
}


void constr_type::get_dchi2(const double twoJ[], const double delta,
			    const double twoJ_delta[], double *bn,
			    double *df) const
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
    df[k+1] =
      lat_constr.get_chi2(twoJ, delta, twoJ_delta, lat_prms, bn, true, false);

    constr_dparam(lat_prms.Fnum[k], lat_prms.n[k], -2e0*eps);
    eps_x = get_lin_opt(lat_constr);
    df[k+1] -=
      lat_constr.get_chi2(twoJ, delta, twoJ_delta, lat_prms, bn, true, false);
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

  printf("    b_3L    = [");
  for (k = 0; k < constr.n_b3; k++) {
    get_bnL_design_elem(constr.Fnum_b3[k], 1, Sext, b3L, a3L);
    printf("%10.3e", b3L);
    if (k != constr.n_b3-1) printf(", ");
  }
  printf("]\n");
  printf("    b_3     = [");
  for (k = 0; k < constr.n_b3; k++) {
    get_bn_design_elem(constr.Fnum_b3[k], 1, Sext, b3, a3);
    printf("%10.3e", b3);
    if (k != constr.n_b3-1) printf(", ");
  }
  printf("]\n");
}


void prt_high_ord_achr(const constr_type &constr)
{
  int j, n;

  n = constr.high_ord_achr_nu.size();
  printf("\n  Higher-Order-Achromat:\n");
  for (j = 0; j < n; j++)
    printf("    [%7.5f, %7.5f]\n",
	   constr.high_ord_achr_nu[j][X_], constr.high_ord_achr_nu[j][Y_]);
  printf("\n");
  for (j = 0; j < n; j++)
    printf("    [%7.5f, %7.5f]\n",
	   constr.high_ord_achr_dnu[j][X_], constr.high_ord_achr_dnu[j][Y_]);
}


void constr_type::prt_constr(const double chi2)
{
  int    loc, k;
  double phi;

  printf("\n%3d chi2: %12.6e -> %12.6e\n", n_iter, this->chi2_prt, chi2);
  this->chi2_prt = chi2;
  printf("\n  Linear Optics:\n");
  printf("    eps_x   = %5.3f (%5.3f)\n"
	 "    nu      = [%5.3f, %5.3f]\n"
	 "    dnu     = [%5.3f, %5.3f]\n",
	 eps_x, eps0_x, nu[X_], nu[Y_], nu[X_]-nu_ref[X_], nu[Y_]-nu_ref[Y_]);
  if (ring)
    printf("    ksi1    = [%5.3f, %5.3f]\n"
	   "    alpha_c = %9.3e\n",
	   lat_constr.drv_terms.h_c[0], lat_constr.drv_terms.h_s[0],
	   globval.Alphac);
  if (lat_constr.ksi1_ctrl.size() != 0)
    printf("    ksi1_ctrl    = [%5.3f, %5.3f, %5.3f]\n",
	   lat_constr.ksi1_ctrl[0], lat_constr.ksi1_ctrl[1],
	   lat_constr.ksi1_ctrl[2]);
  if (phi_scl != 0e0) {
    loc = Elem_GetPos(Fnum_b1[n_b1-1], 1);
    phi = rad2deg(Cell[loc].Elem.PL*Cell[loc].Elem.M->Pirho);
    printf("    phi     = %7.5f (%7.5f)\n    ", phi_tot, phi0);
    prt_name(stdout, Cell[loc].Elem.PName, "_phi", 3);
    printf(" = %7.5f\n", phi);
  }
  if (L_scl != 0e0)
    printf("    L       = %7.5f (%7.5f)\n", Cell[globval.Cell_nLoc].S, L0);

  prt_h(*this);

  if ((high_ord_achr_scl[X_] != 0e0) || (high_ord_achr_scl[Y_] != 0e0))
    prt_high_ord_achr(*this);

  if ((mI_scl[X_] != 0e0) || (mI_scl[Y_] != 0e0)) {
    printf("\n  -I Transf.:\n");
    printf("    [%8.5f, %8.5f]\n\n", mI0[X_], mI0[Y_]);
    for (k = 0; k < (int)mI.size(); k++)
      printf("    [%8.5f, %8.5f] [%3d, %3d]\n",
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


void fit_ksi1(const std::vector<int> &Fnum_b3,
	      const double ksi_x, const double ksi_y, const double db3L)
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
  // no_sxt();

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
  for (j = 1; j <= 2; j++)
    b[j] = -(globval.Chrom[j-1]-ksi0[j-1]);

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
  if (!Cell.Elem.Reverse)
    fprintf(outf, " bending, l = %11.8f,"
	    "\n    t = %11.8f, t1 = %11.8f, t2 = %11.8f, k = %11.8f,"
	    "\n    n = nbend, method = 4;\n",
	    Cell.Elem.PL, phi, Cell.Elem.M->PTx1,
	    Cell.Elem.M->PTx2, Cell.Elem.M->PBpar[Quad+HOMmax]);
  else {
    fprintf(outf, " bending, l = %11.8f,"
	    "\n    t = %11.8f, t1 = %11.8f, t2 = %11.8f, k = %11.8f,"
	    "\n    n = nbend, method = 4;\n",
	    Cell.Elem.PL, phi, Cell.Elem.M->PTx2,
	    Cell.Elem.M->PTx1, Cell.Elem.M->PBpar[Quad+HOMmax]);
  }
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
  int      k, loc;
  double   b_n_prev, b_n;
  CellType *cellp;

 const bool prt = false;

  loc = Elem_GetPos(abs(lat_prms.Fnum[n-1]), 1);
  cellp = &Cell[loc];
  if (prt) {
    printf("prt_elem:\n  ");
    prt_name(stdout, cellp->Elem.PName, "\n", 0);
  }
  if (cellp->Elem.Pkind == drift) {
    prt_name(outf, cellp->Elem.PName, ":", 8);
    prt_drift(outf, Cell[loc]);
  } else if (cellp->Elem.Pkind == Mpole) {
    if (prt) {
      printf("  n_design:     %1d\n", cellp->Elem.M->n_design);
      printf("  n:            %1d\n", cellp->Elem.M->Porder);
    }
    if (cellp->Elem.M->n_design == Dip) {
      if (lat_prms.Fnum[n-1] > 0) {
	prt_name(outf, cellp->Elem.PName, ":", 8);
	prt_dip(outf, Cell[loc]);
      } else
      	prt_grad_dip(outf, lat_prms.grad_dip_Fnum[n-1]);
    } else if (cellp->Elem.M->n_design == Quad) {
      prt_name(outf, cellp->Elem.PName, ":", 8);
      fprintf(outf, " quadrupole, l = %10.8f, k = %13.8f, n = nquad"
	      ", method = 4;\n",
	      cellp->Elem.PL, cellp->Elem.M->PBpar[Quad+HOMmax]);
    } else if (cellp->Elem.M->Porder == Sext) {
      prt_name(outf, cellp->Elem.PName, ":", 8);
      fprintf(outf, " sextupole,  l = %10.8f, k = %13.8f, n = nsext"
	      ", method = 4;\n",
	      cellp->Elem.PL, cellp->Elem.M->PBpar[Sext+HOMmax]);
    } else {
      prt_name(outf, cellp->Elem.PName, ":", 8);
      fprintf(outf, " multipole,  l = %10.8f,\n          hom = (",
	      cellp->Elem.PL);
      for (k = Sext; k <= cellp->Elem.M->Porder; k++){
	b_n_prev = 0e0;
	b_n = cellp->Elem.M->PBpar[k+HOMmax];
	if (b_n != 0e0) {
	  if ((k > Sext) && (b_n_prev != 0e0))
	    fprintf(outf, "                 ");
	  fprintf(outf, "%1d, %14.8e, 0.0", k, b_n);
	  if (k < cellp->Elem.M->Porder)
	    fprintf(outf, ",\n");
	  else
	    fprintf(outf, "),\n");
	}
	b_n_prev = b_n;
      }
      fprintf(outf, "          n = nsext, method = 4;\n");
    }
  } else {
    printf("\nprt_elem: %s %d\n", cellp->Elem.PName, cellp->Elem.Pkind);
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


void prt_bn(FILE *outf, const param_type &lat_prms)
{
  long int         loc;
  int              k;
  std::vector<int> Fam_prt;

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
}


void prt_b3(FILE *outf, const std::vector<int> &Fnum)
{
  long int loc;
  int      k;

  for (k = 0; k < (int)Fnum.size(); k++) {
    loc = Elem_GetPos(Fnum[k], 1);
    prt_name(outf, Cell[loc].Elem.PName, ":", 8);
    fprintf(outf, " sextupole,  l = %10.8f, k = %13.8f, n = nsext"
	    ", method = 4;\n",
	    Cell[loc].Elem.PL, Cell[loc].Elem.M->PBpar[Sext+HOMmax]);
  }
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
      printf("\nget_nu: unstable in %s plane %10.3e\n",
	     (k == 0)? "hor" : "ver", cosmu);
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

  int          n_bn, i, j, iter;
  double       *bn, **xi, fret, eps_x, w[n_w], x[n_prm];
  ss_vect<tps> A;

  const double ftol = 1e-8;

  n_bn = lat_prms.n_prm;

  bn = dvector(1, n_bn); xi = dmatrix(1, n_bn, 1, n_bn);

  lat_prms.ini_prm(bn);
  f(bn);

  if (false) {
    lat_constr.get_Jacobian(lat_prms);
    lat_constr.prt_Jacobian(n_bn);
    exit(0);
  }

  if (true) {
    // Set initial directions (unit vectors).
    for (i = 1; i <= n_bn; i++)
      for (j = 1; j <= n_bn; j++)
	xi[i][j] = (i == j)? eps : 0e0;

    dpowell(bn, xi, n_bn, ftol, &iter, &fret, f);
  } else {
    // for (i = 1; i <= n_prm; i++)
    //   x[i-1] = bn[i];
    // newuoa_(n_prm, n_pt, x, rho_beg, rho_end, n_prt, max_fun, w);
  }

  printf("\n  iter = %d fret = %12.5e\n", iter, fret);
  printf("bns:\n");
  lat_prms.prt_prm(bn);
  // lat_prms.set_prm(bn);
  // eps_x = get_lin_opt(false);
  // f_prt(bn);

  free_dvector(bn, 1, n_bn); free_dmatrix(xi, 1, n_bn, 1, n_bn);
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

  int          n_bn, i, j, iter;
  double       *bn, **xi, fret, eps_x, w[n_w], x[n_prm];
  ss_vect<tps> A;

  const double ftol = 1e-8;

  n_bn = lat_prms.n_prm;

  bn = dvector(1, n_bn); xi = dmatrix(1, n_bn, 1, n_bn);

  lat_prms.ini_prm(bn);
  f(bn);

  if (false) {
    lat_constr.get_Jacobian(lat_prms);
    lat_constr.prt_Jacobian(n_bn);
    exit(0);
  }

  if (true) {
    // Set initial directions (unit vectors).
    for (i = 1; i <= n_bn; i++)
      for (j = 1; j <= n_bn; j++)
	xi[i][j] = (i == j)? eps : 0e0;

    dpowell(bn, xi, n_bn, ftol, &iter, &fret, f);
  } else {
    // for (i = 1; i <= n_prm; i++)
    //   x[i-1] = bn[i];
    // newuoa_(n_prm, n_pt, x, rho_beg, rho_end, n_prt, max_fun, w);
  }

  printf("\n  iter = %d fret = %12.5e\n", iter, fret);
  printf("bns:\n");
  lat_prms.prt_prm(bn);
  // lat_prms.set_prm(bn);
  // eps_x = get_lin_opt(false);
  // f_prt(bn);

  free_dvector(bn, 1, n_bn); free_dmatrix(xi, 1, n_bn, 1, n_bn);
}


void fit_conj_grad(param_type &lat_prms, double (*f)(double *))
{
  int          n_bn, iter;
  double       *bn, fret, eps_x;
  ss_vect<tps> A;

  const double ftol = 1e-8;

  n_bn = lat_prms.n_prm;

  bn = dvector(1, n_bn);

  lat_prms.ini_prm(bn);
  f(bn);

  dfrprmn(bn, n_bn, ftol, &iter, &fret, f, f_der);

  printf("\n  iter = %d fret = %12.5e\n", iter, fret);
  printf("bns:\n");
  lat_prms.prt_prm(bn);
  // lat_prms.set_prm(bn);
  // eps_x = get_lin_opt(false);
  // f_prt(bn);

  free_dvector(bn, 1, n_bn);
}


void prt_f(double *bn, const double chi2, constr_type &lat_constr,
	   param_type &lat_prms, const bool chrom)
{
  FILE *outf;

  const std::string file_name = "opt_des_2_bn.out";

  outf = file_write(file_name.c_str());

  lat_constr.prt_constr(chi2);

  printf("\n  Parameters:\n");
  lat_prms.prt_prm(bn);
  prtmfile("flat_file.fit");
  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);
  if (chrom) prt_chrom_lat();

  prt_bn(outf, lat_prms);
  prt_b3(outf, lat_constr.Fnum_b3);

  fclose(outf);
}


void prt_prms(constr_type &constr)
{
  printf("\n  eps_x_scl         =  %9.3e\n"
	 "  mI_nu_ref         = [%7.5f,   %7.5f]\n"
	 "  mI_scl            = [%9.3e, %9.3e]\n"
	 "  high_ord_achr_scl = [%9.3e, %9.3e]\n"
	 "  phi_scl           =  %9.3e\n",
	 constr.eps_x_scl,
	 constr.mI0[X_], constr.mI0[Y_], constr.mI_scl[X_], constr.mI_scl[Y_],
	 constr.high_ord_achr_scl[X_], constr.high_ord_achr_scl[Y_],
	 constr.phi_scl);
}
