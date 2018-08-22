#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


double rad2deg(const double a) { return a*180e0/M_PI; }

double deg2rad(const double a) { return a*M_PI/180e0; }


// Standard Cell.
const double ic[][2] =
  {{0.0, 0.0}, {4.97103, 5.52181}, {0.0, 0.0}, {0.0, 0.0}};


double eps_x;


double bn_internal(const double bn_bounded,
		   const double bn_min, const double bn_max);
double bn_bounded(const double bn_internal,
		  const double bn_min, const double bn_max);
void set_ds(const int Fnum, const double ds);
void set_bend(const int Fnum, const double phi);
void get_S(void);


struct param_type {
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

  void add_prm(const std::string Fname, const int n,
	       const double bn_min, const double bn_max,
	       const double bn_scl);
  void ini_prm(double *bn);
  void set_prm(double *bn);
  void prt_prm(double *bn) const;
};


struct constr_type {
private:

public:
  int
    n_loc, n_b3, n_b1,
    n_iter;
  double
    chi2,
    eps_x_scl,
    eps0_x,       // Hor. emittance [nm.rad].
    drv_terms_scl,
    nu[2],
    phi,
    phi0,         // Total bend angle.
    drv_terms[2];
  std::vector< std::vector<double> >
    value,
    scl;
  std::vector<int>
    Fnum,
    Fnum_b3,
    Fnum_b1,
    loc,
    type;

  void add_constr(const int loc,
		  const double scl1, const double scl2, const double scl3,
		  const double scl4, const double scl5, const double scl6,
		  const double v1, const double v2, const double v3,
		  const double v4, const double v5, const double v6);
  void ini_constr(const double eps_x_scl, const double eps0_x,
		  const double drv_terms_scl, const double phi0);
  double get_chi2(void) const;
  void prt_Jacobian(void) const;
  void prt_constr(const double chi2) const;
};


param_type  lat_prms;
constr_type lat_constr;


void prt_name(FILE *outf, const char *name, const int cut)
{
  int j, k, len;

  len = strlen(name);
  j = 0;
  do {
    fprintf(outf, "%c", name[j]);
    j++;
  } while ((j < len) && (name[j] != ' '));
  fprintf(outf, ":");
  if (cut == -1) len = cut;
  for (k = j; k < len; k++)
    fprintf(outf, "%c", name[k]);
}


void param_type::add_prm(const std::string Fname, const int n,
			 const double bn_min, const double bn_max,
			 const double bn_scl)
{
  Fnum.push_back(ElemIndex(Fname.c_str()));
  this->n.push_back(n);
  this->bn_min.push_back(bn_min);
  this->bn_max.push_back(bn_max);
  this->bn_scl.push_back(bn_scl);
  this->l0.push_back(0e0);
  n_prm = Fnum.size();
}


void param_type::ini_prm(double *bn)
{
  int    i, loc;
  double an;

  for (i = 1; i <= n_prm; i++) {
    if (n[i-1] > 0)
      // Multipole.
      get_bn_design_elem(Fnum[i-1], 1, n[i-1], bn[i], an);
    else if (n[i-1] == -1) {
      // Displacement.
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
      loc = Elem_GetPos(Fnum[i-1], 1);
      bn[i] = rad2deg(Cell[loc].Elem.PL*Cell[loc].Elem.M->Pirho);
    }
    // Bounded.
    bn[i] = bn_internal(bn[i], bn_min[i-1], bn_max[i-1]);
  }
}


void param_type::set_prm(double *bn)
{
  int    i;
  double bn_ext;

  const bool prt   = false;
  const int  n_prt = 5;

  if (prt) printf("\nset_prm:\n  ");
  for (i = 1; i <= n_prm; i++) {
    // Bounded.
    bn_ext = bn_bounded(bn[i], bn_min[i-1], bn_max[i-1]);
    if (n[i-1] > 0)
      // Multipole strength.
      set_bn_design_fam(Fnum[i-1], n[i-1], bn_ext, 0e0);
    else if (n[i-1] == -1) {
      // Displacement.
      set_ds(Fnum[i-1], bn_ext-l0[i-1]);
      l0[i-1] = bn_ext;
    } else if (n[i-1] == -2) {
      // Length.
      set_L(Fnum[i-1], bn_ext); get_S();
    } else if (n[i-1] == -3)
      // Bend angle; L is fixed.
      set_bend(Fnum[i-1], bn_ext);
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

  string labels[] = {"phi", "L  ", "s  ", "   ", "   ", "b_2"};

  for (i = 1; i <= n_prm; i++) {
    // Bounded.
    bn_ext = bn_bounded(bn[i], bn_min[i-1], bn_max[i-1]);
    loc = Elem_GetPos(Fnum[i-1], 1);
    printf("    ");
    prt_name(stdout, Cell[loc].Elem.PName, 4);
    printf("   %3s %9.5f [%9.5f, %9.5f]",
	   labels[n[i-1]+3].c_str(), bn_ext, bn_min[i-1], bn_max[i-1]);
    if (n[i-1] == -1) {
      printf(" ");
      prt_name(stdout, Cell[loc-1].Elem.PName, -1);
      printf(" %7.5f, ", Cell[loc-1].Elem.PL);
      prt_name(stdout, Cell[loc+1].Elem.PName, -1);
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

  if (Knum % 2 == 1) {
    set_dL(Cell[loc+1].Fnum, Knum, -ds);
    set_dL(Cell[loc-1].Fnum, Knum, ds);
  } else {
    set_dL(Cell[loc+1].Fnum, Knum, ds);
    set_dL(Cell[loc-1].Fnum, Knum, -ds);
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


void set_bend(const int Fnum, const double phi)
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
  this->scl.push_back(vec1);
  this->value.push_back(vec2);
  this->n_loc = this->loc.size();
}


void constr_type::ini_constr(const double eps_x_scl, const double eps0_x,
			     const double drv_terms_scl, const double phi0)
{
  n_b3 = Fnum_b3.size();
  n_b1 = Fnum_b1.size();
  this->eps_x_scl = eps_x_scl;
  this->eps0_x = eps0_x;
  this->drv_terms_scl = drv_terms_scl;
  this->phi0 = phi0;

  n_iter = 0; chi2 = 1e30;
}


double constr_type::get_chi2(void) const
{
  int    k;
  double chi2;

  chi2 = 0e0;

  chi2 += 1e2*sqr(eps_x-eps0_x);

  for (k = 0; k < n_loc; k++) {
    chi2 += scl[k][0]*sqr(Cell[loc[k]].Alpha[X_]-value[k][0]);
    chi2 += scl[k][1]*sqr(Cell[loc[k]].Alpha[Y_]-value[k][1]);
    chi2 += scl[k][2]*sqr(Cell[loc[k]].Beta[X_]-value[k][2]);
    chi2 += scl[k][3]*sqr(Cell[loc[k]].Beta[Y_]-value[k][3]);
    chi2 += scl[k][4]*sqr(Cell[loc[k]].Eta[X_]-value[k][4]);
    chi2 += scl[k][5]*sqr(Cell[loc[k]].Etap[X_]-value[k][5]);
  }

  chi2 += 1e-7*drv_terms[X_];
  chi2 += 1e-7*drv_terms[Y_];

  return chi2;
}


void constr_type::prt_constr(const double chi2) const
{
  int           k, loc;
  double        b3L, a3L;
  static double chi2_ref;

  printf("\n%3d chi2: %11.5e -> %11.5e\n", n_iter, chi2_ref, chi2);
  chi2_ref = chi2;

  printf("\n  Linear Optics:\n");
  printf("    eps_x     = %5.3f\n"
	 "    phi       = %7.5f\n"
	 "    nu        = [%5.3f, %5.3f]\n"
	 "    drv terms = [%10.3e, %10.3e]\n",
	 eps_x, phi, nu[X_], nu[Y_],
	 sqrt(drv_terms[X_]), sqrt(drv_terms[Y_]));

  printf("    b_3       = [");
  for (k = 0; k < n_b3; k++) {
    get_bnL_design_elem(Fnum_b3[k], 1, Sext, b3L, a3L);
    printf("%10.3e", b3L);
    if (k != n_b3-1) printf(", ");
  }
  printf("]\n");

  printf("\n     alpha_x   alpha_y  beta_x   beta_y    eta_x    eta'_x\n");
  for (k = 0; k < n_loc; k++) {
    loc = this->loc[k];
    printf("    %8.5f  %8.5f %8.5f %8.5f %8.5f %8.5f\n",
	   Cell[loc].Alpha[X_], Cell[loc].Alpha[Y_],
	   Cell[loc].Beta[X_], Cell[loc].Beta[Y_],
	   Cell[loc].Eta[X_], Cell[loc].Etap[X_]);
  }
}


void get_S(void)
{
  int k;

  Cell[0].S = 0e0;
  for (k = 1; k <= globval.Cell_nLoc; k++)
    Cell[k].S = Cell[k-1].S + Cell[k].Elem.PL;
}


double get_eps_x1(const bool track)
{
  // eps_x [pm.rad].
  long int     lastpos;
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

  if (prt) {
    printf("\neps_x = %5.3f nm.rad\n", eps_x);
    printf("J_x   = %5.3f, J_z = %5.3f\n", 1.0-I4/I2, 2.0+I4/I2);
  }

  return 1470e0*sqr(globval.Energy)*I5/(I2-I4);
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


void prt_system(std::vector<string> &drv_terms, const int m, const int n,
		double **A, double *b)
{
  int i, j;

  printf("\n Ax = b:\n          ");
  for (j = 1; j <= n; j++)
    printf("%11d", j);
  printf("\n");
  for (i = 1; i <= m; i++) {
    printf("%4d %10s", i, drv_terms[i-1].c_str());
    for (j = 1; j <= n; j++)
      printf("%11.3e", A[i][j]);
    printf("%11.3e\n", b[i]);
  }
}


void fit_ksi1(const std::vector<int> &Fnum_b3,
	      const double ksi_x, const double ksi_y, const double db3L)
{
  int                 n_b3, j, k;
  double              **A, **U, **V, *w, *b, *x, b3L, a3L;
  std::vector<string> drv_terms;

  const bool   prt = false;
  const int    m   = 2;
  const double
    ksi0[]  = {ksi_x, ksi_y},
    svd_cut = 1e-10;

  n_b3 = Fnum_b3.size();

  A = dmatrix(1, m, 1, n_b3); U = dmatrix(1, m, 1, n_b3);
  V = dmatrix(1, n_b3, 1, n_b3);
  w = dvector(1, n_b3); b = dvector(1, m); x = dvector(1, n_b3);

  for (k = 1; k <= n_b3; k++) {
    set_dbnL_design_fam(Fnum_b3[k-1], Sext, db3L, 0e0);
    Ring_Getchrom(0e0);
    for (j = 1; j <= m; j++)
      A[j][k] = globval.Chrom[j-1];
    set_dbnL_design_fam(Fnum_b3[k-1], Sext, -2e0*db3L, 0e0);
    Ring_Getchrom(0e0);
    for (j = 1; j <= 2; j++) {
      A[j][k] -= globval.Chrom[j-1]; A[j][k] /= 2e0*db3L;
    }
    set_dbnL_design_fam(Fnum_b3[k-1], Sext, db3L, 0e0);
  }
  Ring_Getchrom(0e0);
  for (j = 1; j <= 2; j++)
    b[j] = -(globval.Chrom[j-1]-ksi0[j-1]);
  drv_terms.push_back("ksi1_x"); drv_terms.push_back("ksi1_y");

  if (prt) prt_system(drv_terms, m, n_b3, A, b);

  dmcopy(A, m, n_b3, U); dsvdcmp(U, m, n_b3, w, V);
  if (prt) {
    printf("\nfit_ksi1:\n  singular values:\n  ");
    for (j = 1; j <= n_b3; j++) {
      printf("%11.3e", w[j]);
      if (w[j] < svd_cut) {
	w[j] = 0e0;
	printf(" (zeroed)");
      }
    }
    printf("\n");
  }
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


void prt_elem(FILE *outf, CellType &Cell)
{
  double phi;

  if (Cell.Elem.Pkind == drift)
    prt_drift(outf, Cell);
  else if (Cell.Elem.Pkind == Mpole) {
    if (Cell.Elem.M->n_design == Dip) {
      phi = rad2deg(Cell.Elem.PL*Cell.Elem.M->Pirho);
      fprintf(outf, " bending, l = %11.8f,"
	      "\n    t = %11.8f, t1 = %11.8f, t2 = %11.8f, k = %11.8f,"
	      "\n    n = nbend, method = 4;\n",
	      Cell.Elem.PL, phi, Cell.Elem.M->PTx1,
	      Cell.Elem.M->PTx2, Cell.Elem.M->PBpar[Quad+HOMmax]);
    } else if (Cell.Elem.M->n_design == Quad)
      fprintf(outf, " quadrupole, l = %11.8f, k = %11.8f, n = nquad"
	      ", method = 4;\n",
	      Cell.Elem.PL, Cell.Elem.M->PBpar[Quad+HOMmax]);
  } else {
    printf("\nprt_elem: %s %d\n", Cell.Elem.PName, Cell.Elem.Pkind);
    exit(1);
  }
}


void prt_b2(const param_type &lat_prms)
{
  long int loc;
  int      k;
  FILE     *outf;

  std::string file_name = "opt_des_1_b2.out";

  outf = file_write(file_name.c_str());

  for (k = 0; k < lat_prms.n_prm; k++) {
    loc = Elem_GetPos(lat_prms.Fnum[k], 1);
    prt_name(outf, Cell[loc].Elem.PName, -1);
    if (lat_prms.n[k] != -1)
      prt_elem(outf, Cell[loc]);
    else {
      fprintf(outf, "\n  ");
      prt_name(outf, Cell[loc+1].Elem.PName, -1);
      prt_drift(outf, Cell[loc+1]);
      fprintf(outf, "  ");
      prt_name(outf, Cell[loc-1].Elem.PName, -1);
      prt_drift(outf, Cell[loc-1]);
    }
  }

  fclose(outf);
}


void prt_b3(const constr_type &constr)
{
  long int loc;
  int      k;
  FILE     *outf;

  std::string file_name = "opt_des_1_b3.out";

  outf = file_write(file_name.c_str());

  for (k = 0; k < constr.n_b3; k++) {
    loc = Elem_GetPos(constr.Fnum_b3[k], 1);
    prt_name(outf, Cell[loc].Elem.PName, -1);
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
    if (cosmu < 1e0) {
      nu[k] = acos(cosmu)/(2e0*M_PI);
      if (ps[2*k][2*k+1] < 0e0) nu[k] = 1e0 - nu[k];
    } else {
      printf("\nget_nu: unstable in plane %d %7.5f\n", k, cosmu);
      return false;
    }
  }
  if (prt) printf("\nget_nu: nu = [%7.5f, %7.5f]\n", nu[X_], nu[Y_]);
  return true;
}


double get_lin_opt(const bool periodic)
{
  double       eps_x;
  ss_vect<tps> A;

  if (periodic) {
    Ring_GetTwiss(true, 0e0);
    eps_x = get_eps_x1(true);
  } else {
    globval.emittance = true;
    // ttwiss(ic[0], ic[1], ic[2], ic[3], 0e0);
    A = get_A(ic[0], ic[1], ic[2], ic[3]);
    Cell_Twiss(0, globval.Cell_nLoc, A, false, false, 0e0);
    eps_x = get_eps_x1(false);
    globval.emittance = false;
  }

  return eps_x;
}


void get_drv_terms(constr_type &constr)
{
  int    j, k, n_kid;
  double b3L, a3L;

  for (k = 0; k < 2; k++)
    constr.drv_terms[k] = 0e0;
  for (k = 0; k < constr.n_b3; k++) {
    get_bnL_design_elem(constr.Fnum_b3[k], 1, Sext, b3L, a3L);
    for (j = 0; j < 2; j++) {
      n_kid = GetnKid(constr.Fnum_b3[k]);
      constr.drv_terms[j] +=
	sqr(n_kid*b3L*Cell[Elem_GetPos(constr.Fnum_b3[k], 1)].Beta[j]);
    }
  }
}


void fit_powell(param_type &lat_prms, const double eps, double (*f)(double *))
{
  int          n_b2, i, j, iter;
  double       *b2, **xi, fret, eps_x;
  ss_vect<tps> A;

  n_b2 = lat_prms.n_prm;

  b2 = dvector(1, n_b2); xi = dmatrix(1, n_b2, 1, n_b2);

  lat_prms.ini_prm(b2);
  f(b2);

  // Set initial directions (unit vectors).
  for (i = 1; i <= n_b2; i++)
    for (j = 1; j <= n_b2; j++)
      xi[i][j] = (i == j)? eps : 0e0;

  dpowell(b2, xi, n_b2, 1e-16, &iter, &fret, f);

  printf("\n  fret = %12.5e\n", fret);
  printf("b2s:\n");
  lat_prms.prt_prm(b2);
  // lat_prms.set_prm(b2);
  // eps_x = get_lin_opt(false);
  // f_prt(b2);

  free_dvector(b2, 1, n_b2); free_dmatrix(xi, 1, n_b2, 1, n_b2);
}


void prt_achrom(double *b2, const double chi2, constr_type &lat_constr,
		param_type &lat_prms)
{
  lat_constr.prt_constr(chi2);

  printf("\n  Parameters:\n");
  lat_prms.prt_prm(b2);
  prtmfile("flat_file.fit");
  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);

  prt_b2(lat_prms);
  prt_b3(lat_constr);
}


void phi_corr(constr_type &constr)
{
  // Correct total bend angle.
  int    k, loc;
  double phi1;

  constr.phi = 0e0;
  for (k = 0; k < constr.n_b1; k++) {
    loc = Elem_GetPos(constr.Fnum_b1[k], 1);
    constr.phi +=
      GetnKid(constr.Fnum_b1[k])
      *rad2deg(Cell[loc].Elem.PL*Cell[loc].Elem.M->Pirho);
  }
  loc = Elem_GetPos(constr.Fnum_b1[constr.n_b1-1], 1);
  phi1 =
    rad2deg(Cell[loc].Elem.PL*Cell[loc].Elem.M->Pirho)
    - (constr.phi-constr.phi0)
    /GetnKid(constr.Fnum_b1[constr.n_b1-1]);
  set_bend(constr.Fnum_b1[constr.n_b1-1], phi1);
 
}


double f_achrom(double *b2)
{
  bool   stable;
  double chi2;

  const int n_prt = 10;

  lat_prms.set_prm(b2);
  phi_corr(lat_constr);

  stable = get_nu(lat_constr.nu);
  if (stable) {
    fit_ksi1(lat_constr.Fnum_b3, 0e0, 0e0, 1e-1);
    eps_x = get_lin_opt(true);
    get_drv_terms(lat_constr);

    chi2 = lat_constr.get_chi2();

    if (chi2 < lat_constr.chi2) {
      if (lat_constr.n_iter % n_prt == 0)
	prt_achrom(b2, chi2, lat_constr, lat_prms);
      lat_constr.n_iter++;
      lat_constr.chi2 = chi2;
    }
  } else
    chi2 = 1e30;

  return chi2;
}


void set_lat_prms_1(param_type &prms)
{
    // lat_prms.add_prm("sd",  -1,   0.07,   0.1,  1.0);

    prms.add_prm("qf1",  2, -20.0,   20.0,  1.0);
    prms.add_prm("qd2",  2, -20.0,   20.0,  1.0);
    prms.add_prm("d_3", -2,   0.2,    0.45, 1.0);
    prms.add_prm("qf3",  2, -20.0,   20.0,  1.0);

    prms.add_prm("b1",  -3, -20.0,   20.0,  1.0);
    prms.add_prm("b1",  -2,   0.1,    0.8,  1.0);
    prms.add_prm("b1",  -1,   0.075,  0.07, 1.0);
    prms.add_prm("b1",   2, -20.0,   20.0,  1.0);

    prms.add_prm("b2",  -2,   0.1,    0.5,  1.0);
    prms.add_prm("b2",   2, -20.0,   20.0,  1.0);

    prms.add_prm("qf4",  2, -20.0,   20.0,  1.0);
    prms.add_prm("qf4", -1,   0.1,    0.1,  1.0);

    // Parameters are initialized in optimizer.
}


void set_lat_constr_1(constr_type &constr)
{
  // Lattice constraints are: alpha_x,y, beta_x,y, eta_x, eta'_x.
  constr.add_constr(Elem_GetPos(ElemIndex("sfh"), 1),
			1e1, 1e1, 0e0, 0e0, 1e0,  0e0,
			0e0, 0e0, 0e0, 0e0, 6e-2, 0e0);
  constr.add_constr(Elem_GetPos(ElemIndex("b2"), 1),
			1e1, 1e1, 0e0, 0e0, 0e0, 1e1,
			0e0, 0e0, 0e0, 0e0, 0e0, 0e0);
  constr.add_constr(Elem_GetPos(ElemIndex("b1"), 2),
			0e0, 0e0, 0e0, 0e0, 1e5, 1e5,
			0e0, 0e0, 0e0, 0e0, 0e0, 0e0);
  constr.add_constr(globval.Cell_nLoc,
			1e1, 1e1, 1e-1, 1e-1, 0e0, 0e0,
			0e0, 0e0, 4e0,  2.5e0,  0e0, 0e0);

  lat_prms.bn_tol = 1e-6; lat_prms.step = 1.0;

  lat_constr.Fnum_b3.push_back(ElemIndex("sfh"));
  lat_constr.Fnum_b3.push_back(ElemIndex("sd"));

  lat_constr.Fnum_b1.push_back(ElemIndex("b1"));
  lat_constr.Fnum_b1.push_back(ElemIndex("b2"));

  lat_constr.ini_constr(1e2, 0.190, 1e-7, 7.5);
}


void set_lat_prms_2(param_type &prms)
{
    // lat_prms.add_prm("sd",  -1,   0.07,   0.1,  1.0);

    prms.add_prm("qf1",  2, -20.0,   20.0,  1.0);
    prms.add_prm("qd2",  2, -20.0,   20.0,  1.0);
    prms.add_prm("d_3", -2,   0.2,    0.45, 1.0);
    prms.add_prm("qf3",  2, -20.0,   20.0,  1.0);

    prms.add_prm("b1",  -3, -20.0,   20.0,  1.0);
    prms.add_prm("b1",  -2,   0.1,    1.2,  1.0);
    prms.add_prm("b1",  -1,   0.075,  0.07, 1.0);
    prms.add_prm("b1",   2, -20.0,   20.0,  1.0);

    prms.add_prm("b2",  -2,   0.1,    0.7,  1.0);
    prms.add_prm("b2",   2, -20.0,   20.0,  1.0);

    prms.add_prm("qf4",  2, -20.0,   20.0,  1.0);
    prms.add_prm("qf4", -1,   0.1,    0.1,  1.0);

    // Parameters are initialized in optimizer.
}


void set_lat_constr_2(constr_type &constr)
{
  // Lattice constraints are: alpha_x,y, beta_x,y, eta_x, eta'_x.
  constr.add_constr(Elem_GetPos(ElemIndex("sfh"), 1),
			1e1, 1e1, 0e0, 0e0, 1e3,  0e0,
			0e0, 0e0, 0e0, 0e0, 7e-2, 0e0);
  constr.add_constr(Elem_GetPos(ElemIndex("b2"), 1),
			1e1, 1e1, 0e0, 0e0, 0e0, 1e1,
			0e0, 0e0, 0e0, 0e0, 0e0, 0e0);
  constr.add_constr(Elem_GetPos(ElemIndex("b1"), 2),
			0e0, 0e0, 0e0, 0e0, 1e5, 1e5,
			0e0, 0e0, 0e0, 0e0, 0e0, 0e0);
  constr.add_constr(globval.Cell_nLoc,
			1e1, 1e1, 1e-1, 1e-1, 0e0, 0e0,
			0e0, 0e0, 4e0,  2.5e0,  0e0, 0e0);

  lat_prms.bn_tol = 1e-6; lat_prms.step = 1.0;

  lat_constr.Fnum_b3.push_back(ElemIndex("sfh"));
  lat_constr.Fnum_b3.push_back(ElemIndex("sd1a"));
  lat_constr.Fnum_b3.push_back(ElemIndex("sd1b"));
  lat_constr.Fnum_b3.push_back(ElemIndex("sd2"));

  lat_constr.Fnum_b1.push_back(ElemIndex("b1"));
  lat_constr.Fnum_b1.push_back(ElemIndex("b2"));

  lat_constr.ini_constr(1e2, 0.190, 1e-7, 7.5);
}


int main(int argc, char *argv[])
{
  double eps_x, dnu[2];

  reverse_elem = !false;

  trace = false;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  globval.H_exact        = false;  globval.quad_fringe = false;
  globval.Cavity_on      = false;  globval.radiation   = false;
  globval.emittance      = false;  globval.IBS         = false;
  globval.pathlength     = false;  globval.bpm         = 0;
  globval.dip_edge_fudge = true;

  // set_map_reversal(ElemIndex("line_inv"));

  if (!false) {
    Ring_GetTwiss(true, 0e0); printglob();
    dnu[X_] = 0.0; dnu[Y_] = 0.1;
    set_map(ElemIndex("ps_rot"), dnu);
  }

  if (true) {
    Ring_GetTwiss(true, 0e0); printglob();
  } else
    ttwiss(ic[0], ic[1], ic[2], ic[3], 0e0);

  eps_x = get_eps_x1(true);
  printf("\neps_x = %5.3f nm.rad\n\n", eps_x);

  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);

  if (!false) {
    // set_lat_prms_1(lat_prms);
    // set_lat_constr_1(lat_constr);

    set_lat_prms_2(lat_prms);
    set_lat_constr_2(lat_constr);

    no_sxt();
    fit_powell(lat_prms, 1e-3, f_achrom);
  }
}
