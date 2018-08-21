#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


double rad2deg(const double a) { return a*180e0/M_PI; }

double deg2rad(const double a) { return a*M_PI/180e0; }


// Standard Cell.
const double ic[][2] =
  {{0.0, 0.0}, {4.97103, 5.52181}, {0.0, 0.0}, {0.0, 0.0}};


int                 n_iter, n_b1, n_b3;
double              chi2_ref, chi2_prt, chi2, eps_x, beta_b3L_2[2], nu[2];
double              phi;
std::vector<int>    Fnum_b1, Fnum_b3;
std::vector<string> drv_terms;


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
  int                 n_prm;
  double              bn_tol, step;
  std::vector<double> bn_min, bn_max, bn_scl, l0;
  std::vector<int>    Fnum, n;

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
  int                                n_loc;
  std::vector< std::vector<double> > value, scl;
  std::vector<int>                   Fnum, Fnum_b3, loc, type;

  void add_constr(const int loc,
		  const double scl1, const double scl2, const double scl3,
		  const double scl4, const double scl5, const double scl6,
		  const double v1, const double v2, const double v3,
		  const double v4, const double v5, const double v6);
  double get_chi2(void);
  void prt_Jacobian(void) const;
  void prt_constr(void) const;
};


param_type  b2_prms;
constr_type achr_constr;


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

  n_prm = Fnum.size();
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
  }
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
}


double constr_type::get_chi2(void)
{
  int    j, k;
  double chi2;

  chi2 = 0e0;
  for (k = 0; k < (int)loc.size(); k++) {
    chi2 += scl[k][0]*sqr(Cell[loc[k]].Alpha[X_]-value[k][0]);
    chi2 += scl[k][1]*sqr(Cell[loc[k]].Alpha[Y_]-value[k][1]);
    chi2 += scl[k][2]*sqr(Cell[loc[k]].Beta[X_]-value[k][2]);
    chi2 += scl[k][3]*sqr(Cell[loc[k]].Beta[Y_]-value[k][3]);
    chi2 += scl[k][4]*sqr(Cell[loc[k]].Eta[X_]-value[k][4]);
    chi2 += scl[k][5]*sqr(Cell[loc[k]].Etap[X_]-value[k][5]);
  }
  return chi2;
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


void prt_system(const int m, const int n, double **A, double *b)
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


void fit_ksi1(const double ksi_x, const double ksi_y, const double db3L)
{
  int    j, k, n;
  double **A, **U, **V, *w, *b, *x, b3L, a3L;

  const bool   prt = false;
  const int    m   = 2;
  const double
    ksi0[]  = {ksi_x, ksi_y},
    svd_cut = 1e-10;

  n = Fnum_b3.size();

  A = dmatrix(1, m, 1, n); U = dmatrix(1, m, 1, n); V = dmatrix(1, n, 1, n);
  w = dvector(1, n); b = dvector(1, m); x = dvector(1, n);

  for (k = 1; k <= n; k++) {
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

  if (prt) prt_system(m, n, A, b);

  dmcopy(A, m, n, U); dsvdcmp(U, m, n, w, V);
  if (prt) {
    printf("\nfit_ksi1:\n  singular values:\n  ");
    for (j = 1; j <= n; j++) {
      printf("%11.3e", w[j]);
      if (w[j] < svd_cut) {
	w[j] = 0e0;
	printf(" (zeroed)");
      }
    }
    printf("\n");
  }
  dsvbksb(U, w, V, m, n, b, x);

  for (k = 1; k <= n; k++)
    set_dbnL_design_fam(Fnum_b3[k-1], Sext, x[k], 0e0);

  if (prt) {
    printf("  b3:\n  ");
    for (k = 0; k < (int)Fnum_b3.size(); k++) {
      get_bn_design_elem(Fnum_b3[k], 1, Sext, b3L, a3L);
      printf(" %9.5f", b3L);
    }
    printf("\n");
  }

  free_dmatrix(A, 1, m, 1, n); free_dmatrix(U, 1, m, 1, n);
  free_dmatrix(V, 1, n, 1, n);
  free_dvector(w, 1, n); free_dvector(b, 1, m); free_dvector(x, 1, n);
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


void prt_b2(const param_type &b2_prms, const double *b2)
{
  long int loc;
  int      k;
  FILE     *outf;

  std::string file_name = "opt_des_1_b2.out";

  outf = file_write(file_name.c_str());

  for (k = 0; k < b2_prms.n_prm; k++) {
    loc = Elem_GetPos(b2_prms.Fnum[k], 1);
    prt_name(outf, Cell[loc].Elem.PName, -1);
    if (b2_prms.n[k] != -1)
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


void prt_b3(void)
{
  long int loc;
  int      k;
  FILE     *outf;

  std::string file_name = "opt_des_1_b3.out";

  outf = file_write(file_name.c_str());

  for (k = 0; k < (int)Fnum_b3.size(); k++) {
    loc = Elem_GetPos(Fnum_b3[k], 1);
    prt_name(outf, Cell[loc].Elem.PName, 0);
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


void fit_powell(param_type &b2_prms, const double eps, double (*f)(double *),
		void (*f_prt)(double *))
{
  int          n_b2, i, j, iter;
  double       *b2, **xi, fret, eps_x;
  ss_vect<tps> A;

  n_b2 = b2_prms.n_prm;

  b2 = dvector(1, n_b2); xi = dmatrix(1, n_b2, 1, n_b2);

  b2_prms.ini_prm(b2);
  f(b2);

  // Set initial directions (unit vectors).
  for (i = 1; i <= n_b2; i++)
    for (j = 1; j <= n_b2; j++)
      xi[i][j] = (i == j)? eps : 0e0;

  n_iter = 0; chi2_ref = 1e30; chi2_prt = 1e30;
  dpowell(b2, xi, n_b2, 1e-16, &iter, &fret, f);

  printf("\n  fret = %12.5e\n", fret);
  printf("b2s:\n");
  b2_prms.prt_prm(b2);
  // b2_prms.set_prm(b2);
  // eps_x = get_lin_opt(false);
  // f_prt(b2);

  free_dvector(b2, 1, n_b2); free_dmatrix(xi, 1, n_b2, 1, n_b2);
}


void prt_achrom(double *b2)
{
  int    k;
  double b3L, a3L;

  chi2_prt = chi2;
  printf("\n%3d chi2: %12.5e -> %12.5e\n", n_iter, chi2_ref, chi2_prt);
  chi2_ref = chi2;

  printf("\n  Linear Optics:\n"
	 "     %5.3f"
	 "\n    %8.5f"
	 "\n    %6.3f       %6.3f"
	 "\n    %6.3f       %6.3f"
	 "\n    %8.5f     %8.5f     %8.5f     %8.5f"
	 "\n    %8.5f     %8.5f     %8.5f     %8.5f"
	 "\n    %8.5f     %8.5f"
	 "\n    %8.5f     %8.5f"
	 "\n    %12.5e %12.5e\n",
	 eps_x, phi,
	 Cell[achr_constr.loc[3]].Beta[X_], Cell[achr_constr.loc[3]].Beta[Y_],
	 nu[X_], nu[Y_],
	 Cell[achr_constr.loc[0]].Alpha[X_], Cell[achr_constr.loc[0]].Alpha[Y_],
	 Cell[achr_constr.loc[0]].Eta[X_], Cell[achr_constr.loc[0]].Etap[X_],
	 Cell[achr_constr.loc[1]].Alpha[X_], Cell[achr_constr.loc[1]].Alpha[Y_],
	 Cell[achr_constr.loc[1]].Eta[X_], Cell[achr_constr.loc[1]].Etap[X_],
	 Cell[achr_constr.loc[2]].Eta[X_], Cell[achr_constr.loc[2]].Etap[X_],
	 Cell[achr_constr.loc[3]].Alpha[X_], Cell[achr_constr.loc[3]].Alpha[Y_],
	 sqrt(beta_b3L_2[X_]), sqrt(beta_b3L_2[Y_]));
  printf("\n  b2s:\n");
  b2_prms.prt_prm(b2);
  printf("\n  b3s:\n");
  for (k = 0; k < (int)Fnum_b3.size(); k++) {
    get_bn_design_elem(Fnum_b3[k], 1, Sext, b3L, a3L);
    if (k == 0)
      // Half sextupole.
      printf("   %9.5f", 2e0*b3L);
    else
      printf("   %9.5f", b3L);
  }
  printf("\n");

  prtmfile("flat_file.fit");
  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);

  prt_b2(b2_prms, b2);
  prt_b3();
}


double f_achrom(double *b2)
{
  bool   stable;
  int    j, k, loc;
  double b3L, a3L, phi1;

  const int n_prt = 10;
  const double
    eps0_x = 0.190, // [nm.rad].
    phi0   = 7.5;

  b2_prms.set_prm(b2);

  // Correct total bend angle.
  phi = 0e0;
  for (k = 0; k < n_b1; k++) {
    loc = Elem_GetPos(Fnum_b1[k], 1);
    phi +=
      GetnKid(Fnum_b1[k])*rad2deg(Cell[loc].Elem.PL*Cell[loc].Elem.M->Pirho);
  }
  loc = Elem_GetPos(Fnum_b1[n_b1-1], 1);
  phi1 =
    rad2deg(Cell[loc].Elem.PL*Cell[loc].Elem.M->Pirho)
    - (phi-phi0)/GetnKid(Fnum_b1[n_b1-1]);
  set_bend(Fnum_b1[n_b1-1], phi1);
 
  stable = get_nu(nu);
  if (stable) {
    fit_ksi1(0e0, 0e0, 1e-1);
    eps_x = get_lin_opt(true);

    for (k = 0; k < 2; k++)
      beta_b3L_2[k] = 0e0;
    for (k = 0; k < n_b3; k++) {
      get_bnL_design_elem(Fnum_b3[k], 1, Sext, b3L, a3L);
      for (j = 0; j < 2; j++)
	beta_b3L_2[j] +=
	  sqr(GetnKid(Fnum_b3[k])*b3L*Cell[Elem_GetPos(Fnum_b3[k], 1)].Beta[j]);
    }

    chi2 = 0e0;
    chi2 += 1e2*sqr(eps_x-eps0_x);

    chi2 += achr_constr.get_chi2();

    chi2 += 1e-7*beta_b3L_2[X_];
    chi2 += 1e-7*beta_b3L_2[Y_];

    if (chi2 < chi2_ref) {
      n_iter++;
      if (n_iter % n_prt == 0) prt_achrom(b2);
    }
  } else
    chi2 = 1e30;

  return chi2;
}


int main(int argc, char *argv[])
{
  double eps_x, dnu[2];

  reverse_elem = !false;

  trace = false;

  if (!true)
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
    dnu[X_] = 0.0; dnu[Y_] = -0.2;
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
    // b2_prms.add_prm("sd",  -1,   0.07,   0.1,  1.0);

    b2_prms.add_prm("qf1",  2, -20.0,   20.0,  1.0);
    b2_prms.add_prm("qd2",  2, -20.0,   20.0,  1.0);
    b2_prms.add_prm("d_3", -2,   0.2,    0.45, 1.0);
    b2_prms.add_prm("qf3",  2, -20.0,   20.0,  1.0);

    b2_prms.add_prm("b1",  -3, -20.0,   20.0,  1.0);
    b2_prms.add_prm("b1",  -2,   0.1,    0.8,  1.0);
    b2_prms.add_prm("b1",  -1,   0.075,  0.07, 1.0);
    b2_prms.add_prm("b1",   2, -20.0,   20.0,  1.0);

    b2_prms.add_prm("b2",  -2,   0.1,    0.5,  1.0);
    b2_prms.add_prm("b2",   2, -20.0,   20.0,  1.0);

    b2_prms.add_prm("qf4",  2, -20.0,   20.0,  1.0);
    b2_prms.add_prm("qf4", -1,   0.1,    0.1,  1.0);

    // Parameters are: alpha_x,y, beta_x,y, eta_x, eta'_x.
    achr_constr.add_constr(Elem_GetPos(ElemIndex("sfh"), 1),
			   1e1, 1e1, 0e0, 0e0, 1e0,  0e0,
			   0e0, 0e0, 0e0, 0e0, 6e-2, 0e0);
    achr_constr.add_constr(Elem_GetPos(ElemIndex("b2"), 1),
			   1e1, 1e1, 0e0, 0e0, 0e0, 1e1,
			   0e0, 0e0, 0e0, 0e0, 0e0, 0e0);
    achr_constr.add_constr(Elem_GetPos(ElemIndex("b1"), 2),
			   0e0, 0e0, 0e0, 0e0, 1e1, 1e1,
			   0e0, 0e0, 0e0, 0e0, 0e0, 0e0);
    achr_constr.add_constr(globval.Cell_nLoc,
			   1e1, 1e1, 1e-1, 1e-1, 0e0, 0e0,
			   0e0, 0e0, 5e0,  3e0,  0e0, 0e0);

   b2_prms.bn_tol = 1e-6; b2_prms.step = 1.0;

    Fnum_b3.push_back(ElemIndex("sfh"));
    Fnum_b3.push_back(ElemIndex("sd"));
    n_b3 = (int)Fnum_b3.size();

    Fnum_b1.push_back(ElemIndex("b1"));
    Fnum_b1.push_back(ElemIndex("b2"));
    n_b1 = (int)Fnum_b1.size();

    no_sxt();
    fit_powell(b2_prms, 1e-3, f_achrom, prt_achrom);
  }
}
