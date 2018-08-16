#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


double rad2deg(const double a) { return a*180e0/M_PI; }

double deg2rad(const double a) { return a*M_PI/180e0; }


// Standard Cell.
const double ic[][2] =
  {{0.0, 0.0}, {6.89632, 2.63882}, {0.0, 0.0}, {0.0, 0.0}};
// Dipole Cel.
// const double ic[][2] =
//   {{0.91959, 3.19580}, {1.27047, 17.77859}, {0.0, 0.0}, {0.0, 0.0}};


int                 n_iter;
double              chi2_ref, chi2_prt, beta_b3L_2[2], phi;
std::vector<int>    Fnum_b3, Fnum_b1;
std::vector<string> drv_terms;


double bn_internal(const double bn_bounded,
		   const double bn_min, const double bn_max);
double bn_bounded(const double bn_internal,
		  const double bn_min, const double bn_max);
double get_bn_s(const int Fnum, const int Knum);
void set_bn_s(const int Fnum, const double dbn);
void get_S(void);


struct param_type {
private:

public:
  int                 n_prm;
  double              bn_tol, step;
  std::vector<double> bn_min, bn_max, bn_scl;
  std::vector<int>    Fnum, n;

  void add_prm(const std::string Fname, const int n,
	       const double bn_min, const double bn_max,
	       const double bn_scl);
  void ini_prm(double *bn);
  void set_prm(double *bn) const;
  void prt_prm(double *bn) const;
};


param_type b2_prms;


void param_type::add_prm(const std::string Fname, const int n,
			 const double bn_min, const double bn_max,
			 const double bn_scl)
{
  Fnum.push_back(ElemIndex(Fname.c_str()));
  this->n.push_back(n);
  this->bn_min.push_back(bn_min);
  this->bn_max.push_back(bn_max);
  this->bn_scl.push_back(bn_scl);
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
    else if (n[i-1] == -1)
      // Length.
      bn[i] = get_L(Fnum[i-1], 1);
    else if (n[i-1] == -2) {
      // Bend angle; L is fixed.
      loc = Elem_GetPos(Fnum[i-1], 1);
      bn[i] = rad2deg(Cell[loc].Elem.PL*Cell[loc].Elem.M->Pirho);
    } else if (n[i-1] == -3)
      // Location.
      bn[i] = get_bn_s(Fnum[i-1], 1);
    // Bounded.
    bn[i] = bn_internal(bn[i], bn_min[i-1], bn_max[i-1]);
  }
}


void param_type::set_prm(double *bn) const
{
  int    i, k, loc;
  double bn_ext, i_rho;

  for (i = 1; i <= n_prm; i++) {
    // Bounded.
    bn_ext = bn_bounded(bn[i], bn_min[i-1], bn_max[i-1]);
    if (n[i-1] > 0)
      // Multipole strength.
      set_bn_design_fam(Fnum[i-1], n[i-1], bn_ext, 0e0);
    else if (n[i-1] == -1) {
      // Length.
      set_L(Fnum[i-1], bn_ext); get_S();
    } else if (n[i-1] == -2) {
      // Bend angle; L is fixed.
      loc = Elem_GetPos(Fnum[i-1], 1);
      i_rho = deg2rad(bn_ext)/Cell[loc].Elem.PL;
      for (k = 1; k <= GetnKid(Fnum[i-1]); k++) {
	loc = Elem_GetPos(Fnum[i-1], k);
	Cell[loc].Elem.M->Pirho = i_rho;
      }
    } else if (n[i-1] == -3)
      // Location.
      set_bn_s(Fnum[i-1], bn_ext);
  }
}


void param_type::prt_prm(double *bn) const
{
  int    i;
  double bn_ext;

  const int n_prt = 8;

  for (i = 1; i <= n_prm; i++) {
    // Bounded.
    bn_ext = bn_bounded(bn[i], bn_min[i-1], bn_max[i-1]);
    printf(" %9.5f", bn_ext);
    if (i % n_prt == 0) printf("\n");
  }
  if (n_prm % n_prt != 0) printf("\n");
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


void get_s_loc(const int Fnum, const int Knum, int loc[])
{
  partsName name;

  const bool prt = false;

  // Point to multipole.
  loc[1] = Elem_GetPos(Fnum, Knum);
  if (Cell[loc[1]-1].Elem.PName[1] == 'u') {
    loc[0] = loc[1] - 1;
    strcpy(name, Cell[loc[1]-1].Elem.PName); name[1] = 'd';
    loc[2] = Elem_GetPos(ElemIndex(name), Knum);
  } else if (Cell[loc[1]-1].Elem.PName[1] == 'd') {
    loc[2] = loc[1] - 1;
    strcpy(name, Cell[loc[1]-1].Elem.PName); name[1] = 'u';
    loc[0] = Elem_GetPos(ElemIndex(name), Knum);
  } else if (Cell[loc[1]+1].Elem.PName[1] == 'd') {
    loc[2] = loc[1] + 1;
    strcpy(name, Cell[loc[1]+1].Elem.PName); name[1] = 'u';
    loc[0] = Elem_GetPos(ElemIndex(name), Knum);
  } else if (Cell[loc[1]+1].Elem.PName[1] == 'u') {
    loc[0] = loc[1] + 1;
    strcpy(name, Cell[loc[1]+1].Elem.PName); name[1] = 'd';
    loc[2] = Elem_GetPos(ElemIndex(name), Knum);
  } else {
    printf("\nget_s_loc: configuration error %s (%d)\n",
	   Cell[loc[1]].Elem.PName, loc[1]);
    exit(1);
  }

  if (prt)
    printf("\nget_s_loc: %s %s %s\n",
	   Cell[loc[0]].Elem.PName, Cell[loc[1]].Elem.PName,
	   Cell[loc[2]].Elem.PName);
}


double get_bn_s(const int Fnum, const int Knum)
{
  int    loc[3];
  double ds;

  const bool prt = false;

  get_s_loc(Fnum, Knum, loc);
  ds = Cell[loc[0]].Elem.PL;

  if (prt)
    printf("\nget_bn_s:  %s %s(%d) %s %10.3e %10.3e\n",
	   Cell[loc[0]].Elem.PName, Cell[loc[1]].Elem.PName, Knum,
	   Cell[loc[2]].Elem.PName, Cell[loc[0]].Elem.PL, Cell[loc[2]].Elem.PL);

  return ds;
}


void set_bn_s(const int Fnum, const int Knum, const double ds)
{
  int loc[3];

  const bool prt = false;

  get_s_loc(Fnum, Knum, loc);
  set_L(Cell[loc[0]].Fnum, Knum, ds);
  set_L(Cell[loc[2]].Fnum, Knum, -ds);

  if (prt)
    printf("\nset_bn_s:  %s %s(%d) %s %10.3e %10.3e\n",
	   Cell[loc[0]].Elem.PName, Cell[loc[1]].Elem.PName, Knum,
	   Cell[loc[2]].Elem.PName, ds, -ds);
}


void set_bn_s(const int Fnum, const double ds)
{
  int k;

  for (k = 1; k <= GetnKid(Fnum); k++)
    set_bn_s(Fnum, k, ds);
}


void get_S(void)
{
  int    j;
  double S;

  S = 0e0;
  for (j = 0; j <= globval.Cell_nLoc; j++) {
    S += Cell[j].Elem.PL; Cell[j].S = S;
  }
}


double get_eps_x1(const bool track)
{
  // eps_x [pm.rad].
  ss_vect<tps> A;

  if (track) {
    globval.emittance = true;
    A = get_A(ic[0], ic[1], ic[2], ic[3]);
    Cell_Twiss(0, globval.Cell_nLoc, A, false, false, 0e0);
    globval.emittance = false;
  }

  return 1470e0*sqr(globval.Energy)*I5/(I2-I4);
}


void get_nu(ss_vect<tps> &A, const double delta, double nu[])
{
  long int     lastpos;
  ss_vect<tps> A_delta;

  A_delta = get_A_CS(2, A, nu); A_delta[delta_] += delta;
  Cell_Pass(0, globval.Cell_nLoc, A_delta, lastpos);
  get_A_CS(2, A_delta, nu);
  printf("\n[%18.16f, %18.16f]\n", nu[X_], nu[Y_]);
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
  double **A, **U, **V, *w, *b, *x;

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
    for (j = 1; j <= 2; j++)
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
    printf("\nfit_ksi1 singular values:\n");
    for (j = 1; j <= n; j++) {
      printf("%11.3e", w[j]);
      if (w[j] < svd_cut) {
	w[j] = 0e0;
	printf(" (zeroed)");
      }
      printf("\n");
    }
  }
  dsvbksb(U, w, V, m, n, b, x);

  for (k = 1; k <= n; k++)
    set_dbnL_design_fam(Fnum_b3[k-1], Sext, x[k], 0e0);

  free_dmatrix(A, 1, m, 1, n); free_dmatrix(U, 1, m, 1, n);
  free_dmatrix(V, 1, n, 1, n);
  free_dvector(w, 1, n); free_dvector(b, 1, m); free_dvector(x, 1, n);
}


void prt_name(FILE *outf, const char *name)
{
  int j, k, len;

  len = strlen(name);
  j = 0;
  do {
    fprintf(outf, "%c", name[j]);
    j++;
  } while ((j < len) && (name[j] != ' '));
  fprintf(outf, ":");
  // for (k = j; k < len; k++)
  //   fprintf(outf, "%c", name[k]);
}


void prt_b2(const param_type &b2_prms, const double *b2)
{
  long int loc;
  int      k;
  double   phi;
  FILE     *outf;

  std::string file_name = "opt_des_1.out";

  outf = file_write(file_name.c_str());

  for (k = 0; k < b2_prms.n_prm; k++) {
    loc = Elem_GetPos(b2_prms.Fnum[k], 1);
    prt_name(outf, Cell[loc].Elem.PName);
    if (b2_prms.n[k] == -2) {
      phi = rad2deg(Cell[loc].Elem.PL*Cell[loc].Elem.M->Pirho);
      fprintf(outf, " bending, l = %11.8f, t = %11.8f, t1 = %11.8f"
	      ", t2 = %11.8f,  k = %11.8f, n = nbend, method = 4;\n",
	      Cell[loc].Elem.PL, phi, Cell[loc].Elem.M->PTx1,
	      Cell[loc].Elem.M->PTx2, Cell[loc].Elem.M->PBpar[Quad+HOMmax]);
    } else if (b2_prms.n[k] == -1)
      fprintf(outf, " drift, l = %11.8f;\n", Cell[loc].Elem.PL);
    else if (b2_prms.n[k] == 2)
      fprintf(outf, " quadrupole, l = %11.8f, k = %11.8f, n = nquad"
	      ", method = 4;\n",
	      Cell[loc].Elem.PL, Cell[loc].Elem.M->PBpar[Quad+HOMmax]);
  }

  fclose(outf);
}


double get_lin_opt(const bool periodic)
{
  double       eps_x;
  ss_vect<tps> A;

  if (periodic) {
    eps_x = get_eps_x1(true);
    Ring_GetTwiss(true, 0e0);
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
		void (*f_prt)(double *, const double, const double))
{
  int          n_b2, i, j, iter;
  double       *b2, **xi, fret, eps_x;
  ss_vect<tps> A;

  n_b2 = b2_prms.n_prm;

  b2 = dvector(1, n_b2); xi = dmatrix(1, n_b2, 1, n_b2);

  b2_prms.ini_prm(b2);
  eps_x = get_lin_opt(true);

  printf("\neps_x = %5.3f nm.rad\n\n", eps_x);
  f_prt(b2, eps_x, 0e0);

  // Set initial directions (unit vectors).
  for (i = 1; i <= n_b2; i++)
    for (j = 1; j <= n_b2; j++)
      xi[i][j] = (i == j)? eps : 0e0;

  n_iter = 0; chi2_ref = 1e30; chi2_prt = 1e30;
  dpowell(b2, xi, n_b2, 1e-16, &iter, &fret, f);

  b2_prms.set_prm(b2);
  eps_x = get_lin_opt(true);

  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);

  free_dvector(b2, 1, n_b2); free_dmatrix(xi, 1, n_b2, 1, n_b2);
}


void prt_achrom(double *b2, const double eps_x, const double chi2)
{
  int    loc1, loc2, k;
  double b3L, a3L;

  chi2_prt = chi2;
  printf("\n%3d chi2: %12.5e -> %12.5e\n", n_iter, chi2_ref, chi2_prt);
  chi2_ref = chi2;

  loc1 = Elem_GetPos(ElemIndex("b1"), 2);
  loc2 = globval.Cell_nLoc;

  printf("Linear Optics:\n %5.3f %6.3f %6.3f %6.3f %8.5f %8.5f %8.5f %8.5f\n",
	 eps_x, 2e0*phi, Cell[loc2].Beta[X_], Cell[loc2].Beta[Y_],
	 Cell[loc1].Eta[X_], Cell[loc1].Etap[X_],
	 sqrt(beta_b3L_2[X_]), sqrt(beta_b3L_2[Y_]));
  printf("b2s:\n");
  b2_prms.prt_prm(b2);
  printf("b3s:\n");
  for (k = 0; k < (int)Fnum_b3.size(); k++) {
    get_bn_design_elem(Fnum_b3[k], 1, Sext, b3L, a3L);
    if (k == 0)
      // Half sextupole.
      printf(" %9.5f", 1e0*b3L);
    else
      printf(" %9.5f", b3L);
  }
  printf("\n");


  prtmfile("flat_file.fit");
  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);

  prt_b2(b2_prms, b2);
}


double f_achrom(double *b2)
{
  int    n_b3, n_b1, j, k, loc;
  double eps_x, chi2, b3L, a3L;

  const int    n_prt = 1;
  const double
    eps0_x = 0.150,    // [nm.rad].
    phi0   = 15.0/2e0; // Cell bend angle [deg].

  n_b3 = (int)Fnum_b3.size();
  n_b1 = (int)Fnum_b1.size();

  b2_prms.set_prm(b2);
  fit_ksi1(0e0, 0e0, 1e-3);
  eps_x = get_lin_opt(true);

  for (k = 0; k < 2; k++)
    beta_b3L_2[k] = 0e0;
  for (k = 0; k < n_b3; k++) {
    get_bnL_design_elem(Fnum_b3[k], 1, Sext, b3L, a3L);
    loc =  Elem_GetPos(Fnum_b3[k], 1);
    for (j = 0; j < 2; j++)
      beta_b3L_2[j] += sqr(b3L*Cell[loc].Beta[j]);
  }

  phi = 0e0;
  for (k = 0; k < n_b1; k++) {
    loc = Elem_GetPos(Fnum_b1[k], 1);
    phi += rad2deg(Cell[loc].Elem.PL*Cell[loc].Elem.M->Pirho);
  }

  chi2 = 0e0;
  if (globval.stable) {
    loc = Elem_GetPos(ElemIndex("b1"), 2);

    chi2 += 1e0*sqr(eps_x-eps0_x);
    chi2 += 1e0*sqr(phi-phi0);
    chi2 += 1e1*sqr(Cell[loc].Eta[X_]);
    chi2 += 1e1*sqr(Cell[loc].Etap[X_]);
    chi2 += 1e-8*beta_b3L_2[X_];
    chi2 += 1e-8*beta_b3L_2[Y_];
 
    if (chi2 < chi2_ref) {
      n_iter++;
      if (n_iter % n_prt == 0) prt_achrom(b2, eps_x, chi2);
    }
  } else
    chi2 = 1e30;

  return chi2;
}


void prt_match(double *b2, const double eps_x, const double chi2)
{
  long int loc;

  chi2_prt = chi2;
  printf("\n%3d chi2: %12.5e -> %12.5e\n", n_iter, chi2_ref, chi2_prt);
  chi2_ref = chi2;

  loc = globval.Cell_nLoc;

  printf("Linear Optics:\n %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
	 eps_x, Cell[loc].Alpha[X_], Cell[loc].Alpha[Y_],
	 Cell[loc].Beta[X_], Cell[loc].Beta[Y_],
	 Cell[loc].Etap[X_]);
  printf("b2s:\n");
  b2_prms.prt_prm(b2);

  prtmfile("flat_file.fit");
  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);
}


double f_match(double *b2)
{
  int    loc;
  double eps_x, chi2;

  const int n_prt = 10;

  b2_prms.set_prm(b2);
  eps_x = get_lin_opt(true);

  loc = globval.Cell_nLoc;

  chi2 = 0e0;
  chi2 += 1e0*sqr(Cell[loc].Alpha[X_]);
  chi2 += 1e0*sqr(Cell[loc].Alpha[Y_]);
  chi2 += 1e0*sqr(Cell[loc].Etap[X_]); 

  if (chi2 < chi2_ref) {
    n_iter++;
    if (n_iter % n_prt == 0) prt_match(b2, eps_x, chi2);
  }

  return chi2;
}


void prt_emit(double *b2, const double eps_x, const double chi2)
{
  long int loc;

  chi2_prt = chi2;
  printf("\n%3d chi2: %12.5e -> %12.5e\n", n_iter, chi2_ref, chi2_prt);
  chi2_ref = chi2;

  loc = globval.Cell_nLoc;

  printf("Linear Optics:\n %12.5e %12.5e %12.5e\n",
	 eps_x, Cell[loc].Eta[X_], Cell[loc].Etap[X_]);
  printf("b2s:\n");
  b2_prms.prt_prm(b2);

  prtmfile("flat_file.fit");
  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);
}


double f_emit(double *b2)
{
  int    loc;
  double eps_x, chi2;

  const int n_prt = 1;

  b2_prms.set_prm(b2);
  eps_x = get_lin_opt(true);

  loc = globval.Cell_nLoc;

  chi2 = 0e0;
  
  chi2 += 1e0*sqr(Cell[loc].Eta[X_]);
  chi2 += 1e0*sqr(Cell[loc].Etap[X_]);
  // chi2 += 1e-10*sqr(eps_x);

  if (chi2 < chi2_ref) {
    n_iter++;
    if (n_iter % n_prt == 0) prt_emit(b2, eps_x, chi2);
  }

  return chi2;
}


void prt_ksi1(double *b2, const double eps_x, const double chi2)
{
  long int loc;

  chi2_prt = chi2;
  printf("\n%3d chi2: %12.5e -> %12.5e\n", n_iter, chi2_ref, chi2_prt);
  chi2_ref = chi2;

  loc = Elem_GetPos(ElemIndex("ms"), 1);

  printf("Linear Optics:\n %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
	 eps_x, globval.Chrom[X_], globval.Chrom[Y_],
	 Cell[loc].Eta[X_], Cell[loc].Etap[X_],
	 Cell[loc].Alpha[X_], Cell[loc].Alpha[Y_]);
  loc = Elem_GetPos(ElemIndex("b1_5"), 2);
  printf(" %12.5e %12.5e\n", Cell[loc].Eta[X_], Cell[loc].Etap[X_]);
  printf("b2s:\n");
  b2_prms.prt_prm(b2);

  prtmfile("flat_file.fit");
  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);
}


double f_ksi1(double *b2)
{
  int    loc;
  double eps_x, chi2;

  const int    n_prt = 10;
  const double
    eta_x_ref  = 2.44e-2,
    eps_x_ref  = 0.188, // nm.rad.
    ksi1_ref[] = {-2.9, -3.0};

  b2_prms.set_prm(b2);
  eps_x = get_lin_opt(true);

  chi2 = 0e0;
  if (globval.stable) {
    loc = Elem_GetPos(ElemIndex("ms"), 1);
    chi2 += 1e0*sqr(Cell[loc].Eta[X_]-eta_x_ref);
    chi2 += 1e4*sqr(Cell[loc].Etap[X_]);
    chi2 += 1e4*sqr(Cell[loc].Alpha[X_]);
    chi2 += 1e4*sqr(Cell[loc].Alpha[Y_]);

    loc = Elem_GetPos(ElemIndex("b1_5"), 2);
    chi2 += 1e4*sqr(Cell[loc].Etap[X_]);
    chi2 += 1e4*sqr(Cell[loc].Etap[X_]);

    chi2 += 1e1*sqr(eps_x-eps_x_ref);

    chi2 += 1e-1*sqr(globval.Chrom[X_]-ksi1_ref[X_]);
    chi2 += 1e0*sqr(globval.Chrom[Y_]-ksi1_ref[Y_]);

    if (chi2 < chi2_ref) {
      n_iter++;
      if (n_iter % n_prt == 0) prt_ksi1(b2, eps_x, chi2);
    }
  } else
    chi2 = 1e30;

  return chi2;
}


int main(int argc, char *argv[])
{
  double       eps_x;
  ss_vect<tps> A;

  reverse_elem = !false;

  if (!true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  globval.H_exact        = false;  globval.quad_fringe = false;
  globval.Cavity_on      = false;  globval.radiation   = false;
  globval.emittance      = false;  globval.IBS         = false;
  globval.pathlength     = false;  globval.bpm         = 0;
  globval.dip_edge_fudge = true;

  if (true) {
    Ring_GetTwiss(true, 0e0); printglob();
  } else
    ttwiss(ic[0], ic[1], ic[2], ic[3], 0e0);

  eps_x = get_eps_x1(true);
  printf("\neps_x = %5.3f nm.rad\n\n", eps_x);

  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);

  if (!false) {
    b2_prms.add_prm("b1", -2, -20.0, 20.0, 1.0);
    b2_prms.add_prm("b2", -2, -20.0, 20.0, 1.0);
    b2_prms.add_prm("b3", -2, -20.0, 20.0, 1.0);

    b2_prms.add_prm("b1", 2, -20.0, 20.0, 1.0);
    b2_prms.add_prm("b2", 2, -20.0, 20.0, 1.0);
    b2_prms.add_prm("b3", 2, -20.0, 20.0, 1.0);

    b2_prms.add_prm("qf1", 2, -20.0, 20.0, 1.0);
    b2_prms.add_prm("qd2", 2, -20.0, 20.0, 1.0);
    b2_prms.add_prm("qd3", 2, -20.0, 20.0, 1.0);
    b2_prms.add_prm("qf4", 2, -20.0, 20.0, 1.0);
    b2_prms.add_prm("qd5", 2, -20.0, 20.0, 1.0);
    b2_prms.add_prm("qf6", 2, -20.0, 20.0, 1.0);

    b2_prms.bn_tol = 1e-6; b2_prms.step = 1.0;

    no_sxt();

    Fnum_b3.push_back(ElemIndex("sfh"));
    Fnum_b3.push_back(ElemIndex("sda"));
    Fnum_b3.push_back(ElemIndex("sdb"));

    Fnum_b1.push_back(ElemIndex("b1"));
    Fnum_b1.push_back(ElemIndex("b2"));
    Fnum_b1.push_back(ElemIndex("b3"));
    Fnum_b1.push_back(ElemIndex("b4_1"));
    Fnum_b1.push_back(ElemIndex("b4_2"));

    fit_powell(b2_prms, 1e-3, f_achrom, prt_achrom);
  }

  if (false) {
    b2_prms.add_prm("q_b1",   2, -20.0, 20.0, 1.0);
    b2_prms.add_prm("qf",     2, -20.0, 20.0, 1.0);
    b2_prms.add_prm("d_qd3", -1, -20.0, 20.0, 1.0);

    b2_prms.bn_tol = 1e-6; b2_prms.step = 1.0;

    fit_powell(b2_prms, 1e-3, f_match, prt_match);
  }

  if (false) {
    b2_prms.add_prm("b3",  2, -20.0, 20.0, 1.0);
    b2_prms.add_prm("qf5", 2, -20.0, 20.0, 1.0);
    b2_prms.add_prm("q7",  2, -20.0, 20.0, 1.0);

    b2_prms.bn_tol = 1e-6; b2_prms.step = 1.0;

    fit_powell(b2_prms, 1e-3, f_emit, prt_emit);
  }

  if (false) {
    b2_prms.add_prm("qf1ss", 2, -20.0, 20.0, 1.0);
    b2_prms.add_prm("qd2ss", 2, -20.0, 20.0, 1.0);

    b2_prms.add_prm("qd3l",  2, -20.0, 20.0, 1.0);
    b2_prms.add_prm("qf4l",  2, -20.0, 20.0, 1.0);
    b2_prms.add_prm("qd3r",  2, -20.0, 20.0, 1.0);
    b2_prms.add_prm("qf4r",  2, -20.0, 20.0, 1.0);
    b2_prms.add_prm("qf5",   2, -20.0, 20.0, 1.0);
    b2_prms.add_prm("qf6",   2, -20.0, 20.0, 1.0);

    b2_prms.bn_tol = 1e-6; b2_prms.step = 1.0;

    b2_prms.add_prm("q_b1",  2, -20.0, 20.0, 1.0);
    b2_prms.add_prm("q_b2",  2, -20.0, 20.0, 1.0);

    fit_powell(b2_prms, 1e-3, f_ksi1, prt_ksi1);
  }
}
