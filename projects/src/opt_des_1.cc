#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


const double ic[][2] =
  {{0.91959, 3.19580}, {1.27047, 17.77859}, {0.0, 0.0}, {0.0, 0.0}};


int loc[10], n, n_strength;

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
  int    i;
  double an;

  n_prm = Fnum.size();
  for (i = 1; i <= n_prm; i++) {
    if (n[i-1] > 0)
      // Multipole.
      get_bn_design_elem(Fnum[i-1], 1, n[i-1], bn[i], an);
    else if (n[i-1] == -1)
      // Drift.
      bn[i] = get_L(Fnum[i-1], 1);
    else if (n[i-1] == -2)
      // Placement.
      bn[i] = get_bn_s(Fnum[i-1], 1);
    // Bounded.
    bn[i] = bn_internal(bn[i], bn_min[i-1], bn_max[i-1]);
  }
}


void param_type::set_prm(double *bn) const
{
  int    i;
  double bn_ext;

  for (i = 1; i <= n_prm; i++) {
    // Bounded.
    bn_ext = bn_bounded(bn[i], bn_min[i-1], bn_max[i-1]);
    if (n[i-1] > 0)
	set_bn_design_fam(Fnum[i-1], n[i-1], bn_ext, 0e0);
    else if (n[i-1] == -1) {
      set_L(Fnum[i-1], bn_ext); get_S();
    } else if (n[i-1] == -2)
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
  }

  return 1470e0*sqr(globval.Energy)*I5/(I2-I4);
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
  FILE *outf;

  std::string file_name = "match.out";

  outf = file_write(file_name.c_str());

  for (k = 0; k < b2_prms.n_prm; k++) {
    loc = Elem_GetPos(b2_prms.Fnum[k], 1);
    prt_name(outf, Cell[loc].Elem.PName);
    if (b2_prms.n[k] == -1)
      fprintf(outf, " drift, l = %f;\n", Cell[loc].Elem.PL);
    else if (b2_prms.n[k] == 2)
      fprintf(outf, " quadrupole, l = %f, k = %11.8f, n = nquad"
	      ", method = meth;\n",
	      Cell[loc].Elem.PL,
	      bn_bounded(b2[k+1], b2_prms.bn_min[k], b2_prms.bn_max[k]));
  }

  fclose(outf);
}


void  prt_match(double *b2, const double eps_x)
{
  long int loc;

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
  static double chi2_ref = 1e30, chi2_prt = 1e30;

  int          loc;
  double       chi2, eps_x;
  ss_vect<tps> A;

  const int n_prt = 10;

  b2_prms.set_prm(b2);

  loc = globval.Cell_nLoc;

  globval.emittance = true;
  A = get_A(ic[0], ic[1], ic[2], ic[3]);
  Cell_Twiss(0, loc, A, false, false, 0e0);
  eps_x = get_eps_x1(false);
  globval.emittance = false;

  chi2 = 0e0;
  
  chi2 += 1e0*sqr(Cell[loc].Alpha[X_]);
  chi2 += 1e0*sqr(Cell[loc].Alpha[Y_]);
  chi2 += 1e0*sqr(Cell[loc].Etap[X_]); 

  if (chi2 < chi2_ref) {
    n++;

    if (n % n_prt == 0) {
      printf("\n%3d chi2: %12.5e -> %12.5e\n", n, chi2_prt, chi2);
      prt_match(b2, eps_x);
      prt_b2(b2_prms, b2);

      chi2_prt = min(chi2, chi2_prt);
    }

    chi2_ref = min(chi2, chi2_ref);
  }

  return chi2;
}


void fit_match(param_type &b2_prms)
{
  int          n_b2, i, j, iter;
  double       *b2, **xi, fret, eps_x;
  ss_vect<tps> A;

  n_b2 = b2_prms.n_prm;

  b2 = dvector(1, n_b2); xi = dmatrix(1, n_b2, 1, n_b2);

  b2_prms.ini_prm(b2);

  globval.emittance = true;
  A = get_A(ic[0], ic[1], ic[2], ic[3]);
  Cell_Twiss(0, globval.Cell_nLoc, A, false, false, 0e0);
  eps_x = get_eps_x1(false);
  globval.emittance = false;

  // Set initial directions (unit vectors).
  for (i = 1; i <= n_b2; i++)
    for (j = 1; j <= n_b2; j++)
      xi[i][j] = (i == j)? 1e0 : 0e0;

  n = 0;
  dpowell(b2, xi, n_b2, 1e-16, &iter, &fret, f_match);

  b2_prms.set_prm(b2);
  A = get_A(ic[0], ic[1], ic[2], ic[3]);
  Cell_Twiss(0, globval.Cell_nLoc, A, false, false, 0e0);

  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);

  free_dvector(b2, 1, n_b2); free_dmatrix(xi, 1, n_b2, 1, n_b2);
}


void  prt_emit(const double eps_x, double *b2)
{
  long int loc;

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
  static double chi2_ref = 1e30, chi2_prt = 1e30;

  int          loc;
  double       chi2, eps_x;
  ss_vect<tps> A;

  const int n_prt = 1;

  b2_prms.set_prm(b2);

  loc = globval.Cell_nLoc;

  globval.emittance = true;
  A = get_A(ic[0], ic[1], ic[2], ic[3]);
  Cell_Twiss(0, loc, A, false, false, 0e0);
  eps_x = get_eps_x1(false);
  globval.emittance = false;

  chi2 = 0e0;
  
  chi2 += 1e0*sqr(Cell[loc].Eta[X_]);
  chi2 += 1e0*sqr(Cell[loc].Etap[X_]);
  // chi2 += 1e-10*sqr(eps_x);

  if (chi2 < chi2_ref) {
    n++;

    if (n % n_prt == 0) {
      printf("\n%3d chi2: %12.5e -> %12.5e\n", n, chi2_prt, chi2);
      prt_emit(eps_x, b2);
      prt_b2(b2_prms, b2);

      chi2_prt = min(chi2, chi2_prt);
    }

    chi2_ref = min(chi2, chi2_ref);
  }

  return chi2;
}


void fit_emit(param_type &b2_prms)
{
  int          n_b2, i, j, iter;
  double       *b2, **xi, fret, eps_x;
  ss_vect<tps> A;

  n_b2 = b2_prms.n_prm;

  b2 = dvector(1, n_b2); xi = dmatrix(1, n_b2, 1, n_b2);

  b2_prms.ini_prm(b2);

  globval.emittance = true;
  A = get_A(ic[0], ic[1], ic[2], ic[3]);
  Cell_Twiss(0, globval.Cell_nLoc, A, false, false, 0e0);
  eps_x = get_eps_x1(false);
  prt_emit(eps_x, b2);
  globval.emittance = false;

  // Set initial directions (unit vectors).
  for (i = 1; i <= n_b2; i++)
    for (j = 1; j <= n_b2; j++)
      xi[i][j] = (i == j)? 1e0 : 0e0;

  n = 0;
  dpowell(b2, xi, n_b2, 1e-16, &iter, &fret, f_emit);

  b2_prms.set_prm(b2);
  A = get_A(ic[0], ic[1], ic[2], ic[3]);
  Cell_Twiss(0, globval.Cell_nLoc, A, false, false, 0e0);
 
  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);

  free_dvector(b2, 1, n_b2); free_dmatrix(xi, 1, n_b2, 1, n_b2);
}


int main(int argc, char *argv[])
{

  reverse_elem = false;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  globval.H_exact        = false;  globval.quad_fringe = false;
  globval.Cavity_on      = false;  globval.radiation   = false;
  globval.emittance      = false;  globval.IBS         = false;
  globval.pathlength     = false;  globval.bpm         = 0;
  globval.dip_edge_fudge = true;

  // Ring_GetTwiss(true, 0e0); printglob();
  ttwiss(ic[0], ic[1], ic[2], ic[3], 0e0);

  get_eps_x1(true);

  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);

  if (!false) {
    b2_prms.add_prm("q_b1",  2, -20.0, 20.0, 1.0);
    b2_prms.add_prm("qf4",   2, -20.0, 20.0, 1.0);
    b2_prms.add_prm("d_qd3",  -1, -20.0, 20.0, 1.0);

    b2_prms.bn_tol = 1e-6; b2_prms.step = 1.0;

    fit_match(b2_prms);
  }

  if (false) {
    b2_prms.add_prm("b3",  2, -20.0, 20.0, 1.0);
    b2_prms.add_prm("qf5", 2, -20.0, 20.0, 1.0);
    b2_prms.add_prm("q7",  2, -20.0, 20.0, 1.0);

    b2_prms.bn_tol = 1e-6; b2_prms.step = 1.0;

    fit_emit(b2_prms);
  }
}
