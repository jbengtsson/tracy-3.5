#define NO 1

#include "tracy_lib.h"

int no_tps = NO;

// Initial conditions: alpha, beta, eta, etap.
const double ic[][2] =
  {{1.05266, -0.25384}, {0.62733, 5.60502}, {0.06552, 0.0}, {-0.10478, 0.0}};

int loc[10], n;

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
  void ini_prm(double *bn, double *bn_lim);
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


void param_type::ini_prm(double *bn, double *bn_lim)
{
  int    i;
  double an;

  n_prm = Fnum.size();
  for (i = 1; i <= n_prm; i++) {
    bn_lim[i] = bn_max[i-1];
    if (n[i-1] > 0)
      // Multipole.
      get_bn_design_elem(Fnum[i-1], 1, n[i-1], bn[i], an);
    else if (n[i-1] == -1)
      // Drift.
      bn[i] = get_L(Fnum[i-1], 1);
    else if (n[i-1] == -2)
      // Location.
      // bn[i] = get_bn_s(-Fnum[i-1], 1, n[i-1]);
      ;
  }
}


void get_S(void);


void param_type::set_prm(double *bn) const
{
  int i, j;

  for (i = 1; i <= n_prm; i++) {
    if (n[i-1] > 0)
      for (j = 1; j <= GetnKid(Fnum[i-1]); j++)
	set_bn_design_elem(Fnum[i-1], j, n[i-1], bn[i], 0e0);
    else if (n[i-1] == -1) {
      set_L(Fnum[i-1], bn[i]); get_S();
    } else if (n[i-1] == -2)
      // set_bn_s(-Fnum[i-1], n[i-1], bn[i]);
      ;
  }
}


void param_type::prt_prm(double *bn) const
{
  int i;

  const int n_prt = 8;

  for (i = 1; i <= n_prm; i++) {
    printf(" %9.5f", bn[i]);
    if (i % n_prt == 0) printf("\n");
  }
  if (n_prm % n_prt != 0) printf("\n");
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


void prt_lin_opt(const int loc[])
{
  printf("\n      s    alpha_x  beta_x  eta_x  etap_x  alpha_y  beta_y\n");
  printf("  %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
	 Cell[loc[0]].S,
	 Cell[loc[0]].Alpha[X_], Cell[loc[0]].Beta[X_],
	 Cell[loc[0]].Eta[X_], Cell[loc[0]].Etap[X_],
	 Cell[loc[0]].Alpha[Y_], Cell[loc[0]].Beta[Y_]);
  printf("  %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
	 Cell[loc[1]].S,
	 Cell[loc[1]].Alpha[X_], Cell[loc[1]].Beta[X_],
	 Cell[loc[1]].Eta[X_], Cell[loc[1]].Etap[X_],
	 Cell[loc[1]].Alpha[Y_], Cell[loc[1]].Beta[Y_]);
  printf("  %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
	 Cell[loc[2]].S,
	 Cell[loc[2]].Alpha[X_], Cell[loc[2]].Beta[X_],
	 Cell[loc[2]].Eta[X_], Cell[loc[2]].Etap[X_],
	 Cell[loc[2]].Alpha[Y_], Cell[loc[2]].Beta[Y_]);
  printf("  %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
	 Cell[loc[3]].S,
	 Cell[loc[3]].Alpha[X_], Cell[loc[3]].Beta[X_],
	 Cell[loc[3]].Eta[X_], Cell[loc[3]].Etap[X_],
	 Cell[loc[3]].Alpha[Y_], Cell[loc[3]].Beta[Y_]);
  printf("  %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
	 Cell[loc[4]].S,
	 Cell[loc[4]].Alpha[X_], Cell[loc[4]].Beta[X_],
	 Cell[loc[4]].Eta[X_], Cell[loc[4]].Etap[X_],
	 Cell[loc[4]].Alpha[Y_], Cell[loc[4]].Beta[Y_]);
}


void prt_match(const param_type &b2_prms, const double *b2)
{
  FILE   *outf;

  std::string file_name = "match.out";

  outf = file_write(file_name.c_str());

  fprintf(outf, "l5h: drift, l = %7.5f;\n", get_L(ElemIndex("l5h"), 1));
  fprintf(outf, "l6h: drift, l = %7.5f;\n", get_L(ElemIndex("l6h"), 1));
  fprintf(outf, "l7:  drift, l = %7.5f;\n", get_L(ElemIndex("l7"), 1));
  fprintf(outf, "l8:  drift, l = %7.5f;\n", get_L(ElemIndex("l8"), 1));

  fprintf(outf, "\nbm:  bending, l = 0.14559, t = 0.5, k = %9.5f, t1 = 0.0"
	  ", t2 = 0.0,\n     gap = 0.00, N = Nbend, Method = Meth;\n", b2[1]);
  fprintf(outf, "qm:  quadrupole, l = 0.2, k = %9.5f, N = Nquad"
	  ", Method = Meth;\n", b2[4]);
  fprintf(outf, "qfe: quadrupole, l = 0.1, k = %9.5f, N = Nquad"
	  ", Method = Meth;\n", b2[2]);
  fprintf(outf, "qde: quadrupole, l = 0.1, k = %9.5f, N = Nquad"
	  ", Method = Meth;\n", b2[3]);

  fclose(outf);
}


double f_match(double *b2)
{
  static double chi2_ref = 1e30;

  int          i;
  double       chi2;
  ss_vect<tps> Ascr;

  const int n_prt = 50;

  b2_prms.set_prm(b2);

  Ascr = get_A(ic[0], ic[1], ic[2], ic[3]);
  Cell_Twiss(loc[0]+1, loc[4], Ascr, false, false, 0e0);

  chi2 = 0e0;
  // Downstream of 10 degree dipole.
  chi2 += sqr(1e5*(Cell[loc[1]].Eta[X_]));
  chi2 += sqr(1e5*(Cell[loc[1]].Etap[X_]));
  chi2 += sqr(5e1*(Cell[loc[1]].Beta[Y_]-25e0));

  // Entrance of 1st straight.
  chi2 += sqr(5e1*(Cell[loc[2]].Beta[X_]-8e0));

  // Exit of 1st straight.
  chi2 += sqr(5e1*(Cell[loc[3]].Beta[X_]-8e0));
  chi2 += sqr(1e-1*(Cell[loc[3]].Beta[Y_]-10e0));

  // Center of 2nd straight.
  chi2 += sqr(1e4*(Cell[loc[4]].Alpha[X_]));
  chi2 += sqr(1e4*(Cell[loc[4]].Alpha[Y_]));
  chi2 += sqr(1e0*(Cell[loc[4]].Beta[X_]-5e0));
  chi2 += sqr(1e0*(Cell[loc[4]].Beta[Y_]-10e0));

  for (i = 1; i <= b2_prms.n_prm; i++)
    if (fabs(b2[i]) > b2_prms.bn_max[i-1]) chi2 += 1e10;

  if (chi2 < chi2_ref) {
    n++;

    if (n % n_prt == 0) {
      printf("\n%3d chi2: %12.5e -> %12.5e\n", n, chi2_ref, chi2);
      printf("b2s: ");
      b2_prms.prt_prm(b2);

      // Downstream of 10 degree dipole.
      printf("\neta_x   = %8.5f etap_x  = %8.5f\n",
	     Cell[loc[1]].Eta[X_], Cell[loc[1]].Etap[X_]);
      printf("\nbeta_x  = %8.5f beta_y  = %8.5f\n",
	     Cell[loc[1]].Beta[X_], Cell[loc[1]].Beta[Y_]);
      // Entrance of 1st straight.
      printf("\nbeta_x  = %8.5f beta_y  = %8.5f\n",
	     Cell[loc[2]].Beta[X_], Cell[loc[2]].Beta[Y_]);
      // Exit of 1st straight.
      printf("\nbeta_x  = %8.5f beta_y  = %8.5f\n",
	     Cell[loc[3]].Beta[X_], Cell[loc[3]].Beta[Y_]);

      // Center of 2nd straight.
      printf("\nalpha_x = %8.5f alpha_y = %8.5f\n",
	     Cell[loc[4]].Alpha[X_], Cell[loc[4]].Alpha[Y_]);
      printf("beta_x  = %8.5f beta_y  = %8.5f\n",
	     Cell[loc[4]].Beta[X_], Cell[loc[4]].Beta[Y_]);

      // prt_match(b2_prms, b2);

      prtmfile("flat_file.fit");
      prt_lat(loc[0]+1, loc[4], "linlat1.out", globval.bpm, true);
      prt_lat(loc[0]+1, loc[4], "linlat.out", globval.bpm, true, 10);
    }
  }

  chi2_ref = min(chi2, chi2_ref);

  return chi2;
}


void fit_match(param_type &b2_prms)
{
  int          n_b2, i, j, iter;
  double       *b2, *b2_lim, **xi, fret;
  ss_vect<tps> Ascr;

  n_b2 = b2_prms.n_prm;

  b2 = dvector(1, n_b2); b2_lim = dvector(1, n_b2);
  xi = dmatrix(1, n_b2, 1, n_b2);

  // Upstream of 20 degree dipole.
  loc[0] = Elem_GetPos(ElemIndex("sb"),  7);
  // Downstream of 10 degree dipole.
  loc[1] = Elem_GetPos(ElemIndex("b10"), 1);
  // Entrance of 1st straight.
  loc[2] = Elem_GetPos(ElemIndex("e04"), 1);
  // Exit of 1st straight.
  loc[3] = Elem_GetPos(ElemIndex("ef2"), 8);
  // Center of 2nd straight.
  loc[4] = Elem_GetPos(ElemIndex("ef2"), 16);
  // Downstream of 20 degree dipole.
  loc[5] = Elem_GetPos(ElemIndex("b20"), 5);

  prt_lin_opt(loc);

  Ascr = get_A(ic[0], ic[1], ic[2], ic[3]);
  Cell_Twiss(loc[0]+1, loc[3], Ascr, false, false, 0e0);

  prt_lat(loc[0]+1, loc[5], "linlat1.out", globval.bpm, true);
  prt_lat(loc[0]+1, loc[5], "linlat.out", globval.bpm, true, 10);

  b2_prms.ini_prm(b2, b2_lim);

  // Set initial directions (unit vectors).
  for (i = 1; i <= n_b2; i++)
    for (j = 1; j <= n_b2; j++)
      xi[i][j] = (i == j)? 1e0 : 0e0;

  n = 0;
  dpowell(b2, xi, n_b2, 1e-16, &iter, &fret, f_match);

  Ascr = get_A(ic[0], ic[1], ic[2], ic[3]);
  Cell_Twiss(loc[0]+1, loc[3], Ascr, false, false, 0e0);

  prt_lat(loc[0]+1, loc[5], "linlat1.out", globval.bpm, true);
  prt_lat(loc[0]+1, loc[5], "linlat.out", globval.bpm, true, 10);

  free_dvector(b2, 1, n_b2);  free_dvector(b2_lim, 1, n_b2);
  free_dmatrix(xi, 1, n_b2, 1, n_b2);
}


int main(int argc, char *argv[])
{
  ss_vect<tps> Ascr;

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  no_sxt();

  Ring_GetTwiss(true, 0e0); printglob();

  // prt_lat("linlat1.out", globval.bpm, true);
  // prt_lat("linlat.out", globval.bpm, true, 10);

  if (true) {
    b2_prms.add_prm("q01",  2, 0.0, 25.0, 1.0);
    b2_prms.add_prm("q02",  2, 0.0, 25.0, 1.0);
    b2_prms.add_prm("q03",  2, 0.0, 25.0, 1.0);
    b2_prms.add_prm("eq01", 2, 0.0, 25.0, 1.0);
    b2_prms.add_prm("eq02", 2, 0.0, 25.0, 1.0);
    b2_prms.add_prm("eq03", 2, 0.0, 25.0, 1.0);
    b2_prms.add_prm("eq04", 2, 0.0, 25.0, 1.0);
    b2_prms.add_prm("eq05", 2, 0.0, 25.0, 1.0);

    b2_prms.bn_tol = 1e-6; b2_prms.step = 1.0;

    fit_match(b2_prms);
  }
}
