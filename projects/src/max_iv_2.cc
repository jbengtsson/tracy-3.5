#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


const double eps_x       = 15e-3,
             nu_uc[]     = {4.0/15.0, 3.0/15.0}, // Cell tune.
             nu_sc[]     = {1.410, 0.628},       // Super period tune.
             L_uc        = 1.25,                 // Unit cell length.
             L_ss        = 10.5,                 /* Super period length;
                                                    with one unit cell. */
             eta_cuc[]   = {0.00920897, 0.0},    // Center of unit cell.
             etap_cuc[]  = {0.0, 0.0},           // Center of unit cell.
             alpha_cuc[] = {0.0, 0.0},           // Center of unit cell.
             beta_cuc[]  = {1.59916, 0.57961},
             beta_cs[]   = {3.0, 3.0};           // Center of straight.

const string sf_sd_name[] = {"sfh", "sd"};


struct param_type {
private:

public:
  int                 n_prm;
  double              bn_tol, svd_cut, step;
  std::vector<double> bn_min, bn_max, bn_scl;
  std::vector<int>    Fnum, n;

  void add_prm(const std::string Fname, const int n,
	       const double bn_min, const double bn_max,
	       const double bn_scl);
  void ini_prm(double *bn, double *bn_lim);
  void set_prm(double *bn) const;
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


void get_S(void)
{
  int    j;
  double S;

  S = 0e0;
  for (j = 0; j <= globval.Cell_nLoc; j++) {
    S += Cell[j].Elem.PL; Cell[j].S = S;
  }
}


double get_eps_x1(void)
{
  // Evaluate emittance and damping partition from the synchrotron integrals.
  bool         cav, emit;
  long int     lastpos;
  double       eps_x;
  ss_vect<tps> A;

  const bool prt = false;

  cav = globval.Cavity_on; emit = globval.emittance;

  globval.Cavity_on = false; globval.emittance = false;

  Ring_GetTwiss(false, 0e0);

  putlinmat(6, globval.Ascr, A);

  // prt_lin_map(3, A);

  globval.emittance = true;

  Cell_Pass(0, globval.Cell_nLoc, A, lastpos);

  eps_x = 1470e0*pow(globval.Energy, 2)*I5/(I2-I4);

  if (prt)
    printf("eps_x = %5.3f pm.rad, J_x = %5.3f, J_z = %5.3f \n",
	   1e3*eps_x, 1e0-I4/I2, 2e0+I4/I2);

  globval.Cavity_on = cav; globval.emittance = emit;

  return eps_x;
}


void quad_scan(const int n,
	       const char *qf, const double dk_qf,
	       const char *qd, const double dk_qd)
{
  // Parametric scan of gradients for unit cell.
  int    i, j, qf_num, qd_num;
  double k_qf, k_qd, k_qf_step, k_qd_step, a2, eps_x;
  FILE   *outf;

  const char file_name[] = "quad_scan.out";

  outf = fopen(file_name, "w");

  qf_num = ElemIndex(qf); qd_num = ElemIndex(qd);
  get_bn_design_elem(qf_num, 1, Quad, k_qf, a2);
  get_bn_design_elem(qd_num, 1, Quad, k_qd, a2);
  k_qf_step = dk_qf/n; k_qf -= dk_qf;
  k_qd_step = dk_qd/n; k_qd += dk_qd;
  for (i = -n; i <= n; i++) {
    set_bn_design_fam(qf_num, Quad, k_qf, a2);
    k_qd -= (2*n+1)*k_qd_step;
    for (j = -n; j <= n; j++) {
      set_bn_design_fam(qd_num, Quad, k_qd, a2);
      Ring_GetTwiss(true, 0e0);
      if (globval.stable)
	eps_x = get_eps_x1();
      else {
	globval.TotalTune[X_] = NAN; globval.TotalTune[Y_] = NAN;
	globval.Chrom[X_] = NAN; globval.Chrom[Y_] = NAN; eps_x = NAN;
      }
      fprintf(outf, "%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n",
	      k_qf, k_qd,
	      globval.TotalTune[X_], globval.TotalTune[Y_],
	      globval.Chrom[X_], globval.Chrom[Y_], 1e3*eps_x);
      k_qd += k_qd_step;
    }
    fprintf(outf, "\n");
    k_qf += k_qf_step;
  }

  fclose(outf);
}


void prt_emit(const param_type &b2_prms, const double *b2)
{
  int    sf_sd[2];
  double b3[2], a3;
  FILE   *outf;

  std::string file_name = "emit.out";

  sf_sd[0] = ElemIndex(sf_sd_name[0].c_str());
  sf_sd[1] = ElemIndex(sf_sd_name[1].c_str());
  get_bn_design_elem(sf_sd[0], 1, Sext, b3[0], a3);
  get_bn_design_elem(sf_sd[1], 1, Sext, b3[1], a3);

  outf = file_write(file_name.c_str());

  fprintf(outf, "l1:  drift, l = %7.5f;\n", b2[3]);
  fprintf(outf, "l2:  drift, l = %7.5f;\n", b2[4]);

  fprintf(outf, "\nbh:  bending, l = %8.5f, t = 0.5, k = %8.5f, t1 = 0.0"
	  ", t2 = 0.0,\n     gap = 0.0, N = Nbend, Method = Meth;\n",
	  b2[5], b2[1]);
  fprintf(outf, "qf:  quadrupole, l = 0.08, k = %8.5f, N = Nquad"
	  ", Method = Meth;\n\n", b2[2]);

  fprintf(outf, "sfh: sextupole, l = 0.05, k = %8.5f, N = Nsext"
	  ", Method = Meth;\n", b3[0]);
  fprintf(outf, "sd:  sextupole, l = 0.1, k = %8.5f, N = Nsext"
	  ", Method = Meth;\n", b3[1]);

  fclose(outf);
}


double f_emit(double *b2)
{
  // Optimize unit cell.
  // Lattice: unit cell.

  static double chi2_ref = 1e30;

  int    i, sf_sd[2];
  double eps1_x, tr[2], b3L[2], a3L, chi2;

  const double scl_eps = 1e2, scl_L = 1e0, scl_ksi = 1e-5;

  b2_prms.set_prm(b2);

  Ring_GetTwiss(false, 0e0); eps1_x = get_eps_x1();

  tr[X_] = globval.OneTurnMat[x_][x_] + globval.OneTurnMat[px_][px_];
  tr[Y_] = globval.OneTurnMat[y_][y_] + globval.OneTurnMat[py_][py_];
  // printf("trace: %6.3f %6.3f\n", tr[X_], tr[Y_]);

  sf_sd[0] = ElemIndex(sf_sd_name[0].c_str());
  sf_sd[1] = ElemIndex(sf_sd_name[1].c_str());
  FitChrom(sf_sd[0], sf_sd[1], 0e0, 0e0);
  get_bnL_design_elem(sf_sd[0], 1, Sext, b3L[0], a3L);
  get_bnL_design_elem(sf_sd[1], 1, Sext, b3L[1], a3L);

  chi2 = 0e0;
  chi2 += sqr(scl_eps*(eps1_x-eps_x));
  chi2 += sqr(globval.TotalTune[X_]-nu_uc[X_]);
  chi2 += sqr(globval.TotalTune[Y_]-nu_uc[Y_]);
  chi2 += sqr(scl_L*(Cell[globval.Cell_nLoc].S-L_uc));
  chi2 += sqr(scl_ksi*b3L[0]);
  chi2 += sqr(scl_ksi*b3L[1]);

  if ((fabs(tr[X_]) > 2e0) || (fabs(tr[Y_]) > 2e0)) chi2 += 1e10;
  for (i = 1; i <= b2_prms.n_prm; i++) {
    if ((b2_prms.n[i-1] == -1) && (b2[i] < b2_prms.bn_min[i-1]))
      chi2 += 1e10;
    if (fabs(b2[i]) > b2_prms.bn_max[i-1]) chi2 += 1e10;
  }

  if (chi2 < chi2_ref) {
    printf("\nchi2: %12.5e, %12.5e\n", chi2, chi2_ref);
    printf("b:    %10.3e %10.5f %10.5f %10.5f %10.5f %10.5f\n",
	   eps1_x,
	   globval.TotalTune[X_], globval.TotalTune[Y_],
	   Cell[globval.Cell_nLoc].S,
	   b3L[0], b3L[1]);
	   // globval.Chrom[X_], globval.Chrom[Y_]);
    printf("b2s: ");
    for (i = 1; i <= b2_prms.n_prm; i++)
      printf("%9.5f", b2[i]);
    printf("\n");

    prtmfile("flat_file.fit");
    prt_emit(b2_prms, b2);

    get_S();
    prt_lat("linlat1.out", globval.bpm, true);
    prt_lat("linlat.out", globval.bpm, true, 10);
  }

  chi2_ref = min(chi2, chi2_ref);

  return chi2;
}


void fit_emit(param_type &b2_prms)
{
  // Optimize unit cell.
  // Lattice: unit cell.

  int    n_b2, i, j, iter;
  double *b2, *b2_lim, **xi, fret;

  n_b2 = b2_prms.n_prm;

  b2 = dvector(1, n_b2); b2_lim = dvector(1, n_b2);
  xi = dmatrix(1, n_b2, 1, n_b2);

  b2_prms.ini_prm(b2, b2_lim);

  // Set initial directions (unit vectors).
  for (i = 1; i <= n_b2; i++)
    for (j = 1; j <= n_b2; j++)
      xi[i][j] = (i == j)? 1e0 : 0e0;

  dpowell(b2, xi, n_b2, 1e-16, &iter, &fret, f_emit);

  free_dvector(b2, 1, n_b2);  free_dvector(b2_lim, 1, n_b2);
  free_dmatrix(xi, 1, n_b2, 1, n_b2);
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


double f_hcell(double *b2)
{
  // Match half super period.
  // Lattice: half super period.

  static double chi2_ref = 1e30;

  int          i, loc1, loc2;
  double       tr[2], chi2;
  ss_vect<tps> A;

  const double scl_eta = 1e4, scl_alpha = 1e1, scl_beta = 1e0, scl_L = 0e1;

  b2_prms.set_prm(b2);

  // End of bm.
  loc1 = Elem_GetPos(ElemIndex("bm"), 1);
  // Center of straight.
  loc2 = globval.Cell_nLoc;

  A = get_A(alpha_cuc, beta_cuc, eta_cuc, etap_cuc);
  Cell_Twiss(0, globval.Cell_nLoc, A, false, false, 0e0);

  chi2 = 0e0;
  chi2 += sqr(scl_eta*Cell[loc1].Eta[X_]);
  chi2 += sqr(scl_eta*Cell[loc1].Etap[X_]);
  chi2 += sqr(scl_alpha*Cell[loc2].Alpha[X_]);
  chi2 += sqr(scl_alpha*Cell[loc2].Alpha[Y_]);
  chi2 += sqr(scl_beta*(Cell[loc2].Beta[X_]-beta_cs[X_]));
  chi2 += sqr(scl_beta*(Cell[loc2].Beta[Y_]-beta_cs[Y_]));
  chi2 += sqr(scl_L*(Cell[globval.Cell_nLoc].S-L_ss));

  if ((fabs(tr[X_]) > 2e0) || (fabs(tr[Y_]) > 2e0)) chi2 += 1e10;
  for (i = 1; i <= b2_prms.n_prm; i++) {
    if ((b2_prms.n[i-1] == -1) && (b2[i] < b2_prms.bn_min[i-1]))
      chi2 += 1e10;
    if (fabs(b2[i]) > b2_prms.bn_max[i-1]) chi2 += 1e10;
  }

  if (chi2 < chi2_ref) {
    printf("\nchi2: %12.5e, %12.5e\n", chi2, chi2_ref);
    printf("b: %10.3e %10.3e %10.3e %10.3e %10.5f %10.5f %10.5f\n",
	   Cell[loc1].Eta[X_],  Cell[loc1].Etap[X_],
	   Cell[loc2].Alpha[X_], Cell[loc2].Alpha[Y_],
	   Cell[loc2].Beta[X_],  Cell[loc2].Beta[Y_],
	   Cell[globval.Cell_nLoc].S);
    printf("b2s: ");
    for (i = 1; i <= b2_prms.n_prm; i++)
      printf(" %9.5f", b2[i]);
    printf("\n");

    prtmfile("flat_file.fit");
    prt_match(b2_prms, b2);

    prt_lat("linlat1.out", globval.bpm, true);
    prt_lat("linlat.out", globval.bpm, true, 10);
  }

  chi2_ref = min(chi2, chi2_ref);

  return chi2;
}


void fit_hcell(param_type &b2_prms)
{
  // Match half super period.
  // Lattice: super period.

  int    n_b2, i, j, iter;
  double *b2, *b2_lim, **xi, fret;

  n_b2 = b2_prms.n_prm;

  b2 = dvector(1, n_b2); b2_lim = dvector(1, n_b2);
  xi = dmatrix(1, n_b2, 1, n_b2);

  b2_prms.ini_prm(b2, b2_lim);

  // Set initial directions (unit vectors).
  for (i = 1; i <= n_b2; i++)
    for (j = 1; j <= n_b2; j++)
      xi[i][j] = (i == j)? 1e0 : 0e0;

  dpowell(b2, xi, n_b2, 1e-16, &iter, &fret, f_hcell);

  free_dvector(b2, 1, n_b2);  free_dvector(b2_lim, 1, n_b2);
  free_dmatrix(xi, 1, n_b2, 1, n_b2);
}


double f_match(double *b2)
{
  // Optimize super period.
  // Lattice: super period.

  static double chi2_ref = 1e30;

  int          i, loc1, loc2, loc3;
  double       tr[2], chi2;

  const double scl_eta = 1e4, scl_alpha = 1e1,  scl_beta = 1e1,
               scl_L   = 1e1, scl_ksi   = 5e-1, scl_nu   = 1e3;

  b2_prms.set_prm(b2);

  // Center of unit cell.
  loc1 = Elem_GetPos(ElemIndex("sfh"), 1);
  // End of 2nd bm.
  loc2 = Elem_GetPos(ElemIndex("bm"), 2);
  // Center of straight.
  loc3 = globval.Cell_nLoc;

  Ring_GetTwiss(true, 0e0);

  tr[X_] = globval.OneTurnMat[x_][x_] + globval.OneTurnMat[px_][px_];
  tr[Y_] = globval.OneTurnMat[y_][y_] + globval.OneTurnMat[py_][py_];
  // printf("trace: %6.3f %6.3f\n", tr[X_], tr[Y_]);

  chi2 = 0e0;
  chi2 += sqr(scl_eta*(Cell[loc1].Eta[X_]-eta_cuc[X_]));
  chi2 += sqr(scl_eta*Cell[loc1].Etap[X_]);
  chi2 += sqr(scl_alpha*(Cell[loc1].Alpha[X_]-alpha_cuc[X_]));
  chi2 += sqr(scl_alpha*(Cell[loc1].Alpha[Y_]-alpha_cuc[Y_]));
  chi2 += sqr(scl_beta*(Cell[loc1].Beta[X_]-beta_cuc[X_]));
  chi2 += sqr(scl_beta*(Cell[loc1].Beta[Y_]-beta_cuc[Y_]));
  chi2 += sqr(scl_eta*Cell[loc2].Eta[X_]);
  chi2 += sqr(scl_eta*Cell[loc2].Etap[X_]);
  chi2 += sqr(scl_beta*(Cell[loc3].Beta[X_]-beta_cs[X_]));
  chi2 += sqr(scl_beta*(Cell[loc3].Beta[Y_]-beta_cs[Y_]));
  chi2 += sqr(scl_L*(Cell[globval.Cell_nLoc].S-L_ss));
  chi2 += sqr(scl_ksi*globval.Chrom[X_]);
  chi2 += sqr(scl_ksi*globval.Chrom[Y_]);
  chi2 += sqr(scl_nu*(globval.TotalTune[X_]-nu_sc[X_]));
  chi2 += sqr(scl_nu*(globval.TotalTune[Y_]-nu_sc[Y_]));

  if ((fabs(tr[X_]) > 2e0) || (fabs(tr[Y_]) > 2e0)) chi2 += 1e10;
  for (i = 1; i <= b2_prms.n_prm; i++) {
    if ((b2_prms.n[i-1] == -1) && (b2[i] < b2_prms.bn_min[i-1]))
      chi2 += 1e10;
    if (fabs(b2[i]) > b2_prms.bn_max[i-1]) chi2 += 1e10;
  }

  if (chi2 < chi2_ref) {
    printf("\nchi2: %12.5e, %12.5e\n", chi2, chi2_ref);
    printf("b: %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e"
	   " %10.5f %10.5f\n %10.5f %10.5f %10.5f %10.5f %10.5f\n",
	   Cell[loc1].Etap[X_],  Cell[loc1].Etap[X_],
	   Cell[loc1].Alpha[X_], Cell[loc1].Alpha[Y_],
	   Cell[loc1].Beta[X_]-beta_cuc[X_], Cell[loc1].Beta[Y_]-beta_cuc[Y_],
	   Cell[loc2].Eta[X_], Cell[loc2].Etap[X_],
	   Cell[loc3].Beta[X_], Cell[loc3].Beta[Y_],
	   Cell[globval.Cell_nLoc].S,
	   globval.Chrom[X_], globval.Chrom[Y_],
	   globval.TotalTune[X_], globval.TotalTune[Y_]);
    printf("b2s: ");
    for (i = 1; i <= b2_prms.n_prm; i++)
      printf(" %9.5f", b2[i]);
    printf("\n");

    prtmfile("flat_file.fit");
    prt_match(b2_prms, b2);

    prt_lat("linlat1.out", globval.bpm, true);
    prt_lat("linlat.out", globval.bpm, true, 10);
  }

  chi2_ref = min(chi2, chi2_ref);

  return chi2;
}


void fit_match(param_type &b2_prms)
{
  // Minimize Optimize super period.
  // Lattice: super period.

  int    n_b2, i, j, iter;
  double *b2, *b2_lim, **xi, fret;

  n_b2 = b2_prms.n_prm;

  b2 = dvector(1, n_b2); b2_lim = dvector(1, n_b2);
  xi = dmatrix(1, n_b2, 1, n_b2);

  b2_prms.ini_prm(b2, b2_lim);

  // Set initial directions (unit vectors).
  for (i = 1; i <= n_b2; i++)
    for (j = 1; j <= n_b2; j++)
      xi[i][j] = (i == j)? 1e0 : 0e0;

  dpowell(b2, xi, n_b2, 1e-16, &iter, &fret, f_match);

  free_dvector(b2, 1, n_b2);  free_dvector(b2_lim, 1, n_b2);
  free_dmatrix(xi, 1, n_b2, 1, n_b2);
}


int main(int argc, char *argv[])
{

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  if (false)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  Ring_GetTwiss(true, 0e0); printglob();

  if (false) {
    get_eps_x1();
    GetEmittance(ElemIndex("cav"), true);
  }

  if (false) quad_scan(10, "qf", 3e0, "bh", 2e0);

  if (false) {
    b2_prms.add_prm("bh",  2, 0.0, 25.0,  1.0);
    b2_prms.add_prm("qf",  2, 0.0, 25.0,  1.0);
    b2_prms.add_prm("l1", -1, 0.05, 0.35, 0.01);
    b2_prms.add_prm("l2", -1, 0.05, 0.35, 0.01);
    b2_prms.add_prm("bh", -1, 0.05, 0.3,  0.01);

    b2_prms.bn_tol = 1e-6; b2_prms.svd_cut = 1e-8; b2_prms.step = 1.0;

    no_sxt();
    fit_emit(b2_prms);
  }

  if (false) {
    b2_prms.add_prm("bm",   2, 0.0, 25.0, 1.0);
    b2_prms.add_prm("qfe",  2, 0.0, 25.0, 1.0);
    b2_prms.add_prm("qde",  2, 0.0, 25.0, 1.0);
    b2_prms.add_prm("qm",   2, 0.0, 25.0, 1.0);

    b2_prms.add_prm("l5h", -1, 0.05, 1.0,  0.01);
    b2_prms.add_prm("l6h", -1, 0.05, 1.0,  0.01);
    b2_prms.add_prm("l7",  -1, 0.05, 1.0,  0.01);
    b2_prms.add_prm("l8",  -1, 0.05, 1.0,  0.01);

    b2_prms.bn_tol = 1e-6; b2_prms.svd_cut = 1e-8; b2_prms.step = 1.0;

    no_sxt();
    fit_hcell(b2_prms);
  }

  if (true) {
    b2_prms.add_prm("bm",   2, 0.0, 25.0, 1.0);
    b2_prms.add_prm("qm",   2, 0.0, 25.0, 1.0);
    b2_prms.add_prm("qfe",  2, 0.0, 25.0, 1.0);
    b2_prms.add_prm("qde",  2, 0.0, 25.0, 1.0);

    b2_prms.add_prm("l5h", -1, 0.05, 1.0,  0.01);
    b2_prms.add_prm("l6h", -1, 0.05, 1.0,  0.01);
    b2_prms.add_prm("l7",  -1, 0.05, 1.0,  0.01);
    b2_prms.add_prm("l8",  -1, 0.05, 1.0,  0.01);

    // b2_prms.add_prm("bm",  -2,  0.05, 0.5, 1.0);
    // b2_prms.add_prm("qde", -2,  0.05, 0.5, 1.0);
    // b2_prms.add_prm("qm",  -2,  0.05, 0.5, 1.0);

    b2_prms.bn_tol = 1e-6; b2_prms.svd_cut = 1e-8; b2_prms.step = 1.0;

    no_sxt();
    fit_match(b2_prms);
  }
}
