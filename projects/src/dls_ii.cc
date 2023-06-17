#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


const double eps_x       = 150e-3,               // Hor. emittance [nm.rad].
             nu_uc[]     = {2.43, 0.89},         // Cell tune.
             nu_sc[]     = {2.43, 0.89},         // Super period tune.
             L_uc        = 1.25,                 // Unit Cell length.
             L_ss        = 10.5,                 /* Super Period length;
                                                    with one unit cell. */
             eta_cuc[]   = {0.04, 0.0},          /* Linear Optics
						    Center of Unit Cell. */
             etap_cuc[]  = {0.0, 0.0},
             alpha_cuc[] = {0.0, 0.0},
             beta_cuc[]  = {5.6, 1.9},
             beta_cs[]   = {3.0, 3.0};          /* Linear Optics
						    Center of Straight. */

int              n_b3;
std::vector<int> b3_Fnum, loc;

const int n_prt = 9;

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

  const bool prt = false;

  if (prt) printf("set_prm:\n");
  for (i = 1; i <= n_prm; i++) {
    if (n[i-1] > 0)
      for (j = 1; j <= GetnKid(Fnum[i-1]); j++)
	set_bn_design_elem(Fnum[i-1], j, n[i-1], bn[i], 0e0);
    else if (n[i-1] == -1) {
      set_L(Fnum[i-1], bn[i]); get_S();
    } else if (n[i-1] == -2)
      // set_bn_s(-Fnum[i-1], n[i-1], bn[i]);
      ;
    if (prt) {
      printf(" %12.5e", bn_scl[i-1]*bn[i]);
      if (i % n_prt == 0) printf("\n");
    }
  }
  if (prt &&(n_prm % n_prt != 0)) printf("\n");
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


void nu_cell_tweak(const double dnu_diff_min)
{
  bool          stable[2];
  long int      lastpos;
  int           k;
  double        dnu_diff[2], beta0[2], beta1[2], dnu[2];
  ss_vect<tps>  M;

  const bool prt = false;

  for (k = 0; k < 2; k++)
    // Avoid half integer resonances.
    dnu_diff[k] = (2e0*globval.TotalTune[k]-nint(2e0*globval.TotalTune[k]))/2e0;
  if (prt) printf("\nnu = [%8.5f, %8.5f] dnu_diff = [%8.5f, %8.5f]\n",
		  globval.TotalTune[X_], globval.TotalTune[Y_],
		  dnu_diff[X_], dnu_diff[Y_]);

  for (k = 0; k < 2; k++)
    if (fabs(dnu_diff[k]) < dnu_diff_min)
      dnu[k] = sgn(dnu_diff[k])*dnu_diff_min - dnu_diff[k];
  set_map(ElemIndex("M"), dnu);

  M.identity();
  Cell_Pass(0, globval.Cell_nLoc, M, lastpos);
  get_map_twiss(M, beta0, beta1, globval.TotalTune, stable);

  for (k = 0; k < 2; k++)
    dnu_diff[k] =
      (2e0*globval.TotalTune[k]-nint(2e0*globval.TotalTune[k]))/2e0;
  if (prt) printf("nu = [%8.5f, %8.5f] dnu_diff = [%8.5f, %8.5f]"
		  " dnu = [%8.5f, %8.5f]\n",
		  globval.TotalTune[X_], globval.TotalTune[Y_],
		  dnu_diff[X_], dnu_diff[Y_], dnu[X_], dnu[Y_]);
}


double get_eps_x1(void)
{
  // Evaluate emittance [nm.rad] and damping partition from synchrotron
  // integrals.
  bool         cav, emit;
  long int     lastpos;
  double       I[6], eps_x;
  ss_vect<tps> A;

  const bool prt = false;

  cav = globval.Cavity_on; emit = globval.emittance;

  globval.Cavity_on = false; globval.emittance = false;

  Ring_GetTwiss(false, 0e0);
  if (!globval.stable) {
    printf("\nget_eps_x1: unstable\n");
    return -1e0;
  }

  A = putlinmat(6, globval.Ascr);

  // prt_lin_map(3, A);

  globval.emittance = true;

  Cell_Pass(0, globval.Cell_nLoc, A, lastpos);

  get_I(I, false);

  eps_x = 1470e0*pow(globval.Energy, 2)*I[5]/(I[2]-I[4]);

  if (prt)
    printf("eps_x = %5.3f pm.rad, J_x = %5.3f, J_z = %5.3f \n",
	   1e3*eps_x, 1e0-I[4]/I[2], 2e0+I[4]/I[2]);

  globval.Cavity_on = cav; globval.emittance = emit;

  return eps_x;
}


void fit_ksi1(const double ksi_x, const double ksi_y,
	      const std::vector<int> &Fnum,
	      const double db3L, const double eps, const int iter_max)
{
  int      j, k, n, n_iter;
  double   b3L, a3L, chi, ksi0[2], **A, **U, **V, *w, *b, *x;

  const bool   prt = false;
  const int    m = 2, n_prt = 8;
  const double ksi_ref[] = {ksi_x, ksi_y}, svd_eps = 1e-10;

  n = Fnum.size();

  A = dmatrix(1, m, 1, n); U = dmatrix(1, m, 1, n); V = dmatrix(1, n, 1, n);
  w = dvector(1, n); b = dvector(1, m); x = dvector(1, n);

  Ring_Getchrom(0e0);
  for (j = 1; j <= 2; j++) {
    ksi0[j-1] = globval.Chrom[j-1]; b[j] = -(ksi0[j-1]-ksi_ref[j-1]);
  }
  chi = sqrt(sqr(b[1])+sqr(b[2]));

  n_iter = 0; chi = 1e30;
  while ((chi > eps) && (n_iter < iter_max)) {
    n_iter++;
    for (k = 1; k <= n; k++) {
      Ring_Getchrom(0e0);
      set_dbnL_design_fam(Fnum[k-1], Sext, db3L, 0e0);
      Ring_Getchrom(0e0);
      for (j = 1; j <= 2; j++)
	A[j][k] = (globval.Chrom[j-1]-ksi0[j-1])/db3L;
      set_dbnL_design_fam(Fnum[k-1], Sext, -db3L, 0e0);
      Ring_Getchrom(0e0);
    }

    dmcopy(A, m, n, U); dsvdcmp(U, m, n, w, V);

  if (prt) printf("\nfit_ksi1 singular values:\n");
    for (j = 1; j <= n; j++) {
      if (prt) printf("%11.3e", w[j]);
      if (w[j] < svd_eps) {
	w[j] = 0e0;
	if (prt) printf(" (zeroed)");
      }
      if (prt && (j % n_prt == 0)) printf("\n");
    }
    if (prt && (n % n_prt != 0)) printf("\n");

    dsvbksb(U, w, V, m, n, b, x);

    for (k = 1; k <= n; k++) {
      set_dbnL_design_fam(Fnum[k-1], Sext, x[k], 0e0);
      get_bnL_design_elem(Fnum[k-1], 1, Sext, b3L, a3L);
    }

    Ring_Getchrom(0e0);
    for (j = 1; j <= 2; j++) {
      ksi0[j-1] = globval.Chrom[j-1]; b[j] = -(ksi0[j-1]-ksi_ref[j-1]);
    }
    chi = sqrt(sqr(b[1])+sqr(b[2]));
  }
  if (prt) printf("\nfit_ksi1: ksi = [%12.5e, %12.5e]\n",
		  globval.Chrom[0], globval.Chrom[1]);

  free_dmatrix(A, 1, m, 1, n); free_dmatrix(U, 1, m, 1, n);
  free_dmatrix(V, 1, n, 1, n);
  free_dvector(w, 1, n); free_dvector(b, 1, m); free_dvector(x, 1, n);
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


void prt_b2(FILE *outf, const param_type &b2_prms)
{
  long int loc;
  int      k;
  double   b2, a2, phi;

  fprintf(outf, "\n");
  for (k = 0; k < (int)b2_prms.n_prm; k++) {
    loc = Elem_GetPos(b2_prms.Fnum[k], 1);
    get_bn_design_elem(b2_prms.Fnum[k], 1, Quad, b2, a2);
    if (Cell[loc].Elem.M->n_design == Dip) {
      phi = Cell[loc].Elem.M->Pirho*Cell[loc].Elem.PL*180e0/M_PI;
      fprintf(outf,
	      "%-8s: Bending, L = %7.5f, T = %7.5f, T1 = %7.5f, T2 = %7.5f"
	      ", K = %12.5e, N = Ndip, Method = Meth;\n",
	      Cell[loc].Elem.PName, Cell[loc].Elem.PL, phi,
	      Cell[loc].Elem.M->PTx1, Cell[loc].Elem.M->PTx2, b2);
    } else
      fprintf(outf,
	      "%-8s: Quadrupole, L = %7.5f, K = %12.5e"
	      ", N = Nquad, Method = Meth;\n",
	      Cell[loc].Elem.PName, Cell[loc].Elem.PL, b2);
  }
}


void prt_emit(const string &file_name,
	      const param_type &b2_prms, const double *b2)
{
  int    k, loc;
  double b3, a3;
  FILE   *outf;

  outf = file_write(file_name.c_str());

  prt_b2(outf, b2_prms);

  fprintf(outf, "\n");
  for (k = 0; k < n_b3; k++) {
    loc = Elem_GetPos(b3_Fnum[k], 1);
    get_bn_design_elem(b3_Fnum[k], 1, Sext, b3, a3);
    fprintf(outf,
	    "%-8s: sextupole, l = %7.5f"
	    ", k = %12.5e, n = nquad, Method = Meth;\n",
	    Cell[loc].Elem.PName, Cell[loc].Elem.PL, b3);
  }

  fclose(outf);
}


double f_emit(double *b2)
{
  // Optimize unit cell.
  // Lattice: unit cell.
  // Nota Bene: singular for crossing integer Cell Tunes.

  static double chi2_ref = 1e30;

  bool         stable[2];
  long int     lastpos;
  int          k;
  double       eps1_x, tr[2], b3L[n_b3], beta_b3L[n_b3], a3L, chi2;
  double       beta0[2], beta1[2];
  ss_vect<tps> M;

  const std::string file_name = "emit.out";
  const int         ksi_i_max = 10;
  const double      db3L = 1e-2, ksi_eps = 1e-3,
                    scl_eps = 1e0, scl_eta = 1e1, scl_ksi = 1e-6,
                    dnu_diff_min = 0.05;

  b2_prms.set_prm(b2);

  M.identity();
  Cell_Pass(0, Elem_GetPos(ElemIndex("M"), 1)-1, M, lastpos);
  // prt_lin_map(3, M);

  for (k = 0; k < 2; k++)
    tr[k] = M[2*k][2*k] + M[2*k+1][2*k+1];
  if ((fabs(tr[X_]) > 2e0) || (fabs(tr[Y_]) > 2e0)) {
    printf("\nf_emit: unstable [%5.3f, %5.3f]\n", tr[X_], tr[Y_]);
    chi2 = 1e10;
    return chi2;
  }

  get_map_twiss(M, beta0, beta1, globval.TotalTune, stable);

  // Tweak Cell Tune as needed. Update for each iteration; since Linear Optics
  // is changing.
  nu_cell_tweak(dnu_diff_min);

  eps1_x = get_eps_x1();

  fit_ksi1(0e0, 0e0, b3_Fnum, db3L, ksi_eps, ksi_i_max);

  for (k = 0; k < n_b3; k++)
    get_bnL_design_elem(b3_Fnum[k], 1, Sext, b3L[k], a3L);
  beta_b3L[0] = Cell[loc[1]].Beta[Y_]*b3L[0];
  beta_b3L[1] = Cell[loc[2]].Beta[X_]*b3L[1];
  beta_b3L[2] = Cell[loc[3]].Beta[Y_]*b3L[0];

  chi2 = 0e0;

  chi2 += sqr(scl_eta*Cell[loc[0]].Eta[X_]);
  chi2 += sqr(scl_eta*Cell[loc[0]].Etap[X_]);

  chi2 += sqr(scl_eps*(eps1_x-eps_x));

  for (k = 0; k < 3; k++)
    chi2 += sqr(scl_ksi*sqr(beta_b3L[k]));

  // for (k = 1; k <= b2_prms.n_prm; k++) {
  //   if (fabs(b2[k]) > b2_prms.bn_max[k-1]) chi2 += 1e10;
  // }

  if (chi2 < chi2_ref) {
    printf("\nchi2: %11.5e -> %11.5e\n", chi2_ref, chi2);
    printf("b:    %5.1f %10.3e %10.3e %7.5f %7.5f %10.3e %10.3e %10.3e\n",
	   1e3*eps1_x,
	   Cell[loc[0]].Eta[X_], Cell[loc[0]].Etap[X_],
	   globval.TotalTune[X_], globval.TotalTune[Y_],
	   beta_b3L[0], beta_b3L[1], beta_b3L[2]);
    printf("      %6.3f %6.3f %6.3f\n",
	   Cell[loc[1]].Beta[Y_], Cell[loc[2]].Beta[X_], Cell[loc[3]].Beta[Y_]);
    printf("      %9.5f %9.5f\n", b3L[0], b3L[1]);

    for (k = 1; k <= b2_prms.n_prm; k++) {
      printf(" %12.5e", b2[k]);
      if (k % n_prt == 0) printf("\n");
    }
    if (b2_prms.n_prm % n_prt != 0) printf("\n");

    prt_emit(file_name, b2_prms, b2);
    prtmfile("flat_file.fit");
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

  // Set initial directions and magnitude (unit vectors).
  for (i = 1; i <= n_b2; i++)
    for (j = 1; j <= n_b2; j++)
      xi[i][j] = (i == j)? 1e-2 : 0e0;

  dpowell(b2, xi, n_b2, 1e-8, &iter, &fret, f_emit);

  free_dvector(b2, 1, n_b2);  free_dvector(b2_lim, 1, n_b2);
  free_dmatrix(xi, 1, n_b2, 1, n_b2);
}


double f_stable(double *b2)
{
  // Optimize unit cell.
  // Lattice: unit cell.
  // Nota Bene: singular for crossing integer Cell Tunes.

  static double chi2_ref = 1e30;

  bool         stable[2];
  long int     lastpos;
  int          k;
  double       eps1_x, tr[2], chi2;
  double       beta0[2], beta1[2];
  ss_vect<tps> M;

  const std::string file_name = "stable.out";
  const double      scl_eps = 1e0, scl_eta = 1e1, scl_tr = 1e0;

  b2_prms.set_prm(b2);

  M.identity();
  Cell_Pass(0, Elem_GetPos(ElemIndex("M"), 1)-1, M, lastpos);
  // prt_lin_map(3, M);

  for (k = 0; k < 2; k++)
    tr[k] = M[2*k][2*k] + M[2*k+1][2*k+1];
  if ((fabs(tr[X_]) > 2e0) || (fabs(tr[Y_]) > 2e0))
    printf("\nf_emit: unstable [%5.3f, %5.3f]\n", tr[X_], tr[Y_]);

  get_map_twiss(M, beta0, beta1, globval.TotalTune, stable);

  eps1_x = get_eps_x1();

  chi2 = 0e0;

  chi2 += sqr(scl_eta*Cell[loc[0]].Eta[X_]);
  chi2 += sqr(scl_eta*Cell[loc[0]].Etap[X_]);

  chi2 += sqr(scl_tr*tr[X_]);
  chi2 += sqr(scl_tr*tr[Y_]);

  chi2 += sqr(scl_eps*(eps1_x-eps_x));

  if (chi2 < chi2_ref) {
    printf("\nchi2: %11.5e -> %11.5e\n", chi2_ref, chi2);
    printf("b:    %5.1f %10.3e %10.3e %7.5f %7.5f\n",
	   1e3*eps1_x,
	   Cell[loc[0]].Eta[X_], Cell[loc[0]].Etap[X_],
	   globval.TotalTune[X_], globval.TotalTune[Y_]);

    for (k = 1; k <= b2_prms.n_prm; k++) {
      printf(" %12.5e", b2[k]);
      if (k % n_prt == 0) printf("\n");
    }
    if (b2_prms.n_prm % n_prt != 0) printf("\n");

    prt_emit(file_name, b2_prms, b2);
    prtmfile("flat_file.fit");
    prt_lat("linlat1.out", globval.bpm, true);
    prt_lat("linlat.out", globval.bpm, true, 10);
  }

  chi2_ref = min(chi2, chi2_ref);

  return chi2;
}


void fit_stable(param_type &b2_prms)
{
  // Optimize unit cell.
  // Lattice: unit cell.

  int    n_b2, i, j, iter;
  double *b2, *b2_lim, **xi, fret;

  n_b2 = b2_prms.n_prm;

  b2 = dvector(1, n_b2); b2_lim = dvector(1, n_b2);
  xi = dmatrix(1, n_b2, 1, n_b2);

  b2_prms.ini_prm(b2, b2_lim);

  // Set initial directions and magnitude (unit vectors).
  for (i = 1; i <= n_b2; i++)
    for (j = 1; j <= n_b2; j++)
      xi[i][j] = (i == j)? 1e-2 : 0e0;

  dpowell(b2, xi, n_b2, 1e-8, &iter, &fret, f_stable);

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

  tr[X_] = globval.OneTurnMat[x_][x_] + globval.OneTurnMat[px_][px_];
  tr[Y_] = globval.OneTurnMat[y_][y_] + globval.OneTurnMat[py_][py_];
  // printf("trace: %6.3f %6.3f\n", tr[X_], tr[Y_]);

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
    printf("\nchi2: %11.5e, %11.5e\n", chi2, chi2_ref);
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

  int    i, loc1, loc2, loc3;
  double tr[2], chi2;

  const double scl_eta = 1e4, scl_alpha = 1e1,  scl_beta = 1e1,
               scl_L   = 1e1, scl_ksi   = 5e-1, scl_nu   = 1e3;

  b2_prms.set_prm(b2);

  // Center of Unit Cell.
  loc1 = Elem_GetPos(ElemIndex("ms"), 1);
  // End of 3rd Bend.
  loc2 = Elem_GetPos(ElemIndex("dq1"), 1);
  // Center of Straight.
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

  chi2 += sqr(scl_eta*Cell[loc2].Etap[X_]);

  chi2 += sqr(scl_eta*Cell[loc2].Eta[X_]);
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
    printf("\nchi2: %12.5e -> %12.5e\n", chi2_ref, chi2);
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
  int k;

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  reverse_elem = !false;

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

  if (!false) {
    b2_prms.add_prm("b1",    2, 0.0, 25.0, 1.0);

    b2_prms.add_prm("qf1",   2, 0.0, 25.0, 1.0);
    b2_prms.add_prm("qd2",   2, 0.0, 25.0, 1.0);
    b2_prms.add_prm("qd3",   2, 0.0, 25.0, 1.0);
    b2_prms.add_prm("qf421", 2, 0.0, 25.0, 1.0);
    b2_prms.add_prm("qf422", 2, 0.0, 25.0, 1.0);
    b2_prms.add_prm("qd5",   2, 0.0, 25.0, 1.0);
    b2_prms.add_prm("qf6",   2, 0.0, 25.0, 1.0);
    b2_prms.add_prm("qf8",   2, 0.0, 25.0, 1.0);

    b2_prms.bn_tol = 1e-6; b2_prms.svd_cut = 1e-8; b2_prms.step = 1.0;

    b3_Fnum.push_back(ElemIndex("sd2"));
    b3_Fnum.push_back(ElemIndex("sf2"));
    n_b3 = b3_Fnum.size();

    // Upstream of B2.
    loc.push_back(Elem_GetPos(ElemIndex("b2_5"), 1)-1);
    // Center of SD2 (SD3).
    loc.push_back(Elem_GetPos(b3_Fnum[0], 1));
    // Center of SF2.
    loc.push_back(Elem_GetPos(b3_Fnum[1], 1));
    // Center of SD2.
    loc.push_back(Elem_GetPos(b3_Fnum[0], 3));
    printf("\nloc:\n");
    for (k = 0; k < (int)loc.size(); k++)
      printf("\n  %10s %5.3f", Cell[loc[k]].Elem.PName, Cell[loc[k]].S);

    no_sxt();
    // fit_stable(b2_prms);
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

  if (false) {
    b2_prms.add_prm("qf1",   2, 0.0, 25.0, 1.0);
    b2_prms.add_prm("qd2",   2, 0.0, 25.0, 1.0);
    b2_prms.add_prm("qd3",   2, 0.0, 25.0, 1.0);
    b2_prms.add_prm("qf421", 2, 0.0, 25.0, 1.0);
    b2_prms.add_prm("qf422", 2, 0.0, 25.0, 1.0);
    b2_prms.add_prm("qd5",   2, 0.0, 25.0, 1.0);
    b2_prms.add_prm("qf6",   2, 0.0, 25.0, 1.0);
    b2_prms.add_prm("qf8",   2, 0.0, 25.0, 1.0);

    b2_prms.bn_tol = 1e-6; b2_prms.svd_cut = 1e-8; b2_prms.step = 1.0;

    no_sxt();
    fit_match(b2_prms);
  }
}
