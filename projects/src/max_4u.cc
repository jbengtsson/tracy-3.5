#define NO 1

#include <assert.h>

#include "tracy_lib.h"


int no_tps = NO;


void set_ps_rot(const string &fam_name, const double dnu_x, const double dnu_y)
{
  const int
    Fnum = ElemIndex(fam_name.c_str());
  const double
    dnu_0[] = {0e0, 0e0},
    dnu[]   = {dnu_x, dnu_y};

  set_map(Fnum, dnu_0);
  printf("\ntweak_nu:\n");
  printf("  dnu = [%8.5f, %8.5f]\n", dnu[X_], dnu[Y_]);
  set_map(Fnum, dnu);
  Ring_GetTwiss(true, 0e0);
  printglob();
 }


void get_dnu_straight(const int loc)
{
  const double
    dnu[] = {
      Cell[globval.Cell_nLoc].Nu[X_]-Cell[loc].Nu[X_],
      Cell[globval.Cell_nLoc].Nu[Y_]-Cell[loc].Nu[Y_]
    };

  printf("\nget_dnu_straight:\n");
  printf(" nu               = [%7.5f, %7.5f]\n",
	 Cell[globval.Cell_nLoc].Nu[X_], Cell[globval.Cell_nLoc].Nu[Y_]);
  printf(" dnu              = [%7.5f, %7.5f]\n",
	 Cell[loc].Nu[X_],Cell[loc].Nu[Y_] ),
  printf(" dnu_1/2_straight = [%7.5f, %7.5f]\n", dnu[X_], dnu[Y_]);
}


void fit_ksi_jb(const std::vector<int> &Fnum_b3,
	      const double ksi_x, const double ksi_y, const double db3L)
{
  int    n_b3, j, k;
  double **A, **U, **V, *w, *b, *x, b3, a3;

  const bool   prt = !false;
  const int    m   = 2;
  const double
    ksi0[]  = {ksi_x, ksi_y},
    svd_cut = 1e-10;

  n_b3 = Fnum_b3.size();

  A = dmatrix(1, m, 1, n_b3); U = dmatrix(1, m, 1, n_b3);
  V = dmatrix(1, n_b3, 1, n_b3);
  w = dvector(1, n_b3); b = dvector(1, m); x = dvector(1, n_b3);

  // Zero sextupoles to track linear chromaticity.
  if (false) no_sxt();

  for (k = 1; k <= n_b3; k++) {
    set_dbnL_design_fam(Fnum_b3[k-1], Sext, db3L, 0e0);
    Ring_Getchrom(0e0);
    if (prt)
      printf("\nfit_ksi_jb: ksi1+ = [%9.5f, %9.5f]\n",
	     globval.Chrom[X_], globval.Chrom[Y_]);

    for (j = 1; j <= m; j++)
      A[j][k] = globval.Chrom[j-1];
    set_dbnL_design_fam(Fnum_b3[k-1], Sext, -2e0*db3L, 0e0);
    Ring_Getchrom(0e0);
    if (prt)
      printf("fit_ksi_jb: ksi1- = [%9.5f, %9.5f]\n",
	 globval.Chrom[X_], globval.Chrom[Y_]);
    for (j = 1; j <= 2; j++) {
      A[j][k] -= globval.Chrom[j-1]; A[j][k] /= 2e0*db3L;
    }

    set_dbnL_design_fam(Fnum_b3[k-1], Sext, db3L, 0e0);
  }

  Ring_Getchrom(0e0);
  if (prt)
    printf("\nfit_ksi_jb: ksi1  = [%9.5f, %9.5f]\n",
	   globval.Chrom[X_], globval.Chrom[Y_]);
  for (j = 1; j <= 2; j++)
    b[j] = -(globval.Chrom[j-1]-ksi0[j-1]);

  dmcopy(A, m, n_b3, U); dsvdcmp(U, m, n_b3, w, V);

  printf("\nfit_ksi_jb:\n  singular values:\n");
  for (j = 1; j <= n_b3; j++) {
    printf("    %9.3e", w[j]);
    if (w[j] < svd_cut) {
      w[j] = 0e0;
      printf(" (zeroed)");
    }
    printf("\n");
  }

  dsvbksb(U, w, V, m, n_b3, b, x);

  for (k = 1; k <= n_b3; k++)
    set_dbnL_design_fam(Fnum_b3[k-1], Sext, x[k], 0e0);

  if (prt) {
    printf("\n  b3:\n");
    for (k = 0; k < n_b3; k++) {
      get_bn_design_elem(Fnum_b3[k], 1, Sext, b3, a3);
      printf("    %-8s %10.5f\n", ElemFam[Fnum_b3[k]-1].ElemF.PName, b3);
    }
    printf("\n");
  }

  free_dmatrix(A, 1, m, 1, n_b3); free_dmatrix(U, 1, m, 1, n_b3);
  free_dmatrix(V, 1, n_b3, 1, n_b3);
  free_dvector(w, 1, n_b3); free_dvector(b, 1, m); free_dvector(x, 1, n_b3);
}


void fit_ksi_jb(const int lat_case, const double ksi_x, const double ksi_y)
{
  std::vector<int> Fnum;

  Fnum.push_back(ElemIndex("sfoh"));
  Fnum.push_back(ElemIndex("sdqd"));
  // Fnum.push_back(ElemIndex("sd2"));
  
  fit_ksi_jb(Fnum, 0e0, 0e0, 1e0);
}


void track(const int n_turn, const double Ax, const double Ay)
{
  const string
    file_name = "track.dat";
  const double
    A[] = {Ax, Ay};

  long int        lastpos;
  ss_vect<double> ps;
  ofstream        outf;

  file_wr(outf, file_name.c_str());
  ps.zero();
  for (int k = 0; k < 2; k++)
    ps[2*k] = A[k];
  outf << "# turn          x                     p_x"
       << "                     y                     p_y"
       << "                   delta                   c*t\n"
       << "#              [m]                   [rad]"
       << "                   [m]                   [rad]"
       << "                                          [m] \n";
    outf << scientific << setprecision(14)
	 << setw(4) << 0 << setw(23) << ps << "\n";
  for (int k = 1; k <= n_turn; k++) {
    Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
    outf << scientific << setprecision(14)
	 << setw(4) << k << setw(23) << ps << "\n";
  }
  outf.close();
}


void chk_optics(const double alpha_x, const double beta_x,
		const double eta_x, const double etap_x,
		const double alpha_y, const double beta_y,
		const double eta_y, const double etap_y)
{
  Vector2 alpha, beta, eta, etap;

  alpha[X_] = alpha_x; alpha[Y_] = alpha_y;
  beta[X_]  = beta_x;  beta[Y_]  = beta_y;
  eta[X_]   = eta_x;   eta[Y_]   = eta_y;
  etap[X_]  = etap_x;  etap[Y_]  = etap_y;

  ttwiss(alpha, beta, eta, etap, 0e0);
}


void set_state(void)
{
  globval.H_exact        = false;
  globval.quad_fringe    = false;
  globval.Cavity_on      = false;
  globval.radiation      = false;
  globval.emittance      = false;
  globval.IBS            = false;
  globval.pathlength     = false;
  globval.Aperture_on    = false;
  globval.Cart_Bend      = false;
  globval.dip_edge_fudge = true;
}


int main(int argc, char *argv[])
{
  int loc;

  globval.mat_meth = false;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  set_state();

  if (false) {
    no_mult(Sext);
    no_mult(Oct);
  }

  if (false) {
    fit_ksi_jb(1, 0e0, 0e0);
    assert(false);
  }

  Ring_GetTwiss(true, 0e0);
  printglob();

  if (false) {
    // A 1/2 ps_rot at the entrance & exit of the super period for a symmetric
    // approach.
    // set_ps_rot("ps_rot1", 0.1491/2e0, -0.0413/2e0);
    // set_ps_rot("ps_rot2", 0.3991/2e0, -0.2914/2e0);

    set_ps_rot("ps_rot", 0.06/2.0, 0.0/2.0);
  }

  prtmfile("flat_file.dat");
  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);
  prt_chrom_lat("chromlat.out");

  if (false)
    GetEmittance(ElemIndex("cav"), false, true);

  if (false) {
    loc = Elem_GetPos(ElemIndex("sd2"), 2);
    printf("\nloc = %d\n", loc);
    get_dnu_straight(loc);
  }

  if (false) {
    chk_optics(0.0, 13.24261, 0.0, 0.0, 0.0, 2.35728, 0.0, 0.0);
    prt_lat("linlat1.out", globval.bpm, true);
    prt_lat("linlat.out", globval.bpm, true, 10);
    prtmfile("flat_file.dat");
  }

  if (!false) {
    globval.Cavity_on = false;
    track(100, -6e-3, 0e0);
  }
}
