#include <assert.h>

#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


void cod_stat
(const std::vector<int> &bpm, double mean[], double sigma[], double peak[])
{
  int    n;
  double sum[2], sum2[2];

  for (int k = 0; k < 2; k++) {
    sum[k]  = 0e0;
    sum2[k] = 0e0;
    peak[k] = 0e0;
  }

  n = 0;
  for (int j = 0; j < bpm.size(); j++) {
    n++;
    for (int k = 0; k < 2; k++) {
      sum[k] += Cell[bpm[j]].BeamPos[k*2];
      sum2[k] += sqr(Cell[bpm[j]].BeamPos[k*2]);
      peak[k] = max(peak[k], fabs(Cell[bpm[j]].BeamPos[k*2]));
    }
  }

  for (int k = 0; k < 2; k++) {
    mean[k]  = (n != 0)? sum[k]/n : -1e0;
    sigma[k] = (n != 0 && n != 1)? (n*sum2[k]-sqr(sum[k]))/(n*(n-1e0)) : 0e0;
    sigma[k] = (sigma[k] >= 0e0)? sqrt(sigma[k]) : -1e0;
  }
}


void print_cod(void)
{
  printf("\nClosed Orbit:\n [");
  for (int k = 0; k < 2*nd_tps; k++) {
    printf("%13.6e", globval.CODvect[k]);
    printf((k < 2*nd_tps-1)? ", ":"]\n");
  }
}


std::vector<int> get_bpm(void)
{
  std::vector<int> bpm;

  for (int k = 0; k <= globval.Cell_nLoc; k++)
    if (strncmp("bpm", Cell[k].Elem.PName, 3) == 0)
      bpm.push_back(k);
  printf("\nget_bpm: number of bpms %lu\n", bpm.size());
  return bpm;
}


void compute_Ddisp(const double delta, const std::vector<int> &bpm)
{
  const string file_name = "Ddisp.out";
  const int n = bpm.size();

  std::vector<double> disp(n);
  FILE                *fp;

  fp = file_write(file_name.c_str());
  Ring_GetTwiss(true, delta);
  for (int k = 0; k < n; k++)
    disp[k] = Cell[bpm[k]].Eta[X_];
  Ring_GetTwiss(true, 0e0);
  for (int k = 0; k < n; k++) {
    disp[k] -= Cell[bpm[k]].Eta[X_];
    fprintf
      (fp, "%4d %15s %9.5f %4.1f %8.5f\n",
       bpm[k], Cell[bpm[k]].Elem.PName, Cell[bpm[k]].S, get_code(Cell[k]),
       disp[k]);
  }
  fclose(fp);
}


double get_V_RF(const int Fnum)
{
  const int loc = Elem_GetPos(Fnum, 1);
  return Cell[loc].Elem.C->V_RF;
}


double get_f_RF(const int Fnum)
{
  const int loc = Elem_GetPos(Fnum, 1);
  return Cell[loc].Elem.C->f_RF;
}


int get_h_RF(const int Fnum)
{
  const int loc = Elem_GetPos(Fnum, 1);
  return Cell[loc].Elem.C->harm_num;
}


void set_f_RF(const int Fnum, const double f_RF)
{
  const int loc = Elem_GetPos(Fnum, 1);
  Cell[loc].Elem.C->f_RF = f_RF;
}


double nu_s_track_fft(const int n, const double delta_ampl)
{
  // Obtain synchrotron tune from tracking & FFT analysis.
  const bool
    prt = false;
  const string
    file_name_1 = "track.out",
    file_name_2 = "fft.out";
  const int
    ps_dim = 6,
    n_peak = 3;

  long int        lastpos;
  double          ps_buf[ps_dim][n], nu[n_peak], A[n_peak];
  ss_vect<double> ps;
  ofstream        outf;

  if (prt)
    file_wr(outf, file_name_1.c_str());
  ps.zero();
  ps += globval.CODvect;
  ps[delta_] += delta_ampl;
  for (int k = 0; k < ps_dim; k++)
    ps_buf[k][0] = ps[k];
  if (prt)
    outf << scientific << setprecision(3) << setw(4) << 0 << setw(12) << ps
	 << "\n";
  for (int j = 1; j < n; j++) {
    for (int k = 0; k < ps_dim; k++)
      ps_buf[k][j] = ps[k];
    Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
    if (prt)
      outf << scientific << setprecision(3) << setw(4) << j << setw(12) << ps
	   << "\n";
  }
  if (prt)
    outf.close();

  rm_mean(n, ps_buf[ct_]);
  sin_FFT(n, ps_buf[ct_]);
  GetPeaks(n, ps_buf[ct_], 1, nu, A);

  if (prt) {
    file_wr(outf, file_name_2.c_str());
    for (int k = 0; k < n/2+1; k++)
      outf << setprecision(5)
	   << fixed << setw(8) << (double)k/(double)n
	   << scientific << setw(13) << ps_buf[ct_][k] << "\n";
    outf.close();
  }

  return nu[0];
}

double compute_nu_s(const int Fnum)
{
  // Compute synchrotron tune from analytic formula.
  const int
    loc     = Elem_GetPos(Fnum, 1),
    h       = Cell[loc].Elem.C->harm_num;
  const double
    E_0     = 1e9*globval.Energy,
    U_0     = globval.U0,
    V_RF    = Cell[loc].Elem.C->V_RF,
    alpha_c = globval.Alphac,
    phi_0   = fabs(asin(U_0/V_RF));

  return sqrt(h*alpha_c*V_RF*cos(phi_0*M_PI/180e0)/(2e0*M_PI*E_0));
}


void rf_gymnastics
(const int Fnum, const std::vector<int> &bpm,
 const std::vector<double> &alpha_c)
{
  // Remark: Globval.pathlength must be set to true in GetEmittance in
  //         nsls-ii_lib.

  const string
    file_name = "f_RF_scan.out";
  const int
    n[]       = {-15, 10};
  const double
    f_RF_step = 1e3;

  long int
    lastpos;
  int
    h_RF;
  double
    f_0_RF, df_RF, f_RF, df_RF_est, f_rev, mean[2], sigma[2], peak[2],
    dct, delta;
  FILE
    *fp;

  trace = false;

  // fp = file_write(file_name.c_str());
  fp = stdout;

  globval.Cavity_on = globval.radiation = globval.pathlength = true;

  f_0_RF = get_f_RF(Fnum);
  h_RF = get_h_RF(Fnum);

  fprintf(fp, "\n");
  fprintf(fp,
	  "# df_RF    df_RF    delta   ct      nu_s       f_s   hor orbit"
	  "  sigma_ct   nu_x   nu_y\n");
  fprintf(fp,
	  "# [kHz]  estimated   [%%]    [m]               [kHz]   rms [mm]"
	  "   [psec]\n");
  for (int j = n[0]; j <= n[1]; j++) {
    df_RF = j*f_RF_step;
    f_RF = f_0_RF + df_RF;
    set_f_RF(Fnum, f_RF);

    getcod(0e0, lastpos);

    cod_stat(bpm, mean, sigma, peak);
    dct = globval.CODvect[ct_];
    delta = globval.CODvect[delta_];

    df_RF_est = 0e0;
    for (int k = 0; k < alpha_c.size(); k++)
      df_RF_est -= alpha_c[k]*pow(delta, k+1);
    df_RF_est *= f_0_RF;
    f_rev = f_RF/h_RF;

    GetEmittance(Fnum, true, false);

    fprintf
      (fp,
       "  %5.1f    %5.1f   %6.3f  %5.3f  %9.7f   %5.3f    %5.3f     %5.3f"
       "   %6.3f  %5.3f\n",
       1e-3*df_RF, 1e-3*df_RF_est, 1e2*delta, dct,
       -globval.Omega, -1e-3*globval.Omega*f_rev, 1e3*sigma[X_],
       1e12*sqrt(Cell[0].sigma[ct_][ct_])/c0,
       globval.TotalTune[X_], globval.TotalTune[Y_]);
  }
}


void compare_nu_s(const int Fnum)
{
  long int lastpos;

  globval.Cavity_on = true;
  Ring_GetTwiss(true, 0e0);
  printglob();
  getcod(0e0, lastpos);
  cout << scientific << setprecision(5)
       << "\n" << setw(13) << globval.CODvect << "\n";
  printf("\n  nu_s = %9.7f\n", compute_nu_s(Fnum));

  globval.radiation = true;
  Ring_GetTwiss(true, 0e0);
  printglob();
  getcod(0e0, lastpos);
  cout << scientific << setprecision(5)
       << "\n" << setw(13) << globval.CODvect << "\n";
  printf("\n  nu_s = %9.7f\n", compute_nu_s(Fnum));
}


void prt_lin_opt(void)
{
  const int
    n_bpm = 7;
  const string
    bpm[n_bpm] = {
      "BPMZ5D8R", "BPMZ6D8R", "BPMZ7D8R", "BPMZ1T8R", "BPMZ2T8R", "BPMZ3T8R",
      "BPMZ4T8R"
    };

  bool   first = true;
  int    loc;
  double beta_j[2], beta_rel[2], nu_j[2], dnu[2];

  printf("\n  BPM     Rel. beta_x  Rel. beta_y  dnu_x  dnu_y\n");
  for (int j = 0; j < n_bpm; j++) {
    loc = Elem_GetPos(ElemIndex(bpm[j]), 1);
    if (first) {
      for (int k = 0; k < 2; k++) {
	beta_j[k] = Cell[loc].Beta[k];
	nu_j[k] = Cell[loc].Nu[k];
	dnu[k] = 0e0;
      }
      first = false;
    }
    for (int k = 0; k < 2; k++) {
      beta_rel[k] = Cell[loc].Beta[k]/beta_j[k];
      beta_j[k] = Cell[loc].Beta[k];
      dnu[k] = Cell[loc].Nu[k] - nu_j[k];
      nu_j[k] = Cell[loc].Nu[k];
    }
    printf("%-10s   %5.3f        %5.3f     %5.3f  %5.3f\n",
	   bpm[j].c_str(), beta_rel[X_], beta_rel[Y_], dnu[X_], dnu[Y_]);
  }
}


void fit_ksi_1
(const std::vector<int> &Fnum_b3, const double ksi_x, const double ksi_y,
 const double db_3_L)
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
    set_dbnL_design_fam(Fnum_b3[k-1], Sext, db_3_L, 0e0);
    Ring_Getchrom(0e0);
    if (prt)
      printf("\nfit_ksi_1: ksi1+ = [%9.5f, %9.5f]\n",
	     globval.Chrom[X_], globval.Chrom[Y_]);

    for (j = 1; j <= m; j++)
      A[j][k] = globval.Chrom[j-1];
    set_dbnL_design_fam(Fnum_b3[k-1], Sext, -2e0*db_3_L, 0e0);
    Ring_Getchrom(0e0);
    if (prt)
      printf("fit_ksi_1: ksi1- = [%9.5f, %9.5f]\n",
	     globval.Chrom[X_], globval.Chrom[Y_]);
    for (j = 1; j <= 2; j++) {
      A[j][k] -= globval.Chrom[j-1]; A[j][k] /= 2e0*db_3_L;
    }

    set_dbnL_design_fam(Fnum_b3[k-1], Sext, db_3_L, 0e0);
  }

  Ring_Getchrom(0e0);
  if (prt)
    printf("\nfit_ksi_1: ksi1  = [%9.5f, %9.5f]\n",
	   globval.Chrom[X_], globval.Chrom[Y_]);
  for (j = 1; j <= 2; j++)
    b[j] = -(globval.Chrom[j-1]-ksi0[j-1]);

  dmcopy(A, m, n_b3, U); dsvdcmp(U, m, n_b3, w, V);

  printf("\nfit_ksi_1:\n  singular values:\n");
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


void fit_ksi_1(const double ksi_x, const double ksi_y)
{
  const double db_3_L = 1.0;

  std::vector<int> Fnum;

  Fnum.push_back(ElemIndex("S4M2D1R"));
  Fnum.push_back(ElemIndex("S3M2D1R"));
  Fnum.push_back(ElemIndex("S2M2D1R"));
  Fnum.push_back(ElemIndex("S1MT1R"));
  Fnum.push_back(ElemIndex("S2M1T1R"));
  Fnum.push_back(ElemIndex("S3M1T1R"));
  Fnum.push_back(ElemIndex("S4M1T1R"));
  Fnum.push_back(ElemIndex("S4M2T1R"));
  Fnum.push_back(ElemIndex("S3M2T1R"));
  Fnum.push_back(ElemIndex("S2M2T1R"));
  Fnum.push_back(ElemIndex("S1MD2R"));
  Fnum.push_back(ElemIndex("S2M1D2R"));
  Fnum.push_back(ElemIndex("S3M1D2R"));
  Fnum.push_back(ElemIndex("S4M1D2R"));
  Fnum.push_back(ElemIndex("S4M2D2R"));
  Fnum.push_back(ElemIndex("S3M2D2R"));
  Fnum.push_back(ElemIndex("S2M2D2R"));
  Fnum.push_back(ElemIndex("S1MT2R"));
  Fnum.push_back(ElemIndex("S2M1T2R"));
  Fnum.push_back(ElemIndex("S3M1T2R"));
  Fnum.push_back(ElemIndex("S4M1T2R"));
  Fnum.push_back(ElemIndex("S4M2T2R"));
  Fnum.push_back(ElemIndex("S3M2T2R"));
  Fnum.push_back(ElemIndex("S2M2T2R"));
  Fnum.push_back(ElemIndex("S1MD3R"));
  Fnum.push_back(ElemIndex("S2M1D3R"));
  Fnum.push_back(ElemIndex("S3M1D3R"));
  Fnum.push_back(ElemIndex("S4M1D3R"));
  Fnum.push_back(ElemIndex("S4M2D3R"));
  Fnum.push_back(ElemIndex("S3M2D3R"));
  Fnum.push_back(ElemIndex("S2M2D3R"));
  Fnum.push_back(ElemIndex("S1MT3R"));
  Fnum.push_back(ElemIndex("S2M1T3R"));
  Fnum.push_back(ElemIndex("S3M1T3R"));
  Fnum.push_back(ElemIndex("S4M1T3R"));
  Fnum.push_back(ElemIndex("S4M2T3R"));
  Fnum.push_back(ElemIndex("S3M2T3R"));
  Fnum.push_back(ElemIndex("S2M2T3R"));
  Fnum.push_back(ElemIndex("S1MD4R"));
  Fnum.push_back(ElemIndex("S2M1D4R"));
  Fnum.push_back(ElemIndex("S3M1D4R"));
  Fnum.push_back(ElemIndex("S4M1D4R"));
  Fnum.push_back(ElemIndex("S4M2D4R"));
  Fnum.push_back(ElemIndex("S3M2D4R"));
  Fnum.push_back(ElemIndex("S2M2D4R"));
  Fnum.push_back(ElemIndex("S1MT4R"));
  Fnum.push_back(ElemIndex("S2M1T4R"));
  Fnum.push_back(ElemIndex("S3M1T4R"));
  Fnum.push_back(ElemIndex("S4M1T4R"));
  Fnum.push_back(ElemIndex("S4M2T4R"));
  Fnum.push_back(ElemIndex("S3M2T4R"));
  Fnum.push_back(ElemIndex("S2M2T4R"));
  Fnum.push_back(ElemIndex("S1MD5R"));
  Fnum.push_back(ElemIndex("S2M1D5R"));
  Fnum.push_back(ElemIndex("S3M1D5R"));
  Fnum.push_back(ElemIndex("S4M1D5R"));
  Fnum.push_back(ElemIndex("S4M2D5R"));
  Fnum.push_back(ElemIndex("S3M2D5R"));
  Fnum.push_back(ElemIndex("S2M2D5R"));
  Fnum.push_back(ElemIndex("S1MT5R"));
  Fnum.push_back(ElemIndex("S2M1T5R"));
  Fnum.push_back(ElemIndex("S3M1T5R"));
  Fnum.push_back(ElemIndex("S4M1T5R"));
  Fnum.push_back(ElemIndex("S4M2T5R"));
  Fnum.push_back(ElemIndex("S3M2T5R"));
  Fnum.push_back(ElemIndex("S2M2T5R"));
  Fnum.push_back(ElemIndex("S1MD6R"));
  Fnum.push_back(ElemIndex("S2M1D6R"));
  Fnum.push_back(ElemIndex("S3M1D6R"));
  Fnum.push_back(ElemIndex("S4M1D6R"));
  Fnum.push_back(ElemIndex("S4M2D6R"));
  Fnum.push_back(ElemIndex("S3M2D6R"));
  Fnum.push_back(ElemIndex("S2M2D6R"));
  Fnum.push_back(ElemIndex("S1MT6R"));
  Fnum.push_back(ElemIndex("S2M1T6R"));
  Fnum.push_back(ElemIndex("S3M1T6R"));
  Fnum.push_back(ElemIndex("S4M1T6R"));
  Fnum.push_back(ElemIndex("S4M2T6R"));
  Fnum.push_back(ElemIndex("S3M2T6R"));
  Fnum.push_back(ElemIndex("S2M2T6R"));
  Fnum.push_back(ElemIndex("S1MD7R"));
  Fnum.push_back(ElemIndex("S2M1D7R"));
  Fnum.push_back(ElemIndex("S3M1D7R"));
  Fnum.push_back(ElemIndex("S4M1D7R"));
  Fnum.push_back(ElemIndex("S4M2D7R"));
  Fnum.push_back(ElemIndex("S3M2D7R"));
  Fnum.push_back(ElemIndex("S2M2D7R"));
  Fnum.push_back(ElemIndex("S1MT7R"));
  Fnum.push_back(ElemIndex("S2M1T7R"));
  Fnum.push_back(ElemIndex("S3M1T7R"));
  Fnum.push_back(ElemIndex("S4M1T7R"));
  Fnum.push_back(ElemIndex("S4M2T7R"));
  Fnum.push_back(ElemIndex("S3M2T7R"));
  Fnum.push_back(ElemIndex("S2M2T7R"));
  Fnum.push_back(ElemIndex("S1MD8R"));
  Fnum.push_back(ElemIndex("S4M2D1R"));

  fit_ksi_1(Fnum, ksi_x, ksi_y, db_3_L);
}


void set_lat_state(void)
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
  const std::vector<string>
    cav = {"cavh1t8r", "cavh2t8r", "cavh3t8r", "cavh4t8r"};
  const std::vector<double>
    alpha_c = {7.038e-04, 5.164e-04, 3.293e-02, -2.415e-01};

  std::vector<int>    bpm;
  std::vector<double> disp;

  reverse_elem     = true;
  globval.mat_meth = false;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  set_lat_state();

  Ring_GetTwiss(true, 0e0);
  printglob();

  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);
  prtmfile("flat_file.dat");

  if (!false)
    prt_lin_opt();

  if (false) {
    fit_ksi_1(0.0, 0.0);
    get_ksi2(3e-2, 20);
  }

  if (false)
    compare_nu_s(ElemIndex(cav[0]));

  if (false)
    GetEmittance(ElemIndex(cav[0]), false, true);

  if (false) {
    bpm = get_bpm();
    compute_Ddisp(2.0e-2, bpm);
  }

  if (false) {
    bpm = get_bpm();
    rf_gymnastics(ElemIndex(cav[0]), bpm, alpha_c);
  }
}
