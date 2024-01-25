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


double get_f_RF(const int Fnum)
{
  int loc = Elem_GetPos(Fnum, 1);
  return Cell[loc].Elem.C->f_RF;
}


void set_f_RF(const int Fnum, const double f_RF)
{
  int loc = Elem_GetPos(Fnum, 1);
  Cell[loc].Elem.C->f_RF = f_RF;
}


void rf_gymnastics(const int Fnum, const std::vector<int> &bpm)
{
  // Remark: Globval.pathlength must be set to true in GetEmittance in
  //         nsls-ii_lib.

  const string
    file_name = "f_RF_scan.out";
  const int
    n[] = {-15, 10};
  const double
    f_RF_step = 1e3,
    C         = Cell[globval.Cell_nLoc].S;

  long int lastpos;
  double   f_RF, Df_RF, f_0, mean[2], sigma[2], peak[2];
  FILE     *fp;

  trace = false;

  // fp = file_write(file_name.c_str());
  fp = stdout;

  globval.Cavity_on = globval.radiation = globval.pathlength = true;

  f_RF = get_f_RF(Fnum);

  fprintf(fp, "\n");
  fprintf(fp, "# Df_RF   f_s    rms hor    sigma_ct    nu_x   nu_y\n");
  fprintf(fp, "# [kHz]  [kHz]  orbit [mm]  [picosec]\n");
  for (int k = n[0]; k <= n[1]; k++) {
    Df_RF = k*f_RF_step;
    set_f_RF(Fnum, f_RF+Df_RF);
    get_f_RF(Fnum);

    getcod(0e0, lastpos);
    f_0 = c0/(C+globval.CODvect[ct_]);
    cod_stat(bpm, mean, sigma, peak);
    Ring_GetTwiss(true, 0e0);
    GetEmittance(ElemIndex("cavh1t8r"), true, !true);
    fprintf(fp, "  %5.1f  %5.3f    %5.3f      %5.3f    %6.3f  %5.3f\n",
	    1e-3*Df_RF, -1e-3*globval.Omega*f_0,
	    1e3*sigma[X_],
	    1e12*sqrt(Cell[0].sigma[ct_][ct_])/c0,
	    globval.TotalTune[X_], globval.TotalTune[Y_]);
  }
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

  if (false)
    GetEmittance(ElemIndex(cav[0]), false, true);

  if (false) {
    bpm = get_bpm();
    compute_Ddisp(2.0e-2, bpm);
  }

  if (!false) {
    bpm = get_bpm();
    rf_gymnastics(ElemIndex(cav[0]), bpm);
  }
}
