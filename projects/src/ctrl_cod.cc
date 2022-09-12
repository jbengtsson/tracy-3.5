#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


void get_cod_rms(const double dxy_rms[], const int n_seed, const bool all)
{
  bool                cod;
  int                 i, j, k, n, n_cod;
  std::vector<double> x1[6], x2[6], x_mean[6], x_sigma[6];
  FILE                *fp;

  const int n_cod_corr = 5;

  globval.Cavity_on = false;

  for (j = 0; j <= globval.Cell_nLoc; j++)
    for (k = 0; k < 6; k++) {
      x1[k].push_back(0e0); x2[k].push_back(0e0);
    }
  
  fp = file_write("cod_rms.out");
  
  n_cod = 0;
  for (i = 0; i < n_seed; i++) {
    printf("\norb_corr: seed no %d\n", i+1);

    misalign_rms_type(Dip,  dxy_rms[X_], dxy_rms[Y_], 0e0, true);
    misalign_rms_type(Quad, dxy_rms[X_], dxy_rms[Y_], 0e0, true);
    
    cod = orb_corr(n_cod_corr);

    if (cod) {
      n_cod++;

      n = 0;
      for (j = 0; j <= globval.Cell_nLoc; j++)
	if (all || ((Cell[j].Elem.Pkind == Mpole) &&
		    (Cell[j].Elem.M->n_design == Sext))) {
	  n++;
	  for (k = 0; k < 6; k++) {
	    x1[k][n-1] += Cell[j].BeamPos[k];
	    x2[k][n-1] += sqr(Cell[j].BeamPos[k]);
	  }
	}
    } else
      printf("orb_corr: failed\n");

    zero_trims();
  }

  printf("\nget_cod_rms: no of seeds %d, no of cods %d\n", n_seed, n_cod);

  n = 0;
  for (j = 0; j <= globval.Cell_nLoc; j++)
    if (all || ((Cell[j].Elem.Pkind == Mpole) &&
		(Cell[j].Elem.M->n_design == Sext))) {
      n++;
      for (k = 0; k < 6; k++) {
	x_mean[k].push_back(x1[k][n-1]/n_cod);
	x_sigma[k].push_back(sqrt((n_cod*x2[k][n-1]-sqr(x1[k][n-1]))
				  /(n_cod*(n_cod-1.0))));
      }
      fprintf(fp, "%8.3f %6.2f %10.3e +/- %10.3e %10.3e +/- %10.3e\n",
	      Cell[j].S, get_code(Cell[j]),
	      1e3*x_mean[x_][n-1], 1e3*x_sigma[x_][n-1],
	      1e3*x_mean[y_][n-1], 1e3*x_sigma[y_][n-1]);
    } else
      fprintf(fp, "%8.3f %6.2f\n", Cell[j].S, get_code(Cell[j]));
  
  fclose(fp);
}


void config_cod(param_data_type &prms)
{
  const int n_bpm_Fam = 1, n_hcorr_Fam = 7, n_vcorr_Fam = 6;

  const std::string
    bpm_names[n_bpm_Fam] =
    {"bpm"},
    hcorr_names[n_hcorr_Fam] =
    {"sd1h", "sd2h", "sd3ah", "sd3bh", "oxx", "oxy", "oyy"},
    vcorr_names[n_vcorr_Fam] =
    {"sf1",  "sf2",  "sf3", "oxx", "oxy", "oyy"};

  prms.ini_COD_corr(n_bpm_Fam, bpm_names, n_hcorr_Fam, hcorr_names, n_vcorr_Fam,
		    vcorr_names, true);

  prt_gcmat(1); prt_gcmat(2);
}


void chk_cod_corr(const double dx_rms, const double dy_rms, const int seed,
		  const int n_seed)
{
  param_data_type prms;

  const double dxy_rms[] = {dx_rms, dy_rms};

  iniranf(seed);
  setrancut(1e0);

  config_cod(prms);

  get_cod_rms(dxy_rms, n_seed, true);
}


int main(int argc, char *argv[])
{
  const long seed = 1121;
  
  trace            = false;
  reverse_elem     = true;
  globval.mat_meth = false;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  globval.H_exact    = false; globval.quad_fringe    = false;
  globval.Cavity_on  = false; globval.radiation      = false;
  globval.emittance  = false; globval.IBS            = false;
  globval.pathlength = false; globval.Aperture_on    = false;
  globval.Cart_Bend  = false; globval.dip_edge_fudge = true;

  Ring_GetTwiss(true, 0e-3); printglob();

  GetEmittance(ElemIndex("cav"), true);

  chk_cod_corr(100e-6, 100e-6, seed, 100);
}
