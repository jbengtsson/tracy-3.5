#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


void err_and_corr(const string &param_file, const int mode)
{
  bool            cod;
  long int        lastpos;
  double          m_dbeta[2], s_dbeta[2], m_dnu[2], s_dnu[2];
  param_data_type params;
  orb_corr_type   orb_corr[2];

  params.err_and_corr_init(param_file, orb_corr);

  globval.CODeps = 1e-10;

  if (params.fe_file != "") params.LoadFieldErr(false, 1e0, true);
  if (params.ae_file != "") {
    // Load misalignments; set seed, no scaling of rms errors.
    params.LoadAlignTol(false, 1e0, true, 1);
    // Beam based alignment.
    if (params.bba) params.Align_BPMs(Quad);

    cod = params.cod_corr(params.n_cell, 1e0, orb_corr);
  } else
    cod = getcod(0e0, lastpos);

  params.Orb_and_Trim_Stat();

  if (params.N_calls > 0) {
    params.ID_corr(params.N_calls, params.N_steps, false);
    cod = params.cod_corr(params.n_cell, 1e0, orb_corr);
  }

  prtmfile("flat_file.dat");

  if (cod) {
    if (mode == 1) {
      printf("fmap\n");
      fmap(params.n_x, params.n_y, params.n_tr, params.x_max_FMA,
	   params.y_max_FMA, 0e0, true, false);
    } else if (mode == 2) {
      printf("fmapdp\n");
      fmapdp(params.n_x, params.n_dp, params.n_tr, params.x_max_FMA,
	     -params.delta_FMA, 0.1e-3, true, false);
    } else {
      printf("err_and_corr: bad param. mode = %d\n", mode);
      exit(1);
    }
  } else {
    printf("err_and_corr: orbit correction failed\n");
    exit(1);
  }

  params.err_and_corr_exit(orb_corr);
}


int main(int argc, char *argv[])
{
  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.Aperture_on = false;

  if (argc < 2) {
    printf("*** bad command line\n");
    exit(1);
  }

  err_and_corr(argv[1], atoi(argv[2]));
}
