#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


const bool   set_dnu = false;
const double dnu[]   = {0.03, 0.02};


void err_and_corr(const string &param_file)
{
  param_data_type params;
  orb_corr_type   orb_corr[2];
  DA_data_type    DA;

  params.get_param(param_file);

  globval.dPcommon = 1e-10;

  Ring_GetTwiss(true, 0e0); printglob();

  if (set_dnu) {
    // Do not use for Real Lattice.
    set_map(ElemIndex("ps_rot"), dnu);
    Ring_GetTwiss(true, 0e0); printglob();
  }

  params.err_and_corr_init(param_file, orb_corr);

  globval.CODeps = 1e-10;

  globval.Cavity_on = true;

  if (params.DA_bare) DA.get_DA_bare(params);

  DA.get_DA_real(params, orb_corr);

  params.err_and_corr_exit(orb_corr);
}


int main(int argc, char *argv[])
{
  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.Aperture_on = false;

  if (argc == 2)
    err_and_corr(argv[1]);
  else {
    printf("*** bad command line\n");
    exit(1);
  }
}
