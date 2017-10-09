#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


void err_and_corr(const string &param_file)
{
  param_data_type params;
  orb_corr_type   orb_corr[2];
  DA_data_type    DA;

  params.err_and_corr_init(param_file, orb_corr);

  Lattice.param.CODeps = 1e-10;

  Lattice.param.Cavity_on = true;

  if (params.DA_bare) DA.get_DA_bare(params);

  DA.get_DA_real(params, orb_corr);

  params.err_and_corr_exit(orb_corr);
}


int main(int argc, char *argv[])
{
  Lattice.param.H_exact    = false; Lattice.param.quad_fringe = false;
  Lattice.param.Cavity_on  = false; Lattice.param.radiation   = false;
  Lattice.param.emittance  = false; Lattice.param.IBS         = false;
  Lattice.param.pathlength = false; Lattice.param.Aperture_on = false;

  if (argc < 1) {
    printf("*** bad command line\n");
    exit(1);
  }

  err_and_corr(argv[1]);
}
