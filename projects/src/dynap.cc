#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


void err_and_corr(const string &param_file)
{
  param_data_type params;
  DA_data_type    DA;

  params.get_param(param_file);

  globval.dPcommon = 1e-10;

  Ring_GetTwiss(true, 0e0);
  printglob();

  params.err_and_corr_init(param_file);

  globval.CODeps = 1e-10;

  globval.Cavity_on = true;

  if (params.DA_bare) DA.get_DA_bare(params);

  DA.get_DA_real(params);

  params.err_and_corr_exit();
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
  set_state();

  globval.mat_meth = false;

  trace = false;

  if (argc == 2)
    err_and_corr(argv[1]);
  else {
    printf("*** bad command line\n");
    exit(1);
  }
}
