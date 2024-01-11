#include <assert.h>

#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


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

  reverse_elem     = true;
  globval.mat_meth = false;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  set_state();

  globval.Cavity_on = globval.radiation = false;
  Ring_GetTwiss(true, 0e0);
  printglob();

  globval.Cavity_on = globval.radiation = true;
  Ring_GetTwiss(true, 0e0);
  printglob();
}
