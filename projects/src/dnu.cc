#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


const double
#if 1
  A_max[]   = {2e-3, 2e-3},
  delta_max = 3e-2;
#else
  A_max[]   = {8e-3, 8e-3},
  delta_max = 4e-2;
#endif


int main(int argc, char *argv[])
{

  globval.H_exact    = false; globval.quad_fringe    = false;
  globval.Cavity_on  = false; globval.radiation      = false;
  globval.emittance  = false; globval.IBS            = false;
  globval.pathlength = false; globval.bpm            = 0;
  globval.Cart_Bend  = false; globval.dip_edge_fudge = true;

  reverse_elem = !false;

  globval.mat_meth = false;

  if (false)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  if (false) {
    Ring_GetTwiss(true, 0e0); printglob();
  }

  dnu_dA(A_max[X_], A_max[Y_], 0e0, 25);
  get_ksi2(delta_max);
}
