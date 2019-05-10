#define NO 1

#include "tracy_lib.h"

int  no_tps = NO;


const int  lat_case = 3;
const double
  A_max[]   = {6e-3, 2e-3},
  delta_max = 2.5e-2;
  // ALS-U.
  // A_max[]   = {4e-3, 2.5e-3},
  // delta_max = 4e-2;


int main(int argc, char *argv[])
{

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  reverse_elem = !false;

  if (false)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  globval.EPU = false;

  dnu_dA(A_max[X_], A_max[Y_], 0e0, 25);
  get_ksi2(delta_max);
}
