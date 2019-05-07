#define NO 1

#include "tracy_lib.h"

int  no_tps = NO;


const bool set_dnu  = !false;
const int  lat_case = 3;
const double
  A_max[]   = {6e-3, 2e-3},
  delta_max = 2.5e-2,
  // ALS-U.
  // A_max[]   = {4e-3, 2.5e-3},
  // delta_max = 4e-2,
  dnu[]     = {-0.0/6.0, -0.02/6.0};


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

  if (set_dnu) {
    Ring_GetTwiss(true, 0e0); printglob();
    set_map(ElemIndex("ps_rot"), dnu);
    Ring_GetTwiss(true, 0e0); printglob();
  }

  dnu_dA(A_max[X_], A_max[Y_], 0e0, 25);
  get_ksi2(delta_max);
}
