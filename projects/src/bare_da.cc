#define NO 1

#include "tracy_lib.h"

int no_tps = NO;

// MAX-IV              1,
// SLS-2               2,
// DIAMOND             3,
// DIAMOND-II 4-BA     4,
// DIAMOND-II 6-BA     5.
// DIAMOND-II 8-BA     6.
// DIAMOND-II RB-6-BA  7.
// DIAMOND-II D-TBA    8.
// DELTA               9.
// ALS-U              10.

const int lat_case = 5;

const double
  delta_max[] =
    {3e-2, 3e-2, 1.5e-2, 3e-2, 3e-2, 1e-2, 3e-2, 2e-2, 1e-2, 3e-2},
  dnu[]       = {0.1, 0.0};


int main(int argc, char *argv[])
{

  const int    n_turn = 2064;

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  if (false)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  Ring_GetTwiss(true, 0e0); printglob();

  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);

  if (true) GetEmittance(ElemIndex("cav"), true);

  if (false) {
    Ring_GetTwiss(true, 0e0); printglob();
    set_map("ps_rot", dnu[X_], dnu[Y_]);
    Ring_GetTwiss(true, 0e0); printglob();
  }

  if (true) {
    globval.Cavity_on = true;
    get_dynap(delta_max[lat_case-1], 25, n_turn, false);
  }

}
