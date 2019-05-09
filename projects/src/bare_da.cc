#define NO 1

#include "tracy_lib.h"

int no_tps = NO;

// MAX-IV   1,
// SLS-2    2,
// M-H6BAi  3.
// M-H6BA   4.
// M-H8BA   5.
// RB-6BA   6.
// M-6BA    7.
// DIAMOND  8,
// DELTA    9.
// ALS-U    10.

const bool set_dnu  = false;

const double
  delta_max = 2.5e-2,
  dnu[]     = {0.03, 0.02};


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

  if (set_dnu) {
    Ring_GetTwiss(true, 0e0); printglob();
    set_map(ElemIndex("ps_rot"), dnu);
    Ring_GetTwiss(true, 0e0); printglob();
  }

  if (true) {
    globval.Cavity_on = true;
    get_dynap(delta_max, 25, n_turn, false);
  }

}
