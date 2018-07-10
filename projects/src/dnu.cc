#define NO 1

#include "tracy_lib.h"

int  no_tps = NO;


// MAX-IV              1,
// SLS-2               2,
// DIAMOND             3,
// DIAMOND-II DTBA     4,
// DIAMOND-II H-6-BA   5.
// DIAMOND-II H-8-BA   6.
// DIAMOND-II RB-6-BA  7.
// DIAMOND-II D-DBA    8.
// DELTA               9.
// ALS-U              10.

const int lat_case = 2;

const double
  A_max[][2] =
    {{1.5e-3, 1.5e-3}, {7e-3, 5e-3}, {15e-3, 8e-3}, {10e-3, 4e-3},
     {  7e-3,   4e-3}, {3e-3, 2e-3},  {3e-3, 2e-3}, {3e-3, 2e-3},
      {35e-3,   6e-3}, {4e-3, 4e-3}},
  delta_max[] = {3e-2, 5e-2, 1.5e-2, 3e-2, 3e-2, 3e-2, 3e-2, 2e-2, 3e-2, 4e-2},
  dnu[]       = {0.1, 0.0};


int main(int argc, char *argv[])
{
  
  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  if (false)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  globval.EPU = false;

  if (false) {
    Ring_GetTwiss(true, 0e0); printglob();
    set_map("ps_rot", dnu[X_], dnu[Y_]);
    Ring_GetTwiss(true, 0e0); printglob();
  }

  dnu_dA(A_max[lat_case-1][X_], A_max[lat_case-1][Y_], 0e0, 25);
  get_ksi2(delta_max[lat_case-1]);
}
