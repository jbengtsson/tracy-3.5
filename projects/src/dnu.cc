#define NO 1

#include "tracy_lib.h"

int  no_tps = NO;


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
const int  lat_case = 3;
const double
  A_max[][2]  =
    {{1.5e-3, 1.5e-3}, { 7e-3, 5e-3}, {6e-3, 3e-3}, {10e-3, 4e-3},
     {  5e-3,   3e-3}, { 3e-3, 2e-3},  {3e-3, 2e-3}, {3e-3, 2e-3},
     {  3e-3,   2e-3}, {35e-3, 6e-3}, {4e-3, 4e-3}},
  delta_max[] = {3e-2, 5e-2, 2e-2, 3e-2,
		 3e-2, 3e-2, 3e-2, 2e-2,
		 3e-2, 3e-2, 4e-2},
  dnu[]       = {0.03, 0.02},
  eta_x[]     = {0.0,  0.0};


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

  if (set_dnu) {
    Ring_GetTwiss(true, 0e0); printglob();
    set_map("ps_rot", dnu[X_], dnu[Y_], eta_x[0], eta_x[1]);
    Ring_GetTwiss(true, 0e0); printglob();
  }

  dnu_dA(A_max[lat_case-1][X_], A_max[lat_case-1][Y_], 0e0, 25);
  get_ksi2(delta_max[lat_case-1]);
}
