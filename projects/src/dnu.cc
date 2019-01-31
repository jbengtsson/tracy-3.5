#define NO 1

#include "tracy_lib.h"

int  no_tps = NO;


// MAX-IV       1,
// SLS-2        2,
// M-H6BAi      3,
// M-H6BA-0-.-. 4,
// DIAMOND      5,
// DELTA        6,
// ALS-U        7.

const bool set_dnu  = !false;
const int  lat_case = 4;
const double
A_max[][2]  =
  {{1.5e-3, 1.5e-3}, { 7e-3, 5e-3}, {6e-3, 3e-3}, {8e-3, 5e-3},
   {  5e-3,   3e-3}, { 3e-3, 2e-3}, {3e-3, 2e-3}},
  delta_max[] = {3e-2, 5e-2, 2e-2, 3e-2,
		 3e-2, 3e-2, 3e-2},
  dnu[]       = {0.05/6.0, 0.0};


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

  dnu_dA(A_max[lat_case-1][X_], A_max[lat_case-1][Y_], 0e0, 25);
  get_ksi2(delta_max[lat_case-1]);
}
