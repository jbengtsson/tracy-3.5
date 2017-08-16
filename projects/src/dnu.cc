#define NO 1

#include "tracy_lib.h"

int  no_tps = NO;


// MAX-IV           1,
// SLS-2            2,
// DIAMOND          3,
// DIAMOND-II 4-BA  4,
// DIAMOND-II 6-BA  5.
const int lat_case = 3;

const double
  A_max[][2] =
    {{1.5e-3, 1.5e-3}, {7.0e-3, 5.0e-3}, {15.0e-3, 8e-3}, {5.0e-3, 3.0e-3},
     {6.0e-3, 4.0e-3}},
  delta_max[] = {3.0e-2, 3.0e-2, 1.5e-2, 3e-2, 3e-2};


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

  dnu_dA(A_max[lat_case-1][X_], A_max[lat_case-1][Y_], 0e0, 25);
  get_ksi2(delta_max[lat_case-1]);
}
