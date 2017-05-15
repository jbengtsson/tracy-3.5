#define NO 1

#include "tracy_lib.h"

int  no_tps = NO;


int main(int argc, char *argv[])
{
  
  // const double A_max[] = {1.5e-3, 1.5e-3}, delta_max = 4.5e-2;
  const double A_max[] = {6.0e-3, 6.0e-3}, delta_max = 5.0e-2;

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  if (false)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  globval.EPU = false;

  dnu_dA(A_max[X_], A_max[Y_], 0e0, 25); get_ksi2(delta_max);
}
