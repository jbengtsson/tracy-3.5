#define NO 1

#include "tracy_lib.h"

int no_tps = NO;

// MAX-IV           1,
// SLS-2            2,
// DIAMOND          3,
// DIAMOND-II 4-BA  4,
// DIAMOND-II 6-BA  5.
// DIAMOND-II 8-BA  6.
const int lat_case = 3;

const double delta_max[] = {3.0e-2, 3.0e-2, 1.5e-2, 3e-2, 3e-2, 3e-2};

int main(int argc, char *argv[])
{

  const int    n_turn = 2064;
  const double delta  = 1.5e-2;

  Lattice.param.H_exact    = false; Lattice.param.quad_fringe = false;
  Lattice.param.Cavity_on  = false; Lattice.param.radiation   = false;
  Lattice.param.emittance  = false; Lattice.param.IBS         = false;
  Lattice.param.pathlength = false; Lattice.param.bpm         = 0;

  if (false)
    Lattice.Read_Lattice(argv[1]);
  else
    Lattice.rdmfile(argv[1]);

  Lattice.Ring_GetTwiss(true, 0e0); printglob();

  Lattice.prt_lat("linlat1.out", Lattice.param.bpm, true);
  Lattice.prt_lat("linlat.out", Lattice.param.bpm, true, 10);

  if (true) Lattice.GetEmittance(Lattice.Elem_Index("cav"), true);

  if (true) {
    Lattice.param.Cavity_on = true;
    Lattice.get_dynap(delta_max[lat_case-1], 25, n_turn, false);
  }

}
