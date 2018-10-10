#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


int main(int argc, char *argv[])
{

  // 1: DIAMOND, 2: NSLS-II, 3: Oleg I, 4: Oleg II, 5: SRW.
  FieldMap_filetype = 5; sympl = !false;

  reverse_elem = !false;

  trace = false;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  globval.H_exact        = false;  globval.quad_fringe = false;
  globval.Cavity_on      = false;  globval.radiation   = false;
  globval.emittance      = false;  globval.IBS         = false;
  globval.pathlength     = false;  globval.bpm         = 0;
  globval.dip_edge_fudge = true;

  Ring_GetTwiss(true, 0e0); printglob();

  GetEmittance(ElemIndex("cav"), true);
}
