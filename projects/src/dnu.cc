#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


const int
  n_step = 25;

const double
#if 1
  A_max[]   = {4e-3, 4e-3},
  delta_max = 6e-2;
#else
  A_max[]   = {6e-3, 2.5e-3},
  delta_max = 3.5e-2;
#endif


void set_state(void)
{
  globval.H_exact        = false;
  globval.quad_fringe    = false;
  globval.Cavity_on      = false;
  globval.radiation      = false;
  globval.emittance      = false;
  globval.IBS            = false;
  globval.pathlength     = false;
  globval.Aperture_on    = false;
  globval.Cart_Bend      = false;
  globval.dip_edge_fudge = true;
  globval.EPU            = !false;
}


int main(int argc, char *argv[])
{
  reverse_elem = !false;

  globval.mat_meth = false;

  trace = false;

  if (false)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  set_state();

  // globval.CODeps   = 1e-10;
  // globval.dPcommon = 1e-6;

  if (!false) {
    Ring_GetTwiss(true, 0e0); printglob();
    prt_lat("linlat1.out", globval.bpm, true);
    prt_lat("linlat.out", globval.bpm, true, 10);
  }

  dnu_dA(A_max[X_], A_max[Y_], 0e0, n_step);
  get_ksi2(delta_max, n_step);
}
