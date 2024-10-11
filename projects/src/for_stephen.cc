#define NO 1

#include <cstdio>
#include <assert.h>

#include "tracy_lib.h"


int no_tps = NO;


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
}


void compute_mat(void)
{
  const string file_name = "for_stephen.txt";

  long int     lastpos;
  ss_vect<tps> M;

  // Redirect stdout to text file.
  freopen(file_name.c_str(), "a+", stdout); 

  M.identity();
  for (auto k = 0; k <= globval.Cell_nLoc; k++) {
    Cell_Pass(k, k, M, lastpos);
    printf("\n  %2d %10s %7.3f", k, Cell[k].Elem.PName, Cell[k].S);
    prt_lin_map(3, M);
  }

  // Restore stdout.
  freopen("/dev/tty", "w", stdout);
}


int main(int argc, char *argv[])
{
  globval.mat_meth = !false;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  set_state();

  Ring_GetTwiss(true, 0e0);
  printglob();

  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);

  compute_mat();
}
