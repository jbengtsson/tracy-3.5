#define NO 1

#include <assert.h>

#include "tracy_lib.h"


int no_tps = NO;


void prt_kick_map(CellType &Cell)
{
  for (int j = 0; j < Cell.ID->nx; j++) {
    for (int k = 0; k < Cell.ID->ny; k++)
      printf("", Cell.ID->tabx[j], Cell.ID->tabz[k], thetax[j][k]);
    printf("\n");
  }
}


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


int main(int argc, char *argv[])
{
  globval.mat_meth = false;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  set_state();

  Ring_GetTwiss(true, 0e0);
  printglob();

  auto loc = Elem_GetPos(ElemIndex("cpmu", 1)
  prt_kick_map(Cell[loc]);
}
