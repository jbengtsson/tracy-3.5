#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


void compute_map(void)
{
  long int lastpos;
  ss_vect<tps> map;

  map.identity();
  Cell_Pass(0, globval.Cell_nLoc, map, lastpos);
  prt_lin_map(3, map);
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
  globval.EPU            = !true;
}


int main(int argc, char *argv[])
{

  globval.mat_meth = false;

  FieldMap_filetype = 6;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  set_state();

  if (false) {
    compute_map();
    assert(false);
  }
  

  Ring_GetTwiss(true, 0e0);
  printglob();

  prtmfile("flat_file.dat");
  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);
}
