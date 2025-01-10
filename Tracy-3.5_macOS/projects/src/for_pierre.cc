#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


void track_one_turn(const double A_x, const double A_y, const double delta)
{
  long int        lastpos;
  ss_vect<double> ps;
  
  globval.Cavity_on = false; globval.radiation = false;
  globval.pathlength = false;

  cout  << "\ntrack_one_turn:\n";
  ps.zero();
  ps[x_] = A_x; ps[y_] = A_y; ps[delta_] = delta;
  cout << scientific << setprecision(5) << setw(13) << ps << "\n";
  Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
  cout << scientific << setprecision(5) << setw(13) << ps
       << " lastpos = " << lastpos << "(" << globval.Cell_nLoc << ")" << "\n";
}


ss_vect<tps> get_map()
{
  long int     lastpos;
  ss_vect<tps> map;
  
  globval.Cavity_on = false; globval.radiation = false;
  globval.pathlength = false;

  map.identity();
  Cell_Pass(0, globval.Cell_nLoc, map, lastpos);
  prt_lin_map(3, map);
  return map;
}


int main(int argc, char *argv[])
{

  reverse_elem = true;

  trace = false;

  globval.mat_meth = !false;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  globval.H_exact    = false; globval.quad_fringe    = false;
  globval.Cavity_on  = false; globval.radiation      = false;
  globval.emittance  = false; globval.IBS            = false;
  globval.pathlength = false; globval.Aperture_on    = false;
  globval.Cart_Bend  = false; globval.dip_edge_fudge = true;

  // Check lattice.
  Ring_GetTwiss(true, 0e0); printglob();

  track_one_turn(1e-3, -1e-3, 0e0);

  // "Jelly Bean" (a treat for you, for fun):
  // How to compute the Poincaré/one-turn map.
  get_map();
}
