#define NO 1

#include "tracy_lib.h"

int no_tps = NO; // arbitrary TPSA order is defined locally


void track(const double x, const double p_x,
	   const double y, const double p_y, const double delta)
{
  long int        lastpos;
  ss_vect<double> ps;

  ps.zero();
  ps[x_] = x; ps[px_] = p_x; ps[y_] = y; ps[py_] = p_y; ps[delta_] = delta;

  std::cout << std::scientific << std::setprecision(3)
	    << std::setw(11) << ps << "\n"; 
  Lattice.Cell_Pass(0, Lattice.param.Cell_nLoc, ps, lastpos);  
  std::cout << std::scientific << std::setprecision(3)
	    << std::setw(11) << ps << "\n"; 
}



int main(int argc, char *argv[])
{

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  reverse_elem = !false;

  if (true)
    Lattice.Read_Lattice(argv[1]); 
  else
    Lattice.rdmfile("flat_file.dat");

  trace = true;

  track(1e-3, 0e0, 1e-3, 0e0, 0e0);

  Lattice.Ring_GetTwiss(true, 0e-2); printglob();
}
