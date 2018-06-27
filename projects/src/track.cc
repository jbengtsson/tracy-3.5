#define NO 1

#include "tracy_lib.h"

int no_tps = NO; // arbitrary TPSA order is defined locally


void prt_Beam_Pos(const string &file_name)
{
  int      k;
  ofstream outf;

  file_wr(outf, file_name.c_str());

  for (k = 0; k < globval.Cell_nLoc; k++)
    outf << std::scientific << std::setprecision(15)
	 << std::setw(4) << k << " " << std::setw(10) << Cell[k].Elem.PName
	 << std::setw(23) << Cell[k].BeamPos << "\n";

  outf.close();
}


void track(const double x, const double p_x,
	   const double y, const double p_y, const double delta)
{
  long int        lastpos;
  ss_vect<double> ps;

  ps.zero();
  ps[x_] = x; ps[px_] = p_x; ps[y_] = y; ps[py_] = p_y; ps[delta_] = delta;

  std::cout << std::scientific << std::setprecision(3)
	    << "\n" << std::setw(11) << ps << "\n"; 
  Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);  
  std::cout << std::scientific << std::setprecision(3)
	    << std::setw(11) << ps << "\n"; 

  prt_Beam_Pos("track.out");
}



int main(int argc, char *argv[])
{

  reverse_elem = !false;

  trace = !true;

  if (true)
    Read_Lattice(argv[1]); 
  else
    rdmfile("flat_file.dat");

  globval.H_exact        = false;  globval.quad_fringe = false;
  globval.Cavity_on      = false;  globval.radiation   = false;
  globval.emittance      = false;  globval.IBS         = false;
  globval.pathlength     = false;  globval.bpm         = 0;
  globval.dip_edge_fudge = true;

  Ring_GetTwiss(true, 0e-2); printglob();

  track(1e-3, 0e0, 1e-3, 0e0, 0e0);
}
