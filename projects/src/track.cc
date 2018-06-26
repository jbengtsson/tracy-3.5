#define NO 1

#include "tracy_lib.h"

int no_tps = NO; // arbitrary TPSA order is defined locally


void prt_Beam_Pos(const string &file_name)
{
  int      k;
  ofstream outf;

  file_wr(outf, file_name.c_str());

  for (k = 0; k < Lattice.param.Cell_nLoc; k++)
    outf << std::scientific << std::setprecision(15)
	 << std::setw(4) << k << " " << std::setw(10) << Lattice.Cell[k]->Name
	 << std::setw(23) << Lattice.Cell[k]->BeamPos << "\n";

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
  Lattice.Cell_Pass(0, Lattice.param.Cell_nLoc, ps, lastpos);  
  std::cout << std::scientific << std::setprecision(3)
	    << std::setw(11) << ps << "\n"; 

  prt_Beam_Pos("track.out");

  Lattice.prtmfile("flat_file.dat");
}



int main(int argc, char *argv[])
{

  Lattice.param.H_exact           = false; Lattice.param.quad_fringe  = false;
  Lattice.param.Cavity_on         = false; Lattice.param.radiation    = false;
  Lattice.param.emittance         = false; Lattice.param.IBS          = false;
  Lattice.param.pathlength        = false; Lattice.param.bpm          = 0;
  Lattice.param.dip_edge_fudge    = true;  Lattice.param.reverse_elem = !false;
  // 1: DIAMOND, 3: Oleg I, 4: Oleg II.
  Lattice.param.FieldMap_filetype = 1;     Lattice.param.sympl        = false;

  if (true)
    Lattice.Read_Lattice(argv[1]); 
  else
    Lattice.rdmfile("flat_file.dat");

  trace = !true;

  Lattice.Ring_GetTwiss(true, 0e-2); printglob();

  track(1e-3, 0e0, 1e-3, 0e0, 0e0);
}
