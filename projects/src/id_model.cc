#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


template<typename T>
void prt_ps(const ss_vect<T> &ps)
{
  cout << "\n" << scientific << setprecision(3) << setw(11)
       << is_double< ss_vect<T> >::cst(ps) << "\n";
}


int main(int argc, char *argv[])
{
  long int        lastpos;
  int             loc;
  ss_vect<double> ps;

  // DIAMOND 1, NSLS-II 2, Oleg I 3, Oleg II 4, SRW 5.
  FieldMap_filetype = 5; sympl = false;

  reverse_elem = !false;

  trace = !false;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  globval.H_exact        = false;  globval.quad_fringe = false;
  globval.Cavity_on      = false;  globval.radiation   = false;
  globval.emittance      = false;  globval.IBS         = false;
  globval.pathlength     = false;  globval.bpm         = 0;
  globval.dip_edge_fudge = true;


  loc = Elem_GetPos(ElemIndex("tpw"), 1);
  ps.zero(); ps[x_] = 0e-3; ps[y_] = 1e-3;
  prt_ps(ps);
  Cell_Pass(loc, globval.Cell_nLoc, ps, lastpos);
  prt_ps(ps);

  // Ring_GetTwiss(true, 0e0); printglob();

  // GetEmittance(ElemIndex("cav"), true);

  // prt_lat("linlat1.out", globval.bpm, true);
  // prt_lat("linlat.out", globval.bpm, true, 10);
}
