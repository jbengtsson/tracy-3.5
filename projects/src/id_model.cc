#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


template<typename T>
void prt_ps(const ss_vect<T> &ps)
{
  cout << "\n" << scientific << setprecision(3) << setw(11)
       << is_double< ss_vect<T> >::cst(ps) << "\n";
}


void prt_lin_map1(const int n_DOF, const ss_vect<tps> &map)
{
  int i, j;

  std::cout << std::endl;
  for (i = 1; i <= 2*n_DOF; i++) {
    for (j = 1; j <= 2*n_DOF; j++)
      if (true)
	std::cout << std::scientific << std::setprecision(10)
	     << std::setw(18) << getmat(map, i, j);
      else
	std::cout << std::scientific << std::setprecision(16)
	     << std::setw(24) << getmat(map, i, j);
    std::cout << std::endl;
  }
}


int main(int argc, char *argv[])
{
  long int        lastpos;
  int             loc1, loc2;
  ss_vect<double> ps;
  ss_vect<tps>    map;
  Matrix          M;

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


  if (false) {
    globval.radiation = true;

    loc1 = Elem_GetPos(ElemIndex("d_dw"), 1);
    loc2 = Elem_GetPos(ElemIndex("d_dw"), 2);
    // loc1 = Elem_GetPos(ElemIndex("tpw"), 1);
    // loc2 = Elem_GetPos(ElemIndex("tpw"), 1);
    ps.zero();
    prt_ps(ps);
    Cell_Pass(loc1, loc2, ps, lastpos);
    prt_ps(ps);
    exit(0);
  }

  if (false) {
    globval.radiation = !true;

    loc1 = Elem_GetPos(ElemIndex("d_dw"), 1);
    loc2 = Elem_GetPos(ElemIndex("d_dw"), 2);
    // loc1 = Elem_GetPos(ElemIndex("tpw"), 1);
    // loc2 = Elem_GetPos(ElemIndex("tpw"), 1);
    map.identity();
    Cell_Pass(loc1, loc2, map, lastpos);
    prt_lin_map1(3, map);

    if (!false) cout << scientific << setprecision(3) << setw(10) << map;

    getlinmat(6, map, M);
    printf("\n  1-Det{} = %9.3e\n", 1e0-DetMat(6, M));

    exit(0);
  }

  Ring_GetTwiss(true, 0e0); printglob();

  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);

  GetEmittance(ElemIndex("cav"), false, true);
}
