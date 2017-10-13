#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


const bool NSLS_II = false;


ss_vect<tps> get_fix_point(const int Fnum)
{
  long int     lastpos, loc;
  int          j;
  FieldMapType *FM;
  
  const int  n_iter = 3;

  
  loc = Lattice.Elem_GetPos(Fnum, 1); FM = Lattice.Cell[loc].Elem.FM;

  FM->Ld = 0.0; FM->L1 = 0.0;
  if (!NSLS_II) {
    // DIAMOND
    FM->cut = 0; FM->x0 = -30e-3;
  } else {
    // NSLS-II
    FM->cut = 50; FM->x0 = 55e-3;
  }

  trace = true;
  for (j = 1; j <= n_iter; j++) {
    map.identity(); Cell_Pass(loc, loc, map, lastpos);

    FM->phi -= map[px_].cst();
    FM->Lm = map[ct_].cst();
    FM->Ld += 2.0*map[x_].cst()/FM->phi;
    FM->L1 += Lattice.Cell[loc].L - FM->Lm;
  }

  map.identity(); Cell_Pass(loc, loc, map, lastpos);

  cout << endl;
  cout << scientific << setprecision(3)
       << "get_fix_pont:" << setw(11) << map.cst();
  cout << fixed << setprecision(5)
       << "phi [deg] = " << setw(7) << FM->phi*180.0/M_PI
       << ", L [m] = " << setw(7) << Lattice.Cell[loc].L << endl;
  cout << fixed << setprecision(5)
       << "Lr [m] = " << setw(7) << FM->Lr
       << ", Lm [m] = " << setw(7) << FM->Lm
       << ", Ld [m] = " << setw(7) << FM->Ld
       << ", L1 [m] = " << setw(7) << FM->L1 << endl;

  for (j = 1; j <= Lattice.GetnKid(Fnum); j++) {
    loc = Lattice.Elem_GetPos(Fnum, j);
    Lattice.Cell[loc].Elem.FM->cut = FM->cut;
    Lattice.Cell[loc].Elem.FM->x0 = FM->x0;
    Lattice.Cell[loc].Elem.FM->phi = FM->phi;
    Lattice.Cell[loc].Elem.FM->Lr = FM->Lr;
    Lattice.Cell[loc].Elem.FM->Lm = FM->Lm;
    Lattice.Cell[loc].Elem.FM->Ld = FM->Ld;
    Lattice.Cell[loc].Elem.FM->L1 = FM->L1;
  }

  return map;
}


int main(int argc, char *argv[])
{
  bool         tweak;
  long int     lastpos;
  double       dx;
  Matrix       M;
  ss_vect<tps> map;

    
  Lattice.param.H_exact    = false; Lattice.param.quad_fringe = false;
  Lattice.param.Cavity_on  = false; Lattice.param.radiation   = false;
  Lattice.param.emittance  = false; Lattice.param.IBS         = false;
  Lattice.param.pathlength = false; Lattice.param.bpm         = 0;

  // 1: DIAMOND, 3: Oleg I, 4: Oleg II.
  FieldMap_filetype = 1; sympl = false;

  Lattice.Read_Lattice(argv[1]);

  // no_sxt();

//  Ring_GetTwiss(true, 0.0); printglob();

  if (true) {
    trace = false;

    Lattice.param.H_exact    = true;
    Lattice.param.dip_fringe = true;

    map.identity();
    // Tweak to remain within field map range at entrance.
    tweak = true;
    if (tweak) {
      dx = -1.4e-3; map[x_] += dx;
    }
    Cell_Pass(Lattice.Elem_GetPos(Lattice.Elem_Index("bb"), 1),
	      Lattice.Elem_GetPos(Lattice.Elem_Index("bb"), 1), map, lastpos);
    if (tweak) map[x_] -= dx;

    getlinmat(6, map, M);
    prt_lin_map(3, map);
    cout << "\n" << scientific << setprecision(3)
	 << setw(11) << map.cst() << "\n";
    cout << scientific << setprecision(3)
	  << endl << "1-Det: " << setw(9) << 1-DetMat(6, M) << endl;

    exit(0);
  }

  if (false) {
    map = get_fix_point(Lattice.Elem_Index("bb"));

    prt_lin_map(3, map);

    getlinmat(6, map, M);
    cout << endl;
    cout << scientific << setprecision(3)
	 << "1-Det: " << setw(9) << 1-DetMat(6, M) << endl;
    exit(0);
  }

  if (true) {
    trace = true;

    Lattice.Ring_GetTwiss(true, 0.0); printglob();

    Lattice.prt_lat("linlat.out", Lattice.param.bpm, true);  
  }
}
