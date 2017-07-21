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

  
  loc = Elem_GetPos(Fnum, 1); FM = Cell[loc].Elem.FM;

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
    FM->L1 += Cell[loc].Elem.PL - FM->Lm;
  }

  map.identity(); Cell_Pass(loc, loc, map, lastpos);

  cout << endl;
  cout << scientific << setprecision(3)
       << "get_fix_pont:" << setw(11) << map.cst();
  cout << fixed << setprecision(5)
       << "phi [deg] = " << setw(7) << FM->phi*180.0/M_PI
       << ", L [m] = " << setw(7) << Cell[loc].Elem.PL << endl;
  cout << fixed << setprecision(5)
       << "Lr [m] = " << setw(7) << FM->Lr
       << ", Lm [m] = " << setw(7) << FM->Lm
       << ", Ld [m] = " << setw(7) << FM->Ld
       << ", L1 [m] = " << setw(7) << FM->L1 << endl;

  for (j = 1; j <= GetnKid(Fnum); j++) {
    loc = Elem_GetPos(Fnum, j);
    Cell[loc].Elem.FM->cut = FM->cut;
    Cell[loc].Elem.FM->x0 = FM->x0;
    Cell[loc].Elem.FM->phi = FM->phi;
    Cell[loc].Elem.FM->Lr = FM->Lr; Cell[loc].Elem.FM->Lm = FM->Lm;
    Cell[loc].Elem.FM->Ld = FM->Ld; Cell[loc].Elem.FM->L1 = FM->L1;
  }

  return map;
}


int main(int argc, char *argv[])
{
  long int     lastpos;
  double       dx;
  Matrix       M;
  ss_vect<tps> map;

    
  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  // 1: DIAMOND, 3: Oleg I, 4: Oleg II.
  FieldMap_filetype = 1; sympl = false;

  Read_Lattice(argv[1]);

  // no_sxt();

//  Ring_GetTwiss(true, 0.0); printglob();

  if (true) {
    trace = true;

    map.identity();
    // Tweak to remain within field map range at entrance.
    dx = -1.3e-3;
    map[x_] += dx;
      Cell_Pass(Elem_GetPos(ElemIndex("bb"), 1),
		Elem_GetPos(ElemIndex("bb"), 1), map, lastpos);
    // Tweak to remain within field map range at entrance.
    map[x_] -= dx;

    getlinmat(6, map, M);
    prt_lin_map(3, map);
    cout << "\n" << scientific << setprecision(3)
	 << setw(11) << map.cst() << "\n";
    cout << scientific << setprecision(3)
	  << endl << "1-Det: " << setw(9) << 1-DetMat(6, M) << endl;

    exit(0);
  }

  if (false) {
    map = get_fix_point(ElemIndex("bb"));

    prt_lin_map(3, map);

    getlinmat(6, map, M);
    cout << endl;
    cout << scientific << setprecision(3)
	 << "1-Det: " << setw(9) << 1-DetMat(6, M) << endl;
    exit(0);
  }

  if (true) {
    trace = true;

    Ring_GetTwiss(true, 0.0); printglob();

    prt_lat("linlat.out", globval.bpm, true);  
  }
}
