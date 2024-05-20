#include <assert.h>

#define NO 1

#include "tracy_lib.h"


int no_tps = NO;


double get_f_RF(const int Fnum)
{
  const int loc = Elem_GetPos(Fnum, 1);
  return Cell[loc].Elem.C->f_RF;
}


void set_f_RF(const int Fnum, const double f_RF)
{
  const int loc = Elem_GetPos(Fnum, 1);
  Cell[loc].Elem.C->f_RF = f_RF;
}


void scan_f_RF(const string &Fam_name, const double df_RF_max, const int n_step)
{
  const int
    Fnum = ElemIndex(Fam_name.c_str());
  const double
    f_RF_0    = get_f_RF(Fnum),
    f_RF_step = df_RF_max/n_step;

  long int
    lastpos;
  double
    df_RF = 0e0;

  globval.Cavity_on = globval.radiation = globval.pathlength = true;

  printf("\n   f_Rf   eps_x    delta\n");
  printf("  [kHz]  [pm.rad]   [%%]\n");
  for (int k = 0; k <= n_step; k++) {
    GetEmittance(ElemIndex("cav"), true, false);
    printf("  %5.3f   %5.1f   %5.2f\n",
	   1e-3*df_RF, 1e12*globval.eps[X_], 1e2*globval.CODvect[delta_]);
    df_RF += f_RF_step;
    set_f_RF(Fnum, f_RF_0+df_RF);
  }

  prt_cod("cod.out", 0, true);
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

  globval.mat_meth       = false;
}


int main(int argc, char *argv[])
{
  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  set_state();

  Ring_GetTwiss(true, 0e-3);
  printglob();

  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);

  GetEmittance(ElemIndex("cav"), true, true);

  scan_f_RF("cav", 500.0, 10);
}
