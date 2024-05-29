#define NO 1

#include "tracy_lib.h"


int no_tps = NO;


void set_ps_rot(const double dnu_x, const double dnu_y)
{
  const double
    dnu_0[] = {0e0, 0e0},
    dnu[]   = {dnu_x, dnu_y};

  set_map(ElemIndex("ps_rot"), dnu_0);
  printf("\ntweak_nu:\n");
  printf("  dnu = [%8.5f, %8.5f]\n", dnu[X_], dnu[Y_]);
  set_map(ElemIndex("ps_rot"), dnu);
  Ring_GetTwiss(true, 0e0);
  printglob();
 }


void get_dnu_straight(const int loc_1, const int loc_2)
{
  const double
    dnu[] = {
      Cell[globval.Cell_nLoc].Nu[X_]-Cell[loc_2].Nu[X_]+Cell[loc_1].Nu[X_],
      Cell[globval.Cell_nLoc].Nu[Y_]-Cell[loc_2].Nu[Y_]+Cell[loc_1].Nu[Y_]
    };

  printf("\nget_dnu_straight:\n");
  printf(" dnu = [%7.5f, %7.5f]\n", dnu[X_], dnu[Y_]);
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
}


int main(int argc, char *argv[])
{
  int loc_1, loc_2;

  globval.mat_meth = false;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  set_state();

  if (false) {
    no_mult(Sext);
    no_mult(Oct);
  }

  Ring_GetTwiss(true, 0e0);
  printglob();

  if (false)
    set_ps_rot(-0.0351-0.29, 0.35);

  prtmfile("flat_file.dat");
  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);
  prt_chrom_lat("chromlat.out");

  if (!false)
    GetEmittance(ElemIndex("cav"), false, true);

  if (!false) {
    loc_1 = Elem_GetPos(ElemIndex("sd2"), 1);
    loc_2 = Elem_GetPos(ElemIndex("sd2"), 2);
    get_dnu_straight(loc_1, loc_2);
  }
}
