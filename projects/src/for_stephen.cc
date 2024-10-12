#define NO 1

#include <cstdio>
#include <assert.h>

#include "tracy_lib.h"


int no_tps = NO;


void chk_sympl(void)
{
  int long     lastpos;
  ss_vect<tps> M, S, diff;

  M.identity();
  Cell_Pass(0, globval.Cell_nLoc, M, lastpos);

  S = get_S(3);
  diff = tp_map(3, M)*S*M;
  printf("\ndiff:");
  prt_lin_map(3, diff);
}


void compute_mat(void)
{
  const bool   incremental = false;
  const string file_name   = "for_stephen.txt";

  long int     lastpos;
  ss_vect<tps> M;

  // Redirect stdout to text file.
  freopen(file_name.c_str(), "w", stdout); 

  M.identity();
  for (auto k = 0; k <= globval.Cell_nLoc; k++) {
    if (!incremental)
      M.identity();
    Cell_Pass(k, k, M, lastpos);
    printf("\n  %2d %10s S = %7.3f", k, Cell[k].Elem.PName, Cell[k].S);
    if (Cell[k].Elem.Pkind == Mpole)
      printf("  h = %21.16e", Cell[k].Elem.M->Pirho);
    prt_lin_map(3, M);
  }

  // Restore stdout.
  freopen("/dev/tty", "w", stdout);
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
  double I[6], eps_x, sigma_delta, U_0, J[3], tau[3];

  globval.mat_meth = !false;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  set_state();

  Ring_GetTwiss(true, 0e0);
  printglob();

  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);

  if (!false) {
    if (!globval.mat_meth)
      GetEmittance(ElemIndex("cav"), false, true);
    else
      get_eps_x(eps_x, sigma_delta, U_0, J, tau, I, true);
  }

  if (!false)
    compute_mat();

  if (false)
    chk_sympl();
}
