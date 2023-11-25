#include <assert.h>

#define NO 1

#include "tracy_lib.h"

#include "prt_ZAP.cc"

#define PM 1
#if PM
#include "PoincareMap.cc"
#else
#include "get_Poincare_Map.cc"
#endif

#include "prt_lat_param.cc"

int no_tps = NO;


void tst_D(void)
{
  long int        lastpos;
  ss_vect<double> ps, cod;
  ss_vect<tps>    D, A, M;

  globval.rad_D = globval.Cavity_on = globval.radiation = true;

  getcod(0e0, lastpos);
  cod = globval.CODvect;

  ps.zero();
  ps += cod;
  globval.Diff_mat.zero();
  Cell_Pass(0, globval.Cell_nLoc-1, ps, lastpos);

  cout << scientific << setprecision(3)
       << "\ncod:\n" << setw(11) << cod << "\n";
  cout << scientific << setprecision(3)
       << "ps:\n" << setw(11) << ps << "\n";

  globval.rad_D = false;

  Ring_GetTwiss(true, 0e0);
  printglob();

  A = putlinmat(6, globval.Ascr);

  M = putlinmat(6, globval.OneTurnMat);

  printf("\nA.R.A^-1:\n");
  prt_lin_map(3, Inv(A)*M*A);

  GetEmittance(ElemIndex("cav"), true);

  D.zero();
  for (int k = 0; k < 2*nd_tps; k++)
    D[k] = globval.D_rad[k/2]*tps(0e0, k+1);

  printf("\nDiffusion Matrix - Floquet Space:\n");
  prt_lin_map(3, D);

  D =  A*D*tp_map(3, A);

  printf("\nDiffusion Matrix - Phase Space:\n");
  prt_lin_map(3, D);

  printf("\nDiffusion Matrix:\n");
#if 0
  for (int k = 0; k < 6; k++)
    globval.Diff_mat[k] /= 2e0;
#endif
  prt_lin_map(3, globval.Diff_mat);
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
  trace            = false;
  reverse_elem     = true;
  globval.mat_meth = false;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  set_state();

  Ring_GetTwiss(true, 0e0);
  printglob();

  tst_D();
}
