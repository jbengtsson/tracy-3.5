#include <assert.h>

#define NO 1

#include "tracy_lib.h"


int no_tps = NO;


void tst_D(void)
{
  // Average over one half synchrotron oscillation.
  const int n_turn = 215;

  long int        lastpos;
  ss_vect<double> ps, cod;
  ss_vect<tps>    A, M, D, Sigma;

  globval.rad_D = globval.Cavity_on = globval.radiation = true;

  getcod(0e0, lastpos);
  cod = globval.CODvect;

  ps.zero();
  ps += cod;
  globval.Diff_mat.zero();
  printf("\n  n_turn = %d\n", n_turn);
  for (int k = 0; k < n_turn; k++)
    Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
  for (int k = 0; k < 2*nd_tps; k++)
    globval.Diff_mat[k] /= n_turn;

  cout << scientific << setprecision(3)
       << "\ncod:\n" << setw(11) << cod << "\n";
  cout << scientific << setprecision(3)
       << "ps:\n" << setw(11) << ps << "\n";

  globval.rad_D = false;

  Ring_GetTwiss(true, 0e0);
  printglob();

  A = putlinmat(6, globval.Ascr);
  M = putlinmat(6, globval.OneTurnMat);

  printf("\n|M|-1 = %9.3e\n", det_map(3, M)-1e0);


  printf("\nA.R.A^-1:");
  prt_lin_map(3, Inv(A)*M*A);

  GetEmittance(ElemIndex("cav"), true);

  D.zero();
  for (int k = 0; k < 2*nd_tps; k++)
    D[k] = globval.D_rad[k/2]*tps(0e0, k+1);

  printf("\nDiffusion Matrix - Floquet Space:\n");
  prt_lin_map(3, D);

  D = A*D*tp_map(3, A);

  printf("\nDiffusion Matrix - Phase Space:");
  prt_lin_map(3, D);

  printf("\nDiffusion Matrix - Beam Envelope Theory:");
  printf("\nD - Phase Space:");
  prt_lin_map(3, globval.Diff_mat);

  printf("\nD - Floquet Space:");
  prt_lin_map(3, Inv(A)*globval.Diff_mat*Inv(tp_map(3, A)));

  Sigma.zero();
  for (int k = 0; k < 2*nd_tps; k++)
    Sigma[k] = globval.eps[k/2]*tps(0e0, k+1);

  printf("\nSigma - Floquet Space:");
  prt_lin_map(3, Sigma);

  Sigma = A*Sigma*tp_map(3, A);

  printf("\nSigma - Phase Space:");
  prt_lin_map(3, Sigma);

  printf("\nM.Sigma.M^tp + D - Sigma:");
  prt_lin_map(3, M*Sigma*tp_map(3, M)+D-Sigma);
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
  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);


  tst_D();
}
