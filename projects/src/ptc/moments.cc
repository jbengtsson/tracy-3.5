#define NO 2

#include "tracy_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;


void track_twoJ()
{
  long int     lastpos;
  int          k;
  tps          twoJ;
  ss_vect<tps> Id, M;

  const double
    alpha[] = {Cell[0].Alpha[X_], Cell[0].Alpha[Y_]},
    beta[]  = {Cell[0].Beta[X_],  Cell[0].Beta[Y_]},
    gamma[] = {(1e0+sqr(alpha[X_]))/beta[X_], (1e0+sqr(alpha[Y_]))/beta[Y_]};

  Id.identity();
  printf("\n  alpha_x = %5.3f beta_x = %5.3f gamma_x = %5.3f\n",
	 alpha[X_], beta[X_], gamma[X_]);

  M.identity();
  Cell_Pass(0, globval.Cell_nLoc, M, lastpos);
  printf("\nM:\n");
  prt_lin_map(3, M);

  twoJ = 0e0;
  for (k = 0; k < 2; k++)
    twoJ +=
      gamma[k]*sqr(Id[2*k]) + 2e0*alpha[k]*Id[2*k]*Id[2*k+1]
      + beta[k]*sqr(Id[2*k+1]);

  cout << scientific << setprecision(3) <<  "\ntwoJ:\n" << twoJ << "\n";
  twoJ = twoJ*M;
  cout << scientific << setprecision(3) <<  "\ntwoJ:\n" << twoJ << "\n";
}


int main(int argc, char *argv[])
{
  globval.H_exact    = false; globval.quad_fringe    = false;
  globval.Cavity_on  = false; globval.radiation      = false;
  globval.emittance  = false; globval.IBS            = false;
  globval.pathlength = false; globval.bpm            = 0;
  globval.Cart_Bend  = false; globval.dip_edge_fudge = true;
  globval.mat_meth   = false;

  // disable from TPSALib- and LieLib log messages
  idprset_(-1);

  if (false)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  Ring_GetTwiss(true, 0e0); printglob();

  track_twoJ();
}

