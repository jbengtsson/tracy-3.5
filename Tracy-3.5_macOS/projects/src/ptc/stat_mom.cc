#define NO 2

#include "tracy_lib.h"

int
  no_tps   = NO,
  ndpt_tps = 5;


void get_map(void)
{
  long int     lastpos;
  int          k;
  tps          sigma;
  ss_vect<tps> Id, M, A_A_t;

  danot_(1);

  Id.identity();
  M.identity();
  Cell_Pass(0, globval.Cell_nLoc, M, lastpos);
  MNF = MapNorm(M, 1);
  A_A_t = MNF.A1*tp_S(3, MNF.A1);

  prt_lin_map(3, M);
  prt_lin_map(3, A_A_t);

  danot_(NO);

  sigma = 0e0;
  for (k = 0; k < 2; k++)
    sigma += sqr(Id[2*k]) + sqr(Id[2*k+1]);
  sigma = sigma*Inv(MNF.A1);
  cout << "\n  sigma = " << sigma << "\n";

  sigma = sigma*Inv(M);
  cout << "\n  sigma = " << sigma << "\n";

  M.identity();
  Cell_Pass(0, 10, M, lastpos);
  sigma = sigma*Inv(M);
  cout << "\n  sigma = " << sigma << "\n";
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
  reverse_elem     = true;
  globval.mat_meth = false;

  // disable from TPSALib- and LieLib log messages
  idprset_(-1);

  if (!true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  set_state();

  get_map();
}
