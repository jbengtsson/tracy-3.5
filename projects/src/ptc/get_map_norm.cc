#include <assert.h>

#define NO 4

#include "tracy_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;


const double
  beta_inj[] = {6.0, 3.0},
  A_max[]    = {6e-3, 3e-3},
  delta_max  = 6e-2,
  twoJ[]     = {sqr(A_max[X_])/beta_inj[X_], sqr(A_max[Y_])/beta_inj[Y_]};


void get_map_norm(void)
{
  tps          h_re, h_im, g, g_re, g_im;
  ss_vect<tps> R, Id_scl;


  daeps_(eps_tps);

  Id_scl.identity();
  for (auto k = 0; k < 4; k++)
    Id_scl[k] *= sqrt(twoJ[k/2]);
  Id_scl[delta_] *= delta_max;

  danot_(no_tps-1);

  get_map(false);

  printf("\nM:");
  prt_lin_map(3, map);

  danot_(no_tps);

  MNF = MapNorm(map, no_tps);

  auto M_Fl = Inv(MNF.A0*MNF.A1)*map*MNF.A0*MNF.A1;
  auto h = LieFact_DF(M_Fl, R);

  CtoR(h*Id_scl, h_re, h_im);
  cout << "\nh_re:" << h_re;
  cout << "\nh_im:" << h_im;
}


int main(int argc, char *argv[])
{

  // disable from TPSALib- and LieLib log messages
  idprset_(-1);

  daeps_(1e-30);

  if (false)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  danot_(1);

   // Ring_GetTwiss(true, 0.0);
   // printglob();

  get_map_norm();
}
