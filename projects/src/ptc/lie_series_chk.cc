#define NO 3

#include "tracy_lib.h"


int no_tps   = NO,
    ndpt_tps = 5;


tps get_h_local(const ss_vect<tps> &map, const bool dragt_finn)
{
  ss_vect<tps> map1, R;

  if (dragt_finn)
    // Dragt-Finn factorization.
    return LieFact_DF(map, R);
  else {
    // Single Lie exponent.
    danot_(1);
    map1 = map;
    danot_(no_tps);
    return LieFact(map*Inv(map1));
  }
}


void get_K(void)
{
  tps g_re, g_im;

  danot_(no_tps-1);
  get_map(false);
  danot_(no_tps);
  MNF = MapNorm(map, 1);

  CtoR(MNF.g, g_re, g_im);
  daeps_(1e-8);
  cout << scientific << setprecision(5) << 1e0*g_im;
}


void compute_h(void)
{
  tps h, h_re, h_im;

  danot_(no_tps-1);
  get_map(false);

  prt_lin_map(3, map);

  danot_(no_tps);
  h = get_h_local(map, true);
  CtoR(h, h_re, h_im);

  cout << scientific << setprecision(5) << h_re << h_im;
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
  set_state();

  // Disable TPSALib & LieLib log messages.
  idprset(-1);

  daeps_(1e-30);

  if (false)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  compute_h();
}
