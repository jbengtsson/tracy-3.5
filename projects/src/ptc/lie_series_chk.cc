#define NO 4

#include "tracy_lib.h"


int no_tps   = NO,
    ndpt_tps = 5;


void compute_map_norm(void)
{
  tps          h, h_re, h_im, g_re, g_im, k_re, k_im;
  ss_vect<tps> M_lin;

  danot_(no_tps-1);
  get_map(false);
  danot_(no_tps);

  prt_lin_map(3, map);

  MNF = MapNorm(map, 1);
  CtoR(MNF.g, g_re, g_im);
  CtoR(MNF.K, k_re, k_im);

  daeps_(1e-8);
  cout << scientific << setprecision(5) << 1e0*g_re << 1e0*g_im;
  cout << scientific << setprecision(5) << 1e0*k_re << 1e0*k_im;
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

  if (!false)
    compute_map_norm();
}
