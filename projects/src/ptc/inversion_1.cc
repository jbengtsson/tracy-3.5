#define NO 4

#include "tracy_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;


ss_vect<tps> compute_map(void)
{
  long int     lastpos;
  tps          h;
  ss_vect<tps> map, R;

  map.identity();
  Cell_Pass(0, globval.Cell_nLoc, map, lastpos);
  return map;
}


void chk_sympl(ss_vect<tps> &map)
{
  const int dof = 3;

  Matrix       Omega_1;
  ss_vect<tps> Omega;

  Omega = get_S(dof);
  getlinmat(2*dof, map*Omega*tp_S(dof, map), Omega_1);
  printf("\nM^T*Omega*M:\n");
  prtmat(2*dof, Omega_1);
  for (int k = 0; k < dof; k++) {
    Omega_1[2*k][2*k+1] -= 1e0;
    Omega_1[2*k+1][2*k] += 1e0;
  }
  printf("\nM^T*Omega*M - Omega:\n");
  prtmat(2*dof, Omega_1);
}


void analyse_nl_dyn(ss_vect<tps> &map)
{
  tps g_re, g_im, k_re, k_im;

  MNF = MapNorm(map, 1);
  CtoR(MNF.g, g_re, g_im);
  CtoR(MNF.K, k_re, k_im);
  daeps_(1e0);
  cout << scientific << setprecision(5) << setw(13) << 1e0*g_im << "\n";
  cout << scientific << setprecision(5) << setw(13) << 1e0*k_re << "\n";
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
  globval.EPU            = true;
}


int main(int argc, char *argv[])
{
  ss_vect<tps> map;
  
  trace = false;

  globval.mat_meth = false;

  if (!true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  set_state();

  if (!false)
    no_sxt();

  // Disable from TPSALib and LieLib log messages.
  idprset(-1);

  daeps_(1e-30);

  if (!true) {
    Ring_GetTwiss(true, 0e0);
    printglob();
  }

  if (false) {
    globval.Cavity_on = true;

    map = compute_map();
    prt_lin_map(3, map);
    chk_sympl(map);
  }
  
  if (!false) {
    map = compute_map();
    analyse_nl_dyn(map);
  }
}
