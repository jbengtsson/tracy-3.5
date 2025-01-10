#define NO 1

#include "tracy_lib.h"

int no_tps = NO;

const int    n_turn    = 2064;
const double delta_max = 2e-2;


void tst_lsoc(void)
{

  const long   seed   = 1121;
  const int    n_corr = 5;
  const double dx[]   = {100e-6, 100e-6};

  const int
    bpm   = ElemIndex("bpm"),
    hcorr = ElemIndex("chv"),
    vcorr = ElemIndex("chv");

  trace = true;

  iniranf(seed); setrancut(1e0);

  gcmat(bpm, hcorr, 1); gcmat(bpm, vcorr, 2);

  misalign_rms_type(Quad, dx[X_], dx[Y_], 0e0, true);
    
  orb_corr(n_corr);

  prt_cod("orb_corr.out", globval.bpm, true);
}


int main(int argc, char *argv[])
{
  double eps_x, sigma_delta, U_0, J[3], tau[3], I[6];

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;
  globval.Cart_Bend  = false; globval.dip_edge_fudge = true;

  reverse_elem = !false;

  trace = false;

  globval.mat_meth = false;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  Ring_GetTwiss(true, 0e0); printglob();

  if (!false) tst_lsoc();

  if (globval.mat_meth) get_eps_x(eps_x, sigma_delta, U_0, J, tau, I, true);

  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);
  prtmfile("flat_file.dat");
  exit(0);

  if (true) GetEmittance(ElemIndex("cav"), false, true);

  if (true) {
    globval.Cavity_on = true;
    get_dynap(delta_max, 25, n_turn, false);
  }
}
