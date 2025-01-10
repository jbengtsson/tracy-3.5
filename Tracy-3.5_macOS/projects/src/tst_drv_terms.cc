#define NO 1

#include "tracy_lib.h"

#include "drv_terms.cc"

int no_tps = NO;


void tst_drv_terms()
{
  const double
    delta_eps       = 1e-10,
    twoJ[]          = {1e0, 1e0},
    delta           = 1e0,
    twoJ_delta[]    = {1e0, 1e0},
    scl_ksi_1       = 1e0,
    scl_h_3         = 1e0,
    scl_h_3_delta   = 1e0,
    scl_h_4         = 1e0,
    scl_ksi_2       = 1e0,
    scl_chi_2       = 1e0,
    scl_chi_delta_2 = 1e0,
    scl_ksi_3       = 1e0;

  drv_terms_type rdt;

  rdt.get_h_scl
    (scl_ksi_1, scl_h_3, scl_h_3_delta, scl_h_4, scl_ksi_2, scl_chi_2,
     scl_chi_delta_2, scl_ksi_3);
  rdt.get_h(delta_eps, twoJ, delta, twoJ_delta);
  rdt.print();
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
  globval.mat_meth = !false;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  set_state();

  Ring_GetTwiss(true, 0e0);
  printglob();

  tst_drv_terms();
}
