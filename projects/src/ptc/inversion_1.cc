#define NO 4

#include "tracy_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;


void compute_map(void)
{
  long int     lastpos;
  tps          h;
  ss_vect<tps> map, R;

  map.identity();
  Cell_Pass(0, globval.Cell_nLoc, map, lastpos);
  prt_lin_map(3, map);

  h = LieFact_DF(map, R);
  daeps_(1e0);
  cout << scientific << setprecision(5) << setw(13) << 1e0*h << "\n";
}


void chk_sympl()
{
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
  globval.EPU            = !true;
}


int main(int argc, char *argv[])
{

  globval.mat_meth = false;

  FieldMap_filetype = 6;

  if (!true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  set_state();

  // Disable from TPSALib and LieLib log messages.
  idprset(-1);

  if (!false) {
    compute_map();
    assert(false);
  }
  

  Ring_GetTwiss(true, 0e0);
  printglob();

  prtmfile("flat_file.dat");
  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);
}
