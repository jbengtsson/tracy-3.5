#define NO 1

#include "tracy_lib.h"

int no_tps = NO;

const int    n_turn    = 2064;
const double delta_max = 2e-2;


void tst_lat(LatticeType &lat)
{
  int          k;
  ss_vect<tps> map;

  printf("\ntst_lat\n");
  map.identity();
  for (k = 0; k <= 10; k++)
    lat.elems[k]->Elem_Pass(map);
  prt_lin_map(3, map);
}


void get_lat(const char *file_name)
{
  LatticeType lat;
  double      eps_x, sigma_delta, U_0, J[3], tau[3], I[6];
  FILE        *inf, *outf;

  const string str = file_name;

  inf  = file_read((str + ".lat").c_str());
  outf = file_write((str + ".lax").c_str());

  lat.Lattice_Read(inf, outf);
  lat.Lat_Init();
  lat.ChamberOff();

  if (false) {
    lat.prt_fam();
    lat.prt_elem();
  }

  lat.Ring_GetTwiss(true, 0e0); printglob(lat.elems[0]);
  if (globval.mat_meth)
    lat.get_eps_x(eps_x, sigma_delta, U_0, J, tau, I, true);

  if (false) tst_lat(lat);

  lat.prt_lat("linlat1.out", globval.bpm, true);
  lat.prt_lat("linlat.out", globval.bpm, true, 10);
  lat.prtmfile("flat_file.dat");
}


int main(int argc, char *argv[])
{
  string file_name;
  double eps_x, sigma_delta, U_0, J[3], tau[3], I[6];

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;
  globval.Cart_Bend  = false; globval.dip_edge_fudge = true;

  reverse_elem = !false;

  trace = false;

  globval.mat_meth = !false;

  if (!true) {
    if (true)
      Read_Lattice(argv[1]);
    else
      rdmfile(argv[1]);
    lat.Ring_GetTwiss(true, 0e0); printglob(lat.elems[0]);
    // if (globval.mat_meth)
    //   lat.get_eps_x(eps_x, sigma_delta, U_0, J, tau, I, true);
  } else
    get_lat(argv[1]);

  exit(0);

  if (true) GetEmittance(ElemIndex("cav"), true);

  if (true) {
    globval.Cavity_on = true;
    get_dynap(delta_max, 25, n_turn, false);
  }
}
