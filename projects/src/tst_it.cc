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
    lat.elems[k]->Elem_Pass(lat.conf, map);
  prt_lin_map(3, map);
}


void get_lat(const char *file_name)
{
  LatticeType lat;
  double      eps_x, sigma_delta, U_0, J[3], tau[3], I[6];
  FILE        *inf, *outf;

  const string str = file_name;

  lat.conf.H_exact    = false; lat.conf.quad_fringe = false;
  lat.conf.Cavity_on  = false; lat.conf.radiation   = false;
  lat.conf.emittance  = false; lat.conf.IBS         = false;
  lat.conf.pathlength = false; lat.conf.bpm         = 0;
  lat.conf.Cart_Bend  = false; lat.conf.dip_edge_fudge = true;

  reverse_elem = !false;

  trace = false;

  lat.conf.mat_meth = !false;

  inf  = file_read((str + ".lat").c_str());
  outf = file_write((str + ".lax").c_str());

  lat.Lattice_Read(inf, outf);
  lat.Lat_Init();
  lat.ChamberOff();

  if (false) {
    lat.prt_fam();
    lat.prt_elem();
  }

  lat.Ring_GetTwiss(true, 0e0); printglob(lat);
  if (lat.conf.mat_meth)
    lat.get_eps_x(eps_x, sigma_delta, U_0, J, tau, I, true);

  if (false) tst_lat(lat);

  lat.prt_lat("linlat1.out", true);
  lat.prt_lat("linlat.out", true, 10);
  lat.prtmfile("flat_file.dat");

  exit(0);

  if (true) GetEmittance(lat, ElemIndex("cav"), true);

  if (true) {
    lat.conf.Cavity_on = true;
    get_dynap(lat, delta_max, 25, n_turn, false);
  }
}


int main(int argc, char *argv[])
{
  string file_name;
  double eps_x, sigma_delta, U_0, J[3], tau[3], I[6];

  get_lat(argv[1]);
}
