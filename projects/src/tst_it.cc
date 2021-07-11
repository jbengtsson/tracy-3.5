#define _GLIBCXX_DEBUG 1

#define NO_TPSA 1

#include "tracy_lib.h"

int no_tps = NO_TPSA;

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


void cod_stat(const LatticeType &lat, double mean[], double sigma[],
	      double xmax[], long lastpos, bool all)
{
  long   i, n, loc;
  int    j;
  double sum[2], sum2[2];

  // printf("  %2d %3d", (int)bpms_[0].size(), bpms_[2][15]);

  for (j = 0; j < 2; j++) {
    sum[j] = 0e0; sum2[j] = 0e0; xmax[j] = 0e0;
  }

  n = 0;
  if (all) {
    for (i = 0; i < lastpos; i++) {
      n++;
      for (j = 0; j < 2; j++) {
	sum[j] += lat.elems[i]->BeamPos[j*2];
	sum2[j] += sqr(lat.elems[i]->BeamPos[j*2]);
	xmax[j] = max(xmax[j], fabs(lat.elems[i]->BeamPos[j*2]));
      }
    }
  } else {
    for (i = 0; i < n_bpm_[X_]; i++) {
      n++;
      for (j = 0; j < 2; j++) {
	loc = bpms_[j][i];
	sum[j] += lat.elems[loc]->BeamPos[j*2];
	sum2[j] += sqr(lat.elems[loc]->BeamPos[j*2]);
	xmax[j] = max(xmax[j], fabs(lat.elems[loc]->BeamPos[j*2]));
      }
    }
  }

  for (j = 0; j < 2; j++) {
    if (n != 0)
      mean[j] = sum[j] / n;
    else
      mean[j] = -1e0;
    if (n != 0 && n != 1) {
      sigma[j] = (n*sum2[j]-sqr(sum[j]))/(n*(n-1e0));
    } else
      sigma[j] = 0e0;
    if (sigma[j] >= 0e0)
      sigma[j] = sqrt(sigma[j]);
    else
      sigma[j] = -1e0;
  }
}


bool orb_corr(LatticeType &lat, const int n_orbit)
{
  bool   cod = false;
  int    i;
  long   lastpos;
  double xmean[2], xsigma[2], xmax[2];

  lat.conf.CODvect.zeros();
  for (i = 1; i <= n_orbit; i++) {
    cod = lat.getcod(0e0, lastpos);
    if (cod) {
      cod_stat(lat, xmean, xsigma, xmax, lat.conf.Cell_nLoc, true);
      printf("\nRMS orbit [mm]: (%8.1e +/- %7.1e, %8.1e +/- %7.1e)\n",
	     1e3*xmean[X_], 1e3*xsigma[X_], 1e3*xmean[Y_], 1e3*xsigma[Y_]);
      lsoc(lat, 1, 1e0); lsoc(lat, 2, 1e0);

      cod = lat.getcod(0e0, lastpos);
      if (cod) {
	cod_stat(lat, xmean, xsigma, xmax, lat.conf.Cell_nLoc, true);
	printf("RMS orbit [mm]: (%8.1e +/- %7.1e, %8.1e +/- %7.1e)\n",
	       1e3*xmean[X_], 1e3*xsigma[X_], 1e3*xmean[Y_], 1e3*xsigma[Y_]);
      } else
	printf("orb_corr: failed\n");
    } else {
      printf("orb_corr: failed\n");
      break;
    }
  }

  lat.prt_cod("orb_corr.out", true);

  return cod;
}


void tst_lsoc(LatticeType &lat)
{
  const long
    seed   = 1121;
  const int
    n_corr = 5,
    bpm    = lat.ElemIndex("bpm"),
    hcorr  = lat.ElemIndex("chv"),
    vcorr  = lat.ElemIndex("chv");
  const
    double dx[]   = {100e-6, 100e-6};

  iniranf(seed); setrancut(1e0);
  gcmat(lat, bpm, hcorr, 1); gcmat(lat, bpm, vcorr, 2);
  misalign_rms_type(lat, Quad, dx[X_], dx[Y_], 0e0, true);
  orb_corr(lat, n_corr);
  lat.prtmfile("flat_file_err.dat");
}


void get_lat(const string &file_name, LatticeType &lat)
{
  double eps_x, sigma_delta, U_0, J[3], tau[3], I[6];

  const string str = file_name;

  lat.conf.trace        = false;
  lat.conf.reverse_elem = !false;
  lat.conf.mat_meth     = !false;

  lat.conf.H_exact    = false; lat.conf.quad_fringe = false;
  lat.conf.Cavity_on  = false; lat.conf.radiation   = false;
  lat.conf.emittance  = false; lat.conf.IBS         = false;
  lat.conf.pathlength = false; lat.conf.bpm         = 0;
  lat.conf.Cart_Bend  = false; lat.conf.dip_edge_fudge = true;

  if (true) {
    lat.Lat_Read(file_name);
    lat.Lat_Init();
  } else
    lat.rdmfile(file_name);

  lat.ChamberOff();

  lat.conf.CODimax = 10;

  if (false) {
    lat.prt_elems();
    lat.prt_fams();
  }

  if (false) {
    long int lastpos;
    double   xmean[2], xsigma[2], xmax[2];

    lat.getcod(0e0, lastpos);
    cod_stat(lat, xmean, xsigma, xmax, lat.conf.Cell_nLoc, true);
    printf("\nRMS orbit [mm]: (%8.1e +/- %7.1e, %8.1e +/- %7.1e)\n",
	   1e3*xmean[X_], 1e3*xsigma[X_], 1e3*xmean[Y_], 1e3*xsigma[Y_]);

    lat.getcod(0e0, lastpos);

    cod_stat(lat, xmean, xsigma, xmax, lat.conf.Cell_nLoc, true);
    printf("RMS orbit [mm]: (%8.1e +/- %7.1e, %8.1e +/- %7.1e)\n",
	   1e3*xmean[X_], 1e3*xsigma[X_], 1e3*xmean[Y_], 1e3*xsigma[Y_]);

    exit(0);
  }

  lat.Ring_GetTwiss(true, 0e-3); printglob(lat);
  if (lat.conf.mat_meth)
    lat.get_eps_x(eps_x, sigma_delta, U_0, J, tau, I, true);

  if (false) tst_lat(lat);

  lat.prt_lat("linlat1.out", true);
  lat.prt_lat("linlat.out", true, 10);
  lat.prtmfile("flat_file.dat");

  if (!lat.conf.mat_meth) GetEmittance(lat, lat.ElemIndex("cav"), true);

  if (false) {
    lat.conf.Cavity_on = true;
    get_dynap(lat, delta_max, 25, n_turn, false);
  }
}


int main(int argc, char *argv[])
{
  string      file_name;
  LatticeType lat;

  get_lat(argv[1], lat);

  if (!false) tst_lsoc(lat);
}
