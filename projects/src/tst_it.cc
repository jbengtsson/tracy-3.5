#define NO_TPSA 1

#include "tracy_lib.h"

int no_tps = NO_TPSA;

const int    n_turn    = 2064;
const double delta_max = 2e-2;


long int rseed0, rseed;
double   normcut_;

const long int k = 19, c = 656329L, m = 100000001;
   

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


void get_lat(const char *file_name, LatticeType &lat)
{
  int    k;
  double eps_x, sigma_delta, U_0, J[3], tau[3], I[6];
  FILE   *inf, *outf;

  const string str = file_name;

  lat.conf.trace = false;
  lat.conf.reverse_elem = !false;
  lat.conf.mat_meth = !false;

  lat.conf.H_exact    = false; lat.conf.quad_fringe = false;
  lat.conf.Cavity_on  = false; lat.conf.radiation   = false;
  lat.conf.emittance  = false; lat.conf.IBS         = false;
  lat.conf.pathlength = false; lat.conf.bpm         = 0;
  lat.conf.Cart_Bend  = false; lat.conf.dip_edge_fudge = true;

  inf  = file_read((str + ".lat").c_str());
  outf = file_write((str + ".lax").c_str());

  lat.Lattice_Read(inf, outf);
  lat.Lat_Init();
  lat.ChamberOff();

  lat.conf.CODimax = 10;

  if (false) {
    lat.prt_fam();
    lat.prt_elem();
  }

  if (!false) {
    // for (k = 0; k <= lat.conf.Cell_nLoc; k++)
    //   if (lat.elems[k]->Pkind == Mpole) {
    //   MpoleType *M = dynamic_cast<MpoleType*>(lat.elems[k]);
    //   char str[20];
    //   sprintf(str, "\n%2d: %8s", k, lat.elems[k]->PName);
    //   M->M_lin.print(str);
    // }

    long int     lastpos;
    ss_vect<tps> ps;
    ps.identity();
    lat.Cell_Pass(0, lat.conf.Cell_nLoc, ps, lastpos);
    prt_lin_map(3, ps);
  }

  lat.Ring_GetTwiss(true, 0e0); printglob(lat);
  exit(0);
  if (lat.conf.mat_meth)
    lat.get_eps_x(eps_x, sigma_delta, U_0, J, tau, I, true);

  if (false) tst_lat(lat);

  lat.prt_lat("linlat1.out", true);
  lat.prt_lat("linlat.out", true, 10);
  lat.prtmfile("flat_file.dat");

  if (true) GetEmittance(lat, ElemIndex("cav"), true);

  if (!true) {
    lat.conf.Cavity_on = true;
    get_dynap(lat, delta_max, 25, n_turn, false);
  }
}


void cod_stat(LatticeType &lat, double *mean, double *sigma, double *xmax,
	     long lastpos, bool all)
{
  long    i, n, loc;
  int     j;
  Vector2 sum, sum2;

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
    for (i = 1; i <= n_bpm_[X_]; i++) {
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
  bool    cod = false;
  int     i;
  long    lastpos;
  Vector2 xmean, xsigma, xmax;

  printf("\n");
  lat.conf.CODvect.zeros();
  for (i = 1; i <= n_orbit; i++) {
    cod = lat.getcod(0e0, lastpos);
    if (cod) {
      cod_stat(lat, xmean, xsigma, xmax, lat.conf.Cell_nLoc, false);
      printf("\n");
      printf("RMS orbit [mm]: (%8.1e +/- %7.1e, %8.1e +/- %7.1e)\n",
	     1e3*xmean[X_], 1e3*xsigma[X_], 1e3*xmean[Y_], 1e3*xsigma[Y_]);
      lsoc(lat, 1, 1e0); lsoc(lat, 2, 1e0);
      cod = lat.getcod(0e0, lastpos);
      if (cod) {
	cod_stat(lat, xmean, xsigma, xmax, lat.conf.Cell_nLoc, false);
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


void iniranf(const long i) { rseed0 = i; rseed = i; }

void newseed(void) { rseed0 = (k*rseed0+c) % m; rseed = (rseed0+54321) % m; }

double ranf_(void)
{
  /* Generate random number with rectangular distribution */
  rseed = (k*rseed+c) % m; return (rseed/1e8);
}


void setrancut(const double cut)
{

  printf("\n");
  printf("setrancut: cut set to %3.1f\n", cut);

  normcut_ = cut;
}


void tst_lsoc(LatticeType &lat)
{

  const long   seed   = 1121;
  const int    n_corr = 5;
  const double dx[]   = {100e-6, 100e-6};

  const int
    bpm   = ElemIndex("bpm"),
    hcorr = ElemIndex("chv"),
    vcorr = ElemIndex("chv");

  lat.conf.trace = true;

  iniranf(seed); setrancut(1e0);

  gcmat(lat, bpm, hcorr, 1); gcmat(lat, bpm, vcorr, 2);

  misalign_rms_type(lat, Quad, dx[X_], dx[Y_], 0e0, true);
    
  orb_corr(lat, n_corr);
}


int main(int argc, char *argv[])
{
  string      file_name;
  double      eps_x, sigma_delta, U_0, J[3], tau[3], I[6];
  LatticeType lat;

  get_lat(argv[1], lat);

  if (!false) tst_lsoc(lat);
}
