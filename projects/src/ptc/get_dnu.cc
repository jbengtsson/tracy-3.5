
#define NO 6

#include "tracy_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;


// MAX-IV       1,
// SLS-2        2,
// M-H6BAi      3,
// M-H6BA-0-.-. 4,
// DIAMOND      5,
// DELTA        6,
// ALS-U        7.

const bool set_dnu  = false;
const double
  beta_inj[] = {8.7, 2.1},
  A_max[]    = {5e-3, 0.3e-3},
  twoJ[]     = {sqr(A_max[X_])/beta_inj[X_], sqr(A_max[Y_])/beta_inj[Y_]},
  delta_max  = 3e-2,
  dnu[]      = {0.1/6.0, 0.0};

const char home_dir[] = "/home/bengtsson";

double       nu0[2];
ss_vect<tps> A_inv, nus;


void get_map_n(const int n)
{
  ss_vect<tps> map2, map4;

  get_map(false);

  switch (n) {
  case 1:
    break;
  case 2:
    map = map*map;
    break;
  case 3:
    map2 = map*map; map = map2*map;
    break;
  case 4:
    map2 = map*map; map = map2*map2;
    break;
  case 5:
    map2 = map*map; map = map2*map2*map;
    break;
  case 6:
    map2 = map*map; map = map2*map2*map2;
    break;
  case 7:
    map2 = map*map; map4 = map2*map2; map = map4*map2*map;
    break;
  case 8:
    map2 = map*map; map4 = map2*map2; map = map4*map4;
    break;
  case 9:
    map2 = map*map; map4 = map2*map2; map = map4*map4*map;
    break;
  case 10:
    map2 = map*map; map4 = map2*map2; map = map4*map4*map2;
    break;
  case 11:
    map2 = map*map; map4 = map2*map2; map = map4*map4*map2*map;
    break;
  case 12:
    map2 = map*map; map4 = map2*map2; map = map4*map4*map4;
    break;
  case 13:
    map2 = map*map; map4 = map2*map2; map = map4*map4*map4*map;
    break;
  case 14:
    map2 = map*map; map4 = map2*map2; map = map4*map4*map4*map2;
    break;
  case 15:
    map2 = map*map; map4 = map2*map2; map = map4*map4*map4*map2*map;
    break;
  default:
    cout << "get_map_n: n not defined " << n << endl;
    exit(1);
    break;
  }
}


tps g_renorm(const double nu0_x, const double nu0_y,
	     const double nu1_x, const double nu1_y,
	     const tps &g)
{
  // Renormalize g: (1-R^-1)^-1 * h 

  int           i, j, k, l, m;
  long int      jj1[ss_dim], jj2[ss_dim];
  double        re, im, cotan0, cotan1, cotan0_sqr;
  tps           h_re, h_im, g_re, g_im, G_re, G_im, mn1, mn2;
  ss_vect<tps>  Id;

  CtoR(g, g_re, g_im);

  for (k = 0; k < ss_dim; k++) {
    jj1[k] = 0; jj2[k] = 0;
  }

  Id.identity(); G_re = 0.0; G_im = 0.0;
  for (i = 0; i <= no_tps; i++) {
    jj1[x_] = i; jj2[px_] = i;
    for (j = 0; j <= i; j++) {
      jj1[px_] = j; jj2[x_] = j;
      for (k = 0; k <= no_tps; k++) {
	jj1[y_] = k; jj2[py_] = k;
	for (l = 0; l <= no_tps; l++) {
	  jj1[py_] = l; jj2[y_] = l;

	  if (i+j+k+l <= no_tps) {
	    cotan0 = 1.0/tan(((i-j)*nu0_x+(k-l)*nu0_y)*M_PI);
	    cotan0_sqr = sqr(cotan0);
	    cotan1 = 1.0/tan(((i-j)*nu1_x+(k-l)*nu1_y)*M_PI);
	    mn1 =
	      pow(Id[x_], i)*pow(Id[px_], j)*pow(Id[y_], k)*pow(Id[py_], l);
	    mn2 =
	      pow(Id[x_], j)*pow(Id[px_], i)*pow(Id[y_], l)*pow(Id[py_], k);

	    for (m = 0; m <= no_tps; m++) {
	      if (i+j+k+l+m <= no_tps) {
		jj1[delta_] = m; jj2[delta_] = m;
		if ((i != j) || (k != l)) {
		  re = g_re[jj1]; im = g_im[jj1];

		  // compute h
		  h_re = (re+cotan0*im)*2.0/(1.0+cotan0_sqr);
		  h_im = (im-cotan0*re)*2.0/(1.0+cotan0_sqr);

		  // renormalize g
		  G_re += (h_re-cotan1*h_im)*(mn1+mn2)*pow(Id[delta_], m)/2.0;
		  G_im += (h_im+cotan1*h_re)*(mn1-mn2)*pow(Id[delta_], m)/2.0;
		  g_re.pook(jj2, 0.0); g_im.pook(jj2, 0.0);
		}
	      }
	    }
	  }
	}
      }
    }
  }

  return RtoC(G_re, G_im);
}


void get_map_normal_form()
{

  danot_(no_tps);

  MNF = MapNorm(map, no_tps);
}


tps get_H(void)
{
  int          i;
  tps          H, gn;
  ss_vect<tps> Id, Mn;

  const bool prt = false;

  // Construct generator.
  // K is in Dragt-Finn form but the generators commute.
  H = MNF.K; Id.identity();
  for (i = no_tps; i >= 3; i--) {
    gn = Take(MNF.g, i); H = H*LieExp(-gn, Id);
  }

  if (prt) cout << (H-H*Inv(MNF.A0*MNF.A1)*map*MNF.A0*MNF.A1);

  if (false) {
    // Normalize map (=> Map_res)
    for (i = 3; i <= no_tps; i++) {
      gn = Take(MNF.g, i); Mn = LieExp(gn, Id); map = Inv(Mn)*map*Mn;
    }
  }

  return H;
}


void get_A(void)
{
  int          j;
  long int     jj[ss_dim];
  tps          gn;
  ss_vect<tps> Id, A;

  Id.identity(); A = MNF.A1;
  for (j = no_tps; j >= 3; j--) {
    gn = Take(MNF.g, j); A = A*LieExp(gn, Id);
  }

  for (j = 0; j < nv_tps; j++)
    jj[j] = (j < 4)? 1 : 0;

  A_inv = PInv(A, jj);
}


void get_twoJ(const ss_vect<double> &ps, double twoJ[])
{
  int             j;
  ss_vect<double> z;
  ss_vect<tps>    Id;

  z = (A_inv*ps).cst();

  for (j = 0; j < 2; j++)
    twoJ[j] = sqr(z[2*j]) + sqr(z[2*j+1]);
}


void get_dnu(const double Ax_max, const double Ay_max, const double delta_max)
{
  char            str[max_str];
  int             i;
  double          twoJ[2], nux, nuy;
  ss_vect<double> ps;
  ss_vect<tps>    Id;
  ifstream        inf;
  ofstream        outf;

  const int    n_ampl = 25, n_delta = 20;
  const double A_min = 1e-6;

  if (false) {
    sprintf(str, "%s%s", home_dir, "/Thor-2.0/thor/wrk");
    file_rd(inf, strcat(str, "/nus.dat"));
    inf >> nus[3] >> nus[4];
    inf.close();
  }

//  sprintf(str, "%s%s", home_dir, "/projects/src/");
//  file_wr(outf, strcat(str, "dnu_dAx_pert.out"));
  file_wr(outf, "dnu_dAx_pert.out");
  Id.zero(); ps.zero();
  for (i = -n_ampl; i <= n_ampl; i++) {
    ps[x_] = i*Ax_max/n_ampl;
    if (ps[x_] == 0.0) ps[x_] = A_min;
    ps[y_] = A_min;
    get_twoJ(ps, twoJ);
    Id[x_] = sqrt(twoJ[X_]); Id[px_] = sqrt(twoJ[X_]);
    Id[y_] = sqrt(twoJ[Y_]); Id[py_] = sqrt(twoJ[Y_]);
    nux = (nus[3]*Id).cst(); nuy = (nus[4]*Id).cst();

    outf << scientific << setprecision(3)
	 << setw(12) << 1e3*ps[x_] << setw(12) << 1e3*ps[y_]
	 << fixed << setprecision(5)
	 << setw(9) << nux << setw(9) << nuy << endl;
  }
  outf.close();

//  sprintf(str, "%s%s", home_dir, "/projects/src/");
//  file_wr(outf, strcat(str, "dnu_dAy_pert.out"));
  file_wr(outf, "dnu_dAy_pert.out");
  Id.zero(); ps.zero();
  for (i = -n_ampl; i <= n_ampl; i++) {
    ps[x_] = A_min;
    ps[y_] = i*Ay_max/n_ampl;
    if (ps[y_] == 0.0) ps[y_] = A_min;
//    get_twoJ(2, ps, MNF.A1, twoJ);
    get_twoJ(ps, twoJ);
    Id[x_] = sqrt(twoJ[X_]); Id[px_] = sqrt(twoJ[X_]);
    Id[y_] = sqrt(twoJ[Y_]); Id[py_] = sqrt(twoJ[Y_]);
    nux = (nus[3]*Id).cst(); nuy = (nus[4]*Id).cst();

    outf << scientific << setprecision(3)
	 << setw(12) << 1e3*ps[x_] << setw(12) << 1e3*ps[y_]
	 << fixed << setprecision(6)
	 << setw(10) << nux << setw(10) << nuy << endl;
  }
  outf.close();

//  sprintf(str, "%s%s", home_dir, "/projects/src/");
//  file_wr(outf, strcat(str, "chrom2_pert.out"));
  file_wr(outf, "chrom2_pert.out");
  Id.zero(); ps.zero();
  for (i = -n_delta; i <= n_delta; i++) {
    ps[delta_] = i*delta_max/n_delta; Id[delta_] = ps[delta_];

    nux = (nus[3]*Id).cst(); nuy = (nus[4]*Id).cst();

    outf << scientific << setprecision(3)
	 << setw(12) << 1e2*ps[delta_]
	 << fixed << setprecision(5)
	 << setw(9) << nux << setw(9) << nuy << endl;
  }
  outf.close();
}


void get_dnu2(const double Ax_max, const double Ay_max, const double delta)
{
  char            str[max_str];
  int             i, j;
  double          twoJ[2], nux, nuy;
  ss_vect<double> ps;
  ss_vect<tps>    Id;
  ifstream        inf;
  ofstream        outf;

  const int n_ampl = 10;

  if (false) {
    sprintf(str, "%s%s", home_dir, "/Thor-2.0/thor/wrk");
    file_rd(inf, strcat(str, "/nus.dat"));
    inf >> nus[3] >> nus[4];
    inf.close();
  }

  file_wr(outf, "dnu_dAxy_pert.out");
  Id.zero(); Id[delta_] = delta;
  ps.zero();
  for (i = -n_ampl; i <= n_ampl; i++) {
    for (j = -n_ampl; j <= n_ampl; j++) {
      ps[x_] = i*Ax_max/n_ampl; ps[y_] = j*Ay_max/n_ampl;
      get_twoJ(ps, twoJ);
      Id[x_] = sqrt(twoJ[X_]); Id[px_] = sqrt(twoJ[X_]);
      Id[y_] = sqrt(twoJ[Y_]); Id[py_] = sqrt(twoJ[Y_]);
      nux = (nus[3]*Id).cst(); nuy = (nus[4]*Id).cst();

      outf << scientific << setprecision(3)
	   << setw(12) << 1e3*ps[x_] << setw(12) << 1e3*ps[y_]
	   << fixed << setprecision(5)
	   << " " << setw(8) << nux << " " << setw(8) << nuy << endl;
    }
    outf << endl;
  }
  outf.close();
}


int main(int argc, char *argv[])
{
  tps          H, H_re, H_im, g_re, g_im, K_re, K_im;
  ss_vect<tps> Id;
  ofstream     outf;

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  // disable from TPSALib- and LieLib log messages
  idprset_(-1);

  if (false)
    Read_Lattice(argv[1]);
  else {
    rdmfile(argv[1]);
  }

  globval.EPU = true;

  danot_(1);

//  Ring_GetTwiss(true, 0.0); printglob();

  prt_lat("linlat.out", globval.bpm, true);

  danot_(no_tps-1);

//  get_map_n(3);
  get_map(false);

  danot_(no_tps);
  get_map_normal_form(); nus = dHdJ(MNF.K);

  CtoR(MNF.K, K_re, K_im); CtoR(MNF.g, g_re, g_im);

  Id.identity();
  Id[x_] *= sqrt(twoJ[X_]); Id[px_] *= sqrt(twoJ[X_]);
  Id[y_] *= sqrt(twoJ[Y_]); Id[py_] *= sqrt(twoJ[Y_]);
  Id[delta_] *= delta_max;

  if (true) {
//    H = get_H(); CtoR(H, H_re, H_im);

//    outf.open("H.dat", ios::out);
//    outf << H_re*Id;
//    outf.close();

    outf.open("K.dat", ios::out);
    outf << K_re;
    outf.close();

    outf.open("g.dat", ios::out);
    outf << g_im;
    outf.close();

//    cout << N*nus[3]*Id << N*nus[4]*Id;

    outf.open("nus.dat", ios::out);
    outf << nus[3] << nus[4];
    outf.close();
  }

  daeps_(eps_tps);

  // Note, nus are in Floquet space.
  get_A();

  get_dnu(A_max[X_], A_max[Y_], delta_max);

//  get_dnu2(Ax, Ay, 0.0);

//  get_dnu2(Ax, Ay, 2.5e-2);
}
