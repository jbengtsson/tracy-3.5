#define NO 6

#include "tracy_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;


const double
  beta_inj[] = {2.8, 2.8},
  A_max[]    = {3e-3, 2.5e-3},
  delta_max  = 5e-2,
  twoJ[]     = {sqr(A_max[X_])/beta_inj[X_], sqr(A_max[Y_])/beta_inj[Y_]};

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

  z = (A_inv*ps).cst();

  for (j = 0; j < 2; j++)
    twoJ[j] = sqr(z[2*j]) + sqr(z[2*j+1]);
}


void get_dnu(const double Ax_max, const double Ay_max, const double delta_max)
{
  char            str[max_str];
  int             i, k;
  double          twoJ[2], nux, nuy;
  ss_vect<double> ps;
  ss_vect<tps>    Id_scl;
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
  Id_scl.zero(); ps.zero();
  for (i = -n_ampl; i <= n_ampl; i++) {
    ps[x_] = i*Ax_max/n_ampl;
    if (ps[x_] == 0.0) ps[x_] = A_min;
    ps[y_] = A_min;
    get_twoJ(ps, twoJ);
    for (k = 0; k < 4; k++)
      Id_scl[k] = sqrt(twoJ[k/2]);
    nux = (nus[3]*Id_scl).cst(); nuy = (nus[4]*Id_scl).cst();

    outf << scientific << setprecision(3)
	 << setw(12) << 1e3*ps[x_] << setw(12) << 1e3*ps[y_]
	 << fixed << setprecision(5)
	 << setw(9) << nux << setw(9) << nuy << endl;
  }
  outf.close();

//  sprintf(str, "%s%s", home_dir, "/projects/src/");
//  file_wr(outf, strcat(str, "dnu_dAy_pert.out"));
  file_wr(outf, "dnu_dAy_pert.out");
  Id_scl.zero(); ps.zero();
  for (i = -n_ampl; i <= n_ampl; i++) {
    ps[x_] = A_min;
    ps[y_] = i*Ay_max/n_ampl;
    if (ps[y_] == 0.0) ps[y_] = A_min;
//    get_twoJ(2, ps, MNF.A1, twoJ);
    get_twoJ(ps, twoJ);
    for (k = 0; k < 4; k++)
      Id_scl[k] = sqrt(twoJ[k/2]);
    nux = (nus[3]*Id_scl).cst(); nuy = (nus[4]*Id_scl).cst();

    outf << scientific << setprecision(3)
	 << setw(12) << 1e3*ps[x_] << setw(12) << 1e3*ps[y_]
	 << fixed << setprecision(6)
	 << setw(10) << nux << setw(10) << nuy << endl;
  }
  outf.close();

//  sprintf(str, "%s%s", home_dir, "/projects/src/");
//  file_wr(outf, strcat(str, "chrom2_pert.out"));
  file_wr(outf, "chrom2_pert.out");
  Id_scl.zero(); ps.zero();
  for (i = -n_delta; i <= n_delta; i++) {
    ps[delta_] = i*delta_max/n_delta; Id_scl[delta_] = ps[delta_];

    nux = (nus[3]*Id_scl).cst(); nuy = (nus[4]*Id_scl).cst();

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
  int             i, j, k;
  double          twoJ[2], nux, nuy;
  ss_vect<double> ps;
  ss_vect<tps>    Id_scl;
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
  Id_scl.zero(); Id_scl[delta_] = delta;
  ps.zero();
  for (i = -n_ampl; i <= n_ampl; i++) {
    for (j = -n_ampl; j <= n_ampl; j++) {
      ps[x_] = i*Ax_max/n_ampl; ps[y_] = j*Ay_max/n_ampl;
      get_twoJ(ps, twoJ);
      for (k = 0; k < 4; k++)
	Id_scl[k] = sqrt(twoJ[k/2]);
      nux = (nus[3]*Id_scl).cst(); nuy = (nus[4]*Id_scl).cst();

      outf << scientific << setprecision(3)
	   << setw(12) << 1e3*ps[x_] << setw(12) << 1e3*ps[y_]
	   << fixed << setprecision(5)
	   << " " << setw(8) << nux << " " << setw(8) << nuy << endl;
    }
    outf << endl;
  }
  outf.close();
}


void wtf()
{
  long int jj[ss_dim];
  int      k;
  tps      g_re, g_im, K_re, K_im;
  
  Ring_GetTwiss(true, 0.0); printglob();

  get_map(false);

  danot_(no_tps);
  get_map_normal_form(); nus = dHdJ(MNF.K);

  CtoR(MNF.K, K_re, K_im); CtoR(MNF.g, g_re, g_im);

  cout << scientific << setprecision(5) << K_re;
  for (k = 0; k < nv_tps; k++)
    jj[k] = 0;
  jj[x_] = 1; jj[px_] = 1; jj[delta_] = 1;
  printf("\nTune Shift Terms:\n  h_11001 = %12.5e\n", -2e0*K_re[jj]);
  jj[x_] = 0; jj[px_] = 0; jj[y_] = 1; jj[py_] = 1;
  printf("  h_00111 = %12.5e\n", -2e0*K_re[jj]);
  jj[x_] = 2; jj[px_] = 2; jj[y_] = 0; jj[py_] = 0; jj[delta_] = 0;
  printf("  h_22000 = %12.5e\n", -4e0*K_re[jj]);
  jj[x_] = 0; jj[px_] = 0; jj[y_] = 2; jj[py_] = 2;
  printf("  h_00220 = %12.5e\n", -4e0*K_re[jj]);
  jj[x_] = 1; jj[px_] = 1; jj[y_] = 1; jj[py_] = 1;
  printf("  h_11110 = %12.5e\n", -4e0*K_re[jj]);
}


int main(int argc, char *argv[])
{
  int          k;
  tps          H, H_re, H_im, g_re, g_im, K_re, K_im;
  ss_vect<tps> Id_scl;
  ofstream     outf;

  globval.H_exact    = false; globval.quad_fringe    = false;
  globval.Cavity_on  = false; globval.radiation      = false;
  globval.emittance  = false; globval.IBS            = false;
  globval.pathlength = false; globval.bpm            = 0;
  globval.Cart_Bend  = false; globval.dip_edge_fudge = true;
  globval.mat_meth   = false;

  // disable from TPSALib- and LieLib log messages
  idprset_(-1);

  if (false)
    Read_Lattice(argv[1]);
  else {
    rdmfile(argv[1]);
  }

  globval.EPU = true;

  if (false) {
    wtf();
    exit(0);
  }

  danot_(1);

//  Ring_GetTwiss(true, 0.0); printglob();

  prt_lat("linlat.out", globval.bpm, true);

  danot_(no_tps-1);

//  get_map_n(3);
  get_map(false);

  danot_(no_tps);
  get_map_normal_form(); nus = dHdJ(MNF.K);

  CtoR(MNF.K, K_re, K_im); CtoR(MNF.g, g_re, g_im);

  Id_scl.identity();
  for (k = 0; k < 4; k++)
    Id_scl[k] *= sqrt(twoJ[k/2]);
  Id_scl[delta_] *= delta_max;

  if (true) {
//    H = get_H(); CtoR(H, H_re, H_im);

//    outf.open("H.dat", ios::out);
//    outf << H_re*Id_scl;
//    outf.close();

    outf.open("K.dat", ios::out);
    outf << K_re*Id_scl;
    outf.close();

    outf.open("g.dat", ios::out);
    outf << g_im*Id_scl;
    outf.close();

//    cout << N*nus[3]*Id_scl << N*nus[4]*Id_scl;

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
