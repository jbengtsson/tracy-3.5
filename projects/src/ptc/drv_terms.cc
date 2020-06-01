#define NO 3

#include "tracy_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;


const bool   set_dnu = false;
const int    n_cell  = 2;
const double
  beta_inj[] = {7.9, 3.1},
  A_max[]    = {3e-3, 1.5e-3},
  delta_max  = 2e-2,
  twoJ[]     = {sqr(A_max[X_])/beta_inj[X_], sqr(A_max[Y_])/beta_inj[Y_]},
  dnu[]      = {0.03, 0.02};


void chk_bend()
{
  int    k;
  double phi;

  phi = 0e0;
  for (k = 0; k <= globval.Cell_nLoc; k++) {
    if ((Cell[k].Elem.Pkind == Mpole) &&
	(Cell[k].Elem.M->n_design == Dip)) {
      phi += Cell[k].Elem.PL*Cell[k].Elem.M->Pirho;
    }
  }
  phi = 180e0*phi/M_PI;
  printf("\nphi = %8.6f\n", phi);
}


double h_abs_ijklm(const tps &h_re, const tps &h_im,
		   const int i, const int j, const int k, const int l,
		   const int m)
{
  return sqrt(sqr(h_ijklm(h_re, i, j, k, l, m))
  	      +sqr(h_ijklm(h_im, i, j, k, l, m)));
}


tps get_h_local(const ss_vect<tps> &map)
{
  ss_vect<tps>  map1, R;

  if (true)
    // Dragt-Finn factorization.
    return LieFact_DF(map, R);
  else {
    // Single Lie exponent.
    danot_(1); map1 = map; danot_(no_tps);
    return LieFact(map*Inv(map1));
  }
}


void prt_drv_terms(ofstream &outf, const int k,
		   const double twoJ[], const double delta,
		   const ss_vect<tps> &map_Fl)
{
  int          i, loc;
  double       s;
  tps          h_re, h_im;
  ss_vect<tps> Id_scl, A0, A1;

  Id_scl.identity();
  for (i = 0; i < 4; i++)
    Id_scl[i] *= sqrt(twoJ[i/2]);
  Id_scl[delta_] *= delta;

  CtoR(get_h_local(map_Fl)*Id_scl, h_re, h_im);

  loc = k % (globval.Cell_nLoc+1);
  s = Cell[loc].S + k/(globval.Cell_nLoc+1)*Cell[globval.Cell_nLoc].S;

  outf << fixed << setw(3) << k
       << setprecision(3) << setw(9) << s
       << setprecision(1) << setw(5) << get_code(Cell[loc])
       << scientific << setprecision(5)

       << setw(13) << h_abs_ijklm(h_re, h_im, 1, 0, 0, 0, 2)

       << setw(13) << h_abs_ijklm(h_re, h_im, 2, 0, 0, 0, 1)
       << setw(13) << h_abs_ijklm(h_re, h_im, 0, 0, 2, 0, 1)

       << setw(13) << h_abs_ijklm(h_re, h_im, 3, 0, 0, 0, 0)
       << setw(13) << h_abs_ijklm(h_re, h_im, 2, 1, 0, 0, 0)
       << setw(13) << h_abs_ijklm(h_re, h_im, 1, 0, 2, 0, 0)
       << setw(13) << h_abs_ijklm(h_re, h_im, 1, 0, 0, 2, 0)
       << setw(13) << h_abs_ijklm(h_re, h_im, 1, 0, 1, 1, 0)

       << setw(13) << h_abs_ijklm(h_re, h_im, 4, 0, 0, 0, 0)
       << setw(13) << h_abs_ijklm(h_re, h_im, 3, 1, 0, 0, 0)
       << setw(13) << h_abs_ijklm(h_re, h_im, 2, 0, 2, 0, 0)
       << setw(13) << h_abs_ijklm(h_re, h_im, 1, 1, 2, 0, 0)
       << setw(13) << h_abs_ijklm(h_re, h_im, 0, 0, 4, 0, 0)
       << setw(13) << h_abs_ijklm(h_re, h_im, 2, 0, 1, 1, 0)
       << setw(13) << h_abs_ijklm(h_re, h_im, 0, 0, 3, 1, 0)
       << setw(13) << h_abs_ijklm(h_re, h_im, 2, 0, 0, 2, 0)

       << setw(13) << h_abs_ijklm(h_re, h_im, 5, 0, 0, 0, 0)
       << setw(13) << h_abs_ijklm(h_re, h_im, 4, 1, 0, 0, 0)
       << setw(13) << h_abs_ijklm(h_re, h_im, 3, 2, 0, 0, 0)
       << setw(13) << h_abs_ijklm(h_re, h_im, 3, 0, 2, 0, 0)
       << setw(13) << h_abs_ijklm(h_re, h_im, 2, 1, 2, 0, 0)
       << setw(13) << h_abs_ijklm(h_re, h_im, 1, 0, 4, 0, 0)
       << setw(13) << h_abs_ijklm(h_re, h_im, 3, 0, 1, 1, 0)
       << setw(13) << h_abs_ijklm(h_re, h_im, 1, 0, 4, 0, 0)
       << setw(13) << h_abs_ijklm(h_re, h_im, 2, 1, 1, 1, 0)
       << setw(13) << h_abs_ijklm(h_re, h_im, 1, 0, 3, 1, 0)
       << setw(13) << h_abs_ijklm(h_re, h_im, 3, 0, 0, 2, 0)
       << setw(13) << h_abs_ijklm(h_re, h_im, 2, 1, 0, 2, 0)
       << setw(13) << h_abs_ijklm(h_re, h_im, 1, 0, 2, 2, 0)
       << setw(13) << h_abs_ijklm(h_re, h_im, 1, 0, 1, 3, 0)
       << setw(13) << h_abs_ijklm(h_re, h_im, 1, 0, 0, 4, 0)
       << setw(13) << h_abs_ijklm(h_re, h_im, 0, 1, 0, 4, 0)
       << "\n";

  outf.flush();
}


void get_drv_terms(const double twoJ[], const double delta)
{
  long int     lastpos;
  int          k, loc;
  ss_vect<tps> Id, map, map_Fl, A_k;
  ofstream     outf;

  const string str = "drv_terms.out";
  
  outf.open(str.c_str(), ios::out);

  Id.identity();

  danot_(2);
  map.identity();
  Cell_Pass(0, globval.Cell_nLoc, map, lastpos);
  MNF = MapNorm(map, 1);

  A_k = MNF.A0*MNF.A1;
  prt_lin_map(3, map);
  prt_lin_map(3, Inv(A_k)*map*A_k);

  map.identity();
  map[x_] += MNF.A0[x_][delta_]*Id[delta_];
  map[px_] += MNF.A0[px_][delta_]*Id[delta_];
  for (k = 0; k < n_cell*(globval.Cell_nLoc+1); k++) {
    loc = k % (globval.Cell_nLoc+1);
    printf("%5d (%3ld)\n", k, n_cell*globval.Cell_nLoc);
    danot_(1);
    Elem_Pass(loc, A_k);
    danot_(no_tps);
    Elem_Pass(loc, map);
    map_Fl = Inv(A_k)*map;
    prt_drv_terms(outf, k, twoJ, delta, map_Fl);
  }

  outf.close();
}


void tst_g()
{
  long int     lastpos;
  ss_vect<tps> Id, Id_scl, dx_fl, M, M_fl, A1;

  const int loc = 30;

  daeps_(1e-10);

  Id.identity();

  Id_scl.identity();
  Id_scl[delta_] = 0e0; Id_scl[ct_] = 0e0;

  danot_(no_tps-1);
  map.identity();
  Cell_Pass(loc+1, globval.Cell_nLoc, map, lastpos);
  Cell_Pass(0, loc, map, lastpos);
  danot_(no_tps);
  MNF = MapNorm(map, 1);
  // dx_fl = LieExp(MNF.g, Id);
  // cout << scientific << setprecision(5) << setw(13) << (dx_fl*Id_scl)[x_];
  cout << scientific << setprecision(5) << setw(13) << (MNF.g)*Id_scl;

  danot_(no_tps-1);
  map.identity(); Cell_Pass(0, globval.Cell_nLoc, map, lastpos);
  M.identity(); Cell_Pass(0, loc, M, lastpos);
  danot_(no_tps);
  MNF = MapNorm(map, 1);
  A1 = get_A(Cell[loc].Alpha, Cell[loc].Beta, Cell[loc].Eta, Cell[loc].Etap);
  M_fl = Inv(A1)*M*MNF.A1;
  // dx_fl = M_fl*LieExp(MNF.g, Id)*Inv(M_fl);
  dx_fl = LieExp(MNF.g*Inv(M_fl), Id);
  // cout << scientific << setprecision(5) << setw(13) << (dx_fl*Id_scl)[x_];
  cout << scientific << setprecision(5) << setw(13) << (MNF.g*Inv(M_fl))*Id_scl;
}


tps dacfu1(const tps &a, double (*func)(const long int []))
{
  char    name[11];
  int     j, n;
  long int jj[ss_dim], ibuf1[bufsize], ibuf2[bufsize];
  double  rbuf[bufsize];
  tps     b;

  a.exprt(rbuf, ibuf1, ibuf2, name); n = (int)rbuf[0];

  for (j = 0; j < n; j++) {
    dehash_(no_tps, ss_dim, ibuf1[j], ibuf2[j], jj);

    rbuf[j+1] *= (*func)(jj);
  }

  b.imprt(n, rbuf, ibuf1, ibuf2);

  // Remove zeroes.
  return 1e0*b;
}


double f_kernel(const long int jj[])
{

  return
    (((jj[x_] != 0) || (jj[y_] != 0)) &&
     (jj[x_] == jj[px_]) && (jj[y_] == jj[py_]))?
    1e0 : 0e0;
}


void get_dx_dJ(const double twoJ[])
{
  long int     lastpos;
  int          j, k;
  double       b3L, a3L, dx[2];
  tps          g_re, g_im;
  ss_vect<tps> Id, Id_scl, dx_fl, dx_re, dx_im, M;
  ofstream     outf;

  outf.open("dx_dJ.out", ios::out);

  Id.identity();

  Id_scl.identity();
  for (j = 0; j < 4; j++)
    Id_scl[j] *= sqrt(twoJ[j/2]);
  Id_scl[delta_] = 0e0;

  danot_(no_tps-1);
  map.identity();
  Cell_Pass(0, globval.Cell_nLoc, map, lastpos);

  M.identity();
  for (j = 0; j <= globval.Cell_nLoc; j++) {
    danot_(no_tps-1);
    Elem_Pass(j, M);
    danot_(no_tps);

    if (!false || (Cell[j].Elem.Pkind == Mpole)) {
      MNF = MapNorm(M*map*Inv(M), 1);
#if 0
      dx_fl = LieExp(MNF.g, Id);
#else
      for (k = 0; k < 4; k++)
	dx_fl[k] = PB(MNF.g, Id[k]);
#endif
      for (k = 0; k < 4; k++)
	CtoR(dx_fl[k], dx_re[k], dx_im[k]);
#if 0
      // Phase space.
      get_bnL_design_elem(Cell[j].Fnum, Cell[j].Knum, Sext, b3L, a3L);
      dx_re = MNF.A1*dx_re; dx_im = MNF.A1*dx_im;
      dx[X_] = b3L*h_ijklm(dx_re[x_]*Id_scl, 1, 1, 0, 0, 0);
      dx[Y_] = b3L*h_ijklm(dx_re[x_]*Id_scl, 0, 0, 1, 1, 0);
#else
      // Floquet space.
      dx[X_] = h_ijklm(dx_re[x_]*Id_scl, 1, 1, 0, 0, 0);
      dx[Y_] = h_ijklm(dx_re[x_]*Id_scl, 0, 0, 1, 1, 0);
#endif
      outf << setw(4) << j << fixed << setprecision(3) << setw(8) << Cell[j].S
	   << " " << setw(8) << Cell[j].Elem.PName;
      for (k = 0; k < 2; k++)
      	outf << scientific << setprecision(5) << setw(13) << dx[k];
      outf << "\n";
    }
  }

  outf.close();
}


ss_vect<tps> get_map(const MNF_struct MNF)
{
  ss_vect<tps> Id, map;

  Id.identity();

  map =
    MNF.A0*MNF.A1*FExpo(MNF.g, Id, 3, no_tps, -1)
    *FExpo(MNF.K, Id, 2, no_tps, -1)
    *Inv(MNF.A0*MNF.A1*FExpo(MNF.g, Id, 3, no_tps, -1));

  return map;
}


void zero_res(tps &g)
{
  long int jj[ss_dim];
  int      k;

  for (k = 0; k < ss_dim; k++)
    jj[k] = 0;
  jj[x_] = 2; jj[px_] = 1;
  g.pook(jj, 0e0);
  jj[x_] = 1; jj[px_] = 2;
  g.pook(jj, 0e0);
  jj[x_] = 1; jj[px_] = 0; jj[y_] = 1; jj[py_] = 1;
  g.pook(jj, 0e0);
  jj[x_] = 0; jj[px_] = 1;
  g.pook(jj, 0e0);
}


void map_gymn(void)
{
  int          k;
  tps          g_re, g_im;
  ss_vect<tps> Id, map1, map_res, dx_fl, dx_re, dx_im;

  Id.identity();

  danot_(no_tps-1);
  get_map(false);
  danot_(no_tps);
  MNF = MapNorm(map, 1);

  CtoR(MNF.g, g_re, g_im);
  daeps_(1e-8);
  cout << scientific << setprecision(5) << 1e0*g_im;
  if (true) zero_res(g_im);
  cout << scientific << setprecision(5) << 1e0*g_im;
  MNF.g = RtoC(g_re, g_im);

#if 0
  dx_fl = LieExp(MNF.g, Id);
#else
  for (k = 0; k < 4; k++)
    dx_fl[k] = PB(MNF.g, Id[k]);
#endif
  for (k = 0; k < 4; k++)
    CtoR(dx_fl[k], dx_re[k], dx_im[k]);
  dx_re = MNF.A1*dx_re; dx_im = MNF.A1*dx_im;
  printf("\ndx(J) = %10.3e %10.3e %10.3e %10.3e\n",
  	 h_ijklm(dx_re[x_], 1, 1, 0, 0, 0),
  	 h_ijklm(dx_re[x_], 0, 0, 1, 1, 0),
  	 h_ijklm(dx_im[x_], 1, 1, 0, 0, 0),
  	 h_ijklm(dx_im[x_], 0, 0, 1, 1, 0));
  exit(0);

#if 1
  daeps_(eps_tps);
  map_res =
    Inv(FExpo(MNF.K, Id, 2, no_tps, -1))
    *Inv(MNF.A0*MNF.A1*FExpo(MNF.g, Id, 3, no_tps, -1))*map
    *MNF.A0*MNF.A1*FExpo(MNF.g, Id, 3, no_tps, -1);
  daeps_(1e-8);
  danot_(no_tps-1);
  cout << scientific << setprecision(5)
       << 1e0*MNF.map_res[x_] << 1e0*map_res[x_];
#else
  map_res =
    Inv(MNF.A0*MNF.A1*FExpo(MNF.g, Id, 3, no_tps, -1))*map
    *MNF.A0*MNF.A1*FExpo(MNF.g, Id, 3, no_tps, -1);
  daeps_(1e-8);
  danot_(no_tps-1);
  cout << scientific << setprecision(5) << (map_res-MNF.map_res)[x_];
#endif
}


int main(int argc, char *argv[])
{

  globval.H_exact    = false; globval.quad_fringe    = false;
  globval.Cavity_on  = false; globval.radiation      = false;
  globval.emittance  = false; globval.IBS            = false;
  globval.pathlength = false; globval.bpm            = 0;
  globval.Cart_Bend  = false; globval.dip_edge_fudge = true;
  globval.mat_meth   = false;

  // disable from TPSALib and LieLib log messages
  idprset(-1);

  std::string home_dir = "";

  daeps_(eps_tps);

  if (!true)
    Read_Lattice((home_dir+argv[1]).c_str());
  else
    rdmfile(argv[1]);

  if (false) chk_bend();

  // Ring_GetTwiss(true, 0e0); printglob();

  if (set_dnu) {
    set_map(ElemIndex("ps_rot"), dnu);
    Ring_GetTwiss(true, 0e0); printglob();
  }

  if (false) {
    map_gymn();
    exit(0);
  }

  if (false) {
    tst_g();
    exit(0);
  }

  if (true) get_drv_terms(twoJ, delta_max);

  if (false) get_dx_dJ(twoJ);
}
