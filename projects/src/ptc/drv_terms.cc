#define NO 3

#include "tracy_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;


const bool   set_dnu = false;
const int    n_cell  = 2;
const double
  beta_inj[] = {7.9, 3.1},
  // A_max[]    = {3e-3, 1.5e-3},
  A_max[]    = {1e-3, 0e-3},
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
  s = Cell[loc].S;
  if (k > globval.Cell_nLoc) s += Cell[globval.Cell_nLoc].S;

  printf("%5d (%3ld)\n", k, n_cell*globval.Cell_nLoc);
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
  double       dnu[2];
  ss_vect<tps> Id, map, map_Fl, A0, A1;
  ofstream     outf;

  const double eta0[]  = {0e0, 0e0},
               etap0[] = {0e0, 0e0};
  
  outf.open("drv_terms.out", ios::out);

  Id.identity();

  // Get linear dispersion.
  danot_(2);
  map.identity();
  Cell_Pass(0, globval.Cell_nLoc, map, lastpos);
  MNF = MapNorm(map, 1);
  danot_(no_tps);

  A0 = MNF.A1; A1 = A0;
  map.identity();
  map[x_] += MNF.A0[x_][delta_]*Id[delta_];
  map[px_] += MNF.A0[px_][delta_]*Id[delta_];
  for (k = 0; k <= n_cell*globval.Cell_nLoc; k++) {
    loc = k % (globval.Cell_nLoc+1);
    Elem_Pass(loc, map);
    danot_(1);
    Elem_Pass(loc, A1);
    danot_(no_tps);
    A1 = get_A_CS(2, A1, dnu);
    map_Fl = map;
    map_Fl[x_] -= map[x_][delta_]*Id[delta_];
    map_Fl[px_] -= map[px_][delta_]*Id[delta_];
    A1 = get_A(Cell[loc].Alpha, Cell[loc].Beta, eta0, etap0);
    map_Fl = Inv(A1)*map_Fl*A0;
    prt_drv_terms(outf, k, twoJ, delta, map_Fl);
  }

  outf.close();

  // cout << scientific << setprecision(3) << MNF.A0 << "\n";
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


void get_ampl_orb(const double twoJ[])
{
  long int        lastpos;
  int             j, k;
  ss_vect<double> Id_scl;
  ss_vect<tps>    Id, dx, dx_fl, dx_re, dx_im, M;
  ofstream        outf;

  outf.open("ampl_orb.out", ios::out);

  Id.identity();

  Id_scl.zero();
  for (j = 0; j < 4; j++)
    Id_scl[j] = sqrt(twoJ[j/2]);

  danot_(no_tps-1);
  map.identity();
  Cell_Pass(0, globval.Cell_nLoc, map, lastpos);

  M.identity();
  for (j = 0; j <= globval.Cell_nLoc; j++) {
    danot_(no_tps-1);
    Elem_Pass(j, M);
    danot_(no_tps);

    if (false || ((Cell[j].Elem.Pkind == Mpole) &&
		   (Cell[j].Elem.M->PBpar[Sext+HOMmax] != 0e0))) {
      MNF = MapNorm(M*map*Inv(M), 1);
#if 0
      dx_fl = LieExp(MNF.g, Id);
#else
      for (k = 0; k < 4; k++)
	dx_fl[k] = PB(MNF.g, Id[k]);
#endif
      for (k = 0; k < 4; k++) {
	CtoR(dx_fl[k], dx_re[k], dx_im[k]);
	dx_re[k] = dacfu1(dx_re[k], f_kernel);
      }
      dx_re = MNF.A1*dx_re; dx_im = MNF.A1*dx_im;
      outf << setw(4) << j << fixed << setprecision(3) << setw(8) << Cell[j].S
	   << " " << setw(8) << Cell[j].Elem.PName;
      for (k = 0; k < 4; k++)
	outf << scientific << setprecision(5)
	     << setw(13) << (dx_re[k]*Id_scl).cst();
      outf << "\n";
    }
  }

  outf.close();
}


int main(int argc, char *argv[])
{

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  // disable from TPSALib and LieLib log messages
  idprset(-1);

  std::string home_dir = "";

  daeps_(eps_tps);

  if (!true)
    Read_Lattice((home_dir+argv[1]).c_str());
  else
    rdmfile(argv[1]);

  if (false) chk_bend();

  Ring_GetTwiss(true, 0e0); printglob();

  if (set_dnu) {
    set_map(ElemIndex("ps_rot"), dnu);
    Ring_GetTwiss(true, 0e0); printglob();
  }

  if (false) {
    tst_g();
    exit(0);
  }

  if (!true) get_drv_terms(twoJ, delta_max);

  if (!false) get_ampl_orb(twoJ);
}
