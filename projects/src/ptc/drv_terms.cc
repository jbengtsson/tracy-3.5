#define NO 3

#include "tracy_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;


const bool   set_dnu = false;
const int    n_cell  = 1;
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


void get_ampl_orb(const double twoJ[])
{
  long int     lastpos;
  int          j, k;
  ss_vect<tps> Id, Id_scl, dx_fl, dx_fl_lin, M, map1;
  ofstream     outf;

  outf.open("ampl_orb.out", ios::out);

  Id.identity();

  Id_scl.identity();
  for (k = 0; k < 4; k++)
    Id_scl[k] *= sqrt(twoJ[k/2]);
  Id_scl[delta_] = 0e0;

  map.identity();
  danot_(no_tps-1);
  Cell_Pass(0, globval.Cell_nLoc, map, lastpos);
  M.identity();
  for (j = 0; j <= globval.Cell_nLoc; j++) {
    danot_(no_tps-1);
    Elem_Pass(j, M);
    if (!false ||
  	((Cell[j].Elem.Pkind == Mpole) &&
  	 (Cell[j].Elem.M->PBpar[Sext+HOMmax] != 0e0))) {
      map1 = M*map*Inv(M);
      danot_(no_tps);
      MNF = MapNorm(map1, 1);
      dx_fl = LieExp(MNF.g, Id);
      // Transpose before removing linear terms.
      dx_fl = tp_S(3, dx_fl);
      // Remove linear terms.
      danot_(1);
      dx_fl_lin = dx_fl;
      danot_(no_tps);
      dx_fl = (dx_fl-dx_fl_lin)*Id_scl;
      outf << setw(4) << j << fixed << setprecision(3) << setw(8) << Cell[j].S
	   << " " << setw(8) << Cell[j].Elem.PName;
      for (k = 0; k < 4; k++)
  	outf << scientific << setprecision(5) << setw(13)
  	     << sqrt(abs2(dx_fl[k]));
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

  // get_drv_terms(twoJ, delta_max);

  if (false) {
    tst_g();
    exit(0);
  }

  get_ampl_orb(twoJ);
}
