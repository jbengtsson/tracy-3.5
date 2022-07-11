#define NO 5

#include "tracy_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;


const double
#if 1
  beta_inj[] = {3.1, 3.0},
#else
  beta_inj[] = {3.1, 2.8},
#endif
  A_max[]    = {3e-3, 1.5e-3},
  delta_max  = 2e-2,
  twoJ[]     = {sqr(A_max[X_])/beta_inj[X_], sqr(A_max[Y_])/beta_inj[Y_]},
  dnu[]      = {0.03, 0.02};


ss_vect<tps> get_map_Fl(ss_vect<tps> &map)
{
  Matrix       R_mat;
  ss_vect<tps> Id, A0_inv, R;

  const int  n_dim = 4;

  Id.identity();

  // Find fixed point.
  GoFix(map, MNF.A0, A0_inv, no_tps);

  // Translate to fix point.
  map = A0_inv*map*MNF.A0;

  getlinmat(n_dim, map, globval.OneTurnMat);
  GDiag(n_dim, 1.0, globval.Ascr, globval.Ascrinv, R_mat,
        globval.OneTurnMat, globval.Omega, globval.Alphac);

  MNF.A1 = putlinmat(n_dim, globval.Ascr);
  MNF.A1[ct_] = Id[ct_];
  MNF.A1[delta_] = Id[delta_];

  R = putlinmat(n_dim, R_mat);
  R[ct_] = Id[ct_];
  R[delta_] = Id[delta_];

  // Transform to Floquet space.
  return Inv(MNF.A1)*map*MNF.A1;
}


ss_vect<tps> get_map_Fl_2(ss_vect<tps> &map)
{
  MNF = MapNorm(map, 1);
  return Inv(MNF.A0*MNF.A1)*map*MNF.A0*MNF.A1;
}


tps get_mns(const tps &a, const int no1, const int no2)
{
  tps b;

  danot_(no1-1);
  b = -a;
  danot_(no2);
  b += a;
  danot_(no_tps);
  return b;
}


ss_vect<tps> get_mns(const ss_vect<tps> &x, const int no1, const int no2)
{
  int          k;
  ss_vect<tps> y;

  for (k = 0; k < nv_tps; k++)
    y[k] = get_mns(x[k], no1, no2);

  return y;
}


void Dragt_Finn_Fact(const ss_vect<tps> &map, ss_vect<tps> &R, tps &h)
{
  /* Dragt-Finn factorization:

       M = exp(:h_no:)*...*exp(:h_4:)*exp(:h_3:)*R                            */

  int          k;
  tps          hn;
  ss_vect<tps> Id, A0_inv, Fn, map1;

  Id.identity();

  danot_(1);
  map1 = map;
  R = get_map_Fl(map1);

  danot_(no_tps);
  map1 = map*Inv(R);
  h = 0e0;
  for (k = 3; k <= no_tps; k++) {
    Fn = get_mns(map1, k-1, k-1);
    hn = Intd(Fn, -1.0);
    h += hn;
    map1 = map1*LieExp(-hn, Id);
  }
}


void Dragt_Finn_Fact_2(ss_vect<tps> &map, ss_vect<tps> &M, tps &h)
{
  h = LieFact_DF(map, M);
}


ss_vect<tps> Dragt_Finn_Map(ss_vect<tps> &M, tps &h)
{
  int           k;
  ss_vect<tps>  Id, map;

  Id.identity();
  map = M;
  for (k = 3; k <= no_tps; k++)
    map = LieExp(Take(h, k), Id)*map;
  return map;
}


ss_vect<tps> Dragt_Finn_Map_2(ss_vect<tps> &M, tps &h)
{
  ss_vect<tps> Id;

  Id.identity();
  return FExpo(h, Id, 3, no_tps, 1)*M;
}


void analyse()
{
  tps          h, h_2;
  ss_vect<tps> map_Fl, map_Fl_2, R, R_2;

  danot_(no_tps-1);
  get_map(false);
  danot_(no_tps);

  printf("\nM:\n");
  prt_lin_map(3, map);

  map_Fl = get_map_Fl(map);

  if (!false) {
    map_Fl_2 = get_map_Fl_2(map);
    daeps_(1e-6);
    cout << scientific << setprecision(3) << setw(11) << map_Fl-map_Fl_2
	 << "\n";
    daeps_(eps_tps);
    printf("\nmap_Fl:\n");
    prt_lin_map(3, map_Fl);
  }

  Dragt_Finn_Fact(map_Fl, R, h);
  if (!false) {
    Dragt_Finn_Fact_2(map_Fl_2, R_2, h_2);
    daeps_(1e-6);
    cout << scientific << setprecision(3) << "\nR-R_2\n"
	 << setw(11) << R-R_2 << "\n";
    cout << scientific << setprecision(3) << "\nh-h_2\n"
	 << setw(11) << h-h_2 << "\n";
    daeps_(eps_tps);
    printf("\nmap_Fl:\n");
    prt_lin_map(3, map_Fl);
  }

  map_Fl = Dragt_Finn_Map(R, h);

  if (!false) {
    map_Fl_2 = Dragt_Finn_Map(R, h);
    danot_(no_tps-1);
    daeps_(1e-6);
    map_Fl_2 = map_Fl_2;
    cout << scientific << setprecision(3) << setw(11) << map_Fl-map_Fl_2
	 << "\n";
    daeps_(eps_tps);
  }
}


int main(int argc, char *argv[])
{
  int          k;
  ss_vect<tps> Id_scl;

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

  Id_scl.identity();
  for (k = 0; k < 4; k++)
    Id_scl[k] *= sqrt(twoJ[k/2]);
  Id_scl[delta_] *= delta_max;

  if (!true)
    Read_Lattice((home_dir+argv[1]).c_str());
  else
    rdmfile(argv[1]);

  analyse();
}
