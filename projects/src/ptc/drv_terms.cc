#define NO 5

#include "tracy_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;


const double A_max[]    = {2e-2, 1e-2},
             delta_max  = 2e-2,
             beta_inj[] = {10.5, 5.2};


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


double h_ijklm_abs(const tps &h_re, const tps &h_im,
		   const int i, const int j, const int k, const int l,
		   const int m)
{
  return sqrt(sqr(h_ijklm(h_re, i, j, k, l, m))
	      +sqr(h_ijklm(h_im, i, j, k, l, m)));
}


tps get_h_local(const ss_vect<tps> &map)
{
  ss_vect<tps>  map1, R;

  // Parallel transport nonlinear kick to start of lattice,
  // assumes left to right evaluation.

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
  int          i;
  tps          h_re, h_im;
  ss_vect<tps> Id_scl, A0, A1;

  Id_scl.identity();
  for (i = 0; i < 4; i++)
    Id_scl[i] *= sqrt(twoJ[i/2]);
  Id_scl[delta_] *= delta;

  CtoR(get_h_local(map_Fl)*Id_scl, h_re, h_im);

  printf("%5d (%3ld)\n", k, globval.Cell_nLoc);
  outf << fixed << setprecision(3)
       << setw(6) << k << setw(9) << Cell[k].S
       << scientific << setprecision(5)

       << setw(13) << h_ijklm_abs(h_re, h_im, 1, 0, 0, 0, 2)

       << setw(13) << h_ijklm_abs(h_re, h_im, 2, 0, 0, 0, 1)
       << setw(13) << h_ijklm_abs(h_re, h_im, 0, 0, 2, 0, 1)

       << setw(13) << h_ijklm_abs(h_re, h_im, 3, 0, 0, 0, 0)
       << setw(13) << h_ijklm_abs(h_re, h_im, 2, 1, 0, 0, 0)
       << setw(13) << h_ijklm_abs(h_re, h_im, 1, 0, 2, 0, 0)
       << setw(13) << h_ijklm_abs(h_re, h_im, 1, 0, 0, 2, 0)
       << setw(13) << h_ijklm_abs(h_re, h_im, 1, 0, 1, 1, 0)

       << setw(13) << h_ijklm_abs(h_re, h_im, 4, 0, 0, 0, 0)
       << setw(13) << h_ijklm_abs(h_re, h_im, 3, 1, 0, 0, 0)
       << setw(13) << h_ijklm_abs(h_re, h_im, 2, 0, 2, 0, 0)
       << setw(13) << h_ijklm_abs(h_re, h_im, 1, 1, 2, 0, 0)
       << setw(13) << h_ijklm_abs(h_re, h_im, 0, 0, 4, 0, 0)
       << setw(13) << h_ijklm_abs(h_re, h_im, 2, 0, 1, 1, 0)
       << setw(13) << h_ijklm_abs(h_re, h_im, 0, 0, 3, 1, 0)
       << setw(13) << h_ijklm_abs(h_re, h_im, 2, 0, 0, 2, 0)

       << setw(13) << h_ijklm_abs(h_re, h_im, 5, 0, 0, 0, 0)
       << setw(13) << h_ijklm_abs(h_re, h_im, 4, 1, 0, 0, 0)
       << setw(13) << h_ijklm_abs(h_re, h_im, 3, 2, 0, 0, 0)
       << setw(13) << h_ijklm_abs(h_re, h_im, 3, 0, 2, 0, 0)
       << setw(13) << h_ijklm_abs(h_re, h_im, 2, 1, 2, 0, 0)
       << setw(13) << h_ijklm_abs(h_re, h_im, 1, 0, 4, 0, 0)
       << setw(13) << h_ijklm_abs(h_re, h_im, 3, 0, 1, 1, 0)
       << setw(13) << h_ijklm_abs(h_re, h_im, 2, 1, 1, 1, 0)
       << setw(13) << h_ijklm_abs(h_re, h_im, 1, 0, 3, 1, 0)
       << setw(13) << h_ijklm_abs(h_re, h_im, 3, 0, 0, 2, 0)
       << setw(13) << h_ijklm_abs(h_re, h_im, 2, 1, 0, 2, 0)
       << setw(13) << h_ijklm_abs(h_re, h_im, 1, 0, 2, 2, 0)
       << setw(13) << h_ijklm_abs(h_re, h_im, 1, 0, 1, 3, 0)
       << setw(13) << h_ijklm_abs(h_re, h_im, 1, 0, 0, 4, 0)
       << setw(13) << h_ijklm_abs(h_re, h_im, 0, 1, 0, 4, 0)

       << "\n";

  outf.flush();
}


void get_drv_terms(const double twoJ[], const double delta)
{
  long int     lastpos;
  int          k;
  double       alpha[2], beta[2], eta[2], etap[2];
  ss_vect<tps> map, A0, A1, A_Atp;
  ofstream     outf;

  outf.open("drv_terms.out", ios::out);
  printf("\n");
  map.identity(); putlinmat(6, globval.Ascr, A0); A1 = A0;
  for (k = 0; k <= globval.Cell_nLoc; k++) {
    Cell_Pass(k, k, map, lastpos); Cell_Pass(k, k, A1, lastpos);
    A1 = get_A(Cell[k].Alpha, Cell[k].Beta, Cell[k].Eta, Cell[k].Etap);
    prt_drv_terms(outf, k, twoJ, delta, Inv(A1)*map*A0);
  }
  outf.close();
}


int main(int argc, char *argv[])
{
  int    j;
  double twoJ[2];

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

  for (j = 0; j < 2; j++)
    twoJ[j] = sqr(A_max[j])/beta_inj[j];
  get_drv_terms(twoJ, delta_max);
}
