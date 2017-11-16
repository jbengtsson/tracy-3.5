#define NO 5

#include "tracy_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;


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


void prt_drv_terms(ofstream &outf, const int k,
		   const double twoJx, const double twoJy, const double delta)
{
  tps          h_re, h_im;
  ss_vect<tps> Id_scl, nus;

  Id_scl.identity();
  Id_scl[x_] *= sqrt(twoJx); Id_scl[px_] *= sqrt(twoJx);
  Id_scl[y_] *= sqrt(twoJy); Id_scl[py_] *= sqrt(twoJy);
  Id_scl[delta_] *= delta;

  MNF = MapNorm(map, no_tps);
  CtoR(get_h(), h_re, h_im); h_re = h_re*Id_scl; h_im = h_im*Id_scl;
  nus = dHdJ(MNF.K);
  printf("%5d (%3ld) nu = [%7.5f %7.5f]\n",
	 k, globval.Cell_nLoc, nus[0].cst(), nus[1].cst());
  outf << fixed << setprecision(3)
       << setw(6) << k << setw(9) << Cell[k].S
       << scientific << setprecision(5)

       << setw(13) << h_ijklm(h_re, 2, 0, 0, 0, 1)
       << setw(13) << h_ijklm(h_re, 0, 0, 2, 0, 1)

       << setw(13) << h_ijklm(h_re, 3, 0, 0, 0, 0)
       << setw(13) << h_ijklm(h_re, 2, 1, 0, 0, 0)
       << setw(13) << h_ijklm(h_re, 1, 0, 2, 0, 0)
       << setw(13) << h_ijklm(h_re, 1, 0, 0, 2, 0)
       << setw(13) << h_ijklm(h_re, 1, 0, 1, 1, 0)

       << setw(13) << h_ijklm(h_re, 4, 0, 0, 0, 0)
       << setw(13) << h_ijklm(h_re, 3, 1, 0, 0, 0)
       << setw(13) << h_ijklm(h_re, 2, 0, 2, 0, 0)
       << setw(13) << h_ijklm(h_re, 1, 1, 2, 0, 0)
       << setw(13) << h_ijklm(h_re, 0, 0, 4, 0, 0)
       << setw(13) << h_ijklm(h_re, 2, 0, 1, 1, 0)
       << setw(13) << h_ijklm(h_re, 0, 0, 3, 1, 0)
       << setw(13) << h_ijklm(h_re, 2, 0, 0, 2, 0)

       << setw(13) << h_ijklm(h_re, 5, 0, 0, 0, 0)
       << setw(13) << h_ijklm(h_re, 4, 1, 0, 0, 0)
       << setw(13) << h_ijklm(h_re, 3, 2, 0, 0, 0)
       << setw(13) << h_ijklm(h_re, 3, 0, 2, 0, 0)
       << setw(13) << h_ijklm(h_re, 2, 1, 2, 0, 0)
       << setw(13) << h_ijklm(h_re, 1, 0, 4, 0, 0)
       << setw(13) << h_ijklm(h_re, 3, 0, 1, 1, 0)
       << setw(13) << h_ijklm(h_re, 2, 1, 1, 1, 0)
       << setw(13) << h_ijklm(h_re, 1, 0, 3, 1, 0)
       << setw(13) << h_ijklm(h_re, 3, 0, 0, 2, 0)
       << setw(13) << h_ijklm(h_re, 2, 1, 0, 2, 0)
       << setw(13) << h_ijklm(h_re, 1, 0, 2, 2, 0)
       << setw(13) << h_ijklm(h_re, 1, 0, 1, 3, 0)
       << setw(13) << h_ijklm(h_re, 1, 0, 0, 4, 0)
       << setw(13) << h_ijklm(h_re, 0, 1, 0, 4, 0)

       << "\n";
}


void get_drv_terms(const double twoJx, const double twoJy, const double delta)
{
  long int lastpos;
  int      k;
  ofstream outf;

  const int n_step = 1;

  outf.open("drv_terms.out", ios::out);
  printf("\n");
  k = 0;
  while (k <= globval.Cell_nLoc) {
    map.identity();
    Cell_Pass(k, globval.Cell_nLoc, map, lastpos);
    if (k != 0) Cell_Pass(0, k-1, map, lastpos);
    prt_drv_terms(outf, k, twoJx, twoJy, delta);
    k += n_step;
  }
  if (k != globval.Cell_nLoc)
    prt_drv_terms(outf, globval.Cell_nLoc, twoJx, twoJy, delta);
  outf.close();
}


int main(int argc, char *argv[])
{
  int    j;
  double twoJ[2];

  const double A_max[]    = {2e-2, 1e-2},
               delta_max  = 2e-2,
               beta_inj[] = {10.5, 5.2};

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  // disable from TPSALib and LieLib log messages
  idprset(-1);

  std::string home_dir = "";

  if (true)
    Read_Lattice((home_dir+argv[1]).c_str());
  else
    rdmfile(argv[1]);

  if (false) chk_bend();

  Ring_GetTwiss(true, 0e0); printglob();

  for (j = 0; j < 2; j++)
    twoJ[j] = sqr(A_max[j])/beta_inj[j];
  get_drv_terms(twoJ[X_], twoJ[Y_], delta_max);
}
