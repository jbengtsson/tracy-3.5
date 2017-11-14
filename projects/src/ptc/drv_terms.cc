#define NO 4

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


void get_drv_terms()
{
  long int     lastpos;
  int          k;
  tps          h_re, h_im;
  ss_vect<tps> nus;
  ofstream     outf;

  outf.open("drv_terms.out", ios::out);
  printf("\n");
  for (k = 0; k <= globval.Cell_nLoc; k += 10) {
    map.identity();
    if (k == 0)
      Cell_Pass(0, globval.Cell_nLoc, map, lastpos);
    else {
      Cell_Pass(k, globval.Cell_nLoc, map, lastpos);
      Cell_Pass(0, k-1, map, lastpos);
    }
    MNF = MapNorm(map, no_tps);
    CtoR(get_h(), h_re, h_im); nus = dHdJ(MNF.K);
    printf("%5d (%3d) nu = [%7.5f %7.5f]\n",
	   k, globval.Cell_nLoc, nus[0].cst(), nus[1].cst());
    outf << fixed << setprecision(3)
	 << setw(6) << k+1 << setw(9) << Cell[k].S
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
  outf.close();
}


int main(int argc, char *argv[])
{

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false;  globval.bpm         = 0;

  // disable from TPSALib and LieLib log messages
  idprset(-1);

  std::string home_dir = "";

  if (true)
    Read_Lattice((home_dir+argv[1]).c_str());
  else
    rdmfile(argv[1]);

  if (false) chk_bend();

  Ring_GetTwiss(true, 0e0); printglob();

  get_drv_terms();
}
