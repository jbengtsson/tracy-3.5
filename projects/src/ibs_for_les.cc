#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


void prt_ZAP(const int n)
{
  // Generate linear optics file for ZAP:
  //   M. Zisman, S. Chattopadhyay, J. Bisognano "ZAP User's Manual"
  //   LBL-21270 (1986).
  long int k;
  FILE     *outf;

  outf = file_write("ZAPLAT.DAT");

  fprintf(outf, "%ld %7.5f\n",
	  Lattice.param.Cell_nLoc+1, n*Lattice.Cell[Lattice.param.Cell_nLoc].S);
  fprintf(outf, "One super period\n");

  for (k = 0; k <= Lattice.param.Cell_nLoc; k++)
    fprintf(outf, "%10.5f %8.5f %9.6f %8.5f %7.3f %8.5f %8.5f %7.5f\n",
	    Lattice.Cell[k].S,
	    Lattice.Cell[k].Beta[X_], Lattice.Cell[k].Alpha[X_],
	    Lattice.Cell[k].Beta[Y_], Lattice.Cell[k].Alpha[Y_],
	    Lattice.Cell[k].Eta[X_], Lattice.Cell[k].Etap[X_],
	    Lattice.Cell[k].maxampl[X_][1]);

  fprintf(outf, "0\n");

  fclose(outf);
}


void get_IBS(const double Qb, const double eps[])
{
  int    j;
  double eps1[3], sigma_s, sigma_delta;

  const int n_iter = 20;

  for (j = 0; j < 3; j++)
    eps1[j] = eps[j];

  printf("\n");
  for (j = 1; j <= n_iter; j++) {
    if (j == 1) {
      printf("\nIBS %d:\n", j);
      Lattice.IBS_BM(Qb, eps, eps1, true, true);
    } else if ((j == n_iter-1) || (j == n_iter)) {
      printf("\nIBS %d:\n", j);
      Lattice.IBS_BM(Qb, eps, eps1, false, true);
    } else
      Lattice.IBS_BM(Qb, eps, eps1, false, false);
  }

  sigma_s = sqrt(Lattice.param.beta_z*eps1[Z_]);
  sigma_delta = eps1[Z_]/sigma_s;
  printf("eps_x = %9.3e, eps_y = %9.3e, sigma_s = %9.3e, sigma_delta = %9.3e\n",
	 eps1[X_], eps1[Y_], sigma_s, sigma_delta);
}


int main(int argc, char *argv[])
{
  double sigma_s, sigma_delta, gamma_z;

  Lattice.param.H_exact    = false; Lattice.param.quad_fringe = false;
  Lattice.param.Cavity_on  = false; Lattice.param.radiation   = false;
  Lattice.param.emittance  = false; Lattice.param.IBS         = false;
  Lattice.param.pathlength = false; Lattice.param.bpm         = 0;

  // Set bunch charge and vertical emittance here.
  const double Qb = 5e-9, eps_y = 0.008e-9;

  if (true)
    Lattice.Read_Lattice(argv[1]);
  else {
    Lattice.param.Energy = 3e0;
    Lattice.rdmfile(argv[1]);
  }

  Lattice.Ring_GetTwiss(true, 0.0); printglob();

  if (true) prt_ZAP(20);

  Lattice.GetEmittance(Lattice.Elem_Index("cav"), true);

  gamma_z = (1e0+sqr(Lattice.param.alpha_z))/Lattice.param.beta_z;

  // Tweak bunch length here.
  sigma_s     = sqrt(Lattice.param.beta_z*Lattice.param.eps[Z_]);
  sigma_delta = sqrt(gamma_z*Lattice.param.eps[Z_]);
  printf("\nsigma_s [m] = %9.3e, sigma_delta = %9.3e\n", sigma_s, sigma_delta);

  Lattice.param.eps[Y_] = eps_y;

  // alpha_z << 1 => eps_z ~ sigma_s * sigma_delta.
  Lattice.param.eps[Z_] = sigma_s*sigma_delta;
  Lattice.param.beta_z  = sqr(sigma_s)/Lattice.param.eps[Z_];

  printf("\nbeta_z = %5.3f, eps_z = %10.3e\n",
	 Lattice.param.beta_z, Lattice.param.eps[Z_]);

  get_IBS(Qb, Lattice.param.eps);
}
