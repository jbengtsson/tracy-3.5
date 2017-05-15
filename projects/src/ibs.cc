#define NO 1

#include "tracy_lib.h"

int no_tps = NO;

extern bool freq_map;
extern int  N_Fam, Q_Fam[];

// dynamic aperture
const double delta  = 3e-2; // delta for off-momentum aperture


double get_eps_x1(void)
{
  bool         cav, emit;
  long int     lastpos;
  double       eps_x;
  ss_vect<tps> A;

  cav = globval.Cavity_on; emit = globval.emittance;

  globval.Cavity_on = false; globval.emittance = false;

  Ring_GetTwiss(false, 0e0);

  putlinmat(6, globval.Ascr, A);

  prt_lin_map(3, A);

  globval.emittance = true;

  Cell_Pass(0, globval.Cell_nLoc, A, lastpos);

  eps_x = 1470e0*pow(globval.Energy, 2)*I5/(I2-I4);

  printf("\neps_x = %5.3f pm.rad, J_x = %5.3f, J_z = %5.3f \n",
	 1e3*eps_x, 1e0-I4/I2, 2e0+I4/I2);

  globval.Cavity_on = cav; globval.emittance = emit;

  return eps_x;
}


void prt_ZAP(const int n)
{
  long int  k;
  FILE      *outf;

  outf = file_write("ZAPLAT.DAT");

  fprintf(outf, "%ld %7.5f\n",
	  globval.Cell_nLoc+1, n*Cell[globval.Cell_nLoc].S);
  fprintf(outf, "One super period\n");

  for (k = 0; k <= globval.Cell_nLoc; k++)
    fprintf(outf, "%10.5f %8.5f %9.6f %8.5f %7.3f %8.5f %8.5f %7.5f\n",
	    Cell[k].S,
	    Cell[k].Beta[X_], Cell[k].Alpha[X_],
	    Cell[k].Beta[Y_], Cell[k].Alpha[Y_],
	    Cell[k].Eta[X_], Cell[k].Etap[X_], Cell[k].maxampl[X_][1]);

  fprintf(outf, "0\n");

  fclose(outf);
}


void get_IBS(const int n, const double ds, const double Qb, const double eps[])
{
  int    i, j;
  double eps0[3], eps1[3], sigma0_s, sigma0_delta, sigma_s, sigma_delta;
  FILE   *outf;

  const char file_name[] = "ibs.out";
  const int n_iter = 20;

  outf = file_write(file_name);

  for (j = 0; j < 3; j++)
    eps0[j] = eps[j];

  // alpha_z << 1 => eps_z ~ sigma_s * sigma_delta.
  sigma0_s = sqrt(globval.beta_z*eps0[Z_]); sigma0_delta = eps0[Z_]/sigma0_s;

  fprintf(outf, "# sigma0_s eps_x    ratio_x  eps_y    ratio_y"
	  "  sigma_s  sigma_delta\n");
  fprintf(outf, "#   [cm]  [pm.rad]          [pm.rad]   [cm]\n");
  for (i = 1; i <= n; i++) {
    for (j = 0; j < 3; j++)
      eps1[j] = eps0[j];

    printf("\n");
    for (j = 1; j <= n_iter; j++) {
      if (j == 1) {
	printf("\nIBS %d:\n", j);
	IBS_BM(Qb, eps, eps1, true, true);
      } else if ((j == n_iter-1) || (j == n_iter)) {
	printf("\nIBS %d:\n", j);
	IBS_BM(Qb, eps0, eps1, false, true);
      } else
	IBS_BM(Qb, eps0, eps1, false, false);
    }

    sigma_s = sqrt(globval.beta_z*eps1[Z_]);
    sigma_delta = eps1[Z_]/sigma_s;
    fprintf(outf, "%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %12.3e\n",
	    1e2*sigma0_s, 1e12*eps1[X_], eps1[X_]/eps0[X_], 1e12*eps1[Y_],
	    eps1[Y_]/eps0[Y_], 1e2*sigma_s, sigma_delta);
    fflush(outf);

    // alpha_z << 1 => eps_z ~ sigma_s * sigma_delta.
    sigma0_s += ds; eps0[Z_] = sigma0_s*sigma0_delta;
    globval.beta_z = sqr(sigma0_s)/eps0[Z_];
  }

  fclose(outf);
}


int main(int argc, char *argv[])
{
  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  const double Qb   = 5e-9, sigma_s = 1e-2, sigma_delta = 1e-3;

  if (true)
    Read_Lattice(argv[1]);
  else {
    globval.Energy = 3e0;
    rdmfile(argv[1]);
  }

  Ring_GetTwiss(true, 0.0); printglob();

 if (false) prt_ZAP(20);

  get_eps_x1();
  GetEmittance(ElemIndex("cav"), true);

  printf("\nalpha_z = %11.3e, beta_z = %10.3e\n",
	 globval.alpha_z,  globval.beta_z);

  globval.eps[Y_] = globval.eps[X_];

  // alpha_z << 1 => eps_z ~ sigma_s * sigma_delta.
  globval.eps[Z_] = sigma_s*sigma_delta;
  globval.beta_z  = sqr(sigma_s)/globval.eps[Z_];

  printf("\nbeta_z = %5.3f, eps_z = %10.3e\n",
	 globval.beta_z, globval.eps[Z_]);

  get_IBS(21, 1e-2, Qb, globval.eps);
}
