#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


double get_eps_x1(void)
{
  bool             cav, emit;
  long int         lastpos;
  double           eps_x;
  ss_vect<tps>     A;

  cav = Lattice.param.Cavity_on; emit = Lattice.param.emittance;

  Lattice.param.Cavity_on = false; Lattice.param.emittance = false;

  Lattice.Ring_GetTwiss(false, 0e0);

  putlinmat(6, Lattice.param.Ascr, A);

  prt_lin_map(3, A);

  Lattice.param.emittance = true;

  Lattice.Cell_Pass(0, Lattice.param.Cell_nLoc, A, lastpos);

  eps_x = 1470e0*pow(Lattice.param.Energy, 2)*I5/(I2-I4);

  printf("\neps_x = %5.3f pm.rad, J_x = %5.3f, J_z = %5.3f \n",
	 1e3*eps_x, 1e0-I4/I2, 2e0+I4/I2);

  Lattice.param.Cavity_on = cav; Lattice.param.emittance = emit;

  return eps_x;
}


void prt_ZAP(const int n)
{
  long int  k;
  FILE      *outf;

  outf = file_write("ZAPLAT.DAT");

  fprintf(outf, "%ld %7.5f\n",
	  Lattice.param.Cell_nLoc+1, n*Lattice.Cell[Lattice.param.Cell_nLoc].S);
  fprintf(outf, "One super period\n");

  for (k = 0; k <= Lattice.param.Cell_nLoc; k++)
    fprintf(outf, "%10.5f %8.5f %9.6f %8.5f %7.3f %8.5f %8.5f %7.5f\n",
	    Lattice.Cell[k].S,
	    Lattice.Cell[k].Beta[X_], Lattice.Cell[k].Alpha[X_],
	    Lattice.Cell[k].Beta[Y_], Lattice.Cell[k].Alpha[Y_],
	    Lattice.Cell[k].Eta[X_], Lattice.Cell[k].Etap[X_], Lattice.Cell[k].maxampl[X_][1]);

  fprintf(outf, "0\n");

  fclose(outf);
}


void get_Touschek(void)
{
  long int      k;
  double        gamma_z, eps[3];
  double        sum_delta[Lattice.param.Cell_nLoc+1][2];
  double        sum2_delta[Lattice.param.Cell_nLoc+1][2];
  FILE          *fp;

  const double  Qb = 1.3e-9;

  Lattice.param.eps[X_] = 2.040e-9;
  Lattice.param.eps[Y_] = 8e-12;
  Lattice.param.eps[Z_] = 1.516e-6;

  Lattice.param.alpha_z = 2.502e-02; Lattice.param.beta_z = 5.733;

  if (false) {
    Lattice.param.delta_RF = 3.0e-2;
//     Touschek(Qb, Lattice.param.delta_RF, 0.85e-9, 8e-12, 0.84e-3, 4.4e-3);
    gamma_z = (1.0+sqr(Lattice.param.alpha_z))/Lattice.param.beta_z;
//     Touschek(Qb, Lattice.param.delta_RF, eps[X_], eps[Y_],
// 	     sqrt(gamma_z*eps[Z_]), sqrt(Lattice.param.beta_z*eps[Z_]));
    Lattice.Touschek(Qb, 3.03e-2, 2.446e-9, 9.595e-12, 0.580e-3, 3.33e-3);

    if (false) {
      // initialize momentum aperture arrays
      for(k = 0; k <= Lattice.param.Cell_nLoc; k++){
	sum_delta[k][0] = 0.0; sum_delta[k][1] = 0.0;
	sum2_delta[k][0] = 0.0; sum2_delta[k][1] = 0.0;
      }

      Lattice.Touschek(1.3e-9, Lattice.param.delta_RF, false,
		       0.85e-9, 0.008e-9, 0.84e-3, 4.4e-3,
		       512, true, sum_delta, sum2_delta);

      fp = file_write("mom_aper.out");
      for(k = 0; k <= Lattice.param.Cell_nLoc; k++)
	fprintf(fp, "%4ld %7.2f %5.3f %6.3f\n",
		k, Lattice.Cell[k].S, 1e2*sum_delta[k][0], 1e2*sum_delta[k][1]);
    }
  }
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
  sigma0_s = sqrt(Lattice.param.beta_z*eps0[Z_]); sigma0_delta = eps0[Z_]/sigma0_s;

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
	Lattice.IBS_BM(Qb, eps, eps1, true, true);
      } else if ((j == n_iter-1) || (j == n_iter)) {
	printf("\nIBS %d:\n", j);
	Lattice.IBS_BM(Qb, eps0, eps1, false, true);
      } else
	Lattice.IBS_BM(Qb, eps0, eps1, false, false);
    }

    sigma_s = sqrt(Lattice.param.beta_z*eps1[Z_]);
    sigma_delta = eps1[Z_]/sigma_s;
    fprintf(outf, "%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %12.3e\n",
	    1e2*sigma0_s, 1e12*eps1[X_], eps1[X_]/eps0[X_], 1e12*eps1[Y_],
	    eps1[Y_]/eps0[Y_], 1e2*sigma_s, sigma_delta);
    fflush(outf);

    // alpha_z << 1 => eps_z ~ sigma_s * sigma_delta.
    sigma0_s += ds; eps0[Z_] = sigma0_s*sigma0_delta;
    Lattice.param.beta_z = sqr(sigma0_s)/eps0[Z_];
  }

  fclose(outf);
}


int main(int argc, char *argv[])
{
  int           qf, qd, sf, sd;
  double        b2[2], a2, b3[2], a3;
  ss_vect<tps>  map;

  Lattice.param.H_exact    = false; Lattice.param.quad_fringe = false;
  Lattice.param.Cavity_on  = false; Lattice.param.radiation   = false;
  Lattice.param.emittance  = false; Lattice.param.IBS         = false;
  Lattice.param.pathlength = false; Lattice.param.bpm         = 0;

  const double delta = 3e-2;
  const double Qb    = 5e-9, sigma_s = 1e-2, sigma_delta = 1e-3;
  const double nu[]  = {101.9/20.0, 27.6/20.0};

  if (true)
    Lattice.Read_Lattice(argv[1]);
  else {
    Lattice.param.Energy = 3e0;
    Lattice.rdmfile(argv[1]);
  }

  if (false) no_sxt();

  Lattice.Ring_GetTwiss(true, 0e0); printglob();

  if (false) {
    qf = Lattice.Elem_Index("qf"); qd = Lattice.Elem_Index("bh");
    Lattice.FitTune(qf, qd, nu[X_], nu[Y_]);
    get_bn_design_elem(qf, 1, Quad, b2[0], a2);
    get_bn_design_elem(qd, 1, Quad, b2[1], a2);

    printf("\nnu_x = %7.5f nu_y = %7.5f\n",
	   Lattice.param.TotalTune[X_], Lattice.param.TotalTune[Y_]);
    printf("  qf = %8.5f   qd = %8.5f\n", b2[0], b2[1]);

    Lattice.Ring_GetTwiss(true, 0.0); printglob();
  }

  if (false) {
    sf = Lattice.Elem_Index("sfh"); sd = Lattice.Elem_Index("sd");
    Lattice.FitChrom(sf, sd, 0e0, 0e0);
    get_bn_design_elem(sf, 1, Sext, b3[0], a3);
    get_bn_design_elem(sd, 1, Sext, b3[1], a3);

    printf("\nsf = %10.5f, sd = %10.5f", b3[0], b3[1]);

    Lattice.Ring_GetTwiss(true, 0.0); printglob();
  }

  Lattice.prt_lat("linlat1.out", Lattice.param.bpm, true);
  Lattice.prt_lat("linlat.out", Lattice.param.bpm, true, 10);
  Lattice.prt_chrom_lat();

  Lattice.prtmfile("flat_file.dat");

  if (false) prt_ZAP(20);

  get_eps_x1();
  Lattice.GetEmittance(Lattice.Elem_Index("cav"), true);

  if (false) {
    printf("\nalpha_z = %11.3e, beta_z = %10.3e\n",
	   Lattice.param.alpha_z,  Lattice.param.beta_z);

    Lattice.param.eps[Y_] = Lattice.param.eps[X_];

    // alpha_z << 1 => eps_z ~ sigma_s * sigma_delta.
    Lattice.param.eps[Z_] = sigma_s*sigma_delta;
    Lattice.param.beta_z  = sqr(sigma_s)/Lattice.param.eps[Z_];

    printf("\nbeta_z = %5.3f, eps_z = %10.3e\n",
	   Lattice.param.beta_z, Lattice.param.eps[Z_]);

    get_IBS(21, 1e-2, Qb, Lattice.param.eps);
  }

//   get_Touschek();

  if (false) {
    Lattice.param.Cavity_on = true;
    Lattice.get_dynap(delta, 25, 512, true);
  }
}
