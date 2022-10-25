#define NO 4

#include "tracy_lib.h"


// F. Klein's Erlangen Program.

int no_tps   = NO,
    ndpt_tps = 5;


tps get_H(const MNF_struct &MNF)
{
  int          i;
  tps          H, g_n;
  ss_vect<tps> Id, M_n;

  // Compute Lie generator in Floquet space.
  // K is in Dragt-Finn form but the generators commute.
  Id.identity();
  H = MNF.K;
  for (i = no_tps; i >= 3; i--) {
    g_n = Take(MNF.g, i);
    H = H*LieExp(-g_n, Id);
  }
  return H;
}


void compute_invariant()
{
  long int     lastpos;
  int          k;
  double       nu_int;
  tps          H_2, DH_2, H, DH;
  ss_vect<tps> Id, M, B;

  const double
    alpha[] = {Cell[0].Alpha[X_], Cell[0].Alpha[Y_]},
    beta[]  = {Cell[0].Beta[X_],  Cell[0].Beta[Y_]},
    gamma[] = {(1e0+sqr(alpha[X_]))/beta[X_], (1e0+sqr(alpha[Y_]))/beta[Y_]},
    eta_x   = Cell[0].Eta[X_],
    etap_x  = Cell[0].Etap[X_];

  Id.identity();

  printf("\n  alpha_x = %6.3f beta_x = %5.3f gamma_x = %5.3f"
	 " eta_x = %10.3e eta'_x = %10.3e\n"
	 "  alpha_y = %6.3f beta_y = %5.3f gamma_y = %5.3f\n",
	 alpha[X_], beta[X_], gamma[X_], eta_x, etap_x,
	 alpha[Y_], beta[Y_], gamma[Y_]);

  danot_(no_tps-1);

  M.identity();
  Cell_Pass(0, globval.Cell_nLoc, M, lastpos);
  printf("\nM:\n");
  prt_lin_map(3, M);
  // cout << scientific << setprecision(3) << "\nM:\n" << M << "\n";

  danot_(2);

  H_2 = 0e0;
  for (k = 0; k < 2; k++) {
    H_2 +=
      -M_PI*modf(globval.TotalTune[k], &nu_int)
      *(gamma[k]*sqr(Id[2*k])+2e0*alpha[k]*Id[2*k]*Id[2*k+1]
	+beta[k]*sqr(Id[2*k+1]));
  }
  H_2 += M[ct_][delta_]*sqr(Id[delta_])/2e0;

  B.identity();
  B[x_]  += eta_x*Id[delta_];
  B[px_] += etap_x*Id[delta_];
  B[ct_] += -etap_x*Id[x_] - eta_x*Id[px_];

  printf("\nLinear dispersion computed by numerical differentiation.\n");
  printf("\nB:\n");
  prt_lin_map(3, B);

  danot_(no_tps);

  printf("\nLinear dispersion computed by TPSA.\n");
  MNF = MapNorm(M, 1);

  printf("\nA0:\n");
  prt_lin_map(3, MNF.A0);

  H_2 = H_2*Inv(MNF.A0);

  daeps_(1e-10);
  H_2 = 1e0*H_2;
  cout << scientific << setprecision(3) << "\nH_2:\n" << H_2;
  daeps_(eps_tps);

  printf("\ne^-H_2:\n");
  prt_lin_map(3, LieExp(H_2, Id));

  // Restore max order after call to LieExp.
  danot_(2);

  DH_2 = H_2*M - H_2;
  daeps_(1e-13);
  DH_2 = 1e0*DH_2;
  cout << scientific << setprecision(3) << "\nH_2*M - H_2:\n" << DH_2;
  daeps_(eps_tps);

  danot_(no_tps);

  MNF = MapNorm(M, no_tps);
  H = get_H(MNF)*Inv(MNF.A0*MNF.A1);

  daeps_(1e-7);
  H = 1e0*H;
  cout << scientific << setprecision(3) << "\nH:\n" << H;
  daeps_(eps_tps);

  DH = H*M - H;
  daeps_(1e-7);
  DH = 1e0*DH;
  cout << scientific << setprecision(3) << "\nDH:\n" << DH;
  daeps_(eps_tps);
}


int main(int argc, char *argv[])
{
  globval.H_exact    = false; globval.quad_fringe    = false;
  globval.Cavity_on  = false; globval.radiation      = false;
  globval.emittance  = false; globval.IBS            = false;
  globval.pathlength = false; globval.bpm            = 0;
  globval.Cart_Bend  = false; globval.dip_edge_fudge = true;
  globval.mat_meth   = false;

  // disable from TPSALib- and LieLib log messages
  idprset_(-1);

  if (false)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  Ring_GetTwiss(true, 0e0); printglob();

  compute_invariant();
}

