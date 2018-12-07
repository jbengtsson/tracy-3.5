#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


ss_vect<tps> get_sympl_form(const int dof)
{
  int          k;
  ss_vect<tps> Id, omega;

  Id.identity(); omega.zero();
  for (k = 0; k < dof; k++) {
    omega[2*k] = Id[2*k+1]; omega[2*k+1] = -Id[2*k];
  }
  return omega;
}


void A_At_pass(void)
{
  long int     lastpos;
  int          i;
  ss_vect<tps> omega, omega_tp, A, A_Atp;

  omega = get_sympl_form(2); omega_tp = tp_S(2, omega);

  A.identity(); putlinmat(4, globval.Ascr, A);
  A_Atp = A*tp_S(2, A);

  printf("\n   alpha_x  beta_x   alpha_y  beta_y:\n"
	 "  %8.5f %8.5f %8.5f %8.5f\n",
	 -A_Atp[x_][px_], A_Atp[x_][x_], -A_Atp[y_][py_], A_Atp[y_][y_]);

  Cell_Pass(0, globval.Cell_nLoc, A_Atp, lastpos); A_Atp = tp_S(2, A_Atp);
  Cell_Pass(0, globval.Cell_nLoc, A_Atp, lastpos);

  printf("  %8.5f %8.5f %8.5f %8.5f\n",
	 -A_Atp[x_][px_], A_Atp[x_][x_], -A_Atp[y_][py_], A_Atp[y_][y_]);
}


int main(int argc, char *argv[])
{
  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  globval.H_exact        = false;  globval.quad_fringe = false;
  globval.Cavity_on      = false;  globval.radiation   = false;
  globval.emittance      = false;  globval.IBS         = false;
  globval.pathlength     = false;  globval.bpm         = 0;
  globval.dip_edge_fudge = true;

  Ring_GetTwiss(true, 0e0); printglob();

  A_At_pass();
}
