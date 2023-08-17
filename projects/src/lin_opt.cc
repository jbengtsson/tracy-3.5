#include <assert.h>

#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


void set_lat_state(const bool mat_meth)
{
  trace                  = false;
  reverse_elem           = true;
  globval.mat_meth       = mat_meth;

  globval.H_exact        = false;
  globval.quad_fringe    = false;
  globval.Cavity_on      = false;
  globval.radiation      = false;
  globval.emittance      = false;
  globval.IBS            = false;
  globval.pathlength     = false;
  globval.Aperture_on    = false;
  globval.Cart_Bend      = false;
  globval.dip_edge_fudge = true;
}


int main(int argc, char *argv[])
{
  double       eps_x, sigma_delta, U_0, J[3], tau[3], I[6], dnu[3];
  Matrix       A_mat_tp;
  ss_vect<tps> A, A_tp, Sigma;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  set_lat_state(false);

  Ring_GetTwiss(true, 0e0);

  printglob();

  if (!globval.mat_meth) {
    GetEmittance(ElemIndex("cav"), true);

    A = putlinmat(6, globval.Ascr);
    CopyMat(6, globval.Ascr, A_mat_tp);
    TpMat(6, A_mat_tp);
    A_tp = putlinmat(6, A_mat_tp);

    printf("\nA_CS:\n");
    prt_lin_map(3, get_A_CS(2, A, dnu));

    printf("\nA:\n");
    prt_lin_map(3, A);

    printf("\nA^T:\n");
    prt_lin_map(3, A_tp);

    printf("\nsigma [");
    for (auto k = 0; k < 2*nd_tps; k++) {
      printf("%9.3e", sqrt(Cell[globval.Cell_nLoc].sigma[k][k]));
      if (k != 5)
	printf(" ");

    }
    printf("]\n");

    printf("\nCell[{end}].Sigma:\n");
    prtmat(6, Cell[globval.Cell_nLoc].sigma);

    Sigma = putlinmat(6, Cell[globval.Cell_nLoc].sigma);

    printf("\nA^-1*Sigma*A\n");
    prt_lin_map(3, A*Sigma*Inv(A_tp));
  } else
    get_eps_x(eps_x, sigma_delta, U_0, J, tau, I, true);

  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);
  prt_chrom_lat();
  prtmfile("flat_file.dat");
}
