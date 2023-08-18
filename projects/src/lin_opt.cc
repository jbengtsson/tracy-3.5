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
  globval.Cavity_on      = true;
  globval.radiation      = true;
  globval.emittance      = false;
  globval.IBS            = false;
  globval.pathlength     = false;
  globval.Aperture_on    = false;
  globval.Cart_Bend      = false;
  globval.dip_edge_fudge = true;
}


void compute_sigma(Matrix &Sigma)
{
   Matrix A_tp;
   
   for (auto j = 0; j < 6; j++)
     for (auto k = 0; k < 6; k++)
       Sigma[j][k] = (j == k)? globval.eps[k/2] : 0e0;
   
   MulLMat(6, globval.Ascr, Sigma);
   CopyMat(6, globval.Ascr, A_tp);
   TpMat(6, A_tp);
   MulRMat(6, Sigma, A_tp);
}


void check_A()
{
  long int     lastpos;
  double       dnu[3];
  ss_vect<tps> A1, A0;

  A1 = A0 = get_A_CS(2, putlinmat(6, globval.Ascr), dnu);
  
  printf("\nA before:\n");
  prt_lin_map(3, A0);
  
  Cell_Pass(0, globval.Cell_nLoc, A1, lastpos);
  
  A1 = get_A_CS(2, A1, dnu);
    
  printf("\nA after:\n");
  prt_lin_map(3, A1);
  
  printf("\nA1 - A0:\n");
  prt_lin_map(3, A1-A0);
}


ss_vect<tps> compute_map()
{

  long int     lastpos;
  ss_vect<tps> M;

  M.identity(); 
  getcod(0e0, lastpos);
  M += globval.CODvect;
  Cell_Pass(0, globval.Cell_nLoc, M, lastpos);
  M -= globval.CODvect;
  
  return M;
  
}


void track_sigma(void)
{
  Matrix       M, M_tp, Sigma0, Sigma1, DSigma;
  ss_vect<tps> map;
  
  map = compute_map();
  
  getlinmat(6, map, M);
  CopyMat(6, M, M_tp);
  TpMat(6, M_tp);
  
  compute_sigma(Sigma0);

  CopyMat(6, Sigma0, Sigma1);
  MulLMat(6, M, Sigma1);
  MulRMat(6, Sigma1, M_tp);
  
  CopyMat(6, Sigma1, DSigma);
  SubMat (6, Sigma0, DSigma);

  printf("\nSigma before:\n");
  prtmat(6, Sigma0);
  printf("\nSigma after:\n");
  prtmat(6, Sigma1);
  printf("\nSigma1 - Sigma0:\n");
  prtmat(6, DSigma);
 
}


void benchmark(void)
{
  double
    eps_x, sigma_delta, U_0, J[3], tau[3], I[6], dnu[3];
  Matrix
    A_tp_mat;
  ss_vect<tps>
    A, A_tp, Sigma;
  
  if (!globval.mat_meth) {
    A = putlinmat(6, globval.Ascr);
    CopyMat(6, globval.Ascr, A_tp_mat);
    TpMat(6, A_tp_mat);
    A_tp = putlinmat(6, A_tp_mat);

    printf("\nA_CS:\n");
    prt_lin_map(3, get_A_CS(2, A, dnu));

    printf("\nA:\n");
    prt_lin_map(3, A);

    printf("\nA^T:\n");
    prt_lin_map(3, A_tp);

    printf("\nsigma [");
    for (auto k = 0; k < 6; k++) {
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


int main(int argc, char *argv[])
{
  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  set_lat_state(false);

  Ring_GetTwiss(true, 0e0);

  printglob();
  
  GetEmittance(ElemIndex("cav"), true);
  
  if (false)
    check_A();
  
  if (false)
    benchmark();
    
  if (!false)
    track_sigma();
}
