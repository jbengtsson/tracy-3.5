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
   double dnu[nd_tps];
   Matrix A_tp;
   
   for (auto j = 0; j < 2*nd_tps; j++)
     for (auto k = 0; k < 2*nd_tps; k++)
       Sigma[j][k] = (j == k)? globval.eps[k/2] : 0e0;
   
   MulLMat(2*nd_tps, globval.Ascr, Sigma);
   CopyMat(2*nd_tps, globval.Ascr, A_tp);
   TpMat(2*nd_tps, A_tp);
   MulRMat(2*nd_tps, Sigma, A_tp);
}


void check_A()
{
  long int     lastpos;
  double       dnu[nd_tps];
  ss_vect<tps> A1,A0;

  A1 = A0 = get_A_CS(2, putlinmat(2*nd_tps, globval.Ascr), dnu);
  
  printf("\nA before:\n");
  prt_lin_map(nd_tps, A0);
  
  Cell_Pass(0, globval.Cell_nLoc, A1, lastpos);
  
  A1 = get_A_CS(2, A1, dnu);
    
  printf("\nA after:\n");
  prt_lin_map(nd_tps, A1);
  
  printf("\nA1 - A0:\n");
  prt_lin_map(nd_tps, A1-A0);
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
  long int     lastpos;
  Matrix       M, M_tp, Sigma0, Sigma1, DSigma;
  ss_vect<tps> map;
  
  map = compute_map();
  
  getlinmat(2*nd_tps, map, M);
  CopyMat(2*nd_tps, M, M_tp);
  TpMat(2*nd_tps, M_tp);
  
  compute_sigma(Sigma0);

  CopyMat(2*nd_tps, Sigma0, Sigma1);
  MulLMat(6, M, Sigma1);
  MulRMat(6, Sigma1, M_tp);
  
  CopyMat(2*nd_tps, Sigma1, DSigma);
  SubMat (2*nd_tps, Sigma0, DSigma);

  printf("\nSigma before:\n");
  prtmat(2*nd_tps, Sigma0);
  printf("\nSigma after:\n");
  prtmat(2*nd_tps, Sigma1);
  printf("\nSigma1 - Sigma0:\n");
  prtmat(2*nd_tps, DSigma);
 
}


void benchmark(void)
{
  double
    eps_x, sigma_delta, U_0, J[nd_tps], tau[nd_tps], I[2*nd_tps], dnu[nd_tps];
  Matrix
    A_tp_mat;
  ss_vect<tps>
    A, A_tp, Sigma;
  
  if (!globval.mat_meth) {
    A = putlinmat(2*nd_tps, globval.Ascr);
    CopyMat(2*nd_tps, globval.Ascr, A_tp_mat);
    TpMat(2*nd_tps, A_tp_mat);
    A_tp = putlinmat(2*nd_tps, A_tp_mat);

    printf("\nA_CS:\n");
    prt_lin_map(nd_tps, get_A_CS(2, A, dnu));

    printf("\nA:\n");
    prt_lin_map(nd_tps, A);

    printf("\nA^T:\n");
    prt_lin_map(nd_tps, A_tp);

    printf("\nsigma [");
    for (auto k = 0; k < 2*nd_tps; k++) {
      printf("%9.3e", sqrt(Cell[globval.Cell_nLoc].sigma[k][k]));
      if (k != 5)
	printf(" ");
    }
    printf("]\n");

    printf("\nCell[{end}].Sigma:\n");
    prtmat(2*nd_tps, Cell[globval.Cell_nLoc].sigma);

    Sigma = putlinmat(2*nd_tps, Cell[globval.Cell_nLoc].sigma);

    printf("\nA^-1*Sigma*A\n");
    prt_lin_map(nd_tps, A*Sigma*Inv(A_tp));
  } else
    get_eps_x(eps_x, sigma_delta, U_0, J, tau, I, true);

  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);
  prt_chrom_lat();
  prtmfile("flat_file.dat");
}


int main(int argc, char *argv[])
{
  double dnu[nd_tps];
  
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
