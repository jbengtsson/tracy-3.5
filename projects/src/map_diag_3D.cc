#include <assert.h>

#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


ss_vect<tps> compute_R(const int n_dof)
{
  ss_vect<tps> Id, R;

  Id.identity();
  R.identity();
  for (int k = 0; k < n_dof; k++) {
    if (k < 2) {
      R[2*k] =
	cos(2*M_PI*globval.TotalTune[k])*Id[2*k]
	+ sin(2*M_PI*globval.TotalTune[k])*Id[2*k+1];
      R[2*k+1] =
	-sin(2*M_PI*globval.TotalTune[k])*Id[2*k]
	+ cos(2*M_PI*globval.TotalTune[k])*Id[2*k+1];
    } else {
      R[ct_] =
	cos(2*M_PI*globval.TotalTune[k])*Id[ct_]
	- sin(2*M_PI*globval.TotalTune[k])*Id[delta_];
      R[delta_] =
	sin(2*M_PI*globval.TotalTune[k])*Id[ct_]
	+ cos(2*M_PI*globval.TotalTune[k])*Id[delta_];
    }
  }
  return R;
}


ss_vect<tps> compute_A0(void)
{
  ss_vect<tps> Id, A0;

  Id.identity();
  A0.identity();
  A0[x_]  += Cell[0].Eta[X_]*Id[delta_];
  A0[px_] += Cell[0].Etap[X_]*Id[delta_];
  A0[ct_] +=
    Cell[0].Etap[X_]*Id[x_] - Cell[0].Eta[X_]*Id[px_]
    + Cell[0].Eta[X_]*Cell[0].Etap[X_]*Id[delta_];
  return A0;
}


ss_vect<tps> compute_A1(const int n_dof)
{
  ss_vect<tps> Id, A1;

  Id.identity();
  A1.identity();
  for (int k = 0; k < n_dof; k++) {
    if (k < 2) {
      A1[2*k] = sqrt(Cell[0].Beta[k])*Id[2*k];
      A1[2*k+1] =
	-Cell[0].Alpha[k]/sqrt(Cell[0].Beta[k])*Id[2*k]
	+ 1e0/sqrt(Cell[0].Beta[k])*Id[2*k+1];
    } else {
      A1[ct_] =
	sqrt(globval.beta_z)*Id[ct_]
	- globval.alpha_z/sqrt(globval.beta_z)*Id[delta_];
      A1[delta_] =
	
	+ 1e0/sqrt(globval.beta_z)*Id[delta_];
    }
  }
  return A1;
}


void compute_M(const int n_dof, ss_vect<tps> A0)
{
  ss_vect<tps> R, A1, M;

  R = compute_R(n_dof);
  A1 = compute_A1(n_dof);
  M = A0*A1*R*Inv(A0*A1);

  printf("\nA0:");
  prt_lin_map(nd_tps, A0);

  printf("\nA1:");
  prt_lin_map(nd_tps, A1);

  printf("\nA_scr:");
  prt_lin_map(nd_tps, Inv(A0)*putlinmat(2*nd_tps, globval.Ascr));

  if (n_dof == 2) {
    double C = Cell[globval.Cell_nLoc].S;
    M[ct_] += globval.Alphac*C*tps(0e0, delta_+1);
  }

  printf("\nM:");
  prt_lin_map(nd_tps, M);
  // prt_lin_map(nd_tps, M-putlinmat(6, globval.OneTurnMat));
}


void set_state(void)
{
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
  ss_vect<tps> A0;

  trace            = false;
  reverse_elem     = true;
  globval.mat_meth = !false;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  set_state();

  globval.Cavity_on = globval.radiation = false;
  Ring_GetTwiss(true, 0e0);
  A0 = compute_A0();
  compute_M(2, A0);
  printglob();

  globval.Cavity_on = globval.radiation = true;
  Ring_GetTwiss(true, 0e0);
  globval.TotalTune[Z_] = globval.Omega;
  globval.alpha_z =
    -globval.Ascr[ct_][ct_]*globval.Ascr[delta_][ct_]
    - globval.Ascr[ct_][delta_]*globval.Ascr[delta_][delta_];
  globval.beta_z = sqr(globval.Ascr[ct_][ct_]) + sqr(globval.Ascr[ct_][delta_]);
  printglob();
  compute_M(3, A0);
}
