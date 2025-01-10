/* Author: Johan Bengtsson.

   Module:

   Control of dynamic aperture.
   Pascal module developed for the SLS conceptual design in 1995.
   Machine translated to C - by utilising P2C - in 2004 for the conceptual
   design of NSLS-II.

   Description:

   Procedures to compute and minimize the sextupolar driving terms to second
   order, amplitude dependent tune shifts(, and second order chromaticity).

   References:

   [1] J. Bengssson 𝐶𝑜𝑛𝑡𝑟𝑜𝑙 𝑜𝑓 𝐷𝑦𝑛𝑎𝑚𝑖𝑐 𝐴𝑝𝑒𝑟𝑡𝑢𝑟𝑒 𝑓𝑜𝑟 𝑆𝑦𝑛𝑐ℎ𝑟𝑜𝑡𝑟𝑜𝑛 𝐿𝑖𝑔ℎ𝑡 𝑆𝑜𝑢𝑟𝑐𝑒𝑠 PAC 2005.

       https://accelconf.web.cern.ch/P05/PAPERS/MPPE020.PDF

   [2] J. Bengtsson 𝑇ℎ𝑒 𝑆𝑒𝑥𝑡𝑢𝑝𝑜𝑙𝑒 𝑆𝑐ℎ𝑒𝑚𝑒 𝑓𝑜𝑟 𝑡ℎ𝑒 𝑆𝑤𝑖𝑠𝑠 𝐿𝑖𝑔ℎ𝑡 𝑆𝑜𝑢𝑟𝑐𝑒 (𝑆𝐿𝑆):
       𝐴𝑛 𝐴𝑛𝑎𝑙𝑦𝑡𝑖𝑐 𝐴𝑝𝑝𝑟𝑜𝑎𝑐ℎ SLS Note 9/97 (1997).

       https://ados.web.psi.ch/slsnotes/sls0997.pdf

*/

#define NO 1

#include "tracy_lib.h"

int no_tps = NO;

// spatial components
//enum spatial_index { X_ = 0, Y_ = 1, Z_ = 2 };

// phase space components
//enum ps_index { x_ = 0, px_ = 1, y_ = 2, py_ = 3, delta_ = 4, ct_ = 5 };

// 2.5 degrees of freedom, i.e. delta is treated as a parameter
const int ps_dim = 4+1;

typedef double sp_vec[2], ps_vec[ps_dim], mat[ps_dim][ps_dim];

int iter, n_b3, b3s[n_b3_max];

const double
  beta_inj[] = {2.8,  2.8},
  A_max[]    = {3e-3, 1.5e-3},
  delta_max  = 2e-2,
  // twoJ[] = {sqr(A_max[X_])/beta_inj[X_], sqr(A_max[Y_])/beta_inj[Y_]};
  twoJ[] = {1e0, 1e0};


double get_bn(CellType &Cell, long int Order)
{
  double bn;

  if (Cell.Elem.Pkind == Mpole)
    bn = Cell.Elem.M->PBpar[Order+HOMmax];
  else
    bn = 0e0;

  return(bn);
}


double get_bnL(CellType &Cell, long int Order)
{
  double bnL;

  if (Cell.Elem.Pkind == Mpole)
    if (Cell.Elem.M->Pthick == thick)
      bnL = Cell.Elem.M->PBpar[Order+HOMmax]*Cell.Elem.PL;
    else
      bnL = Cell.Elem.M->PBpar[Order+HOMmax];
  else
    bnL = 0e0;

  return(bnL);
}


long int ind(long int k)
{
  long int n = 0;

  if ((0 <= k) && (k <= globval.Cell_nLoc))
    n = k;
  else if (k == -1)
    n = globval.Cell_nLoc;
  else
    printf("ind: incorrect index %ld\n", k);

  return n;
}


void K(const double nu_x, const double nu_y, double a[])
{
  /* Amplitude dependent tune shifts. */

  long int k, n1, n2;
  double   b3L1, b3L2, pi_1x, pi_3x, pi_1xm2y, pi_1xp2y;
  double   s_1x, s_3x, s_1xm2y, s_1xp2y, c_1x, c_3x, c_1xm2y, c_1xp2y;
  double   dmu_x, dmu_y, A;
  Vector2  alpha0, beta0, nu0, eta0, etap0;
  Vector2  beta1, nu1, beta2, nu2;

  pi_1x = M_PI*nu_x;
  pi_3x = 3e0*M_PI*nu_x;
  pi_1xm2y = M_PI*(nu_x-2e0*nu_y);
  pi_1xp2y = M_PI*(nu_x+2e0*nu_y);

  s_1x = sin(pi_1x);
  s_3x = sin(pi_3x);
  s_1xm2y = sin(pi_1xm2y);
  s_1xp2y = sin(pi_1xp2y);

  for (k = 0; k <= 2; k++)
    a[k] = 0e0;

  for (n1 = 0; n1 <= globval.Cell_nLoc; n1++) {
    if ((Cell[n1].Elem.Pkind == Mpole) && (Cell[n1].Elem.M->Porder >= Sext)) {
      for (k = 0; k <= 1; k++) {
	alpha0[k] = Cell[ind(n1-1)].Alpha[k];
	beta0[k] = Cell[ind(n1-1)].Beta[k];
	nu0[k] = Cell[ind(n1-1)].Nu[k];
	eta0[k] = Cell[ind(n1-1)].Eta[k];
	etap0[k] = Cell[ind(n1-1)].Etap[k];
      }
      for (k = 0; k <= 1; k++) {
	beta1[k] = beta0[k];
	nu1[k] = nu0[k];
      }
      b3L1 = get_bnL(Cell[n1], Sext);

      for (n2 = 0; n2 <= globval.Cell_nLoc; n2++) {
	if ((Cell[n2].Elem.Pkind == Mpole) &&
	    (Cell[n2].Elem.M->Porder >= Sext)) {
	  for (k = 0; k <= 1; k++) {
	    alpha0[k] = Cell[n2-1].Alpha[k];
	    beta0[k] = Cell[n2-1].Beta[k];
	    nu0[k] = Cell[n2-1].Nu[k];
	    eta0[k] = Cell[n2-1].Eta[k];
	    etap0[k] = Cell[n2-1].Etap[k];
	  }
	  for (k = 0; k <= 1; k++) {
	    beta2[k] = beta0[k];
	    nu2[k] = nu0[k];
	  }
	  b3L2 = get_bnL(Cell[n2], Sext);
	  dmu_x = 2e0*M_PI*fabs(nu1[X_]-nu2[X_]);
	  dmu_y = 2e0*M_PI*fabs(nu1[Y_]-nu2[Y_]);
	  A = b3L1*b3L2*sqrt(beta1[X_]*beta2[X_]);
	  c_1x = cos(dmu_x-pi_1x)/s_1x;
	  c_3x = cos(3e0*dmu_x-pi_3x)/s_3x;
	  c_1xm2y = cos(dmu_x-2e0*dmu_y-pi_1xm2y)/s_1xm2y;
	  c_1xp2y = cos(dmu_x+2e0*dmu_y-pi_1xp2y)/s_1xp2y;
	  a[0] += A*beta1[X_]*beta2[X_]*(3*c_1x+c_3x)/2e0;
	  a[1] -=
	    A*beta1[Y_]*(4e0*twoJ[X_]*beta2[X_]*c_1x+2e0*twoJ[Y_]*beta2[Y_]
			 *(c_1xm2y-c_1xp2y));
	  a[2] += A*beta1[Y_]*beta2[Y_]*(4e0*c_1x+c_1xm2y+c_1xp2y)/2e0;
	}
      }
    }
  }

  for (k = 0; k < 3; k++)
    a[k] /= 32e0*M_PI;
}


void dnu_dJ(const double twoJ[])
{
  double a[3];

  K(globval.TotalTune[X_], globval.TotalTune[Y_], a);
  printf("\n");
  printf("Amplitude dependent tune shifts\n");
  printf("  a_xx = %23.16e\n", a[0]);
  printf("  a_xy = %23.16e\n", a[1]);
  printf("  a_yy = %23.16e\n", a[2]);
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
  const int n_aper = 25;

  trace            = false;
  reverse_elem     = !false;
  globval.mat_meth = false;

  Read_Lattice(argv[1]);

  set_state();

  Ring_GetTwiss(true, 0e0);
  printglob();

  prtmfile("flat_file.dat");
  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true);

  dnu_dJ(twoJ);
}
