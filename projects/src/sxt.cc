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


void cpy_mat(const mat A, mat B)
{
  int i, j;

  for (i = 0; i < ps_dim; i++)
    for (j = 0; j < ps_dim; j++)
      B[i][j] = A[i][j];
}


void get_Ascr(const sp_vec alpha, const sp_vec beta,
	      const sp_vec eta, const sp_vec etap, mat A)
{
  /* Transformation to Floquet space. */

  int i, j;

  for (i = 0; i < ps_dim; i++)
    for (j = 0; j < ps_dim; j++)
      A[i][j] = 0e0;

  A[delta_][delta_] = 1e0;

  for (i = 0; i <= 1; i++) {
    A[2*i][2*i] = sqrt(beta[i]);
    A[2*i][2*i+1] = 0e0;
    A[2*i+1][2*i] = -alpha[i]/sqrt(beta[i]);
    A[2*i+1][2*i+1] = 1e0/sqrt(beta[i]);
    A[2*i][delta_] = eta[i];
    A[2*i+1][delta_] = etap[i];
  }
}


void get_Twiss(const mat A0, const mat A1,
	       sp_vec alpha, sp_vec beta, sp_vec nu, sp_vec eta, sp_vec etap)
{
  /* Compute Twiss parameters. */

  int k;

  for (k = 0; k <= 1; k++) {
    alpha[k] = -A1[2*k][2*k]*A1[2*k+1][2*k] - A1[2*k][2*k+1]*A1[2*k+1][2*k+1];
    beta[k] = sqr(A1[2*k][2*k]) + sqr(A1[2*k][2*k+1]);
    nu[k] += (atan2(A1[2*k][2*k+1], A1[2*k][2*k])
             -atan2(A0[2*k][2*k+1], A0[2*k][2*k]))/(2e0*M_PI);
    eta[k] = A1[2*k][delta_];
    etap[k] = A1[2*k+1][delta_];
  }
}


void propagate_drift(const double L, mat A)
{
  /* Drift propagator. */

  int i, j, k;
  mat M, A1;

  for (j = 0; j < ps_dim; j++)
    for (k = 0; k < ps_dim; k++)
      if (j == k)
	M[j][k] = 1e0;
      else
	M[j][k] = 0e0;

  M[x_][px_] = L;
  M[y_][py_] = L;

  for (i = 0; i < ps_dim; i++)
    for (j = 0; j < ps_dim; j++) {
      A1[i][j] = 0e0;
      for (k = 0; k < ps_dim; k++)
	A1[i][j] += M[i][k]*A[k][j];
    }

  cpy_mat(A1, A);
}


void propagate_thin_kick(const double L, const double rho_inv,
			 const double b2, mat A)
{
  /* Thin-kick propagator. */

  int i, j, k;
  mat M, A1;

  for (j = 0; j < ps_dim; j++)
    for (k = 0; k < ps_dim; k++)
      if (j == k)
	M[j][k] = 1e0;
      else
	M[j][k] = 0e0;

  M[px_][x_] = -(b2+sqr(rho_inv))*L;
  M[py_][y_] = b2*L;
  M[px_][delta_] = rho_inv*L;

  for (i = 0; i < ps_dim; i++)
    for (j = 0; j < ps_dim; j++) {
      A1[i][j] = 0e0;
      for (k = 0; k < ps_dim; k++)
	A1[i][j] += M[i][k]*A[k][j];
    }

  cpy_mat(A1, A);
}


void propagate_fringe_field(const double L,
			    const double rho_inv, const double phi, mat A)
{
  /* Dipole harde-edge fringe field propagator. */

  int k;

  if (phi != 0e0)
    for (k = 0; k < ps_dim; k++) {
      A[px_][k] += rho_inv*tan(phi*M_PI/180e0)*A[x_][k];
      A[py_][k] -= rho_inv*tan(phi*M_PI/180e0)*A[y_][k];
    }
}


void get_quad(const double L,
	      const double rho_inv, const double phi1, const double phi2,
	      const double b2,
	      const sp_vec alpha0, const sp_vec beta0, const sp_vec nu0,
	      const sp_vec eta0, const sp_vec etap0,
	      const int m_x, const int m_y, const int n_x, const int n_y,
	      const int m, double &h_c, double &h_s)
{
  /* Quad propagator: second order symplectic integrator.
     Use Simpson's rule to integrate for thick magnets. */

  const bool prt    = false;
  const int  n_step = 25;

  int    j, k;
  double h, r, phi;
  sp_vec alpha, beta, nu, eta, etap;
  mat    A0, A1;

  for (k = 0; k <= 1; k++) {
    alpha[k] = alpha0[k];
    beta[k] = beta0[k];
    nu[k] = nu0[k];
    eta[k] = eta0[k];
    etap[k] = etap0[k];
  }

  if (prt) {
    printf("\n");
    printf("%6.3f %6.3f %5.3f %6.3f %6.3f %6.3f %6.3f %5.3f\n",
	   alpha[X_], beta[X_], nu[X_], eta[X_], etap[X_],
	   alpha[Y_], beta[Y_], nu[Y_]);
  }

  get_Ascr(alpha, beta, eta, etap, A0);
  cpy_mat(A0, A1);

  propagate_fringe_field(L, rho_inv, phi1, A1);

  h = L/n_step;
  h_c = 0e0;
  h_s = 0e0;
  for (j = 1; j <= n_step; j++) {
    get_Twiss(A0, A1, alpha, beta, nu, eta, etap);
    cpy_mat(A1, A0);

    phi = 2e0*M_PI*(n_x*nu[X_]+n_y*nu[Y_]);
    r = pow(beta[X_], m_x/2e0)*pow(beta[Y_], m_y/2e0)*pow(eta[X_], m-1);
    h_c += r*cos(phi);
    h_s += -r*sin(phi);

    propagate_drift(h/2e0, A1);
    
    if ((rho_inv != 0e0) || (b2 != 0e0)) {
      if (prt) printf("1/rho^2 = %5.3f, b2 = %5.3f\n", sqr(rho_inv), b2);
      propagate_thin_kick(h, rho_inv, b2, A1);
    }

    get_Twiss(A0, A1, alpha, beta, nu, eta, etap);
    cpy_mat(A1, A0);
    
    if (prt) {
      printf("%6.3f %6.3f %5.3f %6.3f %6.3f %6.3f %6.3f %5.3f\n",
	     alpha[X_], beta[X_], nu[X_], eta[X_], etap[X_],
	     alpha[Y_], beta[Y_], nu[Y_]);
      }

    phi = 2e0*M_PI*(n_x*nu[X_]+n_y*nu[Y_]);
    r = pow(beta[X_], m_x/2e0)*pow(beta[Y_], m_y/2e0)*pow(eta[X_], m-1);
    h_c += 4e0*r*cos(phi);
    h_s += -4e0*r*sin(phi);

    propagate_drift(h/2e0, A1);

    get_Twiss(A0, A1, alpha, beta, nu, eta, etap);
    cpy_mat(A1, A0);
    
    if (prt) {
      printf("%6.3f %6.3f %5.3f %6.3f %6.3f %6.3f %6.3f %5.3f\n",
	     alpha[X_], beta[X_], nu[X_], eta[X_], etap[X_],
	     alpha[Y_], beta[Y_], nu[Y_]);
      }

    phi = 2e0*M_PI*(n_x*nu[X_]+n_y*nu[Y_]);
    r = pow(beta[X_], m_x/2e0)*pow(beta[Y_], m_y/2e0)*pow(eta[X_], m-1);
    h_c += r*cos(phi);
    h_s += -r*sin(phi);

  }

  h_c *= b2*L/(6e0*n_step);
  h_s *= b2*L/(6e0*n_step);

  propagate_fringe_field(L, rho_inv, phi2, A1);

  get_Twiss(A0, A1, alpha, beta, nu, eta, etap);
  if (prt) {
    printf("%6.3f %6.3f %5.3f %6.3f %6.3f %6.3f %6.3f %5.3f\n",
	   alpha[X_], beta[X_], nu[X_], eta[X_], etap[X_],
	   alpha[Y_], beta[Y_], nu[Y_]);
  }
}


void get_quad1(const double L,
	      const double rho_inv, const double phi1, const double phi2, 
	      const double b2,
	      sp_vec alpha3[], sp_vec beta3[], sp_vec nu3[],
	      sp_vec eta3[], sp_vec etap3[],
	      const int m_x, const int m_y, const int n_x, const int n_y,
	      const int m, double &h_c, double &h_s)
{
  /* Use Simpson's Rule for quadrupoles */

  int    k;
  double A, phi;

//  get_Twiss3(L, rho_inv, phi1, phi2, b2, alpha3, beta3, nu3, eta3, etap3);
    
  h_c = 0e0;
  h_s = 0e0;
  for (k = 0; k <= 2; k++) {
    phi = 2e0*M_PI*(n_x*nu3[k][X_]+n_y*nu3[k][Y_]);
    A = pow(beta3[k][X_], m_x/2e0)*pow(beta3[k][Y_], m_y/2e0)
        *pow(eta3[k][X_], m-1);
    if (k == 1) A *= 4e0;
    h_c += A*cos(phi);
    h_s += -A*sin(phi);
  }
  if (m_x != 0) {
    h_c *= (b2+sqr(rho_inv))*L/6e0;
    h_s *= (b2+sqr(rho_inv))*L/6e0;
  } else {
    h_c *= b2*L/6e0;
    h_s *= b2*L/6e0;
  }
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


void sxt_1(const double scl,
	   const int i, const int j, const int k, const int l, const int m,
	   double &h_c, double &h_s)
{
  /* First order generators. */

  int     n, k1, m_x, m_y, n_x, n_y;
  double  c, s, b3L, A, phi;
  sp_vec  alpha0, beta0, nu0, eta0, etap0, beta, nu, eta;

  h_c = 0e0;
  h_s = 0e0;
  m_x = i + j;
  m_y = k + l;
  n_x = i - j;
  n_y = k - l;
  for (n = 0; n <= globval.Cell_nLoc; n++) {
    if (Cell[n].Elem.Pkind == Mpole) {
      if (((Cell[n].Elem.M->Pirho != 0e0) ||
	   (Cell[n].Elem.M->Porder == Quad)) && (m >= 1)) {
	for (k1 = 0; k1 <= 1; k1++) {
	  alpha0[k1] = Cell[ind(n-1)].Alpha[k1];
	  beta0[k1] = Cell[ind(n-1)].Beta[k1];
	  nu0[k1] = Cell[ind(n-1)].Nu[k1];
	  eta0[k1] = Cell[ind(n-1)].Eta[k1];
	  etap0[k1] = Cell[ind(n-1)].Etap[k1];
	}
	get_quad(Cell[n].Elem.PL, Cell[n].Elem.M->Pirho,
		 Cell[n].Elem.M->PTx1, Cell[n].Elem.M->PTx2,
		 get_bn(Cell[n], Quad),
		 alpha0, beta0, nu0, eta0, etap0,
		 m_x, m_y, n_x, n_y, m, c, s);
	h_c += scl*c;
	h_s += scl*s;
      } else if (Cell[n].Elem.M->Porder >= Sext) {
	for (k1 = 0; k1 <= 1; k1++) {
	  alpha0[k1] = Cell[ind(n-1)].Alpha[k1];
	  beta0[k1] = Cell[ind(n-1)].Beta[k1];
	  nu0[k1] = Cell[ind(n-1)].Nu[k1];
	  eta0[k1] = Cell[ind(n-1)].Eta[k1];
	  etap0[k1] = Cell[ind(n-1)].Etap[k1];
	}
	for (k1 = 0; k1 <= 1; k1++) {
	  beta[k1] = beta0[k1];
	  nu[k1] = nu0[k1];
	  eta[k1] = eta0[k1];
	}
	b3L = get_bnL(Cell[n], Sext);
	if (b3L != 0e0) {
	  A = scl*b3L*pow(beta[X_], m_x/2e0)*pow(beta[Y_], m_y/2e0);
	  if (m >= 1) {
	    A *= -pow(eta[X_], m);
	    if (m == 1) A *= 2e0;
	  }
	  phi = 2e0*M_PI*(n_x*nu[X_]+n_y*nu[Y_]);
	  h_c += A*cos(phi);
	  h_s -= A*sin(phi);
	}
      }
    }
  }
}


void sxt_2(const double scl,
	   const int i1, const int j1, const int k1, const int l1,
	   const int i2, const int j2, const int k2, const int l2,
	   double &h_c, double &h_s)
{
  /* Second order generators. */

  long int n1, n2, k;
  double   b3L1, b3L2, dnu, A1, A, phi;
  Vector2  alpha0, beta0, nu0, eta0, etap0;
  Vector2  beta1, nu1, beta2, nu2;

  h_c = 0e0;
  h_s = 0e0;
  for (n1 = 0; n1 <= globval.Cell_nLoc; n1++) {
    if ((Cell[n1].Elem.Pkind == Mpole) &&
	(Cell[n1].Elem.M->Porder >= Sext)) {
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
      dnu = (i1-j1)*nu1[X_] + (k1-l1)*nu1[Y_];
      A1 =
	b3L1*pow(beta1[X_], (i1+j1)/2e0)*pow(twoJ[Y_]*beta1[Y_], (k1+l1)/2e0);
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
	  phi = 2e0*M_PI*(dnu+(i2-j2)*nu2[X_]+(k2-l2)*nu2[Y_]);
	  A = scl*A1*b3L2*pow(beta2[X_], (i2+j2)/2e0)
		 *pow(beta2[Y_], (k2+l2)/2e0);
	  // [..., ...] -> i*...
	  if (n2 < n1) {
	    h_s += A*cos(phi);
	    h_c += A*sin(phi);
	  } else if (n2 > n1) {
	    h_s -= A*cos(phi);
	    h_c -= A*sin(phi);
	  }
	}
      }
    }
  }
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


void sext_terms(const double twoJ[])
{
  double h_c, h_s, c, s, a[3];

  printf("\nFirst order chromatic terms:\n");
  sxt_1(1e0/2e0, 1, 1, 0, 0, 1, h_c, h_s);
  printf("  h_11001: %23.16e %23.16e\n", h_c, h_s);
  sxt_1(-1e0/2e0, 0, 0, 1, 1, 1, h_c, h_s);
  printf("  h_00111: %23.16e %23.16e\n", h_c, h_s);

  printf("\n");
  sxt_1(1e0/4e0, 2, 0, 0, 0, 1, h_c, h_s);
  h_c *= twoJ[X_];
  h_s *= twoJ[X_];
  printf("  h_20001: %23.16e %23.16e\n", h_c, h_s);
  sxt_1(-1e0/4e0, 0, 0, 2, 0, 1, h_c, h_s);
  h_c *= twoJ[Y_];
  h_s *= twoJ[Y_];
  printf("  h_00201: %23.16e %23.16e\n", h_c, h_s);
  sxt_1(1e0, 1, 0, 0, 0, 2, h_c, h_s);
  h_c *= sqrt(twoJ[X_]);
  h_s *= sqrt(twoJ[X_]);
  printf("  h_10002: %23.16e %23.16e\n", h_c, h_s);

  printf("\nFirst order geometric terms:\n");
  sxt_1(-1e0/4e0, 2, 1, 0, 0, 0, h_c, h_s);
  h_c *= pow(twoJ[X_], 3e0/2e0);
  h_s *= pow(twoJ[X_], 3e0/2e0);
  printf("  h_21000: %23.16e %23.16e\n", h_c, h_s);
  sxt_1(-1e0/12e0, 3, 0, 0, 0, 0, h_c, h_s);
  h_c *= pow(twoJ[X_], 3e0/2e0);
  h_s *= pow(twoJ[X_], 3e0/2e0);
  printf("  h_30000: %23.16e %23.16e\n", h_c, h_s);
  sxt_1(1e0/2e0, 1, 0, 1, 1, 0, h_c, h_s);
  h_c *= sqrt(twoJ[X_])*twoJ[Y_];
  h_s *= sqrt(twoJ[X_])*twoJ[Y_];
  printf("  h_10110: %23.16e %23.16e\n", h_c, h_s);
  sxt_1(1e0/4e0, 1, 0, 2, 0, 0, h_c, h_s);
  h_c *= sqrt(twoJ[X_])*twoJ[Y_];
  h_s *= sqrt(twoJ[X_])*twoJ[Y_];
  printf("  h_10200: %23.16e %23.16e\n", h_c, h_s);
  sxt_1(1e0/4e0, 1, 0, 0, 2, 0, h_c, h_s);
  h_c *= sqrt(twoJ[X_])*twoJ[Y_];
  h_s *= sqrt(twoJ[X_])*twoJ[Y_];
  printf("  h_10020: %23.16e %23.16e\n", h_c, h_s);

  printf("\nSecond order geometric terms:\n");
  sxt_2(-1e0/32e0, 2, 1, 0, 0, 3, 0, 0, 0, h_c, h_s);
  h_c *= sqr(twoJ[X_]);
  h_s *= sqr(twoJ[X_]);
  printf("  h_40000: %23.16e %23.16e\n", h_c, h_s);
  
  sxt_2(-1e0/16e0, 1, 2, 0, 0, 3, 0, 0, 0, h_c, h_s);
  h_c *= sqr(twoJ[X_]);
  h_s *= sqr(twoJ[X_]);
  printf("  h_31000: %23.16e %23.16e\n", h_c, h_s);

  sxt_2(-1e0/16e0, 3, 0, 0, 0, 0, 1, 1, 1, c, s);
  h_c = c;
  h_s = s;
  sxt_2(-1e0/16e0, 1, 0, 1, 1, 2, 1, 0, 0, c, s);
  h_c += c;
  h_s += s;
  sxt_2(-1e0/8e0, 1, 0, 0, 2, 1, 0, 2, 0, c, s);
  h_c += c;
  h_s += s;
  h_c *= twoJ[X_]*twoJ[Y_];
  h_s *= twoJ[X_]*twoJ[Y_];
  printf("  h_20110: %23.16e %23.16e\n", h_c, h_s);
  
  sxt_2(-1e0/32e0, 1, 0, 0, 2, 2, 1, 0, 0, c, s);
  h_c = c;
  h_s = s;
  sxt_2(-1e0/32e0, 3, 0, 0, 0, 0, 1, 0, 2, c, s);
  h_c += c;
  h_s += s;
  sxt_2(-4e0/32e0, 1, 0, 0, 2, 1, 0, 1, 1, c, s);
  h_c += c;
  h_s += s;
  h_c *= twoJ[X_]*twoJ[Y_];
  h_s *= twoJ[X_]*twoJ[Y_];
  printf("  h_20020: %23.16e %23.16e\n", h_c, h_s);

  sxt_2(-1e0/32e0, 3, 0, 0, 0, 0, 1, 2, 0, c, s);
  h_c = c;
  h_s = s;
  sxt_2(-1e0/32e0, 1, 0, 2, 0, 2, 1, 0, 0, c, s);
  h_c += c;
  h_s += s;
  sxt_2(-4e0/32e0, 1, 0, 1, 1, 1, 0, 2, 0, c, s);
  h_c += c;
  h_s += s;
  h_c *= twoJ[X_]*twoJ[Y_];
  h_s *= twoJ[X_]*twoJ[Y_];
  printf("  h_20200: %23.16e %23.16e\n", h_c, h_s);
  
  sxt_2(-1e0/16e0, 1, 0, 2, 0, 1, 2, 0, 0, c, s);
  h_c = c;
  h_s = s;
  sxt_2(-1e0/16e0, 2, 1, 0, 0, 0, 1, 2, 0, c, s);
  h_c += c;
  h_s += s;
  sxt_2(-2e0/16e0, 0, 1, 1, 1, 1, 0, 2, 0, c, s);
  h_c += c;
  h_s += s;
  sxt_2(-2e0/16e0, 1, 0, 1, 1, 0, 1, 2, 0, c, s);
  h_c += c;
  h_s += s;
  h_c *= twoJ[X_]*twoJ[Y_];
  h_s *= twoJ[X_]*twoJ[Y_];
  printf("  h_11200: %23.16e %23.16e\n", h_c, h_s);
    
  sxt_2(-1e0/16e0, 0, 1, 2, 0, 1, 0, 1, 1, c, s);
  h_c = c;
  h_s = s;
  sxt_2(-1e0/16e0, 0, 1, 1, 1, 1, 0, 2, 0, c, s);
  h_c += c;
  h_s += s;
  h_c *= sqr(twoJ[Y_]);
  h_s *= sqr(twoJ[Y_]);
  printf("  h_00310: %23.16e %23.16e\n", h_c, h_s);
  
  sxt_2(-1e0/32e0, 0, 1, 2, 0, 1, 0, 2, 0, h_c, h_s);
  h_c *= sqr(twoJ[Y_]);
  h_s *= sqr(twoJ[Y_]);
  printf("  h_00400: %23.16e %23.16e\n", h_c, h_s);

  sxt_2(-1e0/32e0, 0, 3, 0, 0, 3, 0, 0, 0, c, s);
  h_c = c;
  h_s = s;
  sxt_2(-3e0/32e0, 1, 2, 0, 0, 2, 1, 0, 0, c, s);
  h_c += c;
  h_s += s;
  h_c *= sqr(twoJ[X_]);
  h_s *= sqr(twoJ[X_]);
  printf("\n  h_22000: %23.16e %23.16e\n", h_c, h_s);
  
  sxt_2(-1e0/8e0, 1, 0, 1, 1, 1, 2, 0, 0, c, s);
  h_c = c;
  h_s = s;
  sxt_2(-1e0/8e0, 2, 1, 0, 0, 0, 1, 1, 1, c, s);
  h_c += c;
  h_s += s;
  sxt_2(-1e0/8e0, 0, 1, 0, 2, 1, 0, 2, 0, c, s);
  h_c += c;
  h_s += s;
  sxt_2(-1e0/8e0, 1, 0, 0, 2, 0, 1, 2, 0, c, s);
  h_c += c;
  h_s += s;
  h_c *= twoJ[X_]*twoJ[Y_];
  h_s *= twoJ[X_]*twoJ[Y_];
  printf("  h_11110: %23.16e %23.16e\n", h_c, h_s);
  
  sxt_2(-1e0/8e0, 0, 1, 1, 1, 1, 0, 1, 1, c, s);
  h_c = c;
  h_s = s;
  sxt_2(-1e0/32e0, 0, 1, 0, 2, 1, 0, 2, 0, c, s);
  h_c += c;
  h_s += s;
  sxt_2(-1e0/32e0, 0, 1, 2, 0, 1, 0, 0, 2, c, s);
  h_c += c;
  h_s += s;
  h_c *= sqr(twoJ[Y_]);
  h_s *= sqr(twoJ[Y_]);
  printf("  h_00220: %23.16e %23.16e\n", h_c, h_s);

  K(globval.TotalTune[X_], globval.TotalTune[Y_], a);
  printf("\n");
  printf("Amplitude dependent tune shifts\n");
  printf("  a_xx = %23.16e\n", a[0]);
  printf("  a_xy = %23.16e\n", a[1]);
  printf("  a_yy = %23.16e\n", a[2]);
}


double f_sxt(double p[])
{
  const int    n_prt = 100;
  const double scl_ksi1 = 1e3, scl_geom = 1e0, scl_dnu = 1e0, b3L_max = 3.5;

  int    i;
  double f, h_c, h_s, c, s;

  iter++;

  for (i = 1; i <= n_b3; i++) {
    Setbn(b3s[i-1], Sext, p[i]);
    if (fabs(p[i]) > b3L_max) {
      if (iter % n_prt == 0) {
	printf("\n");
	printf("%4d b3L_max exceeded: %9.6f\n", iter, p[i]);
      }
      return 1e30;
    }
  }

  // h_11001
  sxt_1(1e0/4e0, 1, 1, 0, 0, 1, h_c, h_s);
  f = scl_ksi1*sqr(h_c);
  // h_00111
  sxt_1(-1e0/4e0, 0, 0, 1, 1, 1, h_c, h_s);
  f += scl_ksi1*sqr(h_c);

  // h_20001
  sxt_1(1e0/4e0, 2, 0, 0, 0, 1, h_c, h_s);
  // h_00201
  sxt_1(-1e0/4e0, 0, 0, 2, 0, 1, h_c, h_s);
  // h_10002
  sxt_1(1e0, 1, 0, 0, 0, 2, h_c, h_s);

  // h_21000
  sxt_1(-1e0/4e0, 2, 1, 0, 0, 0, h_c, h_s);
  f += scl_geom*pow(twoJ[X_], 3e0/2e0)*sqr(h_c);
  // h_30000
  sxt_1(-1e0/12e0, 3, 0, 0, 0, 0, h_c, h_s);
  f += scl_geom*pow(twoJ[X_], 3e0/2e0)*sqr(h_c);
  // h_10110
  sxt_1(1e0/2e0, 1, 0, 1, 1, 0, h_c, h_s);
  f += scl_geom*sqrt(twoJ[X_])*twoJ[Y_]*sqr(h_c);
  // h_10020
  sxt_1(1e0/4e0, 1, 0, 0, 2, 0, h_c, h_s);
  f += scl_geom*sqrt(twoJ[X_])*twoJ[Y_]*sqr(h_c);
  // h_10200
  sxt_1(1e0/4e0, 1, 0, 2, 0, 0, h_c, h_s);
  f += scl_geom*sqrt(twoJ[X_])*twoJ[Y_]*sqr(h_c);

  // h_40000
  sxt_2(-1e0/32e0, 2, 1, 0, 0, 3, 0, 0, 0, h_c, h_s);
  f += scl_geom*sqr(twoJ[X_])*sqr(h_c);
  
  // h_31000
  sxt_2(-1e0/16e0, 1, 2, 0, 0, 3, 0, 0, 0, h_c, h_s);
  f += scl_geom*sqr(twoJ[X_])*sqr(h_c);

  // h_20110
  sxt_2(-1e0/16e0, 3, 0, 0, 0, 0, 1, 1, 1, c, s);
  h_c = c;
  h_s = s;
  sxt_2(-1e0/16e0, 1, 0, 1, 1, 2, 1, 0, 0, c, s);
  h_c += c;
  h_s += s;
  sxt_2(-1e0/8e0, 1, 0, 0, 2, 1, 0, 2, 0, c, s);
  h_c += c;
  h_s += s;
  f += scl_geom*twoJ[X_]*twoJ[Y_]*sqr(h_c);
  
  // h_20020
  sxt_2(-1e0/32e0, 1, 0, 0, 2, 2, 1, 0, 0, c, s);
  h_c = c;
  h_s = s;
  sxt_2(-1e0/32e0, 3, 0, 0, 0, 0, 1, 0, 2, c, s);
  h_c += c;
  h_s += s;
  sxt_2(-4e0/32e0, 1, 0, 0, 2, 1, 0, 1, 1, c, s);
  h_c += c;
  h_s += s;
  f += scl_geom*twoJ[X_]*twoJ[Y_]*sqr(h_c);

  // h_20200
  sxt_2(-1e0/32e0, 3, 0, 0, 0, 0, 1, 2, 0, c, s);
  h_c = c;
  h_s = s;
  sxt_2(-1e0/32e0, 1, 0, 2, 0, 2, 1, 0, 0, c, s);
  h_c += c;
  h_s += s;
  sxt_2(-4e0/32e0, 1, 0, 1, 1, 1, 0, 2, 0, c, s);
  h_c += c;
  h_s += s;
  f += scl_geom*twoJ[X_]*twoJ[Y_]*sqr(h_c);
  
  // h_11200
  sxt_2(-1e0/16e0, 1, 0, 2, 0, 1, 2, 0, 0, c, s);
  h_c = c;
  h_s = s;
  sxt_2(-1e0/16e0, 2, 1, 0, 0, 0, 1, 2, 0, c, s);
  h_c += c;
  h_s += s;
  sxt_2(-2e0/16e0, 0, 1, 1, 1, 1, 0, 2, 0, c, s);
  h_c += c;
  h_s += s;
  sxt_2(-2e0/16e0, 1, 0, 1, 1, 0, 1, 2, 0, c, s);
  h_c += c;
  h_s += s;
  f += scl_geom*twoJ[X_]*twoJ[Y_]*sqr(h_c);
    
  // h_00310
  sxt_2(-1e0/16e0, 0, 1, 2, 0, 1, 0, 1, 1, c, s);
  h_c = c;
  h_s = s;
  sxt_2(-1e0/16e0, 0, 1, 1, 1, 1, 0, 2, 0, c, s);
  h_c += c;
  h_s += s;
  f += scl_geom*sqr(twoJ[Y_])*sqr(h_c);
  
  // h_00400
  sxt_2(-1e0/32e0, 0, 1, 2, 0, 1, 0, 2, 0, h_c, h_s);
  f += scl_geom*twoJ[X_]*twoJ[Y_]*sqr(h_c);

  // h_22000
  sxt_2(-1e0/32e0, 0, 3, 0, 0, 3, 0, 0, 0, c, s);
  h_c = c;
  h_s = s;
  sxt_2(-3e0/32e0, 1, 2, 0, 0, 2, 1, 0, 0, c, s);
  h_c += c;
  h_s += s;
  f += scl_dnu*sqr(twoJ[X_])*sqr(h_c);
  
  // h_11110
  sxt_2(-1e0/8e0, 1, 0, 1, 1, 1, 2, 0, 0, c, s);
  h_c = c;
  h_s = s;
  sxt_2(-1e0/8e0, 2, 1, 0, 0, 0, 1, 1, 1, c, s);
  h_c += c;
  h_s += s;
  sxt_2(-1e0/8e0, 0, 1, 0, 2, 1, 0, 2, 0, c, s);
  h_c += c;
  h_s += s;
  sxt_2(-1e0/8e0, 1, 0, 0, 2, 0, 1, 2, 0, c, s);
  h_c += c;
  h_s += s;
  f += scl_dnu*twoJ[X_]*twoJ[Y_]*sqr(h_c);
  
  // h_00220
  sxt_2(-1e0/8e0, 0, 1, 1, 1, 1, 0, 1, 1, c, s);
  h_c = c;
  h_s = s;
  sxt_2(-1e0/32e0, 0, 1, 0, 2, 1, 0, 2, 0, c, s);
  h_c += c;
  h_s += s;
  sxt_2(-1e0/32e0, 0, 1, 2, 0, 1, 0, 0, 2, c, s);
  h_c += c;
  h_s += s;
  f += scl_dnu*sqr(twoJ[Y_])*sqr(h_c);

  if (iter % n_prt == 0) {
    printf("\n");
    printf("%4d %10.3e\n", iter, f);
    for (i = 1; i <= n_b3; i++)
      printf(" %9.6f", p[i]);
    printf("\n");
  }

  return f;
}


void h_min(void)
{
  const double ftol = 1e-6;

  int    i, j, n_iter;
  double *p, **xi, fret;

  p = dvector(1, n_b3_max);
  xi = dmatrix(1, n_b3_max, 1, n_b3_max);
  
  n_b3 = 0;
  b3s[n_b3++] = ElemIndex("s4");
  b3s[n_b3++] = ElemIndex("s6");
  b3s[n_b3++] = ElemIndex("s13");
  b3s[n_b3++] = ElemIndex("s19");
  b3s[n_b3++] = ElemIndex("s22");
  b3s[n_b3++] = ElemIndex("s24");

  for (i = 1; i <= n_b3; i++) {
    p[i] = GetKpar(b3s[i-1], 1, Sext);
    for (j = 1; j <= n_b3; j++)
      if (i == j)
	xi[i][j] = 1e-3;
      else
	xi[i][j] = 0e0;
  }

  iter = 0;
  dpowell(p, xi, n_b3, ftol, &n_iter, &fret, f_sxt);
  f_sxt(p);

  free_dvector(p, 1, n_b3_max);
  free_dmatrix(xi, 1, n_b3_max, 1, n_b3_max);
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
  const bool no_oct = !false;
  const int n_aper  = 25;

  double x_aper[n_aper], y_aper[n_aper];
  FILE   *fp;

  trace            = false;
  reverse_elem     = !false;
  globval.mat_meth = false;

  Read_Lattice(argv[1]);

  set_state();

  if (no_oct)
    no_mult(Oct);

  Ring_GetTwiss(true, 0e0);
  printglob();

  prtmfile("flat_file.dat");
  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true);

  sext_terms(twoJ);

  if (no_oct)
    printf("\n*** Octupoles turned off\n");

  if (false) {
    no_sxt();
    h_min();
    sext_terms(twoJ);
  }

  if (false) {
    fp = file_write("dynap.out");
    dynap
      (fp, 5e-3, 0e0, 0.1e-3, n_aper, 512, x_aper, y_aper, false, false, true);
    fclose(fp);
  }
}
