
/* Description:

   Control of on & off-momentum dynamic aperture.
   Evaluates & minimizes the Lie generators to 2nd order in the sextupole
   strengths.
   Initially coded in Pascal and ported to OPA (for the conceptual design of
   SLS). The latter is a derivatvie of K. Wille's OPTIK, DELTA, Dortmund Univ.
   (i.e., an interactive Linear Optics Design Tool; coded in Pascal).

   References:

   [1] J. Bengtsson "The sextupole scheme for the Swiss Light Source (SLS):
       An Analytic Approach" SLS Note 9/97.
   [2] A. Streun OPA.

   Johan Bengtsson                                                            */


#include "tracy_lib.h"

int              iter, n_b3;
double           twoJ_[2];
std::vector<int> b3s;


void propagate_drift(const double L, ss_vect<tps> &A)
{
  // Drift propagator.

  ss_vect<tps> Id, M;

  Id.identity();

  M.identity();
  M[x_] += L*Id[px_]; M[y_] += L*Id[py_];
  A = M*A;
}


void propagate_thin_kick(const double L, const double rho_inv,
			 const double b2, ss_vect<tps> &A)
{
  // Thin-kick propagator.

  ss_vect<tps> Id, M, A1;

  Id.identity();

  M.identity();
  M[px_] += -(b2+sqr(rho_inv))*L*Id[x_]; M[py_] += b2*L*Id[y_];
  M[px_] += rho_inv*L*Id[delta_];

  A = M*A;
}


void propagate_fringe_field(const double L, const double rho_inv,
			    const double phi, ss_vect<tps> &A)
{
  // Dipole harde-edge fringe field propagator.

  int          k;
  ss_vect<tps> Id;

  Id.identity();

  if (phi != 0e0)
    for (k = 0; k < 4; k++) {
      A[px_] += rho_inv*tan(phi*M_PI/180e0)*A[x_][k]*Id[k];
      A[py_] -= rho_inv*tan(phi*M_PI/180e0)*A[y_][k]*Id[k];
    }
}


void get_quad(const double L,
	      const double rho_inv, const double phi1, const double phi2,
	      const double b2,
	      const double alpha0[], const double beta0[], const double nu0[],
	      const double eta0[], const double etap0[],
	      const int m_x, const int m_y, const int n_x, const int n_y,
	      const int m, double &h_c, double &h_s)
{
  // Quad propagator: second order symplectic integrator.
  // Use Simpson's rule to integrate for thick magnets.

  int          j, k;
  double       h, r, phi;
  double       alpha[2], beta[2], nu[2], eta[2], etap[2], dnu[2];
  ss_vect<tps> A1;

  const bool prt    = false;
  const int  n_step = 25;

  for (k = 0; k <= 1; k++) {
    alpha[k] = alpha0[k]; beta[k] = beta0[k]; nu[k] = nu0[k];
    eta[k] = eta0[k]; etap[k] = etap0[k];
  }

  if (prt) {
    printf("\n  %7.3f %6.3f %5.3f %6.3f %6.3f %7.3f %6.3f %5.3f\n",
	   alpha[X_], beta[X_], nu[X_], eta[X_], etap[X_],
	   alpha[Y_], beta[Y_], nu[Y_]);
  }

  A1 = get_A(alpha, beta, eta, etap);
  A1 = get_A_CS(2, A1, dnu);

  propagate_fringe_field(L, rho_inv, phi1, A1);

  h = L/n_step; h_c = 0e0; h_s = 0e0;
  for (j = 1; j <= n_step; j++) {
    get_ab(A1, alpha, beta, dnu, eta, etap);
    for (k = 0; k <= 1; k++)
      nu[k] += dnu[k];
    A1 = get_A_CS(2, A1, dnu);

    phi = 2e0*M_PI*(n_x*nu[X_]+n_y*nu[Y_]);
    r = pow(beta[X_], m_x/2e0)*pow(beta[Y_], m_y/2e0)*pow(eta[X_], m-1);
    h_c += r*cos(phi); h_s += -r*sin(phi);

    propagate_drift(h/2e0, A1);
    
    if ((rho_inv != 0e0) || (b2 != 0e0)) {
      if (!true && prt)
	printf("     1/rho^2 = %5.3f, b2 = %5.3f\n", sqr(rho_inv), b2);
      propagate_thin_kick(h, rho_inv, b2, A1);
    }

    get_ab(A1, alpha, beta, dnu, eta, etap);
    for (k = 0; k <= 1; k++)
      nu[k] += dnu[k];
    A1 = get_A_CS(2, A1, dnu);
    
    if (prt) {
      printf("  %7.3f %6.3f %5.3f %6.3f %6.3f %7.3f %6.3f %5.3f\n",
	     alpha[X_], beta[X_], nu[X_], eta[X_], etap[X_],
	     alpha[Y_], beta[Y_], nu[Y_]);
      }

    phi = 2e0*M_PI*(n_x*nu[X_]+n_y*nu[Y_]);
    r = pow(beta[X_], m_x/2e0)*pow(beta[Y_], m_y/2e0)*pow(eta[X_], m-1);
    h_c += 4e0*r*cos(phi); h_s += -4e0*r*sin(phi);

    propagate_drift(h/2e0, A1);

    get_ab(A1, alpha, beta, dnu, eta, etap);
    for (k = 0; k <= 1; k++)
      nu[k] += dnu[k];
    A1 = get_A_CS(2, A1, dnu);
    
    if (prt) {
      printf("  %7.3f %6.3f %5.3f %6.3f %6.3f %7.3f %6.3f %5.3f\n",
	     alpha[X_], beta[X_], nu[X_], eta[X_], etap[X_],
	     alpha[Y_], beta[Y_], nu[Y_]);
      }

    phi = 2e0*M_PI*(n_x*nu[X_]+n_y*nu[Y_]);
    r = pow(beta[X_], m_x/2e0)*pow(beta[Y_], m_y/2e0)*pow(eta[X_], m-1);
    h_c += r*cos(phi); h_s += -r*sin(phi);

  }

  h_c *= b2*L/(6e0*n_step); h_s *= b2*L/(6e0*n_step);

  propagate_fringe_field(L, rho_inv, phi2, A1);

  get_ab(A1, alpha, beta, dnu, eta, etap);
  for (k = 0; k <= 1; k++)
    nu[k] += dnu[k];
  A1 = get_A_CS(2, A1, dnu);
  if (prt) {
    printf("  %7.3f %6.3f %5.3f %6.3f %6.3f %7.3f %6.3f %5.3f\n",
	   alpha[X_], beta[X_], nu[X_], eta[X_], etap[X_],
	   alpha[Y_], beta[Y_], nu[Y_]);
  }
}


void get_quad1(const double L,
	      const double rho_inv, const double phi1, const double phi2, 
	      const double b2,
	      double alpha3[][2], double beta3[][2], double nu3[][2],
	      double eta3[][2], double etap3[][2],
	      const int m_x, const int m_y, const int n_x, const int n_y,
	      const int m, double &h_c, double &h_s)
{
  // Use Simpson's Rule for quadrupoles.

  int    k;
  double A, phi;

//  get_Twiss3(L, rho_inv, phi1, phi2, b2, alpha3, beta3, nu3, eta3, etap3);
    
  h_c = 0e0; h_s = 0e0;
  for (k = 0; k <= 2; k++) {
    phi = 2e0*M_PI*(n_x*nu3[k][X_]+n_y*nu3[k][Y_]);
    A = pow(beta3[k][X_], m_x/2e0)*pow(beta3[k][Y_], m_y/2e0)
        *pow(eta3[k][X_], m-1);
    if (k == 1) A *= 4e0;
    h_c += A*cos(phi); h_s += -A*sin(phi);
  }
  if (m_x != 0) {
    h_c *= (b2+sqr(rho_inv))*L/6e0; h_s *= (b2+sqr(rho_inv))*L/6e0;
  } else {
    h_c *= b2*L/6e0; h_s *= b2*L/6e0;
  }
}


void get_h_ijklm(const int i, const int j, const int k, const int l,
		 const int m, ss_vect<tps> &A, double nu[], double &c,
		 double &s)
{
  int    k1;
  double alpha[2], beta[2], dnu[2], eta[2], etap[2], phi, ampl;

  const bool prt = false;

  get_ab(A, alpha, beta, dnu, eta, etap);
  A = get_A_CS(2, A, dnu);
  for (k1 = 0; k1 < 2; k1++)
    nu[k1] += dnu[k1];
  if (prt)
    printf("  %7.3f %6.3f %5.3f %6.3f %6.3f %7.3f %6.3f %5.3f\n",
	   alpha[X_], beta[X_], nu[X_], eta[X_], etap[X_],
	   alpha[Y_], beta[Y_], nu[Y_]);

  phi = 2e0*M_PI*((i-j)*nu[X_]+(k-l)*nu[Y_]); 
  ampl = pow(beta[X_], (i+j)/2e0)*pow(beta[Y_], (k+l)/2e0);
  if (m >= 1) {
    ampl *= -pow(eta[X_], m);
    if (m == 1) ampl *= 2e0;
  }

  c += ampl*cos(phi); s -= ampl*sin(phi);
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


void get_twiss(const int n, double alpha[], double beta[], double eta[],
	       double etap[], double nu[])
{
  int k;

  for (k = 0; k < 2; k++) {
    alpha[k] = Cell[n].Alpha[k];
    beta[k] = Cell[n].Beta[k];
    nu[k] = Cell[n].Nu[k];
    eta[k] = Cell[n].Eta[k];
    etap[k] = Cell[n].Etap[k];
  }
}


void sxt_1_(const double scl,
	   const int i, const int j, const int k, const int l, const int m,
	   double &h_c, double &h_s)
{
  // First order generators.

  int    n, k1;
  double c, s, b2, a2, b3, a3, b3L, a3L, A, phi;
  double alpha0[2], beta0[2], nu0[2], eta0[2], etap0[2], beta[2], nu[2], eta[2];

  const bool
    prt = false;
  const int  
    m_x = i+j,
    m_y = k+l,
    n_x = i-j,
    n_y = k-l;

  h_c = h_s = 0e0;
  if (prt) printf("sxt_1:\n");
  for (n = 0; n <= globval.Cell_nLoc; n++) {
    if (Cell[n].Elem.Pkind == Mpole) {
      if (((Cell[n].Elem.M->Pirho != 0e0) ||
    	   (Cell[n].Elem.M->Porder == Quad)) && (m >= 1)) {
	if (prt) printf("  %4d %10s quad\n", n, Cell[n].Elem.PName);
	get_bn_design_elem(Cell[n].Fnum, 1, Quad, b2, a2);
	get_twiss(ind(n-1), alpha0, beta0, eta0, etap0, nu0);
    	get_quad(Cell[n].Elem.PL, Cell[n].Elem.M->Pirho,
    		 Cell[n].Elem.M->PTx1, Cell[n].Elem.M->PTx2,
    		 b2, alpha0, beta0, nu0, eta0, etap0,
    		 m_x, m_y, n_x, n_y, m, c, s);
    	h_c += scl*c; h_s += scl*s;
      } else if (Cell[n].Elem.M->Porder >= Sext) {
    	if (prt) printf("  %4d %10s sext\n", n, Cell[n].Elem.PName);
	get_bnL_design_elem(Cell[n].Fnum, 1, Sext, b3L, a3L);
	get_twiss(ind(n-1), alpha0, beta0, eta0, etap0, nu0);
    	for (k1 = 0; k1 <= 1; k1++) {
    	  beta[k1] = beta0[k1]; nu[k1] = nu0[k1];
    	  eta[k1] = eta0[k1];
    	}
	phi = 2e0*M_PI*(n_x*nu[X_]+n_y*nu[Y_]); 
    	A = scl*b3L*pow(beta[X_], m_x/2e0)*pow(beta[Y_], m_y/2e0);
    	if (m >= 1) {
    	  A *= -pow(eta[X_], m);
    	  if (m == 1) A *= 2e0;
    	}
    	h_c += A*cos(phi); h_s -= A*sin(phi);
      }
    }
  }
}


void get_lin_map(const int n_step, const double delta, const elemtype &Elem, 
		 ss_vect<tps> &M1, ss_vect<tps> &M2, ss_vect<tps> &M)
{

  if (Elem.PL != 0e0) {
    M1 = get_edge_lin_map(Elem.M->Pirho, Elem.M->PTx1, Elem.M->Pgap, delta);
    M2 = get_edge_lin_map(Elem.M->Pirho, Elem.M->PTx2, Elem.M->Pgap, delta);
    M =
      get_sbend_lin_map(Elem.PL/n_step, Elem.M->Pirho, Elem.M->PB[Quad+HOMmax],
			delta);
  } else {
    M1.identity(); M2.identity();
    M = get_thin_kick_lin_map(Elem.PL*Elem.M->PB[Quad+HOMmax], delta);
  }
}


void get_mult(const int loc, const int n_step, const double delta,
	      const int i, const int j, const int k, const int l, const int m,
	      double &h_c, double &h_s)
{
  int          k1;
  double       nu[2];
  ss_vect<tps> M1, M2, M, A;

  const bool prt = false;

  A =
    get_A(Cell[ind(loc-1)].Alpha, Cell[ind(loc-1)].Beta,
	  Cell[ind(loc-1)].Eta, Cell[ind(loc-1)].Etap);
  for (k1 = 0; k1 < 2; k1++)
    nu[k1] = Cell[ind(loc-1)].Nu[k1];

  get_lin_map(n_step-1, delta, Cell[loc].Elem, M1, M2, M);

  if (prt) printf("\n");
  h_c = h_s = 0e0;
  get_h_ijklm(i, j, k, l, m, A, nu, h_c, h_s);
  A = M1*A;
  for (k1 = 0; k1 < n_step-2; k1++) {
    A = M*A;
    get_h_ijklm(i, j, k, l, m, A, nu, h_c, h_s);
  }
  A = M2*M*A;
  get_h_ijklm(i, j, k, l, m, A, nu, h_c, h_s);
  h_c /= n_step; h_s /= n_step;
}


void sxt_1(const double scl,
	   const int i, const int j, const int k, const int l, const int m,
	   double &h_c, double &h_s)
{
  // First order generators.

  int    n, k1;
  double c, s, b2, a2, b3, a3, b3L, a3L, A, phi;
  double alpha0[2], beta0[2], nu0[2], eta0[2], etap0[2];

  const bool
    prt = false;
  const int
    n_step = 3;
  const int  
    m_x = i+j,
    m_y = k+l,
    n_x = i-j,
    n_y = k-l;

  h_c = h_s = 0e0;
  if (prt) printf("sxt_1:\n");
  for (n = 0; n <= globval.Cell_nLoc; n++) {
    if (Cell[n].Elem.Pkind == Mpole) {
      if (((Cell[n].Elem.M->Pirho != 0e0) ||
    	   (Cell[n].Elem.M->Porder == Quad)) && (m >= 1)) {
	if (prt) printf("  %4d %10s quad\n", n, Cell[n].Elem.PName);
	get_bn_design_elem(Cell[n].Fnum, 1, Quad, b2, a2);
	get_twiss(ind(n-1), alpha0, beta0, eta0, etap0, nu0);
    	get_quad(Cell[n].Elem.PL, Cell[n].Elem.M->Pirho,
    		 Cell[n].Elem.M->PTx1, Cell[n].Elem.M->PTx2,
    		 b2, alpha0, beta0, nu0, eta0, etap0,
    		 m_x, m_y, n_x, n_y, m, c, s);
    	h_c += scl*c; h_s += scl*s;
      } else if (Cell[n].Elem.M->Porder >= Sext) {
    	if (prt) printf("  %4d %10s sext\n", n, Cell[n].Elem.PName);
	get_bnL_design_elem(Cell[n].Fnum, 1, Sext, b3L, a3L);
	get_mult(n, n_step, 0e0, i, j, k, l, m, c, s);
    	h_c += scl*b3L*c; h_s += scl*b3L*s;
      }
    }
  }
}


void sxt_2(const double twoJx, const double twoJy, const double scl,
	   const int i1, const int j1, const int k1, const int l1,
	   const int i2, const int j2, const int k2, const int l2,
	   double &h_c, double &h_s)
{
  // Second order generators.

  long int n1, n2, k;
  double   b3L1, a3L1, b3L2, a3L2, dnu, A1, A, phi;
  double   alpha0[2], beta0[2], nu0[2], eta0[2], etap0[2];
  double   beta1[2], nu1[2], beta2[2], nu2[2];

  h_c = 0e0; h_s = 0e0;
  for (n1 = 0; n1 <= globval.Cell_nLoc; n1++) {
    if ((Cell[n1].Elem.Pkind == Mpole) &&
	(Cell[n1].Elem.M->Porder >= Sext)) {
      for (k = 0; k <= 1; k++) {
	alpha0[k] = Cell[ind(n1-1)].Alpha[k];
	beta0[k] = Cell[ind(n1-1)].Beta[k];
	nu0[k] = Cell[ind(n1-1)].Nu[k];
	eta0[k] = Cell[ind(n1-1)].Eta[k]; etap0[k] = Cell[ind(n1-1)].Etap[k];
      }
      for (k = 0; k <= 1; k++) {
	beta1[k] = beta0[k]; nu1[k] = nu0[k];
      }
      get_bnL_design_elem(Cell[n1].Fnum, 1, Sext, b3L1, a3L1);
      dnu = (i1-j1)*nu1[X_] + (k1-l1)*nu1[Y_]; 
      A1 = b3L1*pow(beta1[X_], (i1+j1)/2e0)*pow(twoJy*beta1[Y_], (k1+l1)/2e0);
      for (n2 = 0; n2 <= globval.Cell_nLoc; n2++) {
	if ((Cell[n2].Elem.Pkind == Mpole) &&
	    (Cell[n2].Elem.M->Porder >= Sext)) {
	  for (k = 0; k <= 1; k++) {
	    alpha0[k] = Cell[n2-1].Alpha[k];
	    beta0[k] = Cell[n2-1].Beta[k]; nu0[k] = Cell[n2-1].Nu[k];
	    eta0[k] = Cell[n2-1].Eta[k]; etap0[k] = Cell[n2-1].Etap[k];
	  }
	  for (k = 0; k <= 1; k++) {
	    beta2[k] = beta0[k]; nu2[k] = nu0[k];
	  }
	  get_bnL_design_elem(Cell[n2].Fnum, 1, Sext, b3L2, a3L2);
	  phi = 2e0*M_PI*(dnu+(i2-j2)*nu2[X_]+(k2-l2)*nu2[Y_]);
	  A = scl*A1*b3L2*pow(beta2[X_], (i2+j2)/2e0)
		 *pow(beta2[Y_], (k2+l2)/2e0);
	  // [..., ...] -> i*...
	  if (n2 < n1) {
	    h_s += A*cos(phi); h_c += A*sin(phi);
	  } else if (n2 > n1) {
	    h_s -= A*cos(phi); h_c -= A*sin(phi);
	  }
	}
      }
    }
  }
}


void K(const double nu_x, const double nu_y,
       const double twoJx, const double twoJy, double a[])
{
  // Amplitude dependent tune shifts.

  long int k, n1, n2;
  double   b3L1, a3L1, b3L2, a3L2, pi_1x, pi_3x, pi_1xm2y, pi_1xp2y;
  double   s_1x, s_3x, s_1xm2y, s_1xp2y, c_1x, c_3x, c_1xm2y, c_1xp2y;
  double   dmu_x, dmu_y, A;
  double   alpha0[2], beta0[2], nu0[2], eta0[2], etap0[2];
  double   beta1[2], nu1[2], beta2[2], nu2[2];

  pi_1x = M_PI*nu_x; pi_3x = 3e0*M_PI*nu_x; pi_1xm2y = M_PI*(nu_x-2e0*nu_y);
  pi_1xp2y = M_PI*(nu_x+2e0*nu_y);

  s_1x = sin(pi_1x); s_3x = sin(pi_3x); s_1xm2y = sin(pi_1xm2y);
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
	beta1[k] = beta0[k]; nu1[k] = nu0[k];
      }
      get_bnL_design_elem(Cell[n1].Fnum, 1, Sext, b3L1, a3L1);

      for (n2 = 0; n2 <= globval.Cell_nLoc; n2++) {
	if ((Cell[n2].Elem.Pkind == Mpole) &&
	    (Cell[n2].Elem.M->Porder >= Sext)) {
	  for (k = 0; k <= 1; k++) {
	    alpha0[k] = Cell[n2-1].Alpha[k];
	    beta0[k] = Cell[n2-1].Beta[k]; nu0[k] = Cell[n2-1].Nu[k];
	    eta0[k] = Cell[n2-1].Eta[k]; etap0[k] = Cell[n2-1].Etap[k];
	  }
	  for (k = 0; k <= 1; k++) {
	    beta2[k] = beta0[k]; nu2[k] = nu0[k];
	  }
	  get_bnL_design_elem(Cell[n2].Fnum, 1, Sext, b3L2, a3L2);
	  dmu_x = 2e0*M_PI*fabs(nu1[X_]-nu2[X_]);
	  dmu_y = 2e0*M_PI*fabs(nu1[Y_]-nu2[Y_]);
	  A = b3L1*b3L2*sqrt(beta1[X_]*beta2[X_]);
	  c_1x = cos(dmu_x-pi_1x)/s_1x; c_3x = cos(3e0*dmu_x-pi_3x)/s_3x;
	  c_1xm2y = cos(dmu_x-2e0*dmu_y-pi_1xm2y)/s_1xm2y;
	  c_1xp2y = cos(dmu_x+2e0*dmu_y-pi_1xp2y)/s_1xp2y;
	  a[0] += A*beta1[X_]*beta2[X_]*(3*c_1x+c_3x)/2e0;
	  a[1] -= A*beta1[Y_]*(4e0*twoJx*beta2[X_]*c_1x+2e0*twoJy*beta2[Y_]
		  *(c_1xm2y-c_1xp2y));
	  a[2] += A*beta1[Y_]*beta2[Y_]*(4e0*c_1x+c_1xm2y+c_1xp2y)/2e0;
	}
      }
    }
  }

  for (k = 0; k < 3; k++)
    a[k] /= 32e0;
}


void get_nu_(Matrix &M, double nu[])
{
  int    k;
  double tr2;

  const bool prt = false;

  for (k = 0; k < 2; k++) {
    tr2 = globval.OneTurnMat[2*k][2*k] + globval.OneTurnMat[2*k+1][2*k+1];
    nu[k] = acos(tr2/2e0)/(2e0*M_PI);
    if (globval.OneTurnMat[2*k][2*k+1] < 0e0) nu[k] = 1e0 - nu[k];
  }

  if (prt) printf("\nget_nu: [%10.8f, %10.8f]\n", nu[X_], nu[Y_]);
}


void get_ksi2_(const double delta, double ksi2[])
{
  long int lastpos;
  int      k;
  double   nu[3][2];

  const bool prt = false;

  getcod(-delta, lastpos);
  get_nu_(globval.OneTurnMat, nu[0]);
  getcod(0e0, lastpos);
  get_nu_(globval.OneTurnMat, nu[1]);
  getcod(delta, lastpos);
  get_nu_(globval.OneTurnMat, nu[2]);
  for (k = 0; k < 2; k++)
    ksi2[k] = (nu[2][k]-2e0*nu[1][k]+nu[0][k])/(2e0*sqr(delta));

  if (prt) printf("\nget_ksi2_: [%10.8f, %10.8f]\n", ksi2[X_], ksi2[Y_]);
}


void get_cross_terms(const double delta, const double twoJx, const double twoJy,
		     double ad[])
{
  int    k;
  double a[2][3];

  Ring_GetTwiss(false, -delta);
  K(globval.TotalTune[X_], globval.TotalTune[Y_], twoJx, twoJy, a[0]);
  Ring_GetTwiss(false, delta);
  K(globval.TotalTune[X_], globval.TotalTune[Y_], twoJx, twoJy, a[1]);

  for (k = 0; k < 3; k++)
    ad[k] = (a[1][k]-a[0][k])/(2e0*delta);
}


void drv_terms(const double delta, const double twoJx, const double twoJy)
{
  double h_c, h_s, c, s, a[3], ksi2[2];

  printf("\nFirst order chromatic terms:\n");
  sxt_1(1e0/2e0, 1, 1, 0, 0, 1, h_c, h_s);
  printf("  h_11001 = %12.5e ksi_x = %7.3f\n", h_c/2e0, -h_c/(2e0*M_PI));
  sxt_1(-1e0/2e0, 0, 0, 1, 1, 1, h_c, h_s);
  printf("  h_00111 = %12.5e ksi_y = %7.3f\n", h_c/2e0, -h_c/(2e0*M_PI));

  sxt_1(1e0/4e0, 2, 0, 0, 0, 1, h_c, h_s);
  h_c *= twoJx; h_s *= twoJx;
  printf("\n  h_20001 = [%12.5e, %12.5e]\n", h_c/2e0, h_s/2e0);
  sxt_1(-1e0/4e0, 0, 0, 2, 0, 1, h_c, h_s);
  h_c *= twoJy; h_s *= twoJy;
  printf("  h_00201 = [%12.5e, %12.5e]\n", h_c/2e0, h_s/2e0);
  sxt_1(1e0, 1, 0, 0, 0, 2, h_c, h_s);
  h_c *= sqrt(twoJx); h_s *= sqrt(twoJx);
  printf("  h_10002 = [%12.5e, %12.5e]\n", h_c/2e0, h_s/2e0);

  printf("\nFirst order geometric terms:\n");
  sxt_1(-1e0/12e0, 3, 0, 0, 0, 0, h_c, h_s);
  h_c *= pow(twoJx, 3e0/2e0); h_s *= pow(twoJx, 3e0/2e0);
  printf("  h_30000 = [%12.5e, %12.5e]\n", h_c/2e0, h_s/2e0);
  sxt_1(-1e0/4e0, 2, 1, 0, 0, 0, h_c, h_s);
  h_c *= pow(twoJx, 3e0/2e0); h_s *= pow(twoJx, 3e0/2e0);
  printf("  h_21000 = [%12.5e, %12.5e]\n", h_c/2e0, h_s/2e0);
  sxt_1(1e0/4e0, 1, 0, 2, 0, 0, h_c, h_s);
  h_c *= sqrt(twoJx)*twoJy; h_s *= sqrt(twoJx)*twoJy;
  printf("  h_10200 = [%12.5e, %12.5e]\n", h_c/2e0, h_s/2e0);
  sxt_1(1e0/2e0, 1, 0, 1, 1, 0, h_c, h_s);
  h_c *= sqrt(twoJx)*twoJy; h_s *= sqrt(twoJx)*twoJy;
  printf("  h_10110 = [%12.5e, %12.5e]\n", h_c/2e0, h_s/2e0);
  sxt_1(1e0/4e0, 1, 0, 0, 2, 0, h_c, h_s);
  h_c *= sqrt(twoJx)*twoJy; h_s *= sqrt(twoJx)*twoJy;
  printf("  h_10020 = [%12.5e, %12.5e]\n", h_c/2e0, h_s/2e0);
 
  printf("\nSecond order geometric terms:\n");
  sxt_2(twoJx, twoJy, -1e0/32e0, 2, 1, 0, 0, 3, 0, 0, 0, h_c, h_s);
  h_c *= sqr(twoJx); h_s *= sqr(twoJx);
  printf("  h_40000 = [%12.5e, %12.5e]\n", h_c/2e0, h_s/2e0);
  
  sxt_2(twoJx, twoJy, -1e0/16e0, 1, 2, 0, 0, 3, 0, 0, 0, h_c, h_s);
  h_c *= sqr(twoJx); h_s *= sqr(twoJx);
  printf("  h_31000 = [%12.5e, %12.5e]\n", h_c/2e0, h_s/2e0);

  sxt_2(twoJx, twoJy, -1e0/32e0, 0, 3, 0, 0, 3, 0, 0, 0, c, s);
  h_c = c; h_s = s;
  sxt_2(twoJx, twoJy, -3e0/32e0, 1, 2, 0, 0, 2, 1, 0, 0, c, s);
  h_c += c; h_s += s;
  h_c *= sqr(twoJx); h_s *= sqr(twoJx);
  printf("  h_22000 = [%12.5e, %12.5e]\n", h_c/2e0, h_s/2e0);

  sxt_2(twoJx, twoJy, -1e0/32e0, 3, 0, 0, 0, 0, 1, 2, 0, c, s);
  h_c = c; h_s = s;
  sxt_2(twoJx, twoJy, -1e0/32e0, 1, 0, 2, 0, 2, 1, 0, 0, c, s);
  h_c += c; h_s += s;
  sxt_2(twoJx, twoJy, -4e0/32e0, 1, 0, 1, 1, 1, 0, 2, 0, c, s);
  h_c += c; h_s += s;
  h_c *= twoJx*twoJy; h_s *= twoJx*twoJy;
  printf("  h_20200 = [%12.5e, %12.5e]\n", h_c/2e0, h_s/2e0);
  
  sxt_2(twoJx, twoJy, -1e0/16e0, 1, 0, 2, 0, 1, 2, 0, 0, c, s);
  h_c = c; h_s = s;
  sxt_2(twoJx, twoJy, -1e0/16e0, 2, 1, 0, 0, 0, 1, 2, 0, c, s);
  h_c += c; h_s += s;
  sxt_2(twoJx, twoJy, -2e0/16e0, 0, 1, 1, 1, 1, 0, 2, 0, c, s);
  h_c += c; h_s += s;
  sxt_2(twoJx, twoJy, -2e0/16e0, 1, 0, 1, 1, 0, 1, 2, 0, c, s);
  h_c += c; h_s += s;
  h_c *= twoJx*twoJy; h_s *= twoJx*twoJy;
  printf("  h_11200 = [%12.5e, %12.5e]\n", h_c/2e0, h_s/2e0);
    
  sxt_2(twoJx, twoJy, -1e0/16e0, 3, 0, 0, 0, 0, 1, 1, 1, c, s);
  h_c = c; h_s = s;
  sxt_2(twoJx, twoJy, -1e0/16e0, 1, 0, 1, 1, 2, 1, 0, 0, c, s);
  h_c += c; h_s += s;
  sxt_2(twoJx, twoJy, -1e0/8e0, 1, 0, 0, 2, 1, 0, 2, 0, c, s);
  h_c += c; h_s += s;
  h_c *= twoJx*twoJy; h_s *= twoJx*twoJy;
  printf("  h_20110 = [%12.5e, %12.5e]\n", h_c/2e0, h_s/2e0);
  
  sxt_2(twoJx, twoJy, -1e0/8e0, 1, 0, 1, 1, 1, 2, 0, 0, c, s);
  h_c = c; h_s = s;
  sxt_2(twoJx, twoJy, -1e0/8e0, 2, 1, 0, 0, 0, 1, 1, 1, c, s);
  h_c += c; h_s += s;
  sxt_2(twoJx, twoJy, -1e0/8e0, 0, 1, 0, 2, 1, 0, 2, 0, c, s);
  h_c += c; h_s += s;
  sxt_2(twoJx, twoJy, -1e0/8e0, 1, 0, 0, 2, 0, 1, 2, 0, c, s);
  h_c += c; h_s += s;
  h_c *= twoJx*twoJy; h_s *= twoJx*twoJy;
  printf("  h_11110 = [%12.5e, %12.5e]\n", h_c/2e0, h_s/2e0);
  
  sxt_2(twoJx, twoJy, -1e0/16e0, 0, 1, 2, 0, 1, 0, 1, 1, c, s);
  h_c = c; h_s = s;
  sxt_2(twoJx, twoJy, -1e0/16e0, 0, 1, 1, 1, 1, 0, 2, 0, c, s);
  h_c += c; h_s += s;
  h_c *= sqr(twoJy); h_s *= sqr(twoJy);
  printf("  h_00310 = [%12.5e, %12.5e]\n", h_c/2e0, h_s/2e0);
  
  sxt_2(twoJx, twoJy, -1e0/32e0, 1, 0, 0, 2, 2, 1, 0, 0, c, s);
  h_c = c; h_s = s;
  sxt_2(twoJx, twoJy, -1e0/32e0, 3, 0, 0, 0, 0, 1, 0, 2, c, s);
  h_c += c; h_s += s;
  sxt_2(twoJx, twoJy, -4e0/32e0, 1, 0, 0, 2, 1, 0, 1, 1, c, s);
  h_c += c; h_s += s;
  h_c *= twoJx*twoJy; h_s *= twoJx*twoJy;
  printf("  h_20020 = [%12.5e, %12.5e]\n", h_c/2e0, h_s/2e0);

  sxt_2(twoJx, twoJy, -1e0/32e0, 0, 1, 2, 0, 1, 0, 2, 0, h_c, h_s);
  h_c *= sqr(twoJy); h_s *= sqr(twoJy);
  printf("  h_00400 = [%12.5e, %12.5e]\n", h_c/2e0, h_s/2e0);

  sxt_2(twoJx, twoJy, -1e0/8e0, 0, 1, 1, 1, 1, 0, 1, 1, c, s);
  h_c = c; h_s = s;
  sxt_2(twoJx, twoJy, -1e0/32e0, 0, 1, 0, 2, 1, 0, 2, 0, c, s);
  h_c += c; h_s += s;
  sxt_2(twoJx, twoJy, -1e0/32e0, 0, 1, 2, 0, 1, 0, 0, 2, c, s);
  h_c += c; h_s += s;
  h_c *= sqr(twoJy); h_s *= sqr(twoJy);
  printf("  h_00220 = [%12.5e, %12.5e]\n", h_c/2e0, h_s/2e0);

  K(globval.TotalTune[X_], globval.TotalTune[Y_], twoJx, twoJy, a);
  printf("\nSecond order anharmonic terms:\n");
  printf("  k_22000 = %12.5e\n", a[0]);
  printf("  k_11110 = %12.5e\n", a[1]);
  printf("  k_00220 = %12.5e\n", a[2]);

  get_ksi2_(delta, ksi2);
  printf("\nSecond order chromaticity:\n");
  printf("  k_11002 = %12.5e\n", -M_PI*ksi2[X_]);
  printf("  k_00112 = %12.5e\n", -M_PI*ksi2[Y_]);

  get_cross_terms(delta, twoJx, twoJy, a);
  printf("\nSecond order cross terms:\n");
  printf("  k_22001 = %12.5e\n", a[0]);
  printf("  k_11111 = %12.5e\n", a[1]);
  printf("  k_00221 = %12.5e\n", a[2]);
}


double f_sxt(double p[])
{
  int    i;
  double f, h_c, h_s, c, s;

  const int    n_prt = 100;
  const double scl_ksi1 = 1e3, scl_geom = 1e0, scl_dnu = 1e0, b3L_max = 3.5;

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
  f += scl_geom*pow(twoJ_[X_], 3e0/2e0)*sqr(h_c);
  // h_30000
  sxt_1(-1e0/12e0, 3, 0, 0, 0, 0, h_c, h_s);
  f += scl_geom*pow(twoJ_[X_], 3e0/2e0)*sqr(h_c);
  // h_10110
  sxt_1(1e0/2e0, 1, 0, 1, 1, 0, h_c, h_s);
  f += scl_geom*sqrt(twoJ_[X_])*twoJ_[Y_]*sqr(h_c);
  // h_10020
  sxt_1(1e0/4e0, 1, 0, 0, 2, 0, h_c, h_s);
  f += scl_geom*sqrt(twoJ_[X_])*twoJ_[Y_]*sqr(h_c);
  // h_10200
  sxt_1(1e0/4e0, 1, 0, 2, 0, 0, h_c, h_s);
  f += scl_geom*sqrt(twoJ_[X_])*twoJ_[Y_]*sqr(h_c);

  // h_40000
  sxt_2(twoJ_[X_], twoJ_[Y_], -1e0/32e0, 2, 1, 0, 0, 3, 0, 0, 0, h_c, h_s);
  f += scl_geom*sqr(twoJ_[X_])*sqr(h_c);
  
  // h_31000
  sxt_2(twoJ_[X_], twoJ_[Y_], -1e0/16e0, 1, 2, 0, 0, 3, 0, 0, 0, h_c, h_s);
  f += scl_geom*sqr(twoJ_[X_])*sqr(h_c);

  // h_20110
  sxt_2(twoJ_[X_], twoJ_[Y_], -1e0/16e0, 3, 0, 0, 0, 0, 1, 1, 1, c, s);
  h_c = c; h_s = s;
  sxt_2(twoJ_[X_], twoJ_[Y_], -1e0/16e0, 1, 0, 1, 1, 2, 1, 0, 0, c, s);
  h_c += c; h_s += s;
  sxt_2(twoJ_[X_], twoJ_[Y_], -1e0/8e0, 1, 0, 0, 2, 1, 0, 2, 0, c, s);
  h_c += c; h_s += s;
  f += scl_geom*twoJ_[X_]*twoJ_[Y_]*sqr(h_c);
  
  // h_20020
  sxt_2(twoJ_[X_], twoJ_[Y_], -1e0/32e0, 1, 0, 0, 2, 2, 1, 0, 0, c, s);
  h_c = c; h_s = s;
  sxt_2(twoJ_[X_], twoJ_[Y_], -1e0/32e0, 3, 0, 0, 0, 0, 1, 0, 2, c, s);
  h_c += c; h_s += s;
  sxt_2(twoJ_[X_], twoJ_[Y_], -4e0/32e0, 1, 0, 0, 2, 1, 0, 1, 1, c, s);
  h_c += c; h_s += s;
  f += scl_geom*twoJ_[X_]*twoJ_[Y_]*sqr(h_c);

  // h_20200
  sxt_2(twoJ_[X_], twoJ_[Y_], -1e0/32e0, 3, 0, 0, 0, 0, 1, 2, 0, c, s);
  h_c = c; h_s = s;
  sxt_2(twoJ_[X_], twoJ_[Y_], -1e0/32e0, 1, 0, 2, 0, 2, 1, 0, 0, c, s);
  h_c += c; h_s += s;
  sxt_2(twoJ_[X_], twoJ_[Y_], -4e0/32e0, 1, 0, 1, 1, 1, 0, 2, 0, c, s);
  h_c += c; h_s += s;
  f += scl_geom*twoJ_[X_]*twoJ_[Y_]*sqr(h_c);
  
  // h_11200
  sxt_2(twoJ_[X_], twoJ_[Y_], -1e0/16e0, 1, 0, 2, 0, 1, 2, 0, 0, c, s);
  h_c = c; h_s = s;
  sxt_2(twoJ_[X_], twoJ_[Y_], -1e0/16e0, 2, 1, 0, 0, 0, 1, 2, 0, c, s);
  h_c += c; h_s += s;
  sxt_2(twoJ_[X_], twoJ_[Y_], -2e0/16e0, 0, 1, 1, 1, 1, 0, 2, 0, c, s);
  h_c += c; h_s += s;
  sxt_2(twoJ_[X_], twoJ_[Y_], -2e0/16e0, 1, 0, 1, 1, 0, 1, 2, 0, c, s);
  h_c += c; h_s += s;
  f += scl_geom*twoJ_[X_]*twoJ_[Y_]*sqr(h_c);
    
  // h_00310
  sxt_2(twoJ_[X_], twoJ_[Y_], -1e0/16e0, 0, 1, 2, 0, 1, 0, 1, 1, c, s);
  h_c = c; h_s = s;
  sxt_2(twoJ_[X_], twoJ_[Y_], -1e0/16e0, 0, 1, 1, 1, 1, 0, 2, 0, c, s);
  h_c += c; h_s += s;
  f += scl_geom*sqr(twoJ_[Y_])*sqr(h_c);
  
  // h_00400
  sxt_2(twoJ_[X_], twoJ_[Y_], -1e0/32e0, 0, 1, 2, 0, 1, 0, 2, 0, h_c, h_s);
  f += scl_geom*twoJ_[X_]*twoJ_[Y_]*sqr(h_c);

  // h_22000
  sxt_2(twoJ_[X_], twoJ_[Y_], -1e0/32e0, 0, 3, 0, 0, 3, 0, 0, 0, c, s);
  h_c = c; h_s = s;
  sxt_2(twoJ_[X_], twoJ_[Y_], -3e0/32e0, 1, 2, 0, 0, 2, 1, 0, 0, c, s);
  h_c += c; h_s += s;
  f += scl_dnu*sqr(twoJ_[X_])*sqr(h_c);
  
  // h_11110
  sxt_2(twoJ_[X_], twoJ_[Y_], -1e0/8e0, 1, 0, 1, 1, 1, 2, 0, 0, c, s);
  h_c = c; h_s = s;
  sxt_2(twoJ_[X_], twoJ_[Y_], -1e0/8e0, 2, 1, 0, 0, 0, 1, 1, 1, c, s);
  h_c += c; h_s += s;
  sxt_2(twoJ_[X_], twoJ_[Y_], -1e0/8e0, 0, 1, 0, 2, 1, 0, 2, 0, c, s);
  h_c += c; h_s += s;
  sxt_2(twoJ_[X_], twoJ_[Y_], -1e0/8e0, 1, 0, 0, 2, 0, 1, 2, 0, c, s);
  h_c += c; h_s += s;
  f += scl_dnu*twoJ_[X_]*twoJ_[Y_]*sqr(h_c);
  
  // h_00220
  sxt_2(twoJ_[X_], twoJ_[Y_], -1e0/8e0, 0, 1, 1, 1, 1, 0, 1, 1, c, s);
  h_c = c; h_s = s;
  sxt_2(twoJ_[X_], twoJ_[Y_], -1e0/32e0, 0, 1, 0, 2, 1, 0, 2, 0, c, s);
  h_c += c; h_s += s;
  sxt_2(twoJ_[X_], twoJ_[Y_], -1e0/32e0, 0, 1, 2, 0, 1, 0, 0, 2, c, s);
  h_c += c; h_s += s;
  f += scl_dnu*sqr(twoJ_[Y_])*sqr(h_c);

  if (iter % n_prt == 0) {
    printf("\n");
    printf("%4d %10.3e\n", iter, f);
    for (i = 1; i <= n_b3; i++)
      printf(" %9.6f", p[i]);
    printf("\n");
  }

  return f;
}


void h_min(const double twoJx, const double twoJy)
{
  int    i, j, n_iter;
  double *p, **xi, fret;

  const double ftol = 1e-6;

  b3s.push_back(ElemIndex("s4"));  b3s.push_back(ElemIndex("s6"));
  b3s.push_back(ElemIndex("s13")); b3s.push_back(ElemIndex("s19"));
  b3s.push_back(ElemIndex("s22")); b3s.push_back(ElemIndex("s24"));
  n_b3 = b3s.size();

  p = dvector(1, n_b3); xi = dmatrix(1, n_b3, 1, n_b3);
  
  twoJ_[X_] = twoJx; twoJ_[Y_] = twoJy;

  for (i = 1; i <= n_b3; i++) {
    p[i] = GetKpar(b3s[i-1], 1, Sext);
    for (j = 1; j <= n_b3; j++)
      xi[i][j] = (i == j)? 1e-3 : 0e0;
  }

  iter = 0;
  dpowell(p, xi, n_b3, ftol, &n_iter, &fret, f_sxt); f_sxt(p);

  free_dvector(p, 1, n_b3); free_dmatrix(xi, 1, n_b3, 1, n_b3);
}
