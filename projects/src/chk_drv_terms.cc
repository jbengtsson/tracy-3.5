#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


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
  const int  n_step = 3;

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


void sxt_1(const double scl, const double twoJ[], const double delta, 
	   const int i, const int j, const int k, const int l, const int m,
	   double &h_c, double &h_s, const bool incl_b3, const int loc)
{
  // First order generators.

  int    n;
  double c, s, b2, a2, b3L, a3L, scl1;
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
  scl1 = scl*pow(twoJ[X_], m_x/2e0)*pow(twoJ[Y_], m_y/2e0)*pow(delta, m);
  if (prt) printf("sxt_1:\n");
  for (n = 0; n <= loc; n++) {
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
    	h_c += scl1*c; h_s += scl1*s;
      } else if (incl_b3 && (Cell[n].Elem.M->Porder >= Sext)) {
    	if (prt) printf("  %4d %10s sext\n", n, Cell[n].Elem.PName);
	get_bnL_design_elem(Cell[n].Fnum, 1, Sext, b3L, a3L);
	get_mult(n, n_step, 0e0, i, j, k, l, m, c, s);
    	h_c += scl1*b3L*c; h_s += scl1*b3L*s;
      }
    }
  }
}


void chk_drv_terms()
{
  int    k;
  double h_c, h_s, h[3];
  FILE   *outf;

  const double twoJ[] = {1e0, 1e0}, delta_max = 1e0;

  outf = file_write("sxt.out");

  for (k = 0; k <= globval.Cell_nLoc; k++) {
    sxt_1(1.0/4.0, twoJ, delta_max, 2, 0, 0, 0, 1, h_c, h_s, true, k);
    h[0] = sqrt(sqr(h_c)+sqr(h_s));
    sxt_1(-1.0/4.0, twoJ, delta_max, 0, 0, 2, 0, 1, h_c, h_s, true, k);
    h[1] = sqrt(sqr(h_c)+sqr(h_s));
    sxt_1(1.0, twoJ, delta_max, 1, 0, 0, 0, 2, h_c, h_s, true, k);
    h[2] = sqrt(sqr(h_c)+sqr(h_s));
    fprintf(outf, "%8.3f %10.3e %10.3e %10.3e\n", Cell[k].S, h[0], h[1], h[2]);
  }

  fclose(outf);
}


int main(int argc, char *argv[])
{

  reverse_elem = !false;

  trace = false;

  globval.mat_meth = false;

  if (!true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  globval.H_exact    = false; globval.quad_fringe    = false;
  globval.Cavity_on  = false; globval.radiation      = false;
  globval.emittance  = false; globval.IBS            = false;
  globval.pathlength = false; globval.bpm            = 0;
  globval.Cart_Bend  = false; globval.dip_edge_fudge = true;

  Ring_GetTwiss(true, 0e0); printglob();

  chk_drv_terms();
}
