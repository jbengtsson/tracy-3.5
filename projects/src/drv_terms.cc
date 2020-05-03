
/* Description:

   Control of on & off-momentum dynamic aperture.
   Evaluates & minimizes the Lie generators to 2nd order in the sextupole
   strengths.
   Initially coded in Pascal and ported to OPA (for the conceptual design of
   SLS). The latter is a derivatvie of K. Wille's OPTIK, DELTA, Dortmund Univ.
   (i.e., an interactive Linear Optics Design Tool; coded in Pascal).

   References:

   [1] J. Bengtsson "The Sextupole Scheme for the Swiss Light Source (SLS):
       An Analytic Approach" SLS Note 9/97.
   [2] A. Streun OPA.

   Johan Bengtsson                                                            */


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
  const int  n_step = 5;

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
	   double &h_c, double &h_s)
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
	scl1 = scl*pow(twoJ[X_], m_x/2e0)*pow(twoJ[Y_], m_y/2e0)*pow(delta, m);
    	h_c += scl1*b3L*c; h_s += scl1*b3L*s;
      }
    }
  }
}


void sxt_2(const double scl, const double twoJ[],
	   const int i1, const int j1, const int k1, const int l1,
	   const int i2, const int j2, const int k2, const int l2,
	   double &h_c, double &h_s)
{
  // Second order generators.

  long int n1, n2;
  double   b3L1, a3L1, b3L2, a3L2, dnu, A1, A, phi;
  CellType *cp;

  const int  
    m_x[] = {i1+j1, i2+j2},
    m_y[] = {k1+l1, k2+l2},
    n_x[] = {i1-j1, i2-j2},
    n_y[] = {k1-l1, k2-l2};

  for (n1 = 0; n1 <= globval.Cell_nLoc; n1++) {
    if ((Cell[n1].Elem.Pkind == Mpole) &&
	(Cell[n1].Elem.M->Porder >= Sext)) {
      get_bnL_design_elem(Cell[n1].Fnum, 1, Sext, b3L1, a3L1);
      cp = &Cell[ind(n1-1)];
      dnu = n_x[0]*cp->Nu[X_] + n_y[0]*cp->Nu[Y_]; 
      A1 = scl*b3L1*pow(cp->Beta[X_], m_x[0]/2e0)*pow(cp->Beta[Y_], m_y[0]/2e0);
      for (n2 = 0; n2 <= globval.Cell_nLoc; n2++) {
	if ((Cell[n2].Elem.Pkind == Mpole)
	    && (Cell[n2].Elem.M->Porder >= Sext)) {
	  get_bnL_design_elem(Cell[n2].Fnum, 1, Sext, b3L2, a3L2);
	  cp = &Cell[ind(n2-1)];
	  phi = 2e0*M_PI*(dnu+n_x[1]*cp->Nu[X_]+n_y[1]*cp->Nu[Y_]);
	  A =
	    A1*b3L2*pow(cp->Beta[X_], m_x[1]/2e0)*pow(cp->Beta[Y_], m_y[1]/2e0)
	    *pow(twoJ[X_], ((m_x[0]+m_x[1])-2e0)/2e0)
	    *pow(twoJ[Y_], ((m_y[0]+m_y[1])-2e0)/2e0);
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


void K(const double nu[], const double twoJ[], double a[])
{
  // Anharmonic terms.

  long int k, n1, n2;
  double   b3L1, a3L1, b3L2, a3L2, pi_1x, pi_3x, pi_1xm2y, pi_1xp2y;
  double   s_1x, s_3x, s_1xm2y, s_1xp2y, c_1x, c_3x, c_1xm2y, c_1xp2y;
  double   dmu_x, dmu_y, A;
  CellType *cp1, *cp2;

  pi_1x = M_PI*nu[X_];
  pi_3x = 3e0*M_PI*nu[X_];
  pi_1xm2y = M_PI*(nu[X_]-2e0*nu[Y_]);
  pi_1xp2y = M_PI*(nu[X_]+2e0*nu[Y_]);

  s_1x = sin(pi_1x); s_3x = sin(pi_3x); s_1xm2y = sin(pi_1xm2y);
  s_1xp2y = sin(pi_1xp2y);

  for (k = 0; k < 3; k++)
    a[k] = 0e0;

  for (n1 = 0; n1 <= globval.Cell_nLoc; n1++) {
    if ((Cell[n1].Elem.Pkind == Mpole) && (Cell[n1].Elem.M->Porder >= Sext)) {
      get_bnL_design_elem(Cell[n1].Fnum, 1, Sext, b3L1, a3L1);
      cp1 = &Cell[ind(n1-1)];
      for (n2 = 0; n2 <= globval.Cell_nLoc; n2++) {
	if ((Cell[n2].Elem.Pkind == Mpole) &&
	    (Cell[n2].Elem.M->Porder >= Sext)) {
	  get_bnL_design_elem(Cell[n2].Fnum, 1, Sext, b3L2, a3L2);
	  cp2 = &Cell[ind(n2-1)];
	  dmu_x = 2e0*M_PI*fabs(cp1->Nu[X_]-cp2->Nu[X_]);
	  dmu_y = 2e0*M_PI*fabs(cp1->Nu[Y_]-cp2->Nu[Y_]);
	  A = b3L1*b3L2*sqrt(cp1->Beta[X_]*cp2->Beta[X_]);
	  c_1x = cos(dmu_x-pi_1x)/s_1x; c_3x = cos(3e0*dmu_x-pi_3x)/s_3x;
	  c_1xm2y = cos(dmu_x-2e0*dmu_y-pi_1xm2y)/s_1xm2y;
	  c_1xp2y = cos(dmu_x+2e0*dmu_y-pi_1xp2y)/s_1xp2y;
	  a[0] += A*cp1->Beta[X_]*cp2->Beta[X_]*(3*c_1x+c_3x)/2e0;
	  a[1] -=
	    A*cp1->Beta[Y_]
	    *(4e0*twoJ[X_]*cp2->Beta[X_]*c_1x+2e0*twoJ[Y_]*cp2->Beta[Y_]
	      *(c_1xm2y-c_1xp2y));
	  a[2] += A*cp1->Beta[Y_]*cp2->Beta[Y_]*(4e0*c_1x+c_1xm2y+c_1xp2y)/2e0;
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


void cross_terms(const double delta_eps, const double nu[], const double twoJ[],
		 double ad[])
{
  int    k;
  double a[2][3];

  Ring_GetTwiss(false, -delta_eps);
  K(nu, twoJ, a[0]);
  Ring_GetTwiss(false, delta_eps);
  K(nu, twoJ, a[1]);

  for (k = 0; k < 3; k++)
    ad[k] = (a[1][k]-a[0][k])/(2e0*delta_eps);
}


void h_1(const double scl, const double twoJ[], const double delta,
	 const int i, const int j, const int k, const int l,
	 const int m, std::vector<string> &h_label,
	 std::vector<double> &h_c, std::vector<double> &h_s, const bool ker)
{
  double       c, s;
  stringstream str;

  sxt_1(scl, twoJ, delta, i, j, k, l, m, c, s);
  if (!ker)
    str << "h_" << i << j << k << l << m;
  else
    str << "k_" << i << j << k << l << m;
  h_label.push_back(str.str()); h_c.push_back(c); h_s.push_back(s);
}


void h_2(const string &str, const double c, const double s,
	 std::vector<string> &h_label, std::vector<double> &h_c,
	 std::vector<double> &h_s)
{
  h_label.push_back(str); h_c.push_back(c); h_s.push_back(s);
}


void first_order(const double twoJ[], const double delta,
		 std::vector<string> &h_label, std::vector<double> &h_c,
		 std::vector<double> &h_s)
{
  h_1( 1e0,      twoJ, delta, 1, 1, 0, 0, 1, h_label, h_c, h_s, true);
  h_1(-1e0,      twoJ, delta, 0, 0, 1, 1, 1, h_label, h_c, h_s, true);

  h_1( 1e0/4e0,  twoJ, delta, 2, 0, 0, 0, 1, h_label, h_c, h_s, false);
  h_1(-1e0/4e0,  twoJ, delta, 0, 0, 2, 0, 1, h_label, h_c, h_s, false);
  h_1( 1e0,      twoJ, delta, 1, 0, 0, 0, 2, h_label, h_c, h_s, false);

  h_1(-1e0/12e0, twoJ, delta, 3, 0, 0, 0, 0, h_label, h_c, h_s, false);
  h_1(-1e0/4e0,  twoJ, delta, 2, 1, 0, 0, 0, h_label, h_c, h_s, false);
  h_1( 1e0/4e0,  twoJ, delta, 1, 0, 2, 0, 0, h_label, h_c, h_s, false);
  h_1( 1e0/2e0,  twoJ, delta, 1, 0, 1, 1, 0, h_label, h_c, h_s, false);
  h_1( 1e0/4e0,  twoJ, delta, 1, 0, 0, 2, 0, h_label, h_c, h_s, false);
}


void second_order(const double twoJ[], std::vector<string> &h_label,
		  std::vector<double> &h_c, std::vector<double> &h_s)
{
  double c, s;

  c = s = 0e0;
  sxt_2(-1e0/32e0, twoJ, 2, 1, 0, 0, 3, 0, 0, 0, c, s);
  h_2("h_40000", c, s, h_label, h_c, h_s);
  
  c = s = 0e0;
  sxt_2(-1e0/16e0, twoJ, 1, 2, 0, 0, 3, 0, 0, 0, c, s);
  h_2("h_31000", c, s, h_label, h_c, h_s);

  c = s = 0e0;
  sxt_2(-1e0/32e0, twoJ, 0, 3, 0, 0, 3, 0, 0, 0, c, s);
  sxt_2(-3e0/32e0, twoJ, 1, 2, 0, 0, 2, 1, 0, 0, c, s);
  h_2("h_22000", c, s, h_label, h_c, h_s);

  c = s = 0e0;
  sxt_2(-1e0/32e0, twoJ, 3, 0, 0, 0, 0, 1, 2, 0, c, s);
  sxt_2(-1e0/32e0, twoJ, 1, 0, 2, 0, 2, 1, 0, 0, c, s);
  sxt_2(-4e0/32e0, twoJ, 1, 0, 1, 1, 1, 0, 2, 0, c, s);
  h_2("h_20200", c, s, h_label, h_c, h_s);
  
  c = s = 0e0;
  sxt_2(-1e0/16e0, twoJ, 1, 0, 2, 0, 1, 2, 0, 0, c, s);
  sxt_2(-1e0/16e0, twoJ, 2, 1, 0, 0, 0, 1, 2, 0, c, s);
  sxt_2(-2e0/16e0, twoJ, 0, 1, 1, 1, 1, 0, 2, 0, c, s);
  sxt_2(-2e0/16e0, twoJ, 1, 0, 1, 1, 0, 1, 2, 0, c, s);
  h_2("h_11200", c, s, h_label, h_c, h_s);
    
  c = s = 0e0;
  sxt_2(-1e0/16e0, twoJ, 3, 0, 0, 0, 0, 1, 1, 1, c, s);
  sxt_2(-1e0/16e0, twoJ, 1, 0, 1, 1, 2, 1, 0, 0, c, s);
  sxt_2(-1e0/8e0, twoJ, 1, 0, 0, 2, 1, 0, 2, 0, c, s);
  h_2("h_20110", c, s, h_label, h_c, h_s);
  
  c = s = 0e0;
  sxt_2(-1e0/8e0, twoJ, 1, 0, 1, 1, 1, 2, 0, 0, c, s);
  sxt_2(-1e0/8e0, twoJ, 2, 1, 0, 0, 0, 1, 1, 1, c, s);
  sxt_2(-1e0/8e0, twoJ, 0, 1, 0, 2, 1, 0, 2, 0, c, s);
  sxt_2(-1e0/8e0, twoJ, 1, 0, 0, 2, 0, 1, 2, 0, c, s);
  h_2("h_11110", c, s, h_label, h_c, h_s);
  
  c = s = 0e0;
  sxt_2(-1e0/16e0, twoJ, 0, 1, 2, 0, 1, 0, 1, 1, c, s);
  sxt_2(-1e0/16e0, twoJ, 0, 1, 1, 1, 1, 0, 2, 0, c, s);
  h_2("h_00310", c, s, h_label, h_c, h_s);
  
  c = s = 0e0;
  sxt_2(-1e0/32e0, twoJ, 1, 0, 0, 2, 2, 1, 0, 0, c, s);
  sxt_2(-1e0/32e0, twoJ, 3, 0, 0, 0, 0, 1, 0, 2, c, s);
  sxt_2(-4e0/32e0, twoJ, 1, 0, 0, 2, 1, 0, 1, 1, c, s);
  h_2("h_20020", c, s, h_label, h_c, h_s);

  c = s = 0e0;
  sxt_2(-1e0/32e0, twoJ, 0, 1, 2, 0, 1, 0, 2, 0, c, s);
  h_2("h_00400", c, s, h_label, h_c, h_s);

  c = s = 0e0;
  sxt_2(-1e0/8e0, twoJ, 0, 1, 1, 1, 1, 0, 1, 1, c, s);
  sxt_2(-1e0/32e0, twoJ, 0, 1, 0, 2, 1, 0, 2, 0, c, s);
  sxt_2(-1e0/32e0, twoJ, 0, 1, 2, 0, 1, 0, 0, 2, c, s);
  h_2("h_00220", c, s, h_label, h_c, h_s);
}


void anharm(const double twoJ[], std::vector<string> &h_label,
	    std::vector<double> &h_c, std::vector<double> &h_s)
{
  double a[3];

  K(globval.TotalTune, twoJ, a);
  h_label.push_back("k_22000"); h_c.push_back(a[0]); h_s.push_back(0e0);
  h_label.push_back("k_11110"); h_c.push_back(a[1]); h_s.push_back(0e0);
  h_label.push_back("k_00220"); h_c.push_back(a[2]); h_s.push_back(0e0);
}


void ksi_2(const double delta_eps, const double twoJ[],
	   std::vector<string> &h_label, std::vector<double> &h_c,
	   std::vector<double> &h_s)
{
  double ksi2[2];

  get_ksi2_(delta_eps, ksi2);
  h_label.push_back("k_11002"); h_c.push_back(ksi2[X_]); h_s.push_back(0e0);
  h_label.push_back("k_00112"); h_c.push_back(ksi2[Y_]); h_s.push_back(0e0);
}


void cross_terms(const double delta_eps, const double twoJ[],
		 std::vector<string> &h_label, std::vector<double> &h_c,
		 std::vector<double> &h_s)
{
  double a[3];

  cross_terms(delta_eps, globval.TotalTune, twoJ, a);
  h_label.push_back("k_22001"); h_c.push_back(a[0]); h_s.push_back(0e0);
  h_label.push_back("k_11111"); h_c.push_back(a[1]); h_s.push_back(0e0);
  h_label.push_back("k_00221"); h_c.push_back(a[2]); h_s.push_back(0e0);
}


void prt_drv_terms(std::vector<string> &h_label, std::vector<double> &h_c,
		   std::vector<double> &h_s)
{
  int k;

  printf("\nFirst order chromatic terms:\n");
  for (k = 0; k < 2; k++)
    printf("  %7s = [%12.5e, %12.5e]\n",
	   h_label[k].c_str(), h_c[k]/2e0, h_s[k]/2e0);
  printf("\n");
  for (k = 2; k < 5; k++)
    printf("  %7s = [%12.5e, %12.5e]\n",
	   h_label[k].c_str(), h_c[k]/2e0, h_s[k]/2e0);
  printf("\nFirst order geometric terms:\n");
  for (k = 5; k < 10; k++)
    printf("  %7s = [%12.5e, %12.5e]\n",
	   h_label[k].c_str(), h_c[k]/2e0, h_s[k]/2e0);
  printf("\nSecond order geometric terms:\n");
  for (k = 10; k < 21; k++)
    printf("  %7s = [%12.5e, %12.5e]\n",
	   h_label[k].c_str(), h_c[k]/2e0, h_s[k]/2e0);
  printf("\nSecond order anharmonic terms:\n");
  for (k = 21; k < 24; k++)
    printf("  %7s = [%12.5e, %12.5e]\n", h_label[k].c_str(), h_c[k], h_s[k]);
  printf("\nSecond order chromaticity:\n");
  for (k = 24; k < 26; k++)
    printf("  %7s = [%12.5e, %12.5e]\n",
	   h_label[k].c_str(), -h_c[k]/M_PI, -h_s[k]/M_PI);
  printf("\nSecond order cross terms:\n");
  for (k = 26; k < 29; k++)
    printf("  %7s = [%12.5e, %12.5e]\n", h_label[k].c_str(), h_c[k], h_s[k]);
}


void drv_terms(const double twoJ[], const double delta, const double delta_eps,
	       double k_J_delta[])
{
  int                 k;
  std::vector<string> h_label;
  std::vector<double> h_c, h_s;

  const bool prt = false;

  first_order(twoJ, delta, h_label, h_c, h_s);
  second_order(twoJ, h_label, h_c, h_s);
  anharm(twoJ, h_label, h_c, h_s);
  ksi_2(delta_eps, twoJ, h_label, h_c, h_s);
  cross_terms(delta_eps, twoJ, h_label, h_c, h_s);

  if (prt) prt_drv_terms(h_label, h_c, h_s);

  for (k = 0; k < 3; k++)
    k_J_delta[k] = h_c[26+k];
}
