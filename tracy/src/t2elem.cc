/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -

   Element propagators.                                                      */

bool          first_FM = true;
double        c_1, d_1, c_2, d_2, cl_rad, q_fluct;
double        I2, I4, I5, dcurly_H, dI4, s_FM;
std::ofstream outf_;


template<typename T>
void spline_(const double x[], const T y[], const int n, const double yp1,
	     const double ypn, T y2[])
{
  int    i,k;
  double sig;
  // Variable length arrays for user defined data types is not supported by
  // ANSI C++.
  // T      p, u[n], qn, un;
  T      p, qn, un;
  T      *u = new T[n];

  if (yp1 > 0.99e30)
    y2[1] = u[1] = 0e0;
  else {
    y2[1] = -0.5;
    u[1] = (3e0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
  }
  for (i = 2; i <= n-1; i++) {
    sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p = sig*y2[i-1]+2e0;
    y2[i] = (sig-1e0)/p;
    u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i] = (6e0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn = un = 0e0;
  else {
    qn = 0.5;
    un = (3e0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
  }
  y2[n] = (un-qn*u[n-1])/(qn*y2[n-1]+1e0);
  for (k = n-1; k >= 1; k--)
    y2[k] = y2[k]*y2[k+1]+u[k];

  delete [] u;
}


void splie2_(const double x1a[], const double x2a[], double **ya,
	     const int m, const int n, double **y2a)
{
  int j;

  for (j = 1; j <= m; j++)
    spline_(x2a,ya[j],n,1e30,1e30,y2a[j]);
}


template<typename T, typename U>
void splint_(const double xa[], const U ya[], const U y2a[],
	     const int n, const T &x, T &y)
{
  int    klo,khi,k;
  double h;
  T      a,b;

  klo = 1;
  khi = n;
  while (khi-klo > 1) {
    k = (khi+klo) >> 1;
    if (xa[k] > x) khi = k;
    else klo = k;
  }
  h = xa[khi]-xa[klo];
  if (h == 0e0) nrerror("Bad xa input to routine splint_");
  a = (xa[khi]-x)/h;
  b = (x-xa[klo])/h;
  y = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6e0;
}


template<typename T>
void splin2_(const double x1a[], const double x2a[], double **ya, double **y2a,
	     const int m, const int n, const T &x1, const T &x2, T &y)
{
  int j;
  // Variable length arrays for user defined data types is not supported by
  // ANSI C++.
  // T   ytmp[m+1], yytmp[m+1];
  T   *ytmp = new T[m+1], *yytmp = new T[m+1];

  if ((x1 < x1a[1]) || (x1 > x1a[m])) {
    std::cout << std::fixed << std::setprecision(8)
	 << "splin2_: x undefined ["
	 << is_double<T>::cst(x1) << ", " << is_double<T>::cst(x2) << "] (["
	 << x1a[1] << ", " << x1a[m] << "])" << std::endl;

    y = NAN;

    return;
  }

  if ((x2 < x2a[1]) || (x2 > x2a[n])) {
    std::cout << std::fixed << std::setprecision(8)
	 << "splin2_: y undefined ["
	 << is_double<T>::cst(x1) << ", " << is_double<T>::cst(x2) << "] (["
	 << x2a[1] << ", " << x2a[n] << "])" << std::endl;

    y = NAN;

    return;
  }

  for (j = 1; j<= m; j++)
    splint_(x2a,ya[j],y2a[j],n,x2,yytmp[j]);
  spline_(x1a,yytmp,m,1.0e30,1.0e30,ytmp);
  splint_(x1a,yytmp,ytmp,m,x1,y);

  delete [] ytmp; delete [] yytmp;
}


template<typename T>
void CellType::GtoL(ss_vect<T> &ps,
		    const double c0, const double c1, const double s1)
{
  ss_vect<T> ps1;

  // Simplified rotated p_rot.
  ps[px_] += c1; ps[py_] += s1;
  // Translate.
  ps[x_] -= dS[X_]; ps[y_] -= dS[Y_];
  // Rotate.
  ps1 = ps;
  ps[x_]  =  dR[X_]*ps1[x_]  + dR[Y_]*ps1[y_];
  ps[px_] =  dR[X_]*ps1[px_] + dR[Y_]*ps1[py_];
  ps[y_]  = -dR[Y_]*ps1[x_]  + dR[X_]*ps1[y_];
  ps[py_] = -dR[Y_]*ps1[px_] + dR[X_]*ps1[py_] ;
  // Simplified p_rot.
  ps[px_] -= c0;
}


template<typename T>
void CellType::LtoG(ss_vect<T> &ps,
		    const double c0, const double c1, const double s1)
{
  ss_vect<T> ps1;

  // Simplified p_rot.
  ps[px_] -= c0;
  // Rotate.
  ps1 = ps;
  ps[x_]  = dR[X_]*ps1[x_]  - dR[Y_]*ps1[y_];
  ps[px_] = dR[X_]*ps1[px_] - dR[Y_]*ps1[py_];
  ps[y_]  = dR[Y_]*ps1[x_]  + dR[X_]*ps1[y_];
  ps[py_] = dR[Y_]*ps1[px_] + dR[X_]*ps1[py_];
  // Translate.
  ps[x_] += dS[X_]; ps[y_] += dS[Y_];
  // p_rot rotated.
  ps[px_] += c1; ps[py_] += s1;
}


template<typename T>
inline T get_p_s(const ss_vect<T> &x)
{
  T p_s, p_s2;

  if (!Lattice.param.H_exact)
    // Small angle axproximation.
    p_s = 1e0 + x[delta_];
  else {
    p_s2 = sqr(1e0+x[delta_]) - sqr(x[px_]) - sqr(x[py_]);
    if (p_s2 >= 0e0)
      p_s = sqrt(p_s2);
    else {
//      printf("get_p_s: *** Speed of light exceeded!\n");
      p_s = NAN;
    }
  }
  return(p_s);
}


template<typename T>
void Drift(const double L, ss_vect<T> &ps)
{
  T u;

  if (!Lattice.param.H_exact) {
    // Small angle axproximation.
    u = L/(1e0+ps[delta_]);
    ps[x_]  += u*ps[px_]; ps[y_] += u*ps[py_];
    ps[ct_] += u*(sqr(ps[px_])+sqr(ps[py_]))/(2e0*(1e0+ps[delta_]));
  } else {
    u = L/get_p_s(ps);
    ps[x_]  += u*ps[px_]; ps[y_] += u*ps[py_];
    ps[ct_] += u*(1e0+ps[delta_]) - L;
  }
  if (Lattice.param.pathlength) ps[ct_] += L;
}


template<typename T>
void DriftType::Elem_Propagate(ss_vect<T> &x) { Drift(this->L, x); }


// partial template-class specialization
// primary version
template<typename T>
class is_tps { };

// partial specialization
template<>
class is_tps<double> {
 public:
  static inline void get_ps(const ss_vect<double> &x, CellType &Cell)
  { Cell.BeamPos = x; }

  static inline double set_prm(const int k) { return 1e0; }

  static inline double get_curly_H(const ss_vect<tps> &x)
  {
    std::cout << "get_curly_H: operation not defined for double" << std::endl;
    exit_(1);
    return 0e0;
  }

  static inline double get_dI4(const double h, const double b2, const double L,
			     const ss_vect<tps> &x)
  {
    std::cout << "get_dI4: operation not defined for double" << std::endl;
    exit_(1);
    return 0e0;
  }

  static inline void emittance(const double B2, const double u,
			       const double ps0, const ss_vect<double> &xp) { }

  static inline void diff_mat(const double B2, const double u,
			      const double ps0, const ss_vect<double> &xp) { }
};


// partial specialization
template<>
class is_tps<tps> {
 public:
  static inline void get_ps(const ss_vect<tps> &x, CellType &Cell)
  {
    Cell.BeamPos = x.cst(); getlinmat(6, x, Cell.A);
  }

  static inline tps set_prm(const int k) { return tps(0e0, k); }

  static inline double get_curly_H(const ss_vect<tps> &A)
  {
    int             j;
    double          curly_H[2];
    ss_vect<double> eta;

    eta.zero();
    for (j = 0; j < 4; j++)
      eta[j] = A[j][delta_];

    get_twoJ(2, eta, A, curly_H);

    return curly_H[X_];
  }

  static inline double get_dI4(const ss_vect<tps> &A) { return A[x_][delta_]; }

  static inline void emittance(const tps &B2_perp, const tps &ds,
			       const tps &ps0, const ss_vect<tps> &A)
  {
    // M. Sands "The hysics of Electron Storage Rings" SLAC-121, p. 118.
    // d<delta^2>/ds = 3*C_U*C_gamma*h_bar*c*E_0^5*(1+delta)^4*(B_perp/(Brho))^3
    //                 /(4*pi*m_e^3)
    // A contains the eigenvectors.
    int          j;
    double       B_66;
    ss_vect<tps> A_inv;

    if (B2_perp > 0e0) {
      B_66 = (q_fluct*pow(B2_perp.cst(), 1.5)*pow(ps0, 4)*ds).cst();
      A_inv = Inv(A);
      // D_11 = D_22 = curly_H_x,y * B_66 / 2,
      // curly_H_x,y = eta_Fl^2 + etap_Fl^2
      for (j = 0; j < 3; j++)
	Lattice.param.D_rad[j] +=
	  (sqr(A_inv[j*2][delta_])+sqr(A_inv[j*2+1][delta_]))*B_66/2e0;
    }
  }

  static inline void diff_mat(const tps &B2_perp, const tps &ds, const tps &ps0,
			      ss_vect<tps> &x) { }
};


template<typename T>
void get_B2(const double h_ref, const T B[], const ss_vect<T> &xp,
	    T &B2_perp, T &B2_par)
{
  // compute B_perp^2 and B_par^2
  T xn, e[3];

  xn = 1e0/sqrt(sqr(1e0+xp[x_]*h_ref)+sqr(xp[px_])+sqr(xp[py_]));
  e[X_] = xp[px_]*xn; e[Y_] = xp[py_]*xn; e[Z_] = (1e0+xp[x_]*h_ref)*xn;

  // left-handed coordinate system
  B2_perp =
    sqr(B[Y_]*e[Z_]-B[Z_]*e[Y_]) + sqr(B[X_]*e[Y_]-B[Y_]*e[X_])
    + sqr(B[Z_]*e[X_]-B[X_]*e[Z_]);

//  B2_par = sqr(B[X_]*e[X_]+B[Y_]*e[Y_]+B[Z_]*e[Z_]);
}


template<typename T>
void radiate(ss_vect<T> &x, const double L, const double h_ref, const T B[])
{
  // M. Sands "The hysics of Electron Storage Rings" SLAC-121, p. 98.
  // ddelta/d(ds) = -C_gamma*E_0^3*(1+delta)^2*(B_perp/(Brho))^2/(2*pi)
  T          ps0, ps1, ds, B2_perp = 0e0, B2_par = 0e0;
  ss_vect<T> xp;

  // large ring: x' and y' unchanged
  xp = x; ps0 = get_p_s(x); xp[px_] /= ps0; xp[py_] /= ps0;

  // H = -p_s => ds = H*L
  ds = (1e0+xp[x_]*h_ref+(sqr(xp[px_])+sqr(xp[py_]))/2e0)*L;
  get_B2(h_ref, B, xp, B2_perp, B2_par);

  if (Lattice.param.radiation) {
    x[delta_] -= cl_rad*sqr(ps0)*B2_perp*ds;
    ps1 = get_p_s(x); x[px_] = xp[px_]*ps1; x[py_] = xp[py_]*ps1;
  }

  if (Lattice.param.emittance) is_tps<T>::emittance(B2_perp, ds, ps0, xp);
}


template<typename T>
void radiate_ID(ss_vect<T> &x, const double L, const T &B2_perp)
{
  T          ps0, ps1, ds;
  ss_vect<T> xp;

  // large ring: x' and y' unchanged
  xp = x; ps0 = get_p_s(x); xp[px_] /= ps0; xp[py_] /= ps0;

  // H = -p_s => ds = H*L
  ds = (1e0+(sqr(xp[px_])+sqr(xp[py_]))/2e0)*L;

  if (Lattice.param.radiation) {
    x[delta_] -= cl_rad*sqr(ps0)*B2_perp*ds;
    ps1 = get_p_s(x); x[px_] = xp[px_]*ps1; x[py_] = xp[py_]*ps1;
  }

  if (Lattice.param.emittance) is_tps<T>::emittance(B2_perp, ds, ps0, xp);
}


static double get_psi(double irho, double phi, double gap)
{
  /* Correction for magnet gap (longitudinal fringe field)

       irho h = 1/rho [1/m]
       phi  edge angle
       gap  full gap between poles

                                    2
                   K1*gap*h*(1 + sin phi)
            psi = ----------------------- * (1 - K2*g*gap*tan phi)
                        cos phi

            K1 is usually 1/2
            K2 is zero here                                                  */

  double psi;

  const double k1 = 0.5e0, k2 = 0e0;

  if (phi == 0e0)
    psi = 0e0;
  else
    psi = k1*gap*irho*(1e0+sqr(sin(dtor(phi))))/cos(dtor(phi))
          *(1e0 - k2*gap*irho*tan(dtor(phi)));

  return psi;
}


template<typename T>
void thin_kick(const int Order, const double MB[], const double L,
	       const double h_bend, const double h_ref, ss_vect<T> &x)
{
  // The vector potential for the combined-function sector bend is from:
  // C. Iselin "Lie Transformations and Transport Equations for Combined-
  // Function Dipoles" Part. Accel. 17, 143-155 (1985).
  int        j;
  T          BxoBrho, ByoBrho, ByoBrho1, B[3], u, p_s;
  ss_vect<T> x0;

  if ((h_bend != 0e0) || ((1 <= Order) && (Order <= HOMmax))) {
    x0 = x;
    // Compute magnetic field with Horner's rule.
    ByoBrho = MB[Order+HOMmax]; BxoBrho = MB[HOMmax-Order];
    for (j = Order-1; j >= 1; j--) {
      ByoBrho1 = x0[x_]*ByoBrho - x0[y_]*BxoBrho + MB[j+HOMmax];
      BxoBrho  = x0[y_]*ByoBrho + x0[x_]*BxoBrho + MB[HOMmax-j];
      ByoBrho  = ByoBrho1;
    }

    if (Lattice.param.radiation || Lattice.param.emittance) {
      B[X_] = BxoBrho; B[Y_] = ByoBrho + h_bend; B[Z_] = 0e0;
      radiate(x, L, h_ref, B);
    }

    if (h_ref != 0e0) {
      // Sector bend.
      if (true) {
	x[px_] -= L*(ByoBrho+(h_bend-h_ref)/2e0+h_ref*h_bend*x0[x_]
		     -h_ref*x0[delta_]);
	x[ct_] += L*h_ref*x0[x_];
      } else {
	// The Hamiltonian is split into: H_d + H_k; with [H_d, H_d] = 0.
	p_s = get_p_s(x0); u = L*h_ref*x0[x_]/p_s;
	x[x_]  += u*x0[px_]; x[y_]  += u*x0[py_]; x[ct_] += u*(1e0+x0[delta_]);
	// x[px_] -= L*(h_bend*(1e0+h_ref*x0[x_])-h_ref*p_s);
	// Field expansion up to sextupole like terms.
	ByoBrho += h_bend - MB[Quad+HOMmax]*h_ref*sqr(x0[y_])/2e0;
	x[px_] -= L*((1e0+h_ref*x0[x_])*ByoBrho-h_ref*p_s);
	x[py_] += L*(1e0+h_ref*x0[x_])*BxoBrho;
      }
    } else
      // Cartesian bend.
      x[px_] -= L*(h_bend+ByoBrho);
    x[py_] += L*BxoBrho;
  }
}


template<typename T>
void EdgeFocus(const double irho, const double phi, const double gap,
	       ss_vect<T> &x)
{
  x[px_] += irho*tan(dtor(phi))*x[x_];
  if (!Lattice.param.dip_edge_fudge) {
    // warning: => diverging Taylor map (see SSC-141)
    // x[py_] -=
    //   irho*tan(dtor(phi)-get_psi(irho, phi, gap))*x[y_]/(1e0+x[delta_]);
    // Leading order correction.
    x[py_] -=
      irho*tan(dtor(phi)-get_psi(irho, phi, gap))*x[y_]*(1e0-x[delta_]);
  } else
    x[py_] -= irho*tan(dtor(phi)-get_psi(irho, phi, gap))*x[y_];
}


template<typename T>
void p_rot(double phi, ss_vect<T> &ps)
{
  T          c, s, t, pz, p, val;
  ss_vect<T> ps1;

  c = cos(dtor(phi)); s = sin(dtor(phi)); t = tan(dtor(phi)); pz = get_p_s(ps);

  if (!Lattice.param.H_exact) {
     ps[px_] = s*pz + c*ps[px_];
  } else {
    // ps1 = ps; p = c*pz - s*ps1[px_];
    // px[x_]   = ps1[x_]*pz/p; px[px_] = s*pz + c*ps1[px_];
    // px[y_]  += ps1[x_]*ps1[py_]*s/p;
    // px[ct_] += (1e0+ps1[delta_])*ps1[x_]*s/p;

    ps1 = ps; val = 1e0 - ps1[px_]*t/pz;
    ps[x_]  = ps1[x_]/(c*val);
    ps[px_] = ps1[px_]*c + s*pz;
    ps[y_]  = ps1[y_] + t*ps1[x_]*ps1[py_]/(pz*val);
    ps[ct_] = ps1[ct_] + ps1[x_]*(1e0+ps1[delta_])*t/(pz*val);
  }
}


template<typename T>
void bend_fringe(const double hb, ss_vect<T> &ps)
{
  T          coeff, u, pz, pz2, pz3;
  ss_vect<T> ps1;

  coeff = -hb/2e0; ps1 = ps; pz = get_p_s(ps); pz2 = sqr(pz); pz3 = pz*pz2;
  u = 1e0 + 4e0*coeff*ps1[px_]*ps1[y_]*ps1[py_]/pz3;
  if (u >= 0e0) {
    ps[y_]  = 2e0*ps1[y_]/(1e0+sqrt(u));
    ps[x_]  = ps1[x_] - coeff*sqr(ps[y_])*(pz2+sqr(ps1[px_]))/pz3;
    ps[py_] = ps1[py_] + 2e0*coeff*ps1[px_]*ps[y_]/pz;
    ps[ct_] = ps1[ct_] - coeff*ps1[px_]*sqr(ps[y_])*(1e0+ps1[delta_])/pz3;
  } else {
    printf("bend_fringe: *** Speed of light exceeded!\n");
    ps[x_] = NAN; ps[px_] = NAN; ps[y_] = NAN; ps[py_] = NAN;
    ps[delta_] = NAN; ps[ct_] = NAN;
  }
}


template<typename T>
void quad_fringe(const double b2, ss_vect<T> &x)
{
  T u, ps;

  u = b2/(12e0*(1e0+x[delta_])); ps = u/(1e0+x[delta_]);
  x[py_] /= 1e0 - 3e0*u*sqr(x[y_]); x[y_] -= u*cube(x[y_]);
  if (Lattice.param.Cavity_on) x[ct_] -= ps*cube(x[y_])*x[py_];
  x[px_] /= 1e0 + 3e0*u*sqr(x[x_]);
  if (Lattice.param.Cavity_on) x[ct_] += ps*cube(x[x_])*x[px_];
  x[x_] += u*cube(x[x_]); u = u*3e0; ps = ps*3e0;
  x[y_] = exp(-u*sqr(x[x_]))*x[y_]; x[py_] = exp(u*sqr(x[x_]))*x[py_];
  x[px_] += 2e0*u*x[x_]*x[y_]*x[py_];
  if (Lattice.param.Cavity_on) x[ct_] -= ps*sqr(x[x_])*x[y_]*x[py_];
  x[x_] = exp(u*sqr(x[y_]))*x[x_]; x[px_] = exp(-u*sqr(x[y_]))*x[px_];
  x[py_] -= 2e0*u*x[y_]*x[x_]*x[px_];
  if (Lattice.param.Cavity_on) x[ct_] += ps*sqr(x[y_])*x[x_]*x[px_];
}


template<typename T>
void MpoleType::Elem_Propagate(ss_vect<T> &x)
{
  int    seg = 0;
  double k = 0e0, dL = 0e0, dL1 = 0e0, dL2 = 0e0;
  double dkL1 = 0e0, dkL2 = 0e0, h_ref = 0e0;

  // Global -> Local.
  // Can not be moved to ElemType level because it needs: c0, c1, and s1.
  this->GtoL(x, this->c0, this->c1, this->s1);

  if ((this->method == Meth_Second) || (this->method == Meth_Fourth)) {
    // Fringe fields.
    if (Lattice.param.quad_fringe && (this->B[Quad+HOMmax] != 0e0))
      quad_fringe(this->B[Quad+HOMmax], x);
    if (!Lattice.param.H_exact) {
      if (this->irho != 0e0) EdgeFocus(this->irho, this->Tx1, this->gap, x);
    } else {
      p_rot(this->Tx1, x);
      if (Lattice.param.dip_fringe) bend_fringe(this->irho, x);
    }
  }

  if (this->thick == thicktype(thick_)) {
    if (!Lattice.param.H_exact) {
      // Polar coordinates.
      h_ref = this->irho; dL = this->L/this->N;
    } else {
      // Cartesian coordinates.
      h_ref = 0e0;
      if (this->irho == 0e0)
	dL = this->L/this->N;
      else
	dL = 2e0/this->irho*sin(this->L*this->irho/2e0)/this->N;
    }
  }

  switch (this->method) {

  case Meth_Linear:

  case Meth_First:
    if (this->thick == thicktype(thick_)) {
      /* First Linear  */
      //      LinTrans(5L, this->AU55, x);
      k = this->B[Quad+HOMmax];
      /* retrieve normal quad component already in AU55 */
      this->B[Quad+HOMmax] = 0e0;
      /* Kick w/o quad component */
      thin_kick(this->order, this->B, this->L, 0e0, 0e0, x);
      /* restore quad component */
      this->B[Quad+HOMmax] = k;
      /* Second Linear */
      //      LinTrans(5L, this->AD55, x);
    } else /* thin kick */
      thin_kick(this->order, this->B, 1e0, 0e0, 0e0, x);
    break;

  case Meth_Second:
    std::cout << "Mpole_Pass: Meth_Second not supported" << std::endl;
    exit_(0);
    break;

  case Meth_Fourth:
    if (this->thick == thicktype(thick_)) {
      dL1 = c_1*dL; dL2 = c_2*dL; dkL1 = d_1*dL; dkL2 = d_2*dL;

      dcurly_H = 0e0; dI4 = 0e0;
      for (seg = 1; seg <= this->N; seg++) {
	if (Lattice.param.emittance && (!Lattice.param.Cavity_on) &&
	    (this->irho != 0e0)) {
	  dcurly_H += is_tps<tps>::get_curly_H(x);
	  dI4 += is_tps<tps>::get_dI4(x);
	}

	Drift(dL1, x);
        thin_kick(this->order, this->B, dkL1, this->irho, h_ref, x);
	Drift(dL2, x);
        thin_kick(this->order, this->B, dkL2, this->irho, h_ref, x);

	if (Lattice.param.emittance && (!Lattice.param.Cavity_on) &&
	    (this->irho != 0e0)) {
	  dcurly_H += 4e0*is_tps<tps>::get_curly_H(x);
	  dI4 += 4e0*is_tps<tps>::get_dI4(x);
	}

	Drift(dL2, x);
        thin_kick(this->order, this->B, dkL1, this->irho, h_ref, x);
	Drift(dL1, x);

	if (Lattice.param.emittance && (!Lattice.param.Cavity_on) &&
	    (this->irho != 0e0)) {
	  dcurly_H += is_tps<tps>::get_curly_H(x);
	  dI4 += is_tps<tps>::get_dI4(x);
	}
      }

      if (Lattice.param.emittance && (!Lattice.param.Cavity_on) &&
	  (this->irho != 0)) {
	dcurly_H /= 6e0*this->N;
	dI4 *=
	  this->irho*(sqr(this->irho)+2e0*this->Bpar[Quad+HOMmax])
	  /(6e0*this->N);
	I2 += this->L*sqr(this->irho); I4 += this->L*dI4;
	I5 += this->L*fabs(cube(this->irho))*dcurly_H;
      }
    } else
      thin_kick(this->order, this->B, 1e0, 0e0, 0e0, x);

    break;
  }

  if ((this->method == Meth_Second) || (this->method == Meth_Fourth)) {
    // Fringe fields.
    if (!Lattice.param.H_exact) {
      if (this->irho != 0e0) EdgeFocus(this->irho, this->Tx2, this->gap, x);
    } else {
      if (Lattice.param.dip_fringe) bend_fringe(-this->irho, x);
      p_rot(this->Tx2, x);
    }
    if (Lattice.param.quad_fringe && (this->B[Quad+HOMmax] != 0e0))
      quad_fringe(-this->B[Quad+HOMmax], x);
  }

  // Local -> Global.
  this->LtoG(x, this->c0, this->c1, this->s1);
}


template<typename T>
void MarkerType::Elem_Propagate(ss_vect<T> &X)
{
  /* Global -> Local */
  // this->GtoL(X, 0e0, 0e0, 0e0);
  /* Local -> Global */
  // this->LtoG(X, 0e0, 0e0, 0e0);
}


#if 1

template<typename T>
void CavityType::Elem_Propagate(ss_vect<T> &X)
{
  T delta;

  if (Lattice.param.Cavity_on && this->volt != 0e0) {
    delta = -this->volt/(Lattice.param.Energy*1e9)
            *sin(2e0*M_PI*this->freq/c0*X[ct_]+this->phi);
    X[delta_] += delta;

    if (Lattice.param.radiation) Lattice.param.dE -= is_double<T>::cst(delta);

    if (Lattice.param.pathlength) X[ct_] -= this->h/this->freq*c0;
  }
}

#else

template<typename T, typename U>
void Cav_Focus(const double L, const T delta, const bool entrance,
	       ss_vect<U> &ps)
{
  double sgn;
 
  sgn = (entrance)? -1e0 : 1e0;

  ps[px_] += sgn*ps[x_]*delta/(2e0*L);
  ps[py_] += sgn*ps[y_]*delta/(2e0*L);
}


template<typename T>
void Cav_Pass1(CellType &Cell, ss_vect<T> &ps)
{
  /* J. Rosenzweig and L. Serafini "Transverse Particle Motion in
     Radio-Frequency Linear Accelerators" hys. Rev. E 49(2),
     1599-1602 (1994).                                                        */

  int    k;
  double L, h, p_t1;
  T      delta_max, ddelta, delta;

  L = this->L;

  h = L/(this->N+1e0);
  // Lattice.param.Energy contains p_0.
  delta_max = this->volt/(1e9*Lattice.param.Energy); ddelta = delta_max/this->N;
  delta = delta_max*sin(2e0*M_PI*this->freq*ps[ct_]/c0+this->phi);
  if (this->entry_focus) Cav_Focus(L, delta, true, ps);
  for (k = 0; k < this->N; k++) {
    Drift(h, ps);

    ps[delta_] -= ddelta*sin(2e0*M_PI*this->freq*(ps[ct_]-k*h)/c0+this->phi);

    if (Lattice.param.radiation) Lattice.param.dE -= is_double<T>::cst(ddelta);
    if (Lattice.param.pathlength) ps[ct_] -= this->h/this->freq*c0;
  }
  Drift(h, ps);
  if (this->exit_focus) Cav_Focus(L, delta, false, ps);

  if (false) {
    // Update p_0.
    p_t1 = is_double<T>::cst(ps[delta_]);
    ps[delta_] -= p_t1;
    // Lattice.param.Energy contains p_0.
    Lattice.param.Energy *= sqrt(1e0+2e0*p_t1/Lattice.param.beta0+sqr(p_t1));
    Lattice.param.gamma0 = sqrt(sqr(m_e)+sqr(1e9*Lattice.param.Energy))/m_e;
    Lattice.param.beta0  = sqrt(1e0-1e0/sqr(Lattice.param.gamma0));
    printf("\np0 = %12.5e, beta0 = %12.5e, gamma0 = %12.5e\n",
	   Lattice.param.Energy, Lattice.param.beta0, Lattice.param.gamma0);
  }
}


template<typename T>
void Cav_Pass(ss_vect<T> &ps)
{
  /* J. Rosenzweig and L. Serafini "Transverse Particle Motion in
     Radio-Frequency Linear Accelerators" hys. Rev. E 49(2),
     1599-1602 (1994).                                                        */

  double     L, lambda, phi;
  double     dgammaMax, dgamma, gamma, gamma1;
  double     sf1, f2, f2s;
  double     f5, sf5, dpr, dct, p_t1;
  double     p0, p_t, delta, alpha, dp;
  ss_vect<T> ps0;

  const bool RandS = false;
 
  L = this->L; phi = this->phi; lambda = c0/this->freq;

  p_t = is_double<T>::cst(ps[delta_]);
  delta = sqrt(1e0+2e0*p_t/Lattice.param.beta0+sqr(p_t)) - 1e0;
  // Lattice.param.Energy contains p_0 [GeV].
  p0 = 1e9*Lattice.param.Energy/m_e;
  gamma = sqrt(1e0+sqr(p0));
  dgammaMax = this->volt/m_e; dgamma = dgammaMax*sin(phi);
  gamma1 = gamma + dgamma;
  dp = sqrt(sqr(gamma1)-1e0) - p0;

  printf("delta = %e, p0 = %e, gamma = %e, gamma1 = %e, dp = %e\n",
	 delta, p0, gamma, gamma1, dp);

  if (this->entry_focus) Cav_Focus(L, dgamma/gamma, true, ps);

  if (!RandS) {
    sf1 = sqrt(1e0+sqr(p0));
    f2 = sf1 + dgammaMax*sin(phi); f2s = sqr(f2);
    f5 = f2s - 1e0; sf5 = sqrt(f5);

    printf("p0 = %e, dgammaMax = %e\n", p0, dgammaMax);
    printf("f2= %e, f5 = %e\n", f2, f5);

    dct = L*p0/sin(phi)*(-log(p0+sf1)+log(f2+sf5))/dgammaMax;
    dpr = p0/sf5;

    ps[x_] += dct*ps[px_]; ps[px_] = dpr*ps[px_];
    ps[y_] += dct*ps[py_]; ps[py_] = dpr*ps[py_];
    ps[delta_] =
      (2e0*dgammaMax*M_PI*cos(phi)*f2)/(lambda*f5)*ps[ct_]
      + sqr(p0)*f2/(sf1*f5)*ps[delta_];
  } else {
    if (fabs(sin(phi)) > 1e-6)
      alpha = log(gamma1/gamma)/(2e0*sqrt(2e0)*sin(phi));
    else
      alpha = dgammaMax/(gamma*2e0*sqrt(2e0));

    ps0 = ps;
    // Matrix formalism is for [x, x', y, y'] vs. canonical variables
    // [x, p_x, y, p_y].
    ps[x_] =
      cos(alpha)*ps0[x_]
      + 2e0*sqrt(2e0)*gamma*L*sin(alpha)
      /(dgammaMax*(1e0+ps0[delta_]))*ps0[px_];
    ps[px_] =
      -dgammaMax/(2e0*sqrt(2e0)*L*gamma1)*sin(alpha)*(1e0+ps0[delta_])*ps0[x_]
      + gamma/gamma1*cos(alpha)*ps0[px_];
    ps[y_] =
      cos(alpha)*ps0[y_]
      + 2e0*sqrt(2e0)*gamma*L*sin(alpha)
      /(dgammaMax*(1e0+ps0[delta_]))*ps0[py_]; 
    ps[py_] =
      -dgammaMax/(2e0*sqrt(2e0)*L*gamma1)*sin(alpha)*(1e0+ps0[delta_])*ps0[y_]
      + gamma/gamma1*cos(alpha)*ps0[py_];

   // Lattice.param.Energy contains p_0 [GeV].
    ps[delta_] =
      2e0*M_PI*this->freq*dgammaMax*cos(phi)/(c0*gamma1)*ps0[ct_]
      + 1e0/(1e0+dp/p0)*ps0[delta_];
  }

  if (this->exit_focus) Cav_Focus(L, dgamma/(gamma+dgamma), false, ps);

  if (false) {
    // Update p_0.
    p_t1 = is_double<T>::cst(ps[delta_]);
    ps[delta_] -= p_t1;
    // Lattice.param.Energy contains p_0.
    Lattice.param.Energy *= sqrt(1e0+2e0*p_t1/Lattice.param.beta0+sqr(p_t1));
    Lattice.param.gamma0 = sqrt(sqr(m_e)+sqr(p0))/m_e;
    Lattice.param.beta0  = sqrt(1e0-1e0/sqr(Lattice.param.gamma0));
    printf("\np0 = %12.5e, beta0 = %12.5e, gamma0 = %12.5e\n",
	   Lattice.param.Energy, Lattice.param.beta0, Lattice.param.gamma0);
  }
}

#endif

template<typename T>
inline void WigglerType::get_Axy(const double z, ss_vect<T> &x,
				 T AxoBrho[], T AyoBrho[])

{
  int    i;
  double ky, kz_n;
  T      cx, cz, sx, sz, chy, shy;

  for (i = 0; i <= 3; ++i) {
    AxoBrho[i] = 0e0; AyoBrho[i] = 0e0;
  }

  for (i = 0; i < n_harm; i ++) {
    kz_n = harm[i]*2e0*M_PI/lambda; ky = sqrt(sqr(kxV[i])+sqr(kz_n));
    cx = cos(kxV[i]*x[x_]); sx = sin(kxV[i]*x[x_]);
    chy = cosh(ky*x[y_]); shy = sinh(ky*x[y_]); sz = sin(kz_n*z);

    AxoBrho[0] += BoBrhoV[i]/kz_n*cx*chy*sz;
    AyoBrho[0] += BoBrhoV[i]*kxV[i]/(ky*kz_n)*sx*shy*sz;

    // derivatives with respect to x
    AxoBrho[1] -= BoBrhoV[i]*kxV[i]/kz_n*sx*chy*sz;
    AyoBrho[1] += BoBrhoV[i]*sqr(kxV[i])/(ky*kz_n)*cx*shy*sz;

    // derivatives with respect to y
    AxoBrho[2] += BoBrhoV[i]*ky/kz_n*cx*shy*sz;
    AyoBrho[2] += BoBrhoV[i]*kxV[i]/kz_n*sx*chy*sz;

    if (Lattice.param.radiation) {
      cz = cos(kz_n*z);
      // derivatives with respect to z
      AxoBrho[3] += BoBrhoV[i]*cx*chy*cz;
      AyoBrho[3] += BoBrhoV[i]*kxV[i]/ky*sx*shy*cz;
    }
  }
}

/*
template<typename T>
inline void get_Axy_map(const FieldMapType *FM, const double z,
			const ss_vect<T> &x, T AxoBrho[], T AyoBrho[])
{
  float  y, ax0, ax1, ax2, ay0, ay1, ay2;

  const  float dy = 1e-3, dz = 1e-3;

  y = is_double<T>::cst(x[y_]);

  if ((z < FM->s_pos[1]) || (z > FM->s_pos[FM->n_s])) {
    std::cout << std::scientific << std::setprecision(3)
	 << "get_Axy_map: s out of range " << z << std::endl;
    exit_(1);
  }

  if ((y < FM->y_pos[1]) || (y > FM->y_pos[FM->m_y])) {
    std::cout << std::scientific << std::setprecision(3)
	 << "get_Axy_map: y out of range " << y << std::endl;
    exit_(1);
  }

  splin2(FM->y_pos, FM->s_pos, FM->AxoBrho, FM->AxoBrho2, FM->m_y, FM->n_s,
	 y, z, &ax1);
  AxoBrho[0] = FM->scl*ax1;

  splin2(FM->y_pos, FM->s_pos, FM->AyoBrho, FM->AyoBrho2, FM->m_y, FM->n_s,
	 y, z, &ay1);
  AyoBrho[0] = FM->scl*ay1;

  // derivatives with respect to x
  AxoBrho[1] = FM->scl*0e0; AyoBrho[1] = FM->scl*0e0;

  // derivatives with respect to y
  splin2(FM->y_pos, FM->s_pos, FM->AxoBrho, FM->AxoBrho2, FM->m_y, FM->n_s,
	 y+dy, z, &ax2);
  splin2(FM->y_pos, FM->s_pos, FM->AxoBrho, FM->AxoBrho2, FM->m_y, FM->n_s,
	 y-dy, z, &ax1);
  splin2(FM->y_pos, FM->s_pos, FM->AxoBrho, FM->AxoBrho2, FM->m_y, FM->n_s,
	 y, z, &ax0);
  AxoBrho[2] =
    (ax2-ax1)/(2e0*dy) + (ax2+ax1-2e0*ax0)/sqr(dy)*is_tps<T>::set_prm(y_+1);
  AxoBrho[2] *= FM->scl;

  splin2(FM->y_pos, FM->s_pos, FM->AyoBrho, FM->AyoBrho2, FM->m_y, FM->n_s,
	 y+dy, z, &ay2);
  splin2(FM->y_pos, FM->s_pos, FM->AyoBrho, FM->AyoBrho2, FM->m_y, FM->n_s,
	 y-dy, z, &ay1);
  splin2(FM->y_pos, FM->s_pos, FM->AyoBrho, FM->AyoBrho2, FM->m_y, FM->n_s,
	 y, z, &ay0);
  AyoBrho[2] =
    (ay2-ay1)/(2e0*dy) + (ay2+ay1-2e0*ay0)/sqr(dy)*is_tps<T>::set_prm(y_+1);
  AyoBrho[2] *= FM->scl;

  if (Lattice.param.radiation) {
    // derivatives with respect to z
    splin2(FM->y_pos, FM->s_pos, FM->AxoBrho, FM->AxoBrho2, FM->m_y, FM->n_s,
	   y, z+dz, &ax2);
    splin2(FM->y_pos, FM->s_pos, FM->AxoBrho, FM->AxoBrho2, FM->m_y, FM->n_s,
	   y, z-dz, &ax1);
    AxoBrho[3] = (ax2-ax1)/(2e0*dz); AxoBrho[3] *= FM->scl;

    splin2(FM->y_pos, FM->s_pos, FM->AyoBrho, FM->AyoBrho2, FM->m_y, FM->n_s,
	   y, z+dz, &ay2);
    splin2(FM->y_pos, FM->s_pos, FM->AyoBrho, FM->AyoBrho2, FM->m_y, FM->n_s,
	   y, z-dz, &ay1);
    AyoBrho[3] = (ay2-ay1)/(2e0*dz); AyoBrho[3] *= FM->scl;
    if (false)
      std::cout << std::fixed << std::setprecision(5)
	   << std::setw(8) << z << std::setw(9)
	   << is_double<T>::cst(AxoBrho[3]) << std::endl;
  }
}
*/

template<typename T>
void WigglerType::Wiggler_pass_EF(ss_vect<T> &x)
{
  // First order symplectic integrator for wiggler using expanded Hamiltonian

  int          i, nstep = 0;
  double       h, z;
  T            AxoBrho[4] = {0e0, 0e0, 0e0, 0e0},
               AyoBrho[4] = {0e0, 0e0, 0e0, 0e0},
               psi, hodp, a12, a21, a22, det,
               d1, d2, a11, c11, c12, c21, c22, x2, B[3];

  switch (Elem.Kind) {
  case Wigl:
    nstep = N;
    break;
  case FieldMap:
    printf("Wiggler_pass_EF: problem\n");
    exit(1);
    // FM = static_cast<FieldMapType*>(&Cell1);

    // nstep = FM->n_step;
    break;
  default:
    std::cout << "Wiggler_pass_EF: unknown element type" << std::endl;
    exit_(1);
    break;
  }

  h = L/nstep; z = 0e0;
  for (i = 1; i <= nstep; ++i) {
    switch (Elem.Kind) {
    case Wigl:
      get_Axy(z, x, AxoBrho, AyoBrho);
      break;
    case FieldMap:
//      get_Axy_map(Cell.Elem.FM, z, x, AxoBrho, AyoBrho);
      break;
    default:
      std::cout << "Wiggler_pass_EF: unknown element type" << std::endl;
      exit_(1);
      break;
    }

    psi = 1e0 + x[delta_]; hodp = h/psi;
    a11 = hodp*AxoBrho[1]; a12 = hodp*AyoBrho[1];
    a21 = hodp*AxoBrho[2]; a22 = hodp*AyoBrho[2];
    det = 1e0 - a11 - a22 + a11*a22 - a12*a21;
    d1 = hodp*AxoBrho[0]*AxoBrho[1]; d2 = hodp*AxoBrho[0]*AxoBrho[2];
    c11 = (1e0-a22)/det; c12 = a12/det; c21 = a21/det; c22 = (1e0-a11)/det;
    x2 = c11*(x[px_]-d1) + c12*(x[py_]-d2);

    x[py_] = c21*(x[px_]-d1) + c22*(x[py_]-d2); x[px_] = x2;
    x[x_] += hodp*(x[px_]-AxoBrho[0]); x[y_] += hodp*x[py_];
    x[ct_] += h*(sqr((x[px_]-AxoBrho[0])/psi)
	      + sqr((x[py_]-AyoBrho[0])/psi))/2e0;

    if (false)
      std::cout << std::scientific << std::setprecision(3)
	   << std::setw(8) << z
	   << std::setw(11) << is_double<T>::cst(x[x_])
	   << std::setw(11) << is_double<T>::cst(x[px_])
	   << std::setw(11) << is_double<T>::cst(x[y_])
	   << std::setw(11) << is_double<T>::cst(x[py_])
	   << std::endl;

    if (Lattice.param.pathlength) x[ct_] += h;

    if (Lattice.param.radiation || Lattice.param.emittance) {
      B[X_] = -AyoBrho[3]; B[Y_] = AxoBrho[3]; B[Z_] = AyoBrho[1] - AxoBrho[2];
      radiate(x, h, 0e0, B);
    }

    z += h;
  }
}


template<typename T>
inline void WigglerType::get_Axy2(const double z, const double kz,
				  ss_vect<T> &x, T AxoBrho[], T AyoBrho[])
{
  int i;
  T   cx, sx, cz1, cz2, sz1, sz2, chy, shy, kyH, kyV, chx, shx, cy, sy;

  for (i = 0; i <= 3; ++i) {
    AxoBrho[i] = 0e0; AyoBrho[i] = 0e0;
  }

  kyV = sqrt(sqr(kz)+sqr(kxV[0])); kyH = sqrt(sqr(kz)+sqr(kxH[0]));
  cx = cos(kxV[0]*x[x_]); sx = sin(kxV[0]*x[x_]);
  cy = cos(kxH[0]*x[y_]); sy = sin(kxH[0]*x[y_]);
  chx = cosh(kyH*x[x_]); shx = sinh(kyH*x[x_]);
  chy = cosh(kyV*x[y_]); shy = sinh(kyV*x[y_]);
  sz1 = sin(kz*z); sz2 = sin(kz*z+phi[0]);

  AxoBrho[0] += BoBrhoV[0]/kz*cx*chy*sz1;
  AxoBrho[0] -= BoBrhoH[0]*kxH[0]/(kyH*kz)*shx*sy*sz2;
  AyoBrho[0] += BoBrhoV[0]*kxV[0]/(kyV*kz)*sx*shy*sz1;
  AyoBrho[0] -= BoBrhoH[0]/kz*chx*cy*sz2;

  /* derivatives with respect to x */
  AxoBrho[1] -= BoBrhoV[0]*kxV[0]/kz*sx*chy*sz1;
  AxoBrho[1] -= BoBrhoH[0]*kxH[0]/kz*chx*sy*sz2;
  AyoBrho[1] += BoBrhoV[0]*sqr(kxV[0])/(kyV*kz)*cx*shy*sz1;
  AyoBrho[1] -= BoBrhoH[0]*kyH/kz*shx*cy*sz2;

  /* derivatives with respect to y */
  AxoBrho[2] += BoBrhoV[0]*kyV/kz*cx*shy*sz1;
  AxoBrho[2] -= BoBrhoH[0]*sqr(kxH[0])/(kyH*kz)*shx*cy*sz2;
  AyoBrho[2] += BoBrhoV[0]*kxV[0]/kz*sx*chy*sz1;
  AyoBrho[2] += BoBrhoH[0]*kxH[0]/kz*chx*sy*sz2;

  if (Lattice.param.radiation) {
    cz1 = cos(kz*z); cz2=cos(kz*z+phi[0]);
    /* derivatives with respect to z */
    AxoBrho[3] += BoBrhoV[0]*cx*chy*cz1;
    AxoBrho[3] -= BoBrhoH[0]*kxH[0]/kyH*shx*sy*cz2;
    AyoBrho[3] += BoBrhoV[0]*kxV[0]/kyV*sx*shy*cz1;
    AyoBrho[3] -= BoBrhoH[0]*chx*cy*cz2;
  }
}


template<typename T>
void WigglerType::Wiggler_pass_EF2(ss_vect<T> &x)
{
  // First order symplectic integrator for wiggler using expanded Hamiltonian

  int    i;
  double h, z, kz;
  T      hodp, B[3], px1, px2, px3, py1, py2, AxoBrho[4], AyoBrho[4], psi;
  T      px = 0e0, py = 0e0;

  kz = 2e0*M_PI/this->lambda;
  h = L/N; z = 0e0;
  for (i = 1; i <= N; ++i) {
    get_Axy2(z, kz, x, AxoBrho, AyoBrho);

    psi = 1e0 + x[delta_]; hodp = h/psi;

    px1 = (x[px_]-(AxoBrho[0]*AxoBrho[1]+AyoBrho[0]*AyoBrho[1])*hodp)
          *(1-AyoBrho[2]*hodp);
    px2 = (x[py_]-(AxoBrho[0]*AxoBrho[2]+AyoBrho[0]*AyoBrho[2])*hodp)
          *AyoBrho[1]*hodp;
    px3 = (1-AxoBrho[1]*hodp)*(1-AyoBrho[2]*hodp)
          - AxoBrho[2]*AyoBrho[1]*hodp*hodp;

    py1 = (x[py_]-(AxoBrho[0]*AxoBrho[2]+AyoBrho[0]*AyoBrho[2])*hodp)
          *(1-AxoBrho[1]*hodp);
    py2 = (x[px_]-(AxoBrho[0]*AxoBrho[1]+AyoBrho[0]*AyoBrho[1])*hodp)
          *AxoBrho[2]*hodp;

    py = (py1+py2)/px3; px = (px1+px2)/px3;
    x[x_] += hodp*(px-AxoBrho[0]); x[y_] += hodp*(py-AyoBrho[0]);
    x[ct_] += h*(sqr((px-AxoBrho[0])/psi) + sqr((py-AyoBrho[0])/psi))/2e0;

    if (Lattice.param.pathlength) x[ct_] += h;

    if (Lattice.param.radiation || Lattice.param.emittance) {
      B[X_] = -AyoBrho[3]; B[Y_] = AxoBrho[3]; B[Z_] = AyoBrho[1] - AxoBrho[2];
      radiate(x, h, 0e0, B);
    }

    z += h;
  }

  x[px_] = px; x[py_] = py;
}


template<typename T>
inline void WigglerType::get_Axy_EF3(const double z, const ss_vect<T> &x,
				     T &AoBrho, T dAoBrho[], T &dp,
				     const bool hor)
{
  int    i;
  double ky, kz_n;
  T      cx, sx, sz, chy, shy, cz;

  AoBrho = 0e0; dp = 0e0;

  for (i = 0; i < 3; i++)
    dAoBrho[i] = 0e0;

  for (i = 0; i < n_harm; i++) {
    kz_n = harm[i]*2e0*M_PI/lambda; ky = sqrt(sqr(kxV[i])+sqr(kz_n));

    cx  = cos(kxV[i]*x[x_]); sx = sin(kxV[i]*x[x_]);
    chy = cosh(ky*x[y_]); shy = sinh(ky*x[y_]); sz = sin(kz_n*z);

    if (hor) {
      // A_x/Brho
      AoBrho += BoBrhoV[i]/kz_n*cx*chy*sz;

      if (Lattice.param.radiation) {
	cz = cos(kz_n*z);
	dAoBrho[X_] -= BoBrhoV[i]*kxV[i]/kz_n*sx*chy*sz;
	dAoBrho[Y_] += BoBrhoV[i]*ky/kz_n*cx*shy*sz;
	dAoBrho[Z_] += BoBrhoV[i]*cx*chy*cz;
      }

      // dp_y
      if (kxV[i] == 0e0)
	dp += BoBrhoV[i]/kz_n*ky*x[x_]*shy*sz;
      else
	dp += BoBrhoV[i]/(kxV[i]*kz_n)*ky*sx*shy*sz;
    } else {
      // A_y/Brho
      AoBrho += BoBrhoV[i]*kxV[i]/(ky*kz_n)*sx*shy*sz;

      if (Lattice.param.radiation) {
	cz = cos(kz_n*z);
	dAoBrho[X_] += BoBrhoV[i]*sqr(kxV[i])/(ky*kz_n)*cx*shy*sz;
	dAoBrho[Y_] += BoBrhoV[i]*kxV[i]/kz_n*sx*chy*sz;
	dAoBrho[Z_] += BoBrhoV[i]*kxV[i]/ky*sx*shy*cz;
      }

      // dp_x
      dp += BoBrhoV[i]/kz_n*sqr(kxV[i]/ky)*cx*chy*sz;
    }
  }
}


template<typename T>
void WigglerType::Wiggler_pass_EF3(ss_vect<T> &x)
{
  /* Second order symplectic integrator for insertion devices based on:

       E. Forest, et al "Explicit Symplectic Integrator for s-dependent
       Static Magnetic Field"                                                */

  int    i;
  double h, z;
  T      hd, AxoBrho, AyoBrho, dAxoBrho[3], dAyoBrho[3], dpy, dpx, B[3];

  h = L/N; z = 0e0;
  for (i = 1; i <= N; i++) {
    hd = h/(1e0+x[delta_]);

    // 1: half step in z
    z += 0.5*h;

    // 2: half drift in y
    get_Axy_EF3(z, x, AyoBrho, dAyoBrho, dpx, false);

    x[px_] -= dpx; x[py_] -= AyoBrho;
    x[y_] += 0.5*hd*x[py_];
    x[ct_] += sqr(0.5)*hd*sqr(x[py_])/(1e0+x[delta_]);

    get_Axy_EF3(z, x, AyoBrho, dAyoBrho, dpx, false);

    x[px_] += dpx; x[py_] += AyoBrho;

    // 3: full drift in x
    get_Axy_EF3(z, x, AxoBrho, dAxoBrho, dpy, true);

    x[px_] -= AxoBrho; x[py_] -= dpy; x[x_] += hd*x[px_];
    x[ct_] += 0.5*hd*sqr(x[px_])/(1e0+x[delta_]);

    if (Lattice.param.pathlength) x[ct_] += h;

    get_Axy_EF3(z, x, AxoBrho, dAxoBrho, dpy, true);

    x[px_] += AxoBrho; x[py_] += dpy;

    // 4: a half drift in y
    get_Axy_EF3(z, x, AyoBrho, dAyoBrho, dpx, false);

    x[px_] -= dpx; x[py_] -= AyoBrho;
    x[y_] += 0.5*hd*x[py_];
    x[ct_] += sqr(0.5)*hd*sqr(x[py_])/(1e0+x[delta_]);

    get_Axy_EF3(z, x, AyoBrho, dAyoBrho, dpx, false);

    x[px_] += dpx; x[py_] += AyoBrho;

    // 5: half step in z
    z += 0.5*h;

    if (Lattice.param.radiation || Lattice.param.emittance) {
      get_Axy_EF3(z, x, AyoBrho, dAyoBrho, dpx, false);
      get_Axy_EF3(z, x, AxoBrho, dAxoBrho, dpy, true);
      B[X_] = -dAyoBrho[Z_]; B[Y_] = dAxoBrho[Z_];
      B[Z_] = dAyoBrho[X_] - dAxoBrho[Y_];
      radiate(x, h, 0e0, B);
    }
  }
}


template<typename T>
void WigglerType::Elem_Propagate(ss_vect<T> &X)
{
  int         seg;
  double      L, L1, L2, K1, K2;
  ss_vect<T>  X1;

  // Global -> Local
  this->GtoL(X, 0e0, 0e0, 0e0);
  switch (this->method) {

  case Meth_Linear:
    std::cout << "Wiggler_Pass: Meth_Linear not supported" << std::endl;
    exit_(1);
    break;

  case Meth_First:
    if ((this->BoBrhoV[0] != 0e0) || (this->BoBrhoH[0] != 0e0)) {
      if (!Lattice.param.EPU)
	Wiggler_pass_EF(X);
      else {
	Wiggler_pass_EF2(X);
      }
    } else
      // drift if field = 0
      Drift(this->L, X);
    break;

  case Meth_Second:
    if ((this->BoBrhoV[0] != 0e0) || (this->BoBrhoH[0] != 0e0)) {
      Wiggler_pass_EF3(X);
    } else
      // drift if field = 0
      Drift(this->L, X);
    break;

  case Meth_Fourth:  /* 4-th order integrator */
    L = this->L/this->N;
    L1 = c_1*L; L2 = c_2*L; K1 = d_1*L; K2 = d_2*L;
    for (seg = 1; seg <= this->N; seg++) {
      Drift(L1, X); X1 = X;
      thin_kick(this->order, this->BW, K1, 0e0, 0e0, X1);
      X[py_] = X1[py_]; Drift(L2, X); X1 = X;
      thin_kick(this->order, this->BW, K2, 0e0, 0e0, X1);
      X[py_] = X1[py_]; Drift(L2, X); X1 = X;
      thin_kick(this->order, this->BW, K1, 0e0, 0e0, X1);
      X[py_] = X1[py_]; Drift(L1, X);
    }
    break;
  }
  // Local -> Global
  this->LtoG(X, 0e0, 0e0, 0e0);
}

#undef eps
#undef kx


template<typename T>
void FieldMapType::rk4_(const ss_vect<T> &y, const ss_vect<T> &dydx,
			const double x, const double h, ss_vect<T> &yout,
			void (*derivs)(const double, const ss_vect<T> &,
				       ss_vect<T> &))
{
  double     xh, hh, h6;
  ss_vect<T> dym, dyt, yt;

  hh = h*0.5; h6 = h/6e0;
  xh = x + hh; yt = y + hh*dydx;
  (*derivs)(xh, yt, dyt); yt = y + hh*dyt;
  (*derivs)(xh, yt, dym); yt = y + h*dym; dym += dyt;
  (*derivs)(x+h, yt, dyt);
  yout = y + h6*(dydx+dyt+2e0*dym);
}


template<typename T>
inline T get_p_s_ps(const ss_vect<T> &x, const T &qop_Ax, const T &qop_Ay)
{
  // Compute p_s in phase space (with vector potential in axial gauge).

  return sqrt(sqr(1e0+x[delta_])-sqr(x[px_]-qop_Ax)-sqr(x[py_]-qop_Ay));
}


template<typename T>
void FieldMapType::f_FM(const double z, const ss_vect<T> &ps, ss_vect<T> &Dps)
{
  // Coordinates are: [x, x', y, y', -ct, delta].
  int          j, kz;
  T            BoBrho1[3], p_s;

  const double eps = 1e-5;

  kz = 0;
  for (j = 1; j <= n[Z_]; j++)
    if (fabs(z-x[Z_][j]) < eps) {
      kz = j;
      break;
    }

  if (kz == 0) {
    std::cout << std::fixed << std::setprecision(10)
	 << "z = " << std::setw(12) << z << " undefined" << std::endl;

    for (j = 0; j < ss_dim; j++)
      Dps[j] = NAN;

    return;
  }

  splin2_(x[X_], x[Y_], BoBrho[X_][kz], BoBrho2[X_][kz],
	   n[X_], n[Y_], ps[x_], ps[y_], BoBrho1[X_]);

  if (BoBrho1[X_] == NAN) {
    for (j = 0; j < ss_dim; j++)
      Dps[j] = NAN;
    return;
  }

  splin2_(x[X_], x[Y_], BoBrho[Y_][kz], BoBrho2[Y_][kz],
	   n[X_], n[Y_], ps[x_], ps[y_], BoBrho1[Y_]);

  if (BoBrho1[Y_] == NAN) {
    for (j = 0; j < ss_dim; j++)
      Dps[j] = NAN;
    return;
  }

  splin2_(x[X_], x[Y_], BoBrho[Z_][kz], BoBrho2[Z_][kz],
	   n[X_], n[Y_], ps[x_], ps[y_], BoBrho1[Z_]);

  if (BoBrho1[Z_] == NAN) {
    for (j = 0; j < ss_dim; j++)
      Dps[j] = NAN;
    return;
  }

  p_s = get_p_s_cs(ps);

  Dps[x_] = ps[px_];
  Dps[px_] =
    -(ps[px_]*ps[py_]*BoBrho1[X_]-(1e0+sqr(ps[px_]))*BoBrho1[Y_]
      +ps[py_]*BoBrho1[Z_])/p_s;

  Dps[y_] = ps[py_];
  Dps[py_] =
    -((1e0+sqr(ps[py_]))*BoBrho1[X_]-ps[px_]*ps[py_]*BoBrho1[Y_]
      -ps[px_]*BoBrho1[Z_])/p_s;

  Dps[ct_] = (1e0+ps[delta_])/p_s - 1e0;

  if (Lattice.param.pathlength) Dps[ct_] += 1e0;

  Dps[delta_] = 0e0;
}


template<typename T>
inline T get_p_s_cs(const ss_vect<T> &x)
{
  // Compute p_s in configuration space.

  return (1e0+x[delta_])/sqrt(1e0+sqr(x[px_])+sqr(x[py_]));
}


template<typename T>
void FieldMapType::FieldMap_pass_RK(ss_vect<T> &ps, int k)
{
  int          i;
  double       h, z;
  T            p_s;
  ss_vect<T>   Dps;

  const int n_step = 2; // Each step needs: f(z_n), f(z_n+h), f(z_n+2h)

  switch (Lattice.param.FieldMap_filetype) {
  case 2:
  case 3:
  case 4:
    // Transform to right handed system
    ps[x_] = -ps[x_]; ps[px_] = -ps[px_];
    break;
  }

  // transform to [x, x', y, y']
  p_s = get_p_s(ps); ps[px_] /= p_s; ps[py_] /= p_s;

  h = n_step*dx[Z_]; z = x[Z_][1]; Lr = 0e0;
  if (trace)
    outf_ << std::scientific << std::setprecision(3)
	  << std::setw(11) << s_FM
	  << std::setw(11) << is_double<T>::cst(ps[x_])
	  << std::setw(11) << is_double<T>::cst(ps[px_])
	  << std::setw(11) << is_double<T>::cst(ps[y_])
	  << std::setw(11) << is_double<T>::cst(ps[py_])
	  << std::setw(11) << is_double<T>::cst(ps[delta_])
	  << std::setw(11) << is_double<T>::cst(ps[ct_])
	  << std::endl;
  for(i = 1+cut; i < n[Z_]-cut; i += n_step) {
    if (i <= n[Z_]-cut-2) {
      f_FM(z, ps, Dps);

      if (Dps[x_] == NAN) {
	std::cout << "FieldMap_pass_RK: particle lost" << std::endl;
	std::cout << ps;
	return;
      }

      printf("FieldMap_pass_RK: problem\n");
      exit(1);
      // rk4_(ps, Dps, x[Z_][i], h, ps, f_FM);

      z += h; Lr += h; s_FM += h;
    } else {
      // Use 2nd order Runge-Kutta (aka Midpoint method)
      f_FM(z, ps, Dps);

      if (Dps[x_] == NAN) {
	std::cout << "FieldMap_pass_RK: particle lost" << std::endl;
	std::cout << ps;
	return;
      }

      ps += h/2e0*Dps;

      z += h/2e0; Lr += h/2e0; s_FM += h/2e0;
    }

    if (trace)
      outf_ << std::scientific << std::setprecision(3)
	    << std::setw(11) << s_FM
	    << std::setw(11) << is_double<T>::cst(ps[x_])
	    << std::setw(11) << is_double<T>::cst(ps[px_])
	    << std::setw(11) << is_double<T>::cst(ps[y_])
	    << std::setw(11) << is_double<T>::cst(ps[py_])
	    << std::setw(11) << is_double<T>::cst(ps[delta_])
	    << std::setw(11) << is_double<T>::cst(ps[ct_])
	    << std::endl;
  }

  // transform back to [x, px, y, py] (A_x,y,z = 0)
  p_s = get_p_s_cs(ps); ps[px_] *= p_s; ps[py_] *= p_s;

  switch (Lattice.param.FieldMap_filetype) {
  case 2:
  case 3:
  case 4:
    // Transform back to left handed system
    ps[x_] = -ps[x_]; ps[px_] = -ps[px_];
    break;
  }
}


template<typename T>
void FieldMapType::FieldMap_pass_SI(ss_vect<T> &ps, int k)
{
  /* E. Chacon-Golcher, F. Neri "A Symplectic Integrator with Arbitrary
     Vector and Scalar Potentials" hys. Lett. A 372 p. 4661-4666 (2008).    */

  int          i, j = 0;
  double       h, z;
  T            hd, AoBrho1[2], dAoBrho[2], AoBrho_int, ByoBrho;
  ss_vect<T>   ps1;

  const int    n_step = 2;
  const double d_diff = 1e0;


  switch (Lattice.param.FieldMap_filetype) {
  case 2:
  case 3:
  case 4:
    // Transform to right handed system
    ps[x_] = -ps[x_]; ps[px_] = -ps[px_];
    break;
  }

  h = n_step*dx[Z_]; z = 0e0; Lr = 0e0;
  if (false)
    outf_ << std::scientific << std::setprecision(3)
	  << std::setw(11) << s_FM
	  << std::setw(11) << is_double<T>::cst(ps[x_])
	  << std::setw(11) << is_double<T>::cst(ps[px_])
	  << std::setw(11) << is_double<T>::cst(ps[y_])
	  << std::setw(11) << is_double<T>::cst(ps[py_])
	  << std::setw(11) << is_double<T>::cst(ps[delta_])
	  << std::setw(11) << is_double<T>::cst(ps[ct_])
	  << std::endl;
  for (i = 1+cut; i < n[Z_]-cut; i += n_step) {
    hd = h/(1e0+ps[delta_]);

    // 1. Half step in z.
    z += 0.5*h; j = i + 1; s_FM += 0.5*h;

    // 2. Half drift in y.
    ps1 = ps;

    splin2_(x[X_], x[Y_], AoBrho[Y_][j], AoBrho2[Y_][j],
	    n[X_], n[Y_], ps[x_], ps[y_], AoBrho1[0]);

    if (AoBrho1[0] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    ps1[y_] += (ps[py_]-AoBrho1[0])*0.5*hd;
    ps1[ct_] += sqr(0.5)*hd*sqr(ps[py_]-AoBrho1[0])/(1e0+ps[delta_]);

    splin2_(x[X_], x[Y_], AoBrho[Y_][j], AoBrho2[Y_][j],
	    n[X_], n[Y_], ps[x_]+d_diff*dx[X_], ps[y_],
	    dAoBrho[1]);

    if (dAoBrho[1] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    splin2_(x[X_], x[Y_], AoBrho[Y_][j], AoBrho2[Y_][j],
	    n[X_], n[Y_], ps[x_]-d_diff*dx[X_], ps[y_],
	    dAoBrho[0]);

    if (dAoBrho[0] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    AoBrho_int = (dAoBrho[1]-dAoBrho[0])/(2e0*d_diff*dx[X_]);

    splin2_(x[X_], x[Y_], AoBrho[Y_][j], AoBrho2[Y_][j],
	    n[X_], n[Y_], ps1[x_]+d_diff*dx[X_], ps1[y_],
	    dAoBrho[1]);

    if (dAoBrho[1] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    splin2_(x[X_], x[Y_], AoBrho[Y_][j], AoBrho2[Y_][j],
	    n[X_], n[Y_], ps1[x_]-d_diff*dx[X_], ps1[y_],
	    dAoBrho[0]);

    if (dAoBrho[0] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    AoBrho_int += (dAoBrho[1]-dAoBrho[0])/(2e0*d_diff*dx[X_]);

    // Trapezoidal rule
    ps1[px_] +=
      AoBrho_int*(is_double<T>::cst(ps1[y_])-is_double<T>::cst(ps[y_]))/2e0;

    splin2_(x[X_], x[Y_], AoBrho[Y_][j], AoBrho2[Y_][j],
	    n[X_], n[Y_], ps1[x_], ps1[y_], AoBrho1[1]);

    if (AoBrho1[1] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    ps1[py_] += AoBrho1[1] - AoBrho1[0];

    ps = ps1;

    // 3. Full drift in x.
    ps1 = ps;

    splin2_(x[X_], x[Y_], AoBrho[X_][j], AoBrho2[X_][j],
	    n[X_], n[Y_], ps[x_], ps[y_], AoBrho1[0]);

    if (AoBrho1[0] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    ps1[x_] += (ps[px_]-AoBrho1[0])*hd;
    ps1[ct_] += 0.5*hd*sqr(ps[px_]-AoBrho1[0])/(1e0+ps[delta_]);

    splin2_(x[X_], x[Y_], AoBrho[X_][j], AoBrho2[X_][j],
	    n[X_], n[Y_], ps1[x_], ps1[y_], AoBrho1[1]);

    if (AoBrho1[1] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    ps1[px_] += AoBrho1[1] - AoBrho1[0];

    splin2_(x[X_], x[Y_], AoBrho[X_][j], AoBrho2[X_][j],
	    n[X_], n[Y_], ps[x_], ps[y_]+d_diff*dx[Y_],
	    dAoBrho[1]);

    if (dAoBrho[1] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    splin2_(x[X_], x[Y_], AoBrho[X_][j], AoBrho2[X_][j],
	    n[X_], n[Y_], ps[x_], ps[y_]-d_diff*dx[Y_],
	    dAoBrho[0]);

    if (dAoBrho[0] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    AoBrho_int = (dAoBrho[1]-dAoBrho[0])/(2e0*d_diff*dx[Y_]);

    splin2_(x[X_], x[Y_], AoBrho[X_][j], AoBrho2[X_][j],
	    n[X_], n[Y_], ps1[x_], ps1[y_]+d_diff*dx[Y_],
	    dAoBrho[1]);

    if (dAoBrho[1] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    splin2_(x[X_], x[Y_], AoBrho[X_][j], AoBrho2[X_][j],
	    n[X_], n[Y_], ps1[x_], ps1[y_]-d_diff*dx[Y_],
	    dAoBrho[0]);

    if (dAoBrho[0] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    AoBrho_int += (dAoBrho[1]-dAoBrho[0])/(2e0*d_diff*dx[Y_]);

    // Trapezoidal rule
    ps1[py_] +=
      AoBrho_int*(is_double<T>::cst(ps1[x_])-is_double<T>::cst(ps[x_]))/2e0;

    ps = ps1;

    // 4. Half drift in y.
    ps1 = ps;

    splin2_(x[X_], x[Y_], AoBrho[Y_][j], AoBrho2[Y_][j],
	    n[X_], n[Y_], ps[x_], ps[y_], AoBrho1[0]);

    if (AoBrho1[0] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    ps1[y_] += (ps[py_]-AoBrho1[0])*0.5*hd;
    ps1[ct_] += sqr(0.5)*hd*sqr(ps[py_]-AoBrho1[0])/(1e0+ps[delta_]);

    splin2_(x[X_], x[Y_], AoBrho[Y_][j], AoBrho2[Y_][j],
	    n[X_], n[Y_], ps[x_]+d_diff*dx[X_], ps[y_],
	    dAoBrho[1]);

    if (dAoBrho[1] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    splin2_(x[X_], x[Y_], AoBrho[Y_][j], AoBrho2[Y_][j],
	    n[X_], n[Y_], ps[x_]-d_diff*dx[X_], ps[y_],
	    dAoBrho[0]);

    if (dAoBrho[0] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    AoBrho_int = (dAoBrho[1]-dAoBrho[0])/(2e0*d_diff*dx[X_]);

    splin2_(x[X_], x[Y_], AoBrho[Y_][j], AoBrho2[Y_][j],
	    n[X_], n[Y_], ps1[x_]+d_diff*dx[X_], ps1[y_],
	    dAoBrho[1]);

    if (dAoBrho[1] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    splin2_(x[X_], x[Y_], AoBrho[Y_][j], AoBrho2[Y_][j],
	    n[X_], n[Y_], ps1[x_]-d_diff*dx[X_], ps1[y_],
	    dAoBrho[0]);

    if (dAoBrho[0] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    AoBrho_int += (dAoBrho[1]-dAoBrho[0])/(2e0*d_diff*dx[X_]);

    // Trapezoidal rule
    ps1[px_] +=
      AoBrho_int*(is_double<T>::cst(ps1[y_])-is_double<T>::cst(ps[y_]))/2e0;

    splin2_(x[X_], x[Y_], AoBrho[Y_][j], AoBrho2[Y_][j],
	    n[X_], n[Y_], ps1[x_], ps1[y_], AoBrho1[1]);

    if (AoBrho1[1] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    ps1[py_] += AoBrho[1] - AoBrho[0];

    ps = ps1;

    // 5. Half step in z.
    z += 0.5*h; j = i + 2; s_FM += 0.5*h;

    if (Lattice.param.pathlength) ps[ct_] += h;

    Lr += h;

    if (Lattice.param.radiation || Lattice.param.emittance) {
//      B[X_] = -AoBrhoy[3]; B[Y_] = AoBrho[X_][3];
//      B[Z_] = AoBrhoy[1] - AoBrho[X_][2];
//      radiate(ps, h, 0e0, B);
    }

    if (trace) {
      splin2_(x[X_], x[Y_], AoBrho[X_][j], AoBrho2[X_][j],
	      n[X_], n[Y_], ps[x_], ps[y_], AoBrho1[0]);
      splin2_(x[X_], x[Y_], AoBrho[Y_][j], AoBrho2[Y_][j],
	      n[X_], n[Y_], ps[x_], ps[y_], AoBrho1[1]);
      splin2_(x[X_], x[Y_], BoBrho[Y_][j], BoBrho2[Y_][j],
	      n[X_], n[Y_], ps[x_], ps[y_], ByoBrho);

      outf_ << std::scientific << std::setprecision(3)
	    << std::setw(11) << s_FM
	    << std::setw(11) << is_double<T>::cst(ps[x_])
	    << std::setw(11) << is_double<T>::cst(ps[px_]-AoBrho1[0])
	    << std::setw(11) << is_double<T>::cst(ps[y_])
	    << std::setw(11) << is_double<T>::cst(ps[py_]-AoBrho1[1])
	    << std::setw(11) << is_double<T>::cst(ps[delta_])
	    << std::setw(11) << is_double<T>::cst(ps[ct_])
	    << std::setw(11) << is_double<T>::cst(AoBrho1[0])
	    << std::setw(11) << is_double<T>::cst(AoBrho1[1])
	    << std::setw(11) << is_double<T>::cst(ByoBrho)
	    << std::endl;
    }
  }

  // Change of gauge
  splin2_(x[X_], x[Y_], AoBrho[X_][j], AoBrho2[X_][j],
	  n[X_], n[Y_], ps[x_], ps[y_], AoBrho1[0]);

  if (AoBrho1[0] == NAN) {
    for (j = 0; j < ss_dim; j++)
      ps[j] = NAN;
    return;
  }

  ps[px_] -= AoBrho1[0];

  splin2_(x[X_], x[Y_], AoBrho[Y_][j], AoBrho2[Y_][j],
	  n[X_], n[Y_], ps[x_], ps[y_], AoBrho1[0]);

  if (AoBrho1[0] == NAN) {
    for (j = 0; j < ss_dim; j++)
      ps[j] = NAN;
    return;
  }

  ps[py_] -= AoBrho1[0];

  switch (Lattice.param.FieldMap_filetype) {
  case 2:
  case 3:
  case 4:
    // Transform back to left handed system
    ps[x_] = -ps[x_]; ps[px_] = -ps[px_];
    break;
  }
}


// Instantiate
// template void f_FM(const CellType &, const double, const ss_vect<double> &,
// 		   ss_vect<double> &);
// template void f_FM(const CellType &, const double, const ss_vect<tps> &,
// 		   ss_vect<tps> &);
// template void rk4_(const CellType &, const ss_vect<double> &,
// 		   const ss_vect<double>, const double, const double,
// 		   ss_vect<double> &,
// 		   void (*derivs)(const CellType &, const double,
// 				  const ss_vect<double> &, ss_vect<double> &));
// template void rk4_(const CellType &, const ss_vect<tps> &, const ss_vect<tps>,
// 		   const double, const double, ss_vect<tps> &,
// 		   void (*derivs)(const CellType &, const double,
// 				  const ss_vect<tps> &, ss_vect<tps> &));
// template void FieldMap_pass_RK(CellType &, ss_vect<double> &, int k);
// template void FieldMap_pass_RK(CellType &, ss_vect<tps> &, int k);
// template void FieldMap_pass_SI(CellType &, ss_vect<double> &, int k);
// template void FieldMap_pass_SI(CellType &, ss_vect<tps> &, int k);


template<typename T>
void FieldMapType::Elem_Propagate(ss_vect<T> &ps)
{
  int          k;
  double       Ld;

  if (trace & first_FM) {
    file_wr(outf_, "FieldMap_pass.dat");
    s_FM = 0e0;
    first_FM = false;
  }

//  this->GtoL(ps, 0e0, 0e0, 0e0);

  Ld = (this->Lr-this->L)/2e0;
  p_rot(this->phi/2e0*180e0/M_PI, ps);
  printf("\nFieldMap_Pass: entrance negative drift [m] %12.5e", Ld);
  Drift(-Ld, ps);

  for (k = 1; k <= this->n_step; k++) {
    if (Lattice.param.sympl)
      FieldMap_pass_SI(ps, k);
    else
      FieldMap_pass_RK(ps, k);
  }

  printf("\nFieldMap_Pass: exit negative drift [m]     %12.5e", Ld);
  Drift(-Ld, ps);
  p_rot(this->phi/2e0*180e0/M_PI, ps);

//  this->LtoG(ps, 0e0, 0e0, 0e0);

//  outf_.close();
}


template<typename T>
void InsertionType::Elem_Propagate(ss_vect<T> &x)
{
  /* Purpose:
       Track vector x through a insertion
       If radiation or cavity on insertion is like a drift

   Input:
       Cell element to track through
       x initial coordinates vector

   Output:
       x final coordinates vector

   Return:
       none

   Global variables:
       none

   Specific functions:
       LinearInterpolation2
       Drft
       CopyVec

   Comments:
       Outside of interpolation table simulated by putting 1 in x[4]
       01/07/03 6D tracking activated
       10/01/05 First order kick part added                                  */

  double LN = 0e0;
  T      tx2, tz2;      /* thetax and thetaz retrieved from
			      interpolation routine */
  T      d, B2_perp;
  double alpha0 = 0e0;  // 1/ brh0
  double alpha02= 0e0;  // alpha square
  int    Nslice = 0;
  int    i = 0;
  bool   outoftable = false;

  Nslice = this->N;

  if (this->linear) {
    alpha0 = c0/Lattice.param.Energy*1E-9*this->scaling;
    alpha02 = sgn(this->scaling)*alpha0*alpha0;
  } else
    alpha02 = 1e-6*this->scaling;

//  /* Global -> Local */
//  this->GtoL(X, 0e0, 0e0, 0e0);

  p_rot(this->phi/2e0*180e0/M_PI, x);

  // (Nslice+1) drifts, nslice kicks
  // LN = this->L/(Nslice+1);

  // Nslice drifts and kicks.
  LN = this->L/Nslice;
  Drift(LN/2e0, x);

  for (i = 1; i <= Nslice; i++) {
    // printf("%3d %2d %2d %5.3f %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e\n",
    // 	   i, this->ID->linear, this->ID->secondorder, this->ID->scaling,
    // 	   is_double<T>::cst(x[x_]), is_double<T>::cst(x[px_]),
    // 	   is_double<T>::cst(x[y_]), is_double<T>::cst(x[py_]),
    // 	   is_double<T>::cst(x[delta_]), is_double<T>::cst(x[ct_]));
    // Second order kick
    if (this->secondorder){
      // if (!this->ID->linear)
      //   SplineInterpolation2(x[x_], x[y_], tx2, tz2, outoftable);
      // else {
        LinearInterpolation2(x[x_], x[y_], tx2, tz2, B2_perp, outoftable, 2);

	// Scale locally with (Brho) (as above) instead of when the file
	// is read; since the beam energy might not be known at that time.
	if (Lattice.param.radiation || Lattice.param.emittance)
	  radiate_ID(x, LN, this->scaling*B2_perp);
      // }

      if (outoftable) {
	x[x_] = NAN;
        return;
      }

      d = alpha02/Nslice/(1e0+x[delta_]); x[px_] += d*tx2; x[py_] += d*tz2;
    }
    if (i != Nslice) Drift(LN, x);
  }

  Drift(LN/2e0, x);

  p_rot(this->phi/2e0*180e0/M_PI, x);

//  CopyVec(6L, x, Cell->BeamPos);

//  /* Local -> Global */
//  this->LtoG(X, 0e0, 0e0, 0e0);
}


template<typename T>
void SolenoidType::sol_pass(ss_vect<T> &x)
{
  int          i;
  double       h, z;
  T            hd, AxoBrho, AyoBrho, dAxoBrho[3], dAyoBrho[3], dpy, dpx, B[3];

  h = L/N; z = 0e0;

  for (i = 1; i <= N; i++) {
    hd = h/(1e0+x[delta_]);

    // 1: half step in z
    z += 0.5*h;

    // 2: half drift in y
    AyoBrho = BoBrho*x[x_]/2e0;
    dpx = BoBrho*x[y_]/2e0;
//    get_Axy_EF3(W, z, x, AyoBrho, dAyoBrho, dpx, false);

    x[px_] -= dpx; x[py_] -= AyoBrho;
    x[y_] += 0.5*hd*x[py_];
    x[ct_] += sqr(0.5)*hd*sqr(x[py_])/(1e0+x[delta_]);

    AyoBrho = BoBrho*x[x_]/2e0;
    dpx = BoBrho*x[y_]/2e0;
//    get_Axy_EF3(W, z, x, AyoBrho, dAyoBrho, dpx, false);

    x[px_] += dpx; x[py_] += AyoBrho;

    // 3: full drift in x
    AxoBrho = -BoBrho*x[y_]/2e0;
    dpy = -BoBrho*x[x_]/2e0;
//    get_Axy_EF3(W, z, x, AxoBrho, dAxoBrho, dpy, true);

    x[px_] -= AxoBrho; x[py_] -= dpy; x[x_] += hd*x[px_];
    x[ct_] += 0.5*hd*sqr(x[px_])/(1e0+x[delta_]);

    if (Lattice.param.pathlength) x[ct_] += h;

    AxoBrho = -BoBrho*x[y_]/2e0;
    dpy = -BoBrho*x[x_]/2e0;
//    get_Axy_EF3(W, z, x, AxoBrho, dAxoBrho, dpy, true);

    x[px_] += AxoBrho; x[py_] += dpy;

    // 4: a half drift in y
    AyoBrho = BoBrho*x[x_]/2e0;
    dpx = BoBrho*x[y_]/2e0;
//    get_Axy_EF3(W, z, x, AyoBrho, dAyoBrho, dpx, false);

    x[px_] -= dpx; x[py_] -= AyoBrho;
    x[y_] += 0.5*hd*x[py_];
    x[ct_] += sqr(0.5)*hd*sqr(x[py_])/(1e0+x[delta_]);

    AyoBrho = BoBrho*x[x_]/2e0;
    dpx = BoBrho*x[y_]/2e0;
//    get_Axy_EF3(W, z, x, AyoBrho, dAyoBrho, dpx, false);

    x[px_] += dpx; x[py_] += AyoBrho;

    // 5: half step in z
    z += 0.5*h;

    if (Lattice.param.radiation || Lattice.param.emittance) {
      dAxoBrho[X_] = 0e0;
      dAxoBrho[Y_] = -BoBrho/2e0;
      dAxoBrho[Z_] = 0e0;
      dAyoBrho[X_] = BoBrho/2e0;
      dAyoBrho[Y_] = 0e0;
      dAyoBrho[Z_] = 0e0;
//      get_Axy_EF3(W, z, x, AyoBrho, dAyoBrho, dpx, false);
//      get_Axy_EF3(W, z, x, AxoBrho, dAxoBrho, dpy, true);
      B[X_] = -dAyoBrho[Z_]; B[Y_] = dAxoBrho[Z_];
      B[Z_] = dAyoBrho[X_] - dAxoBrho[Y_];
      radiate(x, h, 0e0, B);
    }
  }
}


template<typename T>
void SolenoidType::Elem_Propagate(ss_vect<T> &ps)
{

  this->GtoL(ps, 0e0, 0e0, 0e0);
  sol_pass(ps);
  this->LtoG(ps, 0e0, 0e0, 0e0);
}


void MapType::Elem_Propagate(ss_vect<double> &ps) { ps = (this->M*ps).cst(); }

void MapType::Elem_Propagate(ss_vect<tps> &ps) { ps = this->M*ps; }


int LatticeType::GetnKid(const int Fnum)
{
  return (Lattice.ElemFam[Fnum-1].nKid);
}


long LatticeType::Elem_GetPos(const int Fnum, const int Knum)
{
  long int loc;

  if (Lattice.ElemFam[Fnum-1].nKid != 0)
    loc = Lattice.ElemFam[Fnum-1].KidList[Knum-1];
  else {
    loc = -1;
    printf("Elem_GetPos: there are no kids in family %d (%s)\n",
	   Fnum, Lattice.ElemFam[Fnum-1].CellF->Name);
    exit_(0);
  }

  return loc;
}


static double thirdroot(double a)
{
  /* By substitution method */
  int i;
  double x;

  x = 1e0; i = 0;
  do {
    i++; x = (x+a)/(x*x+1e0);
  } while (i != 250);
  return x;
}


void SI_init(void)
{
  // SI units are used internally
  // apart from Lattice.param.energy [GeV]
  /*  c_1 = 1/(2*(2-2^(1/3))),    c_2 = (1-2^(1/3))/(2*(2-2^(1/3)))
      d_1 = 1/(2-2^(1/3)),        d_2 = -2^(1/3)/(2-2^(1/3))                 */

  double C_gamma, C_u;

  c_1 = 1e0/(2e0*(2e0-thirdroot(2e0))); c_2 = 0.5e0 - c_1;
  d_1 = 2e0*c_1; d_2 = 1e0 - 2e0*d_1;

  // classical radiation
  // C_gamma = 4*pi*r_e [m]/(3*(m_e [GeV/c^2] *c^2)^3)
  C_gamma = 4e0*M_PI*r_e/(3e0*cube(1e-9*m_e));
  // P_gamma = e^2*c^3/(2*pi)*C_gamma*(E [GeV])^2*(B [T])^2
  // p_s = P_s/P, E = P*c, B/(Brho) = p/e
  cl_rad = C_gamma*cube(Lattice.param.Energy)/(2e0*M_PI);

  // eletron rest mass [GeV]: slightly off???
//  m_e_ = 0.5110034e-03;
  // quantum fluctuations
  C_u = 55e0/(24e0*sqrt(3e0));
  q_fluct =
    3e0*C_u*C_gamma*1e-9*h_bar*c0/(4e0*M_PI*cube(1e-9*m_e))
    *pow(Lattice.param.Energy, 5e0);
}


void Mpole_Print(FILE *f, int Fnum)
{
  ElemFamType *elemfamp;
  MpoleType *M;

  elemfamp = &Lattice.ElemFam[Fnum-1];
  M = static_cast<MpoleType*>(elemfamp->CellF);
  fprintf(f, "Element[%3d ] \n", Fnum);
  fprintf(f, "   Name: %.*s,  Kind:   mpole,  L=% .8E\n",
          SymbolLength, elemfamp->CellF->Name, elemfamp->CellF->L);
  fprintf(f, "   Method: %d, N=%4d\n", M->method, M->N);
}


void Drift_Print(FILE *f, int Fnum)
{
  ElemFamType *elemfamp;

  elemfamp = &Lattice.ElemFam[Fnum-1];
  fprintf(f, "Element[%3d ] \n", Fnum);
  fprintf(f, "   Name: %.*s,  Kind:   drift,  L=% .8E\n",
          SymbolLength, elemfamp->CellF->Name, elemfamp->CellF->L);
  fprintf(f, "   nKid:%3d\n\n", elemfamp->nKid);
}


void Wiggler_Print(FILE *f, int Fnum)
{
  ElemFamType *elemfamp;

  elemfamp = &Lattice.ElemFam[Fnum-1];
  fprintf(f, "Element[%3d ] \n", Fnum);
  fprintf(f, "   Name: %.*s,  Kind:   wiggler,  L=% .8E\n\n",
          NameLength, elemfamp->CellF->Name, elemfamp->CellF->L);
}


void Insertion_Print(FILE *f, int Fnum)
{
  ElemFamType *elemfamp;

  elemfamp = &Lattice.ElemFam[Fnum-1];
  fprintf(f, "Element[%3d ] \n", Fnum);
  fprintf(f, "   Name: %.*s,  Kind:   wiggler,  L=% .8E\n\n",
          SymbolLength, elemfamp->CellF->Name, elemfamp->CellF->L);
}


void Elem_Print(FILE *f, int Fnum)
{
  int i;

  if (Fnum == 0) {
    // print all elements
    for (i = 1; i <= Lattice.param.Elem_nFam; i++)
      Elem_Print(f, i);
    return;
  }

  switch (Lattice.ElemFam[Fnum-1].CellF->Elem.Kind) {
  case PartsKind(drift):
    Drift_Print(f, Fnum);
    break;

  case PartsKind(Mpole):
    Mpole_Print(f, Fnum);
    break;
  case PartsKind(Wigl):
    Wiggler_Print(f, Fnum);
    break;
  case PartsKind(FieldMap):
    break;
  case PartsKind(Insertion):
    Insertion_Print(f, Fnum);
    break;
  case PartsKind(Cavity):
    break;
  case PartsKind(marker):
    break;
  case PartsKind(Spreader):
    break;
  case PartsKind(Recombiner):
    break;
  case PartsKind(Solenoid):
    break;
  case PartsKind(Map):
    break;
  case PartsKind(undef):
    break;
  }
}


void MpoleType::Updateorder(void)
{
  int i;

  if (this->irho != 0e0)
    this->order = 1;
  else
    this->order = 0;
  memcpy(this->B, this->Bpar, sizeof(mpolArray));
  for (i = -HOMmax; i <= HOMmax; i++)
    if (this->B[i+HOMmax] != 0e0 && abs(i) > this->order) this->order = abs(i);
}


// instantiate
template void spline_(const double [], const double [], const int,
		      const double, const double, double []);

template void spline_(const double [], const tps [], const int,
		      const double, const double, tps []);

template void splint_(const double [], const double [], const double [],
		      const int, const double &, double &);

template void splint_(const double [], const tps [], const tps [], const int,
		      const tps &, tps &);

template void splint_(const double [], const double [], const double [],
		      const int, const tps &, tps &);

template void splin2_(const double [], const double [], double **, double **,
		      const int, const int, const double &, const double &,
		      double &);

template void splin2_(const double [], const double [], double **, double **,
		      const int, const int, const tps &, const tps &, tps &);


void FieldMapType::get_B_DIAMOND(const char *filename)
{
  char          line[max_str];
  int           i, j, n1, ny;
  double        x0, y0, z0;
  std::ifstream inf;
  std::ofstream outf;

  const int    skip = 8;
  const double Brho = Lattice.param.Energy*1e9/c0;

  std::cout << std::endl;
  std::cout << "get_B_DIAMOND: loading field map: " << filename << std::endl;

  file_rd(inf, filename);

  for (i = 1; i <= skip; i++)
    inf.getline(line, max_str);

  inf.getline(line, max_str);
  sscanf(line, "Starting point[cm] : %lf %lf %lf", &x0, &y0, &z0);
  inf.getline(line, max_str);
  sscanf(line, "Step size[cm]      : %lf %lf %lf",
	 &dx[X_], &dx[Y_], &dx[Z_]);
  inf.getline(line, max_str);
  sscanf(line, "Number of points   : %d %d %d", &n[X_], &ny, &n[Z_]);

  // Convert from [cm] to [m].
  x0 *= 1e-2; y0 *= 1e-2; z0 *= 1e-2;
  dx[X_] *= 1e-2; dx[Y_] *= 1e-2; dx[Z_] *= 1e-2;
  Lr = dx[Z_]*(n[Z_]-1);

  n[Y_] = 2*ny - 1;

  x[X_] = dvector(1, n[X_]); x[Y_] = dvector(1, n[Y_]);
  x[Z_] = dvector(1, n[Z_]);

  BoBrho[X_]  = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  BoBrho[Y_]  = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  BoBrho[Z_]  = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  BoBrho2[X_] = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  BoBrho2[Y_] = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  BoBrho2[Z_] = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);

  AoBrho[X_]  = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  AoBrho[Y_]  = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  AoBrho2[X_] = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  AoBrho2[Y_] = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);

  for (i = 1; i <= 2; i++)
    inf.getline(line, max_str);

  for (n1 = 1; n1 <= n[Z_]; n1++)
    for (j = 1; j <= ny; j++)
      for (i = 1; i <= n[X_]; i++) {
	inf.getline(line, max_str);
	sscanf(line, "%lf %lf %lf %lf %lf %lf",
	       &x[X_][i], &x[Y_][ny-1+j], &x[Z_][n1],
	       &BoBrho[X_][n1][i][ny-1+j],
	       &BoBrho[Y_][n1][i][ny-1+j],
	       &BoBrho[Z_][n1][i][ny-1+j]);

	// Convert from [cm] to [m].
	x[X_][i] *= 1e-2; x[Y_][ny-1+j] *= 1e-2; x[Z_][n1] *= 1e-2;
	// Convert from [Gauss] to [Tesla].
	BoBrho[X_][n1][i][ny-1+j] /= 1e+4*Brho;
	BoBrho[Y_][n1][i][ny-1+j] /= 1e+4*Brho;
	BoBrho[Z_][n1][i][ny-1+j] /= 1e+4*Brho;
	// Scale.
	BoBrho[X_][n1][i][ny-1+j] *= scl;
	BoBrho[Y_][n1][i][ny-1+j] *= scl;
	BoBrho[Z_][n1][i][ny-1+j] *= scl;

	// Compute vector potential (axial gauge) by extended trapezodial rule.
	if (n1 == 1) {
	  AoBrho[X_][n1][i][ny-1+j] =
	    -BoBrho[Y_][n1][i][ny-1+j]*dx[Z_]/2e0;
	  AoBrho[Y_][n1][i][ny-1+j] =
	    BoBrho[X_][n1][i][ny-1+j]*dx[Z_]/2e0;
	} else if (n1 == n[Z_]) {
	  AoBrho[X_][n1][i][ny-1+j] =
	    AoBrho[X_][n1-1][i][ny-1+j]
	    - BoBrho[Y_][n1][i][ny-1+j]*dx[Z_]/2e0;
	  AoBrho[Y_][n1][i][ny-1+j] =
	    AoBrho[Y_][n1-1][i][ny-1+j]
	    + BoBrho[X_][n1][i][ny-1+j]*dx[Z_]/2e0;
	} else {
	  AoBrho[X_][n1][i][ny-1+j] =
	    AoBrho[X_][n1-1][i][ny-1+j]
	    - BoBrho[Y_][n1][i][ny-1+j]*dx[Z_];
	  AoBrho[Y_][n1][i][ny-1+j] =
	    AoBrho[Y_][n1-1][i][ny-1+j]
	    + BoBrho[X_][n1][i][ny-1+j]*dx[Z_];
	}
      }

  inf.close();

  printf("\n%10.5f %10.5f %10.5f\n", x0, y0, z0);
  printf("%10.5f %10.5f %10.5f\n", dx[X_], dx[Y_], dx[Z_]);
  printf("%10d %10d %10d\n", n[X_], n[Y_], n[Z_]);
  printf("%10.3f -> %10.3f %10.3f -> %10.3f %10.3f -> %10.3f\n",
	 x[X_][1], x[X_][n[X_]],
	 x[Y_][1], x[Y_][n[Y_]],
	 x[Z_][1], x[Z_][n[Z_]]);
  printf("Magnet length [m]: %10.5f\n", Lr);

  for (j = 1; j <= ny-1; j++) {
    x[Y_][j] = -x[Y_][2*ny-j];
    for (i = 1; i <= n[X_]; i++) {
      for (n1 = 1; n1 <= n[Z_]; n1++) {
	// B[X_] is antisymmetric in y (rot(A) = 0)
	BoBrho[X_][n1][i][j] = -BoBrho[X_][n1][i][2*ny-j];
	BoBrho[Y_][n1][i][j] =  BoBrho[Y_][n1][i][2*ny-j];
	// Bz is antisymmetric in y
	BoBrho[Z_][n1][i][j] = -BoBrho[Z_][n1][i][2*ny-j];

	AoBrho[X_][n1][i][j] = AoBrho[X_][n1][i][2*ny-j];
	// Ay is antisymmetric in y
	AoBrho[Y_][n1][i][j] = -AoBrho[Y_][n1][i][2*ny-j];
      }
    }
  }

  if (true) {
    file_wr(outf, "field_map.dat");
    std::cout << std::scientific << std::setprecision(3)
	      << std::setw(11) << x[X_][7] << std::setw(11) << x[Y_][7]
	      << std::endl;
    for (i = 1; i <= n[X_]; i++)
      outf << std::scientific << std::setprecision(3)
	   << std::setw(11) << x[X_][i]
	   << std::setw(11) << BoBrho[Y_][7][i][7] << std::endl;
    outf.close();
  }

  for (n1 = 1; n1 <= n[Z_]; n1++) {
    splie2_(x[X_], x[Y_], BoBrho[X_][n1],
	    n[X_], n[Y_], BoBrho2[X_][n1]);
    splie2_(x[X_], x[Y_], BoBrho[Y_][n1],
	    n[X_], n[Y_], BoBrho2[Y_][n1]);
    splie2_(x[X_], x[Y_], BoBrho[Z_][n1],
	    n[X_], n[Y_], BoBrho2[Z_][n1]);

    splie2_(x[X_], x[Y_], AoBrho[X_][n1],
	    n[X_], n[Y_], AoBrho2[X_][n1]);
    splie2_(x[X_], x[Y_], AoBrho[Y_][n1],
	    n[X_], n[Y_], AoBrho2[Y_][n1]);
  }

  std::cout << "field map loaded: " << filename << std::endl;

/*  free_dvector(x[X_], 1, n[X_]);
  free_dvector(x[Y_], 1, n[Y_]);
  free_dvector(x[Z_], 1, n[Z_]);

  free_df3tensor(BoBrho[X_],  1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(BoBrho[Y_],  1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(BoBrho[Z_],  1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(BoBrho2[X_], 1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(BoBrho2[Y_], 1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(BoBrho2[Z_], 1, n[Z_], 1, n[X_], 1, n[Y_]);

  free_df3tensor(AoBrho[X_],  1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(AoBrho[Y_],  1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(AoBrho2[X_], 1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(AoBrho2[Y_], 1, n[Z_], 1, n[X_], 1, n[Y_]);*/
}


void FieldMapType::get_B_NSLS_II(const char *filename)
{
  char          line[max_str];
  int           i, j, n1;
  double        x_min[3], x_max[3];
  std::ifstream inf;

  const double Brho = Lattice.param.Energy*1e9/c0;

  std::cout << std::endl;
  std::cout << "get_B_NSLS_II: loading field map: " << filename << std::endl;

  file_rd(inf, filename);

  inf.getline(line, max_str);
  inf.getline(line, max_str);
  sscanf(line, "#%lf <= x <= %lf, dx = %lf, nx = %d",
	 &x_min[X_], &x_max[X_], &dx[X_], &n[X_]);
  x_min[X_] *= 1e-2; x_max[X_] *= 1e-2; dx[X_] *= 1e-2;
  inf.getline(line, max_str);
  sscanf(line, "#%lf <= y <= %lf, dy = %lf, ny = %d",
	 &x_min[Y_], &x_max[Y_], &dx[Y_], &n[Y_]);
  x_min[Y_] *= 1e-2; x_max[Y_] *= 1e-2; dx[Y_] *= 1e-2;
  inf.getline(line, max_str);
  sscanf(line, "#%lf <= z <= %lf, dz = %lf, nz = %d",
	 &x_min[Z_], &x_max[Z_], &dx[Z_], &n[Z_]);
  x_min[Z_] *= 1e-2; x_max[Z_] *= 1e-2; dx[Z_] *= 1e-2;

  x[X_] = dvector(1, n[X_]); x[Y_] = dvector(1, n[Y_]);
  x[Z_] = dvector(1, n[Z_]);

  BoBrho[X_]  = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  BoBrho[Y_]  = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  BoBrho[Z_]  = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  BoBrho2[X_] = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  BoBrho2[Y_] = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  BoBrho2[Z_] = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);

  AoBrho[X_]  = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  AoBrho[Y_]  = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  AoBrho2[X_] = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  AoBrho2[Y_] = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);

  for (i = 1; i <= n[X_]; i++)
    for (j = 1; j <= n[Y_]; j++)
      for (n1 = 1; n1 <= n[Z_]; n1++) {
	inf.getline(line, max_str);
	sscanf(line, "%lf %lf %lf %lf %lf %lf",
	       &x[X_][i], &x[Y_][j], &x[Z_][n1],
	       &BoBrho[X_][n1][i][j],
	       &BoBrho[Y_][n1][i][j],
	       &BoBrho[Z_][n1][i][j]);

	// convert from cm to m
	x[X_][i] *= 1e-2; x[Y_][j] *= 1e-2; x[Z_][n1] *= 1e-2;

	BoBrho[X_][n1][i][j] = Brho;
	BoBrho[Y_][n1][i][j] = Brho;
	BoBrho[Z_][n1][i][j] = Brho;

	// Compute vector potential (axial gauge) by extended trapezodial rule
 	if (n1 == 1) {
	  AoBrho[X_][n1][i][j] = -BoBrho[Y_][n1][i][j]*dx[Z_]/2e0;
	  AoBrho[Y_][n1][i][j] =  BoBrho[X_][n1][i][j]*dx[Z_]/2e0;
	} else if (n1 == n[Z_]) {
	  AoBrho[X_][n1][i][j] =
	    AoBrho[X_][n1-1][i][j] - BoBrho[Y_][n1][i][j]*dx[Z_]/2e0;
	  AoBrho[Y_][n1][i][j] =
	    AoBrho[Y_][n1-1][i][j] + BoBrho[X_][n1][i][j]*dx[Z_]/2e0;
	} else {
	  AoBrho[X_][n1][i][j] =
	    AoBrho[X_][n1-1][i][j] - BoBrho[Y_][n1][i][j]*dx[Z_];
	  AoBrho[Y_][n1][i][j] =
	    AoBrho[Y_][n1-1][i][j] + BoBrho[X_][n1][i][j]*dx[Z_];
	}
     }

  inf.close();

  Lr = dx[Z_]*(n[Z_]-1);

  std::cout << std::fixed << std::setprecision(5)
	    << std::setw(10) << 1e3*dx[X_]
	    << std::setw(10) << 1e3*dx[Y_]
	    << std::setw(10) << 1e3*dx[Z_] << std::endl;
  std::cout << std::setw(10) << n[X_] << std::setw(10) << n[Y_]
	    << std::setw(10) << n[Z_] << std::endl;
  std::cout << std::fixed << std::setprecision(3)
	    << std::setw(10) << x[X_][1]
	    << std::setw(10) << x[X_][n[X_]]
	    << std::setw(10) << x[Y_][1]
	    << std::setw(10) << x[Y_][n[Y_]]
	    << std::setw(10) << x[Z_][1]
	    << std::setw(10) << x[Z_][n[Z_]] << std::endl;
  std::cout << std::fixed << std::setprecision(5)
	    << "Magnet length [m]:" << std::setw(10) << Lr << std::endl;

  for (n1 = 1; n1 <= n[Z_]; n1++) {
    splie2_(x[X_], x[Y_], BoBrho[X_][n1],
	    n[X_], n[Y_], BoBrho2[X_][n1]);
    splie2_(x[X_], x[Y_], BoBrho[Y_][n1],
	    n[X_], n[Y_], BoBrho2[Y_][n1]);
    splie2_(x[X_], x[Y_], BoBrho[Z_][n1],
	    n[X_], n[Y_], BoBrho2[Z_][n1]);

    splie2_(x[X_], x[Y_], AoBrho[X_][n1],
	    n[X_], n[Y_], AoBrho2[X_][n1]);
    splie2_(x[X_], x[Y_], AoBrho[Y_][n1],
	    n[X_], n[Y_], AoBrho2[Y_][n1]);
  }

  std::cout << "field map loaded: " << filename << std::endl;

/*  free_dvector(x[X_], 1, n[X_]);
  free_dvector(x[Y_], 1, n[Y_]);
  free_dvector(x[Z_], 1, n[Z_]);

  free_df3tensor(BoBrho[X_],  1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(BoBrho[Y_],  1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(BoBrho[Z_],  1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(BoBrho2[X_], 1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(BoBrho2[Y_], 1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(BoBrho2[Z_], 1, n[Z_], 1, n[X_], 1, n[Y_]);

  free_df3tensor(AoBrho[X_],  1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(AoBrho[Y_],  1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(AoBrho2[X_], 1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(AoBrho2[Y_], 1, n[Z_], 1, n[X_], 1, n[Y_]);*/
}


void FieldMapType::get_B_Oleg1(const char *filename)
{
  char          line[max_str];
  int           i, j, n1;
  double        x_min[3], x_max[3];
  std::ifstream inf;

  const double Brho = Lattice.param.Energy*1e9/c0;

  std::cout << std::endl;
  std::cout << "get_B_Oleg1: loading field map: " << filename << std::endl;

  file_rd(inf, filename);

  inf.getline(line, max_str);
  inf.getline(line, max_str);
  sscanf(line, "#%lf <= x <= %lf, dx = %lf, nx = %d",
	 &x_min[X_], &x_max[X_], &dx[X_], &n[X_]);
  dx[X_] *= 1e-3;
  inf.getline(line, max_str);
  sscanf(line, "#%lf <= y <= %lf, dy = %lf, ny = %d",
	 &x_min[Y_], &x_max[Y_], &dx[Y_], &n[Y_]);
  dx[Y_] *= 1e-3;
  inf.getline(line, max_str);
  sscanf(line, "#%lf <= z <= %lf, dz = %lf, nz = %d",
	 &x_min[Z_], &x_max[Z_], &dx[Z_], &n[Z_]);
  dx[Z_] *= 1e-3;

  x[X_] = dvector(1, n[X_]); x[Y_] = dvector(1, n[Y_]);
  x[Z_] = dvector(1, n[Z_]);

  BoBrho[X_]  = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  BoBrho[Y_]  = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  BoBrho[Z_]  = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  BoBrho2[X_] = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  BoBrho2[Y_] = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  BoBrho2[Z_] = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);

  AoBrho[X_]  = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  AoBrho[Y_]  = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  AoBrho2[X_] = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  AoBrho2[Y_] = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);

  for (n1 = 1; n1 <= n[Z_]; n1++)
    for (j = 1; j <= n[Y_]; j++)
      for (i = 1; i <= n[X_]; i++) {
	inf.getline(line, max_str);
	sscanf(line, "%lf %lf %lf %lf %lf %lf",
	       &x[X_][i], &x[Y_][j], &x[Z_][n1],
	       &BoBrho[X_][n1][i][j],
	       &BoBrho[Y_][n1][i][j],
	       &BoBrho[Z_][n1][i][j]);

	// convert from mm to m
	x[X_][i] *= 1e-3; x[Y_][j] *= 1e-3; x[Z_][n1] *= 1e-3;

	BoBrho[X_][n1][i][j] /= Brho;
	BoBrho[Y_][n1][i][j] /= Brho;
	BoBrho[Z_][n1][i][j] /= Brho;

	// Compute vector potential (axial gauge) by extended trapezodial rule
 	if (n1 == 1) {
	  AoBrho[X_][n1][i][j] = -BoBrho[Y_][n1][i][j]*dx[Z_]/2e0;
	  AoBrho[Y_][n1][i][j] =  BoBrho[X_][n1][i][j]*dx[Z_]/2e0;
	} else if (n1 == n[Z_]) {
	  AoBrho[X_][n1][i][j] =
	    AoBrho[X_][n1-1][i][j] - BoBrho[Y_][n1][i][j]*dx[Z_]/2e0;
	  AoBrho[Y_][n1][i][j] =
	    AoBrho[Y_][n1-1][i][j] + BoBrho[X_][n1][i][j]*dx[Z_]/2e0;
	} else {
	  AoBrho[X_][n1][i][j] =
	    AoBrho[X_][n1-1][i][j] - BoBrho[Y_][n1][i][j]*dx[Z_];
	  AoBrho[Y_][n1][i][j] =
	    AoBrho[Y_][n1-1][i][j] + BoBrho[X_][n1][i][j]*dx[Z_];
	}
     }

  inf.close();

  Lr = dx[Z_]*(n[Z_]-1);

  std::cout << std::fixed << std::setprecision(5)
	    << std::setw(10) << 1e3*dx[X_]
	    << std::setw(10) << 1e3*dx[Y_]
	    << std::setw(10) << 1e3*dx[Z_] << std::endl;
  std::cout << std::setw(10) << n[X_] << std::setw(10) << n[Y_]
	    << std::setw(10) << n[Z_] << std::endl;
  std::cout << std::fixed << std::setprecision(3)
	    << std::setw(10) << x[X_][1]
	    << std::setw(10) << x[X_][n[X_]]
	    << std::setw(10) << x[Y_][1]
	    << std::setw(10) << x[Y_][n[Y_]]
	    << std::setw(10) << x[Z_][1]
	    << std::setw(10) << x[Z_][n[Z_]] << std::endl;
  std::cout << std::fixed << std::setprecision(5)
	    << "Magnet length [m]:" << std::setw(10) << Lr << std::endl;

  for (n1 = 1; n1 <= n[Z_]; n1++) {
    splie2_(x[X_], x[Y_], BoBrho[X_][n1],
	    n[X_], n[Y_], BoBrho2[X_][n1]);
    splie2_(x[X_], x[Y_], BoBrho[Y_][n1],
	    n[X_], n[Y_], BoBrho2[Y_][n1]);
    splie2_(x[X_], x[Y_], BoBrho[Z_][n1],
	    n[X_], n[Y_], BoBrho2[Z_][n1]);

    splie2_(x[X_], x[Y_], AoBrho[X_][n1],
	    n[X_], n[Y_], AoBrho2[X_][n1]);
    splie2_(x[X_], x[Y_], AoBrho[Y_][n1],
	    n[X_], n[Y_], AoBrho2[Y_][n1]);
  }

  std::cout << "field map loaded: " << filename << std::endl;

/*  free_dvector(x[X_], 1, n[X_]);
  free_dvector(x[Y_], 1, n[Y_]);
  free_dvector(x[Z_], 1, n[Z_]);

  free_df3tensor(BoBrho[X_],  1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(BoBrho[Y_],  1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(BoBrho[Z_],  1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(BoBrho2[X_], 1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(BoBrho2[Y_], 1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(BoBrho2[Z_], 1, n[Z_], 1, n[X_], 1, n[Y_]);

  free_df3tensor(AoBrho[X_],  1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(AoBrho[Y_],  1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(AoBrho2[X_], 1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(AoBrho2[Y_], 1, n[Z_], 1, n[X_], 1, n[Y_]);*/
}


void FieldMapType::get_B_Oleg2(const char *filename)
{
  char          line[max_str];
  int           i, j, n1;
  double        x_min[3];
  std::ifstream inf;

  const double  Brho = Lattice.param.Energy*1e9/c0;

  std::cout << std::endl;
  std::cout << "get_B_Oleg2: loading field map: " << filename << std::endl;

  file_rd(inf, filename);

  inf.getline(line, max_str);

  inf.getline(line, max_str); sscanf(line, "#%lf", &x_min[X_]);
  inf.getline(line, max_str); sscanf(line, "#%lf", &dx[X_]);
  inf.getline(line, max_str); sscanf(line, "#%d", &n[X_]);

  inf.getline(line, max_str); sscanf(line, "#%lf", &x_min[Y_]);
  inf.getline(line, max_str); sscanf(line, "#%lf", &dx[Y_]);
  inf.getline(line, max_str); sscanf(line, "#%d", &n[Y_]);

  inf.getline(line, max_str); sscanf(line, "#%lf", &x_min[Z_]);
  inf.getline(line, max_str); sscanf(line, "#%lf", &dx[Z_]);
  inf.getline(line, max_str); sscanf(line, "#%d", &n[Z_]);

  std::cout << std::fixed << std::setprecision(5)
	    << std::setw(10) << 1e3*dx[X_]
	    << std::setw(10) << 1e3*dx[Y_]
	    << std::setw(10) << 1e3*dx[Z_] << std::endl;
  std::cout << std::setw(10) << n[X_] << std::setw(10) << n[Y_]
	    << std::setw(10) << n[Z_] << std::endl;
  std::cout << std::fixed << std::setprecision(3)
	    << std::setw(10) << x_min[X_] << std::setw(10) << x_min[Y_]
	    << std::setw(10) << x_min[Z_] << std::endl;

  x[X_] = dvector(1, n[X_]); x[Y_] = dvector(1, n[Y_]);
  x[Z_] = dvector(1, n[Z_]);

  BoBrho[X_]  = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  BoBrho[Y_]  = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  BoBrho[Z_]  = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  BoBrho2[X_] = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  BoBrho2[Y_] = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  BoBrho2[Z_] = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);

  AoBrho[X_]  = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  AoBrho[Y_]  = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  AoBrho2[X_] = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);
  AoBrho2[Y_] = df3tensor(1, n[Z_], 1, n[X_], 1, n[Y_]);

  for (n1 = 1; n1 <= n[Z_]; n1++) {
    x[Z_][n1] = (n1 == 1)? x_min[Z_] : x[Z_][n1-1] + dx[Z_];

    for (j = 1; j <= n[Y_]; j++) {
      x[Y_][j] = (j == 1)? x_min[Y_] : x[Y_][j-1] + dx[Y_];

      for (i = 1; i <= n[X_]; i++) {
	x[X_][i] = (i == 1)? x_min[X_] : x[X_][i-1] + dx[X_];

	inf.getline(line, max_str);
	sscanf(line, "%lf %lf %lf",
	       &BoBrho[X_][n1][i][j], &BoBrho[Y_][n1][i][j],
	       &BoBrho[Z_][n1][i][j]);

	BoBrho[X_][n1][i][j] /= Brho;
	BoBrho[Y_][n1][i][j] /= Brho;
	BoBrho[Z_][n1][i][j] /= Brho;

	// Compute vector potential (axial gauge) by extended trapezodial rule
 	if (n1 == 1) {
	  AoBrho[X_][n1][i][j] = -BoBrho[Y_][n1][i][j]*dx[Z_]/2e0;
	  AoBrho[Y_][n1][i][j] =  BoBrho[X_][n1][i][j]*dx[Z_]/2e0;
	} else if (n1 == n[Z_]) {
	  AoBrho[X_][n1][i][j] =
	    AoBrho[X_][n1-1][i][j] - BoBrho[Y_][n1][i][j]*dx[Z_]/2e0;
	  AoBrho[Y_][n1][i][j] =
	    AoBrho[Y_][n1-1][i][j] + BoBrho[X_][n1][i][j]*dx[Z_]/2e0;
	} else {
	  AoBrho[X_][n1][i][j] =
	    AoBrho[X_][n1-1][i][j] - BoBrho[Y_][n1][i][j]*dx[Z_];
	  AoBrho[Y_][n1][i][j] =
	    AoBrho[Y_][n1-1][i][j] + BoBrho[X_][n1][i][j]*dx[Z_];
	}
      }
    }
  }

  inf.close();

  Lr = dx[Z_]*(n[Z_]-1);

  std::cout << std::fixed << std::setprecision(5)
	    << std::setw(10) << 1e3*dx[X_]
	    << std::setw(10) << 1e3*dx[Y_]
	    << std::setw(10) << 1e3*dx[Z_] << std::endl;
  std::cout << std::setw(10) << n[X_] << std::setw(10) << n[Y_]
	    << std::setw(10) << n[Z_] << std::endl;
  std::cout << std::fixed << std::setprecision(3)
	    << std::setw(10) << x[X_][1]
	    << std::setw(10) << x[X_][n[X_]]
	    << std::setw(10) << x[Y_][1]
	    << std::setw(10) << x[Y_][n[Y_]]
	    << std::setw(10) << x[Z_][1]
	    << std::setw(10) << x[Z_][n[Z_]] << std::endl;
  std::cout << std::fixed << std::setprecision(5)
	    << "Magnet length [m]:" << std::setw(10) << Lr << std::endl;

  for (n1 = 1; n1 <= n[Z_]; n1++) {
    splie2_(x[X_], x[Y_], BoBrho[X_][n1],
	    n[X_], n[Y_], BoBrho2[X_][n1]);
    splie2_(x[X_], x[Y_], BoBrho[Y_][n1],
	    n[X_], n[Y_], BoBrho2[Y_][n1]);
    splie2_(x[X_], x[Y_], BoBrho[Z_][n1],
	    n[X_], n[Y_], BoBrho2[Z_][n1]);
  }

  std::cout << "field map loaded: " << filename << std::endl;

/*  free_dvector(x[X_], 1, n[X_]);
  free_dvector(x[Y_], 1, n[Y_]);
  free_dvector(x[Z_], 1, n[Z_]);

  free_df3tensor(BoBrho[X_],  1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(BoBrho[Y_],  1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(BoBrho[Z_],  1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(BoBrho2[X_], 1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(BoBrho2[Y_], 1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(BoBrho2[Z_], 1, n[Z_], 1, n[X_], 1, n[Y_]);

  free_df3tensor(AoBrho[X_],  1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(AoBrho[Y_],  1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(AoBrho2[X_], 1, n[Z_], 1, n[X_], 1, n[Y_]);
  free_df3tensor(AoBrho2[Y_], 1, n[Z_], 1, n[X_], 1, n[Y_]);*/
}


void FieldMapType::get_B(const char *filename)
{
  // Do not scale fieldmaps only Hamiltonians, i.e., the kick.  Note that RADIA
  // (2nd order) kick maps are quadratic in the field, and 1st order linear.

  switch (Lattice.param.FieldMap_filetype) {
  case 1:
    get_B_DIAMOND(filename);
    break;
  case 2:
    get_B_NSLS_II(filename);
    break;
  case 3:
    get_B_Oleg1(filename);
    break;
  case 4:
    get_B_Oleg2(filename);
    break;
  default:
    std::cout << "get_B: unknown FieldMap type "
	      << Lattice.param.FieldMap_filetype << std::endl;
    break;
  }
}


double LatticeType::Elem_GetKval(int Fnum, int Knum, int Order)
{
  double      Result = 0e0;
  CellType    *cellp;
  MpoleType   *M;
  WigglerType *W;

  if (Fnum > 0) {
    cellp = Lattice.Cell[Lattice.ElemFam[Fnum-1].KidList[Knum-1]];
    switch (cellp->Elem.Kind) {
    case drift:
      Result = 0e0;
      break;
    case marker:
      Result = 0e0;
      break;
    case Cavity:
      Result = 0e0;
      break;
    case Mpole: /* KL*/
      M = static_cast<MpoleType*>(cellp);
      if (M->thick == thicktype(thick_))
	Result = cellp->L*Lattice.Mpole_GetB(Fnum, Knum, Order);
      else
	Result = Lattice.Mpole_GetB(Fnum, Knum, Order);
      break;
    case Wigl:
      W = static_cast<WigglerType*>(cellp);
      Result =
	cellp->L*sqrt(2e0*W->BW[Order+HOMmax]);
      break;
    case FieldMap:
      Result = 0e0;
      break;
    case Insertion:
      Result = 0e0;
      break;
    case Spreader:
      Result = 0e0;
      break;
    case Recombiner:
      Result = 0e0;
      break;
    case Solenoid:
      Result = 0e0;
      break;
    case Map:
      Result = 0e0;
      break;
    case undef:
      break;
    }
  } else
    Result = 0e0;

  return Result;
}


void Mpole_SetB(int Fnum, int Knum, int Order)
{
  /*  called by Cell_SetdP
      Compute full multipole composent as sum of design, systematic
      and random part
      Compute transport matrix if quadrupole (Order=2)
      Set multipole order to Order if multipole (Order >2)                    */

  CellType  *cellp;  /* pointer on the Cell */
  MpoleType *M;/* Pointer on the Multipole */

  cellp  = Lattice.Cell[Lattice.ElemFam[Fnum-1].KidList[Knum-1]];
  M = static_cast<MpoleType*>(cellp);
  M->B[Order+HOMmax] =
    M->Bpar[Order+HOMmax] + M->Bsys[Order+HOMmax]
    + M->Brms[Order+HOMmax]*M->Brnd[Order+HOMmax];
  if (abs(Order) > M->order && M->B[Order+HOMmax] != 0e0)
    M->order = abs(Order);
}


double LatticeType::Mpole_GetB(int Fnum, int Knum, int Order)
{
  /* Return multipole strength (of order Order) for Knum element of
     family Fnum
       Order =  2 for normal quadrupole
             = -2 for skew quadrupole                                         */

  CellType  *cellp;
  MpoleType *M;

  cellp = Lattice.Cell[Lattice.ElemFam[Fnum-1].KidList[Knum-1]];
  M = static_cast<MpoleType*>(cellp);
  return (M->B[Order+HOMmax]);
}


void Mpole_DefBpar(int Fnum, int Knum, int Order, double Bpar)
{
  CellType  *cellp;
  MpoleType *M;

  cellp = Lattice.Cell[Lattice.ElemFam[Fnum-1].KidList[Knum-1]];
  M = static_cast<MpoleType*>(cellp);
  M->Bpar[Order+HOMmax] = Bpar;
}


void Mpole_DefBsys(int Fnum, int Knum, int Order, double Bsys)
{
  CellType  *cellp;
  MpoleType *M;

  cellp = Lattice.Cell[Lattice.ElemFam[Fnum-1].KidList[Knum-1]];
  M = static_cast<MpoleType*>(cellp);
  M->Bsys[Order+HOMmax] = Bsys;
}


void Mpole_SetdS(int Fnum, int Knum)
{
  int       j;
  CellType  *cellp;
  MpoleType *M;

  cellp = Lattice.Cell[Lattice.ElemFam[Fnum-1].KidList[Knum-1]];
  M = static_cast<MpoleType*>(cellp);
  for (j = 0; j <= 1; j++)
    cellp->dS[j] = M->dSsys[j] + M->dSrms[j]*M->dSrnd[j];
}

void Mpole_SetdR(int Fnum, int Knum)
{
  CellType  *cellp;
  MpoleType *M;

  cellp = Lattice.Cell[Lattice.ElemFam[Fnum-1].KidList[Knum-1]];
  M = static_cast<MpoleType*>(cellp);
  cellp->dR[0] = cos(dtor(M->dRpar + M->dRsys + M->dRrms*M->dRrnd));
  cellp->dR[1] = sin(dtor(M->dRpar + M->dRsys + M->dRrms*M->dRrnd));
  /* Calculate simplified p_rots */
  M->c0 = sin(cellp->L*M->irho/2e0); M->c1 = cos(dtor(M->dRpar))*M->c0;
  M->s1 = sin(dtor(M->dRpar))*M->c0;
}


double Mpole_GetdR(int Fnum, int Knum)
{
  CellType  *cellp;
  MpoleType *M;

  cellp = Lattice.Cell[Lattice.ElemFam[Fnum-1].KidList[Knum-1]];
  M = static_cast<MpoleType*>(cellp);
  return(M->dRpar + M->dRsys + M->dRrms*M->dRrnd);
}


void Mpole_DefdRpar(int Fnum, int Knum, double dRpar)
{
  CellType  *cellp;
  MpoleType *M;

  cellp = Lattice.Cell[Lattice.ElemFam[Fnum-1].KidList[Knum-1]];
  M = static_cast<MpoleType*>(cellp);
  M->dRpar = dRpar;
}


void Mpole_DefdRsys(int Fnum, int Knum, double dRsys)
{
  CellType  *cellp;
  MpoleType *M;

  cellp = Lattice.Cell[Lattice.ElemFam[Fnum-1].KidList[Knum-1]];
  M = static_cast<MpoleType*>(cellp);
  M->dRsys = dRsys;
}


void Wiggler_SetB(int Fnum, int Knum, int Order)
{
  CellType    *cellp;
  WigglerType *W;

  cellp = Lattice.Cell[Lattice.ElemFam[Fnum-1].KidList[Knum-1]];
  W = static_cast<WigglerType*>(cellp);
  if (abs(Order) > W->order) W->order = abs(Order);
}


void Wiggler_SetdS(int Fnum, int Knum)
{
  int         j;
  CellType    *cellp;
  WigglerType *W;

  cellp = Lattice.Cell[Lattice.ElemFam[Fnum-1].KidList[Knum-1]];
  W = static_cast<WigglerType*>(cellp);
  for (j = 0; j <= 1; j++)
    cellp->dS[j] = W->dSsys[j] + W->dSrms[j]*W->dSrnd[j];
}

void Wiggler_SetdR(int Fnum, int Knum)
{
  CellType    *cellp;
  WigglerType *W;

  cellp = Lattice.Cell[Lattice.ElemFam[Fnum-1].KidList[Knum-1]];
  W = static_cast<WigglerType*>(cellp);
  cellp->dR[0] = cos(dtor(W->dRpar+W->dRsys+W->dRrms*W->dRrnd));
  cellp->dR[1] = sin(dtor(W->dRpar+W->dRsys+W->dRrms*W->dRrnd));
}
