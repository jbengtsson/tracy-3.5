/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -

   Element propagators.                                                      */

bool          first_FM = true;
double        C_u, C_gamma, C_q, cl_rad, q_fluct, I[6];
double        c_1, d_1, c_2, d_2;
double        s_FM;
std::ofstream outf_;

// for FieldMap
bool  sympl             = true;
int   FieldMap_filetype = 2;


#if 0

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

#endif


template<typename T>
void GtoL(ss_vect<T> &ps, Vector2 &S, Vector2 &R,
	  const double c0, const double c1, const double s1)
{
  ss_vect<T> ps1;

  // Simplified rotated p_rot: R^-1(theta_des) prot(phi/2) R(theta_des).
  ps[px_] += c1; ps[py_] += s1;
  // Eucluclidian transformation:
  //   first Translate,
  ps[x_] -= S[X_]; ps[y_] -= S[Y_];
  //   then Rotate.
  ps1 = ps;
  ps[x_]  =  R[X_]*ps1[x_]  + R[Y_]*ps1[y_];
  ps[px_] =  R[X_]*ps1[px_] + R[Y_]*ps1[py_];
  ps[y_]  = -R[Y_]*ps1[x_]  + R[X_]*ps1[y_];
  ps[py_] = -R[Y_]*ps1[px_] + R[X_]*ps1[py_] ;
  // Simplified p_rot.
  ps[px_] -= c0;
  // Phase space vector is now in magnet's local coordinates.
}


template<typename T>
void LtoG(ss_vect<T> &ps, Vector2 &S, Vector2 &R,
	  double c0, double c1, double s1)
{
  ss_vect<T> ps1;

  // Reverse of GtoL, with inverted Euclidian.
  // Simplified p_rot.
  ps[px_] -= c0;
  // Inverted Eucluclidian transformation:
  //   first Rotate.
  ps1 = ps;
  ps[x_]  = R[X_]*ps1[x_]  - R[Y_]*ps1[y_];
  ps[px_] = R[X_]*ps1[px_] - R[Y_]*ps1[py_];
  ps[y_]  = R[Y_]*ps1[x_]  + R[X_]*ps1[y_];
  ps[py_] = R[Y_]*ps1[px_] + R[X_]*ps1[py_];
  //   then Translate.
  ps[x_] += S[X_]; ps[y_] += S[Y_];
  // Rotated p_rot.
  ps[px_] += c1; ps[py_] += s1;
}


template<typename T>
inline T get_p_s(ConfigType &conf, const ss_vect<T> &ps)
{
  T p_s, p_s2;

  if (!conf.H_exact)
    // Small angle axproximation.
    p_s = 1e0 + ps[delta_];
  else {
    p_s2 = sqr(1e0+ps[delta_]) - sqr(ps[px_]) - sqr(ps[py_]);
    if (p_s2 >= 0e0)
      p_s = sqrt(p_s2);
    else {
//      printf("get_p_s: *** Speed of light exceeded!\n");
      p_s = NAN;
    }
  }
  return(p_s);
}


// partial template-class specialization
// primary version
template<typename T>
class is_tps { };

// partial specialization
template<>
class is_tps<double> {
 public:
  static inline void get_ps(const ss_vect<double> &x, CellType *Cell)
  { Cell->BeamPos = x; }

  static inline double set_prm(const int k) { return 1e0; }

  static inline double get_curly_H(const ss_vect<tps> &x)
    {
      std::cout << "get_curly_H: operation not defined for double" << std::endl;
      exit_(1);
      return 0e0;
    }

  static inline double get_dI_eta(const ss_vect<tps> &A)
    {
      std::cout << "get_dI_eta: operation not defined for double" << std::endl;
      exit_(1);
      return 0e0;
    }

  static inline void emittance(ConfigType &conf, const double B2,
			       const double u, const double ps0,
			       const ss_vect<double> &xp) { }

  static inline void diff_mat(const double B2, const double u,
			      const double ps0, const ss_vect<double> &xp) { }

};


// partial specialization
template<>
class is_tps<tps> {
 public:
  static inline void get_ps(const ss_vect<tps> &x, CellType *Cell)
  {
    Cell->BeamPos = x.cst(); getlinmat(ss_dim, x, Cell->A);
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

  static inline double get_dI_eta(const ss_vect<tps> &A)
  {
    return A[x_][delta_];
  }

  static inline void emittance(ConfigType &conf, const tps &B2_perp,
			       const tps &ds, const tps &p_s0,
			       const ss_vect<tps> &A)
  {
    // M. Sands "The Physics of Electron Storage Rings" SLAC-121, Eq. (5.20),
    // p. 118:
    //   dN<u^2>/E^2 =
    //     3*C_U*C_gamma*h_bar*c*E_0^5*(1+delta)^4*(B_perp/(Brho))^3
    //     /(4*pi*m_e^3 [eV/c^2])
    // A contains the eigenvectors.
    int          j;
    double       B_66;
    ss_vect<tps> A_inv;

    if (B2_perp > 0e0) {
      B_66 = (q_fluct*pow(B2_perp.cst(), 1.5)*pow(p_s0, 4)*ds).cst();
      A_inv = Inv(A);
      // D_11 = D_22 = curly_H_x,y * B_66 / 2,
      // curly_H_x,y = eta_Fl^2 + etap_Fl^2
      for (j = 0; j < 3; j++)
	conf.D_rad[j] +=
	  (sqr(A_inv[j*2][delta_])+sqr(A_inv[j*2+1][delta_]))*B_66/2e0;
    }
  }

  static inline void diff_mat(const tps &B2_perp, const tps &ds,
			      const tps &p_s0, ss_vect<tps> &x)
  { }

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
void radiate(ConfigType &conf, ss_vect<T> &ps, const double L,
	     const double h_ref, const T B[])
{
  // M. Sands "The Physics of Electron Storage Rings" SLAC-121, p. 98.
  // ddelta/d(ds) = -C_gamma*E_0^3*(1+delta)^2*(B_perp/(Brho))^2/(2*pi)
  T          p_s0, p_s1, ds, B2_perp = 0e0, B2_par = 0e0;
  ss_vect<T> cs;

  // Large ring: x' and y' unchanged.
  p_s0 = get_p_s(conf, ps); cs = ps; cs[px_] /= p_s0; cs[py_] /= p_s0;

  // H = -p_s => ds = H*L.
  ds = (1e0+cs[x_]*h_ref+(sqr(cs[px_])+sqr(cs[py_]))/2e0)*L;
  get_B2(h_ref, B, cs, B2_perp, B2_par);

  if (conf.radiation) {
    ps[delta_] -= cl_rad*sqr(p_s0)*B2_perp*ds;
    p_s1 = get_p_s(conf, ps); ps[px_] = cs[px_]*p_s1; ps[py_] = cs[py_]*p_s1;
  }

  if (conf.emittance) is_tps<T>::emittance(conf, B2_perp, ds, p_s0, cs);
}


template<typename T>
void radiate_ID(ConfigType &conf, ss_vect<T> &ps, const double L,
		const T &B2_perp)
{
  T          p_s0, p_s1, ds;
  ss_vect<T> cs;

  // Large ring: x' and y' unchanged.
  cs = ps; p_s0 = get_p_s(conf, ps); cs[px_] /= p_s0; cs[py_] /= p_s0;

  // H = -p_s => ds = H*L.
  ds = (1e0+(sqr(cs[px_])+sqr(cs[py_]))/2e0)*L;

  if (conf.radiation) {
    ps[delta_] -= cl_rad*sqr(p_s0)*B2_perp*ds;
    p_s1 = get_p_s(conf, ps); ps[px_] = cs[px_]*p_s1; ps[py_] = cs[py_]*p_s1;
  }

  if (conf.emittance) is_tps<T>::emittance(conf, B2_perp, ds, p_s0, cs);
}


template<typename T>
void Drift(ConfigType &conf, const double L, ss_vect<T> &ps)
{
  T u;

  if (!conf.H_exact) {
    // Small angle axproximation.
    u = L/(1e0+ps[delta_]);
    ps[x_]  += u*ps[px_]; ps[y_] += u*ps[py_];
    ps[ct_] += u*(sqr(ps[px_])+sqr(ps[py_]))/(2e0*(1e0+ps[delta_]));
  } else {
    u = L/get_p_s(conf, ps);
    ps[x_]  += u*ps[px_]; ps[y_] += u*ps[py_];
    ps[ct_] += u*(1e0+ps[delta_]) - L;
  }
  if (conf.pathlength) ps[ct_] += L;
}


template<typename T>
void DriftType::Drift_Pass(ConfigType &conf, ss_vect<T> &ps)
{
  Drift(conf, PL, ps);

  if (conf.emittance && !conf.Cavity_on)
    // Needs A^-1.
    curly_dH_x = is_tps<tps>::get_curly_H(ps);
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
    psi = k1*gap*irho*(1e0+sqr(sin(degtorad(phi))))/cos(degtorad(phi))
          *(1e0 - k2*gap*irho*tan(degtorad(phi)));

  return psi;
}


template<typename T>
void thin_kick(ConfigType &conf, const int Order, const MpoleArray &MB,
	       const double L, const double h_bend, const double h_ref,
	       ss_vect<T> &ps)
{
  // The vector potential for the combined-function sector bend is from:
  // C. Iselin "Lie Transformations and Transport Equations for Combined-
  // Function Dipoles" Part. Accel. 17, 143-155 (1985).
  int        j;
  T          BxoBrho, ByoBrho, ByoBrho1, B[3], u, p_s;
  ss_vect<T> ps0;

  if ((h_bend != 0e0) || ((1 <= Order) && (Order <= HOMmax))) {
    ps0 = ps;
    // Compute magnetic field with Horner's rule.
    ByoBrho = MB[Order+HOMmax]; BxoBrho = MB[HOMmax-Order];
    for (j = Order-1; j >= 1; j--) {
      ByoBrho1 = ps0[x_]*ByoBrho - ps0[y_]*BxoBrho + MB[j+HOMmax];
      BxoBrho  = ps0[y_]*ByoBrho + ps0[x_]*BxoBrho + MB[HOMmax-j];
      ByoBrho  = ByoBrho1;
    }

    if (conf.radiation || conf.emittance) {
      B[X_] = BxoBrho; B[Y_] = ByoBrho + h_bend; B[Z_] = 0e0;
      radiate(conf, ps, L, h_ref, B);
    }

    if (h_ref != 0e0) {
      // Sector bend.
      if (true) {
	ps[px_] -= L*(ByoBrho+(h_bend-h_ref)/2e0+h_ref*h_bend*ps0[x_]
		     -h_ref*ps0[delta_]);
	ps[ct_] += L*h_ref*ps0[x_];
      } else {
	// The Hamiltonian is split into: H_d + H_k; with [H_d, H_d] = 0.
	p_s = get_p_s(conf, ps0); u = L*h_ref*ps0[x_]/p_s;
	ps[x_]  += u*ps0[px_];
	ps[y_]  += u*ps0[py_];
	ps[ct_] += u*(1e0+ps0[delta_]);
	// ps[px_] -= L*(h_bend*(1e0+h_ref*ps0[x_])-h_ref*p_s);
	// Field expansion up to sextupole like terms.
	ByoBrho += h_bend - MB[Quad+HOMmax]*h_ref*sqr(ps0[y_])/2e0;
	ps[px_] -= L*((1e0+h_ref*ps0[x_])*ByoBrho-h_ref*p_s);
	ps[py_] += L*(1e0+h_ref*ps0[x_])*BxoBrho;
      }
    } else
      // Cartesian bend.
      ps[px_] -= L*(h_bend+ByoBrho);
    ps[py_] += L*BxoBrho;
  }
}


template<typename T>
void EdgeFocus(ConfigType &conf, const double irho, const double phi,
	       const double gap, ss_vect<T> &ps)
{
  ps[px_] += irho*tan(degtorad(phi))*ps[x_];
  if (!conf.dip_edge_fudge) {
    // warning: => diverging Taylor map (see SSC-141)
    // ps[py_] -=
    //   irho*tan(degtorad(phi)-get_psi(irho, phi, gap))*ps[y_]/(1e0+ps[delta_]);
    // Leading order correction.
    ps[py_] -=
      irho*tan(degtorad(phi)-get_psi(irho, phi, gap))*ps[y_]*(1e0-ps[delta_]);
  } else
    ps[py_] -= irho*tan(degtorad(phi)-get_psi(irho, phi, gap))*ps[y_];
}


template<typename T>
void p_rot(ConfigType &conf, double phi, ss_vect<T> &ps)
{
  T          c, s, t, pz, p, val;
  ss_vect<T> ps1;

  c = cos(degtorad(phi));
  s = sin(degtorad(phi));
  t = tan(degtorad(phi));
  pz = get_p_s(conf, ps);

  if (!conf.H_exact && !conf.Cart_Bend) {
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
void bend_fringe(ConfigType &conf, const double hb, ss_vect<T> &ps)
{
  T          coeff, u, pz, pz2, pz3;
  ss_vect<T> ps1;

  coeff = -hb/2e0;
  ps1 = ps;
  pz = get_p_s(conf, ps);
  pz2 = sqr(pz); pz3 = pz*pz2;
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
void quad_fringe(ConfigType &conf, const double b2, ss_vect<T> &ps)
{
  T u, p_s;

  u = b2/(12e0*(1e0+ps[delta_])); p_s = u/(1e0+ps[delta_]);
  ps[py_] /= 1e0 - 3e0*u*sqr(ps[y_]); ps[y_] -= u*cube(ps[y_]);
  if (conf.Cavity_on) ps[ct_] -= p_s*cube(ps[y_])*ps[py_];
  ps[px_] /= 1e0 + 3e0*u*sqr(ps[x_]);
  if (conf.Cavity_on) ps[ct_] += p_s*cube(ps[x_])*ps[px_];
  ps[x_] += u*cube(ps[x_]); u = u*3e0; p_s = p_s*3e0;
  ps[y_] = exp(-u*sqr(ps[x_]))*ps[y_]; ps[py_] = exp(u*sqr(ps[x_]))*ps[py_];
  ps[px_] += 2e0*u*ps[x_]*ps[y_]*ps[py_];
  if (conf.Cavity_on) ps[ct_] -= p_s*sqr(ps[x_])*ps[y_]*ps[py_];
  ps[x_] = exp(u*sqr(ps[y_]))*ps[x_]; ps[px_] = exp(-u*sqr(ps[y_]))*ps[px_];
  ps[py_] -= 2e0*u*ps[y_]*ps[x_]*ps[px_];
  if (conf.Cavity_on) ps[ct_] += p_s*sqr(ps[y_])*ps[x_]*ps[px_];
}


inline ss_vect<double> mat_pass(arma::mat &M, ss_vect<double> &ps)
{
  ss_vect<double> ps1;
  arma::vec       ps_vec(ss_dim);

  ps1 = is_double< ss_vect<double> >::ps(ps);
  ps_vec = M*pstovec(ps1);
  return vectops(ps_vec);
}


inline ss_vect<tps> mat_pass(arma::mat &M, ss_vect<tps> &ps)
{
  ss_vect<tps> ps1;
  arma::mat    ps_mat(ss_dim, ss_dim);

  ps1 = is_double< ss_vect<tps> >::ps(ps);
  getlinmat(ss_dim, ps1, ps_mat);
  ps_mat = M*ps_mat;
  return putlinmat(ss_dim, ps_mat);
}


template<typename T>
void MpoleType::Mpole_Pass(ConfigType &conf, ss_vect<T> &ps)
{
  int          seg = 0, i;
  double       dL = 0e0, dL1 = 0e0, dL2 = 0e0,
               dkL1 = 0e0, dkL2 = 0e0, h_ref = 0e0;

  GtoL(ps, dS, dT, Pc0, Pc1, Ps1);

  if (conf.emittance && !conf.Cavity_on) {
    // Needs A^-1.
    curly_dH_x = 0e0;
    for (i = 0; i <= 5; i++)
      dI[i] = 0e0;
  }

  switch (Pmethod) {

  case Meth_Fourth:
    if (conf.mat_meth && (Porder <= Quad)) {
      ps = mat_pass(M_lin, ps);

      // if (emittance && !Cavity_on)
      // 	if ((PL != 0e0) && (Pirho != 0e0))
      // 	  get_dI_eta_5(this);
    } else {
      // Fringe fields.
      if (conf.quad_fringe && (PB[Quad+HOMmax] != 0e0))
	quad_fringe(conf, PB[Quad+HOMmax], ps);
      if (!conf.Cart_Bend) {
	if (Pirho != 0e0)
	  EdgeFocus(conf, Pirho, PTx1, Pgap, ps);
      } else {
	p_rot(conf, PTx1, ps); bend_fringe(conf, Pirho, ps);
      }

      if (Pthick == thick) {
	if (!conf.Cart_Bend) {
	  // Polar coordinates.
	  h_ref = Pirho; dL = PL/PN;
	} else {
	  // Cartesian coordinates.
	  h_ref = 0e0;
	  if (Pirho == 0e0)
	    dL = PL/PN;
	  else
	    dL = 2e0/Pirho*sin(PL*Pirho/2e0)/PN;
	}

	dL1 = c_1*dL; dL2 = c_2*dL; dkL1 = d_1*dL; dkL2 = d_2*dL;

	for (seg = 1; seg <= PN; seg++) {
	  if (conf.emittance && !conf.Cavity_on) {
	    // Needs A^-1.
	    curly_dH_x += is_tps<tps>::get_curly_H(ps);
	    dI[4] += is_tps<tps>::get_dI_eta(ps);
	  }

	  Drift(conf, dL1, ps);
	  thin_kick(conf, Porder, PB, dkL1, Pirho, h_ref, ps);
	  Drift(conf, dL2, ps);
	  thin_kick(conf, Porder, PB, dkL2, Pirho, h_ref, ps);

	  if (conf.emittance && !conf.Cavity_on) {
	    // Needs A^-1.
	    curly_dH_x += 4e0*is_tps<tps>::get_curly_H(ps);
	    dI[4] += 4e0*is_tps<tps>::get_dI_eta(ps);
	  }

	  Drift(conf, dL2, ps);
	  thin_kick(conf, Porder, PB, dkL1, Pirho, h_ref, ps);
	  Drift(conf, dL1, ps);

	  if (conf.emittance && !conf.Cavity_on) {
	    // Needs A^-1.
	    curly_dH_x += is_tps<tps>::get_curly_H(ps);
	    dI[4] += is_tps<tps>::get_dI_eta(ps);
	  }
	}

	if (conf.emittance && !conf.Cavity_on) {
	  // Needs A^-1.
	  curly_dH_x /= 6e0*PN;
	  dI[1] += PL*is_tps<tps>::get_dI_eta(ps)*Pirho;
	  dI[2] += PL*sqr(Pirho);
	  dI[3] += PL*fabs(cube(Pirho));
	  dI[4] *=
	    PL*Pirho*(sqr(Pirho)+2e0*PBpar[Quad+HOMmax])
	    /(6e0*PN);
	  dI[5] += PL*fabs(cube(Pirho))*curly_dH_x;
	}
      } else
	thin_kick(conf, Porder, PB, 1e0, 0e0, 0e0, ps);

      // Fringe fields.
      if (!conf.Cart_Bend) {
	if (Pirho != 0e0)
	  EdgeFocus(conf, Pirho, PTx2, Pgap, ps);
      } else {
	bend_fringe(conf, -Pirho, ps); p_rot(conf, PTx2, ps);
      }
      if (conf.quad_fringe && (PB[Quad+HOMmax] != 0e0))
	quad_fringe(conf, -PB[Quad+HOMmax], ps);
    }
    break;

  default:
    printf("Mpole_Pass: Method not supported %10s %d\n",
	   PName, Pmethod);
    exit_(0);
    break;
  }

  LtoG(ps, dS, dT, Pc0, Pc1, Ps1);
}


template<typename T>
void MarkerType::Marker_Pass(ConfigType &conf, ss_vect<T> &ps)
{
  GtoL(ps, dS, dT, 0e0, 0e0, 0e0);

  if (conf.emittance && !conf.Cavity_on)
    // Needs A^-1.
    curly_dH_x = is_tps<tps>::get_curly_H(ps);

  LtoG(ps, dS, dT, 0e0, 0e0, 0e0);
}


template<typename T, typename U>
void Cav_Focus(const double L, const T delta, const bool entrance,
           ss_vect<U> &ps)
{
  double sgn;
 
  sgn = (entrance)? -1e0 : 1e0;

  ps[px_] += sgn*ps[x_]*delta/(2e0*L);
  ps[py_] += sgn*ps[y_]*delta/(2e0*L);
}

#if 1

template<typename T>
void CavityType::Cavity_Pass(ConfigType &conf, ss_vect<T> &ps)
{
  double     L;
  T          delta;

  L = PL;
  Drift(conf, L/2e0, ps);
  if (conf.Cavity_on && Pvolt != 0e0) {
    delta = -Pvolt/(conf.Energy*1e9)
            *sin(2e0*M_PI*Pfreq/c0*ps[ct_]+phi);
    ps[delta_] += delta;

    if (conf.radiation) conf.dE -= is_double<T>::cst(delta);

    if (conf.pathlength) ps[ct_] -= Ph/Pfreq*c0;
  }
  Drift(conf, L/2e0, ps);
}

#else

template<typename T>
void Cav_Pass1(ConfigType &conf, CellType &Cell, ss_vect<T> &ps)
{
  /* J. Rosenzweig and L. Serafini "Transverse Particle Motion in
     Radio-Frequency Linear Accelerators" Phys. Rev. E 49(2),
     1599-1602 (1994).                                                        */

  int        k;
  double     L, h, p_t1;
  T          delta_max, ddelta, delta;

  L = PL;
  h = L/(PN+1e0);
  // Energy contains p_0.
  delta_max = Pvolt/(1e9*Energy); ddelta = delta_max/PN;
  delta = delta_max*sin(2e0*M_PI*Pfreq*ps[ct_]/c0+phi);
  if (entry_focus) Cav_Focus(L, delta, true, ps);
  for (k = 0; k < PN; k++) {
    Drift(conf, h, ps);

    ps[delta_] -= ddelta*sin(2e0*M_PI*Pfreq*(ps[ct_]-k*h)/c0+phi);

    if (conf.radiation) conf.dE -= is_double<T>::cst(ddelta);
    if (conf.pathlength) ps[ct_] -= Ph/Pfreq*c0;
  }
  Drift(conf, h, ps);
  if (exit_focus) Cav_Focus(L, delta, false, ps);

  if (false) {
    // Update p_0.
    p_t1 = is_double<T>::cst(ps[delta_]);
    ps[delta_] -= p_t1;
    // Energy contains p_0.
    Energy *= sqrt(1e0+2e0*p_t1/beta0+sqr(p_t1));
    gamma0 = sqrt(sqr(m_e)+sqr(1e9*Energy))/m_e;
    beta0  = sqrt(1e0-1e0/sqr(gamma0));
    printf("\np0 = %12.5e, beta0 = %12.5e, gamma0 = %12.5e\n",
	   Energy, beta0, gamma0);
  }
}


template<typename T>
void CavityType::Cavity_Pass(ss_vect<T> &ps)
{
  /* J. Rosenzweig and L. Serafini "Transverse Particle Motion in
     Radio-Frequency Linear Accelerators" Phys. Rev. E 49(2),
     1599-1602 (1994).                                                        */

  double     L, Lambda, phi;
  double     dgammaMax, dgamma, gamma, gamma1;
  double     sf1, f2, f2s;
  double     f5, sf5, dpr, dct, p_t1;
  double     p0, p_t, delta, alpha, dp;
  ss_vect<T> ps0;

  const bool RandS = false;
 
  L = PL;
  phi = phi;
  Lambda = c0/Pfreq;

  p_t = is_double<T>::cst(ps[delta_]);
  delta = sqrt(1e0+2e0*p_t/beta0+sqr(p_t)) - 1e0;
  // Energy contains p_0 [GeV].
  p0 = 1e9*Energy/m_e;
  gamma = sqrt(1e0+sqr(p0));
  dgammaMax = Pvolt/m_e; dgamma = dgammaMax*sin(phi);
  gamma1 = gamma + dgamma;
  dp = sqrt(sqr(gamma1)-1e0) - p0;

  printf("delta = %e, p0 = %e, gamma = %e, gamma1 = %e, dp = %e\n",
	 delta, p0, gamma, gamma1, dp);

  if (entry_focus) Cav_Focus(L, dgamma/gamma, true, ps);

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
      (2e0*dgammaMax*M_PI*cos(phi)*f2)/(Lambda*f5)*ps[ct_]
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

   // Energy contains p_0 [GeV].
    ps[delta_] =
      2e0*M_PI*Pfreq*dgammaMax*cos(phi)/(c0*gamma1)*ps0[ct_]
      + 1e0/(1e0+dp/p0)*ps0[delta_];
  }

  if (exit_focus) Cav_Focus(L, dgamma/(gamma+dgamma), false, ps);

  if (false) {
    // Update p_0.
    p_t1 = is_double<T>::cst(ps[delta_]);
    ps[delta_] -= p_t1;
    // Energy contains p_0.
    Energy *= sqrt(1e0+2e0*p_t1/beta0+sqr(p_t1));
    gamma0 = sqrt(sqr(m_e)+sqr(p0))/m_e;
    beta0  = sqrt(1e0-1e0/sqr(gamma0));
    printf("\np0 = %12.5e, beta0 = %12.5e, gamma0 = %12.5e\n",
	   Energy, beta0, gamma0);
  }
}

#endif


void get_dI_eta_5_ID(MpoleType *Cell)
{
  double       L, K, h, b2, alpha, beta, gamma, psi, eta, etap;
  ss_vect<tps> Id;
  CellType     *Cellp;


  Id.identity();

  L = Cell->PL;
  h = Cell->Pirho;
  b2 = Cell->PBpar[Quad+HOMmax];
  K = b2 + sqr(Cell->Pirho);
  psi = sqrt(fabs(K))*L;
  // *** No pointer arithmetic for C++ polymorphic class.
  Cellp = Cell - 1;
  alpha = Cellp->Alpha[X_]; beta = Cellp->Beta[X_];
  gamma = (1e0+sqr(alpha))/beta;
  eta = Cellp->Eta[X_]; etap = Cellp->Etap[X_];

  Cell->dI[1] += L*eta*h;
  Cell->dI[2] += L*sqr(h);
  Cell->dI[3] += L*fabs(cube(h));

  if (K > 0e0) {
    Cell->dI[4] +=
      h/K*(2e0*b2+sqr(h))
      *((eta*sqrt(K)*sin(psi)+etap*(1e0-cos(psi)))
	+ h/sqrt(K)*(psi-sin(psi)));

    Cell->dI[5] +=
      L*fabs(cube(h))
      *(gamma*sqr(eta)+2e0*alpha*eta*etap+beta*sqr(etap))
      - 2e0*pow(h, 4)/(pow(K, 3e0/2e0))
      *(sqrt(K)*(alpha*eta+beta*etap)*(cos(psi)-1e0)
	+(gamma*eta+alpha*etap)*(psi-sin(psi)))
      + fabs(pow(h, 5))/(4e0*pow(K, 5e0/2e0))
      *(2e0*alpha*sqrt(K)*(4e0*cos(psi)-cos(2e0*psi)-3e0)
	+beta*K*(2e0*psi-sin(2e0*psi))
	+gamma*(6e0*psi-8e0*sin(psi)+sin(2e0*psi)));
  } else {
    K = fabs(K);

    Cell->dI[4] +=
      h/K*(2e0*b2+sqr(h))
      *((eta*sqrt(K)*sinh(psi)-etap*(1e0-cosh(psi)))
	- h/sqrt(K)*(psi-sinh(psi)));

    Cell->dI[5] +=
      L*fabs(cube(h))*
      (gamma*sqr(eta)+2e0*alpha*eta*etap+beta*sqr(etap))
      + 2e0*pow(h, 4)/(pow(K, 3e0/2e0))
      *(sqrt(K)*(alpha*eta+beta*etap)*(cosh(psi)-1e0)
	+(gamma*eta+alpha*etap)*(psi-sinh(psi)))
      + fabs(pow(h, 5))/(4e0*pow(K, 5e0/2e0))
      *(2e0*alpha*sqrt(K)*(4e0*cosh(psi)-cosh(2e0*psi)-3e0)
	-beta*K*(2e0*psi-sinh(2e0*psi))
	+gamma*(6e0*psi-8e0*sinh(psi)+sinh(2e0*psi)));
  }
}


template<typename T>
inline void get_Axy(ConfigType &conf, const WigglerType *W,
		    const double z, ss_vect<T> &x, T AxoBrho[], T AyoBrho[])
{
  int     i;
  double  ky, kz_n;
  T       cx, cz, sx, sz, chy, shy;

  for (i = 0; i <= 3; ++i) {
    AxoBrho[i] = 0e0; AyoBrho[i] = 0e0;
  }

  for (i = 0; i < W->n_harm; i ++) {
    kz_n = W->harm[i]*2e0*M_PI/W->Lambda; ky = sqrt(sqr(W->kxV[i])+sqr(kz_n));
    cx = cos(W->kxV[i]*x[x_]); sx = sin(W->kxV[i]*x[x_]);
    chy = cosh(ky*x[y_]); shy = sinh(ky*x[y_]); sz = sin(kz_n*z);

    AxoBrho[0] += W->BoBrhoV[i]/kz_n*cx*chy*sz;
    AyoBrho[0] += W->BoBrhoV[i]*W->kxV[i]/(ky*kz_n)*sx*shy*sz;

    // derivatives with respect to x
    AxoBrho[1] -= W->BoBrhoV[i]*W->kxV[i]/kz_n*sx*chy*sz;
    AyoBrho[1] += W->BoBrhoV[i]*sqr(W->kxV[i])/(ky*kz_n)*cx*shy*sz;

    // derivatives with respect to y
    AxoBrho[2] += W->BoBrhoV[i]*ky/kz_n*cx*shy*sz;
    AyoBrho[2] += W->BoBrhoV[i]*W->kxV[i]/kz_n*sx*chy*sz;

    if (conf.radiation) {
      cz = cos(kz_n*z);
      // derivatives with respect to z
      AxoBrho[3] += W->BoBrhoV[i]*cx*chy*cz;
      AyoBrho[3] += W->BoBrhoV[i]*W->kxV[i]/ky*sx*shy*cz;
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

  if (radiation) {
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
void Wiggler_pass_EF(ConfigType &conf, const ElemType *elem,
		     ss_vect<T> &ps)
{
  // First order symplectic integrator for wiggler using expanded Hamiltonian

  int          i, nstep = 0;
  double       h, z;
  T            AxoBrho[4] = {0e0, 0e0, 0e0, 0e0};
  T            AyoBrho[4] = {0e0, 0e0, 0e0, 0e0};
  T            psi, hodp, a12, a21, a22, det;
  T            d1, d2, a11, c11, c12, c21, c22, x2, B[3];

  switch (elem->Pkind) {
  case Wigl:
    nstep = dynamic_cast<const WigglerType*>(elem)->PN;
    break;
  case FieldMap:
    nstep = dynamic_cast<const FieldMapType*>(elem)->n_step;
    break;
  default:
    std::cout << "Wiggler_pass_EF: unknown element type" << std::endl;
    exit_(1);
    break;
  }

  h = elem->PL/nstep; z = 0e0;
  for (i = 1; i <= nstep; ++i) {
    switch (elem->Pkind) {
    case Wigl:
      get_Axy(conf, dynamic_cast<const WigglerType*>(elem), z, ps, AxoBrho,
	      AyoBrho);
      break;
    // case FieldMap:
    //   get_Axy_map(FM, z, ps, AxoBrho, AyoBrho);
    //   break;
    default:
      std::cout << "Wiggler_pass_EF: unknown element type" << std::endl;
      exit_(1);
      break;
    }

    psi = 1e0 + ps[delta_]; hodp = h/psi;
    a11 = hodp*AxoBrho[1]; a12 = hodp*AyoBrho[1];
    a21 = hodp*AxoBrho[2]; a22 = hodp*AyoBrho[2];
    det = 1e0 - a11 - a22 + a11*a22 - a12*a21;
    d1 = hodp*AxoBrho[0]*AxoBrho[1]; d2 = hodp*AxoBrho[0]*AxoBrho[2];
    c11 = (1e0-a22)/det; c12 = a12/det; c21 = a21/det; c22 = (1e0-a11)/det;
    x2 = c11*(ps[px_]-d1) + c12*(ps[py_]-d2);

    ps[py_] = c21*(ps[px_]-d1) + c22*(ps[py_]-d2); ps[px_] = x2;
    ps[x_] += hodp*(ps[px_]-AxoBrho[0]); ps[y_] += hodp*ps[py_];
    ps[ct_] += h*(sqr((ps[px_]-AxoBrho[0])/psi)
	      + sqr((ps[py_]-AyoBrho[0])/psi))/2e0;

    if (false)
      std::cout << std::scientific << std::setprecision(3)
	   << std::setw(8) << z
	   << std::setw(11) << is_double<T>::cst(ps[x_])
	   << std::setw(11) << is_double<T>::cst(ps[px_])
	   << std::setw(11) << is_double<T>::cst(ps[y_])
	   << std::setw(11) << is_double<T>::cst(ps[py_])
	   << std::endl;

    if (conf.pathlength) ps[ct_] += h;

    if (conf.radiation || conf.emittance) {
      B[X_] = -AyoBrho[3]; B[Y_] = AxoBrho[3]; B[Z_] = AyoBrho[1] - AxoBrho[2];
      radiate(conf, ps, h, 0e0, B);
    }

    z += h;
  }
}


template<typename T>
inline void get_Axy2(ConfigType &conf, const double z,const double kxV,
		     const double kxH, const double kz, const double BoBrhoV,
		     const double BoBrhoH, const double phi, ss_vect<T> &x,
		     T AxoBrho[], T AyoBrho[])
{
  int i;
  T   cx, sx, cz1, cz2, sz1, sz2, chy, shy, kyH, kyV, chx, shx, cy, sy;

  for (i = 0; i <= 3; ++i) {
    AxoBrho[i] = 0e0; AyoBrho[i] = 0e0;
  }

  kyV = sqrt(sqr(kz)+sqr(kxV)); kyH = sqrt(sqr(kz)+sqr(kxH));
  cx = cos(kxV*x[x_]); sx = sin(kxV*x[x_]);
  cy = cos(kxH*x[y_]); sy = sin(kxH*x[y_]);
  chx = cosh(kyH*x[x_]); shx = sinh(kyH*x[x_]);
  chy = cosh(kyV*x[y_]); shy = sinh(kyV*x[y_]);
  sz1 = sin(kz*z); sz2 = sin(kz*z+phi);

  AxoBrho[0] += BoBrhoV/kz*cx*chy*sz1;
  AxoBrho[0] -= BoBrhoH*kxH/(kyH*kz)*shx*sy*sz2;
  AyoBrho[0] += BoBrhoV*kxV/(kyV*kz)*sx*shy*sz1;
  AyoBrho[0] -= BoBrhoH/kz*chx*cy*sz2;

  /* derivatives with respect to x */
  AxoBrho[1] -= BoBrhoV*kxV/kz*sx*chy*sz1;
  AxoBrho[1] -= BoBrhoH*kxH/kz*chx*sy*sz2;
  AyoBrho[1] += BoBrhoV*sqr(kxV)/(kyV*kz)*cx*shy*sz1;
  AyoBrho[1] -= BoBrhoH*kyH/kz*shx*cy*sz2;

  /* derivatives with respect to y */
  AxoBrho[2] += BoBrhoV*kyV/kz*cx*shy*sz1;
  AxoBrho[2] -= BoBrhoH*sqr(kxH)/(kyH*kz)*shx*cy*sz2;
  AyoBrho[2] += BoBrhoV*kxV/kz*sx*chy*sz1;
  AyoBrho[2] += BoBrhoH*kxH/kz*chx*sy*sz2;

  if (conf.radiation) {
    cz1 = cos(kz*z); cz2=cos(kz*z+phi);
    /* derivatives with respect to z */
    AxoBrho[3] += BoBrhoV*cx*chy*cz1;
    AxoBrho[3] -= BoBrhoH*kxH/kyH*shx*sy*cz2;
    AyoBrho[3] += BoBrhoV*kxV/kyV*sx*shy*cz1;
    AyoBrho[3] -= BoBrhoH*chx*cy*cz2;
  }
}


template<typename T>
void Wiggler_pass_EF2(ConfigType &conf, int nstep, double L, double kxV,
		      double kxH, double kz, double BoBrhoV, double BoBrhoH,
		      double phi, ss_vect<T> &x)
{
  // First order symplectic integrator for wiggler using expanded Hamiltonian

  int    i;
  double h, z;
  T      hodp, B[3], px1, px2, px3, py1, py2, AxoBrho[4], AyoBrho[4], psi;
  T      px = 0e0, py = 0e0;

  h = L/nstep; z = 0e0;
  for (i = 1; i <= nstep; ++i) {
    get_Axy2(conf, z, kxV, kxH, kz, BoBrhoV, BoBrhoH, phi, x, AxoBrho, AyoBrho);

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

    if (conf.pathlength) x[ct_] += h;

    if (conf.radiation || conf.emittance) {
      B[X_] = -AyoBrho[3]; B[Y_] = AxoBrho[3]; B[Z_] = AyoBrho[1] - AxoBrho[2];
      radiate(conf, x, h, 0e0, B);
    }

    z += h;
  }

  x[px_] = px; x[py_] = py;
}


template<typename T>
inline void get_Axy_EF3(ConfigType &conf, const WigglerType *W,
			const double z, const ss_vect<T> &ps, T &AoBrho,
			T dAoBrho[], T &dp, const bool hor)
{
  int    i;
  double ky, kz_n;
  T      cx, sx, sz, chy, shy, cz;

  AoBrho = 0e0; dp = 0e0;

  for (i = 0; i < 3; i++)
    dAoBrho[i] = 0e0;

  for (i = 0; i < W->n_harm; i++) {
    kz_n = W->harm[i]*2e0*M_PI/W->Lambda; ky = sqrt(sqr(W->kxV[i])+sqr(kz_n));

    cx  = cos(W->kxV[i]*ps[x_]); sx = sin(W->kxV[i]*ps[x_]);
    chy = cosh(ky*ps[y_]); shy = sinh(ky*ps[y_]); sz = sin(kz_n*z);

    if (hor) {
      // A_x/Brho
      AoBrho += W->BoBrhoV[i]/kz_n*cx*chy*sz;

      if (conf.radiation || (conf.emittance && !conf.Cavity_on)) {
	cz = cos(kz_n*z);
	dAoBrho[X_] -= W->BoBrhoV[i]*W->kxV[i]/kz_n*sx*chy*sz;
	dAoBrho[Y_] += W->BoBrhoV[i]*ky/kz_n*cx*shy*sz;
	dAoBrho[Z_] += W->BoBrhoV[i]*cx*chy*cz;
      }

      // dp_y
      if (W->kxV[i] == 0e0)
	dp += W->BoBrhoV[i]/kz_n*ky*ps[x_]*shy*sz;
      else
	dp += W->BoBrhoV[i]/(W->kxV[i]*kz_n)*ky*sx*shy*sz;
    } else {
      // A_y/Brho
      AoBrho += W->BoBrhoV[i]*W->kxV[i]/(ky*kz_n)*sx*shy*sz;

      if (conf.radiation || (conf.emittance && !conf.Cavity_on)) {
	cz = cos(kz_n*z);
	dAoBrho[X_] +=
	  W->BoBrhoV[i]*sqr(W->kxV[i])/(ky*kz_n)*cx*shy*sz;
	dAoBrho[Y_] += W->BoBrhoV[i]*W->kxV[i]/kz_n*sx*chy*sz;
	dAoBrho[Z_] += W->BoBrhoV[i]*W->kxV[i]/ky*sx*shy*cz;
      }

      // dp_x
      dp += W->BoBrhoV[i]/kz_n*sqr(W->kxV[i]/ky)*cx*chy*sz;
    }
  }
}


template<typename T>
void Wiggler_pass_EF3(ConfigType &conf, ElemType *elem, ss_vect<T> &ps)
{
  /* Symplectic integrator (2nd order) for Insertion Devices based on:

       E. Forest, et al "Explicit Symplectic Integrator for s-dependent
       Static Magnetic Field" Phys. Rev. E 68,  046502 (2003)                 */

  int        i;
  double     h, z, irho, curly_dH_x;
  T          hd, AxoBrho, AyoBrho, dAxoBrho[3], dAyoBrho[3], dpy, dpx, B[3];
  ss_vect<T> ps1;
 
  const WigglerType *W = dynamic_cast<const WigglerType*>(elem);

  h = elem->PL/W->PN; z = 0e0;

  if (conf.emittance && !conf.Cavity_on) {
    // Needs A^-1.
    elem->curly_dH_x = 0e0;
    for (i = 0; i <= 5; i++)
      elem->dI[i] = 0e0;
  }

  for (i = 1; i <= W->PN; i++) {
    hd = h/(1e0+ps[delta_]);

    // 1: half step in z
    z += 0.5*h;

    // 2: half drift in y
    get_Axy_EF3(conf, W, z, ps, AyoBrho, dAyoBrho, dpx, false);

    ps[px_] -= dpx; ps[py_] -= AyoBrho;
    ps[y_] += 0.5*hd*ps[py_];
    ps[ct_] += sqr(0.5)*hd*sqr(ps[py_])/(1e0+ps[delta_]);

    get_Axy_EF3(conf, W, z, ps, AyoBrho, dAyoBrho, dpx, false);

    ps[px_] += dpx; ps[py_] += AyoBrho;

    // 3: full drift in x
    get_Axy_EF3(conf, W, z, ps, AxoBrho, dAxoBrho, dpy, true);

    ps[px_] -= AxoBrho; ps[py_] -= dpy; ps[x_] += hd*ps[px_];
    ps[ct_] += 0.5*hd*sqr(ps[px_])/(1e0+ps[delta_]);

    if (conf.pathlength) ps[ct_] += h;

    get_Axy_EF3(conf, W, z, ps, AxoBrho, dAxoBrho, dpy, true);

    ps[px_] += AxoBrho; ps[py_] += dpy;

    // 4: a half drift in y
    get_Axy_EF3(conf, W, z, ps, AyoBrho, dAyoBrho, dpx, false);

    ps[px_] -= dpx; ps[py_] -= AyoBrho;
    ps[y_] += 0.5*hd*ps[py_];
    ps[ct_] += sqr(0.5)*hd*sqr(ps[py_])/(1e0+ps[delta_]);

    get_Axy_EF3(conf, W, z, ps, AyoBrho, dAyoBrho, dpx, false);

    ps[px_] += dpx; ps[py_] += AyoBrho;

    // 5: half step in z
    z += 0.5*h;

    if (conf.radiation || conf.emittance) {
      get_Axy_EF3(conf, W, z, ps, AyoBrho, dAyoBrho, dpx, false);
      get_Axy_EF3(conf, W, z, ps, AxoBrho, dAxoBrho, dpy, true);
      B[X_] = -dAyoBrho[Z_]; B[Y_] = dAxoBrho[Z_];
      B[Z_] = dAyoBrho[X_] - dAxoBrho[Y_];
      // Tranform from Conjugate to Kinematic Momenta.
      ps[px_] -= AxoBrho; ps[py_] -= AyoBrho;
      radiate(conf, ps, h, 0e0, B);
      // Tranform from Kinematic to Conjugate Momenta.
      ps[px_] += AxoBrho; ps[py_] += AyoBrho;
    }

    if (conf.emittance && !conf.Cavity_on) {
      // Needs A^-1.
      get_Axy_EF3(conf, W, z, ps, AxoBrho, dAxoBrho, dpy, true);
      irho = is_double<T>::cst(dAxoBrho[Z_]);
      // Tranform to Configuration Space.
      ps1 = ps;
      ps1[px_] = (ps[px_]-AxoBrho)/(1e0+ps[delta_]);
      curly_dH_x = is_tps<tps>::get_curly_H(ps1);
      elem->curly_dH_x += curly_dH_x;
      elem->dI[1] += is_tps<tps>::get_dI_eta(ps)*irho;
      elem->dI[2] += sqr(irho);
      elem->dI[3] += fabs(cube(irho));
      elem->dI[4] += is_tps<tps>::get_dI_eta(ps)*cube(irho);
      elem->dI[5] += curly_dH_x*fabs(cube(irho));
    }
  }

  if (conf.emittance && !conf.Cavity_on) {
    // Needs A^-1.
    elem->curly_dH_x /= W->PN;
    for (i = 0; i <= 5; i++)
      elem->dI[i] *= elem->PL/W->PN;
  }
}


template<typename T>
void WigglerType::Wiggler_Pass(ConfigType &conf, ss_vect<T> &ps)
{
  int        seg;
  double     L, L1, L2, K1, K2;
  ss_vect<T> ps1;

  const WigglerType *W = dynamic_cast<const WigglerType*>(this);

  // Global -> Local
  GtoL(ps, dS, dT, 0e0, 0e0, 0e0);
  switch (W->Pmethod) {

  case Meth_Linear:
    std::cout << "Wiggler_Pass: Meth_Linear not supported" << std::endl;
    exit_(1);
    break;

  case Meth_First:
    if ((W->BoBrhoV[0] != 0e0) || (W->BoBrhoH[0] != 0e0)) {
      if (!conf.EPU)
	Wiggler_pass_EF(conf, this, ps);
      else {
	Wiggler_pass_EF2(conf, W->PN, PL, W->kxV[0], W->kxH[0],
		2e0*M_PI/W->Lambda, W->BoBrhoV[0], W->BoBrhoH[0],
		W->phi[0], ps);
      }
    } else
      // drift if field = 0
      Drift(conf, PL, ps);
    break;

  case Meth_Second:
    if ((W->BoBrhoV[0] != 0e0) || (W->BoBrhoH[0] != 0e0)) {
      Wiggler_pass_EF3(conf, this, ps);
    } else
      // drift if field = 0
      Drift(conf, PL, ps);
    break;

  case Meth_Fourth:  /* 4-th order integrator */
    L = PL/W->PN;
    L1 = c_1*L; L2 = c_2*L; K1 = d_1*L; K2 = d_2*L;
    for (seg = 1; seg <= W->PN; seg++) {
      Drift(conf, L1, ps); ps1 = ps;
      thin_kick(conf, W->Porder, W->PBW, K1, 0e0, 0e0, ps1);
      ps[py_] = ps1[py_];
      Drift(conf, L2, ps);
      ps1 = ps;
      thin_kick(conf, W->Porder, W->PBW, K2, 0e0, 0e0, ps1);
      ps[py_] = ps1[py_];
      Drift(conf, L2, ps);
      ps1 = ps;
      thin_kick(conf, W->Porder, W->PBW, K1, 0e0, 0e0, ps1);
      ps[py_] = ps1[py_];
      Drift(conf, L1, ps);
    }
    break;
  }
  // Local -> Global
  LtoG(ps, dS, dT, 0e0, 0e0, 0e0);
}

#undef eps
#undef kx


template<typename T>
inline T get_p_s_ps(const ss_vect<T> &ps, const T &qop_Ax, const T &qop_Ay)
{
  // p_s for Phase Space: [x, p_x, y, p_y, delta, ct], for a Vector Potential
  // with Axial Gauge (A_z = 0); Eq. (53), CERN 88-05).

  return sqrt(sqr(1e0+ps[delta_])-sqr(ps[px_]-qop_Ax)-sqr(ps[py_]-qop_Ay));
}


template<typename T>
inline T get_p_s_cs(const ss_vect<T> &cs)
{
  // p_s for Configuration Space: [x, x', y, y', delta, ct];
  // Eq. (39), CERN 88-05) and:
  //   p_s = gamma*s^dot/p0 = (1+delta)*s^dot/v.

  return (1e0+cs[delta_])/sqrt(1e0+sqr(cs[px_])+sqr(cs[py_]));
}


template<typename T>
void radiate_cs(ConfigType &conf, ss_vect<T> &cs, const double L,
		const T B[])
{
  // M. Sands "The Physics of Electron Storage Rings" SLAC-121, p. 98.
  // ddelta/d(ds) = -C_gamma*E_0^3*(1+delta)^2*(B_perp/(Brho))^2/(2*pi)
  T  p_s, ds, B2_perp = 0e0, B2_par = 0e0;

  // Large ring: x' and y' unchanged, ds = -p_s*L.
  ds = (1e0+(sqr(cs[px_])+sqr(cs[py_]))/2e0)*L;
  get_B2(0e0, B, cs, B2_perp, B2_par);

  p_s = get_p_s_cs(cs);

  if (conf.radiation) cs[delta_] -= cl_rad*sqr(p_s)*B2_perp*ds;

  if (conf.emittance) is_tps<T>::emittance(conf, B2_perp, ds, p_s, cs);
}


#if 0

template<typename T>
bool get_BoBrho(const FieldMapType *FM, const double z, const ss_vect<T> &cs,
		T BoBrho[])
{
  int j, kz;

  const double eps = 1e-5;

  kz = 0;
  for (j = 1; j <= FM->n[Z_]; j++)
    if (fabs(z-FM->x[Z_][j]) < eps) {
      kz = j;
      break;
    }

  if (kz == 0) {
    printf("\nget_BoBrho: undefined z %12.10f\n", z);
    return false;
  }

  splin2_(FM->x[X_], FM->x[Y_], FM->BoBrho[X_][kz], FM->BoBrho2[X_][kz],
	  FM->n[X_], FM->n[Y_], cs[x_], cs[y_], BoBrho[X_]);
  if (BoBrho[X_] == NAN) return false;

  splin2_(FM->x[X_], FM->x[Y_], FM->BoBrho[Y_][kz], FM->BoBrho2[Y_][kz],
	   FM->n[X_], FM->n[Y_], cs[x_], cs[y_], BoBrho[Y_]);
  if (BoBrho[Y_] == NAN) return false;

  splin2_(FM->x[X_], FM->x[Y_], FM->BoBrho[Z_][kz], FM->BoBrho2[Z_][kz],
	   FM->n[X_], FM->n[Y_], cs[x_], cs[y_], BoBrho[Z_]);
  if (BoBrho[Z_] == NAN) return false;

  for (j = 0; j < 3; j++)
    BoBrho[j] *= FM->scl;

  return true;
}


template<typename T>
void rk4_(ConfigType &conf, const ElemType *elem, const ss_vect<T> &y,
	  const ss_vect<T> &dydx, const double x, const double h,
	  ss_vect<T> &cs, const double z,
	  void (*derivs)(ConfigType &conf, const ElemType *,
			 const double, const ss_vect<T> &, ss_vect<T> &))
{
  int        j;
  double     xh, hh, h6;
  T          BoBrho[3];
  ss_vect<T> dym, dyt, yt;

  const FieldMapType *FM = dynamic_cast<const FieldMapType*>(elem);

  hh = h*0.5; h6 = h/6e0;
  xh = x + hh; yt = y + hh*dydx;
  (*derivs)(conf, elem, xh, yt, dyt); yt = y + hh*dyt;
  (*derivs)(conf, elem, xh, yt, dym); yt = y + h*dym; dym += dyt;
  (*derivs)(conf, elem, x+h, yt, dyt);
  cs = y + h6*(dydx+dyt+2e0*dym);

  if (conf.radiation || conf.emittance) {
    if (!get_BoBrho(FM, z, cs, BoBrho)) {
      for (j = 0; j < ss_dim; j++)
	cs[j] = NAN;
      return;
    }

    radiate_cs(conf, cs, h, BoBrho);
  }
}


template<typename T>
void f_FM(ConfigType &conf, const ElemType *elem, const double z,
	  const ss_vect<T> &cs, ss_vect<T> &Dcs)
{
  // Coordinates are: [x, x', y, y', -ct, delta].

  int j;
  T   BoBrho[3], p_s;

  const FieldMapType *FM = dynamic_cast<const FieldMapType*>(elem);

  if (!get_BoBrho(FM, z, cs, BoBrho)) {
    for (j = 0; j < ss_dim; j++)
      Dcs[j] = NAN;
    return;
  }

  p_s = get_p_s_cs(cs);

  Dcs[x_]  = cs[px_];

  Dcs[px_] = -(cs[px_]*cs[py_]*BoBrho[X_]-(1e0+sqr(cs[px_]))*BoBrho[Y_]
             + cs[py_]*BoBrho[Z_])/p_s;

  Dcs[y_]  = cs[py_];

  Dcs[py_] = -((1e0+sqr(cs[py_]))*BoBrho[X_]-cs[px_]*cs[py_]*BoBrho[Y_]
             - cs[px_]*BoBrho[Z_])/p_s;

  Dcs[ct_] = (1e0+cs[delta_])/p_s - ((!conf.pathlength)? 1e0 : 0e0);

  Dcs[delta_] = 0e0;
}


template<typename T>
void FieldMap_pass_RK(ConfigType &conf, ElemType *elem, ss_vect<T> &ps)
{
  int        i;
  double     h, z;
  T          p_s;
  ss_vect<T> Dps;

  const int n_step = 2; // Each step needs: f(z_n), f(z_n+h), f(z_n+2h)

  FieldMapType *FM = dynamic_cast<FieldMapType*>(elem);

  switch (FieldMap_filetype) {
  case 1:
    break;
  case 2 ... 5:
    // Transform to right handed system
    ps[x_] = -ps[x_]; ps[px_] = -ps[px_];
    break;
  default:
    printf("\nFieldMap_pass_RK: unknown Fieldmap type: %d\n",
	   FieldMap_filetype);
    exit(1);
    break;
  }

  // [x, px, y, py, -ct, delta] -> [x, x', y, y', -ct, delta], A_x,y,z = 0.
  p_s = get_p_s(conf, ps); ps[px_] /= p_s; ps[py_] /= p_s;

  h = n_step*FM->dx[Z_]; z = FM->x[Z_][1]; FM->Lr = 0e0;
  if (conf.trace)
    outf_ << std::scientific << std::setprecision(3)
	  << std::setw(5) << 0 << std::setw(11) << s_FM
	  << std::setw(11) << is_double< ss_vect<T> >::cst(ps) << "\n";
  for(i = 1+FM->cut; i < FM->n[Z_]-FM->cut; i += n_step) {
    if (i <= FM->n[Z_]-FM->cut-2) {
      f_FM(conf, elem, z, ps, Dps);

      if (Dps[x_] == NAN) {
	std::cout << "FieldMap_pass_RK: particle lost" << std::endl;
	std::cout << ps;
	return;
      }

      rk4_(conf, elem, ps, Dps, FM->x[Z_][i], h, ps, z, f_FM);

      z += h; FM->Lr += h; s_FM += h;
    } else {
      // Use 2nd order Runge-Kutta (aka Midpoint Method).
      f_FM(conf, elem, z, ps, Dps);

      if (Dps[x_] == NAN) {
	std::cout << "FieldMap_pass_RK: particle lost" << std::endl;
	std::cout << ps;
	return;
      }

      ps += h/2e0*Dps;

      z += h/2e0; FM->Lr += h/2e0; s_FM += h/2e0;
    }

    if (conf.trace)
      outf_ << std::scientific << std::setprecision(3)
	    << std::setw(5) << i << std::setw(11) << s_FM
	    << std::setw(11) << is_double< ss_vect<T> >::cst(ps) << "\n";
  }

  // [x, x', y, y', -ct, delta] -> [x, px, y, py, -ct, delta], A_x,y,z = 0.
  p_s = get_p_s_cs(ps); ps[px_] *= p_s; ps[py_] *= p_s;

  switch (FieldMap_filetype) {
  case 1:
    break;
  case 2 ... 5:
    // Transform back to left handed system.
    ps[x_] = -ps[x_]; ps[px_] = -ps[px_];
    break;
  default:
    printf("\nFieldMap_pass_RK: unknown Fieldmap type: %d\n",
	   FieldMap_filetype);
    exit(1);
    break;
  }
}


template<typename T>
void FieldMap_pass_SI(ConfigType &conf, ElemType *elem, ss_vect<T> &ps)
{
  /* E. Chacon-Golcher, F. Neri "A Symplectic Integrator with Arbitrary
     Vector and Scalar Potentials" Phys. Lett. A 372 p. 4661-4666 (2008).    */

  int          i, j = 0;
  double       h, z;
  T            hd, AoBrho[2], dAoBrho[2], AoBrho_int;
  ss_vect<T>   ps1;

  const int    n_step = 2;
  const double d_diff = 1e0;

  FieldMapType *FM = dynamic_cast<FieldMapType*>(elem);

  switch (FieldMap_filetype) {
  case 1:
    break;
  case 2 ... 5:
    // Transform to right handed system
    ps[x_] = -ps[x_]; ps[px_] = -ps[px_];
    break;
  default:
    printf("\nFieldMap_pass_SI: unknown Fieldmap type: %d\n",
	   FieldMap_filetype);
    exit(1);
    break;
  }

  h = n_step*FM->dx[Z_]; z = 0e0; FM->Lr = 0e0;
  if (conf.trace)
    outf_ << std::scientific << std::setprecision(3)
	  << std::setw(5) << 0 << std::setw(11) << s_FM
	  << std::setw(11) << is_double< ss_vect<T> >::cst(ps) << "\n";
  for (i = 1+FM->cut; i < FM->n[Z_]-FM->cut; i += n_step) {
    hd = h/(1e0+ps[delta_]);

    // 1. Half step in z.
    z += 0.5*h; j = i + 1; s_FM += 0.5*h;

    // 2. Half drift in y.
    ps1 = ps;

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][j], FM->AoBrho2[Y_][j],
	    FM->n[X_], FM->n[Y_], ps[x_], ps[y_], AoBrho[0]);

    if (AoBrho[0] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    ps1[y_] += (ps[py_]-AoBrho[0])*0.5*hd;
    ps1[ct_] += sqr(0.5)*hd*sqr(ps[py_]-AoBrho[0])/(1e0+ps[delta_]);

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][j], FM->AoBrho2[Y_][j],
	    FM->n[X_], FM->n[Y_], ps[x_]+d_diff*FM->dx[X_], ps[y_],
	    dAoBrho[1]);

    if (dAoBrho[1] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][j], FM->AoBrho2[Y_][j],
	    FM->n[X_], FM->n[Y_], ps[x_]-d_diff*FM->dx[X_], ps[y_],
	    dAoBrho[0]);

    if (dAoBrho[0] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    AoBrho_int = (dAoBrho[1]-dAoBrho[0])/(2e0*d_diff*FM->dx[X_]);

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][j], FM->AoBrho2[Y_][j],
	    FM->n[X_], FM->n[Y_], ps1[x_]+d_diff*FM->dx[X_], ps1[y_],
	    dAoBrho[1]);

    if (dAoBrho[1] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][j], FM->AoBrho2[Y_][j],
	    FM->n[X_], FM->n[Y_], ps1[x_]-d_diff*FM->dx[X_], ps1[y_],
	    dAoBrho[0]);

    if (dAoBrho[0] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    AoBrho_int += (dAoBrho[1]-dAoBrho[0])/(2e0*d_diff*FM->dx[X_]);

    // Trapezoidal rule
    ps1[px_] +=
      AoBrho_int*(is_double<T>::cst(ps1[y_])-is_double<T>::cst(ps[y_]))/2e0;

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][j], FM->AoBrho2[Y_][j],
	    FM->n[X_], FM->n[Y_], ps1[x_], ps1[y_], AoBrho[1]);

    if (AoBrho[1] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    ps1[py_] += AoBrho[1] - AoBrho[0];

    ps = ps1;

    // 3. Full drift in x.
    ps1 = ps;

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[X_][j], FM->AoBrho2[X_][j],
	    FM->n[X_], FM->n[Y_], ps[x_], ps[y_], AoBrho[0]);

    if (AoBrho[0] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    ps1[x_] += (ps[px_]-AoBrho[0])*hd;
    ps1[ct_] += 0.5*hd*sqr(ps[px_]-AoBrho[0])/(1e0+ps[delta_]);

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[X_][j], FM->AoBrho2[X_][j],
	    FM->n[X_], FM->n[Y_], ps1[x_], ps1[y_], AoBrho[1]);

    if (AoBrho[1] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    ps1[px_] += AoBrho[1] - AoBrho[0];

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[X_][j], FM->AoBrho2[X_][j],
	    FM->n[X_], FM->n[Y_], ps[x_], ps[y_]+d_diff*FM->dx[Y_],
	    dAoBrho[1]);

    if (dAoBrho[1] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[X_][j], FM->AoBrho2[X_][j],
	    FM->n[X_], FM->n[Y_], ps[x_], ps[y_]-d_diff*FM->dx[Y_],
	    dAoBrho[0]);

    if (dAoBrho[0] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    AoBrho_int = (dAoBrho[1]-dAoBrho[0])/(2e0*d_diff*FM->dx[Y_]);

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[X_][j], FM->AoBrho2[X_][j],
	    FM->n[X_], FM->n[Y_], ps1[x_], ps1[y_]+d_diff*FM->dx[Y_],
	    dAoBrho[1]);

    if (dAoBrho[1] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[X_][j], FM->AoBrho2[X_][j],
	    FM->n[X_], FM->n[Y_], ps1[x_], ps1[y_]-d_diff*FM->dx[Y_],
	    dAoBrho[0]);

    if (dAoBrho[0] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    AoBrho_int += (dAoBrho[1]-dAoBrho[0])/(2e0*d_diff*FM->dx[Y_]);

    // Trapezoidal rule
    ps1[py_] +=
      AoBrho_int*(is_double<T>::cst(ps1[x_])-is_double<T>::cst(ps[x_]))/2e0;

    ps = ps1;

    // 4. Half drift in y.
    ps1 = ps;

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][j], FM->AoBrho2[Y_][j],
	    FM->n[X_], FM->n[Y_], ps[x_], ps[y_], AoBrho[0]);

    if (AoBrho[0] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    ps1[y_] += (ps[py_]-AoBrho[0])*0.5*hd;
    ps1[ct_] += sqr(0.5)*hd*sqr(ps[py_]-AoBrho[0])/(1e0+ps[delta_]);

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][j], FM->AoBrho2[Y_][j],
	    FM->n[X_], FM->n[Y_], ps[x_]+d_diff*FM->dx[X_], ps[y_],
	    dAoBrho[1]);

    if (dAoBrho[1] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][j], FM->AoBrho2[Y_][j],
	    FM->n[X_], FM->n[Y_], ps[x_]-d_diff*FM->dx[X_], ps[y_],
	    dAoBrho[0]);

    if (dAoBrho[0] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    AoBrho_int = (dAoBrho[1]-dAoBrho[0])/(2e0*d_diff*FM->dx[X_]);

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][j], FM->AoBrho2[Y_][j],
	    FM->n[X_], FM->n[Y_], ps1[x_]+d_diff*FM->dx[X_], ps1[y_],
	    dAoBrho[1]);

    if (dAoBrho[1] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][j], FM->AoBrho2[Y_][j],
	    FM->n[X_], FM->n[Y_], ps1[x_]-d_diff*FM->dx[X_], ps1[y_],
	    dAoBrho[0]);

    if (dAoBrho[0] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    AoBrho_int += (dAoBrho[1]-dAoBrho[0])/(2e0*d_diff*FM->dx[X_]);

    // Trapezoidal rule
    ps1[px_] +=
      AoBrho_int*(is_double<T>::cst(ps1[y_])-is_double<T>::cst(ps[y_]))/2e0;

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][j], FM->AoBrho2[Y_][j],
	    FM->n[X_], FM->n[Y_], ps1[x_], ps1[y_], AoBrho[1]);

    if (AoBrho[1] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    ps1[py_] += AoBrho[1] - AoBrho[0];

    ps = ps1;

    // 5. Half step in z.
    z += 0.5*h; j = i + 2; s_FM += 0.5*h;

    if (conf.pathlength) ps[ct_] += h;

    FM->Lr += h;

    if (conf.radiation || conf.emittance) {
//      B[X_] = -AoBrhoy[3]; B[Y_] = AoBrho[X_][3];
//      B[Z_] = AoBrhoy[1] - AoBrho[X_][2];
//      radiate(conf, ps, h, 0e0, B);
    }

    if (conf.trace)
      outf_ << std::scientific << std::setprecision(3)
	    << std::setw(5) << 0 << std::setw(11) << s_FM
	    << std::setw(11) << is_double< ss_vect<T> >::cst(ps) << "\n";
  }

  // Change of gauge
  splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[X_][j], FM->AoBrho2[X_][j],
	  FM->n[X_], FM->n[Y_], ps[x_], ps[y_], AoBrho[0]);

  if (AoBrho[0] == NAN) {
    for (j = 0; j < ss_dim; j++)
      ps[j] = NAN;
    return;
  }

  ps[px_] -= AoBrho[0];

  splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][j], FM->AoBrho2[Y_][j],
	  FM->n[X_], FM->n[Y_], ps[x_], ps[y_], AoBrho[0]);

  if (AoBrho[0] == NAN) {
    for (j = 0; j < ss_dim; j++)
      ps[j] = NAN;
    return;
  }

  ps[py_] -= AoBrho[0];

  switch (FieldMap_filetype) {
  case 1:
    break;
  case 2 ... 5:
    // Transform back to left handed system.
    ps[x_] = -ps[x_]; ps[px_] = -ps[px_];
    break;
  default:
    printf("\nFieldMap_pass_SI: unknown Fieldmap type: %d\n",
	   FieldMap_filetype);
    exit(1);
    break;
  }
}


// Instantiate
template void f_FM(ConfigType &, const ElemType *, const double,
		   const ss_vect<double> &, ss_vect<double> &);
template void f_FM(ConfigType &, const ElemType *, const double,
		   const ss_vect<tps> &, ss_vect<tps> &);
template void rk4_(ConfigType &, const ElemType *,
		   const ss_vect<double> &, const ss_vect<double> &,
		   const double, const double, ss_vect<double> &, const double,
		   void (*derivs)(ConfigType &, const ElemType *,
				  const double, const ss_vect<double> &,
				  ss_vect<double> &));
template void rk4_(ConfigType &, const ElemType *,
		   const ss_vect<tps> &, const ss_vect<tps> &, const double,
		   const double,ss_vect<tps> &, const double,
		   void (*derivs)(ConfigType &, const ElemType *,
				  const double, const ss_vect<tps> &,
				  ss_vect<tps> &));
template void FieldMap_pass_RK(ConfigType &conf, ElemType *,
			       ss_vect<double> &);
template void FieldMap_pass_RK(ConfigType &conf, ElemType *,
			       ss_vect<tps> &);
template void FieldMap_pass_SI(ConfigType &conf, ElemType *, ss_vect<double> &);
template void FieldMap_pass_SI(ConfigType &conf, ElemType *, ss_vect<tps> &);

#endif


template<typename T>
void FieldMapType::FieldMap_Pass(ConfigType &conf, ss_vect<T> &ps)
{
  int          k;
  double       Ld;

  const FieldMapType *FM = dynamic_cast<const FieldMapType*>(this);

  if (conf.trace & first_FM) {
    file_wr(outf_, "FieldMap_pass.dat");
    s_FM = 0e0;
    first_FM = false;
  }

  Ld = (FM->Lr-PL)/2e0;
  p_rot(conf, radtodeg(FM->phi/2e0), ps);
  printf("\nFieldMap_Pass:\n");
  printf("  phi = %12.5e\n  cut = %12d\n", FM->phi, FM->cut);
  printf("  entrance negative drift [m] %12.5e\n", -Ld);
  Drift(conf, -Ld, ps);

  // n_step: number of Field Map repetitions.
  for (k = 1; k <= FM->n_step; k++) {
    // if (sympl)
    //   FieldMap_pass_SI(conf, this, ps);
    // else
    //   FieldMap_pass_RK(conf, this, ps);
  }

  printf("  exit negative drift [m]     %12.5e\n", -Ld);
  Drift(conf, -Ld, ps);
  p_rot(conf, radtodeg(FM->phi/2e0), ps);

//  outf_.close();
}


template<typename T>
void InsertionType::Insertion_Pass(ConfigType &conf, ss_vect<T> &x)
{
  double        LN = 0e0;
  T             tx2, tz2;      /* thetax and thetaz retrieved from
				  interpolation routine */
  T             d, B2_perp;
  double        alpha0 = 0e0;  // 1/ brh0
  double        alpha02= 0e0;  // alpha square
  int           Nslice = 0;
  int           i = 0;
  bool          outoftable = false;

  const InsertionType *ID = dynamic_cast<const InsertionType*>(this);

  Nslice = ID->PN;

  if (ID->linear) {
    alpha0 = c0/conf.Energy*1E-9*ID->scaling;
    alpha02 = sgn(ID->scaling)*alpha0*alpha0;
  } else
    alpha02 = 1e-6*ID->scaling;

  p_rot(conf, radtodeg(ID->phi/2e0), x);

  // (Nslice+1) drifts, nslice kicks
  // LN = PL/(Nslice+1);

  // Nslice drifts and kicks.
  LN = PL/Nslice;
  Drift(conf, LN/2e0, x);

  for (i = 1; i <= Nslice; i++) {
    // printf("%3d %2d %2d %5.3f %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e\n",
    // 	   i, ID->linear, ID->secondorder, ID->scaling,
    // 	   is_double<T>::cst(x[x_]), is_double<T>::cst(x[px_]),
    // 	   is_double<T>::cst(x[y_]), is_double<T>::cst(x[py_]),
    // 	   is_double<T>::cst(x[delta_]), is_double<T>::cst(x[ct_]));
    // Second order kick
    if (ID->secondorder){
      // if (!ID->linear)
      //   SplineInterpolation2(x[x_], x[y_], tx2, tz2, *this, outoftable);
      // else {
        LinearInterpolation2(x[x_], x[y_], tx2, tz2, B2_perp, this,
			     outoftable, 2);

	// Scale locally with (Brho) (as above) instead of when the file
	// is read; since the beam energy might not be known at that time.
	if (conf.radiation || conf.emittance)
	  radiate_ID(conf, x, LN, ID->scaling*B2_perp);
      // }

      if (outoftable) {
	x[x_] = NAN;
        return;
      }

      d = alpha02/Nslice/(1e0+x[delta_]); x[px_] += d*tx2; x[py_] += d*tz2;
    }
    if (i != Nslice) Drift(conf, LN, x);
  }

  Drift(conf, LN/2e0, x);

  p_rot(conf, radtodeg(ID->phi/2e0), x);

//  CopyVec(6L, x, BeamPos);
}

template<typename T>
void SpreaderType::Spreader_Pass(ConfigType &conf, ss_vect<T> &ps) { }

template<typename T>
void RecombinerType::Recombiner_Pass(ConfigType &conf, ss_vect<T> &ps) { }

template<typename T>
void sol_pass(ConfigType &conf, const ElemType *elem, ss_vect<T> &ps)
{
  int          i;
  double       h, z;
  T            hd, AxoBrho, AyoBrho, dAxoBrho[3], dAyoBrho[3], dpy, dpx, B[3];

  const SolenoidType *Sol = dynamic_cast<const SolenoidType*>(elem);

  h = elem->PL/Sol->N; z = 0e0;

  for (i = 1; i <= Sol->N; i++) {
    hd = h/(1e0+ps[delta_]);

    // 1: half step in z
    z += 0.5*h;

    // 2: half drift in y
    AyoBrho = Sol->BoBrho*ps[x_]/2e0; dpx = Sol->BoBrho*ps[y_]/2e0;
//    get_Axy_EF3(elem->W, z, x, AyoBrho, dAyoBrho, dpx, false);

    ps[px_] -= dpx; ps[py_] -= AyoBrho;
    ps[y_] += 0.5*hd*ps[py_];
    ps[ct_] += sqr(0.5)*hd*sqr(ps[py_])/(1e0+ps[delta_]);

    AyoBrho = Sol->BoBrho*ps[x_]/2e0; dpx = Sol->BoBrho*ps[y_]/2e0;
//    get_Axy_EF3(elem->W, z, x, AyoBrho, dAyoBrho, dpx, false);

    ps[px_] += dpx; ps[py_] += AyoBrho;

    // 3: full drift in x
    AxoBrho = -Sol->BoBrho*ps[y_]/2e0; dpy = -Sol->BoBrho*ps[x_]/2e0;
//    get_Axy_EF3(elem->W, z, x, AxoBrho, dAxoBrho, dpy, true);

    ps[px_] -= AxoBrho; ps[py_] -= dpy; ps[x_] += hd*ps[px_];
    ps[ct_] += 0.5*hd*sqr(ps[px_])/(1e0+ps[delta_]);

    if (conf.pathlength) ps[ct_] += h;

    AxoBrho = -Sol->BoBrho*ps[y_]/2e0; dpy = -Sol->BoBrho*ps[x_]/2e0;
//    get_Axy_EF3(elem->W, z, x, AxoBrho, dAxoBrho, dpy, true);

    ps[px_] += AxoBrho; ps[py_] += dpy;

    // 4: a half drift in y
    AyoBrho = Sol->BoBrho*ps[x_]/2e0; dpx = Sol->BoBrho*ps[y_]/2e0;
//    get_Axy_EF3(elem->W, z, x, AyoBrho, dAyoBrho, dpx, false);

    ps[px_] -= dpx; ps[py_] -= AyoBrho;
    ps[y_] += 0.5*hd*ps[py_];
    ps[ct_] += sqr(0.5)*hd*sqr(ps[py_])/(1e0+ps[delta_]);

    AyoBrho = Sol->BoBrho*ps[x_]/2e0; dpx = Sol->BoBrho*ps[y_]/2e0;
//    get_Axy_EF3(elem->W, z, x, AyoBrho, dAyoBrho, dpx, false);

    ps[px_] += dpx; ps[py_] += AyoBrho;

    // 5: half step in z
    z += 0.5*h;

    if (conf.radiation || conf.emittance) {
      dAxoBrho[X_] = 0e0;
      dAxoBrho[Y_] = -Sol->BoBrho/2e0;
      dAxoBrho[Z_] = 0e0;
      dAyoBrho[X_] = Sol->BoBrho/2e0;
      dAyoBrho[Y_] = 0e0;
      dAyoBrho[Z_] = 0e0;
//      get_Axy_EF3(elem->W, z, x, AyoBrho, dAyoBrho, dpx, false);
//      get_Axy_EF3(elem->W, z, x, AxoBrho, dAxoBrho, dpy, true);
      B[X_] = -dAyoBrho[Z_]; B[Y_] = dAxoBrho[Z_];
      B[Z_] = dAyoBrho[X_] - dAxoBrho[Y_];
      radiate(conf, ps, h, 0e0, B);
    }
  }
}


template<typename T>
void SolenoidType::Solenoid_Pass(ConfigType &conf, ss_vect<T> &ps)
{
  GtoL(ps, dS, dT, 0e0, 0e0, 0e0);

  sol_pass(conf, this, ps);

  LtoG(ps, dS, dT, 0e0, 0e0, 0e0);
}


template<typename T>
void MapType::Map_Pass(ConfigType &conf, ss_vect<T> &ps)
{ ps = is_double<ss_vect<T>>::ps(M*ps); }


void LatticeType::getelem(long i, ElemType *cellrec) { cellrec = elems[i]; }

void LatticeType::putelem(long i, ElemType *cellrec) { elems[i] = cellrec; }

int LatticeType::GetnKid(const int Fnum) { return (elemf[Fnum-1].nKid); }


inline long LatticeType::Elem_GetPos(const int Fnum, const int Knum)
{
  if (elemf[Fnum-1].nKid > 0)
    return elemf[Fnum-1].KidList[Knum-1];
  else {
    printf("Elem_GetPos: there are no kids in family %d %s (%d)\n",
	   Fnum, elemf[Fnum-1].ElemF->PName, elemf[Fnum-1].nKid);
    exit_(1);
    return -1;
  }
}


static double thirdroot(double a)
{
  /* By substitution method */
  int    i;
  double x;

  x = 1e0; i = 0;
  do {
    i++; x = (x+a)/(x*x+1e0);
  } while (i != 250);
  return x;
}


void LatticeType::SI_init()
{
  // SI units are used internally
  // apart from energy [GeV]
  /*  c_1 = 1/(2*(2-2^(1/3))),    c_2 = (1-2^(1/3))/(2*(2-2^(1/3)))
      d_1 = 1/(2-2^(1/3)),        d_2 = -2^(1/3)/(2-2^(1/3))                 */

  c_1 = 1e0/(2e0*(2e0-thirdroot(2e0))); c_2 = 0.5e0 - c_1;
  d_1 = 2e0*c_1; d_2 = 1e0 - 2e0*d_1;

  // classical radiation
  C_u = 55e0/(24e0*sqrt(3e0));
  // C_gamma = 4*pi*r_e [m]/(3*(m_e [GeV/c^2] *c^2)^3)
  C_gamma = 4e0*M_PI*r_e/(3e0*cube(1e-9*m_e));
  // P_gamma = e^2*c^3/(2*pi)*C_gamma*(E [GeV])^2*(B [T])^2
  // p_s = P_s/P, E = P*c, B/(Brho) = p/e
  cl_rad = C_gamma*cube(conf.Energy)/(2e0*M_PI);

  // eletron rest mass [GeV]: slightly off???
//  m_e_ = 0.5110034e-03;
  // quantum fluctuations
  C_q = 3e0*C_u*h_bar*c0/(4e0*m_e);
  q_fluct = C_q*C_gamma/(M_PI*sqr(1e-9*m_e))*pow(conf.Energy, 5e0);
}


void prt_name(const char *name, const int len)
{
  int j, k;

  j = 0;
  do {
    printf("%c", name[j]);
    j++;
  } while (name[j] != ' ');
  for (k = j; k < len; k++)
    printf(" ");
}


void prt_elem(ElemType *elem)
{
  printf("  ");
  prt_name(elem->PName, 7);
  printf(" %6.3f %2d %6.3f %2d %2d",
	 elem->S, elem->Pkind, elem->PL, elem->Fnum, elem->Knum);
}


void DriftType::print(void)
{
  prt_elem(this);
  printf(" drift \n");
}


double get_phi(MpoleType *elem)
{
  return (elem->Pirho != 0e0)? radtodeg(elem->Pirho*elem->PL) : 0e0;
}

void MpoleType::print(void)
{
  prt_elem(this);
  printf(" mpole  method %1d n_step %2d n_max %2d, n_design %2d reverse %1d"
	 " phi %10.3e phi_1 %10.3e phi_2 %10.3e b_2 %10.3e b_3 %10.3e\n",
	 Pmethod, PN, Porder, n_design, Reverse,
	 get_phi(this), PTx1, PTx2, PB[Quad+HOMmax], PB[Sext+HOMmax]);
}


void CavityType::print(void)
{
  prt_elem(this);
  printf(" cavity\n");
}


void MarkerType::print(void)
{
  prt_elem(this);
  printf(" marker\n");
}


void WigglerType::print(void)
{
  prt_elem(this);
  printf(" wiggl \n");
}


void InsertionType::print(void)
{
  prt_elem(this);
  printf(" ID    \n");
}


void FieldMapType::print(void)
{
  prt_elem(this);
  printf(" FM    \n");
}


void SpreaderType::print(void)
{
  prt_elem(this);
  printf(" spread\n");
}


void RecombinerType::print(void)
{
  prt_elem(this);
  printf(" recomb\n");
}


void SolenoidType::print(void)
{
  prt_elem(this);
  printf(" sol   \n");
}


void MapType::print(void)
{
  prt_elem(this);
  printf(" map   \n");
}


void LatticeType::prt_fam(void)
{
  int k;

  printf("\nFamilies:\n");
  for (k = 0; k < elemf.size(); k++) {
    printf("  %3d ", k+1);
    prt_name(elemf[k].ElemF->PName, 7);
    printf(" %2d %6.3f\n", elemf[k].ElemF->Pkind, elemf[k].ElemF->PL);
  }
}


void LatticeType::prt_elem(void)
{
  int k;

  printf("\nLattice: %lu\n", elems.size());
  for (k = 0; k < elems.size(); k++) {
    printf("  %3d ", k);
    elems[k]->print();
  }
}


void Init_Euclid(ElemType *elem)
{
  elem->dT[X_] = 1e0; // cos = 1.
  elem->dT[Y_] = 0e0; // sin = 0.
  elem->dS[X_] = 0e0; // no hor displacement.
  elem->dS[Y_] = 0e0; // no ver displacement.
}


DriftType* Drift_Alloc(void)
{
  DriftType *D = new DriftType;
  Init_Euclid(D);
  return D;
}


MpoleType* Mpole_Alloc(void)
{
  int       j;
  MpoleType *M = new MpoleType;

  Init_Euclid(M);
  for (j = -HOMmax; j <= HOMmax; j++) {
    M->PBpar.push_back(0e0);
    M->PBsys.push_back(0e0);
    M->PBrms.push_back(0e0);
    M->PBrnd.push_back(0e0);
    M->PB.push_back(0e0);
  }
  
  M->Pmethod = Meth_Fourth; M->PN = 0;
  /* Displacement errors */
  for (j = 0; j <= 1; j++) {
    M->PdSsys[j] = 0e0; M->PdSrms[j] = 0e0; M->PdSrnd[j] = 0e0;
  }
  M->PdTpar = 0e0; /* Roll angle */
  M->PdTsys = 0e0; /* systematic Roll errors */
  M->PdTrms = 0e0; /* random Roll errors */
  M->PdTrnd = 0e0; /* random seed */
  for (j = -HOMmax; j <= HOMmax; j++) {
    /* Initializes multipoles strengths to zero */
    M->PB[j+HOMmax]    = 0e0; M->PBpar[j+HOMmax] = 0e0;
    M->PBsys[j+HOMmax] = 0e0; M->PBrms[j+HOMmax] = 0e0;
    M->PBrnd[j+HOMmax] = 0e0;
  }
  M->Porder = 0; M->n_design = 0;
  M->Pirho  = 0e0; /* inverse of curvature radius */
  M->PTx1   = 0e0; /* Entrance angle */
  M->PTx2   = 0e0; /* Exit angle */
  M->Pgap   = 0e0; /* Gap for fringe field ??? */

  M->Pc0 = 0e0; M->Pc1 = 0e0; M->Ps1 = 0e0;

  // M_lin is allocated in Mpole_Init.

  return M;
}


CavityType* Cavity_Alloc(void)
{
  CavityType *C = new CavityType;

  Init_Euclid(C);
  C->Pvolt = 0e0; C->Pfreq = 0e0; C->phi = 0e0; C->Ph = 0;
  C->entry_focus = false; C->exit_focus = false;

  return C;
}


MarkerType* Marker_Alloc(void)
{
  MarkerType *Mk = new MarkerType;

  Init_Euclid(Mk);
  return Mk;
}


WigglerType* Wiggler_Alloc(void)
{
  int         j;
  WigglerType *W = new WigglerType;

  Init_Euclid(W);
  W->Pmethod = Meth_Linear; W->PN = 0;
  for (j = 0; j <= 1; j++) {
    W->PdSsys[j] = 0e0; W->PdSrnd[j] = 0e0;
  }
  W->PdTpar = 0e0; W->PdTsys = 0e0; W->PdTrnd = 0e0;
  W->n_harm = 0;
  // 2/21/12 J.B. & J.C.
  W->Lambda = 0e0;
  for (j = 0; j < n_harm_max; j++) {
    W->BoBrhoV[j] = 0e0; W->BoBrhoH[j] = 0e0; W->kxV[j] = 0e0; W->kxH[j] = 0e0;
    W->phi[j] = 0e0;
  }
  for (j = 0; j <= HOMmax; j++)
    W->PBW[j+HOMmax] = 0e0;
  W->Porder = 0;

  return W;
}


InsertionType* Insertion_Alloc(void)
{
  int           i = 0, j = 0;
  InsertionType *ID = new InsertionType;

  Init_Euclid(ID);
  ID->Pmethod = Meth_Linear; ID->PN = 0;
  ID->nx = 0; ID->nz = 0;

  /* Initialisation thetax and thetaz to 0*/

  // first order kick map
  if (ID->firstorder){
    for (i = 0; i < IDZMAX; i++){
      for (j = 0; j < IDXMAX; j++) {
	ID->thetax1[i][j] = 0e0; ID->thetaz1[i][j] = 0e0; ID->B2[i][j] = 0e0;
      }
    }
  }

  // second order kick map
  if (ID->secondorder) {
    for (i = 0; i < IDZMAX; i++) {
      for (j = 0; j < IDXMAX; j++) {
          ID->thetax[i][j] = 0e0; ID->thetaz[i][j] = 0e0; ID->B2[i][j] = 0e0;
      }
    }
  }

  // stuffs for interpolation
  for (j = 0; j < IDXMAX; j++)
    ID->tabx[j] = 0e0;

  for (j = 0; j < IDZMAX; j++)
    ID->tabz[j] = 0e0;

  // filenames
  strcpy(ID->fname1,""); strcpy(ID->fname2,"");

//  ID->kx = 0e0;
  for (j = 0; j <= 1; j++) {
    ID->PdSsys[j] = 0e0; ID->PdSrnd[j] = 0e0;
  }
  ID->PdTpar = 0e0; ID->PdTsys = 0e0; ID->PdTrnd = 0e0;
//  for (j = 0; j <= HOMmax; j++)
//    ID->PBW[j+HOMmax] = 0e0;
  ID->Porder = 0;

  return ID;
}


FieldMapType* FieldMap_Alloc(void)
{
  FieldMapType *FM = new FieldMapType;

  Init_Euclid(FM);
  FM->n_step = 0; FM->n[X_] = 0; FM->n[Y_] = 0; FM->n[Z_] = 0; FM->scl = 1e0;
  FM->phi = 0e0; FM->Ld = 0e0; FM->L1 = 0e0; FM->cut = 0; FM->x0 = 0e0;

  return FM;
}


SpreaderType* Spreader_Alloc(void)
{
  int          k;
  SpreaderType *Spr = new SpreaderType;

  Init_Euclid(Spr);
  for (k = 0; k < Spreader_max; k++)
    Spr->Cell_ptrs[k] = NULL;

  return Spr;
}


RecombinerType* Recombiner_Alloc(void)
{
  RecombinerType *Rec = new RecombinerType;

  Init_Euclid(Rec);
  return Rec;
}


SolenoidType* Solenoid_Alloc(void)
{
  int          j;
  SolenoidType *Sol = new SolenoidType;

  Init_Euclid(Sol);
  Sol->N = 0;
  for (j = 0; j <= 1; j++) {
    Sol->PdSsys[j] = 0e0; Sol->PdSrms[j] = 0e0; Sol->PdSrnd[j] = 0e0;
  }
  Sol->dTpar = 0e0; Sol->dTsys = 0e0; Sol->dTrnd = 0e0;
  return Sol;
}


MapType* Map_Alloc(void)
{
  MapType *Map = new MapType;

  Init_Euclid(Map);
  return Map;
}


ElemType* DriftType::Elem_Init(const ConfigType &conf, const bool reverse)
{
  DriftType *Dp;

  const DriftType* D = dynamic_cast<const DriftType*>(this);

  Dp = Drift_Alloc();
  *Dp = *D;
  return Dp;
}


static int UpdatePorder(MpoleType *M)
{
  int i, order;

  order = (M->Pirho != 0e0)? 1 : 0;
  for (i = -HOMmax; i <= HOMmax; i++)
    if (M->PB[i+HOMmax] != 0e0 && abs(i) > order) order = abs(i);
  return order;
}


arma::mat get_edge_mat(const double h, const double phi, const double gap,
		       const double delta)
{
  arma::mat
    Id = arma::mat(ss_dim, ss_dim),
    M  = arma::mat(ss_dim, ss_dim);

  M.eye(ss_dim, ss_dim);
  M(px_, x_) += h*tan(degtorad(phi));
  M(py_, y_) -= h*tan(degtorad(phi)-get_psi(h, phi, gap));

  return M;
}


arma::mat get_sbend_mat(const double L, const double h, const double b2,
			const double delta)
{
  double K_x, K_y, psi_x, psi_y;
  arma::mat
    Id = arma::mat(ss_dim, ss_dim),
    M  = arma::mat(ss_dim, ss_dim);

  K_x = b2 + sqr(h);
  K_y = fabs(b2);
  psi_x = sqrt(fabs(K_x)/(1e0+delta))*L;
  psi_y = sqrt(K_y/(1e0+delta))*L;

  M.eye(ss_dim, ss_dim);
  if (K_x > 0e0) {
    M(x_, x_)      = cos(psi_x);
    M(x_, px_)     = sin(psi_x)/sqrt(K_x*(1e0+delta));
    M(x_, delta_)  = (1e0-cos(psi_x))*h/K_x;
    M(px_, x_)     = -sqrt(K_x*(1e0+delta))*sin(psi_x);
    M(px_, px_)    = cos(psi_x);
    M(px_, delta_) = sin(psi_x)*sqrt(1e0+delta)*h/sqrt(K_x);

    if (psi_y != 0e0) {
      M(y_, y_)   = cosh(psi_y);
      M(y_, py_)  = sinh(psi_y)/sqrt(K_y*(1e0+delta));
      M(py_, y_)  = sqrt(K_y*(1e0+delta))*sinh(psi_y);
      M(py_, py_) = cosh(psi_y);
    } else
      M(y_, py_)  += L/(1e0+delta);
 
    M(ct_, x_)     = sin(psi_x)*sqrt(1e0+delta)*h/sqrt(K_x);
    M(ct_, px_)    = (1e0-cos(psi_x))*h/K_x;
    M(ct_, delta_) =
      (psi_x-sin(psi_x))*sqrt(1e0+delta)*sqr(h)/pow(K_x, 3e0/2e0);
  } else if (K_x < 0e0) {
    K_x = fabs(K_x);
    M(x_, x_)      = cosh(psi_x);
    M(x_, px_)     = sinh(psi_x)/sqrt(K_x*(1e0+delta));
    M(x_, delta_)  = -(1e0-cosh(psi_x))*h/K_x;
    M(px_, x_)     = sqrt(K_x*(1e0+delta))*sinh(psi_x);
    M(px_, px_)    = cosh(psi_x);
    M(px_, delta_) = sinh(psi_x)*sqrt(1e0+delta)*h/sqrt(K_x);

    if (psi_y != 0e0) {
      M(y_, y_)    = cos(psi_y);
      M(y_, py_)   = sin(psi_y)/sqrt(K_y*(1e0+delta));
      M(py_, y_)   = -sqrt(K_y*(1e0+delta))*sin(psi_y);
      M(py_, py_)  = cos(psi_y);
   } else
      M(y_, py_)   = L/(1e0+delta);

    M(ct_, x_)     = sinh(psi_x)*sqrt(1e0+delta)*h/sqrt(K_x);
    M(ct_, px_)    = -(1e0-cosh(psi_x))*h/K_x;
    M(ct_, delta_) =
      - (psi_x-sinh(psi_x))*sqrt(1e0+delta)*sqr(h)/pow(K_x, 3e0/2e0);
  } else {
    // K_x = 0.
    M(x_, px_) = L/(1e0+delta);
    M(y_, py_) = L/(1e0+delta);
  }

  return M;
}


arma::mat get_thin_kick_mat(const double b2L, const double delta)
{
  arma::mat M(ss_dim, ss_dim);

  M.eye(ss_dim, ss_dim);
  M(x_, px_) = b2L;
  M(y_, py_) = b2L;

  return M;
}


arma::mat get_mat(ElemType *elem, const double delta)
{
  arma::mat M0(ss_dim, ss_dim), M1(ss_dim, ss_dim), M2(ss_dim, ss_dim);

  const MpoleType* M = dynamic_cast<const MpoleType*>(elem);

  if (elem->PL != 0e0) {
    M0 = get_sbend_mat(elem->PL, M->Pirho, M->PB[Quad+HOMmax], delta);
    M1 = get_edge_mat(M->Pirho, M->PTx1, M->Pgap, delta);
    M2 = get_edge_mat(M->Pirho, M->PTx2, M->Pgap, delta);
    M0 = M2*M0*M1;
  } else
    M0 = get_thin_kick_mat(elem->PL*M->PB[Quad+HOMmax], delta);

  return M0;
}


void LatticeType::get_mats(const double delta)
{
  long int  k;
  MpoleType *M;

  for (k = 0; k <= conf.Cell_nLoc; k++) {
    M = dynamic_cast<MpoleType*>(elems[k]);
    if (elems[k]->Pkind == Mpole) M->M_lin = get_mat(elems[k], delta);
  }
}


ElemType* MpoleType::Elem_Init(const ConfigType &conf, const bool reverse)
{
  double     phi;
  MpoleType  *Mp;

  MpoleType* M = dynamic_cast<MpoleType*>(this);

  M->PB = M->PBpar; M->Porder = UpdatePorder(M);

  /* set entrance and exit angles */
  M->dT[X_] = cos(degtorad(M->PdTpar)); M->dT[Y_] = sin(degtorad(M->PdTpar));

  /* set displacement to zero */
  M->dS[X_] = 0e0; M->dS[Y_] = 0e0;

  if (M->PL != 0e0 || M->Pirho != 0e0) {
    /* Thick element or radius non zero element */
    M->Pthick = pthicktype(thick);
    /* sin(L*irho/2) =sin(theta/2) half the angle */
    M->Pc0 = sin(M->PL*M->Pirho/2e0);
    /* cos roll: sin(theta/2)*cos(dT) */
    M->Pc1 = M->dT[X_]*M->Pc0;
    /* sin roll: sin(theta/2)*cos(dT) */
    M->Ps1 = M->dT[Y_]*M->Pc0;
  } else /* element as thin lens */
    M->Pthick = pthicktype(thin);

  // Allocate transport matrix.
  M->M_lin = get_mat(this, 0e0);

  Mp = Mpole_Alloc();
  *Mp = *M;

  Mp->Reverse = reverse;
  if (conf.reverse_elem && (Mp->Reverse == true)) {
    // Swap entrance and exit angles.
    phi = Mp->PTx1; Mp->PTx1 = Mp->PTx2; Mp->PTx2 = phi;
  }

  return Mp;
}


ElemType* CavityType::Elem_Init(const ConfigType &conf, const bool reverse)
{
  CavityType *Cp;

  const CavityType* C =
    dynamic_cast<const CavityType*>(this);

  Cp = Cavity_Alloc();
  *Cp = *C;
  return Cp;
}


ElemType* MarkerType::Elem_Init(const ConfigType &conf, const bool reverse)
{
  MarkerType *Mkp;

  const MarkerType* Mk =
    dynamic_cast<const MarkerType*>(this);

  Mkp = Marker_Alloc();
  *Mkp = *Mk;
  return Mkp;
}


ElemType* WigglerType::Elem_Init(const ConfigType &conf, const bool reverse)
{
  int         i;
  WigglerType *Wp;

  WigglerType* W = dynamic_cast<WigglerType*>(this);

  for (i = -HOMmax; i <= HOMmax; i++)
    W->PBW.push_back(0e0);
  Wp = Wiggler_Alloc();
  *Wp = *W;
  return Wp;
}


ElemType* InsertionType::Elem_Init(const ConfigType &conf, const bool reverse)
{
  InsertionType *IDp;

  const InsertionType* ID =
    dynamic_cast<const InsertionType*>(this);

//  ID->Porder = order;
//  x = ID->PBW[Quad+HOMmax];
  IDp = Insertion_Alloc();
  *IDp = *ID;
  return IDp;
}


ElemType* FieldMapType::Elem_Init(const ConfigType &conf, const bool reverse)
{
  FieldMapType *FMp;

  const FieldMapType *FM =
    dynamic_cast<const FieldMapType*>(this);

  FMp = FieldMap_Alloc();
  *FMp = *FM;
  return FMp;
}


ElemType* SpreaderType::Elem_Init(const ConfigType &conf, const bool reverse)
{
  SpreaderType *Sprp;

  const SpreaderType* Spr =
    dynamic_cast<const SpreaderType*>(this);

  Sprp = Spreader_Alloc();
  *Sprp = *Spr;
  return Sprp;
}


ElemType* RecombinerType::Elem_Init(const ConfigType &conf, const bool reverse)
{
  RecombinerType *Recp;

  const RecombinerType* Rec =
    dynamic_cast<const RecombinerType*>(this);

  Recp = Recombiner_Alloc();
  *Recp = *Rec;
  return Recp;
}


ElemType* SolenoidType::Elem_Init(const ConfigType &conf, const bool reverse)
{
  SolenoidType *Solp;

  const SolenoidType* Sol =
    dynamic_cast<const SolenoidType*>(this);

  Solp = Solenoid_Alloc();
  *Solp = *Sol;
  return Solp;
}


ElemType* MapType::Elem_Init(const ConfigType &conf, const bool reverse)
{
  MapType *Mapp;

  const MapType* Map = dynamic_cast<const MapType*>(this);

  Mapp = Map_Alloc();
  *Mapp = *Map;
  Mapp->M.identity();
  return Mapp;
}


#if 0

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


void get_B_DIAMOND(ConfigType &conf, const char *filename,
		   FieldMapType *FM)
{
  char          line[max_str];
  int           i, j, n, ny;
  double        x0, y0, z0;
  std::ifstream inf;
  std::ofstream outf;

  const int     skip = 8;
  const double  Brho = conf.Energy*1e9/c0;

  std::cout << std::endl;
  std::cout << "get_B_DIAMOND: loading field map: " << filename << std::endl;

  file_rd(inf, filename);

  for (i = 1; i <= skip; i++)
    inf.getline(line, max_str);

  inf.getline(line, max_str);
  sscanf(line, "Starting point[cm] : %lf %lf %lf", &x0, &y0, &z0);
  inf.getline(line, max_str);
  sscanf(line, "Step size[cm]      : %lf %lf %lf",
	 &FM->dx[X_], &FM->dx[Y_], &FM->dx[Z_]);
  inf.getline(line, max_str);
  sscanf(line, "Number of points   : %d %d %d", &FM->n[X_], &ny, &FM->n[Z_]);

  // Convert from [cm] to [m].
  x0 *= 1e-2; y0 *= 1e-2; z0 *= 1e-2;
  FM->dx[X_] *= 1e-2; FM->dx[Y_] *= 1e-2; FM->dx[Z_] *= 1e-2;
  FM->Lr = FM->dx[Z_]*(FM->n[Z_]-1);

  FM->n[Y_] = 2*ny - 1;

  FM->x[X_] = dvector(1, FM->n[X_]); FM->x[Y_] = dvector(1, FM->n[Y_]);
  FM->x[Z_] = dvector(1, FM->n[Z_]);

  FM->BoBrho[X_]  = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho[Y_]  = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho[Z_]  = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho2[X_] = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho2[Y_] = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho2[Z_] = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);

  FM->AoBrho[X_]  = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->AoBrho[Y_]  = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->AoBrho2[X_] = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->AoBrho2[Y_] = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);

  for (i = 1; i <= 2; i++)
    inf.getline(line, max_str);

  for (n = 1; n <= FM->n[Z_]; n++)
    for (j = 1; j <= ny; j++)
      for (i = 1; i <= FM->n[X_]; i++) {
	inf.getline(line, max_str);
	sscanf(line, "%lf %lf %lf %lf %lf %lf",
	       &FM->x[X_][i], &FM->x[Y_][ny-1+j], &FM->x[Z_][n],
	       &FM->BoBrho[X_][n][i][ny-1+j],
	       &FM->BoBrho[Y_][n][i][ny-1+j],
	       &FM->BoBrho[Z_][n][i][ny-1+j]);

	// Convert from [cm] to [m].
	FM->x[X_][i] *= 1e-2; FM->x[Y_][ny-1+j] *= 1e-2; FM->x[Z_][n] *= 1e-2;
	// Convert from [Gauss] to [Tesla].
	FM->BoBrho[X_][n][i][ny-1+j] /= 1e+4*Brho;
	FM->BoBrho[Y_][n][i][ny-1+j] /= 1e+4*Brho;
	FM->BoBrho[Z_][n][i][ny-1+j] /= 1e+4*Brho;
	// Scale.
	FM->BoBrho[X_][n][i][ny-1+j] *= FM->scl;
	FM->BoBrho[Y_][n][i][ny-1+j] *= FM->scl;
	FM->BoBrho[Z_][n][i][ny-1+j] *= FM->scl;

	// Compute vector potential (axial gauge) by extended trapezodial rule.
	if (n == 1) {
	  FM->AoBrho[X_][n][i][ny-1+j] =
	    -FM->BoBrho[Y_][n][i][ny-1+j]*FM->dx[Z_]/2e0;
	  FM->AoBrho[Y_][n][i][ny-1+j] =
	    FM->BoBrho[X_][n][i][ny-1+j]*FM->dx[Z_]/2e0;
	} else if (n == FM->n[Z_]) {
	  FM->AoBrho[X_][n][i][ny-1+j] =
	    FM->AoBrho[X_][n-1][i][ny-1+j]
	    - FM->BoBrho[Y_][n][i][ny-1+j]*FM->dx[Z_]/2e0;
	  FM->AoBrho[Y_][n][i][ny-1+j] =
	    FM->AoBrho[Y_][n-1][i][ny-1+j]
	    + FM->BoBrho[X_][n][i][ny-1+j]*FM->dx[Z_]/2e0;
	} else {
	  FM->AoBrho[X_][n][i][ny-1+j] =
	    FM->AoBrho[X_][n-1][i][ny-1+j]
	    - FM->BoBrho[Y_][n][i][ny-1+j]*FM->dx[Z_];
	  FM->AoBrho[Y_][n][i][ny-1+j] =
	    FM->AoBrho[Y_][n-1][i][ny-1+j]
	    + FM->BoBrho[X_][n][i][ny-1+j]*FM->dx[Z_];
	}
      }

  inf.close();

  printf("\n%10.5f %10.5f %10.5f\n", x0, y0, z0);
  printf("%10.5f %10.5f %10.5f\n", FM->dx[X_], FM->dx[Y_], FM->dx[Z_]);
  printf("%10d %10d %10d\n", FM->n[X_], FM->n[Y_], FM->n[Z_]);
  printf("%10.3f -> %10.3f %10.3f -> %10.3f %10.3f -> %10.3f\n",
	 FM->x[X_][1], FM->x[X_][FM->n[X_]],
	 FM->x[Y_][1], FM->x[Y_][FM->n[Y_]],
	 FM->x[Z_][1], FM->x[Z_][FM->n[Z_]]);
  printf("Magnet length [m]: %10.5f\n", FM->Lr);

  for (j = 1; j <= ny-1; j++) {
    FM->x[Y_][j] = -FM->x[Y_][2*ny-j];
    for (i = 1; i <= FM->n[X_]; i++) {
      for (n = 1; n <= FM->n[Z_]; n++) {
	// B[X_] is antisymmetric in y (rot(A) = 0)
	FM->BoBrho[X_][n][i][j] = -FM->BoBrho[X_][n][i][2*ny-j];
	FM->BoBrho[Y_][n][i][j] =  FM->BoBrho[Y_][n][i][2*ny-j];
	// Bz is antisymmetric in y
	FM->BoBrho[Z_][n][i][j] = -FM->BoBrho[Z_][n][i][2*ny-j];

	FM->AoBrho[X_][n][i][j] = FM->AoBrho[X_][n][i][2*ny-j];
	// Ay is antisymmetric in y
	FM->AoBrho[Y_][n][i][j] = -FM->AoBrho[Y_][n][i][2*ny-j];
      }
    }
  }

  if (true) {
    file_wr(outf, "field_map.dat");
    std::cout << std::scientific << std::setprecision(3)
	      << std::setw(11) << FM->x[X_][7] << std::setw(11) << FM->x[Y_][7]
	      << std::endl;
    for (i = 1; i <= FM->n[X_]; i++)
      outf << std::scientific << std::setprecision(3)
	   << std::setw(11) << FM->x[X_][i]
	   << std::setw(11) << FM->BoBrho[Y_][7][i][7] << std::endl;
    outf.close();
  }

  for (n = 1; n <= FM->n[Z_]; n++) {
    splie2_(FM->x[X_], FM->x[Y_], FM->BoBrho[X_][n],
	    FM->n[X_], FM->n[Y_], FM->BoBrho2[X_][n]);
    splie2_(FM->x[X_], FM->x[Y_], FM->BoBrho[Y_][n],
	    FM->n[X_], FM->n[Y_], FM->BoBrho2[Y_][n]);
    splie2_(FM->x[X_], FM->x[Y_], FM->BoBrho[Z_][n],
	    FM->n[X_], FM->n[Y_], FM->BoBrho2[Z_][n]);

    splie2_(FM->x[X_], FM->x[Y_], FM->AoBrho[X_][n],
	    FM->n[X_], FM->n[Y_], FM->AoBrho2[X_][n]);
    splie2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][n],
	    FM->n[X_], FM->n[Y_], FM->AoBrho2[Y_][n]);
  }

  std::cout << "field map loaded: " << filename << std::endl;

/*  free_dvector(FM->x[X_], 1, FM->n[X_]); free_dvector(FM->x[Y_], 1, FM->n[Y_]);
  free_dvector(FM->x[Z_], 1, FM->n[Z_]);

  free_df3tensor(FM->BoBrho[X_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho[Y_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho[Z_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho2[X_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho2[Y_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho2[Z_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);

  free_df3tensor(FM->AoBrho[X_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->AoBrho[Y_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->AoBrho2[X_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->AoBrho2[Y_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);*/
}


void get_B_NSLS_II(ConfigType &conf, const char *filename,
		   FieldMapType *FM)
{
  char          line[max_str];
  int           i, j, n;
  double        x_min[3], x_max[3];
  std::ifstream inf;

  const double  Brho = conf.Energy*1e9/c0;

  std::cout << std::endl;
  std::cout << "get_B_NSLS_II: loading field map: " << filename << std::endl;

  file_rd(inf, filename);

  inf.getline(line, max_str);
  inf.getline(line, max_str);
  sscanf(line, "#%lf <= x <= %lf, dx = %lf, nx = %d",
	 &x_min[X_], &x_max[X_], &FM->dx[X_], &FM->n[X_]);
  x_min[X_] *= 1e-2; x_max[X_] *= 1e-2; FM->dx[X_] *= 1e-2;
  inf.getline(line, max_str);
  sscanf(line, "#%lf <= y <= %lf, dy = %lf, ny = %d",
	 &x_min[Y_], &x_max[Y_], &FM->dx[Y_], &FM->n[Y_]);
  x_min[Y_] *= 1e-2; x_max[Y_] *= 1e-2; FM->dx[Y_] *= 1e-2;
  inf.getline(line, max_str);
  sscanf(line, "#%lf <= z <= %lf, dz = %lf, nz = %d",
	 &x_min[Z_], &x_max[Z_], &FM->dx[Z_], &FM->n[Z_]);
  x_min[Z_] *= 1e-2; x_max[Z_] *= 1e-2; FM->dx[Z_] *= 1e-2;

  FM->x[X_] = dvector(1, FM->n[X_]); FM->x[Y_] = dvector(1, FM->n[Y_]);
  FM->x[Z_] = dvector(1, FM->n[Z_]);

  FM->BoBrho[X_]  = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho[Y_]  = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho[Z_]  = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho2[X_] = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho2[Y_] = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho2[Z_] = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);

  FM->AoBrho[X_]  = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->AoBrho[Y_]  = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->AoBrho2[X_] = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->AoBrho2[Y_] = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);

  for (i = 1; i <= FM->n[X_]; i++)
    for (j = 1; j <= FM->n[Y_]; j++)
      for (n = 1; n <= FM->n[Z_]; n++) {
	inf.getline(line, max_str);
	sscanf(line, "%lf %lf %lf %lf %lf %lf",
	       &FM->x[X_][i], &FM->x[Y_][j], &FM->x[Z_][n],
	       &FM->BoBrho[X_][n][i][j],
	       &FM->BoBrho[Y_][n][i][j],
	       &FM->BoBrho[Z_][n][i][j]);

	// convert from cm to m
	FM->x[X_][i] *= 1e-2; FM->x[Y_][j] *= 1e-2; FM->x[Z_][n] *= 1e-2;

	FM->BoBrho[X_][n][i][j] = Brho;
	FM->BoBrho[Y_][n][i][j] = Brho;
	FM->BoBrho[Z_][n][i][j] = Brho;

	// Compute vector potential (axial gauge) by extended trapezodial rule
 	if (n == 1) {
	  FM->AoBrho[X_][n][i][j] = -FM->BoBrho[Y_][n][i][j]*FM->dx[Z_]/2e0;
	  FM->AoBrho[Y_][n][i][j] =  FM->BoBrho[X_][n][i][j]*FM->dx[Z_]/2e0;
	} else if (n == FM->n[Z_]) {
	  FM->AoBrho[X_][n][i][j] =
	    FM->AoBrho[X_][n-1][i][j] - FM->BoBrho[Y_][n][i][j]*FM->dx[Z_]/2e0;
	  FM->AoBrho[Y_][n][i][j] =
	    FM->AoBrho[Y_][n-1][i][j] + FM->BoBrho[X_][n][i][j]*FM->dx[Z_]/2e0;
	} else {
	  FM->AoBrho[X_][n][i][j] =
	    FM->AoBrho[X_][n-1][i][j] - FM->BoBrho[Y_][n][i][j]*FM->dx[Z_];
	  FM->AoBrho[Y_][n][i][j] =
	    FM->AoBrho[Y_][n-1][i][j] + FM->BoBrho[X_][n][i][j]*FM->dx[Z_];
	}
     }

  inf.close();

  FM->Lr = FM->dx[Z_]*(FM->n[Z_]-1);

  std::cout << std::fixed << std::setprecision(5)
	    << std::setw(10) << 1e3*FM->dx[X_]
	    << std::setw(10) << 1e3*FM->dx[Y_]
	    << std::setw(10) << 1e3*FM->dx[Z_] << std::endl;
  std::cout << std::setw(10) << FM->n[X_] << std::setw(10) << FM->n[Y_]
	    << std::setw(10) << FM->n[Z_] << std::endl;
  std::cout << std::fixed << std::setprecision(3)
	    << std::setw(10) << FM->x[X_][1]
	    << std::setw(10) << FM->x[X_][FM->n[X_]]
	    << std::setw(10) << FM->x[Y_][1]
	    << std::setw(10) << FM->x[Y_][FM->n[Y_]]
	    << std::setw(10) << FM->x[Z_][1]
	    << std::setw(10) << FM->x[Z_][FM->n[Z_]] << std::endl;
  std::cout << std::fixed << std::setprecision(5)
	    << "Magnet length [m]:" << std::setw(10) << FM->Lr << std::endl;

  for (n = 1; n <= FM->n[Z_]; n++) {
    splie2_(FM->x[X_], FM->x[Y_], FM->BoBrho[X_][n],
	    FM->n[X_], FM->n[Y_], FM->BoBrho2[X_][n]);
    splie2_(FM->x[X_], FM->x[Y_], FM->BoBrho[Y_][n],
	    FM->n[X_], FM->n[Y_], FM->BoBrho2[Y_][n]);
    splie2_(FM->x[X_], FM->x[Y_], FM->BoBrho[Z_][n],
	    FM->n[X_], FM->n[Y_], FM->BoBrho2[Z_][n]);

    splie2_(FM->x[X_], FM->x[Y_], FM->AoBrho[X_][n],
	    FM->n[X_], FM->n[Y_], FM->AoBrho2[X_][n]);
    splie2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][n],
	    FM->n[X_], FM->n[Y_], FM->AoBrho2[Y_][n]);
  }

  std::cout << "field map loaded: " << filename << std::endl;

/*  free_dvector(FM->x[X_], 1, FM->n[X_]);
  free_dvector(FM->x[Y_], 1, FM->n[Y_]);
  free_dvector(FM->x[Z_], 1, FM->n[Z_]);

  free_df3tensor(FM->BoBrho[X_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho[Y_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho[Z_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho2[X_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho2[Y_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho2[Z_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);

  free_df3tensor(FM->AoBrho[X_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->AoBrho[Y_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->AoBrho2[X_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->AoBrho2[Y_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);*/
}


void get_B_Oleg1(ConfigType &conf, const char *filename, FieldMapType *FM)
{
  char          line[max_str];
  int           i, j, n;
  double        x_min[3], x_max[3];
  std::ifstream inf;

  const double Brho = conf.Energy*1e9/c0;

  std::cout << std::endl;
  std::cout << "get_B_Oleg1: loading field map: " << filename << std::endl;

  file_rd(inf, filename);

  inf.getline(line, max_str);
  inf.getline(line, max_str);
  sscanf(line, "#%lf <= x <= %lf, dx = %lf, nx = %d",
	 &x_min[X_], &x_max[X_], &FM->dx[X_], &FM->n[X_]);
  FM->dx[X_] *= 1e-3;
  inf.getline(line, max_str);
  sscanf(line, "#%lf <= y <= %lf, dy = %lf, ny = %d",
	 &x_min[Y_], &x_max[Y_], &FM->dx[Y_], &FM->n[Y_]);
  FM->dx[Y_] *= 1e-3;
  inf.getline(line, max_str);
  sscanf(line, "#%lf <= z <= %lf, dz = %lf, nz = %d",
	 &x_min[Z_], &x_max[Z_], &FM->dx[Z_], &FM->n[Z_]);
  FM->dx[Z_] *= 1e-3;

  FM->x[X_] = dvector(1, FM->n[X_]); FM->x[Y_] = dvector(1, FM->n[Y_]);
  FM->x[Z_] = dvector(1, FM->n[Z_]);

  FM->BoBrho[X_]  = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho[Y_]  = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho[Z_]  = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho2[X_] = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho2[Y_] = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho2[Z_] = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);

  FM->AoBrho[X_]  = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->AoBrho[Y_]  = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->AoBrho2[X_] = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->AoBrho2[Y_] = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);

  for (n = 1; n <= FM->n[Z_]; n++)
    for (j = 1; j <= FM->n[Y_]; j++)
      for (i = 1; i <= FM->n[X_]; i++) {
	inf.getline(line, max_str);
	sscanf(line, "%lf %lf %lf %lf %lf %lf",
	       &FM->x[X_][i], &FM->x[Y_][j], &FM->x[Z_][n],
	       &FM->BoBrho[X_][n][i][j],
	       &FM->BoBrho[Y_][n][i][j],
	       &FM->BoBrho[Z_][n][i][j]);

	// convert from mm to m
	FM->x[X_][i] *= 1e-3; FM->x[Y_][j] *= 1e-3; FM->x[Z_][n] *= 1e-3;

	FM->BoBrho[X_][n][i][j] /= Brho;
	FM->BoBrho[Y_][n][i][j] /= Brho;
	FM->BoBrho[Z_][n][i][j] /= Brho;

	// Compute vector potential (axial gauge) by extended trapezodial rule
 	if (n == 1) {
	  FM->AoBrho[X_][n][i][j] = -FM->BoBrho[Y_][n][i][j]*FM->dx[Z_]/2e0;
	  FM->AoBrho[Y_][n][i][j] =  FM->BoBrho[X_][n][i][j]*FM->dx[Z_]/2e0;
	} else if (n == FM->n[Z_]) {
	  FM->AoBrho[X_][n][i][j] =
	    FM->AoBrho[X_][n-1][i][j] - FM->BoBrho[Y_][n][i][j]*FM->dx[Z_]/2e0;
	  FM->AoBrho[Y_][n][i][j] =
	    FM->AoBrho[Y_][n-1][i][j] + FM->BoBrho[X_][n][i][j]*FM->dx[Z_]/2e0;
	} else {
	  FM->AoBrho[X_][n][i][j] =
	    FM->AoBrho[X_][n-1][i][j] - FM->BoBrho[Y_][n][i][j]*FM->dx[Z_];
	  FM->AoBrho[Y_][n][i][j] =
	    FM->AoBrho[Y_][n-1][i][j] + FM->BoBrho[X_][n][i][j]*FM->dx[Z_];
	}
     }

  inf.close();

  FM->Lr = FM->dx[Z_]*(FM->n[Z_]-1);

  std::cout << std::fixed << std::setprecision(5)
	    << std::setw(10) << 1e3*FM->dx[X_]
	    << std::setw(10) << 1e3*FM->dx[Y_]
	    << std::setw(10) << 1e3*FM->dx[Z_] << std::endl;
  std::cout << std::setw(10) << FM->n[X_] << std::setw(10) << FM->n[Y_]
	    << std::setw(10) << FM->n[Z_] << std::endl;
  std::cout << std::fixed << std::setprecision(3)
	    << std::setw(10) << FM->x[X_][1]
	    << std::setw(10) << FM->x[X_][FM->n[X_]]
	    << std::setw(10) << FM->x[Y_][1]
	    << std::setw(10) << FM->x[Y_][FM->n[Y_]]
	    << std::setw(10) << FM->x[Z_][1]
	    << std::setw(10) << FM->x[Z_][FM->n[Z_]] << std::endl;
  std::cout << std::fixed << std::setprecision(5)
	    << "Magnet length [m]:" << std::setw(10) << FM->Lr << std::endl;

  for (n = 1; n <= FM->n[Z_]; n++) {
    splie2_(FM->x[X_], FM->x[Y_], FM->BoBrho[X_][n],
	    FM->n[X_], FM->n[Y_], FM->BoBrho2[X_][n]);
    splie2_(FM->x[X_], FM->x[Y_], FM->BoBrho[Y_][n],
	    FM->n[X_], FM->n[Y_], FM->BoBrho2[Y_][n]);
    splie2_(FM->x[X_], FM->x[Y_], FM->BoBrho[Z_][n],
	    FM->n[X_], FM->n[Y_], FM->BoBrho2[Z_][n]);

    splie2_(FM->x[X_], FM->x[Y_], FM->AoBrho[X_][n],
	    FM->n[X_], FM->n[Y_], FM->AoBrho2[X_][n]);
    splie2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][n],
	    FM->n[X_], FM->n[Y_], FM->AoBrho2[Y_][n]);
  }

  std::cout << "field map loaded: " << filename << std::endl;

/*  free_dvector(FM->x[X_], 1, FM->n[X_]);
  free_dvector(FM->x[Y_], 1, FM->n[Y_]);
  free_dvector(FM->x[Z_], 1, FM->n[Z_]);

  free_df3tensor(FM->BoBrho[X_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho[Y_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho[Z_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho2[X_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho2[Y_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho2[Z_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);

  free_df3tensor(FM->AoBrho[X_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->AoBrho[Y_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->AoBrho2[X_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->AoBrho2[Y_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);*/
}


void get_B_Oleg2(ConfigType &conf, const char *filename, FieldMapType *FM)
{
  char          line[max_str];
  int           i, j, n;
  double        x_min[3];
  std::ifstream inf;

  const double Brho = conf.Energy*1e9/c0;

  std::cout << std::endl;
  std::cout << "get_B_Oleg2: loading field map: " << filename << std::endl;

  file_rd(inf, filename);

  inf.getline(line, max_str);

  inf.getline(line, max_str); sscanf(line, "#%lf", &x_min[X_]);
  inf.getline(line, max_str); sscanf(line, "#%lf", &FM->dx[X_]);
  inf.getline(line, max_str); sscanf(line, "#%d", &FM->n[X_]);

  inf.getline(line, max_str); sscanf(line, "#%lf", &x_min[Y_]);
  inf.getline(line, max_str); sscanf(line, "#%lf", &FM->dx[Y_]);
  inf.getline(line, max_str); sscanf(line, "#%d", &FM->n[Y_]);

  inf.getline(line, max_str); sscanf(line, "#%lf", &x_min[Z_]);
  inf.getline(line, max_str); sscanf(line, "#%lf", &FM->dx[Z_]);
  inf.getline(line, max_str); sscanf(line, "#%d", &FM->n[Z_]);

  std::cout << std::fixed << std::setprecision(5)
	    << std::setw(10) << 1e3*FM->dx[X_]
	    << std::setw(10) << 1e3*FM->dx[Y_]
	    << std::setw(10) << 1e3*FM->dx[Z_] << std::endl;
  std::cout << std::setw(10) << FM->n[X_] << std::setw(10) << FM->n[Y_]
	    << std::setw(10) << FM->n[Z_] << std::endl;
  std::cout << std::fixed << std::setprecision(3)
	    << std::setw(10) << x_min[X_] << std::setw(10) << x_min[Y_]
	    << std::setw(10) << x_min[Z_] << std::endl;

  FM->x[X_] = dvector(1, FM->n[X_]); FM->x[Y_] = dvector(1, FM->n[Y_]);
  FM->x[Z_] = dvector(1, FM->n[Z_]);

  FM->BoBrho[X_]  = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho[Y_]  = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho[Z_]  = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho2[X_] = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho2[Y_] = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho2[Z_] = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);

  FM->AoBrho[X_]  = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->AoBrho[Y_]  = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->AoBrho2[X_] = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->AoBrho2[Y_] = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);

  for (n = 1; n <= FM->n[Z_]; n++) {
    FM->x[Z_][n] = (n == 1)? x_min[Z_] : FM->x[Z_][n-1] + FM->dx[Z_];

    for (j = 1; j <= FM->n[Y_]; j++) {
      FM->x[Y_][j] = (j == 1)? x_min[Y_] : FM->x[Y_][j-1] + FM->dx[Y_];

      for (i = 1; i <= FM->n[X_]; i++) {
	FM->x[X_][i] = (i == 1)? x_min[X_] : FM->x[X_][i-1] + FM->dx[X_];

	inf.getline(line, max_str);
	sscanf(line, "%lf %lf %lf",
	       &FM->BoBrho[X_][n][i][j], &FM->BoBrho[Y_][n][i][j],
	       &FM->BoBrho[Z_][n][i][j]);

	FM->BoBrho[X_][n][i][j] /= Brho;
	FM->BoBrho[Y_][n][i][j] /= Brho;
	FM->BoBrho[Z_][n][i][j] /= Brho;

	// Compute vector potential (axial gauge) by extended trapezodial rule
 	if (n == 1) {
	  FM->AoBrho[X_][n][i][j] = -FM->BoBrho[Y_][n][i][j]*FM->dx[Z_]/2e0;
	  FM->AoBrho[Y_][n][i][j] =  FM->BoBrho[X_][n][i][j]*FM->dx[Z_]/2e0;
	} else if (n == FM->n[Z_]) {
	  FM->AoBrho[X_][n][i][j] =
	    FM->AoBrho[X_][n-1][i][j] - FM->BoBrho[Y_][n][i][j]*FM->dx[Z_]/2e0;
	  FM->AoBrho[Y_][n][i][j] =
	    FM->AoBrho[Y_][n-1][i][j] + FM->BoBrho[X_][n][i][j]*FM->dx[Z_]/2e0;
	} else {
	  FM->AoBrho[X_][n][i][j] =
	    FM->AoBrho[X_][n-1][i][j] - FM->BoBrho[Y_][n][i][j]*FM->dx[Z_];
	  FM->AoBrho[Y_][n][i][j] =
	    FM->AoBrho[Y_][n-1][i][j] + FM->BoBrho[X_][n][i][j]*FM->dx[Z_];
	}
      }
    }
  }

  inf.close();

  FM->Lr = FM->dx[Z_]*(FM->n[Z_]-1);

  std::cout << std::fixed << std::setprecision(5)
	    << std::setw(10) << 1e3*FM->dx[X_]
	    << std::setw(10) << 1e3*FM->dx[Y_]
	    << std::setw(10) << 1e3*FM->dx[Z_] << std::endl;
  std::cout << std::setw(10) << FM->n[X_] << std::setw(10) << FM->n[Y_]
	    << std::setw(10) << FM->n[Z_] << std::endl;
  std::cout << std::fixed << std::setprecision(3)
	    << std::setw(10) << FM->x[X_][1]
	    << std::setw(10) << FM->x[X_][FM->n[X_]]
	    << std::setw(10) << FM->x[Y_][1]
	    << std::setw(10) << FM->x[Y_][FM->n[Y_]]
	    << std::setw(10) << FM->x[Z_][1]
	    << std::setw(10) << FM->x[Z_][FM->n[Z_]] << std::endl;
  std::cout << std::fixed << std::setprecision(5)
	    << "Magnet length [m]:" << std::setw(10) << FM->Lr << std::endl;

  for (n = 1; n <= FM->n[Z_]; n++) {
    splie2_(FM->x[X_], FM->x[Y_], FM->BoBrho[X_][n],
	    FM->n[X_], FM->n[Y_], FM->BoBrho2[X_][n]);
    splie2_(FM->x[X_], FM->x[Y_], FM->BoBrho[Y_][n],
	    FM->n[X_], FM->n[Y_], FM->BoBrho2[Y_][n]);
    splie2_(FM->x[X_], FM->x[Y_], FM->BoBrho[Z_][n],
	    FM->n[X_], FM->n[Y_], FM->BoBrho2[Z_][n]);
  }

  std::cout << "field map loaded: " << filename << std::endl;

/*  free_dvector(FM->x[X_], 1, FM->n[X_]);
  free_dvector(FM->x[Y_], 1, FM->n[Y_]);
  free_dvector(FM->x[Z_], 1, FM->n[Z_]);

  free_df3tensor(FM->BoBrho[X_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho[Y_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho[Z_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho2[X_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho2[Y_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho2[Z_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);

  free_df3tensor(FM->AoBrho[X_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->AoBrho[Y_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->AoBrho2[X_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->AoBrho2[Y_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);*/
}


void get_B_SRW(ConfigType &conf, const char *filename, FieldMapType *FM)
{
  char          line[max_str];
  int           i, j, n;
  double        x_min[3];
  std::ifstream inf;

  const double Brho = conf.Energy*1e9/c0;

  printf("\nget_B_SRW: loading field map: %s\n", filename);

  file_rd(inf, filename);

  inf.getline(line, max_str);

  inf.getline(line, max_str); sscanf(line, "#%lf", &x_min[X_]);
  inf.getline(line, max_str); sscanf(line, "#%lf", &FM->dx[X_]);
  inf.getline(line, max_str); sscanf(line, "#%d", &FM->n[X_]);

  inf.getline(line, max_str); sscanf(line, "#%lf", &x_min[Y_]);
  inf.getline(line, max_str); sscanf(line, "#%lf", &FM->dx[Y_]);
  inf.getline(line, max_str); sscanf(line, "#%d", &FM->n[Y_]);

  inf.getline(line, max_str); sscanf(line, "#%lf", &x_min[Z_]);
  inf.getline(line, max_str); sscanf(line, "#%lf", &FM->dx[Z_]);
  inf.getline(line, max_str); sscanf(line, "#%d", &FM->n[Z_]);

  printf("\n  dx [mm]   = [%7.5f, %7.5f, %7.5f]\n",
	 1e3*FM->dx[X_], 1e3*FM->dx[Y_], 1e3*FM->dx[Z_]);
  printf("  n         = [%d, %d, %d]\n", FM->n[X_], FM->n[Y_], FM->n[Z_]);
  printf("  x_min [m] = [%7.5f, %7.5f, %7.5f]\n",
	 x_min[X_], x_min[Y_], x_min[Z_]);

  FM->x[X_] = dvector(1, FM->n[X_]); FM->x[Y_] = dvector(1, FM->n[Y_]);
  FM->x[Z_] = dvector(1, FM->n[Z_]);

  FM->BoBrho[X_]  = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho[Y_]  = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho[Z_]  = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho2[X_] = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho2[Y_] = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho2[Z_] = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);

  FM->AoBrho[X_]  = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->AoBrho[Y_]  = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->AoBrho2[X_] = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->AoBrho2[Y_] = df3tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);

  for (n = 1; n <= FM->n[Z_]; n++) {
    FM->x[Z_][n] = (n == 1)? x_min[Z_] : FM->x[Z_][n-1] + FM->dx[Z_];

    for (j = 1; j <= FM->n[Y_]; j++) {
      FM->x[Y_][j] = (j == 1)? x_min[Y_] : FM->x[Y_][j-1] + FM->dx[Y_];

      for (i = 1; i <= FM->n[X_]; i++) {
	FM->x[X_][i] = (i == 1)? x_min[X_] : FM->x[X_][i-1] + FM->dx[X_];

	inf.getline(line, max_str);
	sscanf(line, "%lf %lf %lf",
	       &FM->BoBrho[X_][n][i][j], &FM->BoBrho[Y_][n][i][j],
	       &FM->BoBrho[Z_][n][i][j]);

	FM->BoBrho[X_][n][i][j] /= Brho;
	FM->BoBrho[Y_][n][i][j] /= Brho;
	FM->BoBrho[Z_][n][i][j] /= Brho;

	// Compute vector potential (axial gauge) by extended trapezodial rule
 	if (n == 1) {
	  FM->AoBrho[X_][n][i][j] = -FM->BoBrho[Y_][n][i][j]*FM->dx[Z_]/2e0;
	  FM->AoBrho[Y_][n][i][j] =  FM->BoBrho[X_][n][i][j]*FM->dx[Z_]/2e0;
	} else if (n == FM->n[Z_]) {
	  FM->AoBrho[X_][n][i][j] =
	    FM->AoBrho[X_][n-1][i][j] - FM->BoBrho[Y_][n][i][j]*FM->dx[Z_]/2e0;
	  FM->AoBrho[Y_][n][i][j] =
	    FM->AoBrho[Y_][n-1][i][j] + FM->BoBrho[X_][n][i][j]*FM->dx[Z_]/2e0;
	} else {
	  FM->AoBrho[X_][n][i][j] =
	    FM->AoBrho[X_][n-1][i][j] - FM->BoBrho[Y_][n][i][j]*FM->dx[Z_];
	  FM->AoBrho[Y_][n][i][j] =
	    FM->AoBrho[Y_][n-1][i][j] + FM->BoBrho[X_][n][i][j]*FM->dx[Z_];
	}
      }
    }
  }

  inf.close();

  FM->Lr = FM->dx[Z_]*(FM->n[Z_]-1);

  printf("\n  dx [mm]   = [%7.5f, %7.5f, %7.5f]\n",
	 1e3*FM->dx[X_], 1e3*FM->dx[Y_], 1e3*FM->dx[Z_]);
  printf("  n         = [%d, %d, %d]\n", FM->n[X_], FM->n[Y_], FM->n[Z_]);
  printf("  x [m]     = [%5.3f - %5.3f, %5.3f - %5.3f, %5.3f - %5.3f]\n",
	 FM->x[X_][1], FM->x[X_][FM->n[X_]],
	 FM->x[Y_][1], FM->x[Y_][FM->n[Y_]],
	 FM->x[Z_][1], FM->x[Z_][FM->n[Z_]]);

  printf("\n  Magnet length [m]: %7.5f\n", FM->Lr);

  for (n = 1; n <= FM->n[Z_]; n++) {
    splie2_(FM->x[X_], FM->x[Y_], FM->BoBrho[X_][n],
	    FM->n[X_], FM->n[Y_], FM->BoBrho2[X_][n]);
    splie2_(FM->x[X_], FM->x[Y_], FM->BoBrho[Y_][n],
	    FM->n[X_], FM->n[Y_], FM->BoBrho2[Y_][n]);
    splie2_(FM->x[X_], FM->x[Y_], FM->BoBrho[Z_][n],
	    FM->n[X_], FM->n[Y_], FM->BoBrho2[Z_][n]);
  }

  printf("\n  Field map loaded: %s\n", filename);

/*  free_dvector(FM->x[X_], 1, FM->n[X_]);
  free_dvector(FM->x[Y_], 1, FM->n[Y_]);
  free_dvector(FM->x[Z_], 1, FM->n[Z_]);

  free_df3tensor(FM->BoBrho[X_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho[Y_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho[Z_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho2[X_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho2[Y_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho2[Z_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);

  free_df3tensor(FM->AoBrho[X_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->AoBrho[Y_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->AoBrho2[X_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->AoBrho2[Y_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);*/
}

#endif


void get_B(ConfigType &conf, const char *filename, FieldMapType *FM)
{
  // Do not scale fieldmaps only Hamiltonians, i.e., the kick.  Note that RADIA
  // (2nd order) kick maps are quadratic in the field, and 1st order linear.

  switch (FieldMap_filetype) {
  case 1:
    // get_B_DIAMOND(conf, filename, FM);
    break;
  case 2:
    // get_B_NSLS_II(conf, filename, FM);
    break;
  case 3:
    // get_B_Oleg1(conf, filename, FM);
    break;
  case 4:
    // get_B_Oleg2(conf, filename, FM);
    break;
  case 5:
    // get_B_SRW(conf, filename, FM);
    break;
  default:
    printf("\nget_B: unknown FieldMap type %d", FieldMap_filetype);
    exit(1);
    break;
  }
}

void MpoleType::SetdS(void)
{
  int k;

  for (k = 0; k < 2; k++)
    dS[k] = PdSsys[k] + PdSrms[k]*PdSrnd[k];
}


void MpoleType::SetdT(void)
{
  dT[X_] = cos(degtorad(PdTpar + PdTsys + PdTrms*PdTrnd));
  dT[Y_] = sin(degtorad(PdTpar + PdTsys + PdTrms*PdTrnd));
  // Simplified p_rot.
  Pc0 = sin(PL*Pirho/2e0);
  Pc1 = cos(degtorad(PdTpar))*Pc0;
  Ps1 = sin(degtorad(PdTpar))*Pc0;
}


double MpoleType::GetdT(void)
{
  return PdTpar+PdTsys+PdTrms*PdTrnd;
}


void MpoleType::SetPB(const int n)
{
  // Compute full multipole composent as sum of design, systematic, and
  // random part.

  if ((1 <= abs(n)) && (abs(n) <= HOMmax)) {
    PB[n+HOMmax] =
      PBpar[n+HOMmax] + PBsys[n+HOMmax] + PBrms[n+HOMmax]*PBrnd[n+HOMmax];
    if (abs(n) > Porder && PB[n+HOMmax] != 0e0)
      Porder = abs(n);
  } else {
    printf("SetPB: |n| < 1 (%d)\n", n);
    exit(1);
  }
}


double MpoleType::GetPB(const int n)
{
  //  Return multipole strength of order n:
  //        /  2, normal quadrupole
  //    n = |
  //        \ -2, skew quadrupole                                      

  return PB[n+HOMmax];
}


void WigglerType::SetdS(void)
{
  int k;

  for (k = 0; k <= 1; k++)
    dS[k] = PdSsys[k] + PdSrms[k]*PdSrnd[k];
}


void WigglerType::SetdT(void)
{

  dT[X_] = cos(degtorad(PdTpar+PdTsys+PdTrms*PdTrnd));
  dT[Y_] = sin(degtorad(PdTpar+PdTsys+PdTrms*PdTrnd));
}


double WigglerType::GetdT(void)
{
  return PdTpar+PdTsys+PdTrms*PdTrnd;
}


double WigglerType::GetPB(const int n)
{
  //  Return multipole strength of order n:
  //        /  2, normal quadrupole
  //    n = |
  //        \ -2, skew quadrupole                                      

  return sqrt(2e0*PBW[n+HOMmax]);
}


void LatticeType::SetdS(const int Fnum, const int Knum)
{
  elems[Elem_GetPos(Fnum, Knum)]->SetdS();
}


void LatticeType::SetdT(const int Fnum, const int Knum)
{
  elems[Elem_GetPos(Fnum, Knum)]->SetdT();
}


double LatticeType::GetdT(const int Fnum, const int Knum)
{
  return elems[Elem_GetPos(Fnum, Knum)]->GetdT();
}


void LatticeType::SetPB(const int Fnum, const int Knum, const int n)
{
  elems[Elem_GetPos(Fnum, Knum)]->SetPB(n);
}

double LatticeType::GetPB(const int Fnum, const int Knum, const int n)
{
  return elems[Elem_GetPos(Fnum, Knum)]->GetPB(n);
}



void LatticeType::Mpole_DefPBpar(const int Fnum, const int Knum, const int n,
				 const double PBpar)
{
  MpoleType* M = dynamic_cast<MpoleType*>(elems[Elem_GetPos(Fnum, Knum)]);
  M->PBpar[n+HOMmax] = PBpar;
}


void LatticeType::Mpole_DefPBsys(const int Fnum, const int Knum, const int n,
				 const double PBsys)
{
  MpoleType* M = dynamic_cast<MpoleType*>(elems[Elem_GetPos(Fnum, Knum)]);
  M->PBsys[n+HOMmax] = PBsys;
}


void LatticeType::Mpole_DefdTpar(const int Fnum, const int Knum,
				 const double PdTpar)
{
  MpoleType* M = dynamic_cast<MpoleType*>(elems[Elem_GetPos(Fnum, Knum)]);
  M->PdTpar = PdTpar;
}


void LatticeType::Mpole_DefdTsys(const int Fnum, const int Knum,
				 const double PdTsys)
{
  MpoleType* M = dynamic_cast<MpoleType*>(elems[Elem_GetPos(Fnum, Knum)]);
  M->PdTsys=PdTsys;
}


double LatticeType::Elem_GetKval(const int Fnum, const int Knum, const int n)
{
  ElemType *elemp;
  double   b_2;

  elemp = elems[Elem_GetPos(Fnum, Knum)];
  b_2 = GetPB(Fnum, Knum, n);
  return (elemp->PL != 0e0)? elemp->PL*b_2 : b_2;
}
