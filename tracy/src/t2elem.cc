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
ElemFamType   ElemFam[Elem_nFamMax];
CellType      Cell[Cell_nLocMax+1];
std::ofstream outf_;

// for FieldMap
bool  sympl             = true;
int   FieldMap_filetype = 2;


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
inline T get_p_s(const ss_vect<T> &ps)
{
  T p_s, p_s2;

  if (!globval.H_exact)
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


void zero_mat(const int n, double** A)
{
  int i, j;

  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++)
      A[i][j] = 0e0;
}


void identity_mat(const int n, double** A)
{
  int i, j;

  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++)
      A[i][j] = (i == j)? 1e0 : 0e0;
}


double det_mat(const int n, double **A)
{
  int    i, *indx;
  double **U, d;

  indx = ivector(1, n); U = dmatrix(1, n, 1, n);

  dmcopy(A, n, n, U); dludcmp(U, n, indx, &d);

  for (i = 1; i <= n; i++)
    d *= U[i][i];

  free_dmatrix(U, 1, n, 1, n); free_ivector(indx, 1, n);

  return d;
}


double trace_mat(const int n, double **A)
{
  int    i;
  double d;

  d = 0e0;
  for (i = 1; i <= n; i++)
    d += A[i][i];

  return d;
}

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

  static inline double get_dI_eta(const ss_vect<tps> &A)
    {
      std::cout << "get_dI_eta: operation not defined for double" << std::endl;
      exit_(1);
      return 0e0;
    }

  static inline void emittance
  (CellType &Cell, const double B2, const double u, const double ps0,
   const ss_vect<double> &xp)
  { }

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

  static inline double get_dI_eta(const ss_vect<tps> &A)
  { return A[x_][delta_]; }

  static inline void emittance
  (CellType &Cell, const tps &B2_perp, const tps &ds, const tps &p_s0,
   const ss_vect<tps> &A)
  {
    // M. Sands "The Physics of Electron Storage Rings" SLAC-121, Eq. (5.20),
    // p. 118:
    //   dN<u^2>/E^2 =
    //     3*C_u*C_gamma*h_bar*c*E_0^5*(1+delta)^4*(B_perp/(Brho))^3
    //     /(4*pi*m_e^3 [eV/c^2])
    // A contains the eigenvectors.
    int          j;
    double       B_66, dD[3];
    ss_vect<tps> A_inv;

    bool prt_debug = false;

    if (B2_perp > 0e0) {
      B_66 = (q_fluct*pow(B2_perp.cst(), 1.5)*pow(p_s0, 4)*ds).cst();
      A_inv = Inv(A);
      // D_11 = D_22 = curly_H_x,y * B_66 / 2,
      // curly_H_x,y = eta_Fl^2 + etap_Fl^2

      for (j = 0; j < 3; j++) {
	dD[j] = (sqr(A_inv[j*2][delta_])+sqr(A_inv[j*2+1][delta_]))*B_66/2e0;
	globval.D_rad[j] += dD[j];
	Cell.dD[j] = globval.D_rad[j];
      }
      if (prt_debug)
	printf("emittance:\n  %8s\n  B_66 = %12.5e\n"
	       "  dD   = [%12.5e, %12.5e, %12.5e]\n"
	       "  D    = [%12.5e, %12.5e, %12.5e]\n",
	       Cell.Elem.PName, B_66, dD[X_], dD[Y_], dD[Z_],
	       globval.D_rad[X_], globval.D_rad[Y_], globval.D_rad[Z_]);
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
  const bool prt_debug = false;

  // Compute B_perp^2 and B_par^2.
  T xn, e[3];

  xn = 1e0/sqrt(sqr(1e0+xp[x_]*h_ref)+sqr(xp[px_])+sqr(xp[py_]));
  e[X_] = xp[px_]*xn; e[Y_] = xp[py_]*xn; e[Z_] = (1e0+xp[x_]*h_ref)*xn;

  // Left-handed coordinate system.
  B2_perp =
    sqr(B[Y_]*e[Z_]-B[Z_]*e[Y_]) + sqr(B[X_]*e[Y_]-B[Y_]*e[X_])
    + sqr(B[Z_]*e[X_]-B[X_]*e[Z_]);

  if (prt_debug)
    cout << scientific << setprecision(5)
	 << "\nget_B2:" << setw(13) << B2_perp << "\n";

  //  B2_par = sqr(B[X_]*e[X_]+B[Y_]*e[Y_]+B[Z_]*e[Z_]);
}


template<typename T>
void radiate
(CellType &Cell, ss_vect<T> &ps, const double L, const double h_ref,
 const T B[])
{
  // M. Sands "The Physics of Electron Storage Rings" SLAC-121, p. 98.
  // ddelta/d(ds) = -C_gamma*E_0^3*(1+delta)^2*(B_perp/(Brho))^2/(2*pi)
  T          p_s0, p_s1, ds, B2_perp = 0e0, B2_par = 0e0;
  ss_vect<T> cs;

  const bool prt_debug = false;

  if (prt_debug)
    cout << scientific << setprecision(5)
	 << "\nradiate ->:\n" << setw(13) << ps << "\n";

  // Large ring: x' and y' unchanged.
  p_s0 = get_p_s(ps); cs = ps; cs[px_] /= p_s0; cs[py_] /= p_s0;

  // H = -p_s => ds = H*L.
  ds = (1e0+cs[x_]*h_ref+(sqr(cs[px_])+sqr(cs[py_]))/2e0)*L;
  get_B2(h_ref, B, cs, B2_perp, B2_par);

  if (globval.radiation) {
    ps[delta_] -= cl_rad*sqr(p_s0)*B2_perp*ds;
    p_s1 = get_p_s(ps); ps[px_] = cs[px_]*p_s1; ps[py_] = cs[py_]*p_s1;
  }

  if (globval.emittance)
    is_tps<T>::emittance(Cell, B2_perp, ds, p_s0, cs);

  if (prt_debug)
    cout << scientific << setprecision(5)
	 << "\n<- radiate:\n" << setw(13) << ps << "\n";
}


template<typename T>
void radiate_ID
(CellType &Cell, ss_vect<T> &ps, const double L, const T &B2_perp)
{
  T          p_s0, p_s1, ds;
  ss_vect<T> cs;

  // Large ring: x' and y' unchanged.
  cs = ps; p_s0 = get_p_s(ps); cs[px_] /= p_s0; cs[py_] /= p_s0;

  // H = -p_s => ds = H*L.
  ds = (1e0+(sqr(cs[px_])+sqr(cs[py_]))/2e0)*L;

  if (globval.radiation) {
    ps[delta_] -= cl_rad*sqr(p_s0)*B2_perp*ds;
    p_s1 = get_p_s(ps); ps[px_] = cs[px_]*p_s1; ps[py_] = cs[py_]*p_s1;
  }

  if (globval.emittance)
    is_tps<T>::emittance(Cell, B2_perp, ds, p_s0, cs);
}


template<typename T>
void Drift(const double L, ss_vect<T> &ps)
{
  T u;

  if (!globval.H_exact) {
    // Small angle axproximation.
    u = L/(1e0+ps[delta_]);
    ps[x_]  += u*ps[px_]; ps[y_] += u*ps[py_];
    ps[ct_] += u*(sqr(ps[px_])+sqr(ps[py_]))/(2e0*(1e0+ps[delta_]));
  } else {
    u = L/get_p_s(ps);
    ps[x_]  += u*ps[px_]; ps[y_] += u*ps[py_];
    ps[ct_] += u*(1e0+ps[delta_]) - L;
  }
  if (globval.pathlength) ps[ct_] += L;
}


template<typename T>
void Drift_Pass(CellType &Cell, ss_vect<T> &x)
{
  Drift(Cell.Elem.PL, x);

  if (globval.emittance && !globval.Cavity_on)
    // Needs A^-1.
    Cell.curly_dH_x = is_tps<tps>::get_curly_H(x);
}


static double get_psi(double irho, double phi, double gap)
{
  /* References:
     1. H. Enge 𝐸𝑓𝑓𝑒𝑐𝑡 𝑜𝑓 𝐸𝑥𝑡𝑒𝑛𝑑𝑒𝑑 𝐹𝑟𝑖𝑛𝑔𝑖𝑛𝑔 𝐹𝑖𝑒𝑙𝑑𝑠 𝑜𝑛 𝐼𝑜𝑛-𝐹𝑜𝑐𝑢𝑠𝑖𝑛𝑔 𝑃𝑟𝑜𝑝𝑒𝑟𝑡𝑖𝑒𝑠 𝑜𝑓 𝐷𝑒𝑓𝑙𝑒𝑐𝑡𝑖𝑛𝑔
        𝑀𝑎𝑔𝑛𝑒𝑡𝑠 Rev. Sci. Instr. 35, 278-287 (1964).

        https://doi.org/10.1063/1.1718806

     2. H. Enge 𝐸𝑓𝑓𝑒𝑐𝑡𝑠 𝑜𝑓 𝐸𝑥𝑡𝑒𝑛𝑑𝑒𝑑 𝐹𝑟𝑖𝑛𝑔𝑖𝑛𝑔 𝐹𝑖𝑒𝑙𝑑𝑠 𝐹𝑜𝑐𝑢𝑠𝑖𝑛𝑔 𝑜𝑓 𝐶ℎ𝑎𝑟𝑔𝑒𝑑 𝑃𝑎𝑟𝑡𝑖𝑐𝑙𝑒𝑠, 𝑉𝑜𝑙. 𝐼𝐼,
        239-248, Ed. A. Septier (Academic Press, New York, 1967).

	https://archive.org/details/in.ernet.dli.2015.141778/page/n241/mode/2up

     Magnet gap correction (longitudinal fringe field)

       irho h = 1/rho [1/m]
       phi  edge angle
       gap  full gap between poles

                                     2
                   k_1*gap*h*(1 + sin (phi))
            psi = ----------------------- * (1 - k_2*k_1*gap*h*tan(phi))
                        cos phi

            k_1 is usually 1/2
            k_2 is zero here                                                  */

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
void thin_kick
(CellType &Cell, const int Order, const double MB[], const double L,
 const double h_bend, const double h_ref, ss_vect<T> &ps)
{
  // The vector potential for the combined-function sector bend is from:
  // C. Iselin "Lie Transformations and Transport Equations for Combined-
  // Function Dipoles" Part. Accel. 17, 143-155 (1985).
  int        j;
  T          BxoBrho, ByoBrho, ByoBrho1, B[3], u, p_s;
  ss_vect<T> ps0;

  const bool prt_debug = false;

  if ((h_bend != 0e0) || ((1 <= Order) && (Order <= HOMmax))) {
    ps0 = ps;
    // Compute magnetic field with Horner's rule.
    ByoBrho = MB[Order+HOMmax]; BxoBrho = MB[HOMmax-Order];
    for (j = Order-1; j >= 1; j--) {
      ByoBrho1 = ps0[x_]*ByoBrho - ps0[y_]*BxoBrho + MB[j+HOMmax];
      BxoBrho  = ps0[y_]*ByoBrho + ps0[x_]*BxoBrho + MB[HOMmax-j];
      ByoBrho  = ByoBrho1;
    }

  if (prt_debug)
    cout << scientific << setprecision(5)
	 << "\nthin_kick ->:\n" << "  h_bend = " << h_bend << " h_ref = "
	 << h_ref << "\n  BxoBrho = " << setw(13) << BxoBrho << " ByoBrho = "
	 << setw(13) << ByoBrho << "\n  ps = " << setw(13) << ps << "\n";

  if (globval.radiation || globval.emittance) {
      B[X_] = BxoBrho; B[Y_] = ByoBrho + h_bend; B[Z_] = 0e0;
      radiate(Cell, ps, L, h_ref, B);
    }

    if (h_ref != 0e0) {
      // Sector bend.
      if (true) {
	ps[px_] -= L*(ByoBrho+(h_bend-h_ref)/2e0+h_ref*h_bend*ps0[x_]
		     -h_ref*ps0[delta_]);
	ps[ct_] += L*h_ref*ps0[x_];
      } else {
	// The Hamiltonian is split into: H_d + H_k; with [H_d, H_d] = 0.
	p_s = get_p_s(ps0); u = L*h_ref*ps0[x_]/p_s;
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

  if (prt_debug)
    cout << scientific << setprecision(5)
	 << "\n<- thin_kick:\n  ps = " << setw(13) << ps << "\n";
}


template<typename T>
void EdgeFocus(const double irho, const double phi, const double gap,
		      ss_vect<T> &ps)
{
  ps[px_] += irho*tan(dtor(phi))*ps[x_];
  if (!globval.dip_edge_fudge) {
    // warning: => diverging Taylor map (see SSC-141)
    // ps[py_] -=
    //   irho*tan(dtor(phi)-get_psi(irho, phi, gap))*ps[y_]/(1e0+ps[delta_]);
    // Leading order correction.
    ps[py_] -=
      irho*tan(dtor(phi)-get_psi(irho, phi, gap))*ps[y_]*(1e0-ps[delta_]);
  } else
    ps[py_] -= irho*tan(dtor(phi)-get_psi(irho, phi, gap))*ps[y_];
}


template<typename T>
void p_rot(double phi, ss_vect<T> &ps)
{
  T          c, s, t, pz, p, val;
  ss_vect<T> ps1;

  c = cos(dtor(phi)); s = sin(dtor(phi)); t = tan(dtor(phi)); pz = get_p_s(ps);

  if (!globval.H_exact && !globval.Cart_Bend) {
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
  // Reference:
  //   E. Forest 𝐻𝑒𝑎𝑙𝑦'𝑠 𝑀𝑜𝑑𝑢𝑙𝑎𝑟 𝐴𝑝𝑝𝑟𝑜𝑎𝑐ℎ 𝑡𝑜 𝑡ℎ𝑒 𝐶𝑜𝑚𝑝𝑢𝑡𝑎𝑡𝑖𝑜𝑛 𝑜𝑓 𝑡ℎ𝑒 𝐺𝑒𝑛𝑒𝑟𝑎𝑙 𝐵𝑒𝑛𝑑𝑖𝑛𝑔
  //   𝑀𝑎𝑔𝑛𝑒𝑡 𝑀𝑎𝑝 𝐴𝑝𝑝𝑙𝑖𝑒𝑑 𝑡𝑜 𝑡ℎ𝑒 𝑄𝑢𝑎𝑑𝑟𝑎𝑡𝑖𝑐 𝑃𝑎𝑟𝑡 𝑜𝑓 𝑡ℎ𝑒 𝐻𝑎𝑚𝑖𝑙𝑡𝑜𝑛𝑖𝑎𝑛 𝑤ℎ𝑖𝑐ℎ 𝑖𝑠 𝐸𝑥𝑎𝑐𝑡 𝑖𝑛
  //   𝐷𝑒𝑙𝑡𝑎 𝑝/𝑝_0 SSC-141, 7 (1987).
  //
  //   https://lss.fnal.gov/archive/other/ssc/ssc-142.pdf
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
void quad_fringe(const double b2, ss_vect<T> &ps)
{
  T u, p_s;

  u = b2/(12e0*(1e0+ps[delta_])); p_s = u/(1e0+ps[delta_]);
  ps[py_] /= 1e0 - 3e0*u*sqr(ps[y_]); ps[y_] -= u*cube(ps[y_]);
  if (globval.Cavity_on) ps[ct_] -= p_s*cube(ps[y_])*ps[py_];
  ps[px_] /= 1e0 + 3e0*u*sqr(ps[x_]);
  if (globval.Cavity_on) ps[ct_] += p_s*cube(ps[x_])*ps[px_];
  ps[x_] += u*cube(ps[x_]); u = u*3e0; p_s = p_s*3e0;
  ps[y_] = exp(-u*sqr(ps[x_]))*ps[y_]; ps[py_] = exp(u*sqr(ps[x_]))*ps[py_];
  ps[px_] += 2e0*u*ps[x_]*ps[y_]*ps[py_];
  if (globval.Cavity_on) ps[ct_] -= p_s*sqr(ps[x_])*ps[y_]*ps[py_];
  ps[x_] = exp(u*sqr(ps[y_]))*ps[x_]; ps[px_] = exp(-u*sqr(ps[y_]))*ps[px_];
  ps[py_] -= 2e0*u*ps[y_]*ps[x_]*ps[px_];
  if (globval.Cavity_on) ps[ct_] += p_s*sqr(ps[y_])*ps[x_]*ps[px_];
}


void get_dI_eta_5(CellType &Cell)
{
  double       L, K, h, b2, alpha, beta, gamma, psi, eta, etap;
  ss_vect<tps> Id;
  CellType     *Cellp;

  Id.identity();

  L = Cell.Elem.PL;
  h = Cell.Elem.M->Pirho;
  b2 = Cell.Elem.M->PBpar[Quad+HOMmax];
  K = b2 + sqr(Cell.Elem.M->Pirho);
  psi = sqrt(fabs(K))*L;
  Cellp = &Cell - 1;
  alpha = Cellp->Alpha[X_]; beta = Cellp->Beta[X_];
  gamma = (1e0+sqr(alpha))/beta;
  eta = Cellp->Eta[X_]; etap = Cellp->Etap[X_];

  Cell.dI[1] = L*eta*h;
  Cell.dI[2] = L*sqr(h);
  Cell.dI[3] = L*fabs(cube(h));

  if (K > 0e0) {
    Cell.dI[4] =
      h/K*(2e0*b2+sqr(h))
      *((eta*sqrt(K)*sin(psi)+etap*(1e0-cos(psi)))
	+ h/sqrt(K)*(psi-sin(psi)));

    Cell.dI[5] =
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

    Cell.dI[4] =
      h/K*(2e0*b2+sqr(h))
      *((eta*sqrt(K)*sinh(psi)-etap*(1e0-cosh(psi)))
	- h/sqrt(K)*(psi-sinh(psi)));

    Cell.dI[5] =
      L*fabs(cube(h))
      *(gamma*sqr(eta)+2e0*alpha*eta*etap+beta*sqr(etap))
      + 2e0*pow(h, 4)/(pow(K, 3e0/2e0))
      *(sqrt(K)*(alpha*eta+beta*etap)*(cosh(psi)-1e0)
	+(gamma*eta+alpha*etap)*(psi-sinh(psi)))
      + fabs(pow(h, 5))/(4e0*pow(K, 5e0/2e0))
      *(2e0*alpha*sqrt(K)*(4e0*cosh(psi)-cosh(2e0*psi)-3e0)
	-beta*K*(2e0*psi-sinh(2e0*psi))
	+gamma*(6e0*psi-8e0*sinh(psi)+sinh(2e0*psi)));
  }

  // For completeness.
  Cell.curly_dH_x = Cell.dI[5]/(L*fabs(cube(h)));

}


template<typename T>
void Mpole_Pass(CellType &Cell, ss_vect<T> &ps)
{
  int             seg = 0, i;
  double          dL = 0e0, dL1 = 0e0, dL2 = 0e0;
  double          dkL1 = 0e0, dkL2 = 0e0, h_ref = 0e0;
  elemtype        *elemp;
  MpoleType       *M;

  elemp = &Cell.Elem; M = elemp->M;

  GtoL(ps, Cell.dS, Cell.dT, M->Pc0, M->Pc1, M->Ps1);

  if (globval.emittance) {
    // Needs A^-1.
    Cell.curly_dH_x = 0e0;
    for (i = 0; i <= 5; i++)
      Cell.dI[i] = 0e0;
  }

  switch (M->Pmethod) {

  case Meth_Fourth:
    if (globval.mat_meth && (M->Porder <= Quad)) {
      ps = is_double< ss_vect<T> >::ps(M->M_lin*ps);

      if (globval.emittance)
	if ((Cell.Elem.PL != 0e0) && (Cell.Elem.M->Pirho != 0e0))
	  get_dI_eta_5(Cell);
    } else {
      // Fringe fields.
      if (globval.quad_fringe && (M->PB[Quad+HOMmax] != 0e0))
	quad_fringe(M->PB[Quad+HOMmax], ps);
      if (!globval.Cart_Bend) {
	if (M->Pirho != 0e0) EdgeFocus(M->Pirho, M->PTx1, M->Pgap, ps);
      } else {
	p_rot(M->PTx1, ps); bend_fringe(M->Pirho, ps);
      }

      if (M->Pthick == thick) {
	if (!globval.Cart_Bend) {
	  // Polar coordinates.
	  h_ref = M->Pirho; dL = elemp->PL/M->PN;
	} else {
	  // Cartesian coordinates.
	  h_ref = 0e0;
	  if (M->Pirho == 0e0)
	    dL = elemp->PL/M->PN;
	  else
	    dL = 2e0/M->Pirho*sin(elemp->PL*M->Pirho/2e0)/M->PN;
	}

	dL1 = c_1*dL; dL2 = c_2*dL; dkL1 = d_1*dL; dkL2 = d_2*dL;

	for (seg = 1; seg <= M->PN; seg++) {
	  if (globval.emittance) {
	    // Needs A^-1.
	    Cell.curly_dH_x += is_tps<tps>::get_curly_H(ps);
	    Cell.dI[4] += is_tps<tps>::get_dI_eta(ps);
	  }

	  Drift(dL1, ps);
	  thin_kick(Cell, M->Porder, M->PB, dkL1, M->Pirho, h_ref, ps);
	  Drift(dL2, ps);
	  thin_kick(Cell, M->Porder, M->PB, dkL2, M->Pirho, h_ref, ps);

	  if (globval.emittance) {
	    // Needs A^-1.
	    Cell.curly_dH_x += 4e0*is_tps<tps>::get_curly_H(ps);
	    Cell.dI[4] += 4e0*is_tps<tps>::get_dI_eta(ps);
	  }

	  Drift(dL2, ps);
	  thin_kick(Cell, M->Porder, M->PB, dkL1, M->Pirho, h_ref, ps);
	  Drift(dL1, ps);

	  if (globval.emittance) {
	    // Needs A^-1.
	    Cell.curly_dH_x += is_tps<tps>::get_curly_H(ps);
	    Cell.dI[4] += is_tps<tps>::get_dI_eta(ps);
	  }
	}

	if (globval.emittance) {
	  // Needs A^-1.
	  Cell.curly_dH_x /= 6e0*M->PN;
	  Cell.dI[1] += elemp->PL*is_tps<tps>::get_dI_eta(ps)*M->Pirho;
	  Cell.dI[2] += elemp->PL*sqr(M->Pirho);
	  Cell.dI[3] += elemp->PL*fabs(cube(M->Pirho));
	  Cell.dI[4] *=
	    elemp->PL*M->Pirho*(sqr(M->Pirho)+2e0*M->PBpar[Quad+HOMmax])
	    /(6e0*M->PN);
	  Cell.dI[5] += elemp->PL*fabs(cube(M->Pirho))*Cell.curly_dH_x;
	}
      } else
	thin_kick(Cell, M->Porder, M->PB, 1e0, 0e0, 0e0, ps);

      // Fringe fields.
      if (!globval.Cart_Bend) {
	if (M->Pirho != 0e0) EdgeFocus(M->Pirho, M->PTx2, M->Pgap, ps);
      } else {
	bend_fringe(-M->Pirho, ps); p_rot(M->PTx2, ps);
      }
      if (globval.quad_fringe && (M->PB[Quad+HOMmax] != 0e0))
	quad_fringe(-M->PB[Quad+HOMmax], ps);
    }
    break;

  default:
    printf("Mpole_Pass: Method not supported %10s %d\n",
	   Cell.Elem.PName, M->Pmethod);
    exit_(0);
    break;
  }

  LtoG(ps, Cell.dS, Cell.dT, M->Pc0, M->Pc1, M->Ps1);
}


template<typename T>
void Marker_Pass(CellType &Cell, ss_vect<T> &ps)
{
  GtoL(ps, Cell.dS, Cell.dT, 0e0, 0e0, 0e0);

  if (globval.emittance && !globval.Cavity_on)
    // Needs A^-1.
    Cell.curly_dH_x = is_tps<tps>::get_curly_H(ps);

  LtoG(ps, Cell.dS, Cell.dT, 0e0, 0e0, 0e0);
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
void Cav_Pass(const CellType &Cell, ss_vect<T> &ps)
{
  const elemtype*   elemp = &Cell.Elem;
  const CavityType* C     = elemp->C;
  const double      L     = elemp->PL;

  T delta;

  Drift(L/2e0, ps);
  if (globval.Cavity_on && C->V_RF != 0e0) {
    delta = -C->V_RF/(globval.Energy*1e9)
            *sin(2e0*M_PI*C->f_RF/c0*ps[ct_]-C->phi_RF);
    ps[delta_] += delta;

    if (globval.radiation) globval.dE -= is_double<T>::cst(delta);

    if (globval.pathlength) ps[ct_] -= C->harm_num/C->f_RF*c0;
  }
  Drift(L/2e0, ps);
}

#else

template<typename T>
void Cav_Pass1(const CellType &Cell, ss_vect<T> &ps)
{
  /* J. Rosenzweig and L. Serafini "Transverse Particle Motion in
     Radio-Frequency Linear Accelerators" Phys. Rev. E 49(2),
     1599-1602 (1994).                                                        */

  elemtype   *elemp;
  CavityType *C;
  int        k;
  double     L, h, p_t1;
  T          delta_max, ddelta, delta;

  elemp = &Cell.Elem; C = elemp->C; L = elemp->PL;

  h = L/(C->PN+1e0);
  // globval.Energy contains p_0.
  delta_max = C->Pvolt/(1e9*globval.Energy);
  ddelta = delta_max/C->PN;
  delta = delta_max*sin(2e0*M_PI*C->Pfreq*ps[ct_]/c0+C->phi);
  if (C->entry_focus) Cav_Focus(L, delta, true, ps);
  for (k = 0; k < C->PN; k++) {
    Drift(h, ps);

    ps[delta_] -= ddelta*sin(2e0*M_PI*C->Pfreq*(ps[ct_]-k*h)/c0-C->phi);

    if (globval.radiation) globval.dE -= is_double<T>::cst(ddelta);
    if (globval.pathlength) ps[ct_] -= C->Ph/C->Pfreq*c0;
  }
  Drift(h, ps);
  if (C->exit_focus) Cav_Focus(L, delta, false, ps);

  if (false) {
    // Update p_0.
    p_t1 = is_double<T>::cst(ps[delta_]);
    ps[delta_] -= p_t1;
    // globval.Energy contains p_0.
    globval.Energy *= sqrt(1e0+2e0*p_t1/globval.beta0+sqr(p_t1));
    globval.gamma0 = sqrt(sqr(m_e)+sqr(1e9*globval.Energy))/m_e;
    globval.beta0  = sqrt(1e0-1e0/sqr(globval.gamma0));
    printf("\np0 = %12.5e, beta0 = %12.5e, gamma0 = %12.5e\n",
	   globval.Energy, globval.beta0, globval.gamma0);
  }
}


template<typename T>
void Cav_Pass(const CellType &Cell, ss_vect<T> &ps)
{
  /* J. Rosenzweig and L. Serafini "Transverse Particle Motion in
     Radio-Frequency Linear Accelerators" Phys. Rev. E 49(2),
     1599-1602 (1994).                                                        */

  elemtype   *elemp;
  CavityType *C;
  double     L, Lambda, phi;
  double     dgammaMax, dgamma, gamma, gamma1;
  double     sf1, f2, f2s;
  double     f5, sf5, dpr, dct, p_t1;
  double     p0, p_t, delta, alpha, dp;
  ss_vect<T> ps0;

  const bool RandS = false;
 
  elemp = &Cell.Elem; C = elemp->C; L = elemp->PL; phi = C->phi;
  Lambda = c0/C->Pfreq;

  p_t = is_double<T>::cst(ps[delta_]);
  delta = sqrt(1e0+2e0*p_t/globval.beta0+sqr(p_t)) - 1e0;
  // globval.Energy contains p_0 [GeV].
  p0 = 1e9*globval.Energy/m_e;
  gamma = sqrt(1e0+sqr(p0));
  dgammaMax = C->Pvolt/m_e; dgamma = dgammaMax*sin(phi);
  gamma1 = gamma + dgamma;
  dp = sqrt(sqr(gamma1)-1e0) - p0;

  printf("delta = %e, p0 = %e, gamma = %e, gamma1 = %e, dp = %e\n",
	 delta, p0, gamma, gamma1, dp);

  if (C->entry_focus) Cav_Focus(L, dgamma/gamma, true, ps);

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

   // globval.Energy contains p_0 [GeV].
    ps[delta_] =
      2e0*M_PI*C->Pfreq*dgammaMax*cos(phi)/(c0*gamma1)*ps0[ct_]
      + 1e0/(1e0+dp/p0)*ps0[delta_];
  }

  if (C->exit_focus) Cav_Focus(L, dgamma/(gamma+dgamma), false, ps);

  if (false) {
    // Update p_0.
    p_t1 = is_double<T>::cst(ps[delta_]);
    ps[delta_] -= p_t1;
    // globval.Energy contains p_0.
    globval.Energy *= sqrt(1e0+2e0*p_t1/globval.beta0+sqr(p_t1));
    globval.gamma0 = sqrt(sqr(m_e)+sqr(p0))/m_e;
    globval.beta0  = sqrt(1e0-1e0/sqr(globval.gamma0));
    printf("\np0 = %12.5e, beta0 = %12.5e, gamma0 = %12.5e\n",
	   globval.Energy, globval.beta0, globval.gamma0);
  }
}

#endif


void get_dI_eta_5_ID(CellType &Cell)
{
  double       L, K, h, b2, alpha, beta, gamma, psi, eta, etap;
  ss_vect<tps> Id;
  CellType     *Cellp;


  Id.identity();

  L = Cell.Elem.PL;
  h = Cell.Elem.M->Pirho;
  b2 = Cell.Elem.M->PBpar[Quad+HOMmax];
  K = b2 + sqr(Cell.Elem.M->Pirho);
  psi = sqrt(fabs(K))*L;
  Cellp = &Cell - 1;
  alpha = Cellp->Alpha[X_]; beta = Cellp->Beta[X_];
  gamma = (1e0+sqr(alpha))/beta;
  eta = Cellp->Eta[X_]; etap = Cellp->Etap[X_];

  Cell.dI[1] += L*eta*h;
  Cell.dI[2] += L*sqr(h);
  Cell.dI[3] += L*fabs(cube(h));

  if (K > 0e0) {
    Cell.dI[4] +=
      h/K*(2e0*b2+sqr(h))
      *((eta*sqrt(K)*sin(psi)+etap*(1e0-cos(psi)))
	+ h/sqrt(K)*(psi-sin(psi)));

    Cell.dI[5] +=
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

    Cell.dI[4] +=
      h/K*(2e0*b2+sqr(h))
      *((eta*sqrt(K)*sinh(psi)-etap*(1e0-cosh(psi)))
	- h/sqrt(K)*(psi-sinh(psi)));

    Cell.dI[5] +=
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
inline void get_Axy(const WigglerType *W, const double z,
		    ss_vect<T> &x, T AxoBrho[], T AyoBrho[])

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

    if (globval.radiation) {
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

  if (globval.radiation) {
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
void Wiggler_pass_EF(CellType &Cell, ss_vect<T> &ps)
{
  // First order symplectic integrator for wiggler using expanded Hamiltonian

  int
    i, nstep = 0;
  double
    h, z;
  T
    AxoBrho[4] = {0e0, 0e0, 0e0, 0e0}, AyoBrho[4] = {0e0, 0e0, 0e0, 0e0},
    psi, hodp, a12, a21, a22, det,
    d1, d2, a11, c11, c12, c21, c22, x2, B[3];
  elemtype
    *elemp;

  elemp = &Cell.Elem;

  switch (elemp->Pkind) {
  case Wigl:
    nstep = elemp->W->PN;
    break;
  case FieldMap:
    nstep = elemp->FM->n_step;
    break;
  default:
    std::cout << "Wiggler_pass_EF: unknown element type" << std::endl;
    exit_(1);
    break;
  }

  h = elemp->PL/nstep; z = 0e0;
  for (i = 1; i <= nstep; ++i) {
    switch (elemp->Pkind) {
    case Wigl:
      get_Axy(elemp->W, z, ps, AxoBrho, AyoBrho);
      break;
    case FieldMap:
//      get_Axy_map(elemp->FM, z, ps, AxoBrho, AyoBrho);
      break;
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

    if (globval.pathlength) ps[ct_] += h;

    if (globval.radiation || globval.emittance) {
      B[X_] = -AyoBrho[3]; B[Y_] = AxoBrho[3]; B[Z_] = AyoBrho[1] - AxoBrho[2];
      radiate(Cell, ps, h, 0e0, B);
    }

    z += h;
  }
}


template<typename T>
inline void get_Axy2
(const double z, const double kxV, const double kxH, const double kz,
 const double BoBrhoV, const double BoBrhoH, const double phi,
 ss_vect<T> &x, T AxoBrho[], T AyoBrho[])
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

  if (globval.radiation) {
    cz1 = cos(kz*z); cz2=cos(kz*z+phi);
    /* derivatives with respect to z */
    AxoBrho[3] += BoBrhoV*cx*chy*cz1;
    AxoBrho[3] -= BoBrhoH*kxH/kyH*shx*sy*cz2;
    AyoBrho[3] += BoBrhoV*kxV/kyV*sx*shy*cz1;
    AyoBrho[3] -= BoBrhoH*chx*cy*cz2;
  }
}


template<typename T>
void Wiggler_pass_EF2
(CellType &Cell, int nstep, double L, double kxV, double kxH, double kz,
 double BoBrhoV, double BoBrhoH, double phi, ss_vect<T> &x)
{
  // First order symplectic integrator for wiggler using expanded Hamiltonian

  int    i;
  double h, z;
  T      hodp, B[3], px1, px2, px3, py1, py2, AxoBrho[4], AyoBrho[4], psi;
  T      px = 0e0, py = 0e0;

  h = L/nstep; z = 0e0;
  for (i = 1; i <= nstep; ++i) {
    get_Axy2(z, kxV, kxH, kz, BoBrhoV, BoBrhoH, phi, x, AxoBrho, AyoBrho);

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

    if (globval.pathlength) x[ct_] += h;

    if (globval.radiation || globval.emittance) {
      B[X_] = -AyoBrho[3]; B[Y_] = AxoBrho[3]; B[Z_] = AyoBrho[1] - AxoBrho[2];
      radiate(Cell, x, h, 0e0, B);
    }

    z += h;
  }

  x[px_] = px; x[py_] = py;
}


template<typename T>
inline void get_Axy_EF3
(const WigglerType *W, const double z, const ss_vect<T> &ps, T &AoBrho,
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

      if (globval.radiation || (globval.emittance && !globval.Cavity_on)) {
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

      if (globval.radiation || (globval.emittance && !globval.Cavity_on)) {
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
void Wiggler_pass_EF3(CellType &Cell, ss_vect<T> &ps)
{
  /* Symplectic integrator (2nd order) for Insertion Devices based on:

       E. Forest, et al "Explicit Symplectic Integrator for s-dependent
       Static Magnetic Field" Phys. Rev. E 68,  046502 (2003)                 */

  int         i;
  double      h, z, irho, curly_dH_x;
  T           hd, AxoBrho, AyoBrho, dAxoBrho[3], dAyoBrho[3], dpy, dpx, B[3];
  ss_vect<T>  ps1;
  elemtype    *elemp;
  WigglerType *W;

  elemp = &Cell.Elem; W = elemp->W;

  h = elemp->PL/W->PN; z = 0e0;

  if (globval.emittance && !globval.Cavity_on) {
    // Needs A^-1.
    Cell.curly_dH_x = 0e0;
    for (i = 0; i <= 5; i++)
      Cell.dI[i] = 0e0;
  }

  for (i = 1; i <= W->PN; i++) {
    hd = h/(1e0+ps[delta_]);

    // 1: half step in z
    z += 0.5*h;

    // 2: half drift in y
    get_Axy_EF3(W, z, ps, AyoBrho, dAyoBrho, dpx, false);

    ps[px_] -= dpx; ps[py_] -= AyoBrho;
    ps[y_] += 0.5*hd*ps[py_];
    ps[ct_] += sqr(0.5)*hd*sqr(ps[py_])/(1e0+ps[delta_]);

    get_Axy_EF3(W, z, ps, AyoBrho, dAyoBrho, dpx, false);

    ps[px_] += dpx; ps[py_] += AyoBrho;

    // 3: full drift in x
    get_Axy_EF3(W, z, ps, AxoBrho, dAxoBrho, dpy, true);

    ps[px_] -= AxoBrho; ps[py_] -= dpy; ps[x_] += hd*ps[px_];
    ps[ct_] += 0.5*hd*sqr(ps[px_])/(1e0+ps[delta_]);

    if (globval.pathlength) ps[ct_] += h;

    get_Axy_EF3(W, z, ps, AxoBrho, dAxoBrho, dpy, true);

    ps[px_] += AxoBrho; ps[py_] += dpy;

    // 4: a half drift in y
    get_Axy_EF3(W, z, ps, AyoBrho, dAyoBrho, dpx, false);

    ps[px_] -= dpx; ps[py_] -= AyoBrho;
    ps[y_] += 0.5*hd*ps[py_];
    ps[ct_] += sqr(0.5)*hd*sqr(ps[py_])/(1e0+ps[delta_]);

    get_Axy_EF3(W, z, ps, AyoBrho, dAyoBrho, dpx, false);

    ps[px_] += dpx; ps[py_] += AyoBrho;

    // 5: half step in z
    z += 0.5*h;

    if (globval.radiation || globval.emittance) {
      get_Axy_EF3(W, z, ps, AyoBrho, dAyoBrho, dpx, false);
      get_Axy_EF3(W, z, ps, AxoBrho, dAxoBrho, dpy, true);
      B[X_] = -dAyoBrho[Z_]; B[Y_] = dAxoBrho[Z_];
      B[Z_] = dAyoBrho[X_] - dAxoBrho[Y_];
      // Tranform from Conjugate to Kinematic Momenta.
      ps[px_] -= AxoBrho; ps[py_] -= AyoBrho;
      radiate(Cell, ps, h, 0e0, B);
      // Tranform from Kinematic to Conjugate Momenta.
      ps[px_] += AxoBrho; ps[py_] += AyoBrho;
    }

    if (globval.emittance && !globval.Cavity_on) {
      // Needs A^-1.
      get_Axy_EF3(W, z, ps, AxoBrho, dAxoBrho, dpy, true);
      irho = is_double<T>::cst(dAxoBrho[Z_]);
      // Tranform to Configuration Space.
      ps1 = ps;
      ps1[px_] = (ps[px_]-AxoBrho)/(1e0+ps[delta_]);
      curly_dH_x = is_tps<tps>::get_curly_H(ps1);
      Cell.curly_dH_x += curly_dH_x;
      Cell.dI[1] += is_tps<tps>::get_dI_eta(ps)*irho;
      Cell.dI[2] += sqr(irho);
      Cell.dI[3] += fabs(cube(irho));
      Cell.dI[4] += is_tps<tps>::get_dI_eta(ps)*cube(irho);
      Cell.dI[5] += curly_dH_x*fabs(cube(irho));
    }
  }

  if (globval.emittance && !globval.Cavity_on) {
    // Needs A^-1.
    Cell.curly_dH_x /= W->PN;
    for (i = 0; i <= 5; i++)
      Cell.dI[i] *= elemp->PL/W->PN;
  }
}


template<typename T>
void Wiggler_Pass(CellType &Cell, ss_vect<T> &ps)
{
  int         seg;
  double      L, L1, L2, K1, K2;
  elemtype    *elemp;
  WigglerType *W;
  ss_vect<T>  ps1;

  elemp = &Cell.Elem; W = elemp->W;
  // Global -> Local
  GtoL(ps, Cell.dS, Cell.dT, 0e0, 0e0, 0e0);
  switch (W->Pmethod) {

  case Meth_Linear:
    std::cout << "Wiggler_Pass: Meth_Linear not supported" << std::endl;
    exit_(1);
    break;

  case Meth_First:
    if ((W->BoBrhoV[0] != 0e0) || (W->BoBrhoH[0] != 0e0)) {
      if (!globval.EPU)
	Wiggler_pass_EF(Cell, ps);
      else {
	Wiggler_pass_EF2(Cell, W->PN, elemp->PL, W->kxV[0], W->kxH[0],
		2e0*M_PI/W->Lambda, W->BoBrhoV[0], W->BoBrhoH[0],
		W->phi[0], ps);
      }
    } else
      // drift if field = 0
      Drift(elemp->PL, ps);
    break;

  case Meth_Second:
    if ((W->BoBrhoV[0] != 0e0) || (W->BoBrhoH[0] != 0e0)) {
      Wiggler_pass_EF3(Cell, ps);
    } else
      // drift if field = 0
      Drift(elemp->PL, ps);
    break;

  case Meth_Fourth:  /* 4-th order integrator */
    L = elemp->PL/W->PN;
    L1 = c_1*L; L2 = c_2*L; K1 = d_1*L; K2 = d_2*L;
    for (seg = 1; seg <= W->PN; seg++) {
      Drift(L1, ps); ps1 = ps;
      thin_kick(Cell, W->Porder, W->PBW, K1, 0e0, 0e0, ps1);
      ps[py_] = ps1[py_]; Drift(L2, ps); ps1 = ps;
      thin_kick(Cell, W->Porder, W->PBW, K2, 0e0, 0e0, ps1);
      ps[py_] = ps1[py_]; Drift(L2, ps); ps1 = ps;
      thin_kick(Cell, W->Porder, W->PBW, K1, 0e0, 0e0, ps1);
      ps[py_] = ps1[py_]; Drift(L1, ps);
    }
    break;
  }
  // Local -> Global
  LtoG(ps, Cell.dS, Cell.dT, 0e0, 0e0, 0e0);
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
void radiate_cs(CellType &Cell, ss_vect<T> &cs, const double L, const T B[])
{
  // M. Sands "The Physics of Electron Storage Rings" SLAC-121, p. 98.
  // ddelta/d(ds) = -C_gamma*E_0^3*(1+delta)^2*(B_perp/(Brho))^2/(2*pi)
  T  p_s, ds, B2_perp = 0e0, B2_par = 0e0;

  // Large ring: x' and y' unchanged, ds = -p_s*L.
  ds = (1e0+(sqr(cs[px_])+sqr(cs[py_]))/2e0)*L;
  get_B2(0e0, B, cs, B2_perp, B2_par);

  p_s = get_p_s_cs(cs);

  if (globval.radiation)
    cs[delta_] -= cl_rad*sqr(p_s)*B2_perp*ds;

  if (globval.emittance)
    is_tps<T>::emittance(Cell, B2_perp, ds, p_s, cs);
}


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
void rk4_
(CellType &Cell, const ss_vect<T> &y, const ss_vect<T> &dydx, const double x,
 const double h, ss_vect<T> &cs, const double z,
 void (*derivs)(const CellType &, const double, const ss_vect<T> &,
		ss_vect<T> &))
{
  int        j;
  double     xh, hh, h6;
  T          BoBrho[3], p_s;
  ss_vect<T> dym, dyt, yt;

  hh = h*0.5; h6 = h/6e0;
  xh = x + hh; yt = y + hh*dydx;
  (*derivs)(Cell, xh, yt, dyt); yt = y + hh*dyt;
  (*derivs)(Cell, xh, yt, dym); yt = y + h*dym; dym += dyt;
  (*derivs)(Cell, x+h, yt, dyt);
  cs = y + h6*(dydx+dyt+2e0*dym);

  if (globval.radiation || globval.emittance) {
    if (!get_BoBrho(Cell.Elem.FM, z, cs, BoBrho)) {
      for (j = 0; j < ss_dim; j++)
	cs[j] = NAN;
      return;
    }

    radiate_cs(Cell, cs, h, BoBrho);
  }
}


template<typename T>
void f_FM(const CellType &Cell, const double z, const ss_vect<T> &cs,
	  ss_vect<T> &Dcs)
{
  // Coordinates are: [x, x', y, y', -ct, delta].

  int j;
  T   BoBrho[3], p_s;

  if (!get_BoBrho(Cell.Elem.FM, z, cs, BoBrho)) {
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

  Dcs[ct_] = (1e0+cs[delta_])/p_s - ((!globval.pathlength)? 1e0 : 0e0);

  Dcs[delta_] = 0e0;
}


template<typename T>
void FieldMap_pass_RK(CellType &Cell, ss_vect<T> &ps)
{
  int          i;
  double       h, z;
  T            p_s;
  ss_vect<T>   Dps;
  FieldMapType *FM;

  const int n_step = 2; // Each step needs: f(z_n), f(z_n+h), f(z_n+2h)

  FM = Cell.Elem.FM;

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
  p_s = get_p_s(ps); ps[px_] /= p_s; ps[py_] /= p_s;

  h = n_step*FM->dx[Z_]; z = FM->x[Z_][1]; FM->Lr = 0e0;
  if (trace)
    outf_ << std::scientific << std::setprecision(3)
	  << std::setw(5) << 0 << std::setw(11) << s_FM
	  << std::setw(11) << is_double< ss_vect<T> >::cst(ps) << "\n";
  for(i = 1+FM->cut; i < FM->n[Z_]-FM->cut; i += n_step) {
    if (i <= FM->n[Z_]-FM->cut-2) {
      f_FM(Cell, z, ps, Dps);

      if (Dps[x_] == NAN) {
	std::cout << "FieldMap_pass_RK: particle lost" << std::endl;
	std::cout << ps;
	return;
      }

      rk4_(Cell, ps, Dps, FM->x[Z_][i], h, ps, z, f_FM);

      z += h; FM->Lr += h; s_FM += h;
    } else {
      // Use 2nd order Runge-Kutta (aka Midpoint Method).
      f_FM(Cell, z, ps, Dps);

      if (Dps[x_] == NAN) {
	std::cout << "FieldMap_pass_RK: particle lost" << std::endl;
	std::cout << ps;
	return;
      }

      ps += h/2e0*Dps;

      z += h/2e0; FM->Lr += h/2e0; s_FM += h/2e0;
    }

    if (trace)
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
void FieldMap_pass_SI(CellType &Cell, ss_vect<T> &ps)
{
  /* E. Chacon-Golcher, F. Neri "A Symplectic Integrator with Arbitrary
     Vector and Scalar Potentials" Phys. Lett. A 372 p. 4661-4666 (2008).    */

  int          i, j = 0;
  double       h, z;
  T            hd, AoBrho[2], dAoBrho[2], AoBrho_int;
  ss_vect<T>   ps1;
  FieldMapType *FM;

  const int    n_step = 2;
  const double d_diff = 1e0;


  FM = Cell.Elem.FM;

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
  if (trace)
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

    if (globval.pathlength) ps[ct_] += h;

    FM->Lr += h;

    if (globval.radiation || globval.emittance) {
//      B[X_] = -AoBrhoy[3]; B[Y_] = AoBrho[X_][3];
//      B[Z_] = AoBrhoy[1] - AoBrho[X_][2];
//      radiate(Cell, ps, h, 0e0, B);
    }

    if (trace)
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
template void f_FM(const CellType &, const double, const ss_vect<double> &,
		   ss_vect<double> &);
template void f_FM(const CellType &, const double, const ss_vect<tps> &,
		   ss_vect<tps> &);
template void rk4_
(CellType &, const ss_vect<double> &, const ss_vect<double> &, const double,
 const double, ss_vect<double> &, const double,
 void (*derivs)
 (const CellType &, const double, const ss_vect<double> &, ss_vect<double> &));
template void rk4_
(CellType &, const ss_vect<tps> &, const ss_vect<tps> &, const double,
 const double, ss_vect<tps> &, const double,
 void (*derivs)
 (const CellType &, const double, const ss_vect<tps> &, ss_vect<tps> &));
template void FieldMap_pass_RK(CellType &, ss_vect<double> &);
template void FieldMap_pass_RK(CellType &, ss_vect<tps> &);
template void FieldMap_pass_SI(CellType &, ss_vect<double> &);
template void FieldMap_pass_SI(CellType &, ss_vect<tps> &);


template<typename T>
void FieldMap_Pass(CellType &Cell, ss_vect<T> &ps)
{
  int          k;
  double       Ld;
  FieldMapType *FM;

  if (trace & first_FM) {
    file_wr(outf_, "FieldMap_pass.dat");
    s_FM = 0e0;
    first_FM = false;
  }

  FM = Cell.Elem.FM;

//  GtoL(ps, Cell.dS, Cell.dT, 0e0, 0e0, 0e0);

  Ld = (FM->Lr-Cell.Elem.PL)/2e0;
  p_rot(FM->phi/2e0*180e0/M_PI, ps);
  printf("\nFieldMap_Pass:\n");
  printf("  phi = %12.5e\n  cut = %12d\n", FM->phi, FM->cut);
  printf("  entrance negative drift [m] %12.5e\n", -Ld);
  Drift(-Ld, ps);

  // n_step: number of Field Map repetitions.
  for (k = 1; k <= FM->n_step; k++) {
    if (sympl)
      FieldMap_pass_SI(Cell, ps);
    else
      FieldMap_pass_RK(Cell, ps);
  }

  printf("  exit negative drift [m]     %12.5e\n", -Ld);
  Drift(-Ld, ps);
  p_rot(FM->phi/2e0*180e0/M_PI, ps);

//  LtoG(ps, Cell.dS, Cell.dT, 0e0, 0e0, 0e0);

//  outf_.close();
}


template<typename T>
void Insertion_Pass(CellType &Cell, ss_vect<T> &x)
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

  elemtype *elemp;
  double   LN = 0e0;
  T        tx2, tz2;      /* thetax and thetaz retrieved from
			     interpolation routine */
  T        d, B2_perp;
  double   alpha0 = 0e0;  // 1/ brh0
  double   alpha02= 0e0;  // alpha square
  int      Nslice = 0;
  int      i = 0;
  bool     outoftable = false;

  elemp  = &Cell.Elem; Nslice = elemp->ID->PN;

  if (elemp->ID->linear) {
    alpha0 = c0/globval.Energy*1E-9*elemp->ID->scaling;
    alpha02 = sgn(elemp->ID->scaling)*alpha0*alpha0;
  } else
    alpha02 = 1e-6*elemp->ID->scaling;

//  /* Global -> Local */
//  GtoL(X, Cell->dS, Cell->dT, 0e0, 0e0, 0e0);

  p_rot(elemp->ID->phi/2e0*180e0/M_PI, x);

  // (Nslice+1) drifts, nslice kicks
  // LN = elemp->PL/(Nslice+1);

  // Nslice drifts and kicks.
  LN = elemp->PL/Nslice;
  Drift(LN/2e0, x);

  for (i = 1; i <= Nslice; i++) {
    // printf("%3d %2d %2d %5.3f %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e\n",
    // 	   i, elemp->ID->linear, elemp->ID->secondorder, elemp->ID->scaling,
    // 	   is_double<T>::cst(x[x_]), is_double<T>::cst(x[px_]),
    // 	   is_double<T>::cst(x[y_]), is_double<T>::cst(x[py_]),
    // 	   is_double<T>::cst(x[delta_]), is_double<T>::cst(x[ct_]));
    // Second order kick
    if (elemp->ID->secondorder){
      // if (!elemp->ID->linear)
      //   SplineInterpolation2(x[x_], x[y_], tx2, tz2, Cell, outoftable);
      // else {
        LinearInterpolation2(x[x_], x[y_], tx2, tz2, B2_perp, Cell,
			     outoftable, 2);

	// Scale locally with (Brho) (as above) instead of when the file
	// is read; since the beam energy might not be known at that time.
	if (globval.radiation || globval.emittance)
	  radiate_ID(Cell, x, LN, elemp->ID->scaling*B2_perp);
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

  p_rot(elemp->ID->phi/2e0*180e0/M_PI, x);

//  CopyVec(6L, x, Cell->BeamPos);

//  /* Local -> Global */
//  LtoG(X, Cell->dS, Cell->dT, 0e0, 0e0, 0e0);
}


template<typename T>
void sol_pass(CellType &Cell, ss_vect<T> &x)
{
  int      i;
  double   h, z;
  T        hd, AxoBrho, AyoBrho, dAxoBrho[3], dAyoBrho[3], dpy, dpx, B[3];
  elemtype *elem;

  elem = &Cell.Elem;

  h = elem->PL/elem->Sol->N; z = 0e0;

  for (i = 1; i <= elem->Sol->N; i++) {
    hd = h/(1e0+x[delta_]);

    // 1: half step in z
    z += 0.5*h;

    // 2: half drift in y
    AyoBrho = elem->Sol->BoBrho*x[x_]/2e0; dpx = elem->Sol->BoBrho*x[y_]/2e0;
//    get_Axy_EF3(elem->W, z, x, AyoBrho, dAyoBrho, dpx, false);

    x[px_] -= dpx; x[py_] -= AyoBrho;
    x[y_] += 0.5*hd*x[py_];
    x[ct_] += sqr(0.5)*hd*sqr(x[py_])/(1e0+x[delta_]);

    AyoBrho = elem->Sol->BoBrho*x[x_]/2e0; dpx = elem->Sol->BoBrho*x[y_]/2e0;
//    get_Axy_EF3(elem->W, z, x, AyoBrho, dAyoBrho, dpx, false);

    x[px_] += dpx; x[py_] += AyoBrho;

    // 3: full drift in x
    AxoBrho = -elem->Sol->BoBrho*x[y_]/2e0; dpy = -elem->Sol->BoBrho*x[x_]/2e0;
//    get_Axy_EF3(elem->W, z, x, AxoBrho, dAxoBrho, dpy, true);

    x[px_] -= AxoBrho; x[py_] -= dpy; x[x_] += hd*x[px_];
    x[ct_] += 0.5*hd*sqr(x[px_])/(1e0+x[delta_]);

    if (globval.pathlength) x[ct_] += h;

    AxoBrho = -elem->Sol->BoBrho*x[y_]/2e0; dpy = -elem->Sol->BoBrho*x[x_]/2e0;
//    get_Axy_EF3(elem->W, z, x, AxoBrho, dAxoBrho, dpy, true);

    x[px_] += AxoBrho; x[py_] += dpy;

    // 4: a half drift in y
    AyoBrho = elem->Sol->BoBrho*x[x_]/2e0; dpx = elem->Sol->BoBrho*x[y_]/2e0;
//    get_Axy_EF3(elem->W, z, x, AyoBrho, dAyoBrho, dpx, false);

    x[px_] -= dpx; x[py_] -= AyoBrho;
    x[y_] += 0.5*hd*x[py_];
    x[ct_] += sqr(0.5)*hd*sqr(x[py_])/(1e0+x[delta_]);

    AyoBrho = elem->Sol->BoBrho*x[x_]/2e0; dpx = elem->Sol->BoBrho*x[y_]/2e0;
//    get_Axy_EF3(elem->W, z, x, AyoBrho, dAyoBrho, dpx, false);

    x[px_] += dpx; x[py_] += AyoBrho;

    // 5: half step in z
    z += 0.5*h;

    if (globval.radiation || globval.emittance) {
      dAxoBrho[X_] = 0e0;
      dAxoBrho[Y_] = -elem->Sol->BoBrho/2e0;
      dAxoBrho[Z_] = 0e0;
      dAyoBrho[X_] = elem->Sol->BoBrho/2e0;
      dAyoBrho[Y_] = 0e0;
      dAyoBrho[Z_] = 0e0;
//      get_Axy_EF3(elem->W, z, x, AyoBrho, dAyoBrho, dpx, false);
//      get_Axy_EF3(elem->W, z, x, AxoBrho, dAxoBrho, dpy, true);
      B[X_] = -dAyoBrho[Z_]; B[Y_] = dAxoBrho[Z_];
      B[Z_] = dAyoBrho[X_] - dAxoBrho[Y_];
      radiate(Cell, x, h, 0e0, B);
    }
  }
}


template<typename T>
void Solenoid_Pass(CellType &Cell, ss_vect<T> &ps)
{
  GtoL(ps, Cell.dS, Cell.dT, 0e0, 0e0, 0e0);
  sol_pass(Cell, ps);
  LtoG(ps, Cell.dS, Cell.dT, 0e0, 0e0, 0e0);
}


// template<typename T>
// void Map_Pass(CellType &Cell, ss_vect<T> &ps) { ps = Cell.Elem.Map->M*ps; }

void Map_Pass(const CellType &Cell, ss_vect<double> &ps)
{
  ps = (Cell.Elem.Map->M*ps).cst();
}

void Map_Pass(const CellType &Cell, ss_vect<tps> &ps)
{ ps = Cell.Elem.Map->M*ps; }


void getelem(long i, CellType *cellrec) { *cellrec = Cell[i]; }

void putelem(long i, const CellType *cellrec) { Cell[i] = *cellrec; }


int GetnKid(const int Fnum1) { return (ElemFam[Fnum1-1].nKid); }


long Elem_GetPos(const int Fnum1, const int Knum1)
{
  long int  loc;

  if (ElemFam[Fnum1-1].nKid != 0)
    loc = ElemFam[Fnum1-1].KidList[Knum1-1];
  else {
    loc = -1;
    printf("Elem_GetPos: there are no kids in family %d (%s)\n",
	   Fnum1, ElemFam[Fnum1-1].ElemF.PName);
    exit_(0);
  }

  return loc;
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


void SI_init(void)
{
  // SI units are used internally
  // apart from globval.energy [GeV]
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
  cl_rad = C_gamma*cube(globval.Energy)/(2e0*M_PI);

  // eletron rest mass [GeV]: slightly off???
//  m_e_ = 0.5110034e-03;
  // quantum fluctuations
  C_q = 3e0*C_u*h_bar*c0/(4e0*m_e);
  q_fluct = C_q*C_gamma/(M_PI*sqr(1e-9*m_e))*pow(globval.Energy, 5e0);
}


static void Mpole_Print(FILE *f, int Fnum1)
{
  elemtype  *elemp;
  MpoleType *M;

  elemp  = &ElemFam[Fnum1-1].ElemF; M = elemp->M;
  fprintf(f, "Element[%3d ] \n", Fnum1);
  fprintf(f, "   Name: %.*s,  Kind:   mpole,  L=% .8E\n",
          SymbolLength, elemp->PName, elemp->PL);
  fprintf(f, "   Method: %d, N=%4d\n", M->Pmethod, M->PN);
}


static void Drift_Print(FILE *f, int Fnum1)
{
  ElemFamType *elemfamp;
  elemtype    *elemp;

  elemfamp = &ElemFam[Fnum1-1]; elemp = &elemfamp->ElemF;
  fprintf(f, "Element[%3d ] \n", Fnum1);
  fprintf(f, "   Name: %.*s,  Kind:   drift,  L=% .8E\n",
          SymbolLength, elemp->PName, elemp->PL);
  fprintf(f, "   nKid:%3d\n\n", elemfamp->nKid);
}


static void Wiggler_Print(FILE *f, int Fnum1)
{
  elemtype *elemp;

  elemp = &ElemFam[Fnum1-1].ElemF;
  fprintf(f, "Element[%3d ] \n", Fnum1);
  fprintf(f, "   Name: %.*s,  Kind:   wiggler,  L=% .8E\n\n",
          NameLength, elemp->PName, elemp->PL);
}


static void Insertion_Print(FILE *f, int Fnum1)
{
  elemtype *elemp;

  elemp = &ElemFam[Fnum1-1].ElemF;
  fprintf(f, "Element[%3d ] \n", Fnum1);
  fprintf(f, "   Name: %.*s,  Kind:   wiggler,  L=% .8E\n\n",
          SymbolLength, elemp->PName, elemp->PL);
}


void Elem_Print(FILE *f, int Fnum1)
{
  int i;

  if (Fnum1 == 0) {
    // print all elements
    for (i = 1; i <= globval.Elem_nFam; i++)
      Elem_Print(f, i);
    return;
  }

  switch (ElemFam[Fnum1-1].ElemF.Pkind) {
  case drift:
    Drift_Print(f, Fnum1);
    break;

  case Mpole:
    Mpole_Print(f, Fnum1);
    break;
  case Wigl:
    Wiggler_Print(f, Fnum1);
    break;
  case FieldMap:
    break;
  case Insertion:
    Insertion_Print(f, Fnum1);
    break;
  case Cavity:
    break;
  case marker:
    break;
  case Spreader:
    break;
  case Recombiner:
    break;
  case Solenoid:
    break;
  case Map:
    break;
  case undef:
    break;
  }
}


double Mpole_GetPB(int Fnum1, int Knum1, int Order);


double Elem_GetKval(int Fnum1, int Knum1, int Order)
{
  double   Result = 0e0;
  elemtype *elemp;

  if (Fnum1 > 0) {
    elemp = &Cell[ElemFam[Fnum1-1].KidList[Knum1-1]].Elem;
    switch (elemp->Pkind) {
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
      if (elemp->M->Pthick == thick)
	Result = elemp->PL*Mpole_GetPB(Fnum1, Knum1, Order);
      else
	Result = Mpole_GetPB(Fnum1, Knum1, Order);
      break;
    case Wigl:
      Result =
	elemp->PL*sqrt(2e0
	*Cell[ElemFam[Fnum1-1].KidList[Knum1-1]].Elem.W->PBW[Order+HOMmax]);
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


void Drift_Alloc(elemtype *Elem)
{
  Elem->D = (DriftType *)malloc(sizeof(DriftType));
}


void Mpole_Alloc(elemtype *Elem)
{
  int       j;
  MpoleType *M;

  /* Memory allocation */
  Elem->M = (MpoleType *)malloc(sizeof(MpoleType));
  M = Elem->M; M->Pmethod = Meth_Fourth; M->PN = 0;
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
}


void Cav_Alloc(elemtype *Elem)
{
  CavityType *C;

  Elem->C = (CavityType *)malloc(sizeof(CavityType));
  C = Elem->C;
  C->V_RF = 0e0; C->f_RF = 0e0; C->phi_RF = 0e0; C->harm_num = 0;
  C->entry_focus = false; C->exit_focus = false;
}


void Wiggler_Alloc(elemtype *Elem)
{
  int         j;
  WigglerType *W;

  Elem->W = (WigglerType *)malloc(sizeof(WigglerType)); W = Elem->W;
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
}


void FieldMap_Alloc(elemtype *Elem)
{
  FieldMapType  *FM;

  Elem->FM = (FieldMapType *)malloc(sizeof(FieldMapType)); FM = Elem->FM;
  FM->n_step = 0; FM->n[X_] = 0; FM->n[Y_] = 0; FM->n[Z_] = 0; FM->scl = 1e0;
  FM->phi = 0e0; FM->Ld = 0e0; FM->L1 = 0e0; FM->cut = 0; FM->x0 = 0e0;
}


void Insertion_Alloc(elemtype *Elem)
{
  int           i = 0, j = 0;
  InsertionType *ID;

  Elem->ID = (InsertionType *)malloc(sizeof(InsertionType));
  ID = Elem->ID;

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
}


void Spreader_Alloc(elemtype *Elem)
{
  int k;

  Elem->Spr = (SpreaderType *)malloc(sizeof(SpreaderType));

  for (k = 0; k < Spreader_max; k++)
    Elem->Spr->Cell_ptrs[k] = NULL;
}


void Recombiner_Alloc(elemtype *Elem)
{
  Elem->Rec = (RecombinerType *)malloc(sizeof(RecombinerType));
}


void Solenoid_Alloc(elemtype *Elem)
{
  int          j;
  SolenoidType *Sol;

  Elem->Sol = (SolenoidType *)malloc(sizeof(SolenoidType));
  Sol = Elem->Sol; Sol->N = 0;
  for (j = 0; j <= 1; j++) {
    Sol->PdSsys[j] = 0e0; Sol->PdSrms[j] = 0e0; Sol->PdSrnd[j] = 0e0;
  }
  Sol->dTpar = 0e0; Sol->dTsys = 0e0; Sol->dTrnd = 0e0;
}


void Map_Alloc(elemtype *Elem)
{
  // Elem->Map = (MapType *)malloc(sizeof(MapType));
  // Use new; to allocate TPSA vector.
  Elem->Map = new MapType();
}


void Drift_Init(int Fnum1)
{
  int         i;
  ElemFamType *elemfamp;
  CellType    *cellp;
  elemtype    *elemp;

  elemfamp = &ElemFam[Fnum1-1];
  for (i = 1; i <= elemfamp->nKid; i++) {
    /* Get in Cell kid # i from Family Fnum1 */
    cellp = &Cell[elemfamp->KidList[i-1]]; elemp = &cellp->Elem;
    /* Dynamic memory allocation for element */
    Drift_Alloc(elemp);
    /* copy low level routine */
    memcpy(elemp->PName, elemfamp->ElemF.PName, sizeof(partsName));
    /* Set the drift length */
    elemp->PL = elemfamp->ElemF.PL;
    /* set the kind of element */
    elemp->Pkind = elemfamp->ElemF.Pkind;
    /* set pointer for the D dynamic space */
    *elemp->D = *elemfamp->ElemF.D;
    cellp->dT[0] = 1e0; /* cos = 1 */
    cellp->dT[1] = 0e0; /* sin = 0 */
    cellp->dS[0] = 0e0; /* no H displacement */
    cellp->dS[1] = 0e0; /* no V displacement */
  }
}


static int UpdatePorder(elemtype &Elem)
{
  int       i, order;
  MpoleType *M;

  M = Elem.M;
  if (M->Pirho != 0e0) /* non zero curvature => bend */
    order = 1;
  else /* mutipole */
    order = 0;
  for (i = -HOMmax; i <= HOMmax; i++)
    if (M->PB[i+HOMmax] != 0e0 && abs(i) > order) order = abs(i);

  return order;
}


ss_vect<tps> get_edge_lin_map(const double h, const double phi,
			      const double gap, const double delta)
{
  ss_vect<tps> Id, M;

  Id.identity();

  M.identity();
  M[px_] += h*tan(dtor(phi))*Id[x_];
  M[py_] -= h*tan(dtor(phi)-get_psi(h, phi, gap))*Id[y_];

  return M;
}


ss_vect<tps> get_sbend_lin_map(const double L, const double h, const double b2,
			       const double delta)
{
  double       K_x, K_y, psi_x, psi_y;
  ss_vect<tps> Id, M;

  Id.identity();

  K_x = b2 + sqr(h);
  K_y = fabs(b2);
  psi_x = sqrt(fabs(K_x)/(1e0+delta))*L;
  psi_y = sqrt(K_y/(1e0+delta))*L;

  M.identity();
  if (K_x > 0e0) {
    M[x_] =
      cos(psi_x)*Id[x_] + sin(psi_x)/sqrt(K_x*(1e0+delta))*Id[px_]
      + (1e0-cos(psi_x))*h/K_x*Id[delta_];
    M[px_] =
      -sqrt(K_x*(1e0+delta))*sin(psi_x)*Id[x_] + cos(psi_x)*Id[px_]
      + sin(psi_x)*sqrt(1e0+delta)*h/sqrt(K_x)*Id[delta_];

    if (psi_y != 0e0) {
      M[y_]  = cosh(psi_y)*Id[y_] + sinh(psi_y)/sqrt(K_y*(1e0+delta))*Id[py_];
      M[py_] = sqrt(K_y*(1e0+delta))*sinh(psi_y)*Id[y_] + cosh(psi_y)*Id[py_];
    } else
      M[y_]  += L/(1e0+delta)*Id[py_];
 
    M[ct_] +=
      sin(psi_x)*sqrt(1e0+delta)*h/sqrt(K_x)*Id[x_]
      + (1e0-cos(psi_x))*h/K_x*Id[px_]
      + (psi_x-sin(psi_x))*sqrt(1e0+delta)*sqr(h)/pow(K_x, 3e0/2e0)*Id[delta_];
  } else if (K_x < 0e0) {
    K_x = fabs(K_x);
    M[x_] =
      cosh(psi_x)*Id[x_] + sinh(psi_x)/sqrt(K_x*(1e0+delta))*Id[px_]
      -(1e0-cosh(psi_x))*h/K_x*Id[delta_];
    M[px_] =
      sqrt(K_x*(1e0+delta))*sinh(psi_x)*Id[x_] + cosh(psi_x)*Id[px_]
      + sinh(psi_x)*sqrt(1e0+delta)*h/sqrt(K_x)*Id[delta_];

    if (psi_y != 0e0) {
      M[y_]  = cos(psi_y)*Id[y_] + sin(psi_y)/sqrt(K_y*(1e0+delta))*Id[py_];
      M[py_] = -sqrt(K_y*(1e0+delta))*sin(psi_y)*Id[y_] + cos(psi_y)*Id[py_];
   } else
      M[y_]  += L/(1e0+delta)*Id[py_];

    M[ct_] +=
      sinh(psi_x)*sqrt(1e0+delta)*h/sqrt(K_x)*Id[x_]
      - (1e0-cosh(psi_x))*h/K_x*Id[px_]
      - (psi_x-sinh(psi_x))*sqrt(1e0+delta)*sqr(h)/pow(K_x, 3e0/2e0)*Id[delta_];
  } else {
    // K_x = 0.
    M[x_] += L/(1e0+delta)*Id[px_];
    M[y_] += L/(1e0+delta)*Id[py_];
  }

  return M;
}


ss_vect<tps> get_thin_kick_lin_map(const double b2L, const double delta)
{
  ss_vect<tps> Id, M;

  Id.identity();

  M.identity();
  M[px_] -= b2L*(1e0+delta)*Id[x_];
  M[py_] += b2L*(1e0+delta)*Id[y_];

  return M;
}


ss_vect<tps> get_lin_map(elemtype &Elem, const double delta)
{
  ss_vect<tps> M, M1, M2;

  if (Elem.PL != 0e0) {
    M =
      get_sbend_lin_map(Elem.PL, Elem.M->Pirho, Elem.M->PB[Quad+HOMmax], delta);
    M1 = get_edge_lin_map(Elem.M->Pirho, Elem.M->PTx1, Elem.M->Pgap, delta);
    M2 = get_edge_lin_map(Elem.M->Pirho, Elem.M->PTx2, Elem.M->Pgap, delta);
    M = M2*M*M1;
  } else
    M = get_thin_kick_lin_map(Elem.M->PB[Quad+HOMmax], delta);

  return M;
}


void get_lin_maps(const double delta)
{
  long int k;

  for (k = 0; k <= globval.Cell_nLoc; k++)
    if (Cell[k].Elem.Pkind == Mpole)
      Cell[k].Elem.M->M_lin = get_lin_map(Cell[k].Elem, delta);
}


void Mpole_Init(int Fnum1)
{
  static bool first = true;
  int          i;
  double       phi;
  ElemFamType  *elemfamp;
  CellType     *cellp;
  elemtype     *elemp;

  /* Pointer on element */
  elemfamp = &ElemFam[Fnum1-1];
  memcpy(elemfamp->ElemF.M->PB, elemfamp->ElemF.M->PBpar, sizeof(mpolArray));
  /* Update the right multipole order */
  elemfamp->ElemF.M->Porder = UpdatePorder(elemfamp->ElemF);

  for (i = 1; i <= elemfamp->nKid; i++) {
    cellp = &Cell[elemfamp->KidList[i-1]]; elemp = &cellp->Elem;
    /* Memory allocation and set everything to zero */
    Mpole_Alloc(elemp);
    memcpy(elemp->PName, elemfamp->ElemF.PName, sizeof(partsName));
    /* set length */
    elemp->PL = elemfamp->ElemF.PL;
    /* set element kind (Mpole) */
    elemp->Pkind = elemfamp->ElemF.Pkind;
    *elemp->M = *elemfamp->ElemF.M;

    if (reverse_elem && (elemp->Reverse == true)) {
      // Swap entrance and exit angles.
      if (first) {
	printf("\nSwapping entrance and exit angles for %8s %2d\n",
	       elemp->PName, i);
	printf("...\n");
	first = false;
      }
      phi = elemp->M->PTx1;
      elemp->M->PTx1 = elemp->M->PTx2; elemp->M->PTx2 = phi; 
    }

    /* set entrance and exit angles */
    cellp->dT[0] = cos(dtor(elemp->M->PdTpar));
    cellp->dT[1] = sin(dtor(elemp->M->PdTpar));

    /* set displacement to zero */
    cellp->dS[0] = 0e0; cellp->dS[1] = 0e0;

    if (elemp->PL != 0e0 || elemp->M->Pirho != 0e0) {
      /* Thick element or radius non zero element */
      elemp->M->Pthick = pthicktype(thick);
      /* sin(L*irho/2) =sin(theta/2) half the angle */
      elemp->M->Pc0 = sin(elemp->PL*elemp->M->Pirho/2e0);
      /* cos roll: sin(theta/2)*cos(dT) */
      elemp->M->Pc1 = cellp->dT[0]*elemp->M->Pc0;
      /* sin roll: sin(theta/2)*cos(dT) */
      elemp->M->Ps1 = cellp->dT[1]*elemp->M->Pc0;
    } else /* element as thin lens */
      elemp->M->Pthick = pthicktype(thin);

    // Allocate TPSA vector.
    if (globval.mat_meth)
      elemp->M->M_lin = get_lin_map(*elemp, 0e0);
  }
}


#define order           2
void Wiggler_Init(int Fnum1)
{
  int         i;
  ElemFamType *elemfamp;
  CellType    *cellp;
  elemtype    *elemp;

  elemfamp = &ElemFam[Fnum1-1];
  /* ElemF.M^.PB := ElemF.M^.PBpar; */
  elemfamp->ElemF.W->Porder = order;
  for (i = 1; i <= elemfamp->nKid; i++) {
    cellp = &Cell[elemfamp->KidList[i-1]]; elemp = &cellp->Elem;
    Wiggler_Alloc(elemp);
    memcpy(elemp->PName, elemfamp->ElemF.PName, sizeof(partsName));
    elemp->PL = elemfamp->ElemF.PL;
    elemp->Pkind = elemfamp->ElemF.Pkind;
    *elemp->W = *elemfamp->ElemF.W;

    // 2/21/12 JB & JC
//     cellp->dT[0] = cos(dtor(elemp->M->PdTpar));
//     cellp->dT[1] = sin(dtor(elemp->M->PdTpar));
    cellp->dT[0] = cos(dtor(elemp->W->PdTpar));
    cellp->dT[1] = sin(dtor(elemp->W->PdTpar));

    cellp->dS[0] = 0e0; cellp->dS[1] = 0e0;
 }
}
#undef order


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


void get_B_DIAMOND(const char *filename, FieldMapType *FM)
{
  char          line[max_str];
  int           i, j, n, ny;
  double        x0, y0, z0;
  std::ifstream inf;
  std::ofstream outf;

  const int     skip = 8;
  const double  Brho = globval.Energy*1e9/c0;

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


void get_B_NSLS_II(const char *filename, FieldMapType *FM)
{
  char          line[max_str];
  int           i, j, n;
  double        x_min[3], x_max[3];
  std::ifstream inf;

  const double  Brho = globval.Energy*1e9/c0;

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


void get_B_Oleg1(const char *filename, FieldMapType *FM)
{
  char          line[max_str];
  int           i, j, n;
  double        x_min[3], x_max[3];
  std::ifstream inf;

  const double Brho = globval.Energy*1e9/c0;

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


void get_B_Oleg2(const char *filename, FieldMapType *FM)
{
  char          line[max_str];
  int           i, j, n;
  double        x_min[3];
  std::ifstream inf;

  const double Brho = globval.Energy*1e9/c0;

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


void get_B_SRW(const char *filename, FieldMapType *FM)
{
  char          line[max_str];
  int           i, j, n;
  double        x_min[3];
  std::ifstream inf;

  const double Brho = globval.Energy*1e9/c0;

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


void get_B(const char *filename, FieldMapType *FM)
{
  // Do not scale fieldmaps only Hamiltonians, i.e., the kick.  Note that RADIA
  // (2nd order) kick maps are quadratic in the field, and 1st order linear.

  switch (FieldMap_filetype) {
  case 1:
    get_B_DIAMOND(filename, FM);
    break;
  case 2:
    get_B_NSLS_II(filename, FM);
    break;
  case 3:
    get_B_Oleg1(filename, FM);
    break;
  case 4:
    get_B_Oleg2(filename, FM);
    break;
  case 5:
    get_B_SRW(filename, FM);
    break;
  default:
    printf("\nget_B: unknown FieldMap type %d", FieldMap_filetype);
    exit(1);
    break;
  }
}


void FieldMap_Init(int Fnum1)
{
  int         i;
  ElemFamType *elemfamp;
  CellType    *cellp;
  elemtype    *elemp;

  elemfamp = &ElemFam[Fnum1-1];
  for (i = 1; i <= elemfamp->nKid; i++) {
    cellp = &Cell[elemfamp->KidList[i-1]]; elemp = &cellp->Elem;
    FieldMap_Alloc(elemp);
    memcpy(elemp->PName, elemfamp->ElemF.PName, sizeof(partsName));
    elemp->PL = elemfamp->ElemF.PL;
    elemp->Pkind = elemfamp->ElemF.Pkind;
    *elemp->FM = *elemfamp->ElemF.FM;

    cellp->dT[0] = 1e0; cellp->dT[1] = 0e0;
    cellp->dS[X_] = 0e0; cellp->dS[Y_] = 0e0;
  }
}


void Cav_Init(int Fnum1)
{
  int         i;
  ElemFamType *elemfamp;
  CellType    *cellp;

  elemfamp = &ElemFam[Fnum1-1];
  for (i = 0; i < elemfamp->nKid; i++) {
    cellp = &Cell[elemfamp->KidList[i]];
    cellp->Elem = elemfamp->ElemF;
  }
}


void Marker_Init(int Fnum1)
{
  int         i;
  ElemFamType *elemfamp;
  CellType    *cellp;

  elemfamp = &ElemFam[Fnum1-1];
  for (i = 0; i < elemfamp->nKid; i++) {
    cellp = &Cell[elemfamp->KidList[i]];
    cellp->Elem  = elemfamp->ElemF;
    cellp->dT[0] = 1e0; cellp->dT[1] = 0e0;
    cellp->dS[0] = 0e0; cellp->dS[1] = 0e0;
  }
}


void Insertion_Init(int Fnum1)
{
  int         i;
  ElemFamType *elemfamp;
  CellType    *cellp;
  elemtype    *elemp;

  elemfamp = &ElemFam[Fnum1-1];
//  elemfamp->ElemF.ID->Porder = order;
//  x = elemfamp->ElemF.ID->PBW[Quad + HOMmax];
  for (i = 1; i <= elemfamp->nKid; i++) {
    cellp = &Cell[elemfamp->KidList[i-1]]; elemp = &cellp->Elem;
    Insertion_Alloc(elemp);
    memcpy(elemp->PName, elemfamp->ElemF.PName, sizeof(partsName));
    elemp->PL = elemfamp->ElemF.PL;
    elemp->Pkind = elemfamp->ElemF.Pkind;
    *elemp->ID = *elemfamp->ElemF.ID;

    cellp->dT[0] = cos(dtor(elemp->ID->PdTpar));
    cellp->dT[1] = sin(dtor(elemp->ID->PdTpar));
    cellp->dS[0] = 0e0; cellp->dS[1] = 0e0;
  }
}


void Spreader_Init(int Fnum1)
{
  int         i;
  ElemFamType *elemfamp;
  CellType    *cellp;
  elemtype    *elemp;

  elemfamp = &ElemFam[Fnum1-1];
  for (i = 1; i <= elemfamp->nKid; i++) {
    /* Get in Cell kid # i from Family Fnum1 */
    cellp = &Cell[elemfamp->KidList[i-1]]; elemp = &cellp->Elem;
    /* Dynamic memory allocation for element */
    Spreader_Alloc(elemp);
    /* copy low level routine */
    memcpy(elemp->PName, elemfamp->ElemF.PName, sizeof(partsName));
    /* set the kind of element */
    elemp->Pkind = elemfamp->ElemF.Pkind;
    /* set pointer for the dynamic space */
    *elemp->Spr = *elemfamp->ElemF.Spr;
    cellp->dT[0] = 1e0; /* cos = 1 */
    cellp->dT[1] = 0e0; /* sin = 0 */
    cellp->dS[0] = 0e0; /* no H displacement */
    cellp->dS[1] = 0e0; /* no V displacement */
  }
}


void Recombiner_Init(int Fnum1)
{
  int         i;
  ElemFamType *elemfamp;
  CellType    *cellp;
  elemtype    *elemp;

  elemfamp = &ElemFam[Fnum1-1];
  for (i = 1; i <= elemfamp->nKid; i++) {
    /* Get in Cell kid # i from Family Fnum1 */
    cellp = &Cell[elemfamp->KidList[i-1]]; elemp = &cellp->Elem;
    /* Dynamic memory allocation for element */
    Spreader_Alloc(elemp);
    /* copy low level routine */
    memcpy(elemp->PName, elemfamp->ElemF.PName, sizeof(partsName));
    /* set the kind of element */
    elemp->Pkind = elemfamp->ElemF.Pkind;
    /* set pointer for the dynamic space */
    *elemp->Rec = *elemfamp->ElemF.Rec;
    cellp->dT[0] = 1e0; /* cos = 1 */
    cellp->dT[1] = 0e0; /* sin = 0 */
    cellp->dS[0] = 0e0; /* no H displacement */
    cellp->dS[1] = 0e0; /* no V displacement */
  }
}


void Solenoid_Init(int Fnum1)
{
  int         i;
  ElemFamType *elemfamp;
  CellType    *cellp;
  elemtype    *elemp;

  /* Pointer on element */
  elemfamp = &ElemFam[Fnum1-1];
  for (i = 1; i <= elemfamp->nKid; i++) {
    cellp = &Cell[elemfamp->KidList[i-1]]; elemp = &cellp->Elem;
    /* Memory allocation and set everything to zero */
    Solenoid_Alloc(elemp);
    memcpy(elemp->PName, elemfamp->ElemF.PName, sizeof(partsName));
    /* set length */
    elemp->PL = elemfamp->ElemF.PL;
    /* set element kind */
    elemp->Pkind = elemfamp->ElemF.Pkind;
    *elemp->Sol = *elemfamp->ElemF.Sol;

    /* set entrance and exit angles */
    cellp->dT[0] = 1e0; cellp->dT[1] = 0e0;
    /* set displacement to zero */
    cellp->dS[0] = 0e0; cellp->dS[1] = 0e0;
  }
}


void Map_Init(int Fnum1)
{
  int         i;
  ElemFamType *elemfamp;
  CellType    *cellp;
  elemtype    *elemp;

  /* Pointer on element */
  elemfamp = &ElemFam[Fnum1-1];
  for (i = 1; i <= elemfamp->nKid; i++) {
    cellp = &Cell[elemfamp->KidList[i-1]]; elemp = &cellp->Elem;
    /* Memory allocation and set everything to zero */
    Map_Alloc(elemp);
    memcpy(elemp->PName, elemfamp->ElemF.PName, sizeof(partsName));
    /* set length */
    elemp->PL = 0e0;
    /* set element kind */
    elemp->Pkind = elemfamp->ElemF.Pkind;
    *elemp->Map = *elemfamp->ElemF.Map;

    elemp->Map->M.identity();

    /* set entrance and exit angles */
    cellp->dT[0] = 1e0; cellp->dT[1] = 0e0;
    /* set displacement to zero */
    cellp->dS[0] = 0e0; cellp->dS[1] = 0e0;
  }
}


void Mpole_SetPB(int Fnum1, int Knum1, int Order)
{
  /* Compute full multipole composent as sum of design, systematic
     and random part
     Compute transport matrix if quadrupole (Order=2)
     Set multipole order to Order if multipole (Order >2)                  */

  CellType  *cellp; /* pointer on the Cell */
  elemtype  *elemp; /* pointer on the Elemetype */
  MpoleType *M;     /* Pointer on the Multipole */

  cellp  = &Cell[ElemFam[Fnum1-1].KidList[Knum1-1]]; elemp = &cellp->Elem;
  M = elemp->M;
  M->PB[Order+HOMmax] =
    M->PBpar[Order+HOMmax] + M->PBsys[Order+HOMmax] +
    M->PBrms[Order+HOMmax]*M->PBrnd[Order+HOMmax];
  if (abs(Order) > M->Porder && M->PB[Order+HOMmax] != 0e0)
    M->Porder = abs(Order);
}


double Mpole_GetPB(int Fnum1, int Knum1, int Order)
{
  /*  Return multipole strength (of order Order) for Knum1 element of
      family Fnum1
       Order =  2 for normal quadrupole
             = -2 for skew quadrupole                                        */

  MpoleType *M; /* Pointer on the multipole */

  M = Cell[ElemFam[Fnum1-1].KidList[Knum1-1]].Elem.M;
  return (M->PB[Order+HOMmax]);
}


void Mpole_DefPBpar(int Fnum1, int Knum1, int Order, double PBpar)
{
  elemtype  *elemp;
  MpoleType *M;

  elemp = &Cell[ElemFam[Fnum1-1].KidList[Knum1-1]].Elem; M = elemp->M;

  M->PBpar[Order+HOMmax]=PBpar;
}


void Mpole_DefPBsys(int Fnum1, int Knum1, int Order, double PBsys)
{
  /*Fnum1, Knum1, Order : integer*/
  elemtype  *elemp;
  MpoleType *M;

  elemp = &Cell[ElemFam[Fnum1-1].KidList[Knum1-1]].Elem; M = elemp->M;

  M->PBsys[Order+HOMmax]=PBsys;
}


void Mpole_SetdS(int Fnum1, int Knum1)
{
  int       j;
  CellType  *cellp;
  elemtype  *elemp;
  MpoleType *M;

  cellp = &Cell[ElemFam[Fnum1-1].KidList[Knum1-1]]; elemp = &cellp->Elem;
  M = elemp->M;
  for (j = 0; j <= 1; j++)
    cellp->dS[j] = M->PdSsys[j] + M->PdSrms[j]*M->PdSrnd[j];
}

void Mpole_SetdT(int Fnum1, int Knum1)
{
  CellType  *cellp;
  elemtype  *elemp;
  MpoleType *M;

  cellp = &Cell[ElemFam[Fnum1-1].KidList[Knum1-1]]; elemp = &cellp->Elem;
  M = elemp->M;
  cellp->dT[0] =
    cos(dtor(M->PdTpar + M->PdTsys + M->PdTrms*M->PdTrnd));
  cellp->dT[1] = sin(
      dtor(M->PdTpar + M->PdTsys + M->PdTrms*M->PdTrnd));
  /* Calculate simplified p_rots */
  M->Pc0 = sin(elemp->PL*M->Pirho/2e0);
  M->Pc1 = cos(dtor(M->PdTpar))*M->Pc0;
  M->Ps1 = sin(dtor(M->PdTpar))*M->Pc0;
}


double Mpole_GetdT(int Fnum1, int Knum1)
{
  elemtype  *elemp;
  MpoleType *M;

  elemp = &Cell[ElemFam[Fnum1-1].KidList[Knum1-1]].Elem; M = elemp->M;

  return(M->PdTpar + M->PdTsys + M->PdTrms*M->PdTrnd);
}


void Mpole_DefdTpar(int Fnum1, int Knum1, double PdTpar)
{
  elemtype  *elemp;
  MpoleType *M;

  elemp = &Cell[ElemFam[Fnum1-1].KidList[Knum1-1]].Elem; M = elemp->M;

  M->PdTpar = PdTpar;
}


void Mpole_DefdTsys(int Fnum1, int Knum1, double PdTsys)
{
  elemtype  *elemp;
  MpoleType *M;

  elemp = &Cell[ElemFam[Fnum1-1].KidList[Knum1-1]].Elem; M = elemp->M;

  M->PdTsys=PdTsys;
}


void Wiggler_SetPB(int Fnum1, int Knum1, int Order)
{
  CellType    *cellp;
  elemtype    *elemp;
  WigglerType *W;

  cellp = &Cell[ElemFam[Fnum1-1].KidList[Knum1-1]]; elemp = &cellp->Elem;
  W = elemp->W;
  if (abs(Order) > W->Porder)
    W->Porder = abs(Order);
}


void Wiggler_SetdS(int Fnum1, int Knum1)
{
  int         j;
  CellType    *cellp;
  elemtype    *elemp;
  WigglerType *W;

  cellp = &Cell[ElemFam[Fnum1-1].KidList[Knum1-1]]; elemp = &cellp->Elem;
  W = elemp->W;
  for (j = 0; j <= 1; j++)
    cellp->dS[j] = W->PdSsys[j] + W->PdSrms[j]*W->PdSrnd[j];
}

void Wiggler_SetdT(int Fnum1, int Knum1)
{
  CellType    *cellp;
  elemtype    *elemp;
  WigglerType *W;

  cellp = &Cell[ElemFam[Fnum1-1].KidList[Knum1-1]]; elemp = &cellp->Elem;
  W = elemp->W;
  cellp->dT[0] = cos(dtor(W->PdTpar+W->PdTsys+W->PdTrms*W->PdTrnd));
  cellp->dT[1] = sin(dtor(W->PdTpar+W->PdTsys+W->PdTrms*W->PdTrnd));
}
