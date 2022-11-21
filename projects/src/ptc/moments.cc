#define NO 2

#include "tracy_lib.h"


// F. Klein's Erlangen Program.


int no_tps   = NO,
    ndpt_tps = 5;

//------------------------------------------------------------------------------
// Template class for complex TPSA.

const int ps_dim = 2*nd_tps;

class ctps
{
private:
  tps re, im;

public:
  ctps(void) { this->re = 0e0; this->im = 0e0; }
  template<typename T1, typename T2>
  ctps(const T1 &real, const T2 &imag) { this->re = real; this->im = imag; }
  template<typename T>
  ctps(const T real) { this->re = real; this->im = 0e0; }

  template<typename CharT, class Traits>
  friend std::basic_ostream<CharT, Traits>&
  operator<<(std::basic_ostream<CharT, Traits> &, const ctps &);

  void print(void)
  {
      cout << scientific << setprecision(3) << *this;
  }

  void print(const string &str)
  {
    cout << str;
    this->print();
  }

  tps& real(void) { return this->re; }
  tps& imag(void) { return this->im; }
  const tps& real(void) const { return this->re; }
  const tps& imag(void) const { return this->im; }

  ctps operator+=(const double);
  ctps operator-=(const double);
  ctps operator*=(const double);
  ctps operator/=(const double);

  ctps operator+=(const ctps &);
  ctps operator-=(const ctps &);
  ctps operator*=(const ctps &);
  ctps operator/=(const ctps &);

  friend ctps conjugate(const ctps &);
  friend ctps abs(const ctps &);

  friend ctps operator+(const double a, const ctps &b) { return ctps(a) += b; }
  friend ctps operator-(const double a, const ctps &b) { return ctps(a) -= b; }
  friend ctps operator*(const double a, const ctps &b) { return ctps(a) *= b; }
  friend ctps operator/(const double a, const ctps &b) { return ctps(a) /= b; }

  friend ctps operator+(const ctps &a, const double b) { return ctps(a) += b; }
  friend ctps operator-(const ctps &a, const double b) { return ctps(a) -= b; }
  friend ctps operator*(const ctps &a, const double b) { return ctps(a) *= b; }
  friend ctps operator/(const ctps &a, const double b) { return ctps(a) /= b; }

  friend ctps operator+(const ctps &a, const ctps &b) { return ctps(a) += b; }
  friend ctps operator-(const ctps &a, const ctps &b) { return ctps(a) -= b; }
  friend ctps operator*(const ctps &a, const ctps &b) { return ctps(a) *= b; }
  friend ctps operator/(const ctps &a, const ctps &b) { return ctps(a) /= b; }

  friend ctps operator+(const ctps &x) { return ctps(x); }
  friend ctps operator-(const ctps &x) { return ctps(x) *= ctps(-1e0, 0e0); }
};

template<typename CharT, class Traits>
std::basic_ostream<CharT, Traits>&
operator<<(std::basic_ostream<CharT, Traits> &os, const ctps &a)
{
  std::basic_ostringstream<CharT, Traits> s;

  s.flags(os.flags());
  s.imbue(os.getloc());
  s << std::setprecision(os.precision()) << std::setw(os.width())
    << "\nReal:\n" << a.real() << "\nImag:\n" << a.imag();
  //   s << endl;
  return os << s.str();
}

ctps ctps::operator+=(const double a)
{
  this->real() += a;
  return *this;
}

ctps ctps::operator-=(const double a)
{
  this->real() -= a;
  return *this;
}

ctps ctps::operator*=(const double a)
{
  this->real() *= a;
  this->imag() *= a;
  return *this;
}

ctps ctps::operator/=(const double a)
{
  this->real() /= a;
  this->imag() /= a;
  return *this;
}

ctps ctps::operator+=(const ctps &a)
{
  this->real() += a.real();
  this->imag() += a.imag();
  return *this;
}

ctps ctps::operator-=(const ctps &a)
{
  this->real() -= a.real();
  this->imag() -= a.imag();
  return *this;
}

ctps ctps::operator*=(const ctps &a)
{
  const ctps b = *this;
  this->real() = b.real()*a.real()-b.imag()*a.imag();
  this->imag() = b.real()*a.imag()+b.imag()*a.real();
  return *this;
}

ctps ctps::operator/=(const ctps &a)
{
  return this->operator*=(conjugate(a)/(a*conjugate(a)));
}

inline ctps conjugate(const ctps &a)
{
  return ctps(a.real(), -a.imag());
}

inline ctps abs(const ctps &a)
{
  return a*conjugate(a);
}


class css_vect
{
private:
  ctps css[6];
public:
  css_vect(void)
  {
    for (int k = 0; k < ps_dim; k++) 
      css[k] = ctps(0e0, 0e0);
  }
  css_vect(const ss_vect<tps> &real, const ss_vect<tps> &imag)
  {
    for (int k = 0; k < ps_dim; k++) 
      css[k] = ctps(real[k], imag[k]);
  }
  css_vect(const ss_vect<tps> &real)
  {
    for (int k = 0; k < ps_dim; k++) 
      css[k] = ctps(real[k], 0e0);
  }

  template<typename CharT, class Traits>
  friend std::basic_ostream<CharT, Traits>&
  operator<<(std::basic_ostream<CharT, Traits> &, const css_vect &);

  void print(void)
  {
      cout << scientific << setprecision(3) << *this;
  }

  void print(const string &str)
  {
    cout << str;
    this->print();
  }

  ss_vect<tps> get_real(void) const
  {
    ss_vect<tps> a;
    for (int k = 0; k < ps_dim; k++)
      a[k] = css[k].real();
    return a;
  }
  ss_vect<tps> get_imag(void) const
  {
    ss_vect<tps> a;
    for (int k = 0; k < ps_dim; k++)
      a[k] = css[k].imag();
    return a;
  }

  css_vect& zero(void)
  {
    for (int k = 0; k < ps_dim; k++)
      css[k] = ctps(0e0, 0e0);
    return *this;
  }

  css_vect& identity(void)
  {
    for (int k = 0; k < ps_dim; k++)
      css[k] = ctps(tps(0e0, k+1), 0e0);
    return *this;
  }

  ctps& operator[](const int i) { return css[i]; }
  const ctps& operator[](const int i) const { return css[i]; }

  friend ctps operator*(const ctps &A, const css_vect &B);
  friend css_vect operator*(const css_vect &, const css_vect &);
};

template<typename CharT, class Traits>
std::basic_ostream<CharT, Traits>&
operator<<(std::basic_ostream<CharT, Traits> &os, const css_vect &a)
{
  std::basic_ostringstream<CharT, Traits> s;

  const ss_vect<tps>
    re = a.get_real(),
    im = a.get_imag();
    
  s.flags(os.flags());
  s.imbue(os.getloc());
  s << std::setprecision(os.precision()) << std::setw(os.width())
    << "\nReal:\n" << re << "\nImag:\n" << im;
  //   s << endl;
  return os << s.str();
}

ctps operator*(const ctps &A, const css_vect &B)
{
  const tps
    A_re = A.real(),
    A_im = A.imag();
  const ss_vect<tps>
    B_re = B.get_real(),
    B_im = B.get_imag();

  return ctps(A_re*B_re-A_im*B_im, A_re*B_im+A_im*B_re);
}

css_vect operator*(const css_vect &A, const css_vect &B)
{
  const ss_vect<tps>
    A_re = A.get_real(),
    A_im = A.get_imag(),
    B_re = B.get_real(),
    B_im = B.get_imag();

  return css_vect(A_re*B_re-A_im*B_im, A_re*B_im+A_im*B_re);
}


void prt_cmplx_lin_map(const int n_DOF, const string &str, const css_vect &map)
{
  int          i, j;
  ss_vect<tps> re, im;

  for (i = 0; i < 2*n_DOF; i++) {
    re[i] = map[i].real();
    im[i] = map[i].imag();
  }

  cout << str;
  for (i = 1; i <= 2*n_DOF; i++) {
    for (j = 1; j <= 2*n_DOF; j++)
	cout << scientific << setprecision(3)
	     << setw(11) << getmat(re, i, j)
	     << ((getmat(im, i, j) > 0e0)? " + i":" - i")
	     << setw(0) << fabs(getmat(im, i, j));
    cout << "\n";
  }
}

//------------------------------------------------------------------------------

ss_vect<tps> Id(void)
{
  ss_vect<tps> Id;
  Id.identity();
  return Id;
}


tps get_h_k(const tps &h, const int k1, const int k2)
{
  // Get monomials of order [k1..k2].
  long int no;
  tps      h_k;

  no = getno_();
  danot_(k1-1);
  h_k = -h;
  danot_(k2);
  h_k += h;
  danot_(no);
  return h_k;
}


ss_vect<tps> get_map_k(const ss_vect<tps> &x, const int k1, const int k2)
{
  int          k;
  ss_vect<tps> map_k;

  for (k = 0; k < nv_tps; k++)
    map_k[k] = get_h_k(x[k], k1, k2);
  return map_k;
}


void GoFix(const ss_vect<tps> &map, MNF_struct &MNF, const int no)
{
}


ss_vect<tps> get_map_Fl(ss_vect<tps> &map)
{
  Matrix       R_mat;
  ss_vect<tps> Id, R;

  const int n_dim = 4;

  Id.identity();
  // Find fix point.
  GoFix(map, MNF.A0, MNF.A0_inv, no_tps);
  printf("\nA0:\n");
  prt_lin_map(3, MNF.A0);
  // Translate to fix point.
  map = MNF.A0_inv*map*MNF.A0;
  getlinmat(n_dim, map, globval.OneTurnMat);
  GDiag(n_dim, 1e0, globval.Ascr, globval.Ascrinv, R_mat,
        globval.OneTurnMat, globval.Omega, globval.Alphac);
  MNF.A1 = putlinmat(n_dim, globval.Ascr);
  MNF.A1[ct_]    = Id[ct_];
  MNF.A1[delta_] = Id[delta_];
  R = putlinmat(n_dim, R_mat);
  R[ct_]    = Id[ct_];
  R[delta_] = Id[delta_];
  // Transform to Floquet space.
  printf("\nA1:\n");
  prt_lin_map(3, MNF.A1);
  map = Inv(MNF.A1)*map*MNF.A1;
  return R;
}


ss_vect<tps> compute_Dragt_Finn_Map
(const tps &h, const ss_vect<tps> &map, const int k1, const int k2,
 const int reverse)
{
  // Dragt-Finn map:
  // not reverse:
  //   exp(:h_3:) exp(:h_4:) ... exp(:h_no:)
  // reverse:
  //   exp(:h_no:) exp(:h_no-1:) ... exp(:h_3:)
  int          k;
  tps          h_k;
  ss_vect<tps> map1;

  map1.identity();
  for (k = k2; k >= k1; k--) {
    h_k = get_h_k(h, k, k);
    if (!reverse)
      map1 = map1*LieExp(h_k, map);
    else
      map1 = LieExp(h_k, map)*map1;
  }
  return map1;
}


tps LieFact_JB(const ss_vect<tps> &map)
{
  // Dragt-Finn factorization:
  //   M = M_1 exp(:h_3:) ... exp(:h_4:) exp(:h_no:).
  int          k;
  tps          h, h_k;
  ss_vect<tps> Id, A0_inv, M_1, Fn, map1;

  Id.identity();
  map1 = map;
  M_1 = get_map_Fl(map1);
  map1 = map1*Inv(M_1);
  h = 0e0;
  for (k = 3; k <= no_tps; k++) {
    Fn = get_map_k(map1, k-1, k-1);
    h_k = Intd(Fn, -1e0);
    h += h_k;
    map1 = map1*compute_Dragt_Finn_Map(-h_k, Id, k, k, true);
  }
//  cout << h;
  return h;
}


tps get_Ker(const tps &h)
{
  int          i, j, k;
  long int     jj[ss_dim];
  tps          h_Ke;
  ss_vect<tps> Id;

  for (k = 0; k < nv_tps; k++)
    jj[k] = 0;

  Id.identity();
  h_Ke = 0e0;
  for (i = 0; i <= no_tps; i++) {
    jj[x_] = i; jj[px_] = i;
    for (j = 0; j <= no_tps; j++) {
      jj[y_] = j;
      jj[py_] = j;
      for (k = 0; k <= no_tps; k++) {
	jj[delta_] = k;
	if ((2*i+2*j+k <= no_tps) && ((i != 0) || (j != 0) || (k != 0))) {
	  h_Ke +=
	    h[jj]*pow(Id[x_], i)*pow(Id[px_], i)
	    *pow(Id[y_], j)*pow(Id[py_], j)*pow(Id[delta_], k);
	}
      }
    }
  }

  return h_Ke;
}


tps get_g(const tps nu_x, const tps nu_y, const tps &h)
{
  // Compute g = (1-R)^-1 * h 

  int          i, j, k, l, m;
  long int     jj1[ss_dim], jj2[ss_dim];
  double       re, im;
  tps          h_re, h_im, g_re, g_im, mn1, mn2, cotan;
  ss_vect<tps> Id;

  CtoR(h, h_re, h_im);

  for (k = 0; k < nv_tps; k++) {
    jj1[k] = 0; jj2[k] = 0;
  }

  Id.identity();
  g_re = 0e0;
  g_im = 0e0;
  for (i = 0; i <= no_tps; i++) {
    jj1[x_] = i; jj2[px_] = i;
    for (j = 0; j <= no_tps; j++) {
      jj1[px_] = j; jj2[x_] = j;
      for (k = 0; k <= no_tps; k++) {
	jj1[y_] = k; jj2[py_] = k;
	for (l = 0; l <= no_tps; l++) {
	  jj1[py_] = l; jj2[y_] = l;
	  if ((i+j+k+l <= no_tps) && ((i-j != 0) || (k-l != 0))) {
	    cotan = 1e0/tan(((i-j)*nu_x+(k-l)*nu_y)*M_PI);
	    mn1 =
	      pow(Id[x_], i)*pow(Id[px_], j)*pow(Id[y_], k)*pow(Id[py_], l);
	    mn2 =
	      pow(Id[x_], j)*pow(Id[px_], i)*pow(Id[y_], l)*pow(Id[py_], k);

	    for (m = 0; m <= no_tps-i-j-k-l; m++) {
	      jj1[delta_] = m; jj2[delta_] = m;
	      re = h_re[jj1]; im = h_im[jj1];
	      // compute g
	      g_re += (re-cotan*im)*(mn1+mn2)*pow(Id[delta_], m)/2e0;
	      g_im += (im+cotan*re)*(mn1-mn2)*pow(Id[delta_], m)/2e0;
	      h_re.pook(jj2, 0e0); h_im.pook(jj2, 0e0);
	    }
	  }
	}
      }
    }
  }

  return RtoC(g_re, g_im);
}


ss_vect<tps> compute_map_normal_form(ss_vect<tps> &map)
{
  // Transform map into normal form:
  //   M = (exp(:g_no:) exp(:g_no-1:) ... exp(g_3) A_1 A_0)^-1
  //       exp(:K_2:) exp(:K_3:) ... exp(:K_no:)
  //       exp(:g_no:) exp(:g_no-1:) ... exp(g_3) A_1 A_0
  int          k, n;
  double       nu0[2];
  tps          h_k, h_k_re, h_k_im, h_ke, K, K_k, g, g_k;
  ss_vect<tps> Id, map1, R, A, nus, map2;

  Id.identity();

  danot_(no_tps-1);

  map1 = map;
  R = get_map_Fl(map1);

  danot_(no_tps);

  K = 0e0;
  for (k = 0; k < 2; k++) {
    nu0[k] = atan2(R[2*k][2*k+1], R[2*k][2*k]);
    if (nu0[k] < 0e0) nu0[k] += 2e0*M_PI;
    nu0[k] /= 2e0*M_PI;
    K -= M_PI*nu0[k]*(sqr(Id[2*k])+sqr(Id[2*k+1]));
  }
  cout << endl;
  cout << fixed << setprecision(5)
       << "nu0 = (" << nu0[X_] << ", " << nu0[Y_] << ")" << endl;

  // coasting beam
  K += h_ijklm(map1[ct_], 0, 0, 0, 0, 1)*sqr(Id[delta_])/2e0;
//  CtoR(K, h_k_re, h_k_im);
//  cout << endl << "K:" << h_k_re;

  g = 0e0;
  for (k = 3; k <= no_tps; k++) {
    n = pow(2, k-3);

    map2 = map1*Inv(R*compute_Dragt_Finn_Map(K, Id, 3, k-1, true));
    h_k =
      Intd(get_map_k(map2, k-1, k-1), -1e0);
    g_k = get_g(nu0[X_], nu0[Y_], h_k);
    g += g_k;
    CtoR(h_k, h_k_re, h_k_im);
    K_k = RtoC(get_Ker(h_k_re), get_Ker(h_k_im));
    K += K_k;
//    cout << endl << "k = " << k << h_k_re;

    A = compute_Dragt_Finn_Map(g_k, Id, k, k, true);
    map1 = Inv(A)*map1*A;
  }

#if 1
  MNF = MapNorm(map, no_tps);

  daeps_(1e-8);
  cout << "\nMNF.g-g:\n" << MNF.g-g;
  cout << "\nMNF.K-K:\n" << MNF.K-K;
  daeps_(eps_tps);
#endif

  // MNF.g = g;
  // MNF.K = K;

  return map1;
}


ss_vect<tps> get_A_nl(const tps g)
{
  return compute_Dragt_Finn_Map(g, Id(), 3, no_tps, true);
}


ss_vect<tps> get_A_nl_inv(const tps g)
{
  return compute_Dragt_Finn_Map(-g, Id(), 3, no_tps, false);
}


tps get_H(const tps &g, const tps &K)
{
  return K*get_A_nl_inv(g);
}


tps g_renorm(const double nu0_x, const double nu0_y,
	     const double nu1_x, const double nu1_y,
	     const tps &g)
{
  // Renormalize g: (1-R^-1)^-1 * h.
  long int     jj1[ss_dim], jj2[ss_dim];
  int          i, j, k, l, m;
  double       re, im, cotan0, cotan1, cotan0_sqr;
  tps          h_re, h_im, g_re, g_im, G_re, G_im, mn1, mn2;
  ss_vect<tps> Id;

  CtoR(g, g_re, g_im);

  for (k = 0; k < ss_dim; k++) {
    jj1[k] = 0; jj2[k] = 0;
  }

  Id.identity(); G_re = 0e0; G_im = 0e0;
  for (i = 0; i <= no_tps; i++) {
    jj1[x_] = i; jj2[px_] = i;
    for (j = 0; j <= i; j++) {
      jj1[px_] = j; jj2[x_] = j;
      for (k = 0; k <= no_tps; k++) {
	jj1[y_] = k; jj2[py_] = k;
	for (l = 0; l <= no_tps; l++) {
	  jj1[py_] = l; jj2[y_] = l;

	  if (i+j+k+l <= no_tps) {
	    cotan0 = 1e0/tan(((i-j)*nu0_x+(k-l)*nu0_y)*M_PI);
	    cotan0_sqr = sqr(cotan0);
	    cotan1 = 1e0/tan(((i-j)*nu1_x+(k-l)*nu1_y)*M_PI);
	    mn1 =
	      pow(Id[x_], i)*pow(Id[px_], j)*pow(Id[y_], k)*pow(Id[py_], l);
	    mn2 =
	      pow(Id[x_], j)*pow(Id[px_], i)*pow(Id[y_], l)*pow(Id[py_], k);

	    for (m = 0; m <= no_tps; m++) {
	      if (i+j+k+l+m <= no_tps) {
		jj1[delta_] = m; jj2[delta_] = m;
		if ((i != j) || (k != l)) {
		  re = g_re[jj1]; im = g_im[jj1];

		  // Compute h.
		  h_re = (re+cotan0*im)*2e0/(1e0+cotan0_sqr);
		  h_im = (im-cotan0*re)*2e0/(1e0+cotan0_sqr);

		  // Renormalize g.
		  G_re += (h_re-cotan1*h_im)*(mn1+mn2)*pow(Id[delta_], m)/2e0;
		  G_im += (h_im+cotan1*h_re)*(mn1-mn2)*pow(Id[delta_], m)/2e0;
		  g_re.pook(jj2, 0e0); g_im.pook(jj2, 0e0);
		}
	      }
	    }
	  }
	}
      }
    }
  }
  return RtoC(G_re, G_im);
}


void find_res(const int n, const double nu_x, const double nu_y)
{
  int    j, k;
  double res, int_part;

  const double eps = 0.01;

  printf("\nfind_res:\n"
	 );
  for (j = 0; j <= n; j++)
    for (k = -n; k <= n; k++) {
      res = j*nu_x + k*nu_y;
      if (((j != 0) || (k != 0)) &&
	  fabs(modf(res, &int_part)) < eps)
	printf("  %1d*nu_x %c %1d*nu_y = %8.5f\n", j, (k < 0)? '-':'+',
	       abs(k), res);
    }
}


tps MNF_renorm(ss_vect<double> &ps_Fl, ss_vect<double> &ps_nl, MNF_struct &MNF,
	       const bool prt_res)
{
  int    k;
  double nu[2], scl;
  tps    g_r;

  // Scale action-angle coordinates.
  for (k = 0; k < 2; k++) {
    nu[k] = (MNF.nus[k]*ps_nl).cst();
    scl = sqrt(MNF.nus[k].cst()/nu[k]);
    ps_Fl[2*k] *= scl;
    ps_Fl[2*k+1] *= scl;
  }

  if (prt_res) {
    printf("\ntrack_H - Renormalising:\n  nu = [%7.5f, %7.5f] ->"
	   "\n  nu = [%7.5f, %7.5f]\n",
	   MNF.nus[0].cst(),MNF. nus[1].cst(), nu[X_], nu[Y_]);
    find_res(5, nu[X_], nu[Y_]);
  }

  g_r = g_renorm(MNF.nus[0].cst(), MNF.nus[1].cst(), nu[X_], nu[Y_], MNF.g);

  MNF.A_nl = get_A_nl(g_r);
  MNF.A_nl_inv = get_A_nl_inv(g_r);

  return g_r;
}


tps df_dJ(const int k, const tps &f, const ss_vect<double> ps_Fl)
{
  // d/dJ = (x * df / dx + p_x * df / dp_x) / (x^2 + p_x^2).
  // Care is required for denominator.
  ss_vect<tps> Id;

  Id.identity();
  return
    (Id[2*k]*Der(f, 2*k+1)+Id[2*k+1]*Der(f, 2*k+2))
    /(sqr(ps_Fl[2*k])+sqr(ps_Fl[2*k+1]));
}


tps df_dphi(const int k, const tps &f)
{
  // d/dphi = px*df/dx - x*df/dp_x.
  ss_vect<tps> Id;

  Id.identity();
  return Id[2*k+1]*Der(f, 2*k+1) - Id[2*k]*Der(f, 2*k+2);
}


tps dfg_dphi(const int k, const tps &f, const tps &g,
	     const ss_vect<double> ps_nl)

{
  // phi = -atan(g(x, p_x) / f(x, p_x))
  // is not defined for f(x, p_x) = 0
  // dphi = (g * df - f * dg) / (f^2 + g^2)
  // Care is required for denominator.
  return (g*df_dphi(k, f)-f*df_dphi(k, g))/(sqr(ps_nl[2*k])+sqr(ps_nl[2*k+1]));
}


void track_H(const string &file_name, const double x, const double y,
	     tps &H, const bool prt_res)
{
  long int        lastpos;
  int             j, k;
  double          h, twoJ[2], phi[2], twoJ2[2], phi2[2];
  tps             g_r;
  ss_vect<double> ps, ps_Fl, ps_nl;
  FILE            *outf;

  const int    n   = 500;
  const double A[] = {x, y};

  outf = file_write(file_name.c_str());

  danot_(no_tps);

  ps.zero();
  for (k = 0; k < 2; k++)
    ps[2*k] = A[k];

  ps_Fl = (MNF.A1_inv*MNF.A0_inv*ps).cst();
  ps_nl = (MNF.A_nl_inv*ps_Fl).cst();

  if (true)
    g_r = MNF_renorm(ps_Fl, ps_nl, MNF, prt_res);
  else
    g_r = MNF.g;

  H = get_H(g_r, MNF.K);

  for (j = 0; j < n; j++) {
    Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);

    if (lastpos == globval.Cell_nLoc) {
      ps_Fl = (MNF.A1_inv*ps).cst();
      ps_nl = (MNF.A_nl_inv*ps_Fl).cst();
      h = (H*ps_Fl).cst();

      for (k = 0; k < 2; k++) {
	twoJ[k] = sqr(ps_Fl[2*k]) + sqr(ps_Fl[2*k+1]);
	phi[k] = atan2(ps_Fl[2*k+1], ps_Fl[2*k]);
	twoJ2[k] = sqr(ps_nl[2*k]) + sqr(ps_nl[2*k+1]);
	phi2[k] = atan2(ps_nl[2*k+1], ps_nl[2*k]);
      }

      fprintf(outf,
	      "%3d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e"
	      " %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	      j, 1e3*ps_Fl[x_], 1e3*ps_Fl[px_], 1e3*ps_Fl[y_], 1e3*ps_Fl[py_],
	      1e6*twoJ[X_], phi[X_], 1e6*twoJ[Y_], phi[Y_],
	      1e6*twoJ2[X_], phi2[X_], 1e6*twoJ2[Y_], phi2[Y_], 1e6*fabs(h));
    }
  }

  fclose(outf);
}


void track_H(const MNF_struct &MNF)
{
  tps H;

  H = get_H(MNF.g, MNF.K);

  // One super period: BESSY-III/b3_sf_40Grad_JB.lat.
  track_H("track_H_1.dat", 0.1e-3,     0.1e-3, H, false);
  track_H("track_H_2.dat", 0.5e-3,     0.5e-3, H, false);
  track_H("track_H_3.dat", 1.0e-3,     1.0e-3, H, false);
  track_H("track_H_4.dat", 1.1e-3,     1.0e-3, H, false);
  track_H("track_H_5.dat", 1.2e-3,     1.0e-3, H, false);
  track_H("track_H_6.dat", 1.3e-3,     1.0e-3, H, false);
#if 1
  track_H("track_H_7.dat", 1.47e-3,    1.0e-3, H, false);
#else
  track_H("track_H_7.dat", 1.47005e-3, 1.0e-3, H, true);
#endif
}


void diff_map(const string &file_name, const int n_x, const int n_y,
	      const double x_max, const double y_max, MNF_struct &MNF)
{
  int             i, j, k;
  double          A[2], dphi2_dphi[2];
  tps             g_r;
  ss_vect<double> ps, ps_Fl, ps_nl;
  FILE            *outf;

  const int
    n[] = {n_x, n_y};
  const double
    A_max[] = {x_max, y_max},
    A_min[] = {1e-6, 1e-6},
    cut     = 2e0;

  outf = file_write(file_name.c_str());

  danot_(no_tps);

  for (i = -n[X_]; i <= n[X_]; i++) {
    A[X_] = i*A_max[X_]/n[X_];
    ps.zero();
    ps[x_] = (i != 0)? A[X_] : A_min[X_];
    for (j = -n[Y_]; j <= n[Y_]; j++) {
      A[Y_] = j*A_max[Y_]/n[Y_];
      ps[y_] = (j != 0)? A[Y_] : A_min[Y_];
      ps_Fl = (MNF.A1_inv*ps).cst();
      ps_nl = (MNF.A_nl_inv*ps_Fl).cst();

      if (true)
	g_r = MNF_renorm(ps_Fl, ps_nl, MNF, false);
      else
	g_r = MNF.g;

      for (k = 0; k < 2; k++)
	dphi2_dphi[k] =
	  (dfg_dphi(k, MNF.A_nl_inv[2*k], MNF.A_nl_inv[2*k+1], ps_nl)
	   *ps_Fl).cst();

      fprintf(outf, "\n  %9.5f %9.5f %12.5e %12.5e",
	      1e3*A[X_], 1e3*A[Y_],
	      (fabs(dphi2_dphi[X_]) < cut)? fabs(dphi2_dphi[X_]):NAN,
	      (fabs(dphi2_dphi[Y_]) < cut)? fabs(dphi2_dphi[Y_]):NAN);
    }
    fprintf(outf, "\n");
  }
  fclose(outf);
}


void compute_invariant(ss_vect<tps> &M)
{
  int          k;
  double       nu_int;
  tps          H_2, DH_2, H, DH;
  ss_vect<tps> Id, B;

  const double
    alpha[] = {Cell[0].Alpha[X_], Cell[0].Alpha[Y_]},
    beta[]  = {Cell[0].Beta[X_],  Cell[0].Beta[Y_]},
    gamma[] = {(1e0+sqr(alpha[X_]))/beta[X_], (1e0+sqr(alpha[Y_]))/beta[Y_]},
    eta_x   = Cell[0].Eta[X_],
    etap_x  = Cell[0].Etap[X_];

  Id.identity();

  printf("\n  alpha_x = %6.3f beta_x = %5.3f gamma_x = %5.3f"
	 " eta_x = %10.3e eta'_x = %10.3e\n"
	 "  alpha_y = %6.3f beta_y = %5.3f gamma_y = %5.3f\n",
	 alpha[X_], beta[X_], gamma[X_], eta_x, etap_x,
	 alpha[Y_], beta[Y_], gamma[Y_]);

  danot_(2);

  H_2 = 0e0;
  for (k = 0; k < 2; k++) {
    H_2 +=
      -M_PI*modf(globval.TotalTune[k], &nu_int)
      *(gamma[k]*sqr(Id[2*k])+2e0*alpha[k]*Id[2*k]*Id[2*k+1]
	+beta[k]*sqr(Id[2*k+1]));
  }
  H_2 += M[ct_][delta_]*sqr(Id[delta_])/2e0;

  B.identity();
  B[x_]  += eta_x*Id[delta_];
  B[px_] += etap_x*Id[delta_];
  B[ct_] += -etap_x*Id[x_] - eta_x*Id[px_];

  printf("\nLinear dispersion computed by numerical differentiation.\n");
  printf("\nB:\n");
  prt_lin_map(3, B);

  danot_(no_tps);

  printf("\nLinear dispersion computed by TPSA.\n");
  MNF = MapNorm(M, 1);

  cout << scientific << setprecision(3) << "\ng:\n" << MNF.g;
  printf("\nA0:\n");
  prt_lin_map(3, MNF.A0);

  H_2 = H_2*Inv(MNF.A0);

  daeps_(1e-10);
  H_2 = 1e0*H_2;
  cout << scientific << setprecision(3) << "\nH_2:\n" << H_2;
  daeps_(eps_tps);

  printf("\ne^-H_2:\n");
  prt_lin_map(3, LieExp(H_2, Id));

  // Restore max order after call to LieExp.
  danot_(2);

  DH_2 = H_2*M - H_2;
  daeps_(1e-13);
  DH_2 = 1e0*DH_2;
  cout << scientific << setprecision(3) << "\nH_2*M - H_2:\n" << DH_2;
  daeps_(eps_tps);

  danot_(no_tps);

  MNF = MapNorm(M, no_tps);

  H = get_H(MNF.g, MNF.K)*Inv(MNF.A0*MNF.A1);

  daeps_(1e-7);
  H = 1e0*H;
  cout << scientific << setprecision(3) << "\nH:\n" << H;
  daeps_(eps_tps);

  DH = H*M - H;
  daeps_(1e-7);
  DH = 1e0*DH;
  cout << scientific << setprecision(3) << "\nDH:\n" << DH;
  daeps_(eps_tps);
}


int f_p_k_cmplx_sgn_corr(const int n_dof, const long int jj[])
{
  int ord, k, sgn = 0;

  // Compute the sum of exponents for the momenta for the oscillating planes:
  ord = 0;
  for (k = 0; k < n_dof; k++)
    ord += jj[2*k+1];
  ord = (ord % 4);
  //  Sum_k c_ijkl x^i p_x^j y^k p_y^l
  //  j + l mod 4 = [0, 3: +1; 1, 2: -1]
  switch (ord) {
  case 0:
  case 3:
    sgn = 1;
    break;
  case 1:
  case 2:
    sgn = -1;
    break;
  default:
    printf("\n: undefined case %d\n", ord);
    break;
  }
  return sgn;
}


tps p_k_cmplx_sgn_corr(const int n_dof, const tps &a)
{
  // Correct sign for complex vs. real momenta p_k.
  //   q_k =  (h_q_k^+ + h_q_k^-) / 2
  // i p_k = -(h_q_k^+ - h_q_k^-) / 2
  char     name[name_len_for+1];
  int      j, n;
  long int ibuf1[bufsize], ibuf2[bufsize], jj[nv_tps];
  double   rbuf[bufsize];
  tps      b;

  // Adjust the sign for the momenta for the oscillating planes.
  a.exprt(rbuf, ibuf1, ibuf2, name);
  n = rbuf[0];
  for (j = 1; j <= n; j++) {
    dehash_(no_tps, nv_tps, ibuf1[j-1], ibuf2[j-1], jj);
    rbuf[j] *= f_p_k_cmplx_sgn_corr(n_dof, jj);
  }
  b.imprt(n, rbuf, ibuf1, ibuf2);
  return b;
}

ctps p_k_cmplx_sgn_corr(const int n_dof, const ctps &a)
{
  return
    ctps(p_k_cmplx_sgn_corr(n_dof, a.real()),
	 p_k_cmplx_sgn_corr(n_dof, a.imag()));
}


#if 0
// Obsolete.

void CtoR_JB2(const int n_dof, const tps &a, tps &a_re, tps &a_im)
{
  int          k;
  tps          b, c;
  ss_vect<tps> Id;

  Id.identity();

  b = p_k_cmplx_sgn_corr(n_dof, a);

  // q_k -> (q_k + p_k) / 2
  // p_k -> (q_k - p_k) / 2
  // Complex space:
  // q_k =   (h_q_k^+ + h_q_k^-) / 2
  // p_k = i (h_q_k^+ - h_q_k^-) / 2
  map.identity();
  for (k = 0; k < n_dof; k++) {
    map[2*k]   = (Id[2*k]+Id[2*k+1])/2e0;
    map[2*k+1] = (Id[2*k]-Id[2*k+1])/2e0;
  }
  b = b*map;

  // q_k -> p_k
  // p_k -> q_k
  // Complex space:
  // i (q_k -/+ i p_k) = (i q_k +/- p_k)
  map.identity();
  for (k = 0; k < n_dof; k++) {
    map[2*k]   = Id[2*k+1];
    map[2*k+1] = Id[2*k];
  }
  c = b*map;

  a_re = (b+c)/2e0;
  a_im = (b-c)/2e0;
}


tps RtoC_JB2(const int n_dof, tps &a_re, tps &a_im)
{
  int          k;
  tps          b;
  ss_vect<tps> Id, map;

  Id.identity();

  b = a_re + a_im;

  // q_k -> q_k + p_k
  // p_k -> q_k - p_k
  // Complex space:
  // h_q_k^+ = q_k - i h_p_k
  // h_q_k^- = q_k + i h_p_k
  map.identity();
  for (k = 0; k < n_dof; k++) {
    map[2*k]   = Id[2*k] + Id[2*k+1];
    map[2*k+1] = Id[2*k] - Id[2*k+1];
  }
  b = b*map;
  b = p_k_cmplx_sgn_corr(n_dof, b);
  return b;
}

#endif

ctps CtoR(const int n_dof, const ctps &a)
{
  // Cartesian to resonance basis:
  //   q_k =  (h_q_k^+ + h_q_k^-) / 2
  // i p_k = -(h_q_k^+ - h_q_k^-) / 2
  int          k;
  ss_vect<tps> Id, Zero;

  Id.identity();
  Zero.zero();
  map.identity();
  for (k = 0; k < n_dof; k++) {
    map[2*k]   = (Id[2*k]+Id[2*k+1])/2e0;
    map[2*k+1] = (Id[2*k]-Id[2*k+1])/2e0;
  }
  return p_k_cmplx_sgn_corr(n_dof, a)*css_vect(map, Zero);
}


ctps RtoC(const int n_dof, const ctps &a)
{
  // Resonance to Cartesian basis.
  // h_q_k^+ = q_k - i h_p_k
  // h_q_k^- = q_k + i h_p_k
  int          k;
  ss_vect<tps> Id, Zero;

  Id.identity();
  Zero.zero();
  map.identity();
  for (k = 0; k < n_dof; k++) {
    map[2*k]   = Id[2*k] + Id[2*k+1];
    map[2*k+1] = Id[2*k] - Id[2*k+1];
  }
  return p_k_cmplx_sgn_corr(n_dof, a*css_vect(map, Zero));
}


void tst_ctor(MNF_struct &MNF)
{
  tps          K_re, K_im, K_re_JB, K_im_JB, g_re, g_im, g_re_JB, g_im_JB;
  ss_vect<tps> Id;
  ctps         cK, cg;

  Id.identity();

  CtoR(MNF.K, K_re, K_im);
  CtoR(MNF.g, g_re, g_im);
  cK = CtoR(2, ctps(MNF.K, 0e0));
  cg = CtoR(2, ctps(0e0, MNF.g));

  cout << "\n[K_re-K_re_JB, K_im-K_im_JB]:\n"
       << K_re-cK.real() << K_im-cK.imag();
  daeps_(1e-7);
  cout << "\n[g_re-g_re_JB, g_im-g_im_JB]:\n"
       << g_re-cg.real() << g_im-cg.imag();
  daeps_(eps_tps);

  cK = RtoC(2, ctps(K_re, K_im));
  cg = RtoC(2, ctps(g_re, g_im));

  cout << "\nRtoC_JB(2, cK)-K:\n" << cK.real()-MNF.K << cK.imag();

  daeps_(1e-6);
  cout << "\nRtoC_JB(2, cg)-g:\n" << 1e0*cg.real() << cg.imag()-MNF.g;
  daeps_(eps_tps);
}


void compute_C_S_long(double &alpha_z, double &beta_z)
{
  alpha_z =
    -globval.Ascr[ct_][ct_]*globval.Ascr[delta_][ct_]
    - globval.Ascr[ct_][delta_]*globval.Ascr[delta_][delta_];
  beta_z = sqr(globval.Ascr[ct_][ct_]) + sqr(globval.Ascr[ct_][delta_]);
}


ss_vect<tps> compute_A_A_t(const int n_dof)
{
  int          k;
  double       alpha_z, beta_z;
  ss_vect<tps> Id, A_A_t;

  Id.identity();

  if (n_dof > 2) compute_C_S_long(alpha_z, beta_z);

  A_A_t.zero();
  for (k = 0; k < n_dof; k++) {
    if (k < 2) {
      A_A_t[2*k] = Cell[0].Beta[k]*Id[2*k] - Cell[0].Alpha[k]*Id[2*k+1];
      A_A_t[2*k+1] =
	-Cell[0].Alpha[k]*Id[2*k]
	+ (1e0+sqr(Cell[0].Alpha[k]))/Cell[0].Beta[k]*Id[2*k+1];
    } else {
      A_A_t[ct_] = beta_z*Id[ct_] - alpha_z*Id[delta_];
      A_A_t[delta_] = -alpha_z*Id[ct_]	+ (1e0+sqr(alpha_z))/beta_z*Id[delta_];
    }
  }

  return A_A_t;
}


void tst_moment(void)
{
  long int     lastpos;
  ss_vect<tps> M, A_A_t;

  const int
    n_dof = 3,
#if 0
    loc   = globval.Cell_nLoc;
#else
    loc   = 10;
#endif

  globval.Cavity_on = !false;
  globval.radiation = !false;

  danot_(1);

  Ring_GetTwiss(true, 0e0);
  printglob();

  danot_(no_tps-1);

  M.identity();
  Cell_Pass(0, loc, M, lastpos);
  printf("\nM:");
  prt_lin_map(3, M);

  A_A_t = compute_A_A_t(n_dof);
  printf("\nInitial A_A_t beta = [%5.3f, %5.3f]:",
	 A_A_t[x_][x_], A_A_t[y_][y_]);
  prt_lin_map(3, A_A_t);

  A_A_t = M*tp_S(n_dof, M*A_A_t);

  printf("\nFinal A_A_t beta = [%5.3f, %5.3f]:",
	 A_A_t[x_][x_], A_A_t[y_][y_]);
  prt_lin_map(3, A_A_t);

  A_A_t = compute_A_A_t(n_dof);
  printf("\nInitial A_A_t beta = [%5.3f, %5.3f]:",
	 A_A_t[x_][x_], A_A_t[y_][y_]);
  prt_lin_map(3, A_A_t);

  Cell_Pass(0, loc, A_A_t, lastpos);
  A_A_t = tp_S(n_dof, A_A_t);

  Cell_Pass(0, loc, A_A_t, lastpos);
  printf("\nFinal A_A_t beta = [%5.3f, %5.3f]:",
	 A_A_t[x_][x_], A_A_t[y_][y_]);
  prt_lin_map(3, A_A_t);
}


tps compute_twoJ(const int n_dof, const double eps[], const ss_vect<tps> &A_A_t)
{
  int          j, k;
  tps          twoJ;
  ss_vect<tps> Id, quad_form;

  const ss_vect<tps> omega = get_S(n_dof);

  Id.identity();

  quad_form = tp_S(n_dof, omega)*A_A_t*omega;
  twoJ = 0e0;
  for (j = 0; j < 2*n_dof; j++)
    for (k = 0; k < 2*n_dof; k++)
      twoJ += Id[j]*sqrt(eps[j/2])*quad_form[j][k]*sqrt(eps[k/2])*Id[k];
  return twoJ;
}


void compute_Twiss(const tps &twoJ, double alpha[], double beta[])
{
  long int jj[ss_dim];
  int      k;

  for (k = 0; k < ss_dim; k++)
    jj[k] = 0;
  for (k = 0; k < 2; k++) {
    jj[2*k]   = 1;
    jj[2*k+1] = 1;
    alpha[k] = twoJ[jj]/2e0;
    jj[2*k]   = 0;
    jj[2*k+1] = 0;
    jj[2*k+1] = 2;
    beta[k] = twoJ[jj];
    jj[2*k+1] = 0;
  }
}


void tst_twoJ(void)
{
  long int     lastpos;
  int          k;
  double       alpha[2], beta[2];
  tps          twoJ;
  ss_vect<tps> Id, M, M_inv, A_A_t;

  const int n_dof = 3,
#if 0
    loc   = globval.Cell_nLoc;
#else
  loc   = 10;
#endif

  Id.identity();

  globval.Cavity_on = !false;
  globval.radiation = !false;

  danot_(1);

  Ring_GetTwiss(true, 0e0);
  printglob();

  danot_(no_tps-1);

  M.identity();
  Cell_Pass(0, loc, M, lastpos);
  M_inv = Inv(M);
  printf("\nM:");
  prt_lin_map(3, M);

  danot_(no_tps);

  for (k = 0; k < n_dof; k++)
    globval.eps[k] = 1e0;

  twoJ = compute_twoJ(n_dof, globval.eps, compute_A_A_t(n_dof));
  compute_Twiss(twoJ, alpha, beta);
  printf("\nInitial 2J alpha = [%5.3f, %5.3f] beta = [%5.3f, %5.3f]:\n",
	 alpha[X_], alpha[Y_], beta[X_], beta[Y_]);
  cout << twoJ;

  twoJ = twoJ*Inv(M);

  compute_Twiss(twoJ, alpha, beta);
  printf("\nFinal 2J alpha = [%5.3f, %5.3f] beta = [%5.3f, %5.3f]:\n",
	 alpha[X_], alpha[Y_], beta[X_], beta[Y_]);
  cout << twoJ;
}


int main(int argc, char *argv[])
{
  long int    lastpos;
  ss_vect<tps> M;

  globval.H_exact    = false; globval.quad_fringe    = false;
  globval.Cavity_on  = false; globval.radiation      = false;
  globval.emittance  = false; globval.IBS            = false;
  globval.pathlength = false; globval.bpm            = 0;
  globval.Cart_Bend  = false; globval.dip_edge_fudge = true;
  globval.mat_meth   = false;

  // disable from TPSALib- and LieLib log messages
  idprset_(-1);

  if (false)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  if (!false) {
    danot_(1);
    Ring_GetTwiss(true, 0e0);
    printglob();
  }

  danot_(no_tps-1);

  if (false) {
    M.identity();
    Cell_Pass(0, globval.Cell_nLoc, M, lastpos);
    printf("\nM:");
    prt_lin_map(3, M);
    // cout << scientific << setprecision(3) << "\nM:\n" << M << "\n";

    MNF = MapNorm(M, no_tps);
    MNF.nus = dHdJ(MNF.K);
    MNF.A0_inv = Inv(MNF.A0);
    MNF.A1_inv = Inv(MNF.A1);
    MNF.A_nl = get_A_nl(MNF.g);
    MNF.A_nl_inv = get_A_nl_inv(MNF.g);
  }

  if (false)
    compute_invariant(M);

  if (false)
    track_H(MNF);

  if (false)
    diff_map("diffusion.out", 25, 25, 4e-3, 5e-3, MNF);

  if (false)
    tst_ctor(MNF);

  if (false)
    tst_moment();

  if (!false)
    tst_twoJ();
}
