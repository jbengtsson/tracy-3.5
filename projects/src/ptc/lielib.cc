#define NO 5

#include "tracy_lib.h"


// F. Klein's Erlangen Program.


int no_tps   = NO,
    ndpt_tps = 5;


tps get_tps_k(const tps &h, const int k)
{
  // Take in Forest's F77 LieLib.
  // Get monomials of order k.
  long int no;
  tps      h_k;

  no = getno_();
  danot_(k-1);
  h_k = -h;
  danot_(k);
  h_k += h;
  danot_(no);
  return h_k;
}


ss_vect<tps> get_M_k(const ss_vect<tps> &x, const int k)
{
  // Taked in Forest's F77 LieLib.
  int          i;
  ss_vect<tps> map_k;

  for (i = 0; i < nv_tps; i++)
    map_k[i] = get_tps_k(x[i], k);
  return map_k;
}

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


void prt_cmplx_lin_map(const string &str, const css_vect &map)
{
  int          i, j;
  ss_vect<tps> re, im;

  for (i = 0; i < 2*nd_tps; i++) {
    re[i] = map[i].real();
    im[i] = map[i].imag();
  }

  cout << str;
  for (i = 1; i <= 2*nd_tps; i++) {
    for (j = 1; j <= 2*nd_tps; j++)
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

    map2 = map1*Inv(R*FExpo(K, Id, 3, k-1, true));
    h_k =
      Intd(get_M_k(map2, k-1), -1e0);
    g_k = get_g(nu0[X_], nu0[Y_], h_k);
    g += g_k;
    CtoR(h_k, h_k_re, h_k_im);
    K_k = RtoC(get_Ker(h_k_re), get_Ker(h_k_im));
    K += K_k;
//    cout << endl << "k = " << k << h_k_re;

    A = FExpo(g_k, Id, k, k, true);
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
  return FExpo(g, Id(), 3, no_tps, true);
}


ss_vect<tps> get_A_nl_inv(const tps g)
{
  return FExpo(-g, Id(), 3, no_tps, false);
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


double f_p_k_cmplx_sgn_corr(const long int jj[])
{
  // Adjust the sign for the momenta for the oscillating planes.
  // Correct sign for complex vs. real momenta p_k.
  //   q_k =  (h_q_k^+ + h_q_k^-) / 2
  // i p_k = -(h_q_k^+ - h_q_k^-) / 2
  // Adjust the sign for the momenta for the oscillating planes.
  int ord, k, sgn = 0;

  // Compute the sum of exponents for the momenta for the oscillating planes:
  ord = 0;
  for (k = 0; k < nd_tps; k++)
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


tps tps_compute_function
(const tps &a, std::function<double (const long int [])> fun)
{
  // Dacfu in Forest's LieLib.
  char     name[name_len_for+1];
  int      k, n;
  long int ibuf1[bufsize], ibuf2[bufsize], jj[nv_tps];
  double   rbuf[bufsize];
  tps      b;

  a.exprt(rbuf, ibuf1, ibuf2, name);
  n = rbuf[0];
  for (k = 1; k <= n; k++) {
    dehash_(no_tps, nv_tps, ibuf1[k-1], ibuf2[k-1], jj);
    rbuf[k] *= fun(jj);
  }
  b.imprt(n, rbuf, ibuf1, ibuf2);
  return b;
}


ctps p_k_cmplx_sgn_corr(const ctps &a)
{
  return
    ctps(tps_compute_function(a.real(), f_p_k_cmplx_sgn_corr),
	 tps_compute_function(a.imag(), f_p_k_cmplx_sgn_corr));
}


#if 0
// Obsolete.

void CtoR_JB2(const tps &a, tps &a_re, tps &a_im)
{
  int          k;
  tps          b, c;
  ss_vect<tps> Id;

  Id.identity();

  b = p_k_cmplx_sgn_corr_fun(nd_tps, a);

  // q_k -> (q_k + p_k) / 2
  // p_k -> (q_k - p_k) / 2
  // Complex space:
  // q_k =   (h_q_k^+ + h_q_k^-) / 2
  // p_k = i (h_q_k^+ - h_q_k^-) / 2
  map.identity();
  for (k = 0; k < nd_tps; k++) {
    map[2*k]   = (Id[2*k]+Id[2*k+1])/2e0;
    map[2*k+1] = (Id[2*k]-Id[2*k+1])/2e0;
  }
  b = b*map;

  // q_k -> p_k
  // p_k -> q_k
  // Complex space:
  // i (q_k -/+ i p_k) = (i q_k +/- p_k)
  map.identity();
  for (k = 0; k < nd_tps; k++) {
    map[2*k]   = Id[2*k+1];
    map[2*k+1] = Id[2*k];
  }
  c = b*map;

  a_re = (b+c)/2e0;
  a_im = (b-c)/2e0;
}


tps RtoC_JB2(tps &a_re, tps &a_im)
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
  for (k = 0; k < nd_tps; k++) {
    map[2*k]   = Id[2*k] + Id[2*k+1];
    map[2*k+1] = Id[2*k] - Id[2*k+1];
  }
  b = b*map;
  b = p_k_cmplx_sgn_corr_fun(nd_tps, b);
  return b;
}

#endif

ctps CtoR(const ctps &a)
{
  // Cartesian to resonance basis:
  //   q_k =  (h_q_k^+ + h_q_k^-) / 2
  // i p_k = -(h_q_k^+ - h_q_k^-) / 2
  int          k;
  ss_vect<tps> Id, Zero;

  Id.identity();
  Zero.zero();
  map.identity();
  for (k = 0; k < nd_tps; k++) {
    map[2*k]   = (Id[2*k]+Id[2*k+1])/2e0;
    map[2*k+1] = (Id[2*k]-Id[2*k+1])/2e0;
  }
  return p_k_cmplx_sgn_corr(a)*css_vect(map, Zero);
}


ctps RtoC(const ctps &a)
{
  // Resonance to Cartesian basis.
  // h_q_k^+ = q_k - i h_p_k
  // h_q_k^- = q_k + i h_p_k
  int          k;
  ss_vect<tps> Id, Zero;

  Id.identity();
  Zero.zero();
  map.identity();
  for (k = 0; k < nd_tps; k++) {
    map[2*k]   = Id[2*k] + Id[2*k+1];
    map[2*k+1] = Id[2*k] - Id[2*k+1];
  }
  return p_k_cmplx_sgn_corr(a*css_vect(map, Zero));
}


void tst_ctor(MNF_struct &MNF)
{
  tps          K_re, K_im, K_re_JB, K_im_JB, g_re, g_im, g_re_JB, g_im_JB;
  ss_vect<tps> Id;
  ctps         cK, cg;

  Id.identity();

  CtoR(MNF.K, K_re, K_im);
  CtoR(MNF.g, g_re, g_im);
  cK = CtoR(ctps(MNF.K, 0e0));
  cg = CtoR(ctps(0e0, MNF.g));

  cout << "\n[K_re-K_re_JB, K_im-K_im_JB]:\n"
       << K_re-cK.real() << K_im-cK.imag();
  daeps_(1e-7);
  cout << "\n[g_re-g_re_JB, g_im-g_im_JB]:\n"
       << g_re-cg.real() << g_im-cg.imag();
  daeps_(eps_tps);

  cK = RtoC(ctps(K_re, K_im));
  cg = RtoC(ctps(g_re, g_im));

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


ss_vect<tps> compute_A_A_t(void)
{
  int          k;
  double       alpha_z, beta_z;
  ss_vect<tps> Id, A_A_t;

  Id.identity();

  if (nd_tps > 2) compute_C_S_long(alpha_z, beta_z);

  A_A_t.zero();
  for (k = 0; k < nd_tps; k++) {
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


tps compute_twoJ(const double eps[], const ss_vect<tps> &A_A_t)
{
  int          j, k;
  tps          twoJ;
  ss_vect<tps> Id, quad_form;

  const ss_vect<tps> omega = get_S(nd_tps);

  Id.identity();

  quad_form = tp_S(nd_tps, omega)*A_A_t*omega;
  twoJ = 0e0;
  for (j = 0; j < 2*nd_tps; j++)
    for (k = 0; k < 2*nd_tps; k++)
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
  int          k, loc;
  double       alpha[2], beta[2];
  tps          twoJ;
  ss_vect<tps> Id, M, M_inv, A_A_t;

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

  for (k = 0; k < nd_tps; k++)
    globval.eps[k] = 1e0;

  twoJ = compute_twoJ(globval.eps, compute_A_A_t());
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

//------------------------------------------------------------------------------

tps tps_fun
(const tps &a, std::function<double (const long int [])> fun)
{
  // Dacfu in Forest's F77 LieLib.
  // Multiplies mononials I_vec with function f(I_vec).
  char     name[name_len_for+1];
  int      k, n;
  long int ibuf1[bufsize], ibuf2[bufsize], jj[nv_tps];
  double   rbuf[bufsize];
  tps      b;

  a.exprt(rbuf, ibuf1, ibuf2, name);
  n = rbuf[0];
  for (k = 1; k <= n; k++) {
    dehash_(no_tps, nv_tps, ibuf1[k-1], ibuf2[k-1], jj);
    rbuf[k] *= fun(jj);
  }
  b.imprt(n, rbuf, ibuf1, ibuf2);
  return b;
}


double f_int_mon(const long int jj[])
{
  // Integrate monomials:
  //   scl = 1/(|I_vec|+1)
  int    k;
  double scl;

  scl = 0e0;
  for (k = 0; k < 2*nd_tps; k++)
    scl += jj[k];
  scl += 1e0;
  scl = 1e0/scl;
  return scl;
}


tps M_to_h(const ss_vect<tps> &map)
{
  // Intd in Forest's F77 LieLib.
  // E. Forest, M. Berz, J. Irwin "Normal Form Methods for Complicated
  // Periodic Systems: A Complete Solution Using Differential Algebra and Lie
  // Operators" Part. Accel. 24, 91-107 (1989):
  //   Eqs. (34)-(37).
  // Integrate monomials:
  //   M -> exp(:h:)
  int          k;
  tps          f_x, f_px, h;
  ss_vect<tps> Id;

  Id.identity();
  h = 0e0;
  for (k = 0; k < nd_tps; k++) {
    // Integrate monomials.
    f_x = tps_fun(map[2*k+1], f_int_mon)*Id[2*k];
    f_px = tps_fun(map[2*k], f_int_mon)*Id[2*k+1];
    h += f_x - f_px;
  }
  return h;
}


ss_vect<tps> h_to_v(const tps &h)
{
  // Difd in Forest's F77 LieLib:
  // Compute vector flow operator from Lie operator :h:
  //   v = Omega * [del_x H, del_px H]^T
  int          k;
  ss_vect<tps> v;

  for (k = 0; k < nd_tps; k++) {
    v[2*k+1] = Der(h, 2*k+1);
    v[2*k] = -Der(h, 2*k+2);
  }
  return v;
}


tps v_to_tps(const ss_vect<tps> &v, const tps &x)
{
  // Daflo in Forest's F77 LieLib.
  //   y = v * nabla * x
  int k;
  tps y;

  y = 0e0;
  for (k = 0; k < 2*nd_tps; k++)
    y += v[k]*Der(x, k+1);
  return y;
}


tps exp_v_to_tps(const ss_vect<tps> &v, const tps &x, const double eps,
	      const int n_max)
{
  // Expflo in Forest's F77 LieLib:
  //   y = exp(v*nabla) * x
  int    k;
  double eps1;
  tps    y_k, y;

  y_k = y = x;
  for (k = 1; k <= n_max; k++) {
    y_k = v_to_tps(v, y_k/k);
    y += y_k;
    eps1 = abs(y_k);
    if (eps1 < eps)
      break;
  }
  if (eps1 < eps)
    return y;
  else {
    printf("\n*** exp_v_to_tps: did not converge eps = %9.3e (eps = %9.3e)"
	   " n_max = %1d\n", eps1, eps, n_max);
    return NAN;
  }
}


tps exp_v_to_tps(const ss_vect<tps> &v, const tps &x, const int k1,
	      const int k2, const double scl, const bool reverse)
{
  // Facflo in Forest's F77 LieLib.
  // not reverse:
  //   y = exp(D_k1) * exp(D_k1+1) ...  * exp(D_k2) * x
  // reverse:
  //   y = exp(D_k2) * exp(D_k2-1) ... * exp(D_k1) * x
  int          k;
  tps          y;
  ss_vect<tps> v_k;

  const int n_max = 100; 

  y = x;
  if (!reverse) {
    for (k = k1; k <= k2; k++) {
      v_k = scl*get_M_k(v, k);
      y = exp_v_to_tps(v_k, y, eps_tps, n_max);
    }
  } else {
    for (k = k2; k >= k1; k--) {
      v_k = scl*get_M_k(v, k);
      y = exp_v_to_tps(v_k, y, eps_tps, n_max);
    }
  }
  return y;
}


ss_vect<tps>M_to_M_fact(const ss_vect<tps> &map)
{
  // Flofac in Forest's F77 LieLib.
  // Factor map:
  //   M = M_2 ... * M_n
  int          j, k;
  ss_vect<tps> map_lin, map_res, map_fact;

  map_lin = get_M_k(map, 1);
  map_res = map*Inv(map_lin);
  map_fact.zero();
  for (k = 2; k <= no_tps; k++) {
    map_fact += get_M_k(map_res, k);
    for (j = 0; j < 2*nd_tps; j++)
      map_res[j] = exp_v_to_tps(map_fact, map_res[j], k, k, -1e0, false);
  }
  return map_fact;
}


tps exp_h_to_tps(const tps &h, const tps &x, const double eps,
		   const int n_max)
{
  // Exp1d in Forest's F77 LieLib.
  //   y = exp(:h:) x
  return exp_v_to_tps(h_to_v(h), x, eps, n_max);
}


ss_vect<tps> exp_h_to_M(const tps &h, const ss_vect<tps> &x, const double eps,
		       const int n_max)
{
  // Expnd2 in Forest's F77 LieLib.
  //   Y = exp(:h:) X
  int          k;
  ss_vect<tps> y;

  y = x;
  for (k = 0; k < 2*nd_tps; k++)
    y[k] = exp_h_to_tps(h, y[k], eps, n_max);
  return y;
}


tps M_to_h_DF(const ss_vect<tps> &map)
{
  // Liefact in Forest's F77 LieLib.
  // A. Dragt, J. Finn "Lie Series and Invariant Functions for Analytic
  // Symplectic maps" J. Math. Phys. 17, 2215-2227 (1976).
  // Dragt-Finn factorization:
  //   M ->  M_lin * exp(:h_3:) * exp(:h_4:) ...  * exp(:h_n:)
  return M_to_h(M_to_M_fact(map));
}


ss_vect<tps> h_DF_to_M
(const tps &Lie_DF_gen, const ss_vect<tps> &x, const int k1, const int k2,
 const bool reverse)
{
  // Fexpo in Forest's F77 LieLib.
  // Compute map from Dragt-Finn factorisation:
  // not reverse:
  //   M = exp(:h_3:) * exp(:h_4:) ...  * exp(:h_n:) * X
  // reverse:
  //   M = exp(:h_n:) * exp(:h_no-1:) ... * exp(:h_3:) * X
  int          k;
  tps          h_k;
  ss_vect<tps> map;

  const int n_max = 100;

  map.identity();
  for (k = k2; k >= k1; k--) {
    h_k = get_tps_k(Lie_DF_gen, k);
    if (!reverse)
      map = map*exp_h_to_M(h_k, x, eps_tps, n_max);
    else
      map = exp_h_to_M(h_k, x, eps_tps, n_max)*map;
  }
  return map;
}


int main(int argc, char *argv[])
{
  long int    lastpos;
  tps          h, h_re, h_im, h_DF;
  ss_vect<tps> Id, M, M2, M_lin;

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

  Id.identity();

  if (!false) {
    danot_(1);
    Ring_GetTwiss(true, 0e0);
    printglob();
  }

  if (false) {
    danot_(no_tps);

    M.identity();
    Cell_Pass(0, globval.Cell_nLoc, M, lastpos);
    printf("\nM:");
    prt_lin_map(3, M);

    h = LieFact_DF(M, M_lin);
    cout << exp_h_to_M(h, Id, eps_tps, 100)-LieExp(h, Id) << "\n";
  }

  if(false) {
    danot_(no_tps);

    M.identity();
    Cell_Pass(0, globval.Cell_nLoc, M, lastpos);
    printf("\nM:");
    prt_lin_map(3, M);

    cout << M_to_h_DF(M)-LieFact_DF(M, M_lin) << "\n";
  }

  if(!false) {
    danot_(no_tps-1);

    M.identity();
    Cell_Pass(0, globval.Cell_nLoc, M, lastpos);
    printf("\nM:");
    prt_lin_map(3, M);

    danot_(no_tps);

    M_lin = get_M_k(M, 1);
    prt_lin_map(3, M_lin);

    h_DF = M_to_h_DF(M);
    cout << "\nh_DF:\n" << h_DF << "\n";

    M2 = h_DF_to_M(h_DF, Id, 3, no_tps, false)*M_lin;
    
    danot_(no_tps-1);

    cout << "\nh_M2-M:\n" << M2-M << "\n";
  }

//------------------------------------------------------------------------------

  if (false) {
    danot_(no_tps-1);

    M.identity();
    Cell_Pass(0, globval.Cell_nLoc, M, lastpos);
    printf("\nM:");
    prt_lin_map(3, M);
    // cout << scientific << setprecision(3) << "\nM:\n" << M << "\n";

    danot_(no_tps);

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
    tst_twoJ();

//------------------------------------------------------------------------------

}
