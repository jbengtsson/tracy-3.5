/*

  Author:        Johan Bengtsson

  Definitions:  Truncated Power Series Algebra (to first order).

   Note, the operators "operator*=()", etc. requires
   the use of a local variable

*/


extern const int nv_tps, nd_tps, iref_tps;

extern int    no_tps, ndpt_tps;
extern double eps_tps;

bool ini_tps = false, header = false, res_basis = false, stable = false;

unsigned short int seq_tps = 0;       // sequence no for TPSA vector

const int n_max   = 100;              // max iterations for LieExp
const int bufsize = 250000;           /* Note, max no of monomials is
					 (no+nv)!/(nv!*no!) */
 

long int fact(long int n)
{
  if (n > 0)
    return n*fact(n-1);
  else if (n == 0)
    return 1;
  else {
    std::cout << "fact: neg. argument: " << n << "\n";
    exit_(1);
    // avoid compiler warning
    return -1;
  }
}


long int nok(long int n, long int k)
{
  long int j;
  double   u;

  u = 1.0;
  for (j = 0; j < k; j++)
    u *= (double)(n-j)/(double)(k-j);
  return (long int)(u+0.5);
}


// Interface to Fortran TPSA library

extern "C" {
  // for Fortran compability
//  void f_init(void);

//  void MAIN_() { cout << "call to MAIN_" << "\n"; }
}

void TPSA_Ini(void)
{
  long int n;

  std::cout << std::scientific << "\ninitializing linear TPSA library: no = "
	    << no_tps << ", nv = " << nv_tps << ", eps = " << eps_tps << "\n";
  // Initialize Fortran I/O
//  f_init();
  // Initialize TPSA-lib
  daini_(no_tps, nv_tps, 0); daeps_(eps_tps);
  // Initialize Lie-lib
//  lieinit_(no_tps, nv_tps, nd_tps, ndpt_tps, iref_tps, 0);
  n = nok(no_tps+nv_tps, nv_tps);
  if (n > bufsize) {
    std::cout << "*** bufsize exceeded " << n << " (" << bufsize << ")\n";
    exit_(0);
  }
  ini_tps = true;
}

void TPSAEps(const double eps)
{ daeps_(eps); eps_tps = eps; }

tps::tps(void) {
  char name[11];

  if (!ini_tps) TPSA_Ini();
  seq_tps++; sprintf(name, "tps-%-5hu", seq_tps);
  daall_(ltps, 1, name, no_tps, nv_tps); dacon_(ltps, 0.0);
}

tps::tps(const double r)
{
  char name[11];

  if (!ini_tps) TPSA_Ini();
  seq_tps++; sprintf(name, "tps-%-5hu", seq_tps);
  daall_(ltps, 1, name, no_tps, nv_tps); dacon_(ltps, r);
}

tps::tps(const double r, const int i)
{
  char name[11];

  if (!ini_tps) TPSA_Ini();
  seq_tps++; sprintf(name, "tps-%-5hu", seq_tps);
  daall_(ltps, 1, name, no_tps, nv_tps);
  if (i == 0)
    dacon_(ltps, r);
  else
    davar_(ltps, r, i);
}

tps::tps(const tps &x) {
  char name[11];

  if (!ini_tps) TPSA_Ini();
  seq_tps++; sprintf(name, "tps-%-5hu", seq_tps);
  daall_(ltps, 1, name, no_tps, nv_tps);
  dacop_(x.ltps, ltps);
}

tps::~tps(void) { dadal_(ltps, 1); }


double tps::operator[](const int k) const
{
  int      i;
  long int jj[nv_tps];
  double   r;

  for (i = 0; i < nv_tps; i++)
    jj[i] = 0;
  jj[k] = 1;
  dapek_(ltps, jj, r);
  return(r);
}

double tps::operator[](const long int jj[]) const
{
  double r;

  dapek_(ltps, jj, r);
  return(r);
}

void tps::pook(const long int jj[], const double r) { dapok_(ltps, jj, r); }

tps& tps::operator=(const double r) { dacon_(ltps, r); return *this; }

tps& tps::operator+=(const double x) { dacad_(ltps, x, ltps); return *this; }

tps& tps::operator-=(const double x) { dacad_(ltps, -x, ltps); return *this; }

tps& tps::operator*=(const double x) { dacmu_(ltps, x, ltps); return *this; }

tps& tps::operator/=(const double x)
{ dacmu_(ltps, 1.0/x, ltps); return *this; }


tps& tps::operator=(const tps &x) { dacop_(x.ltps, ltps); return *this; }

tps& tps::operator+=(const tps &x) { daadd_(ltps, x.ltps, ltps); return *this; }

tps& tps::operator-=(const tps &x) { dasub_(ltps, x.ltps, ltps); return *this; }

tps& tps::operator*=(const tps &x) { damul_(ltps, x.ltps, ltps); return *this; }

tps& tps::operator/=(const tps &x) { dadiv_(ltps, x.ltps, ltps); return *this; }


tps sqrt(const tps &a)
{
  tps b;

  dafun_("SQRT", a.ltps, b.ltps);
  return b;
}

tps pow(const tps &a, const int n)
{
  if (n < 0)
    return tps(pow(a, n+1)) /= a;
  else if (n == 0)
    return tps(1.0);
  else if (n == 1)
    return tps(a);
  else if (n > 1)
    return tps(pow(a, n-1)) *= a;
  else {
    std::cout << "pow: should never get here " << n << "\n";
    exit_(0);
    // avoid compiler warning
    return 0.0;
  }
}

tps exp(const tps &a)
{
  tps b;

  dafun_("EXP ", a.ltps, b.ltps);
  return b;
}

tps log(const tps &a)
{
  tps b;

  dafun_("LOG ", a.ltps, b.ltps);
  return b;
}

tps sin(const tps &a)
{
  tps b;

  dafun_("SIN ", a.ltps, b.ltps);
  return b;
}


tps cos(const tps &a)
{
  tps b;

  dafun_("COS ", a.ltps, b.ltps);
  return b;
}

tps tan(const tps &a)
{
  tps b;

  dafun_("TAN ", a.ltps, b.ltps);
  return b;
}

tps asin(const tps &a)
{
  tps b;

  dafun_("ASIN", a.ltps, b.ltps);
  return b;
}

tps atan(const tps &a)
{
  tps b;

  dafun_("ATAN", a.ltps, b.ltps);
  return b;
}


tps atan2(const tps &b, const tps &a) {
  tps c;

  if (a.cst() > 0.0)
    c = atan(b/a);
  else if (a.cst() == 0.0)
    if (b.cst() != 0.0)
      c = sgn(b.cst())*M_PI/2.0;
    else {
      std::cout << "atan2: 0/0 undefined\n";
      exit_(1);
    }
  else
    if (b.cst() >= 0.0)
      c = atan(b/a) + M_PI;
    else
      c = atan(b/a) - M_PI;
  return c;
}


tps sinh(const tps &a)
{
  tps b;

  dafun_("SINH", a.ltps, b.ltps);
  return b;
}

tps cosh(const tps &a)
{
  tps b;

  dafun_("COSH", a.ltps, b.ltps);
  return b;
}

const double tps::cst(void) const
{
  int      i;
  long int jj[nv_tps];
  double   r;

  for (i = 0; i < nv_tps; i++)
    jj[i] = 0;
  dapek_(ltps, jj, r);
  return r;
}

double abs(const tps &a)
{
  double r;

  daabs_(a.ltps, r);
  return r;
}

double abs2(const tps &a)
{
  double r;

  daabs2_(a.ltps, r);
  return r;
}


tps operator*(const tps &x, const ss_vect<tps> &y)
{
  ss_vect<tps> x1, z1;

  x1.zero(); x1[0] = x;
  dacct_(x1, nv_tps, y, nv_tps, z1, 1);
  return z1[0];
}

ss_vect<tps> operator*(const ss_vect<tps> &x, const ss_vect<tps> &y)
{
  ss_vect<tps> z;

  dacct_(x, nv_tps, y, nv_tps, z, nv_tps);
  return z;
}

ss_vect<tps> Inv(const ss_vect<tps> &x)
{
  ss_vect<tps> y;

  dainv_(x, nv_tps, y, nv_tps);
  return y;
}

ss_vect<tps> PInv(const ss_vect<tps> &x, const long int jj[])
{
  int          j, k;
  ss_vect<tps> Id, y, z;

  Id.identity(); y.zero();
  for (j = 0; j < nv_tps; j++)
    if (jj[j] != 0) {
      for (k = 0; k < nv_tps; k++)
	y[j] += jj[k]*x[j][k]*Id[k];
    } else
      y[j] = Id[j];

  dainv_(y, nv_tps, z, nv_tps);

  y.zero();
  for (j = 0; j < nv_tps; j++)
    if (jj[j] != 0) y[j] = z[j];

  return y;
}

std::istream& operator>>(std::istream &is, tps &a)
{
  char	 line[max_str], *token;
  int    i, n, no1, nv1;
//  int    ibuf1[bufsize], ibuf2[bufsize];
  int    jj[nv_tps];
  double rbuf[bufsize];

  const bool prt = false;

  std::cout << "not implemented\n"; exit_(1);

  is.getline(line, max_str); is.getline(line, max_str);
  sscanf(line, "tpsa, NO_TPSA =%d, NV =%d", &no1, &nv1);
  if (prt) std::cout << "no = " << no1 << ", nv = " << nv1 << "\n";
//  ibuf1[0] = no_tps; ibuf2[0] = nv_tps;

  if ((no1 <= no_tps) && (nv1 <= nv_tps)) {
    for (i = 1; i <= 5; i++)
      is.getline(line, max_str);

    n = 0;
    do {
      n++;
      is.getline(line, max_str);
      token = strtok(line, " "); sscanf(token, "%d", &no1);
      token = strtok(NULL, " "); sscanf(token, "%le", &rbuf[n]);
      for (i = 0; i < nv_tps; i++) {
	token = strtok(NULL, " "); sscanf(token, "%d", &jj[i]);
      }
      if (prt) {
	std::cout << std::scientific << std::setprecision(3)
	     << no1 << std::setw(11) << rbuf[n];
	for (i = 0; i < nv_tps; i++)
	  std::cout << std::setw(3) << jj[i];
	std::cout << "\n"; 
      }

//      hash_(no_tps, nv_tps, jj, ibuf1[n-1], ibuf2[n-1]);
    } while (no1 >= 0);

    rbuf[0] = -no1;
//    daimp_(rbuf, ibuf1, ibuf2, a.intptr);
  } else {
    std::cout << "*** illegal no (" << no_tps << ") or nv ("
	      << nv_tps << ")\n";
    exit_(1);
  }

  return is;
}


void prt_coeff(std::ostringstream &s, const tps &a, const long int jj[],
	       const int ord, int &n_coeff)
{
  int j;

  if (fabs(a[jj]) >= eps_tps) {
    s << std::setw(5) << ord << std::scientific << std::setw(24)
      << std::setprecision(16) << a[jj] << " ";
    for (j = 0; j < nv_tps; j++)
      s << std::setw(3) << jj[j];
    s << "\n";
    n_coeff++;
  }
}


void prt_header(std::ostream &s, const bool res_basis)
{
  s << "\n";
  if (!res_basis) {
    s << "                                                        n"
      << "\n      ====     i  i   i  i  i   i  i     i             ====\n"
      << "      \\         1  2   3  4  5   6  7     n            \\   \n"
      << "  P =  |   a  x  p   y  p  d  ct  p ... p  ,    |I| =  |   i\n"
      << "      /     I     x      y         1     n             /     k\n"
      << "      ====                                             ====\n"
      << "       I                                               k=1\n";
  } else {
    s << "                                                          n\n"
      << "      ====      i   i   i   i  i   i  i     i            ====\n"
      << "      \\        + 1 - 2 + 3 - 4  5   6  7     n           \\   \n"
      << "  P =  |   a  h   h   h   h   d  ct  p ... p  ,    |I| =  |   i\n"
      << "      /     I  x   x   y   y          1     n            /     k\n"
      << "      ====                                               ====\n"
      << "       I                                                 k=1\n";
  }
}


std::ostream& operator<<(std::ostream &os, const tps &a)
{
  int                i, j, n;
  long int           jj[nv_tps];
  std::ostringstream s;

  s << "\nNO_TPSA = " << no_tps << ", NV = " << nv_tps << "\n";

  for (i = 1; i <= 66; i++)
    s << "-"; 
  s << "\n";

  if (header) prt_header(s, res_basis);
  
  n = 0;
  for (j = 1; j <= nv_tps; j++)
    if (fabs(a.ltps[j]) >= eps_tps) n++;

  if (n != 0) {
    s << "\n   |I|         a              ";
    for (i = 1; i <= nv_tps; i++)
      s << "  i";
    s << "\n                I              ";
    for (i = 1; i <= nv_tps; i++)
      s << std::setw(3) << i;
    s << "\n\n";
  } else
    s << "   ALL COMPONENTS ZERO " << "\n";

  n = 0;
  for (j = 0; j < nv_tps; j++)
    jj[j] = 0;
  prt_coeff(s, a, jj, 0, n);
  for (i = 0; i < nv_tps; i++) {
    jj[i] = 1;
    prt_coeff(s, a, jj, 1, n);
    jj[i] = 0;
  }

  if (n == 0) n = 1;
  s << std::setw(5) << -n
    << std::scientific << std::setw(24) << std::setprecision(16) << 0e0 << " ";
  for (j = 0; j < nv_tps; j++)
    s << std::setw(3) << 0;
  s << "\n";

  return os << s.str();
}
