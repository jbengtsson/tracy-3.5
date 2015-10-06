/* Author:        Johan Bengtsson

    Definitions:  Interface to Fortran library for Truncated Power
		  Series Algebra.

    Note, linear case is a special case, see e.g. daexp

*/

extern int  no_tps, ndpt_tps;

bool    ini_tps = false, header = false, res_basis = false, stable = false;

unsigned short int  seq_tps = 0;       // sequence no for TPSA vector
//const int           n_max   = 150;     // max iterations for LieExp
const int           n_max   = 200;     // max iterations for LieExp

// Fortran strings are passed from C by: [str, strlen(str)].
const int   name_len_for = 10; // name length in FORTRAN library is 10.

int  bufsize; // Note, max no of monomials is (no+nv)!/(nv!*no!)
 

long int fact(long int n)
{
  if (n > 0)
    return n*fact(n-1);
  else if (n == 0)
    return 1;
  else {
    cout << "fact: neg. argument: " << n << endl;
    exit_(0);
    // avoid compiler warning
    return -1;
  }
}


long int nok(long int n, long int k)
{
  long int  j;
  double    u;

  u = 1.0;
  for (j = 0; j < k; j++)
    u *= (double)(n-j)/(double)(k-j);
  return (long int)(u+0.5);
}


#if NO > 1

double getmat(const ss_vect<tps> &map, const int i, const int j)
{
  int      k;
  double   r;
  iVector  jj;

  for (k = 0; k < nv_tps; k++)
    jj[k] = 0;

  jj[j-1] = 1;

  dapek_(map[i-1].intptr, jj, r);

  return r;
}


void putmat(ss_vect<tps> &map, const int i, const int j, const double r)
{
  int      k;
  iVector  jj;

  for (k = 0; k < nv_tps; k++)
    jj[k] = 0;

  if (j > 0) jj[j-1] = 1;

  dapok_(map[i-1].intptr, jj, r);
}


void getlinmat(const int nv, const ss_vect<tps> &map, Matrix &mat)
{
  int  j, k;

  for (j = 1; j <= nv; j++)
    for (k = 1; k <= nv; k++)
      mat[j-1][k-1] = getmat(map, j, k);
}


void putlinmat(const int nv, const Matrix &mat, ss_vect<tps> &map)
{
  /* Puts zeroes in constant part of da map */
  int j, k;

  for (j = 1; j <= nv; j++) {
    for (k = 0; k <= nv; k++) {
      if (k == 0)
        putmat(map, j, k, 0.0);
      else
        putmat(map, j, k, mat[j-1][k-1]);
    }
  }
}

#endif


// Interface to Fortran TPSA library

extern "C" {
  // for Fortran compability
//  void f_init(void);

//  void MAIN_() { cout << "call to MAIN_" << endl; }
}


void TPSA_Ini(void)
{

  cout << endl;
  cout << scientific << "initializing TPSA library: no = " << no_tps
       << ", nv = " << nv_tps << ", nd = " << nd_tps
       << ", ndpt = " << ndpt_tps << ", eps = " << eps_tps << endl;

  // Initialize Fortran I/O
//  f_init();

#if NO > 1
  // Initialize Lie-lib
  lieinit_(no_tps, nv_tps, nd_tps, ndpt_tps, iref_tps, 0);
#endif

  bufsize = nok(no_tps+nv_tps, nv_tps);

  ini_tps = true;
}


void TPSAEps(const double eps)
{ daeps_(eps); eps_tps = eps; }

tps::tps(void) {
  char  name[11];

  if (!ini_tps) TPSA_Ini();
  seq_tps++; intptr = 0; sprintf(name, "tps-%-5hu", seq_tps);
  daall_(intptr, 1, name, no_tps, nv_tps, name_len_for); dacon_(intptr, 0.0);
}


tps::tps(const double r)
{
  char  name[11];

  if (!ini_tps) TPSA_Ini();
  seq_tps++; intptr = 0; sprintf(name, "tps-%-5hu", seq_tps);
  daall_(intptr, 1, name, no_tps, nv_tps, name_len_for); dacon_(intptr, r);
}


tps::tps(const double r, const int i)
{
    char  name[11];

  if (!ini_tps) TPSA_Ini();
  seq_tps++; intptr = 0; sprintf(name, "tps-%-5hu", seq_tps);
  daall_(intptr, 1, name, no_tps, nv_tps, name_len_for);
  if (i == 0)
    dacon_(intptr, r);
  else
    davar_(intptr, r, i);
}


tps::tps(const double r, const int jj[])
{
    char  name[11];

  if (!ini_tps) TPSA_Ini();
  seq_tps++; intptr = 0; sprintf(name, "tps-%-5hu", seq_tps);
  daall_(intptr, 1, name, no_tps, nv_tps, name_len_for);
  dapok_(intptr, jj, r);
}


tps::tps(const tps &x) {
    char  name[11];

  if (!ini_tps) TPSA_Ini();
  seq_tps++; intptr = 0; sprintf(name, "tps-%-5hu", seq_tps);
  daall_(intptr, 1, name, no_tps, nv_tps, name_len_for);
  dacop_(x.intptr, intptr);
}


tps::~tps(void)
{ dadal_(intptr, 1); }


double tps::operator[](const int k) const
{
  int      i;
  iVector  jj;
  double   r;

  for (i = 0; i < nv_tps; i++)
    jj[i] = 0;
  jj[k] = 1;
  dapek_(intptr, jj, r);
  return(r);
}


double tps::operator[](const int jj[]) const
{
  double  r;

  dapek_(intptr, jj, r);
  return(r);
}


void tps::pook(const int jj[], const double r)
{ dapok_(intptr, jj, r); }

void tps::exprt(double rbuf[], int ibuf1[], int ibuf2[], char *name) const
{ daexp_(intptr, rbuf, ibuf1, ibuf2, name, name_len_for); }

void tps::imprt(const int n, double rbuf[],
		const int ibuf1[], const int ibuf2[])
{ rbuf[0] = n; daimp_(rbuf, ibuf1, ibuf2, intptr); }

tps& tps::operator=(const double r)
{ dacon_(intptr, r); return *this; }

tps& tps::operator+=(const double x)
{ dacad_(intptr, x, intptr); return *this; }

tps& tps::operator-=(const double x)
{ dacad_(intptr, -x, intptr); return *this; }

tps& tps::operator*=(const double x)
{ dacmu_(intptr, x, intptr); return *this; }

tps& tps::operator/=(const double x)
{ dacmu_(intptr, 1.0/x, intptr); return *this; }


tps& tps::operator=(const tps &x)
{ dacop_(x.intptr, intptr); return *this; }

tps& tps::operator+=(const tps &x)
{ daadd_(intptr, x.intptr, intptr); return *this; }

tps& tps::operator-=(const tps &x)
{ dasub_(intptr, x.intptr, intptr); return *this; }

tps& tps::operator*=(const tps &x)
{ damul_(intptr, x.intptr, intptr); return *this; }

tps& tps::operator/=(const tps &x)
{ dadiv_(intptr, x.intptr, intptr); return *this; }


tps sqrt(const tps &a)
{
  tps  b;

  dafun_("SQRT", a.intptr, b.intptr, name_len_for);
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
    cout << "pow: should never get here " << n << endl;
    exit_(0);
    // avoid compiler warning
    return 0.0;
  }
}


tps exp(const tps &a)
{
  tps  b;

  dafun_("EXP ", a.intptr, b.intptr, name_len_for);
  return b;
}


tps log(const tps &a)
{
  tps  b;

  dafun_("LOG ", a.intptr, b.intptr, name_len_for);
  return b;
}


#if false

tps sin(const tps &a)
{
  tps  b;

  dafun_("SIN ", a.intptr, b.intptr, name_len_for);
  return b;
}

tps cos(const tps &a)
{
  tps  b;

  dafun_("COS ", a.intptr, b.intptr, name_len_for);
  return b;
}

tps tan(const tps &a)
{
  tps  b;

  dafun_("TAN ", a.intptr, b.intptr, name_len_for);
  return b;
}

#else

tps sin_tps(const tps &a)
{
  int       n;
  long int  k;
  tps       r, b;

  b = 0.0; k = 1; r = a;
  for (n = 1; n <= no_tps; n += 2) {
    b += r/(double)k; k *= (n+1)*(n+2); r *= -sqr(a);
  }

  return b;
}

tps cos_tps(const tps &a)
{
  int       n;
  long int  k;
  tps       r, b;

  b = 0.0; k = 1; r = 1.0;
  for (n = 0; n <= no_tps; n += 2) {
    b += r/(double)k; k *= (n+1)*(n+2); r *= -sqr(a);
  }

  return b;
}

tps sin(const tps &a)
{
  // sin(a+b) = sin(a)*cos(b) + cos(a)*sin(b)
  double  cst;
  tps     b;

  cst = a.cst(); b = a - cst;

  return sin(cst)*cos_tps(b)+cos(cst)*sin_tps(b);
}

tps cos(const tps &a)
{
  // cos(a+b) = cos(a)*cos(b) - sin(a)*sin(b)
  double  cst;
  tps     b;

  cst = a.cst(); b = a - cst;

  return cos(cst)*cos_tps(b)-sin(cst)*sin_tps(b);
}

tps tan(const tps &a) { return sin(a)/cos(a); }

#endif


tps asin(const tps &a)
{
  tps  b;

  dafun_("ASIN", a.intptr, b.intptr, name_len_for);
  return b;
}


tps acos(const tps &a)
{
  tps  b;

  dafun_("ACOS", a.intptr, b.intptr, name_len_for);
  return b;
}


#if false

tps atan(const tps &a)
{
  tps  b;

  dafun_("ATAN", a.intptr, b.intptr);
  return b;
}

#else

tps atan(const tps &a)
{
  // arctan(a+b) to 8th order
  double  cst;
  tps     b, c;

  if (no_tps <= 8) {
    cst = a.cst(); b = a - cst;

    c = b/(1+sqr(cst)) - (cst*sqr(b))/pow(1.0+sqr(cst), 2) +
      ((-1+3.0*sqr(cst))*pow(b, 3))/(3.0*pow(1.0+sqr(cst), 3)) +
      ((cst-pow(cst, 3))*pow(b, 4))/pow(1+sqr(cst), 4) +
      ((1-10*sqr(cst)+5.0*pow(cst, 4))*pow(b, 5))/(5.0*pow(1+sqr(cst), 5))+
      ((-3*cst+10.0*pow(cst, 3)-3.0*pow(cst, 5))*pow(b, 6))/
      (3.0*pow(1+sqr(cst), 6)) +
      ((-1+21.0*sqr(cst)-35.0*pow(cst, 4)+7.0*pow(cst, 6))*pow(b, 7))/
      (7.0*pow(1+sqr(cst), 7)) +
      ((cst-7.0*pow(cst, 3)+7.0*pow(cst, 5)-pow(cst, 7))*pow(b, 8))/
      pow(1+sqr(cst), 8) + atan(cst);
  } else {
    cout << "atan: only defined to " << no_tps << "th order (" << no_tps
	 << ")" << endl;
    exit_(1);
  }

  return c;
}

#endif


tps atan2(const tps &b,const tps &a) {
  tps  c;

  if (a.cst() > 0.0)
    c = atan(b/a);
  else if (a.cst() == 0.0)
    if (b.cst() != 0.0)
      c = sgn(b.cst())*M_PI/2.0;
    else {
      cout << "atan2: 0/0 undefined" << endl;
      exit_(1);
    }
  else
    if (b.cst() >= 0.0)
      c = atan(b/a) + M_PI;
    else
      c = atan(b/a) - M_PI;
  return c;
}


#if false

tps sinh(const tps &a)
{
  tps  b;

  dafun_("SINH", a.intptr, b.intptr);
  return b;
}

tps cosh(const tps &a)
{
  tps  b;

  dafun_("COSH", a.intptr, b.intptr);
  return b;
}

#else

tps sinh_tps(const tps &a)
{
  int       n;
  long int  k;
  tps       r, b;

  b = 0.0; k = 1; r = a;
  for (n = 1; n <= no_tps; n += 2) {
    b += r/(double)k; k *= (n+1)*(n+2); r *= sqr(a);
  }

  return b;
}

tps cosh_tps(const tps &a)
{
  int       n;
  long int  k;
  tps       r, b;

  b = 0.0; k = 1; r = 1.0;
  for (n = 0; n <= no_tps; n += 2) {
    b += r/(double)k; k *= (n+1)*(n+2); r *= sqr(a);
  }

  return b;
}

tps sinh(const tps &a)
{
  // sinh(a+b) = sinh(a)*cosh(b) + cosh(a)*sinh(b)
  double  cst;
  tps     b;

  cst = a.cst(); b = a - cst;

  return sinh(cst)*cosh_tps(b)+cosh(cst)*sinh_tps(b);
}

tps cosh(const tps &a)
{
  // cos(a+b) = cos(a)*cos(b) - sin(a)*sin(b)
  double  cst;
  tps     b;

  cst = a.cst(); b = a - cst;

  return cosh(cst)*cosh_tps(b)+sinh(cst)*sinh_tps(b);
}

#endif


const double tps::cst(void) const
{
  int      i;
  iVector  jj;
  double   r;

  for (i = 0; i < nv_tps; i++)
    jj[i] = 0;
  dapek_(intptr, jj, r);
  return r;
}


double abs(const tps &a)
{
  double  r;

  daabs_(a.intptr, r);
  return r;
}


double abs2(const tps &a)
{
  double  r;

  daabs2_(a.intptr, r);
  return r;
}


#if NO > 1

void idprset(const int level)
{
  idprset_(level);
}


tps Der(const tps &a, const int k)
{
  tps  b;

  dader_(k, a.intptr, b.intptr);
  return b;
}


tps LieExp(const tps &H, const tps &x)
{
  tps  y;

  exp1d_(H.intptr, x.intptr, y.intptr, eps_tps, n_max);
  return y;
}


tps LieFlo(const ss_vect<tps> &H, const tps &x)
{
  int  i, Hintptrs[nv_tps];
  tps  y;

  for (i = 0; i < nv_tps; i++) {
    Hintptrs[i] = H[i].intptr;
  }

  daflo_(Hintptrs, x.intptr, y.intptr);
  return y;
}


ss_vect<tps> FExpo(const tps &H, const ss_vect<tps> &x,
		   const int k0, const int k1, const int k)
{
  int           i, xintptrs[nv_tps], mapintptrs[nv_tps];
  ss_vect<tps>  map;

  for (i = 0; i < nv_tps; i++) {
    xintptrs[i] = x[i].intptr; mapintptrs[i] = map[i].intptr;
  }
  fexpo_(H.intptr, xintptrs, mapintptrs, k0, k1, 1.0, k);
  for (i = 2*nd_tps; i < nv_tps; i++)
    map[i] = tps(0.0, i+1);
  return map;
}


ss_vect<tps> LieExp(const tps &H, const ss_vect<tps> &x)
{
  int           i, xintptrs[nv_tps], mapintptrs[nv_tps];
  ss_vect<tps>  map;

  for (i = 0; i < nv_tps; i++) {
    xintptrs[i] = x[i].intptr; mapintptrs[i] = map[i].intptr;
  }
  expnd2_(H.intptr, xintptrs, mapintptrs, eps_tps, n_max);
  for (i = 2*nd_tps; i < nv_tps; i++)
    map[i] = tps(0.0, i+1);
  return map;
}


ss_vect<tps> LieFlo(const ss_vect<tps> &H, const ss_vect<tps> &x)
{
  int           i, Hintptrs[nv_tps], xintptrs[nv_tps], mapintptrs[nv_tps];
  ss_vect<tps>  map;

  for (i = 0; i < nv_tps; i++) {
    Hintptrs[i] = H[i].intptr; xintptrs[i] = x[i].intptr;
    mapintptrs[i] = map[i].intptr;
  }
  daflod_(Hintptrs, xintptrs, mapintptrs);
  for (i = 2*nd_tps; i < nv_tps; i++)
    map[i] = tps(0.0, i+1);
  return map;
}


void CCT(const tps x[], const int n_x, const tps y[], const int n_y,
	 tps z[], const int n_z)
{
  int           i, xintptrs[n_x], yintptrs[n_y], zintptrs[n_z];

  for (i = 0; i < n_x; i++)
    xintptrs[i] = x[i].intptr;
  for (i = 0; i < n_y; i++)
    yintptrs[i] = y[i].intptr;
  for (i = 0; i < n_z; i++)
    zintptrs[i] = z[i].intptr;
  dacct_(xintptrs, n_x, yintptrs, n_y, zintptrs, n_z);
}


ss_vect<tps> Inv_Ext(const ss_vect<tps> &x)
{
  int           i, xintptrs[nv_tps], yintptrs[nv_tps];
  ss_vect<tps>  y;

  for (i = 0; i < nv_tps; i++) {
    xintptrs[i] = x[i].intptr; yintptrs[i] = y[i].intptr;
  }
  dainv_(xintptrs, nv_tps, yintptrs, nv_tps);
  return y;
}


ss_vect<tps> MTREE(const ss_vect<tps> &x)
{
  int           i, xintptrs[nv_tps], yintptrs[nv_tps];
  ss_vect<tps>  y;

  for (i = 0; i < nv_tps; i++) {
    xintptrs[i] = x[i].intptr; yintptrs[i] = y[i].intptr;
  }
  etmtree_(xintptrs, yintptrs);
  return y;
}


ss_vect<double> PPUSH(const ss_vect<tps> &x, ss_vect<double> &y)
{
  int              i, xintptrs[nv_tps];
  ss_vect<double>  z;

  for (i = 0; i < nv_tps; i++)
    xintptrs[i] = x[i].intptr;
  etppush2_(xintptrs, y, z);
  return z;
}

#endif


tps operator*(const tps &x, const ss_vect<tps> &y)
{
  int  i, xintptrs[nv_tps], y1intptrs[nv_tps], zintptrs[nv_tps];
  tps           z;
  ss_vect<tps>  y1;

  xintptrs[0] = x.intptr; zintptrs[0] = z.intptr; y1 = y;
  for (i = 2*nd_tps; i < nv_tps; i++)
    y1[i] = tps(0.0, i+1);
  for (i = 0; i < nv_tps; i++)
    y1intptrs[i] = y1[i].intptr;
  dacct_(xintptrs, 1, y1intptrs, nv_tps, zintptrs, 1);
  return z;
}


ss_vect<tps> operator*(const ss_vect<tps> &x, const ss_vect<tps> &y)
{
  int           i, xintptrs[nv_tps], yintptrs[nv_tps], zintptrs[nv_tps];
  ss_vect<tps>  z;

  for (i = 0; i < nv_tps; i++) {
    xintptrs[i] = x[i].intptr; yintptrs[i] = y[i].intptr;
    zintptrs[i] = z[i].intptr;
  }
#if NO == 1
  dacct_(x, nv_tps, y, nv_tps, z, nv_tps);
#else
  etcct_(xintptrs, yintptrs, zintptrs);
#endif
  return z;
}


ss_vect<tps> Inv(const ss_vect<tps> &x)
{
  int           i, xintptrs[nv_tps], yintptrs[nv_tps];
  ss_vect<tps>  y;

  for (i = 0; i < nv_tps; i++) {
    xintptrs[i] = x[i].intptr; yintptrs[i] = y[i].intptr;
  }
#if NO == 1
  dainv_(x, nv_tps, y, nv_tps);
#else
  etinv_(xintptrs, yintptrs);
#endif
  return y;
}


#if NO == 1

ss_vect<tps> PInv(const ss_vect<tps> &x, const iVector &jj)
{
  int           k, n;
  ss_vect<tps>  y, z;

  n = 0;
  for (k = 0; k < ss_dim; k++) {
    if (jj[k] != 0) {
      n++; y[n-1] = x[k];
    }
  }

  dainv_(y, n, z, n);

  n = 0; y.zero();
  for (k = 0; k < ss_dim; k++) {
    if (jj[k] != 0) {
      n++; y[n-1] = z[k];
    }
  }

  return y;
}

#else

ss_vect<tps> PInv(const ss_vect<tps> &x, const iVector &jj)
{
  int           i, xintptrs[nv_tps], yintptrs[nv_tps];
  ss_vect<tps>  y;

  for (i = 0; i < nv_tps; i++) {
    xintptrs[i] = x[i].intptr; yintptrs[i] = y[i].intptr;
  }
  etpin_(xintptrs, yintptrs, jj);
  return y;
}

#endif


void GoFix(const ss_vect<tps> &xy, ss_vect<tps> &a1, ss_vect<tps> &a1inv,
	   const int nord)
{
  int  i, xyintptrs[nv_tps], a1intptrs[nv_tps], a1invintptrs[nv_tps];

  for (i = 0; i < nv_tps; i++) {
    xyintptrs[i] = xy[i].intptr; a1intptrs[i] = a1[i].intptr;
    a1invintptrs[i] = a1inv[i].intptr;
  }
  gofix_(xyintptrs, a1intptrs, a1invintptrs, nord);
}

MNF_struct MapNorm(const ss_vect<tps> &map, const int no)
{
  int         i, mapintptrs[nv_tps], A1intptrs[nv_tps], A0intptrs[nv_tps];
  int         map_resintptrs[nv_tps];
  MNF_struct  MNF;

  for (i = 0; i < nv_tps; i++) {
    mapintptrs[i] = map[i].intptr;
    A0intptrs[i] = MNF.A0[i].intptr; A1intptrs[i] = MNF.A1[i].intptr;
    map_resintptrs[i] = MNF.map_res[i].intptr;
  }
  stable = mapnorm_(mapintptrs, MNF.g.intptr, A1intptrs, A0intptrs,
		    map_resintptrs, MNF.K.intptr, no);
  return MNF;
}


ss_vect<tps> MapNormF(const ss_vect<tps> &x, ss_vect<tps> &g, ss_vect<tps> &a2,
		      ss_vect<tps> &a1, ss_vect<tps> &xy,
		      const int nord, const int kpmax)
{
  int   i, xintptrs[nv_tps], gintptrs[nv_tps], a2intptrs[nv_tps];
  int   a1intptrs[nv_tps], xyintptrs[nv_tps], Kintptrs[nv_tps];
  ss_vect<tps>  K;

  for (i = 0; i < nv_tps; i++) {
    xintptrs[i] = x[i].intptr; gintptrs[i] = g[i].intptr;
    a2intptrs[i] = a2[i].intptr; a1intptrs[i] = a1[i].intptr;
    xyintptrs[i] = xy[i].intptr;  Kintptrs[i] = K[i].intptr;
  }
  stable = mapnormf_(xintptrs, gintptrs, a2intptrs, a1intptrs, xyintptrs,
		     Kintptrs, nord, kpmax);
  return K;
}


ss_vect<tps> dHdJ(const tps &H)
{
  int           i, nuintptrs[nv_tps];
  ss_vect<tps>  nu;

  for (i = 0; i < nv_tps; i++)
    nuintptrs[i] = nu[i].intptr;
  dhdj_(H.intptr, nuintptrs);
  return nu;
}


void CtoR(const tps &a, tps &a_re, tps &a_im)
{
  ctor_(a.intptr, a_re.intptr, a_im.intptr);
}


tps RtoC(const tps &a_re, const tps &a_im)
{
  tps  a;

  rtoc_(a_re.intptr, a_im.intptr, a.intptr);
  return a;
}


tps LieFact_DF(const ss_vect<tps> &xy, ss_vect<tps> &x)
{
  /* Dragt-Finn factorization:

       M = M_lin exp(:h_3:) exp(:h_4:) ... 

  */
  int  i, xyintptrs[nv_tps], xintptrs[nv_tps];
  tps  H;

  for (i = 0; i < nv_tps; i++) {
    xyintptrs[i] = xy[i].intptr; xintptrs[i] = x[i].intptr;
  }
  liefact_(xyintptrs, xintptrs, H.intptr);
  return H;
}


tps LieFact(const ss_vect<tps> &xy)
{
  /* Single exponent Dragt-Finn factorization:

       M = exp(:h_2:) exp(:h_3:) exp(:h_4:) ... 

  */
  return Intd(FlowFact(xy), -1.0); 
}


ss_vect<tps> FlowFact(const ss_vect<tps> &xy)
{
  int           i, xyintptrs[nv_tps], Vintptrs[nv_tps];
  ss_vect<tps>  V;

  for (i = 0; i < nv_tps; i++) {
    xyintptrs[i] = xy[i].intptr; Vintptrs[i] = V[i].intptr;
  }
  flofacg_(xyintptrs, Vintptrs, eps_tps);
  return V;
}


tps Intd(const ss_vect<tps> &V, double scl)
{
  int  i, Vintptrs[nv_tps];
  tps  H;

  for (i = 0; i < nv_tps; i++)
    Vintptrs[i] = V[i].intptr;
  intd_(Vintptrs, H.intptr, scl);
  return H;
}


ss_vect<tps> Difd(const tps &H, double scl)
{
  int           i, Vintptrs[nv_tps];
  ss_vect<tps>  V;

  for (i = 0; i < nv_tps; i++)
    Vintptrs[i] = V[i].intptr;
  difd_(H.intptr, Vintptrs, scl);
  return V;
}


tps PB(const tps &a, const tps &b)
{
  tps  c;

  etpoi_(a.intptr, b.intptr, c.intptr);
  return c;
}


tps Take(const tps &H, const int n)
{
  tps  Hn;

  take_(H.intptr, n, Hn.intptr);
  return Hn;
}


ss_vect<tps> Taked(const ss_vect<tps> &H, const int n)
{
  int           i, Hintptrs[nv_tps], Hnintptrs[nv_tps];
  ss_vect<tps>  Hn;

  for (i = 0; i < nv_tps; i++) {
    Hintptrs[i] = H[i].intptr; Hnintptrs[i] = Hn[i].intptr;
  }
  taked_(Hintptrs, n, Hnintptrs);
  return Hn;
}


istream& operator>>(istream &is, tps &a)
{
  char	  line[max_str], *token;
  int     i, n, no1, nv1;
  int     ibuf1[bufsize], ibuf2[bufsize], jj[ss_dim];
  double  rbuf[bufsize];

  const bool  debug_ = false, prt = false;

  if (debug_) {
//    darea77_(a.intptr, 8);
    darea_(a.intptr, 8);
    return is;
  }

  is.getline(line, max_str); is.getline(line, max_str);
  sscanf(line, "tpsa, NO =%d, NV =%d", &no1, &nv1);
  if (prt) cout << "no = " << no1 << ", nv = " << nv1 << endl;
  ibuf1[0] = no_tps; ibuf2[0] = ss_dim;

  if ((no1 <= no_tps) && (nv1 <= ss_dim)) {
    for (i = 1; i <= 5; i++)
      is.getline(line, max_str);

    n = 0;
    do {
      n++;
      is.getline(line, max_str);
      token = strtok(line, " "); sscanf(token, "%d", &no1);
      token = strtok(NULL, " "); sscanf(token, "%le", &rbuf[n]);
      for (i = 0; i < ss_dim; i++) {
	token = strtok(NULL, " "); sscanf(token, "%d", &jj[i]);
      }
      if (prt) {
	cout << scientific << setprecision(3)
	     << no1 << setw(11) << rbuf[n];
	for (i = 0; i < ss_dim; i++)
	  cout << setw(3) << jj[i];
	cout << endl; 
      }

      hash_(no_tps, ss_dim, jj, ibuf1[n-1], ibuf2[n-1]);
    } while (no1 >= 0);

    rbuf[0] = -no1;
    daimp_(rbuf, ibuf1, ibuf2, a.intptr);
  } else {
    cout << "*** illegal no = " << no1 << " (" << no_tps
	 << ") or nv = " << nv1 << " (" << ss_dim << ")" << endl;
    exit_(1);
  }

  return is;
}


ostream& operator<<(ostream &os, const tps &a)
{
  char           name[11];
  int            i, j, ord, n, no;
  int            ibuf1[bufsize], ibuf2[bufsize], jj[nv_tps];
  double         rbuf[bufsize];
  ostringstream  s;

  const bool  debug_ = false;

  if (debug_) {
    dapri_(a.intptr, 6);
    return os;
  }

  daexp_(a.intptr, rbuf, ibuf1, ibuf2, name, name_len_for);
  s << endl;
  
  name[10] = '\0'; i = 0;
  while ((i <= 9) && (name[i] != ' ')) {
    s << name[i]; i++;
  }
  n = (int) rbuf[0];
  s << ", NO = " << no_tps
    << ", NV = " << nv_tps << ", INA = " << a.intptr << endl;

  for (i = 1; i <= 66; i++)
    s << "-"; 
  s << endl;

  if (header) {
    s << endl;
    if (!res_basis) {
      s << "                                                        n"
	<< endl;
      s << "      ====     i  i   i  i  i   i  i     i             ===="
	<< endl;
      s << "      \\         1  2   3  4  5   6  7     n            \\   "
	<< endl;
      s << "  P =  |   a  x  p   y  p  d  ct  p ... p  ,    |I| =  |   i"
	<< endl;
      s << "      /     I     x      y         1     n             /     k"
	<< endl;
      s << "      ====                                             ===="
	<< endl;
      s << "       I                                               k=1"
	<< endl;
    } else {
      s << "                                                          n"
	<< endl;
      s << "      ====      i   i   i   i  i   i  i     i            ===="
	<< endl;
      s << "      \\        + 1 - 2 + 3 - 4  5   6  7     n           \\   "
	<< endl;
      s << "  P =  |   a  h   h   h   h   d  ct  p ... p  ,    |I| =  |   i"
	<< endl;
      s << "      /     I  x   x   y   y          1     n            /     k"
	<< endl;
      s << "      ====                                               ===="
	<< endl;
      s << "       I                                                 k=1"
	<< endl;
    }
  }
  
  if (n != 0) {
    s << endl;
    s << "   |I|         a              ";
    for (i = 1; i <= nv_tps; i++)
      s << "  i";
    s << endl;
    s << "                I              ";
    for (i = 1; i <= nv_tps; i++)
      s << setw(3) << i;
    s << endl;
    s << endl;
  } else
    s << "   ALL COMPONENTS ZERO " << endl;
  for (no = 0; no <= no_tps; no++) {
    for (i = 1; i <= n; i++) {
      dehash_(no_tps, nv_tps, ibuf1[i-1], ibuf2[i-1], jj);
      ord = 0;
      for (j = 0; j < nv_tps; j++) 
	ord += jj[j];
      if (ord == no)
	if (fabs(rbuf[i]) >= eps_tps) {
          s << setw(5) << ord << scientific << setw(24)
            << setprecision(16) << rbuf[i] << " ";
	  for (j = 0; j < nv_tps; j++)
	    s << setw(3) << jj[j];
	  s << endl;
	}
    }
  }
  if (n == 0) n = 1;
  s << setw(5) << -n
    << scientific << setw(24) << setprecision(16) << 0.0 << " ";
  for (j = 0; j < nv_tps; j++)
    s << setw(3) << 0;
  s << endl;

  return os << s.str();
}
