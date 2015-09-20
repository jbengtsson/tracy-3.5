/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

   Note, the operators "operator*=()", etc. requires
   the use of a local variable

*/


#define danamlen  10
#define dafunlen   4

typedef char    danambuf[danamlen];
typedef char    funnambuf[dafunlen];
typedef int     iVector[ss_dim];

extern const int  nv_tps, nd_tps, iref_tps;
extern int        no_tps, ndpt_tps;
extern double     eps_tps;


void daeps_(const double eps) { eps_tps = eps; }


void danot_(const int no) {

  if (no != 1) {
    printf("danot_: max order exceeded %d (%d)\n", no, 1);
    exit_(0);
  }
}


void daini_(int no, int nv, int fio)
{

  eps_tps = 1e-25;
  if (no != 1) {
    printf("daini_: max order exceeded %d (%d)\n", no, 1);
    exit_(1);
  }
  if ((nv < 1) || (nv > ss_dim))
    printf("daini_: to many dimensions %d(%d)\n", nv, ss_dim);
}


void lieini(const int no, const int nv, const int nd2i)
{

  if ((nv < 1) || (nv > ss_dim))
    printf("lieini: max dim exceeded %d (%d)\n", nv, ss_dim);
}


void daall_(tps_buf &x, const int nd2, const char *daname,
	   const int no, const int nv)
{
  int  j;

  for (j = 0; j <= ss_dim; j++)
    x[j] = 0.0;
}


void dadal_(tps_buf &x, const int ss_dim) { }


void davar_(tps_buf &x, const double r, const int i)
{
  int  j;

  x[0] = r;
  for (j = 1; j <= ss_dim; j++)
    x[j] = 0.0;
  x[i] = 1.0;
}


void dacon_(tps_buf &x, const double r)
{
  int  j;

  x[0] = r;
  for (j = 1; j <= ss_dim; j++)
    x[j] = 0.0;
}


void dapek_(const tps_buf &x, const int jj[], double &r)
{
  int  i, nzero;

  nzero  = 0;
  for (i = 0; i < ss_dim; i++) {
    if (jj[i] == 0)
      nzero++;
    else if (jj[i] == 1)
      r = x[i+1];
    else
      printf("Invalid jj in dapek\n");
  }

  if (nzero == ss_dim)
    r = x[0];
  else if (nzero < ss_dim-1)
    printf("Invalid jj in dapek\n");
}


void dapok_(tps_buf &x, const int jj[], const double r)
{
  int  i, nzero;

  nzero = 0;
  for (i = 0; i < ss_dim; i++) {
    if (jj[i] == 0)
      nzero++;
    else if (jj[i] == 1)
      x[i+1] = r;
    else
      printf("Invalid jj in dapek\n");
  }
  if (nzero == ss_dim)
    x[0] = r;
  else if (nzero < ss_dim-1)
    printf("Invalid jj in dapek\n");
}


double getmat(const ss_vect<tps> &map, const int i, const int j)
{

  return map[i-1].ltps[j];
}


void putmat(ss_vect<tps> &map, const int i, const int j, const double r)
{

  map[i-1].ltps[j] = r;
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


void dacop_(const tps_buf &x, tps_buf &z)
{
  int  i;

  for (i = 0; i <= ss_dim; i++)
    z[i] = x[i];
}


void daadd_(const tps_buf &x, const tps_buf &y, tps_buf &z)
{
  int  i;

  for (i = 0; i <= ss_dim; i++)
    z[i] = x[i] + y[i];
}


void dasub_(const tps_buf &x, const tps_buf &y, tps_buf &z)
{
  int  i;

  for (i = 0; i <= ss_dim; i++)
    z[i] = x[i] - y[i];
}


void damul_(const tps_buf &x, const tps_buf &y, tps_buf &z)
{
  int      i;
  tps_buf  u;

  u[0] = x[0]*y[0];
  for (i = 1; i <= ss_dim; i++)
    u[i] = x[0]*y[i] + x[i]*y[0];
  dacop_(u, z);
}


void dafun_(const char *fun, const tps_buf &x, tps_buf &z);

void dadiv_(const tps_buf &x, const tps_buf &y, tps_buf &z)
{
  tps_buf  yinv;

  dafun_("INV ", y, yinv); damul_(x, yinv, z);
}


void dacad_(const tps_buf &x, const double y, tps_buf &z)
{
  int  i;

  z[0] = x[0] + y;
  for (i = 1; i <= ss_dim; i++)
    z[i] = x[i];
}


void dacsu_(const tps_buf &x, const double y, tps_buf &z)
{

  dacad_(x, -y, z);
}


void dacmu_(const tps_buf &x, const double y, tps_buf &z)
{
  int  i;

  for (i = 0; i <= ss_dim; i++)
    z[i] = x[i]*y;
}


void dasuc_(const tps_buf &x, const double y, tps_buf &z)
{
  tps_buf  negx;

  dacmu_(x, -1.0, negx); dacad_(negx, y, z);
}


void dacdi_(const tps_buf &x, const double y, tps_buf &z)
{

  dacmu_(x, 1.0/y, z);
}


void dadic_(const tps_buf &x, const double y, tps_buf &z)
{
  tps_buf  xinv;

  dafun_("INV ", x, xinv); dacmu_(xinv, y, z);
}


void dapos_(const tps_buf &x, tps_buf &z)
{

  if (x[0] < 0.0) dacmu_(x, -1.0, z);
}


void dasqr_(const tps_buf &x, tps_buf &z) { damul_(x, x, z); }


void dacma_(const tps_buf &x, const tps_buf &y, const double rb, tps_buf &z)
{
  tps_buf  x1;

  dacmu_(y, rb, x1);  /* x1=y*rb */
  daadd_(x, x1, z);   /* z =x+x1 */
}


void dalin_(const tps_buf &x, const double ra, const tps_buf &y,
	    const double rb, tps_buf &z)
{
  tps_buf  x1, x2;

  dacmu_(x, ra, x1); dacmu_(y, rb, x2); daadd_(x1, x2, z);
}


void dainv_(const tps_buf &x, tps_buf &z)
{
  int     i;
  double  a;

  a = -1.0/sqr(x[0]); z[0] = 1.0/x[0];
  for (i = 1; i <= ss_dim; i++)
    z[i] = a*x[i];
}


void dasqrt(const tps_buf &x, tps_buf &z)
{
  int     i;
  double  a;

  a = sqrt(x[0]); z[0] = a; a = 0.5/a;
  for (i = 1; i <= ss_dim; i++)
    z[i] = a*x[i];
}


void daexp(const tps_buf &x, tps_buf &z)
{
  int     i;
  double  a;

  a = exp(x[0]); z[0] = a;
  for (i = 1; i <= ss_dim; i++)
    z[i] = a*x[i];
}


void dalog(const tps_buf &x, tps_buf &z)
{
  int i;

  z[0] = log(x[0]);
  for (i = 1; i <= ss_dim; i++)
    z[i] = x[i]/x[0];
}


void dasin(const tps_buf &x, tps_buf &z)
{
  int     i;
  double  a;

  z[0] = sin(x[0]); a = cos(x[0]);
  for (i = 1; i <= ss_dim; i++)
    z[i] = a*x[i];
}


void dacos(const tps_buf &x, tps_buf &z)
{
  int     i;
  double  a;

  z[0] = cos(x[0]); a = -sin(x[0]);
  for (i = 1; i <= ss_dim; i++)
    z[i] = a*x[i];
}


void dasinh(const tps_buf &x, tps_buf &z)
{
  int     i;
  double  a;

  z[0] = sinh(x[0]); a = cosh(x[0]);
  for (i = 1; i <= ss_dim; i++)
    z[i] = a*x[i];
}


void dacosh(const tps_buf &x, tps_buf &z)
{
  int     i;
  double  a;

  z[0] = cos(x[0]); a = sinh(x[0]);
  for (i = 1; i <= ss_dim; i++)
    z[i] = a*x[i];
}


void datan(const tps_buf &x, tps_buf &z)
{
  tps_buf  c, s;

  dacos(x, c); dasin(x, s); dadiv_(s, c, z);
}


void daarctan(const tps_buf &x, tps_buf &z)
{
  int     i;
  double  a;

  a = x[0]; z[0] = atan(a); a = 1.0/(1.0+sqr(a));
  for (i = 1; i <= ss_dim; i++)
    z[i] = a*x[i];
}


void dafun_(const char *fun, const tps_buf &x, tps_buf &z)
{
  tps_buf  u;

  if (!strncmp(fun, "INV ", dafunlen))
    dainv_(x, u);
  else if (!strncmp(fun, "SQRT", dafunlen))
    dasqrt(x, u);
  else if (!strncmp(fun, "EXP ", dafunlen))
    daexp(x, u);
  else if (!strncmp(fun, "LOG ", dafunlen))
    dalog(x, u);
  else if (!strncmp(fun, "SIN ", dafunlen))
    dasin(x, u);
  else if (!strncmp(fun, "COS ", dafunlen))
    dacos(x, u);
  else if (!strncmp(fun, "SINH ", dafunlen))
    dasinh(x, u);
  else if (!strncmp(fun, "COSH ", dafunlen))
    dacosh(x, u);
  else if (!strncmp(fun, "TAN ", dafunlen))
    datan(x, u);
  else if (!strncmp(fun, "ATAN", dafunlen))
    daarctan(x, u);
  else {
    printf("dafun: illegal function %s\n", fun);
    exit_(0);
  }

  dacop_(u, z);
}


void dacct_(const ss_vect<tps> &x, const int i,
	    const ss_vect<tps> &y, const int j, ss_vect<tps> &z, const int k)
{
  int           l, m, n;
  ss_vect<tps>  u;

  for (l = 0; l < k; l++) {
    u[l].ltps[0] = x[l].ltps[0];
    for (m = 1; m <= j; m++)
      u[l].ltps[m] = 0.0;
    for (m = 0; m <= j; m++)
      for (n = 1; n <= i; n++)
        u[l].ltps[m] += x[l].ltps[n]*y[n-1].ltps[m];
  }
  for (l = 0; l < k; l++)
    z[l] = u[l];
}


void dainv_(const ss_vect<tps> &x, const int i, ss_vect<tps> &z, const int k)
{
  Matrix  mat;

  getlinmat(i, x, mat);  /* gets linear matrix from x */
  if (InvMat(i, mat))    /* inverses matrix of rank i */
    putlinmat(k, mat, z);/* puts linear matrix into z */
  else
    printf("dainv: map is singular\n");
}


void Rotmap(const int n, ss_vect<tps> &map, const Matrix &R)
{
  Matrix  mat;

  getlinmat(n, map, mat); MulRMat(n, mat, R); putlinmat(n, mat, map);
}


void daabs_(const tps_buf &x, double &r)
{
  int     k;

  r = 0.0;
  for (k = 0; k <= ss_dim; k++)
    r += fabs(x[k]);
}


void daabs2_(const tps_buf &x, double &r)
{
  int     k;

  r = 0.0;
  for (k = 0; k <= ss_dim; k++)
    r += sqr(x[k]);
}
