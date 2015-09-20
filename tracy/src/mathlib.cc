/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

*/

// missing in lstdc++
//template double std::__cmath_power<double>(double, unsigned);

double log(const int k) { return log((double)k); }


/* Local variables for DetMat: */
struct LOC_DetMat
{
  const Matrix  *a;
  bool    cross[ss_dim];
};

typedef int iv1[ss_dim];
typedef int iv2[ss_dim][2];
typedef double v1[ss_dim];

/* Local variables for InvMat: */
struct LOC_InvMat
{
  long    n;
  Matrix  *a;
  long    row, column;
  double  determ;
} ;

void iniranf(const long i)
{
  rseed0 = i; rseed = i;
}

#define k               19
#define c               656329L
#define m               100000001
void newseed(void)
{

  rseed0 = (k*rseed0+c) % m; rseed = (rseed0+54321) % m;
}

double ranf(void)
{
  /* Generate random number with rectangular distribution */
  rseed = (k*rseed+c) % m; return (rseed/1e8);
}

#undef k
#undef c
#undef m


void setrancut(const double cut)
{

  printf("\n");
  printf("setrancut: cut set to %3.1f\n", cut);

  normcut_ = cut;
}


#define maxiter         100 
double normranf(void)
{  
  int i, j;
  double f, w;

  j = 0;
  do {
    j++;
    w = 0.0;
    for (i = 1; i <= 12; i++)
      w += ranf();
    f = w - 6.0;
  }
  while (fabs(f) > fabs(normcut_) && j <= maxiter);

  if (j > maxiter)
    fprintf(stdout,"  *** fatal error in normranf\n");
  return f;
}
#undef maxiter


double dtor(const double d) { return (d*M_PI/180.0); }


double GetAngle(const double x, const double y)
{
  double z;
  
//  if (pi == 0e0)
//    fprintf(stdout,"** pi not initialized in GetAngle\n");
  if (x != 0e0)
    z = atan(y/x);
  else
    z = sgn(y)*M_PI/2.0;
  if (x >= 0.0)
    return z;
  if (y >= 0.0)
    z += M_PI;
  else
    z -= M_PI;
  return z;
}


void CopyVec(const int n, const Vector &a, Vector &b)
{
  int i;

  for (i = 0; i < n; i++)
    b[i] = a[i];
}


void AddVec(const int n, const Vector &a, Vector &b)
{
  int i;
  for (i = 0; i < n; i++)
    b[i] = a[i] + b[i];
}


void SubVec(int n, const Vector &a, Vector &b)
{
  int i;

  for (i = 0; i < n; i++)
    b[i] -= a[i];
}


double xabs(long n, Vector &x)
{
  long    i;
  double  sum;

  sum = 0.0;
  for (i = 0; i < n; i++)
    sum += sqr(x[i]);

  return sqrt(sum);
}


void UnitMat(const int n, Matrix &A)
{
  int i, j;

  for (i = 1; i <= n; i++) {
    for (j = 1; j <= n; j++) {
      if (i == j)
        A[i-1][j-1] = 1.0;
      else
        A[i-1][j-1] = 0.0;
    }
  }
}


void ZeroMat(const int n, Matrix &A)
{
  int i, j;

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++)
      A[i][j] = 0.0;
  }
}


void CopyMat(const int n, const Matrix &A, Matrix &B)
{
  int i;

  for (i = 0; i < n; i++)
    CopyVec(n, A[i], B[i]);
}


void AddMat(const int n, const Matrix &a, Matrix &b)
{
  int i;

  for (i = 0; i < n; i++)
    AddVec(n, a[i], b[i]);
}

void SubMat(const int n, const Matrix &a, Matrix &b)
{
  /*n : integer; VAR a, b : matrix*/
  int i;

  for (i = 0; i < n; i++)
    SubVec(n, a[i], b[i]);
}


void LinTrans(const int n, const Matrix &a, Vector &x)
{
  int i, j;
  Vector y;

  for (i = 0; i < n; i++) {
    y[i] = 0e0;
    for (j = 0; j < n; j++)
      y[i] += a[i][j]*x[j];
  }
  CopyVec(n, y, x);
}


void MulcMat(const int n, const double c, Matrix &A)
{
  int i,j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      A[i][j] = A[i][j]*c;
}


void MulLMat(const int n, const Matrix &A, Matrix &B)
{
  int i, j, k;
  double x;
  Matrix C;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      x = 0e0;
      for (k = 0; k < n; k++)
        x += A[i][k]*B[k][j];
      C[i][j] = x;
    }

  CopyMat(n, C, B);
}


void MulRMat(const int n, Matrix &A, const Matrix &B)
{
  int i, j, k;
  double x;
  Matrix C;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      x = 0e0;
      for (k = 0; k < n; k++)
        x += A[i][k]*B[k][j];
      C[i][j] = x;
    }

  CopyMat(n, C, A);
}


double TrMat(const int n, const Matrix &A)
{
  int i;
  double x;

  x = 0e0;
  for (i = 0; i < n; i++)
    x += A[i][i];
  return x;
}


void TpMat(const int n, Matrix &A)
{
  int i, j;
  double x;

  for (i = 1; i <= n; i++) {
    for (j = 0; j <= i-2; j++) {
      x = A[i-1][j];
      A[i-1][j] = A[j][i-1];
      A[j][i-1] = x;
    }
  }
}


void SwapSigmaMat(Matrix &A)
{
  int i, j;
  double x;

  for (i = 0; i <= 2; i++) {
    for (j = 0; j <= 2; j++) {
      x = A[2*i][2*j]; A[2*i][2*j] = A[2*i+1][2*j+1]; A[2*i+1][2*j+1] = x;
      x = A[2*i][2*j+1]; A[2*i][2*j+1] = A[2*i+1][2*j]; A[2*i+1][2*j] = x;
    }
  }
}


double GdetMat(long n, struct LOC_DetMat *LINK)
{
  double Result = 0.0;
  long k, sign;
  double det;

  if (n > 1)
  {
    det = 0e0;
    if ((n & 1) == 1)
      sign = 1;
    else
      sign = -1;
    for (k = 0; k < ss_dim; k++) {
      if (!LINK->cross[k]) {
        LINK->cross[k] = true;
        det += sign * (*LINK->a)[n - 1][k] * GdetMat(n - 1, LINK);
        LINK->cross[k] = false;
        sign = -sign;
      }
    }
    return det;
  }
  for (k = 0; k < ss_dim; k++) {
    if (!LINK->cross[k])
      Result = (*LINK->a)[n - 1][k];
  }
  return Result;
}


double DetMat(const int n, const Matrix &A_)
{
  struct LOC_DetMat V;
  long j;
  double d;

  V.a = &A_;
  if (n == 2)  /* simple case of a matrix of rank 2*/
    return ((*V.a)[0][0] * (*V.a)[1][1] - (*V.a)[0][1] * (*V.a)[1][0]);
  else if (n == 3) { /* simple case of a matrix of rank 3*/
    d = (*V.a)[0][0]
        * ((*V.a)[1][1] * (*V.a)[2][2] - (*V.a)[1][2] * (*V.a)[2][1]);
    d += (*V.a)[0][1]
         * ((*V.a)[1][2] * (*V.a)[2][0] - (*V.a)[1][0] * (*V.a)[2][2]);
    d += (*V.a)[0][2]
         * ((*V.a)[1][0] * (*V.a)[2][1] - (*V.a)[1][1] * (*V.a)[2][0]);
    return d;
  } else {
    for (j = 1; j <= ss_dim; j++) {
      if (j <= n)
        V.cross[j - 1] = false;
      else
        V.cross[j - 1] = true;
    }
    return (GdetMat(n, &V));
  }
}


void swap_(double *x, double *y, struct LOC_InvMat *LINK)
{
  double d;

  d  = *x;
  *x = *y;
  *y = d;
}

void Interchange(struct LOC_InvMat *LINK)
{
  long l, FORLIM;

  if (LINK->row == LINK->column)
    return;
  LINK->determ = -LINK->determ;
  FORLIM = LINK->n;
  for (l = 0; l < FORLIM; l++)
    swap_(&(*LINK->a)[LINK->row - 1][l],
	  &(*LINK->a)[LINK->column - 1][l], LINK);
}


bool InvMat(const int n_, Matrix &A_)
{
  struct LOC_InvMat V;
  bool Result = false;
  long i, j, k, l, l1;
  double amax, t, d;
  Matrix b;
  iv1 ipivot;
  iv2 index;
  v1 pivot;
  long FORLIM, FORLIM1;

  V.n = n_;
  V.a = &A_;

  /* if 2-square matrix */
  if (V.n == 2) {
    d = DetMat(2, *V.a);
    if (d != 0e0) {  /* non zero determinant */
      Result = true;
      b[0][0] = (*V.a)[1][1] / d;
      b[0][1] = -((*V.a)[0][1] / d);
      b[1][0] = -((*V.a)[1][0] / d);
      b[1][1] = (*V.a)[0][0] / d;
      CopyMat(V.n, b, *V.a);
    } else /* non iversible matrix */
      Result = false;
  } else { /* matrix with n greater than 2 */
    V.determ = 1.0;
    FORLIM = V.n;
    for (j = 0; j < FORLIM; j++)
      ipivot[j] = 0;
    i = 1;
    while (i <= V.n && V.determ != 0e0) {
      amax = 0e0;
      FORLIM = V.n;
      for (j = 1; j <= FORLIM; j++) {
        if (ipivot[j - 1] != 1) {
          FORLIM1 = V.n;
          for (k = 1; k <= FORLIM1; k++) {
            if (ipivot[k - 1] > 1) goto _L1;
            if (ipivot[k - 1] < 1) {
              if (fabs(amax) < fabs((*V.a)[j - 1][k - 1])) {
                V.row = j;
                V.column = k;
                amax = (*V.a)[j - 1][k - 1];
              }
            }
          }
        }
      }
      if (amax == 0e0) {
        Result = false;
        V.determ = 0e0;
      } else {
        Result = true;
        ipivot[V.column - 1]++;
        Interchange(&V);
        index[i - 1][0] = V.row;
        index[i - 1][1] = V.column;
        pivot[i - 1] = (*V.a)[V.column - 1][V.column - 1];
        V.determ *= pivot[i - 1];
        (*V.a)[V.column - 1][V.column - 1] = 1.0;
        FORLIM = V.n;
        for (l = 0; l < FORLIM; l++)
          (*V.a)[V.column - 1][l] /= pivot[i - 1];
        FORLIM = V.n;
        for (l1 = 0; l1 < FORLIM; l1++) {
          if (l1 + 1 != V.column) {
            t = (*V.a)[l1][V.column - 1];
            (*V.a)[l1][V.column - 1] = 0e0;
            FORLIM1 = V.n;
            for (l = 0; l < FORLIM1; l++)
              (*V.a)[l1][l] -= (*V.a)[V.column - 1][l] * t;
          }
        }
      }  /*else */
      i++;
    }  /*while*/
    if (V.determ != 0e0) {
      FORLIM = V.n;
      for (i = 1; i <= FORLIM; i++) {
        l = V.n - i + 1;
        if (index[l - 1][0] != index[l - 1][1]) {
          V.row = index[l - 1][0];
          V.column = index[l - 1][1];
          FORLIM1 = V.n;
          for (k = 0; k < FORLIM1; k++)
            swap_(&(*V.a)[k][V.row - 1], &(*V.a)[k][V.column - 1], &V);
        }
      }
    }
  }
_L1:
  return Result;
}


#define SWAP(a,b) {float temp=(a);(a)=(b);(b)=temp;}

bool InvMat2(double a[4][4])
{
  const int n = 4;

  int     indxc[n], indxr[n], ipiv[n];
  int     i, icol = 0, irow = 0, j, k, l, ll;
  double  big, dum, pivinv;

  for (j=0;j<n;j++)
    ipiv[j]=0;
  for (i=0;i<n;i++) {
    big=0.0;
    for (j=0;j<n;j++)
      if (ipiv[j] != 1)
	for (k=0;k<n;k++) {
	  if (ipiv[k] == 0) {
	    if (fabs(a[j][k]) >= big) {
	      big=fabs(a[j][k]);
	      irow=j;
	      icol=k;
	    }
	  } else if (ipiv[k] > 1) fprintf(stdout,"GAUSSJ: Singular Matrix-1");
				}
    ++(ipiv[icol]);
    if (irow != icol) {
      for (l=0;l<n;l++) SWAP(a[irow][l],a[icol][l])
			   //~ for (l=0;l<=n;l++) SWAP(b[irow][l],b[icol][l])
			   }
    indxr[i]=irow;
    indxc[i]=icol;
    if (a[icol][icol] == 0.0) fprintf(stdout,"GAUSSJ: Singular Matrix-2");
    pivinv=1.0/a[icol][icol];
    a[icol][icol]=1.0;
    for (l=0;l<n;l++) a[icol][l] *= pivinv;
    //~ for (l=0;l<n;l++) b[icol][l] *= pivinv;
    for (ll=0;ll<n;ll++)
      if (ll != icol) {
	dum=a[ll][icol];
	a[ll][icol]=0.0;
	for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
	//~ for (l=0;l<n;l++) b[ll][l] -= b[icol][l]*dum;
      }
  }
  for (l=n-1;l>=0;l--) {
    if (indxr[l] != indxc[l])
      for (k=0;k<n;k++)
	SWAP(a[k][indxr[l]],a[k][indxc[l]]);
  }
  return true;
}
#undef SWAP


void prtmat(const int n, const Matrix &A)
{
  int i, j;

  printf("matrix:\n");
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++)
      printf(" %14.6e", A[i][j]);
//      printf(" %24.16e", A[i][j]);
    putchar('\n');
  }
}

