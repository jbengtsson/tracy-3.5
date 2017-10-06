/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

*/


template<typename T>
void spline(const double x[], const T y[], int const n,
	    double const yp1, const double ypn, T y2[])
{
  int    i, k;
  double sig;
  // Variable length arrays for user defined data types is not supported by
  // ANSI C++.
  // T      p, u[n], qn, un;
  T      p, qn, un;
  T      *u = new T[n];

  if (yp1 > 0.99e30)
    y2[1] = u[1] = 0.0;
  else {
    y2[1] = -0.5; u[1] = (3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
  }
  for (i = 2; i <= n-1; i++) {
    sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]); p = sig*y2[i-1]+2.0;
    y2[i] = (sig-1.0)/p;
    u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn = un = 0.0;
  else {
    qn = 0.5; un = (3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
  }
  y2[n] = (un-qn*u[n-1])/(qn*y2[n-1]+1.0);
  for (k = n-1; k >= 1; k--)
    y2[k] = y2[k]*y2[k+1]+u[k];

  delete [] u;
}


template<typename T, typename U>
void splint(const double xa[], const U ya[], const U y2a[],
	    const int n, const T &x, T &y)
{
  int     klo, khi, k;
  double  h;
  T       a, b;
  
  klo = 1; khi = n;
  while (khi-klo > 1) {
    k = (khi+klo) >> 1;
    if (xa[k] > x)
      khi = k;
    else
      klo = k;
  }
  h = xa[khi]-xa[klo];
  if (h == 0.0) nrerror("Bad xa input to routine splint");
  a = (xa[khi]-x)/h; b = (x-xa[klo])/h;
  y = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}


template<typename T>
void splin2(const double x1a[], const double x2a[],
	    double **ya, double **y2a, const int m, const int n,
	    const T &x1, const T &x2, T &y)
{
  int j;
  // Variable length arrays for user defined data types is not supported by
  // ANSI C++.
  // T   ytmp[m+1], yytmp[m+1];
  T   *ytmp = new T[m+1], *yytmp = new T[m+1];
  
  for (j = 1; j <= m; j++)
    splint(x2a, ya[j], y2a[j], n, x2, yytmp[j]);
  spline(x1a, yytmp, m, 1.0e30, 1.0e30, ytmp);
  splint(x1a, yytmp, ytmp, m, x1, y);

  delete [] ytmp; delete [] yytmp;
}


void splie2(double x1a[], double x2a[], double **ya,
	    int m, int n, double **y2a)
{
  int  j;

  for (j = 1; j <= m; j++)
    spline(x2a, ya[j], n, 1.0e30, 1.0e30, y2a[j]);
}


#define ZERO_RADIA 1e-7

void Read_IDfile(char *fic_radia, double &L, int &nx, int &nz,
                 double tabx[IDXMAX],  double tabz[IDZMAX],
                 double thetax[IDZMAX][IDXMAX], double thetaz[IDZMAX][IDXMAX],
		 bool &long_comp, double B2[IDZMAX][IDXMAX], int linear)
{
  FILE *fi;
  char dummy[5000];
  int  i, j;

  traceID = false;

  const double  Brho = globval.Energy*1e9/c0;

  /* open radia text file */
  if ((fi = fopen(fic_radia,"r")) == NULL) {
    printf("Read_IDfile: Error while opening file %s \n", fic_radia);
    exit_(1);
  }
  
  printf("\n");
  printf("Reading ID filename %s \n", fic_radia);
  printf("E      = %6.3f GeV\n", globval.Energy);
  printf("(Brho) = %6.3f\n", Brho);
  printf("linear = %1d\n", linear);
 
  /* first line */
  fscanf(fi, "%[^\n]\n", dummy); /* Read a full line */
  printf("%s\n", dummy);
  /* second line */
  fscanf(fi, "%[^\n]\n", dummy);
  printf("%s\n", dummy);
  /* third line */
  fscanf(fi, "%[^\n]\n", dummy);
  printf("%s\n", dummy);
  /* fourth line : Undulator length */
  fscanf(fi, "%lf\n", &L);
  printf("Insertion de longueur L = %lf m\n",L);
  /* fifth line */
  fscanf(fi, "%[^\n]\n", dummy);
  printf("%s\n", dummy);
  /* sisxth line : Number of Horizontal points */
  fscanf(fi, "%d\n", &nx);
  printf("nx = %d\n", nx);
  /* seventh line */
  fscanf(fi, "%[^\n]\n", dummy);
  printf("%s\n", dummy);
  /* Number of Vertical points */
  fscanf(fi, "%d\n", &nz);
  printf("nz = %d\n", nz);
  
  /* Check dimensions */
  if (nx > IDXMAX || nz > IDZMAX) {
    printf("Read_IDfile:  Increase the size of insertion tables \n");
    printf("nx = % d (IDXmax = %d) and nz = %d (IDZMAX = % d) \n",
	   nx, IDXMAX, nz, IDZMAX);
    exit_(1);
  }
  
  /* ninth line */
  fscanf(fi, "%[^\n]\n", dummy);
  printf("%s\n", dummy);
  /* tenth line */
  fscanf(fi, "%[^\n]\n", dummy);
//   printf("%s\n", dummy);
  
  for (j = 0; j < nx; j++)
    fscanf(fi, "%lf", &tabx[j]);
  fscanf(fi, "%[^\n]\n", dummy);
  
  /* Array creation for thetax */
  for (i = 0; i < nz; i++) {
    fscanf(fi, "%lf", &tabz[i]);
    
    for (j = 0; j < nx; j++) {
      fscanf(fi, "%lf", &thetax[i][j]);
      if (fabs(thetax[i][j]) < ZERO_RADIA)
	thetax[i][j] = 0.0;
      if (traceID) printf("%+12.8lf ", thetax[i][j]);
    }
    fscanf(fi, "\n");
    if (traceID) printf("\n");
  }
  
  fscanf(fi, "%[^\n]\n", dummy);
  printf("%s\n", dummy);
  fscanf(fi, "%[^\n]\n", dummy);
//   printf("%s\n", dummy);
  for (j = 0; j < nx; j++) {
    fscanf(fi, "%lf", &tabx[j]);
  }
  
  /* Array creation for thetaz */
  for (i = 0; i < nz; i++) {
    fscanf(fi, "%*f");
    for (j = 0; j < nx; j++) {
      fscanf(fi, "%lf", &thetaz[i][j]);
      if (fabs(thetaz[i][j]) < ZERO_RADIA)
	thetaz[i][j] = 0.0;
      if (traceID)
	printf("%+12.8lf ", thetaz[i][j]);
    }
    fscanf(fi, "\n");
    if (traceID) printf("\n");
  }
  
  /* Array creation for B2 */
  strcpy(dummy, "");
  fscanf(fi, "%[^\n]\n", dummy);
  printf("%s\n", dummy);
  fscanf(fi, "%[^\n]\n", dummy);
//  printf("%s\n", dummy);
  if (strstr(dummy, "START") == NULL) {
    long_comp = false;
    printf("no longitudinal component\n");
  } else {
    long_comp = true;
    printf("read longitudinal component\n");

    for (j = 0; j < nx; j++)
      fscanf(fi, "%lf", &tabx[j]);
  
    for (i = 0; i < nz; i++) {
      fscanf(fi, "%*f");
      for (j = 0; j < nx; j++) {
	fscanf(fi, "%lf", &B2[i][j]);
	B2[i][j] /= L*sqr(Brho);
	if (fabs(B2[i][j]) < ZERO_RADIA)
	  B2[i][j] = 0e0;
	if (traceID) printf("%+12.8lf ", B2[i][j]);
      }
      fscanf(fi, "\n");
      if (traceID) printf("\n");
    }
  }
  
  if (traceID)
    for (j = 0; j < nx; j++)
      printf("tabx[%d] =% lf\n", j, tabx[j]);
  if (traceID)
    for (j = 0; j < nz; j++)
      printf("tabz[%d] =% lf\n", j, tabz[j]);
  
  fclose(fi);
}



template<typename T>
void LinearInterpolation2(T &X, T &Z, T &TX, T &TZ, T &B2,
			  CellType &Cell, bool &out, int order)
{
  int            i, ix = 0, iz = 0;
  T              T1, U, THX = 0.0, THZ = 0.0;
  double         xstep = 0.0;
  double         zstep = 0.0;
  int            nx = 0, nz = 0;
  InsertionType  *WITH;
  
  WITH = Cell.Elem.ID; nx = WITH->nx; nz = WITH->nz;
  
  xstep = WITH->tabx[1]-WITH->tabx[0]; /* increasing values */
  zstep = WITH->tabz[0]-WITH->tabz[1]; /* decreasing values */
  
  if (traceID) printf("xstep = % f zstep = % f\n", xstep, zstep);
  
  /* test wether X and Z within the transverse map area */
  if (X < WITH->tabx[0] || X > WITH->tabx[nx-1]) {
    printf("LinearInterpolation2: X out of borders \n");
    printf("X = % lf but tabx[0] = % lf and tabx[nx-1] = % lf\n",
	   is_double<T>::cst(X), WITH->tabx[0], WITH->tabx[nx-1]);
    out = true;
    return;
  }
  
  if (Z > WITH->tabz[0] || Z < WITH->tabz[nz-1]) {
    printf("LinearInterpolation2: Z out of borders \n");
    printf("Z = % lf but tabz[0] = % lf and tabz[nz-1] = % lf\n",
	   is_double<T>::cst(Z),  WITH->tabz[0], WITH->tabz[nz-1]);
    out = true;
    return;
  }
  
  out = false;
  
  /* looking for the index for X */
  i = 0;
  while (X >= WITH->tabx[i]  && i <= nx-1) {
    i++;
    if (traceID)
      printf("%2d % lf % lf % lf \n",
	     i, is_double<T>::cst(X), WITH->tabx[i], WITH->tabx[i+1]);
  }
  ix = i - 1;
  
  /* looking for the index for Z */
  i = 0;
  while (Z <= WITH->tabz[i] && i <= nz-1) {
    i++;
    if (traceID)
      printf("%2d % lf % lf % lf \n",
	     i, is_double<T>::cst(Z), WITH->tabz[i], WITH->tabz[i+1]);
  }
  iz = i - 1;
  
  if (traceID) printf("Indices are ix=%d and iz=%d\n", ix, iz);
  
  /** Bilinear Interpolation **/
  U = (X - WITH->tabx[ix])/xstep; T1 = -(Z - WITH->tabz[iz])/zstep;
  
  if (order == 1) { // first order kick map interpolation
    if (traceID) printf("first order kick map interpolation\n");
    if (ix >= 0 && iz >= 0) {
      THX = (1.0-U)*(1.0-T1)*WITH->thetax1[iz][ix]
	    + U*(1.0-T1)*WITH->thetax1[iz][ix+1]
	    + (1.0-U)*T1*WITH->thetax1[iz+1][ix]
	    + U*T1*WITH->thetax1[iz+1][ix+1];
      
      THZ = (1.0-U)*(1.0-T1)*WITH->thetaz1[iz][ix]
	    + U*(1.0-T1)*WITH->thetaz1[iz][ix+1]
	    + (1.0-U)*T1*WITH->thetaz1[iz+1][ix]
	    + U*T1*WITH->thetaz1[iz+1][ix+1];
    }
    
    if (traceID) {
      printf("X=% f interpolation : U= % lf T =% lf\n",
	     is_double<T>::cst(X), is_double<T>::cst(U),
	     is_double<T>::cst(T1));
      printf("THX = % lf 11= % lf 12= %lf 21 = %lf 22 =%lf \n",
	     is_double<T>::cst(THX),
	     WITH->thetax1[iz][ix], WITH->thetax1[iz][ix+1],
	     WITH->thetax1[iz+1][ix], WITH->thetax1[iz+1][ix+1]);
      printf("Z=% f interpolation : U= % lf T =% lf\n",
	     is_double<T>::cst(Z), is_double<T>::cst(U),
	     is_double<T>::cst(T1));
      printf("THZ = % lf 11= % lf 12= %lf 21 = %lf 22 =%lf \n",
	     is_double<T>::cst(THZ),
	     WITH->thetaz1[iz][ix], WITH->thetaz1[iz][ix+1],
	     WITH->thetaz1[iz+1][ix],WITH->thetaz1[iz+1][ix+1]);
    }
  }
  
  if (order == 2) { // second order kick map interpolation
    if (traceID) printf("second order kick map interpolation\n");
    if (ix >= 0 && iz >= 0) {
      THX =
	(1.0-U)*(1.0-T1)*WITH->thetax[iz][ix]
	+ U*(1.0-T1)*WITH->thetax[iz][ix+1]
	+ (1.0-U)*T1*WITH->thetax[iz+1][ix]
	+ U*T1*WITH->thetax[iz+1][ix+1];
      
      THZ =
	(1.0-U)*(1.0-T1)*WITH->thetaz[iz][ix]
	+ U*(1.0-T1)*WITH->thetaz[iz][ix+1]
	+ (1.0-U)*T1*WITH->thetaz[iz+1][ix]
	+ U*T1*WITH->thetaz[iz+1][ix+1];

      if (WITH->long_comp)
	B2 =
	  (1.0-U)*(1.0-T1)*WITH->B2[iz][ix]
	  + U*(1.0-T1)*WITH->B2[iz][ix+1]
	  + (1.0-U)*T1*WITH->B2[iz+1][ix]
	  + U*T1*WITH->B2[iz+1][ix+1];
    }
    
    if (traceID) {
      printf("X=% f interpolation : U= % lf T =% lf\n",
	     is_double<T>::cst(X), is_double<T>::cst(U),
	     is_double<T>::cst(T1));
      printf("THX = % lf 11= % lf 12= %lf 21 = %lf 22 =%lf \n",
	     is_double<T>::cst(THX),
	     WITH->thetax[iz][ix], WITH->thetax[iz][ix+1],
	     WITH->thetax[iz+1][ix], WITH->thetax[iz+1][ix+1]);
      printf("Z=% f interpolation : U= % lf T =% lf\n",
	     is_double<T>::cst(Z), is_double<T>::cst(U),
	     is_double<T>::cst(T1));
      printf("THZ = % lf 11= % lf 12= %lf 21 = %lf 22 =%lf \n",
	     is_double<T>::cst(THZ),
	     WITH->thetaz[iz][ix], WITH->thetaz[iz][ix+1],
	     WITH->thetaz[iz+1][ix], WITH->thetaz[iz+1][ix+1]);
      printf("B2 = % lf 11= % lf 12= %lf 21 = %lf 22 =%lf \n",
	     is_double<T>::cst(THZ),
	     WITH->B2[iz][ix], WITH->B2[iz][ix+1],
	     WITH->B2[iz+1][ix],WITH->B2[iz+1][ix+1]);
    }
  }
  TX = THX; TZ = THZ;
}


/****************************************************************************/
/* void SplineInterpolation2(double X, double Z, double &TX, double &TZ,
                             CellType &Cell, bool &out)
 
   Purpose:
        Computes thx and thz in X and Z values using a bilinear interpolation
        interpolation of the array thetax(x, z) and thetaz(x, z)
 
   Input:
       X, Z location of the interpolation
       Cell elment containing ID device
 
   Output:
       TX, TZ thetax and thetaz interpolated at X and Z
       out true if interpolation out of table
 
   Return:
       none
 
   Global variables:
       none
 
   Specific functions:
 
   Comments:
       none
 
****************************************************************************/
template<typename T>
void SplineInterpolation2(T &X, T &Z, T &thetax, T &thetaz,
			  CellType &Cell, bool &out)
{
    int            nx, nz;
    InsertionType  *WITH;
//    int kx, kz;

    WITH = Cell.Elem.ID; nx = WITH->nx; nz = WITH->nz;

    /* test wether X and Z within the transverse map area */
    if (X < WITH->tabx[0] || X > WITH->tabx[nx-1] ||
	Z > WITH->tabz[0] || Z < WITH->tabz[nz-1]) {
        printf("SplineInterpDeriv2: out of borders in element s= %4.2f %*s\n",
	       Cell.S, 5, Cell.Elem.PName);
        printf("X = % lf but tabx[0] = % lf and tabx[nx-1] = % lf\n",
	       is_double<T>::cst(X), WITH->tabx[0], WITH->tabx[nx-1]);
        printf("Z = % lf but tabz[0] = % lf and tabz[nz-1] = % lf\n",
	       is_double<T>::cst(Z), WITH->tabz[0], WITH->tabz[nz-1]);
        out = true;
        return;
    }

    out = false;
    splin2(WITH->tab2-1, WITH->tab1-1, WITH->tx, WITH->f2x, nz, nx,
	   Z, X, thetax);
/*    if (fabs(temp) > ZERO_RADIA)
      *thetax = (double) temp;
    else
      *thetax = 0.0;*/
    splin2(WITH->tab2-1, WITH->tab1-1, WITH->tz, WITH->f2z, nz, nx,
	   Z, X, thetaz);
/*    if (fabs(temp) > ZERO_RADIA)
      *thetaz = (double) temp;
    else
      *thetaz = 0.0;*/

/*    FILE * fic0;
    char *fic="fit.out";
    fic0 = fopen(fic, "w");
    for (kz = 1; kz <= nz; kz++) {
      for (kx = 1; kx <= nx; kx++)
	fprintf(fic0, "% 12.3e", tz[kz][kx]);
      fprintf(fic0, "\n");
    }
    fclose(fic0);*/
}


void Matrices4Spline(InsertionType *WITH)
{
  int kx, kz;

  for (kx = 0; kx < WITH->nx; kx++) {
    WITH->tab1[kx] = (float) WITH->tabx[kx];
  }

  /** reordering: it has to be in increasing order */
  for (kz = 0; kz < WITH->nz; kz++) {
    WITH->tab2[kz] = (float) WITH->tabz[WITH->nz-kz-1];
  }

  for (kx = 0; kx < WITH->nx; kx++) {
    for (kz = 0; kz <WITH-> nz; kz++) {
      WITH->tx[kz+1][kx+1] = (float) (WITH->thetax[WITH->nz-kz-1][kx]);
      WITH->tz[kz+1][kx+1] = (float) (WITH->thetaz[WITH->nz-kz-1][kx]);
    }
  }

  // computes second derivative matrices
  splie2(WITH->tab2-1,WITH->tab1-1,WITH->tx,WITH->nz,WITH->nx,WITH->f2x);
  splie2(WITH->tab2-1,WITH->tab1-1,WITH->tz,WITH->nz,WITH->nx,WITH->f2z);
}
