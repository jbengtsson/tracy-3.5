/* Tracy-2

   J. Bengtsson  CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

  ety.c --  given a real general matrix, this subroutine
            reduces a submatrix situated in rows and columns
            low through high to upper hessenberg form by
            orthogonal similarity transformations.

*/


/****************************************************************************/
/* double dsign(double x, double y)

   Purpose:
       Returns x*sign(y)

   input:
       x double
       y double

   output:
       none

   return:
       x*sign(y)

   global variables:
       none

   specific functions:
       none

   comments
       none

****************************************************************************/
double dsign(double x, double y)
{
  x = fabs(x);
  if (y >= 0.0)
    return x;
  else
    return (-x);
}

/* Implementation */

/****************************************************************************/
/* void ETY(int n, int low, int high, Matrix &a, Vector &ort)

   Purpose:
        this subroutine is a translation of the algol procedure orthes,
        num. math. 12, 349-368(1968] by Martin and Wilkinson.
        handbook for auto. comp., vol.ii-linear algebra, 339-358(1971].

        given a real general matrix, this subroutine
        reduces a submatrix situated in rows and columns
        low through high to upper hessenberg form by
        orthogonal similarity transformations.

        on input-

           n is the order of the matrix,

           low and high are integers determined by the balancing
             subroutine  balanc.  ifbalanc  has not been used,
             set low:=1, high:=n,

           a contains the input matrix.

        on output-

           a contains the hessenberg matrix.  information about
             the orthogonal transformations used in the reduction
             is stored in the remaining triangle under the
             hessenberg matrix,

           ort contains further information about the transformations.
             only elements low through high are used.

        fortran routine by b. s. garbow
        modified by filippo neri. 


   input:
       none

   output:
       none

   return:
       none

   global variables:
       none

   specific functions:
       none

   comments
       none

****************************************************************************/
void ETY(int n, int low, int high, Matrix &a, Vector &ort)
{
  int    i, j, m, ii, jj, la, mp, kp1;
  double  f, g, h, scale;

  la = high - 1;
  kp1 = low + 1;
  /* if la < kp1 goto 200*/
  if (la < kp1)   /*180*/
    return;
  for (m = kp1; m <= la; m++) {
    h = 0.0;
    ort[m - 1] = 0.0;
    scale = 0.0;   /*90*/

    /*     ********** scale column (algol tol then not needed] ***********/

    for (i = m - 1; i < high; i++)
      scale += fabs(a[i][m - 2]);

    /* if scale = 0.0 goto 180 */

    if (scale != 0) {
      mp = m + high;   /*100*/
      /*     ********** for i:=igh step -1 until m for -- ***********/

      for (ii = m; ii <= high; ii++) {
	i = mp - ii;
	ort[i - 1] = a[i - 1][m - 2] / scale;
	h += ort[i - 1] * ort[i - 1];
	/*100*/
      }

      g = -dsign(sqrt(h), ort[m - 1]);
      h -= ort[m - 1] * g;
      ort[m - 1] -= g;   /*130*/
      /*     ********** form (i-(u*ut]/h]*a ***********/
      for (j = m - 1; j < n; j++) {   /*160*/
	f = 0.0;   /*110*/
	/*     ********** for i:=igh step -1 until m for -- ***********/
	for (ii = m; ii <= high; ii++) {
	  i = mp - ii;
	  f += ort[i - 1] * a[i - 1][j];
	  /*110*/
	}

	f /= h;   /*120*/

	for (i = m - 1; i < high; i++) {
	  /*120*/
	  a[i][j] -= f * ort[i];
	}

	/*130*/
      }

      /*     ********** form (i-(u*ut]/h]*a*(i-(u*ut]/h] ***********/
      for (i = 0; i < high; i++) {
	f = 0.0;   /*140*/
	/*     ********** for j:=igh step -1 until m for -- ***********/
	for (jj = m; jj <= high; jj++) {
	  j = mp - jj;
	  f += ort[j - 1] * a[i][j - 1];
	  /*140*/
	}

	f /= h;

	for (j = m - 1; j < high; j++)
	  a[i][j] -= f * ort[j];

	/*160*/
      }

      ort[m - 1] = scale * ort[m - 1];
      a[m - 1][m - 2] = scale * g;
    }
    /*180*/
  }
}


/****************************************************************************/
/* void ETYT(int n, int low, int high, Matrix &a, Vector &ort, Matrix &z)

   Purpose:
        this subroutine is a translation of the algol procedure ortrans,
        num. math. 16, 181-204[1970] by peters and wilkinson.
        handbook for auto. comp., vol.ii-linear algebra, 372-395[1971].

        this subroutine accumulates the orthogonal similarity
        transformations used in the reduction of a real general
        matrix to upper hessenberg form by  ety.

        on input-

           n is the order of the matrix,

           low and high are integers determined by the balancing
             subroutine  balanc.  if  balanc  has not been used,
             set low:=1, high:=n,

           a contains information about the orthogonal trans-
             formations used in the reduction by  orthes
             in its strict lower triangle,

             ort contains further information about the trans-
             formations used in the reduction by  ety.
             only elements low through high are used.

        on output-

           z contains the transformation matrix produced in the
             reduction by  ety,

           ort has been altered.

        fortran routine by b. s. garbow.
        modified by f. neri.


        ********** initialize z to identity matrix ***********


   input:
       none

   output:
       none

   return:
       none

   global variables:
       none

   specific functions:
       none

   comments
       none

****************************************************************************/
void ETYT(int n, int low, int high, Matrix &a, Vector &ort, Matrix &z)
{
  int i, j, kl, mm, mp, mp1;
  double g;


  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++)
      z[i][j] = 0.0;
    z[i][i] = 1.0;
  }

  kl = high - low - 1;
  if (kl < 1)
    return;
  for (mm = 1; mm <= kl; mm++) {
    mp = high - mm;
    if (a[mp - 1][mp - 2] != 0.0) {
      mp1 = mp + 1;
      for (i = mp1 - 1; i < high; i++)
	ort[i] = a[i][mp - 2];
      for (j = mp - 1; j < high; j++) {
	g = 0.0;
	for (i = mp - 1; i < high; i++)
	  g += ort[i] * z[i][j];
	/*     ********** divisor below is negative of h formed in orthes.
	                 double division avoids possible underflow ***********/
	g = g / ort[mp - 1] / a[mp - 1][mp - 2];
	for (i = mp - 1; i < high; i++)
	  z[i][j] += g * ort[i];
      }
    }
  }
}


#define machep          1e-20

/****************************************************************************/
/* int min0(int i, int j)

   Purpose:
       minimum of integers i and j

   input:
       none

   output:
       none

   return:
       minimum of integers i and j

   global variables:
       none

   specific functions:
       none

   comments
       none

****************************************************************************/
int min0(int i, int j)
{
  if (i >= j)
    return j;
  else
    return i;
}

/****************************************************************************/
/* void etdiv(Vector &a, Vector &b, double c, double d,
                 double e, double f)

   Purpose:
       computes the complex division
          a + ib := (c + id)/(e + if)
       very slow, but tries to be as accurate as
       possible by changing the order of the
       operations, so to avoid under(over)flow
       problems.
       Written by F. Neri Feb. 12 1986


   input:
       none

   output:
       none

   return:
       none

   global variables:
       none

   specific functions:
       none

   comments
       none

****************************************************************************/
void etdiv(double *a, double *b, double c, double d, double e, double f)
{
  double s, t, cc, dd, ee, ff, temp;
  int flip;

  flip = 0;
  cc = c;
  dd = d;
  ee = e;
  ff = f;
  if (fabs(f) >= fabs(e))
  {
    ee = f;
    ff = e;
    cc = d;
    dd = c;
    flip = 1;
  }
  s = 1.0 / ee;
  t = 1.0 / (ee + ff * ff * s);
  if (fabs(ff) >= fabs(s))
  {
    temp = ff;
    ff = s;
    s = temp;
  }

  if (fabs(dd) >= fabs(s))
    *a = t * (cc + s * dd * ff);
  else if (fabs(dd) >= fabs(ff))
    *a = t * (cc + dd * s * ff);
  else
    *a = t * (cc + ff * s * dd);

  if (fabs(cc) >= fabs(s))
    *b = t * (dd - s * cc * ff);
  else if (fabs(cc) >= fabs(ff))
    *b = t * (dd - cc * s * ff);
  else
    *b = t * (dd - ff * s * cc);

  if (flip != 0)
    *b = -*b;

}


/****************************************************************************/
/* void ety2(int n, int low, int high, Matrix &h, Vector &wr,
          Vector &wi, Matrix &z, int *ierr)

   Purpose:
        this subroutine is a translation of the algol procedure hqr2,
        num. math. 16, 181-204[1970] by peters and wilkinson.
        handbookfor{auto. comp., vol.ii-linear algebra, 372-395[1971].

        this subroutine finds the eigenvalues and eigenvectors
        of a real upper hessenberg matrix by the qr method.  the
        eigenvectors of a real general matrix can also be found
        ifelmhes  and  eltran  or  orthes  and  ortran  have
        been used to reduce this general matrix to hessenberg form
        and to accumulate the similarity transformations.

        on input-

           n is the order of the matrix,

           low and high are integers determined by the balancing
             subroutine  balanc.  ifbalanc  has not been used,
             set low:=1, high:=n,

           h contains the upper hessenberg matrix,

           z contains the transformation matrix produced by  eltran
             after the reduction by  elmhes, or by  ortran  after the
             reduction by  orthes, ifperformed.  ifthe eigenvectors
             of the hessenberg matrix are desired, z must contain the
             identity matrix.

        on output-

           h has been destroyed,

           wr and wi contain the real and imaginary parts,
             respectively, of the eigenvalues.  the eigenvalues
             are unordered except that complex conjugate pairs
             of values appear consecutively with the eigenvalue
             having the positive imaginary part first.  ifan
             error exit is made, the eigenvalues should be correct
            for{indices ierr+1,...,n,

           z contains the real and imaginary parts of the eigenvectors.
             ifthe i-th eigenvalue is real, the i-th column of z
             contains its eigenvector.  ifthe i-th eigenvalue is complex
             with positive imaginary part, the i-th and [i+1]-th
             columns of z contain the real and imaginary parts of its
             eigenvector.  the eigenvectors are unnormalized.  ifan
             error exit is made, none of the eigenvectors has been found,

           ierr is set to
             zero      for{normal return,
             j          ifthe j-th eigenvalue has not been
                        determined after 30 iterations.

        arithmetic is double precision. complex division
        is simulated by routin etdiv.

        fortran routine by b. s. garbow.
        modified by f. neri.


                   machep is a machine dependent parameter specifying
                   the relative precision of floating polong arithmetic.



   input:
       none

   output:
       none

   return:
       none

   global variables:
       none

   specific functions:
       none

   comments
       none

****************************************************************************/
void ety2(int n, int low, int high, Matrix &h, Vector &wr,
          Vector &wi, Matrix &z, int &ierr)
{
  int i, j, k, l = 0, m = 0, en, ii, jj, ll, mm, na, nn, its, mp2, enm2;
  double p = 0.0, q = 0.0, r = 0.0, s = 0.0, t, w, x, y, ra, sa, vi, vr, zz = 0.0, norm;
  bool notlas;
  double z3r, z3i;
  int FORLIM;

  ierr = 0;
  norm = 0.0;
  k = 1;   /*50*/
  /*     ********** store roots isolated by balanc
                   and compute matrix norm ********** */
  for (i = 0; i < n; i++) {   /*40*/
    for (j = k - 1; j < n; j++) {
      /*40*/
      norm += fabs(h[i][j]);
    }
    k = i + 1;
    /*   if (i >= low) and (i <= high) goto 50*/
    if (i + 1 < low || i + 1 > high) {
      wr[i] = h[i][i];
      wi[i] = 0.0;
    }
    /*50*/
  }

  en = high;
  t = 0.0;

  /*     ********** searchfor{next eigenvalues ***********/
_L60:
  if (en < low)
    goto _L340;
  its = 0;
  na = en - 1;
  enm2 = na - 1;
  /*     ********** lookfor{single small sub-diagonal element*/

_L70:   /*80*/
  for (ll = low; ll <= en; ll++) {
    l = en + low - ll;
    if (l == low)
      goto _L100;
    s = fabs(h[l - 2][l - 2]) + fabs(h[l - 1][l - 1]);
    if (s == 0.0)
      s = norm;
    if (fabs(h[l - 1][l - 2]) <= machep * s)
      goto _L100;
    /*80*/
  }


  /*     ********** form shift ***********/

_L100:
  x = h[en - 1][en - 1];
  if (l == en)
    goto _L270;
  y = h[na - 1][na - 1];
  w = h[en - 1][na - 1] * h[na - 1][en - 1];
  if (l == na)
    goto _L280;
  if (its == 30) {
    ierr = en;
    goto _L1000;
  }

  /*      if (its <> 10) and (its <> 20)  then goto 130;*/
  if (its == 10 || its == 20) {   /*140*/
    /*     ********** form exceptional shift ***********/
    t += x;   /*120*/

    for (i = low - 1; i < en; i++) {
      /*120*/
      h[i][i] -= x;
    }

    s = fabs(h[en - 1][na - 1]) + fabs(h[na - 1][enm2 - 1]);
    x = 0.75 * s;
    y = x;
    w = -0.4375 * s * s;
    /*130:*/
    its++;
  }

  /*     ********** lookfor{two consecutive small
                   sub-diagonal elements.*/

  for (mm = l; mm <= enm2; mm++) {
    m = enm2 + l - mm;
    zz = h[m - 1][m - 1];
    r = x - zz;
    s = y - zz;
    p = (r * s - w) / h[m][m - 1] + h[m - 1][m];
    q = h[m][m] - zz - r - s;
    r = h[m + 1][m];
    s = fabs(p) + fabs(q) + fabs(r);
    p /= s;
    q /= s;
    r /= s;
    if (m == l)
      goto _L150;
    if (fabs(h[m - 1][m - 2]) * (fabs(q) + fabs(r)) <=
	machep * fabs(p) * (fabs(h[m - 2][m - 2]) + fabs(zz) + fabs(h[m][m])))
      goto _L150;
    /*140*/
  }

_L150:
  mp2 = m + 2;   /*160*/

  for (i = mp2; i <= en; i++) {   /*260*/
    h[i - 1][i - 3] = 0.0;
    /* if i = mp2 goto 160;*/
    if (i != mp2)
      h[i - 1][i - 4] = 0.0;
    /*160*/
  }

  /*     ********** double qr step involving rows l to en and
                   columns m to en ***********/
  FORLIM = na;
  for (k = m; k <= FORLIM; k++) {  /*a*/
    notlas = (k != na);
    /* if k = m goto 170*/
    if (k != m) {  /*b*/
      p = h[k - 1][k - 2];
      q = h[k][k - 2];
      r = 0.0;
      if (notlas)
	r = h[k + 1][k - 2];
      x = fabs(p) + fabs(q) + fabs(r);
      if (x == 0.0)
	goto _L260;
      p /= x;
      q /= x;
      r /= x;
    }  /*b*/

    s = dsign(sqrt(p * p + q * q + r * r), p);
    if (k != m)
      h[k - 1][k - 2] = -s * x;
    else if (l != m)
      h[k - 1][k - 2] = -h[k - 1][k - 2];

    p += s;
    x = p / s;
    y = q / s;
    zz = r / s;
    q /= p;
    r /= p;   /*210*/

    /*     ********** row modification ***********/

    for (j = k - 1; j < n; j++) {  /*d*/
      p = h[k - 1][j] + q * h[k][j];
      if (notlas) {  /*e*/
	p += r * h[k + 1][j];
	h[k + 1][j] -= p * zz;
      }  /*e*/
      h[k][j] -= p * y;
      h[k - 1][j] -= p * x;
      /*210*/
    }  /*d*/

    j = min0(en, k + 3);   /*230*/

    /*     ********** column modification ***********/

    for (i = 0; i < j; i++)   /*250*/
    {  /*d*/
      p = x * h[i][k - 1] + y * h[i][k];
      if (notlas) {  /*e*/
	p += zz * h[i][k + 1];
	h[i][k + 1] -= p * r;
      }  /*e*/
      h[i][k] -= p * q;
      h[i][k - 1] -= p;
    }  /*d*/

    /*     ********** accumulate transformations ***********/

    for (i = low - 1; i < high; i++) {  /*d*/
      p = x * z[i][k - 1] + y * z[i][k];
      if (notlas) {  /*e*/
	p += zz * z[i][k + 1];
	z[i][k + 1] -= p * r;
      }  /*e*/
      z[i][k] -= p * q;
      z[i][k - 1] -= p;
      /*250*/
    }  /*d*/
_L260: ;
  }  /*a*/

  goto _L70;
  /*     ********** one root found ***********/
_L270:
  h[en - 1][en - 1] = x + t;
  wr[en - 1] = h[en - 1][en - 1];
  wi[en - 1] = 0.0;
  en = na;
  goto _L60;
  /*     ********** two roots found ***********/
_L280:
  p = (y - x) / 2.0;
  q = p * p + w;
  zz = sqrt(fabs(q));
  h[en - 1][en - 1] = x + t;
  x = h[en - 1][en - 1];
  h[na - 1][na - 1] = y + t;
  if (q < 0.0)
    goto _L320;
  /*     ********** real pair ***********/
  zz = p + dsign(zz, p);
  wr[na - 1] = x + zz;
  wr[en - 1] = wr[na - 1];
  if (zz != 0.0)
    wr[en - 1] = x - w / zz;
  wi[na - 1] = 0.0;
  wi[en - 1] = 0.0;
  x = h[en - 1][na - 1];
  s = fabs(x) + fabs(zz);
  p = x / s;
  q = zz / s;
  r = sqrt(p * p + q * q);
  p /= r;
  q /= r;   /*290*/
  /*     ********** row modification ***********/
  for (j = na - 1; j < n; j++) {   /*300*/
    zz = h[na - 1][j];
    h[na - 1][j] = q * zz + p * h[en - 1][j];
    h[en - 1][j] = q * h[en - 1][j] - p * zz;
    /*290*/
  }
  /*     ********** column modification ***********/
  for (i = 0; i < en; i++) {   /*310*/
    zz = h[i][na - 1];
    h[i][na - 1] = q * zz + p * h[i][en - 1];
    h[i][en - 1] = q * h[i][en - 1] - p * zz;
    /*300*/
  }
  /*     ********** accumulate transformations ***********/
  for (i = low - 1; i < high; i++) {
    zz = z[i][na - 1];
    z[i][na - 1] = q * zz + p * z[i][en - 1];
    z[i][en - 1] = q * z[i][en - 1] - p * zz;
    /*310*/
  }

  goto _L330;
  /*     ********** complex pair ***********/
_L320:
  wr[na - 1] = x + p;
  wr[en - 1] = x + p;
  wi[na - 1] = zz;
  wi[en - 1] = -zz;
_L330:
  en = enm2;
  goto _L60;
  /*     ********** all roots found.  backsubstitute to find    */
  /*                vectors of upper triangular form ********** */
  /*  340 if norm = 0.0 goto 1001*/
_L340:
  if (norm != 0.0) {  /*0.5*/
    for (nn = 1; nn <= n; nn++) {  /*1*/
      en = n - nn + 1;
      p = wr[en - 1];
      q = wi[en - 1];
      na = en - 1;
      if (q == 0) {   /*2*/
	m = en;
	h[en - 1][en - 1] = 1.0;
	if (na != 0) {  /*3*/
	  for (ii = 1; ii <= na; ii++)   /*4*/
	  {  /*4*/
	    i = en - ii;

	    w = h[i - 1][i - 1] - p;
	    r = h[i - 1][en - 1];

	    if (m <= na) {
	      for (j = m - 1; j < na; j++)
		r += h[i - 1][j] * h[j][en - 1];
	    }

	    if (wi[i - 1] < 0.0)   /*5*/
	    {  /*5*/
	      zz = w;
	      s = r;
	    } else   /*5*/
	    {  /*5*/
	      m = i;
	      if (wi[i - 1] == 0.0)   /*6*/
	      {  /*6*/
		t = w;
		if (w == 0.0)
		  t = machep * norm;
		h[i - 1][en - 1] = -(r / t);
	      } else   /*6*/
	      {  /*6*/
		/* ********** solve double equations ***********/
		x = h[i - 1][i];
		y = h[i][i - 1];
		q = (wr[i - 1] - p) * (wr[i - 1] - p) + wi[i - 1] * wi[i - 1];
		t = (x * s - zz * r) / q;
		h[i - 1][en - 1] = t;
		if (fabs(x) <= fabs(zz))
		  h[i][en - 1] = (-s - y * t) / zz;
		else
		  h[i][en - 1] = (-r - w * t) / x;
	      }
	    }
	  }
	}  /*3*/
      } else {
	if (q < 0)   /* ** Complex ** */
	{  /*2*/
	  m = na;
	  /*  last vector component chosen imaginary so that
	      eigenvector matrix is triangular */
	  /* if Abs(h[en,na))<=Abs(h[na,en)) then  720 */
	  if (na != 0) {  /*2.5*/
	    if (fabs(h[en - 1][na - 1]) > fabs(h[na - 1][en - 1]))   /*3*/
	    {  /*3*/
	      h[na - 1][na - 1] = q / h[en - 1][na - 1];
	      h[na - 1][en - 1] = (p - h[en - 1][en - 1]) / h[en - 1][na - 1];
	    } else {  /*3*/
	      etdiv(&z3r, &z3i, 0.0, -h[na - 1]
		    [en - 1], h[na - 1][na - 1] - p, q);
	      h[na - 1][na - 1] = z3r;
	      h[na - 1][en - 1] = z3i;
	    }  /*3*/
	    h[en - 1][na - 1] = 0.0;
	    h[en - 1][en - 1] = 1.0;
	    enm2 = na - 1;

	    /*if enm2 =0 then  800*/
	    if (enm2 != 0)   /*3*/
	    {  /*3*/
	      for (ii = 1; ii <= enm2; ii++) {  /*4*/
		i = na - ii;
		w = h[i - 1][i - 1] - p;
		ra = 0.0;
		sa = h[i - 1][en - 1];

		for (j = m - 1; j < na; j++) {  /*5*/
		  ra += h[i - 1][j] * h[j][na - 1];
		  sa += h[i - 1][j] * h[j][en - 1];
		}  /*5*/

		if (wi[i - 1] < 0.0)   /*5*/
		{  /*5*/
		  zz = w;
		  r = ra;
		  s = sa;
		} else {  /*5*/
		  m = i;
		  if (wi[i - 1] == 0.0) {  /*6*/
		    etdiv(&z3r, &z3i, -ra, -sa, w, q);
		    h[i - 1][na - 1] = z3r;
		    h[i - 1][en - 1] = z3i;
		  }  /*6*/
		  else {  /*6*/
		    x = h[i - 1][i];
		    y = h[i][i - 1];
		    vr = (wr[i - 1] - p) * (wr[i - 1] - p) +
			 wi[i - 1] * wi[i - 1] - q * q;
		    vi = (wr[i - 1] - p) * 2.0 * q;

		    if (vr == 0.0 && vi == 0.0)
		      vr = machep * norm *
			   (fabs(w) + fabs(q) + fabs(x) + fabs(y) + fabs(zz));

		    etdiv(&z3r, &z3i, x * r - zz * ra + q * sa,
			  x * s - zz * sa - q * ra, vr, vi);

		    h[i - 1][na - 1] = z3r;
		    h[i - 1][en - 1] = z3i;

		    if (fabs(x) <= fabs(zz) + fabs(q))   /*7*/
		    {  /*7*/
		      etdiv(&z3r, &z3i, -r - y * h[i - 1][na - 1],
			    -s - y * h[i - 1][en - 1], zz, q);
		      h[i][na - 1] = z3r;
		      h[i][en - 1] = z3i;
		    } else {  /*7*/
		      h[i]
			[na - 1] = (q * h[i - 1][en - 1] - w * h[i - 1]
							   [na - 1] - ra) / x;
		      h[i]
			[en - 1] = (-sa - w * h[i - 1][en - 1] - q * h[i - 1]
				      [na - 1]) / x;
		    }  /*7*/
		  }  /*6*/
		}  /*4*/
	      }  /*5*/
	    }
	  }  /*2.5*/
	}  /*2*/
      }
      /*** double vector ***/
      /*2: of if q=0 then */
    }  /*1*/


    /**     ********** end back substitution.               **/
    /**               vectors of isolated roots **********  **/
    for (i = 0; i < n; i++) {
      if (i + 1 < low || i + 1 > high) {
	for (j = i; j < n; j++)
	  z[i][j] = h[i][j];
      }
    }

    /**     ********** multiply by transformation matrix to give  **/
    /**                vectors of original full matrix.           **/
    /**                for j:=n step -1 until low for --          **/


    for (jj = low; jj <= n; jj++) {  /*1*/
      j = n + low - jj;
      m = min0(j, high);
      for (i = low - 1; i < high; i++) {  /*2*/
	zz = 0.0;
	for (k = low - 1; k < m; k++)
	  zz += z[i][k] * h[k][j - 1];
	z[i][j - 1] = zz;
      }  /*2*/
    }  /*1*/
  }  /*0.5*/
  en = 0;
  /*     ********** last card of ety2 ***********/
_L1000:
  ierr = en;

  /*  780     ********** solve complex equations ********** */
}

#undef machep

/* End. */
