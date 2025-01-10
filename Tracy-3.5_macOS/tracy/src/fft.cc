/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

   fft.c -- FFT routine

*/


void FFT(long int n, double *xr, double *xi)
{
  /* DFT using FFT algorithm */
  long m, i, j, j1, j2, jj1, jj2, jr, jf, k, l, mh, ifact1, ifact2, nfact1,
       nfact2, nfact3, lcoef1, lcoef2, ncoef1, ncoef2;
  double tcoef, theat, scoef, ccoef, s, c, r1, r2, a, b, sb, cb, TEMP;

  m = (int)(log((double)n) / log(2.0));
  for (l = 1; l <= m; l++) {
    lcoef1 = (int)floor(pow(2.0, l - 1.0) + 0.5);
    lcoef2 = lcoef1 * 2;
    ncoef1 = n / lcoef2;
    ncoef2 = ncoef1 * 2;
    tcoef = 2 * M_PI / n;
    theat = lcoef1 * tcoef;
    scoef = sin(theat);
    TEMP = sin(theat * 0.5);
    ccoef = -2 * (TEMP * TEMP);
    for (k = 0; k < lcoef1; k++) {
      j1 = k * ncoef2;
      j2 = j1 + ncoef1;
      s = 0.0;
      c = 1.0;
      for (j = 1; j <= ncoef1; j++) {
	jj1 = j + j1;
	jj2 = j + j2;
	r1 = xr[jj1 - 1];
	r2 = xi[jj1 - 1];
	xr[jj1 - 1] = r1 + xr[jj2 - 1];
	xi[jj1 - 1] = r2 + xi[jj2 - 1];
	a = r1 - xr[jj2 - 1];
	b = r2 - xi[jj2 - 1];
	xr[jj2 - 1] = a * c + b * s;
	xi[jj2 - 1] = b * c - a * s;
	sb = ccoef * s + scoef * c + s;
	cb = ccoef * c - scoef * s + c;
	s = sb;
	c = cb;
      }
    }
  }
  mh = m / 2;
  for (l = 1; l <= mh; l++) {
    ifact1 = (int)floor(pow(2.0, l - 1.0) + 0.5);
    ifact2 = ifact1 * 4;
    nfact1 = n / ifact2;
    nfact2 = nfact1 * 2;
    nfact3 = nfact2 * 2;
    for (k = 0; k < ifact1; k++) {
      jr = k * nfact3;
      j1 = jr + ifact1;
      j2 = jr + nfact2;
      for (j = 1; j <= nfact1; j++) {
	jf = (j - 1) / ifact1;
	jf = jf * ifact1 + j;
	jj1 = jf + j1;
	jj2 = jf + j2;
	r1 = xr[jj1 - 1];
	r2 = xi[jj1 - 1];
	xr[jj1 - 1] = xr[jj2 - 1];
	xi[jj1 - 1] = xi[jj2 - 1];
	xr[jj2 - 1] = r1;
	xi[jj2 - 1] = r2;
      }
    }
  }
  for (i = 0; i < n; i++) {
    xr[i] = 2.0 / n * xr[i];
    xi[i] = 2.0 / n * xi[i];
  }
}




/* End. */
