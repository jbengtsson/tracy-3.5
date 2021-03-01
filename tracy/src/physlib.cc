/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

*/


double FindRes_eps = 0.5e-6;


/**************************/
/* Routines for printing  */
/**************************/

/** Get time and date **/
void prt_sigma(LatticeType &lat)
{
  long int  i;
  double    code = 0.0;
  MpoleType *M;
  FILE      *outf;

  outf = file_write("../out/sigma.out");

  fprintf(outf, "#  name     s   sqrt(sx)   sqrt(sx')  sqrt(sy)  sqrt(sy')\n");
  fprintf(outf, "#           [m]   [mm]       [mrad]     [mm]     [mrad]\n");
  fprintf(outf, "#\n");

  for (i = 0; i <= lat.conf.Cell_nLoc; i++) {
    switch (lat.elems[i]->Pkind) {
    case drift:
      code = 0.0;
      break;
    case Mpole:
      M = dynamic_cast<MpoleType*>(lat.elems[i]);
      if (M->Pirho != 0)
	code = 0.5;
      else if (M->PBpar[Quad+HOMmax] != 0)
	code = sgn(M->PBpar[Quad+HOMmax]);
      else if (M->PBpar[Sext+HOMmax] != 0)
	code = 1.5*sgn(M->PBpar[Sext+HOMmax]);
      else if (lat.elems[i]->Fnum == lat.conf.bpm)
	code = 2.0;
      else
	code = 0.0;
      break;
    default:
      code = 0.0;
      break;
    }
    fprintf(outf, "%4ld %.*s %6.2f %4.1f %9.3e %9.3e %9.3e %9.3e\n",
            i, SymbolLength, lat.elems[i]->PName, lat.elems[i]->S, code,
            1e3*sqrt(lat.elems[i]->sigma[x_][x_]),
	    1e3*sqrt(fabs(lat.elems[i]->sigma[x_][px_])),
	    1e3*sqrt(lat.elems[i]->sigma[y_][y_]),
	    1e3*sqrt(fabs(lat.elems[i]->sigma[y_][py_])));
  }

  fclose(outf);
}


void recalc_S(LatticeType &lat)
{
  long int  k;
  double    S_tot;

  S_tot = 0.0;
  for (k = 0; k <= lat.conf.Cell_nLoc; k++) {
    S_tot += lat.elems[k]->PL; lat.elems[k]->S = S_tot;
  }
}


void getabn(LatticeType &lat, Vector2 &alpha, Vector2 &beta, Vector2 &nu)
{
  Vector2 gamma;

  Cell_GetABGN(lat.conf.OneTurnMat, alpha, beta, gamma, nu);
}


#define nfloq 4

void getfloqs(LatticeType &lat, psVector &x)
{
  // Transform to Floquet space
  LinTrans(nfloq, lat.conf.Ascrinv, x);
}

#undef nfloq


#define ntrack 4

// 4D tracking in normal or Floquet space over nmax turns

void track(const char *file_name, LatticeType &lat, double ic1, double ic2,
	   double ic3, double ic4, double dp, long int nmax, long int &lastn,
	   long int &lastpos, int floqs, double f_rf)
{
  /* Single particle tracking around closed orbit:

          Output                floqs

        Phase Space               0     [x, px, y, py, delta, ct]  
	Floquet Space             1     [x^, px^, y^, py^, delta, ct]
	Action-Angle Variables    2     [2Jx, phx, 2Jy, phiy, delta, ct]

  */
  long int i;
  double   twoJx, twoJy, phix, phiy, scl_1 = 1.0, scl_2 = 1.0;
  psVector x0, x1, x2, xf;
  FILE     *outf;

  bool  prt = false;

  if (floqs == 0) {
    scl_1 = 1e3; scl_2 = 1e3;
    x0[x_] = ic1; x0[px_] = ic2; x0[y_] = ic3; x0[py_] = ic4;
  } else if (floqs == 1) {
    scl_1 = 1.0; scl_2 = 1.0;
    x0[x_] = ic1; x0[px_] = ic2; x0[y_] = ic3; x0[py_] = ic4;
    LinTrans(4, lat.conf.Ascr, x0);
  } else if (floqs == 2) {
    scl_1 = 1e6; scl_2 = 1.0;
    x0[x_] = sqrt(ic1)*cos(ic2); x0[px_] = -sqrt(ic1)*sin(ic2);
    x0[y_] = sqrt(ic3)*cos(ic4); x0[py_] = -sqrt(ic3)*sin(ic4);
    LinTrans(4, lat.conf.Ascr, x0);
  }

  outf = file_write(file_name);
  fprintf(outf, "# Tracking with TRACY");
  lat.getcod(dp, lastpos);
  if (floqs == 0)
    fprintf(outf, "\n");
  else if (floqs == 1) {
    lat.Ring_GetTwiss(false, dp);
    fprintf(outf, "# (Floquet space)\n");
  } else if (floqs == 2) {
    lat.Ring_GetTwiss(false, dp);
    fprintf(outf, "# (Action-Angle variables)\n");
  }
  fprintf(outf, "#\n");
  fprintf(outf, "#%3d%6ld% .1E% .1E% .1E% .1E% 7.5f% 7.5f\n",
	  1, nmax, 1e0, 1e0, 0e0, 0e0, lat.conf.TotalTune[0],
	  lat.conf.TotalTune[1]);
  if (floqs == 0) {
    fprintf(outf, "#    N       x            p_x            y            p_y");
    fprintf(outf, "          delta          cdt\n");
    fprintf(outf, "#           [mm]         [mrad]"
	    "         [mm]         [mrad]");
  } else if  (floqs == 1) {
    fprintf(outf, "#    N     x^          px^          y^          py^");
    fprintf(outf, "          delta          cdt\n");
    fprintf(outf, "#                              "
	    "                            ");
  } else if  (floqs == 2) {
    fprintf(outf, "#    N     2Jx          phi_x          2Jy          phi_y");
    fprintf(outf, "          delta          cdt\n");
    fprintf(outf, "#                              "
	    "                            ");
  }
  if (f_rf == 0.0){
    fprintf(outf, "         [%%]           [mm]\n");
    fprintf(outf, "#\n");
    fprintf(outf, "%4d %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e\n",
	    0, scl_1*ic1, scl_2*ic2, scl_1*ic3, scl_2*ic4, 1e2*dp, 
	    0*1e3*lat.conf.CODvect[ct_]);
  } else {
    fprintf(outf, "         [%%]           [deg]\n");
    fprintf(outf, "#\n");
    fprintf(outf, "%4d %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e\n",
	    0, scl_1*ic1, scl_2*ic2, scl_1*ic3, scl_2*ic4, 1e2*dp, 
	    2.0*f_rf*180.0*lat.conf.CODvect[ct_]/c0);
  }
  x2[x_] = x0[x_] + lat.conf.CODvect[x_];
  x2[px_] = x0[px_] + lat.conf.CODvect[px_];
  x2[y_] = x0[y_] + lat.conf.CODvect[y_];
  x2[py_] = x0[py_] + lat.conf.CODvect[py_];
  if (lat.conf.Cavity_on) {
    x2[delta_] = dp + lat.conf.CODvect[delta_]; x2[ct_] = lat.conf.CODvect[ct_];
  } else {
    x2[delta_] = dp; x2[ct_] = 0.0;
  }

  lastn = 0;

  if (prt) {
    printf("\n");
    printf("track:\n");
    printf("%4ld %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
	   lastn, 1e3*x2[x_], 1e3*x2[px_], 1e3*x2[y_], 1e3*x2[py_],
	   1e2*x2[delta_], 1e3*x2[ct_]);
  }
    
  do {
    (lastn)++;
    for (i = 0; i < nv_; i++)
      x1[i] = x2[i];

    lat.Cell_Pass(0, lat.conf.Cell_nLoc, x2, lastpos);

    for (i = x_; i <= py_; i++)
      xf[i] = x2[i] - lat.conf.CODvect[i];

    for (i = delta_; i <= ct_; i++)
      if (lat.conf.Cavity_on && (i != ct_)) 
	xf[i] = x2[i] - lat.conf.CODvect[i];
      else
	xf[i] = x2[i];

    if (floqs == 1)
      getfloqs(lat, xf);
    else if (floqs == 2) {
      getfloqs(lat, xf);
      twoJx = pow(xf[x_], 2.0) + pow(xf[px_], 2.0);
      twoJy = pow(xf[y_], 2.0) + pow(xf[py_], 2.0);
      phix = atan2(xf[px_], xf[x_]);
      phiy = atan2(xf[py_], xf[y_]);
      xf[x_] = twoJx; xf[px_] = phix; xf[y_] = twoJy; xf[py_] = phiy;
    }
    if (f_rf == 0.0)
      fprintf(outf,
	      "%4ld %23.16le %23.16le %23.16le %23.16le %23.16le %23.16le\n",
	      lastn, scl_1*xf[0], scl_2*xf[1], scl_1*xf[2], scl_2*xf[3],
	      1e2*xf[4], 1e3*(xf[5]-lat.conf.CODvect[ct_]));
    else
      fprintf(outf,
	      "%4ld %23.16le %23.16le %23.16le %23.16le %23.16le %23.16le\n",
	      lastn, scl_1*xf[0], scl_2*xf[1], scl_1*xf[2], scl_2*xf[3],
	      1e2*xf[4], 2.0*f_rf*180.0*xf[5]/c0);
  } while ((lastn != nmax) && (lastpos == lat.conf.Cell_nLoc));

  fclose(outf);
}

#undef ntrack


#define step            0.1
#define px              0.0
#define py              0.0

void track_(LatticeType &lat, double r, struct LOC_getdynap *LINK)
{
  long i, lastn, lastpos;
  psVector x;

  x[0] = r * cos(LINK->phi);
  x[1] = px;
  x[2] = r * sin(LINK->phi);
  x[3] = py;
  x[4] = LINK->delta;
  x[5] = 0.0;
  /* transform to phase space */
  if (LINK->floqs) {
    LinTrans(5, lat.conf.Ascr, x);
  }
  for (i = 0; i <= 3; i++)
    x[i] += lat.conf.CODvect[i];
  lastn = 0;
  do {
    lastn++;
    lat.Cell_Pass(0, lat.conf.Cell_nLoc, x, lastpos);
  } while (lastn != LINK->nturn && lastpos == lat.conf.Cell_nLoc);
  LINK->lost = (lastn != LINK->nturn);
}

#undef step
#undef px
#undef py



/****************************************************************************/
/* void getcsAscr(void)

   Purpose:
        Get Courant-Snyder Ascr


   Input:
       none

   Output:
       none

   Return:
       none

   Global variables:
       none

   Specific functions:
       none

   Comments:
       none

****************************************************************************/
void getcsAscr(LatticeType &lat)
{
  long i, j;
  double phi;
  Matrix R;

  UnitMat(6, R);
  for (i = 1; i <= 2; i++) {
    phi = -atan2(lat.conf.Ascr[i * 2 - 2][i * 2 - 1], lat.conf.Ascr[i * 2 - 2]
		[i * 2 - 2]);
    R[i * 2 - 2][i * 2 - 2] = cos(phi);
    R[i * 2 - 1][i * 2 - 1] = R[i * 2 - 2][i * 2 - 2];
    R[i * 2 - 2][i * 2 - 1] = sin(phi);
    R[i * 2 - 1][i * 2 - 2] = -R[i * 2 - 2][i * 2 - 1];
  }
  MulRMat(6, lat.conf.Ascr, R);
  for (i = 1; i <= 2; i++) {
    if (lat.conf.Ascr[i * 2 - 2][i * 2 - 2] < 0.0) {
      for (j = 0; j <= 5; j++) {
	lat.conf.Ascr[j][i * 2 - 2] = -lat.conf.Ascr[j][i * 2 - 2];
	lat.conf.Ascr[j][i * 2 - 1] = -lat.conf.Ascr[j][i * 2 - 1];
      }
    }
  }
  if (!InvMat(6, lat.conf.Ascrinv))
    printf("  *** Ascr is singular\n");
}


void GetTrack(const char *file_name,
	      long *n, double x[], double px[], double y[], double py[])
{
  int   k;
  char  line[200];
  FILE  *inf;

  inf = file_read(file_name);

  do {
    fgets(line, 200, inf);
  } while (strstr(line, "#") != NULL);

  // skip initial conditions
  fgets(line, 200, inf);

  do {
     sscanf(line, "%d", &k);
     sscanf(line, "%*d %lf %lf %lf %lf", &x[k-1], &px[k-1], &y[k-1], &py[k-1]);
  } while (fgets(line, 200, inf) != NULL);

  *n = k;

   fclose(inf);
}


/****************************************************************************/
/* void Getj(long n, double *x, double *px, double *y, double *py)

   Purpose:
        Calculates the linear invariant

   Input:
       none

   Output:
       none

   Return:
       none

   Global variables:
       none

   Specific functions:
       none

   Comments:
       none

****************************************************************************/
void Getj(long n, double *x, double *px, double *y, double *py)
{
  long int i;

  for (i = 0; i < n; i++) {
    x[i] = (pow(x[i], 2.0)+pow(px[i], 2.0))/2.0;
    y[i] = (pow(y[i], 2.0)+pow(py[i], 2.0))/2.0;
  }
}

/****************************************************************************/
/* double GetArg(double x, double px, double nu)

   Purpose:
       get argument of x

   Input:
       none

   Output:
       none

   Return:
       none

   Global variables:
       none

   Specific functions:
       none

   Comments:
       17/07/03 use M_PI instead of pi

****************************************************************************/
double GetArg(double x, double px, double nu)
{
  double phi, val;

  phi = GetAngle(x, px);
  if (phi < 0.0)
    phi += 2.0 * M_PI;
  val = phi + Fract(nu) * 2.0 * M_PI;
  if (val < 0.0)
    val += 2.0 * M_PI;
  return val;
}

/****************************************************************************/
/* void GetPhi(long n, double *x, double *px, double *y, double *py)

   Purpose:
       get linear phases

   Input:
       none

   Output:
       none

   Return:
       none

   Global variables:
       none

   Specific functions:
       none

   Comments:
       none

****************************************************************************/
void GetPhi(LatticeType &lat, long n, double *x, double *px, double *y,
	    double *py)
{
  /* Calculates the linear phase */
  long i;

  for (i = 1; i <= n; i++) {
    x[i - 1] = GetArg(x[i - 1], px[i - 1], i * lat.conf.TotalTune[0]);
    y[i - 1] = GetArg(y[i - 1], py[i - 1], i * lat.conf.TotalTune[1]);
  }
}


/*********************************/
/* Routines for Fourier analysis */
/*********************************/

void Sinfft(int n, double xr[])
{
  /* DFT with sine window */
  int     i;
  double  xi[n];

  for (i = 0; i < n; i++) {
    xr[i] = sin((double)i/(double)n*M_PI)*xr[i]; xi[i] = 0.0;
  }
  FFT(n, xr, xi);
  for (i = 0; i < n; i++)
    xr[i] = sqrt(xr[i]*xr[i]+xi[i]*xi[i]);
}

void sin_FFT(int n, double xr[])
{
  /* DFT with sine window */
  int     i;
  double  *xi;

  xi = dvector(1, 2*n);

  for (i = 1; i <= n; i++) {
    xi[2*i-1] = sin((double)i/n*M_PI)*xr[i-1]; xi[2*i] = 0.0;
  }
  dfour1(xi, (unsigned long)n, 1);
  for (i = 1; i <= n; i++)
    xr[i-1] = sqrt(pow(xi[2*i-1], 2.0)+pow(xi[2*i], 2.0))*2.0/n;

  free_dvector(xi, 1, 2*n);
}


void sin_FFT(int n, double xr[], double xi[])
{
  /* DFT with sine window */
  int      i;
  double  *xri;

  xri = dvector(1, 2*n);

  for (i = 1; i <= n; i++) {
    xri[2*i-1] = sin((double)i/n*M_PI)*xr[i-1];
    xri[2*i] = sin((double)i/n*M_PI)*xi[i-1];
  }
  dfour1(xri, (unsigned long)n, 1);
  for (i = 1; i <= n; i++) {
    xr[i-1] = sqrt(pow(xri[2*i-1], 2.0)+pow(xri[2*i], 2.0))*2.0/n;
    xi[i-1] = atan2(xri[2*i], xri[2*i-1]);
  }

  free_dvector(xri, 1, 2*n);
}


void GetInd(int n, int k, int *ind1, int *ind3)
{
  if (k == 1) {
    *ind1 = 2; *ind3 = 2;
  } else if (k == n/2+1) {
    *ind1 = n/2; *ind3 = n/2;
  } else {
    *ind1 = k - 1; *ind3 = k + 1;
  }

  if (trace)
    printf("GetInd: n = %d, k = %d, ind1 = %d, ind3 = %d\n",
	   n, k, *ind1, *ind3);
}


void GetInd1(int n, int k, int *ind1, int *ind3)
{
  if (k == 1) {
    *ind1 = 2; *ind3 = 2;
  } else if (k == n) {
    *ind1 = n - 1; *ind3 = n - 1;
  } else {
    *ind1 = k - 1; *ind3 = k + 1;
  }

  if (trace)
    printf("GetInd1: n = %d, k = %d, ind1 = %d, ind3 = %d\n",
	   n, k, *ind1, *ind3);
}


void GetPeak(int n, double *x, int *k)
{
  /* Locate peak in DFT spectrum */
  int ind1, ind2, ind3;
  double peak;

  peak = 0.0; *k = 1;
  for (ind2 = 1; ind2 <= n/2+1; ind2++) {
    GetInd(n, ind2, &ind1, &ind3);
    if (x[ind2-1] > peak && x[ind1-1] < x[ind2-1] &&
        x[ind3-1] < x[ind2-1])
    {
      peak = x[ind2-1]; *k = ind2;
    }
  }
}


void GetPeak1(int n, double *x, int *k)
{
  /* Locate peak in DFT spectrum */
  int ind1, ind2, ind3;
  double peak;

  peak = 0.0; *k = 1;
  for (ind2 = 1; ind2 <= n; ind2++) {
    GetInd1(n, ind2, &ind1, &ind3);
    if (x[ind2-1] > peak && x[ind1-1] < x[ind2-1] &&
        x[ind3-1] < x[ind2-1])
    {
      peak = x[ind2-1]; *k = ind2;
    }
  }
}


double Int2snu(int n, double *x, int k)
{
  /* Get frequency by nonlinear interpolation with two samples
     for sine window. The interpolation is:

              1              2 A(k)       1
         nu = - [ k - 1 + ------------- - - ] ,      k-1 <= N nu <= k
              N           A(k-1) + A(k)   2
   */
  int ind, ind1, ind3;
  double ampl1, ampl2;

  GetInd(n, k, &ind1, &ind3);
  if (x[ind3 - 1] > x[ind1 - 1]) {
    ampl1 = x[k - 1];
    ampl2 = x[ind3 - 1];
    ind = k;
  } else {
    ampl1 = x[ind1 - 1];
    ampl2 = x[k - 1];
    /* Interpolate in right direction for 0 frequency */
    if (k != 1)
      ind = ind1;
    else
      ind = 0;
  }
  /* Avoid division by zero */
  if (ampl1 + ampl2 != 0.0)
    return ((ind - 1 + 2 * ampl2 / (ampl1 + ampl2) - 0.5) / n);
  else
    return 0.0;
}


double Int2snu1(int n, double *x, int k)
{
  /* Get frequency by nonlinear interpolation with two samples
     for sine window. The interpolation is:

              1              2 A(k)       1
         nu = - [ k - 1 + ------------- - - ] ,      k-1 <= N nu <= k
              N           A(k-1) + A(k)   2
   */
  int ind, ind1, ind3;
  double ampl1, ampl2;

  GetInd1(n, k, &ind1, &ind3);
  if (x[ind3 - 1] > x[ind1 - 1]) {
    ampl1 = x[k - 1];
    ampl2 = x[ind3 - 1];
    ind = k;
  } else {
    ampl1 = x[ind1 - 1];
    ampl2 = x[k - 1];
    /* Interpolate in right direction for 0 frequency */
    if (k != 1)
      ind = ind1;
    else
      ind = 0;
  }
  /* Avoid division by zero */
  if (ampl1 + ampl2 != 0.0)
    return ((ind - 1 + 2 * ampl2 / (ampl1 + ampl2) - 0.5) / n);
  else
    return 0.0;
}


double Sinc(double omega)
{
  /*  Function to calculate:

                        sin( omega )
                        ------------
                           omega
  */
  if (omega != 0.0)
    return (sin(omega) / omega);
  else
    return 1.0;
}


double intsampl(int n, double *x, double nu, int k)
{
  /* Get amplitude by nonlinear interpolation for sine window. The
     distribution is given by:

                   1    sin pi ( k + 1/2 )     sin pi ( k - 1/2 )
           F(k) =  - ( -------------------- + -------------------- )
                   2      pi ( k + 1/2 )          pi ( k - 1/2 )
   */
  double corr;

  corr = (Sinc(M_PI * (k - 1 + 0.5 - nu * n)) +
          Sinc(M_PI * (k - 1 - 0.5 - nu * n))) / 2;
  return (x[k - 1] / corr);
}


double intsampl1(int n, double *x, double nu, int k)
{
  /* Get amplitude by nonlinear interpolation for sine window. The
     distribution is given by:

                   1    sin pi ( k + 1/2 )     sin pi ( k - 1/2 )
           F(k) =  - ( -------------------- + -------------------- )
                   2      pi ( k + 1/2 )          pi ( k - 1/2 )
   */
  double corr;

  corr = (Sinc(M_PI * (k - 1 + 0.5 - nu * n)) +
          Sinc(M_PI * (k - 1 - 0.5 - nu * n))) / 2;
  return (x[k - 1] / corr);
}


double linint(int n, int k, double nu, double *x)
{
  /* Get phase by linear interpolation for rectangular window
     with -pi <= phi <= pi */
  int     i;
  double  phi;
  double  xr[n], xi[n];

  for (i = 0; i < n; i++) {
    xr[i] = x[i]; xi[i] = 0.0;
  }
  FFT(n, xr, xi);
  phi = GetAngle(xr[k-1], xi[k-1]) - (n*nu-k+1)*M_PI;
  if (phi > M_PI)
    phi -= 2.0*M_PI;
  else if (phi < -M_PI)
    phi += 2.0*M_PI;
  return phi;
}

/****************************************************************************/
/* void FndRes(struct LOC_findres *LINK)

   Purpose:


   Input:
       none

   Output:
       none

   Return:
       none

   Global variables:
       none

   Specific functions:
       none

   Comments:
       none

****************************************************************************/
void FndRes(struct LOC_findres *LINK)
{
  int i, j, FORLIM, FORLIM1;
  double delta;

  FORLIM = LINK->n;
  for (i = 0; i <= FORLIM; i++) {
    FORLIM1 = LINK->n;
    for (j = -LINK->n; j <= FORLIM1; j++) {
      delta = fabs(i * LINK->nux + j * LINK->nuy);
      delta -= (int)delta;
      if (delta > 0.5)
	delta = 1 - delta;
      delta = fabs(delta - LINK->f);
      delta -= (int)delta;
      if (delta > 0.5)
	delta = 1 - delta;
      if (delta < LINK->eps) {
	if (abs(i) + abs(j) < LINK->n && (i != 0 || j >= 0)) {
	  LINK->found = true;
	  *LINK->nx = i;
	  *LINK->ny = j;
	}
      }
    }
  }
}


/****************************************************************************/
/* void FindRes(long n_, double nux_, double nuy_, double f_,
             long *nx_, long *ny_)

   Purpose:


   Input:
       none

   Output:
       none

   Return:
       none

   Global variables:
       none

   Specific functions:
       none

   Comments:
       none

****************************************************************************/
void FindRes(int n_, double nux_, double nuy_, double f_, int *nx_, int *ny_)
{
  /* Match f by a linear combination of nux and nuy */
  struct LOC_findres V;

  V.n = n_;
  V.nux = nux_;
  V.nuy = nuy_;
  V.f = f_;
  V.nx = nx_;
  V.ny = ny_;
  V.found = false;
  V.eps = FindRes_eps;
  do {
    V.eps = 10 * V.eps;
    FndRes(&V);
  } while (!V.found);
}


void GetPeaks(int n, double *x, int nf, double *nu, double *A)
{
  int i, k, ind1, ind3;

  for (i = 0; i < nf; i++) {
    GetPeak(n, x, &k);
    nu[i] = Int2snu(n, x, k); A[i] = intsampl(n, x, nu[i], k);
    /* Make peak flat to allow for new call */
    GetInd(n, k, &ind1, &ind3);
    if (x[ind1-1] > x[ind3-1])
      x[k-1] = x[ind1-1];
    else
      x[k-1] = x[ind3-1];
  }
}


void GetPeaks1(int n, double *x, int nf, double *nu, double *A)
{
  int i, k, ind1, ind3;

  for (i = 0; i < nf; i++) {
    GetPeak1(n, x, &k);
    nu[i] = Int2snu1(n, x, k); A[i] = intsampl1(n, x, nu[i], k);
    /* Make peak flat to allow for new call */
    GetInd1(n, k, &ind1, &ind3);
    if (x[ind1-1] > x[ind3-1])
      x[k-1] = x[ind1-1];
    else
      x[k-1] = x[ind3-1];
  }
}

/*******************************/
/* Routines for magnetic error */
/*******************************/

void SetTol(LatticeType &lat, int Fnum, double dxrms, double dyrms,
	    double drrms)
{
  int       i;
  long      k;
  MpoleType *M;

  for (i = 1; i <= lat.GetnKid(Fnum); i++) {
    k = lat.Elem_GetPos(Fnum, i);
    M = dynamic_cast<MpoleType*>(lat.elems[k]);
    M->PdSrms[X_] = dxrms;
    M->PdSrnd[X_] = normranf();
    M->PdSrms[Y_] = dyrms;
    M->PdSrnd[Y_] = normranf();
    M->PdTrms = drrms;
    M->PdTrnd = normranf();
    lat.SetdS(Fnum, i); lat.SetdT(Fnum, i);
  }
}


void Scale_Tol(LatticeType &lat, int Fnum, double dxrms, double dyrms,
	       double drrms)
{
  int       Knum;
  long int  loc;
  MpoleType *M;

  for (Knum = 1; Knum <= lat.GetnKid(Fnum); Knum++) {
    loc = lat.Elem_GetPos(Fnum, Knum);
    M = dynamic_cast<MpoleType*>(lat.elems[loc]);
    M->PdSrms[X_] = dxrms; M->PdSrms[Y_] = dyrms;
    M->PdTrms    = drrms;
    lat.SetdS(Fnum, Knum); lat.SetdT(Fnum, Knum);
  }
}


/****************************************************************************/
/* void SetaTol(int Fnum, int Knum, double dx, double dy, double dr)

   Purpose:
       Set a known random multipole displacement error

   Input:
       none

   Output:
       none

   Return:
       none

   Global variables:
       none

   Specific functions:
       none

   Comments:
       none

****************************************************************************/
void SetaTol(LatticeType &lat, int Fnum, int Knum, double dx, double dy,
	     double dr)
{
  long int  loc;
  MpoleType *M;

  loc = lat.Elem_GetPos(Fnum, Knum);
  M = dynamic_cast<MpoleType*>(lat.elems[loc]);
  M->PdSrms[0] = dx; M->PdSrnd[0] = 1e0;
  M->PdSrms[1] = dy; M->PdSrnd[1] = 1e0;
  M->PdTrms    = dr; M->PdTrnd    = 1e0;
  lat.SetdS(Fnum, Knum); lat.SetdT(Fnum, Knum);
}


void ini_aper(LatticeType &lat, const double Dxmin, const double Dxmax, 
              const double Dymin, const double Dymax) 
{ 
  int  k; 
 
  for (k = 0; k <= lat.conf.Cell_nLoc; k++) { 
    lat.elems[k]->maxampl[X_][0] = Dxmin;
    lat.elems[k]->maxampl[X_][1] = Dxmax; 
    lat.elems[k]->maxampl[Y_][0] = Dymin;
    lat.elems[k]->maxampl[Y_][1] = Dymax; 
  } 
} 
 
void set_aper(LatticeType &lat, const int Fnum, const double Dxmin,
	      const double Dxmax, const double Dymin, const double Dymax)
{
  int       i;
  long int  loc;

  for (i = 1; i <= lat.GetnKid(Fnum); i++) {
    loc = lat.Elem_GetPos(Fnum, i);
    lat.elems[loc]->maxampl[X_][0] = Dxmin;
    lat.elems[loc]->maxampl[X_][1] = Dxmax;
    lat.elems[loc]->maxampl[Y_][0] = Dymin;
    lat.elems[loc]->maxampl[Y_][1] = Dymax;
  }
}


void LoadApertures(LatticeType &lat, const char *ChamberFileName)
{
  char    line[128], FamName[32];
  long    Fnum;
  double  Xmin, Xmax, Ymin, Ymax;
  FILE    *ChamberFile;

  ChamberFile = file_read(ChamberFileName);

  do
    fgets(line, 128, ChamberFile);
  while (strstr(line, "#") != NULL);

  do {
    sscanf(line,"%s %lf %lf %lf %lf", FamName,&Xmin, &Xmax, &Ymin,&Ymax);
      Fnum = ElemIndex(FamName);
      if (Fnum > 0) set_aper(lat, Fnum, Xmin, Xmax, Ymin, Ymax);
  } while (fgets(line, 128, ChamberFile ) != NULL);

  fclose(ChamberFile);
}


// Load tolerances from the file
void LoadTolerances(LatticeType &lat, const char *TolFileName) 
{
  char    line[128], FamName[32];
  int     Fnum;
  double  dx, dy, dr;
  FILE    *tolfile;

  tolfile = file_read(TolFileName);

  do
    fgets(line, 128, tolfile);
  while (strstr(line, "#") != NULL);

  do {
    if (strstr(line, "#") == NULL) {
      sscanf(line,"%s %lf %lf %lf", FamName, &dx, &dy, &dr);
      Fnum = ElemIndex(FamName);
      if (Fnum > 0) {
	SetTol(lat, Fnum, dx, dy, dr);
      } else {
	printf("LoadTolerances: undefined element %s\n", FamName);
	exit_(1);
      }
    }
  } while (fgets(line, 128, tolfile) != NULL);

  fclose(tolfile);
}


// Load tolerances from the file
void ScaleTolerances(LatticeType &lat, const char *TolFileName,
		     const double scl) 
{
  char    line[128], FamName[32];
  int     Fnum;
  double  dx, dy, dr;
  FILE    *tolfile;

  tolfile = file_read(TolFileName);

  do
    fgets(line, 128, tolfile);
  while (strstr(line, "#") != NULL);
  
  do {
    if (strstr(line, "#") == NULL) {
      sscanf(line,"%s %lf %lf %lf", FamName, &dx, &dy, &dr);
      Fnum = ElemIndex(FamName);
      if (Fnum > 0) {
	Scale_Tol(lat, Fnum, scl*dx, scl*dy, scl*dr);
      } else {
	printf("ScaleTolerances: undefined element %s\n", FamName);
	exit_(1);
      }
    }
  } while (fgets(line, 128, tolfile) != NULL);
  fclose(tolfile);
}


void SetKpar(LatticeType &lat, int Fnum, int Knum, int Order, double k)
{
  long int  loc;
  MpoleType *M;

  loc = lat.Elem_GetPos(Fnum, Knum);
  M = dynamic_cast<MpoleType*>(lat.elems[loc]);
  M->PBpar[Order+HOMmax] = k;
  lat.SetPB(Fnum, Knum, Order);
}


void SetL(LatticeType &lat, int Fnum, int Knum, double L)
{

  lat.elems[lat.Elem_GetPos(Fnum, Knum)]->PL = L;
}


void SetL(LatticeType &lat, int Fnum, double L)
{
  int  i;

  for (i = 1; i <= lat.GetnKid(Fnum); i++)
    lat.elems[lat.Elem_GetPos(Fnum, i)]->PL = L;
}


void SetdKpar(LatticeType &lat, int Fnum, int Knum, int Order, double dk)
{
  long int  loc;
  MpoleType *M;

  loc = lat.Elem_GetPos(Fnum, Knum);
  M = dynamic_cast<MpoleType*>(lat.elems[loc]);
  M->PBpar[Order+HOMmax] += dk;
  lat.SetPB(Fnum, Knum, Order);
}


void SetdKLpar(LatticeType &lat, int Fnum, int Knum, int Order, double dkL)
{
  long int  loc;
  MpoleType *M;

  loc = lat.Elem_GetPos(Fnum, Knum);
  M = dynamic_cast<MpoleType*>(lat.elems[loc]);
  if (lat.elems[loc]->PL != 0e0)
    M->PBpar[Order + HOMmax] += dkL/lat.elems[loc]->PL;
  else
    M->PBpar[Order + HOMmax] += dkL;
  lat.SetPB(Fnum, Knum, Order);
}


void SetdKrpar(LatticeType &lat, int Fnum, int Knum, int Order, double dkrel)
{
  long int  loc;
  MpoleType *M;

  loc = lat.Elem_GetPos(Fnum, Knum);
  M = dynamic_cast<MpoleType*>(lat.elems[loc]);
  if (Order == Dip && M->Pthick == thick)
    M->PBpar[Dip+HOMmax] += dkrel*M->Pirho;
  else
    M->PBpar[Order+HOMmax] += dkrel*M->PBpar[Order+HOMmax];
  lat.SetPB(Fnum, Knum, Order);
}


void Setbn(LatticeType &lat, int Fnum, int order, double bn)
{
  int i;

  for (i = 1; i <=  lat.GetnKid(Fnum); i++)
    SetKpar(lat, Fnum, i, order, bn);
}


void SetbnL(LatticeType &lat, int Fnum, int order, double bnL)
{
  int i;

  for (i = 1; i <= lat.GetnKid(Fnum); i++)
    SetKLpar(lat, Fnum, i, order, bnL);
}


void Setdbn(LatticeType &lat, int Fnum, int order, double dbn)
{
  int i;

  for (i = 1; i <= lat.GetnKid(Fnum); i++)
    SetdKpar(lat, Fnum, i, order, dbn);
}


void SetdbnL(LatticeType &lat, int Fnum, int order, double dbnL)
{
  int i;

  for (i = 1; i <= lat.GetnKid(Fnum); i++) {
    SetdKLpar(lat, Fnum, i, order, dbnL);
  }
}


void Setbnr(LatticeType &lat, int Fnum, long order, double bnr)
{
  int  i;

  for (i = 1; i <= lat.GetnKid(Fnum); i++)
    SetdKrpar(lat, Fnum, i, order, bnr);
}


void SetbnL_sys(LatticeType &lat, int Fnum, int Order, double bnL_sys)
{
  int       Knum;
  long int  loc;
  MpoleType *M;

  for (Knum = 1; Knum <= lat.GetnKid(Fnum); Knum++) {
    loc = lat.Elem_GetPos(Fnum, Knum);
    M = dynamic_cast<MpoleType*>(lat.elems[loc]);
    if (lat.elems[loc]->PL != 0.0)
      M->PBsys[Order+HOMmax] = bnL_sys/lat.elems[loc]->PL;
    else
      M->PBsys[Order+HOMmax] = bnL_sys;
    lat.SetPB(Fnum, Knum, Order);
  }
}


void set_dbn_rel(LatticeType &lat, const int type, const int n,
		 const double dbn_rel)
{
  long int  j;
  double    dbn;
  MpoleType *M;

  printf("\n");
  printf("Setting Db_%d/b_%d = %6.1e for:\n", n, type, dbn_rel);
  printf("\n");
  for (j = 0; j <= lat.conf.Cell_nLoc; j++)
    M = dynamic_cast<MpoleType*>(lat.elems[j]);
    if ((lat.elems[j]->Pkind == Mpole) && (M->n_design == type)) {
      printf("%s\n", lat.elems[j]->PName);
      dbn = dbn_rel*M->PBpar[type+HOMmax];
      M->PBrms[n+HOMmax] = dbn;
      M->PBrnd[n+HOMmax] = normranf();
      lat.SetPB(lat.elems[j]->Fnum, lat.elems[j]->Knum, n);
    }
}


double GetL(LatticeType &lat, int Fnum, int Knum)
{
  return (lat.elems[lat.Elem_GetPos(Fnum, Knum)]->PL);
}


double GetKLpar(LatticeType &lat, int Fnum, int Knum, int Order)
{
  long int  loc;
  MpoleType *M;

  loc = lat.Elem_GetPos(Fnum, Knum);
  M = dynamic_cast<MpoleType*>(lat.elems[loc]);
  if (lat.elems[loc]->PL != 0e0)
    return (M->PBpar[Order+HOMmax]*lat.elems[loc]->PL);
  else
    return (M->PBpar[Order+HOMmax]);
}


void SetdKLrms(LatticeType &lat, int Fnum, int Order, double dkLrms)
{
  long int  Knum, loc;
  MpoleType *M;

  for (Knum = 1; Knum <= lat.GetnKid(Fnum); Knum++) {
    loc = lat.Elem_GetPos(Fnum, Knum);
    M = dynamic_cast<MpoleType*>(lat.elems[loc]);
    if (lat.elems[loc]->PL != 0e0)
      M->PBrms[Order+HOMmax] = dkLrms/lat.elems[loc]->PL;
    else
      M->PBrms[Order+HOMmax] = dkLrms;
    M->PBrnd[Order+HOMmax] = normranf();
    lat.SetPB(Fnum, Knum, Order);
  }
}


void Setdkrrms(LatticeType &lat, int Fnum, int Order, double dkrrms)
{
  long int  Knum, loc;
  MpoleType *M;

  for (Knum = 1; Knum <= lat.GetnKid(Fnum); Knum++) {
    loc = lat.Elem_GetPos(Fnum, Knum);
    M = dynamic_cast<MpoleType*>(lat.elems[loc]);
    if (Order == Dip && M->Pthick == thick)
      M->PBrms[Dip+HOMmax] = dkrrms*M->Pirho;
    else
      M->PBrms[Order+HOMmax]
	= dkrrms*M->PBpar[Order+HOMmax];
    M->PBrnd[Order+HOMmax] = normranf();
    lat.SetPB(Fnum, Knum, Order);
  }
}


void SetKL(LatticeType &lat, int Fnum, int Order)
{
  long int  Knum;

  for (Knum = 1; Knum <= lat.GetnKid(Fnum); Knum++)
    lat.SetPB(Fnum, Knum, Order);
}


void set_dx(LatticeType &lat, const int type, const double sigma_x,
	    const double sigma_y)
{
  long int  j;
  MpoleType *M;

  printf("\n");
  printf("Setting sigma_x,y = (%6.1e, %6.1e) for b_%d:\n",
	 sigma_x, sigma_y, type);
  printf("\n");
  for (j = 0; j <= lat.conf.Cell_nLoc; j++)
    M = dynamic_cast<MpoleType*>(lat.elems[j]);
    if ((lat.elems[j]->Pkind == Mpole) && (M->n_design == type)) {
      printf("%s\n", lat.elems[j]->PName);
      M->PdSrms[X_] = sigma_x;
      M->PdSrms[Y_] = sigma_y;
      M->PdSrnd[X_] = normranf();
      M->PdSrnd[Y_] = normranf();
      lat.SetdS(lat.elems[j]->Fnum, lat.elems[j]->Knum);
    }
}


void SetBpmdS(LatticeType &lat, int Fnum, double dxrms, double dyrms)
{
  long int  Knum, loc;

  for (Knum = 1; Knum <= lat.GetnKid(Fnum); Knum++) {
    loc = lat.Elem_GetPos(Fnum, Knum);
    lat.elems[loc]->dS[X_] = normranf()*dxrms;
    lat.elems[loc]->dS[Y_] = normranf()*dyrms;
  }
}


/****************************************/
/* Routines for closed orbit correction */
/****************************************/


/****************************************************************************/
void codstat(LatticeType &lat, double *mean, double *sigma, double *xmax,
	     long lastpos, bool all)
{
  long    i, n, loc;
  int     j;
  Vector2 sum, sum2;

  for (j = 0; j < 2; j++) {
    sum[j] = 0e0; sum2[j] = 0e0; xmax[j] = 0e0;
  }

  n = 0;
  if (all) {
    for (i = 0; i < lastpos; i++) {
      n++;
      for (j = 0; j < 2; j++) {
	sum[j] += lat.elems[i]->BeamPos[j*2];
	sum2[j] += sqr(lat.elems[i]->BeamPos[j*2]);
	xmax[j] = max(xmax[j], fabs(lat.elems[i]->BeamPos[j*2]));
      }
    }
  } else {
    for (i = 1; i <= n_bpm_[X_]; i++) {
      n++;
      for (j = 0; j < 2; j++) {
	loc = bpms_[j][i];
	sum[j] += lat.elems[loc]->BeamPos[j*2];
	sum2[j] += sqr(lat.elems[loc]->BeamPos[j*2]);
	xmax[j] = max(xmax[j], fabs(lat.elems[loc]->BeamPos[j*2]));
      }
    }
  }

  for (j = 0; j < 2; j++) {
    if (n != 0)
      mean[j] = sum[j] / n;
    else
      mean[j] = -1e0;
    if (n != 0 && n != 1) {
      sigma[j] = (n*sum2[j]-sqr(sum[j]))/(n*(n-1e0));
    } else
      sigma[j] = 0e0;
    if (sigma[j] >= 0e0)
      sigma[j] = sqrt(sigma[j]);
    else
      sigma[j] = -1e0;
  }
}

/****************************************************************************/
/* void CodStatBpm(double *mean, double *sigma, double *xmax, long lastpos,
                long bpmdis[])

   Purpose:
       Get statistics for  closed orbit

   Input:
       none

   Output:
       none

   Return:
       none

   Global variables:
       none

   Specific functions:
       none

   Comments:
       none

****************************************************************************/
void CodStatBpm(LatticeType &lat, double *mean, double *sigma, double *xmax,
		long lastpos, long bpmdis[])
{
  long     i, j, m, n;
  Vector2  sum, sum2;
  double   TEMP;

  m= n= 0;
  for (j = 0; j <= 1; j++) {
    sum[j] = 0.0; sum2[j] = 0.0; xmax[j] = 0.0;
  }

  for (i = 0; i <= lastpos; i++) {
    if (lat.elems[i]->Fnum == lat.conf.bpm) {
      if (! bpmdis[m]) {
	for (j = 1; j <= 2; j++) {
	  sum[j - 1] += lat.elems[i]->BeamPos[j * 2 - 2];
	  TEMP = lat.elems[i]->BeamPos[j * 2 - 2];
	  sum2[j - 1] += TEMP * TEMP;
	  xmax[j - 1] =
	    max(xmax[j - 1], fabs(lat.elems[i]->BeamPos[j * 2 - 2]));
	}
        n++;
      }
      m++;
    }
  }
  for (j = 0; j <= 1; j++) {
    if (n != 0)
      mean[j] = sum[j] / n;
    else
      mean[j] = 0.0;
    if (n != 0 && n != 1) {
      TEMP = sum[j];
      sigma[j] = (n * sum2[j] - TEMP * TEMP) / (n * (n - 1.0));
    } else
      sigma[j] = 0.0;
    if (sigma[j] >= 0.0)
      sigma[j] = sqrt(sigma[j]);
    else
      sigma[j] = 0.0;
  }
}



/****************************************************************************/
/* double digitize(double x, double maxkick, double maxsamp)

   Purpose:
       Map x onto the integer interval (-maxsamp ... maxsamp) where maxsamp
       corresponds maxkick.


   Input:
       none

   Output:
       none

   Return:
       none

   Global variables:
       none

   Specific functions:
       none

   Comments:
       none

****************************************************************************/
double digitize(double x, double maxkick, double maxsamp)
{
  if (maxkick>0.)
    if (maxsamp>1.)
      return Sgn(x)*maxkick/maxsamp *min(floor(fabs(x)/maxkick*maxsamp),maxsamp-1.);
    else {
      return Sgn(x)*min(fabs(x),maxkick);
    }
  else
    return x;
}


/****************************************************************************/
/* double digitize2(long plane, long inum, double x, double maxkick, double maxsamp)

   Purpose:


   Input:
       none

   Output:
       none

   Return:
       none

   Global variables:
       none

   Specific functions:
       none

   Comments:
       none

****************************************************************************/
double xmemo[2][2000];

double digitize2(long plane, long inum, double x, double maxkick, double maxsamp)
{
  double xint;

  if (maxkick>0.)
    if (maxsamp>1.)
    {
      xint=min(floor(fabs(x)/maxkick*maxsamp),maxsamp-1.);

      if(fabs(xint-xmemo[inum][plane]) >=1.)
      {
        xmemo[inum][plane]=xint;
      }
      else
      {
        xmemo[inum][plane]+=0.1;
        xint=min(xmemo[inum][plane],maxsamp-1.);
      }
      return Sgn(x)*maxkick/maxsamp*xint;
    }
    else
    {
      return Sgn(x)*min(fabs(x),maxkick);
    }
  else
    return x;
}


// MATH ROUTINE a mettre dans mathlib.c

/****************************************************************************/
/*  void GetMean(n, x)

   Purpose:
       Get out the mean value of vector x

   Input:
       n vector size
       x vector to get out the mean value

   Output:
       none

   Return:
       none

   Global variables:
       none

   Specific functions:
       none

   Comments:
       to be moved in mathlib

****************************************************************************/
void GetMean(long n, double *x)
{
  long i;
  double mean = 0e0;

  if ( n < 1 )
  {
    fprintf(stdout,"GetMean: error wrong vector size n=%ld\n",n);
    exit_(1);
  }
  for (i = 0; i < n; i++)
    mean += x[i];
  mean /= n;
  for (i = 0; i < n; i++)
    x[i] = x[i] - mean;
}

/****************************************************************************/
/* double Fract(double x)

   Purpose:
      Gets fractional part of x 

   Input:
       none

   Output:
       none

   Return:
       none

   Global variables:
       none

   Specific functions:
       none

   Comments:
       none

****************************************************************************/
double Fract(double x)
{
  return (x - (long int)x);
}

/****************************************************************************/
/* double Sgn (double x)

   Purpose:
       Gets sign of x

   Input:
       none

   Output:
       0  if zero
       1  if positive
       -1 if negative

   Return:
       none

   Global variables:
       none

   Specific functions:
       none

   Comments:
       none

****************************************************************************/
double Sgn (double x)
{
  return (x == 0.0 ? 0.0 : (x > 0.0 ? 1.0 : -1.0));
}


/** function from soleilcommon.c **/

void Read_Lattice(LatticeType &lat, const char *fic)
{
  bool     status;
  char     fic_maille[S_SIZE+4] = "", fic_erreur[S_SIZE+4] = "";
  int      i;
  double   dP = 0.0;
  Vector2  beta, alpha, eta, etap;
  psVector codvect;

  const double RFacceptance = 0.060001; // soleil energy acceptance

  strcpy(fic_maille, fic); strcpy(fic_erreur, fic);

  /* generation automatique du nom du fichier maille et erreur */
  strcat(fic_maille, ".lat"); strcat(fic_erreur, ".lax");

  /* Initialisation de Tracy */

  t2init();

  /* open the lattice Input file  */

  if ((fi = fopen(fic_maille, "r")) == NULL) {
    fprintf(stdout, "ReadLattice: Error while opening file %s \n", fic_maille);
    fprintf(stdout, "The lattice file has to end by .lat \n");
    exit_(1);
  }

  /* opens the lattice Output file */
  if ((fo = fopen(fic_erreur, "w")) == NULL) {
    fprintf(stdout, "ReadLattice: Error while opening file %s \n", fic_erreur);
    exit_(1);
  }

  /* Reads lattice and set principle parameters
   * Energy CODeps and energy offset
   * print statistics
   */
  status = lat.Lattice_Read(fi, fo);

  if (status == false) {
    std::cout << "Lattice_Read function has returned false" << std::endl;
    std::cout << "See file " << fic_erreur << std::endl;
    exit_(1);
  }
  std::cout << "Lattice file: " << fic_maille << std::endl;

  /* initializes cell structure: construction of the RING */
  /* Creator of all the matrices for each element         */
  lat.Lat_Init();

  lat.conf.H_exact        = false; // Small Ring Hamiltonian
  lat.conf.Cart_Bend      = false; // Cartesian Bend
  lat.conf.dip_edge_fudge = true;

  if (lat.conf.RingType == 1) { // for a ring
    /* define x/y physical aperture  */
    lat.ChamberOff();

    /* Defines global variables for Tracy code */
    lat.conf.quad_fringe    = false; // quadrupole fringe fields on/off
    lat.conf.EPU            = false; // Elliptically Polarizing Undulator
    lat.conf.Cavity_on      = false; /* Cavity on/off */
    lat.conf.radiation      = false; /* radiation on/off */
    lat.conf.IBS            = false; /* diffusion on/off */
    lat.conf.emittance      = false; /* emittance  on/off */
    lat.conf.pathlength     = false; /* Path lengthening computation */
    lat.conf.CODimax        = 40;    /* maximum number of iterations for COD
				    algo */
    lat.conf.delta_RF = RFacceptance;/* energy acceptance for SOLEIL */
  } else {   
    // for transfer lines
    /* Initial settings : */
    beta[X_] = 8.1; alpha[X_] = 0.0; beta[Y_] = 8.1; alpha[Y_] = 0.0;
    eta[X_] = 0.0; etap[X_] = 0.0; eta[Y_] = 0.0; etap[Y_] = 0.0;

    for (i = 0; i < ss_dim; i++) {
      codvect[i] = 0.0; lat.conf.CODvect[i] = codvect[i];
    }
    dP = codvect[delta_];

    /* Defines global variables for Tracy code */
    lat.conf.Cavity_on    = false; /* Cavity on/off */
    lat.conf.radiation    = false; /* radiation on/off */
    lat.conf.emittance    = false; /* emittance  on/off */
    lat.conf.pathlength   = false; /* Path lengthening computation */
    lat.conf.CODimax      = 10;    /* maximum number of iterations for COD
				     algo */
    lat.conf.delta_RF = RFacceptance; /* 6% + epsilon energy acceptance
                                            for SOLEIL */
    lat.conf.dPparticle = dP;

    lat.ChamberOff();

    lat.ttwiss(alpha, beta, eta, etap, dP);
  }
}


/****************************************************************************/
/* void GetChromTrac(long Nb, long Nbtour, double emax,
                     double *xix, double *xiz)

   Purpose:
       Computes chromaticities by tracking

   Input:
       Nb      point number
       Nbtour  turn number
       emax    energy step

   Output:
       xix horizontal chromaticity
       xiz vertical chromaticity

   Return:
       none

   Global variables:
       trace

   Specific functions:
       Trac_Simple, Get_NAFF

   Comments:
       27/04/03 chromaticities are now output arguments

****************************************************************************/
#define nterm 2
void GetChromTrac(LatticeType &lat, long Nb, long Nbtour, double emax,
		  double *xix, double *xiz)
{
  bool    status = true;
  int     nb_freq[2] = { 0, 0 };  /* frequency number to look for */
  int     i = 0;
  double  Tab[6][NTURN], fx[nterm], fz[nterm], nux1, nux2, nuz1, nuz2;

  double  x = 1e-6, xp = 0.0, z = 1e-6, zp = 0.0;
  double  x0 = 1e-6, xp0 = 0.0, z0 = 1e-6, zp0 = 0.0;

  /* initializations */
  for (i = 0; i < nterm; i++) {
    fx[i] = 0.0; fz[i] = 0.0;
  }
  /* end init */

  /* Tracking for delta = emax and computing tunes */
  x = x0; xp = xp0; z = z0; zp = zp0;

  Trac_Simple(lat, x, xp, z, zp, emax, 0.0, Nbtour, Tab, &status);
  Get_NAFF(nterm, Nbtour, Tab, fx, fz, nb_freq);

  nux1 = (fabs (fx[0]) > 1e-8 ? fx[0] : fx[1]); nuz1 = fz[0];

  if (trace)
    fprintf(stdout,
       "\n Entering routine for chroma using tracking\n");
  if (trace)
    fprintf(stdout, "emax= % 10.6e nux1=% 10.6e nuz1= % 10.6e\n",
       emax, nux1, nuz1);

  /* Tracking for delta = -emax and computing tunes */
  x = x0; xp = xp0; z = z0; zp = zp0;

  Trac_Simple(lat, x, xp, z, zp, -emax, 0.0, Nbtour, Tab, &status);
  Get_NAFF(nterm, Nbtour, Tab, fx, fz, nb_freq);

  if (trace)
    fprintf(stdout, "nturn=%6ld x=% 10.5g xp=% 10.5g z=% 10.5g zp=% 10.5g"
	    " delta=% 10.5g ctau=% 10.5g \n",
	    Nbtour,
	    Tab[0][Nbtour-1], Tab[1][Nbtour-1],
	    Tab[2][Nbtour-1], Tab[3][Nbtour-1],
	    Tab[4][Nbtour-1], Tab[5][Nbtour-1]);

  nux2 = (fabs(fx[0]) > 1e-8 ? fx[0] : fx[1]); nuz2 = fz[0];

  if (trace)
    fprintf(stdout, "emax= % 10.6e nux2= % 10.6e nuz2= % 10.6e\n",
	    -emax, nux2, nuz2);

  /* Computing chromaticities */
  *xix = (nux2-nux1)*0.5/emax; *xiz = (nuz2-nuz1)*0.5/emax;

  if (trace)
    fprintf (stdout, " Exiting routine for chroma using tracking\n\n");
}
#undef nterm

/****************************************************************************/
/* void GetTuneTrac(long Nbtour, double emax, double *nux, double *nuz)

   Purpose:
       Computes chromaticities by tracking

   Input:
       Nb      point number
       Nbtour  turn number
       emax    energy step

   Output:
       none

   Return:
       none

   Global variables:
       trace

   Specific functions:
       Trac_Simple, Get_NAFF

   Comments:
       none

****************************************************************************/
#define nterm  2
void GetTuneTrac(LatticeType &lat, long Nbtour, double emax, double *nux,
		 double *nuz)
{
  double Tab[6][NTURN], fx[nterm], fz[nterm];
  int nb_freq[2];
  bool status;

  double x = 1e-6, xp = 0.0, z = 1e-6, zp = 0.0;

  Trac_Simple(lat, x, xp, z, zp, emax, 0.0, Nbtour, Tab, &status);
  Get_NAFF(nterm, Nbtour, Tab, fx, fz, nb_freq);

  *nux = (fabs (fx[0]) > 1e-8 ? fx[0] : fx[1]);
  *nuz = fz[0];
}
#undef nterm


/****************************************************************************/
/* void findcodS(double dP)

   Purpose: 
       Search for the closed orbit using a numerical method
       Algo: Newton_Raphson method
             Quadratic convergence
             May need a guess starting point
             Simple precision algorithm

   Input:
       dP energy offset

   Output:
       none

   Return:
       none

   Global variables:
       none

   specific functions:
       Newton_Raphson

   Comments:
       Method introduced because of bad convergence of da for ID using RADIA maps

****************************************************************************/
void findcodS(LatticeType &lat, double dP)
{
  double        *vcod;
  psVector       x0;
  const int    ntrial = 40;  // maximum number of trials for closed orbit
  const double  tolx = 1e-8;  // numerical precision
  int          k;
  int          dim;    // 4D or 6D tracking
  long         lastpos;

  vcod = dvector(1, 6);
      
  // starting point
  for (k = 1; k <= 6; k++)
    vcod[k] = 0.0;  
  
  vcod[5] = dP;  // energy offset 
    
  if (lat.conf.Cavity_on){
      dim = 6;   /* 6D tracking*/
    fprintf(stdout,"Error looking for cod in 6D\n");
    exit_(1);
  }
    else{
      dim = 4; /* 4D tracking */
      vcod[1] = lat.elems[0]->Eta[0]*dP; vcod[2] = lat.elems[0]->Etap[0]*dP;
      vcod[3] = lat.elems[0]->Eta[1]*dP; vcod[4] = lat.elems[0]->Etap[1]*dP;
  }
  
  Newton_RaphsonS(lat, ntrial, vcod, dim, tolx);

  if (lat.conf.codflag == false)
    fprintf(stdout, "Error No COD found\n");
  if (trace) {
    for (k = 1; k <= 6; k++)
      x0[k-1] = vcod[k];
    fprintf(stdout, "Before cod % .5e % .5e % .5e % .5e % .5e % .5e \n",
	    x0[0], x0[1], x0[2], x0[3], x0[4], x0[5]);
    lat.Cell_Pass(0, lat.conf.Cell_nLoc, x0, lastpos);
    fprintf(stdout, "After  cod % .5e % .5e % .5e % .5e % .5e % .5e \n",
	    x0[0], x0[1], x0[2], x0[3], x0[4], x0[5]);
    lat.Cell_Pass(0, lat.conf.Cell_nLoc, x0, lastpos);
  }
  free_dvector(vcod,1,6);
}


/****************************************************************************/
/* void computeFandJS(double *x, int n, double **fjac, double *fvect)

   Purpose:
       Simple precision algo
       Tracks x over one turn. And computes the Jacobian matrix of the 
       transformation by numerical differentiation.
       using forward difference formula : faster but less accurate
       using symmetric difference formula

   Input:
       x vector for evaluation
       n dimension 4 or 6

   Output:
      fvect transport of x over one turn
      fjac  Associated jacobian matrix      

   Return:
       none

   Global variables:
       none

   specific functions:
       none

   Comments:
       none

****************************************************************************/

void computeFandJS(LatticeType &lat, double *x, int n, double **fjac,
		   double *fvect)
{
  int     i, k;
  long    lastpos = 0L;
  psVector  x0, fx, fx1, fx2;

  const double deps = 1e-8;  //stepsize for numerical differentiation

  for (i = 1; i <= 6; i++)
    x0[i - 1] = x[i];
  
  lat.Cell_Pass(0, lat.conf.Cell_nLoc, x0, lastpos);

  for (i = 1; i <= n; i++)
  {
    fvect[i] = x0[i - 1];
    fx[i - 1] = x0[i - 1];
  }

  // compute Jacobian matrix by numerical differentiation
  for (k = 0; k < n; k++)
  {
    for (i = 1; i <= 6; i++)
      x0[i - 1] = x[i];
    x0[k] += deps;  // differential step in coordinate k

    lat.Cell_Pass(0, lat.conf.Cell_nLoc, x0, lastpos);  // tracking along the ring
    for (i = 1; i <= 6; i++)
      fx1[i - 1] = x0[i - 1];

    for (i = 1; i <= 6; i++)
      x0[i - 1] = x[i];
    x0[5] = 0.0;
    x0[k] -= deps;  // differential step in coordinate k

    lat.Cell_Pass(0, lat.conf.Cell_nLoc, x0, lastpos);  // tracking along the ring
    for (i = 1; i <= 6; i++)
      fx2[i - 1] = x0[i - 1];

    for (i = 1; i <= n; i++)  // symmetric difference formula
      fjac[i][k + 1] = 0.5 * (fx1[i - 1] - fx2[i - 1]) / deps;
    //~ for (i = 1; i <= n; i++) // forward difference formula
    //~ fjac[i][k + 1] = (float) ((x0[i - 1] - fx[i - 1]) / deps);  
  }
}


/****************************************************************************/
/* void Newton_RaphsonS(int ntrial,double x[],int n,double tolx, double tolf)
 
   Purpose:
       Newton_Rapson algorithm from Numerical Recipes
       single precision algorithm
       Robustess: quadratic convergence
       Hint: for n-dimensional problem, the algo can be stuck on local minimum
             In this case, it should be enough to provide a resonable starting 
             point.

       Method:
         look for closed orbit solution of f(x) = x
         This problems is equivalent to finding the zero of g(x)= f(x) - x
         g(x+h) ~= f(x) - x + (Jacobian(f) -Id) h + O(h*h)
         Then at first order we solve h:
             h = - inverse(Jacobian(f) -Id) * (f(x)-x)
            the new guess is then xnew = x + h
         By iteration, this converges quadratically.
     
     The algo is stopped whenever  |x -xnew| < tolx     

         f(x) is computes by tracking over one turn
     Jacobian(f) is computed numerically by numerical differentiation
     These two operations are provided by the function computeFandJ

   Input:
       ntrial number of iterations for closed zero search
       n number of dimension 4 or 6
     x intial guess for the closed orbit
       tolx tolerance over the solution x
       tolf tolerance over the evalution f(x)  

   Output:
       x closed orbit  

   Return:
       none

   Global variables:
       status

   specific functions:
       computeFandJS
     ludcmp,lubksb 

   Comments:
       none

****************************************************************************/

void Newton_RaphsonS(LatticeType &lat, int ntrial, double x[], int n,
		     double tolx)
{
  int    k, i, *indx;
  double  errx, d, *bet, *fvect, **alpha;

  errx = 0.0;
  // NR arrays start from 1 and not 0 !!!       
  indx = ivector(1, n);
  bet = dvector(1, n);
  fvect = dvector(1, n);
  alpha = dmatrix(1, n, 1, n);

  for (k = 1; k <= ntrial; k++) {      // loop over number of iterations
    // supply function values at x in fvect and Jacobian matrix in fjac
    computeFandJS(lat, x, n, alpha, fvect);

    // Jacobian -Id
    for (i = 1; i <= n; i++)
      alpha[i][i] -= 1.0;
    for (i = 1; i <= n; i++)
      bet[i] = x[i] - fvect[i];  // right side of linear equation
    // solve linear equations using LU decomposition using NR routines
    dludcmp(alpha, n, indx, &d);
    dlubksb(alpha, n, indx, bet);
    errx = 0.0;  // check root convergence
    for (i = 1; i <= n; i++) {    // update solution
      errx += fabs(bet[i]);
      x[i] += bet[i];
    }

    if (trace)
      fprintf(stdout,
         "%02d: cod % .5e % .5e % .5e % .5e % .5e % .5e  errx =% .5e\n",
         k, x[1], x[2], x[3], x[4], x[5], x[6], errx);
    if (errx <= tolx) {
      lat.conf.codflag = true;
      break;
    }
  }
  // check whever closed orbit found out
  if ((k >= ntrial) && (errx >= tolx * 100)) lat.conf.codflag = false;

  free_dmatrix(alpha,1,n,1,n); free_dvector(bet,1,n); free_dvector(fvect,1,n);
  free_ivector(indx,1,n);
}


void rm_mean(long int n, double x[])
{
  long int  i;
  double    mean;

  mean = 0.0;
  for (i = 0; i < n; i++)
    mean += x[i];
  mean /= n;
  for (i = 0; i < n; i++)
    x[i] -= mean;
}
