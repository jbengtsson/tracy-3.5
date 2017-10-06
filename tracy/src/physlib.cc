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

/**** same as asctime in C without the \n at the end****/
char *asctime2(const struct tm *timeptr)
{
    // terminated with \0.
    static char wday_name[7][4] = {
        "Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"
    };
    // terminated with \0.
    static char mon_name[12][4] = {
        "Jan", "Feb", "Mar", "Apr", "May", "Jun",
        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
    };
    static char result[26];

    sprintf(result, "%.3s %.3s%3d %.2d:%.2d:%.2d %d",
        wday_name[timeptr->tm_wday],
        mon_name[timeptr->tm_mon],
        timeptr->tm_mday, timeptr->tm_hour,
        timeptr->tm_min, timeptr->tm_sec,
        1900 + timeptr->tm_year);
    return result;
}

/** Get time and date **/
struct tm* GetTime()
{
  struct tm *whattime;
  /* Get time and date */
  time_t aclock;
  time(&aclock);                 /* Get time in seconds */
  whattime = localtime(&aclock);  /* Convert time to struct */
  return whattime;
}


void printglob(void)
{
  printf("\n***************************************************************"
	 "***************\n");
  printf("\n");
  printf("  dPcommon     =  %9.3e  dPparticle   =  %9.3e"
	 "  Energy [GeV] = %.3f\n",
         globval.dPcommon, globval.dPparticle, globval.Energy);
  printf("  MaxAmplx [m] = %9.3e  MaxAmply [m] = %9.3e"
	 "  RFAccept [%%] = \xB1%4.2f\n",
         Lattice.Cell[0].maxampl[X_][0], Lattice.Cell[0].maxampl[Y_][0],
	 globval.delta_RF*1e2);
  printf(" Cavity_On    =  %s    ", globval.Cavity_on ? "TRUE " : "FALSE");
  printf("  Radiation_On = %s     \n",
	 globval.radiation ? "TRUE " : "FALSE");
  printf("  bpm          =  %3d        qt           = %3d        ",
	 globval.bpm, globval.qt);
  printf(" Chambre_On   = %s     \n", status.chambre ? "TRUE " : "FALSE");
  printf("  hcorr        =  %3d        vcorr        = %3d\n\n",
	 globval.hcorr, globval.vcorr);
  printf("  alphac       =   %22.16e\n", globval.Alphac); 
  printf("  nux          =  %19.16f      nuz  =  %19.16f",
         globval.TotalTune[X_], globval.TotalTune[Y_]);
  if (globval.Cavity_on)
    printf("  omega  = %11.9f\n", globval.Omega);
  else {
    printf("\n");
    printf("  ksix         = %10.6f                ksiz = %10.6f\n",
            globval.Chrom[X_], globval.Chrom[Y_]);
  }
  printf("\n");
  printf("  OneTurn matrix:\n");
  printf("\n");
  prtmat(2*DOF, globval.OneTurnMat);
  fflush(stdout);
}


double int_curly_H1(long int n)
{
  /* Integration with Simpson's Rule */

  double    curly_H;
  Vector2   alpha[3], beta[3], nu[3], eta[3], etap[3];

  // only works for matrix style calculations
  get_twiss3(n, alpha, beta, nu, eta, etap);

  curly_H = (get_curly_H(alpha[0][X_], beta[0][X_], eta[0][X_], etap[0][X_])
            +4.0*get_curly_H(alpha[1][X_], beta[1][X_],
                             eta[1][X_], etap[1][X_])
	    + get_curly_H(alpha[2][X_], beta[2][X_], eta[2][X_], etap[2][X_]))
	    /6.0;

  return curly_H;
}


void prt_sigma(void)
{
  long int  i;
  double    code = 0.0;
  FILE      *outf;

  outf = file_write("../out/sigma.out");

  fprintf(outf, "#  name     s   sqrt(sx)   sqrt(sx')  sqrt(sy)  sqrt(sy')\n");
  fprintf(outf, "#           [m]   [mm]       [mrad]     [mm]     [mrad]\n");
  fprintf(outf, "#\n");

  for (i = 0; i <= globval.Cell_nLoc; i++) {
    switch (Lattice.Cell[i].Elem.Pkind) {
    case drift:
      code = 0.0;
      break;
    case Mpole:
      if (Lattice.Cell[i].Elem.M->Pirho != 0)
	code = 0.5;
      else if (Lattice.Cell[i].Elem.M->PBpar[Quad+HOMmax] != 0)
	code = sgn(Lattice.Cell[i].Elem.M->PBpar[Quad+HOMmax]);
      else if (Lattice.Cell[i].Elem.M->PBpar[Sext+HOMmax] != 0)
	code = 1.5*sgn(Lattice.Cell[i].Elem.M->PBpar[Sext+HOMmax]);
      else if (Lattice.Cell[i].Fnum == globval.bpm)
	code = 2.0;
      else
	code = 0.0;
      break;
    default:
      code = 0.0;
      break;
    }
    fprintf(outf, "%4ld %.*s %6.2f %4.1f %9.3e %9.3e %9.3e %9.3e\n",
            i, SymbolLength, Lattice.Cell[i].Elem.PName, Lattice.Cell[i].S, code,
            1e3*sqrt(Lattice.Cell[i].sigma[x_][x_]),
	    1e3*sqrt(fabs(Lattice.Cell[i].sigma[x_][px_])),
	    1e3*sqrt(Lattice.Cell[i].sigma[y_][y_]),
	    1e3*sqrt(fabs(Lattice.Cell[i].sigma[y_][py_])));
  }

  fclose(outf);
}


void recalc_S(void)
{
  long int  k;
  double    S_tot;

  S_tot = 0.0;
  for (k = 0; k <= globval.Cell_nLoc; k++) {
    S_tot += Lattice.Cell[k].Elem.PL; Lattice.Cell[k].S = S_tot;
  }
}


bool getcod(double dP, long &lastpos)
{
  return GetCOD(globval.CODimax, globval.CODeps, dP, lastpos);
}


void get_twiss3(long int loc,
		Vector2 alpha[], Vector2 beta[], Vector2 nu[],
		Vector2 eta[], Vector2 etap[])
{
  /* Get Twiss functions at magnet entrance-, center-, and exit. */

  long int  j, k;
  Vector2   dnu;
  Matrix    Ascr0, Ascr1;

  // Lattice functions at the magnet entrance
  for (k = 0; k <= 1; k++) {
    alpha[0][k] = Lattice.Cell[loc-1].Alpha[k]; beta[0][k] = Lattice.Cell[loc-1].Beta[k];
    nu[0][k] = Lattice.Cell[loc-1].Nu[k];
    eta[0][k] = Lattice.Cell[loc-1].Eta[k]; etap[0][k] = Lattice.Cell[loc-1].Etap[k];
  }

  UnitMat(5, Ascr0);
  for (k = 0; k <= 1; k++) {
    Ascr0[2*k][2*k] = sqrt(Lattice.Cell[loc-1].Beta[k]); Ascr0[2*k][2*k+1] = 0.0;
    Ascr0[2*k+1][2*k] = -Lattice.Cell[loc-1].Alpha[k]/Ascr0[2*k][2*k];
    Ascr0[2*k+1][2*k+1] = 1/Ascr0[2*k][2*k];
  }
  Ascr0[0][4] = eta[0][X_]; Ascr0[1][4] = etap[0][X_];
  Ascr0[2][4] = eta[1][Y_]; Ascr0[3][4] = etap[1][Y_];
  CopyMat(5, Ascr0, Ascr1); MulLMat(5, Lattice.Cell[loc].Elem.M->AU55, Ascr1);
  dnu[0] = (atan(Ascr1[0][1]/Ascr1[0][0])-atan(Ascr0[0][1]/Ascr0[0][0]))
           /(2.0*M_PI);
  dnu[1] = (atan(Ascr1[2][3]/Ascr1[2][2])-atan(Ascr0[2][3]/Ascr0[2][2]))
           /(2.0*M_PI);

  // Lattice functions at the magnet center
  for (k = 0; k <= 1; k++) {
    alpha[1][k] = -Ascr1[2*k][2*k]*Ascr1[2*k+1][2*k]
                  - Ascr1[2*k][2*k+1]*Ascr1[2*k+1][2*k+1];
    beta[1][k] = pow(Ascr1[2*k][2*k], 2.0) + pow(Ascr1[2*k][2*k+1], 2.0);
    nu[1][k] = nu[0][k] + dnu[k];
    j = 2*k + 1; eta[1][k] = Ascr1[j-1][4]; etap[1][k] = Ascr1[j][4];
  }

  CopyMat(5, Ascr1, Ascr0); MulLMat(5, Lattice.Cell[loc].Elem.M->AD55, Ascr1);
  dnu[0] = (atan(Ascr1[0][1]/Ascr1[0][0])-atan(Ascr0[0][1]/Ascr0[0][0]))
           /(2.0*M_PI);
  dnu[1] = (atan(Ascr1[2][3]/Ascr1[2][2])-atan(Ascr0[2][3]/Ascr0[2][2]))
           /(2.0*M_PI);

  // Lattice functions at the magnet exit
  for (k = 0; k <= 1; k++) {
    alpha[2][k] = -Ascr1[2*k][2*k]*Ascr1[2*k+1][2*k]
                  - Ascr1[2*k][2*k+1]*Ascr1[2*k+1][2*k+1];
    beta[2][k] = pow(Ascr1[2*k][2*k], 2.0) + pow(Ascr1[2*k][2*k+1], 2.0);
    nu[2][k] = nu[1][k] + dnu[k];
    j = 2*k + 1; eta[2][k] = Ascr1[j-1][4]; etap[2][k] = Ascr1[j][4];
  }
}


void getabn(Vector2 &alpha, Vector2 &beta, Vector2 &nu)
{
  Vector2 gamma;
  Cell_GetABGN(globval.OneTurnMat, alpha, beta, gamma, nu);
}


void TraceABN(long i0, long i1, const Vector2 &alpha, const Vector2 &beta,
              const Vector2 &eta, const Vector2 &etap, const double dP)
{
  long          i, j;
  double        sb;
  ss_vect<tps>  Ascr;

  UnitMat(6, globval.Ascr);
  for (i = 1; i <= 2; i++) {
    sb = sqrt(beta[i-1]); j = i*2 - 1;
    globval.Ascr[j-1][j-1] = sb;               globval.Ascr[j-1][j] = 0.0;
    globval.Ascr[j][j - 1] = -(alpha[i-1]/sb); globval.Ascr[j][j] = 1/sb;
  }
  globval.Ascr[0][4] = eta[0]; globval.Ascr[1][4] = etap[0];
  globval.Ascr[2][4] = eta[1]; globval.Ascr[3][4] = etap[1];

  for (i = 0; i < 6; i++)
    globval.CODvect[i] = 0.0;
  globval.CODvect[4] = dP;

  for (i = 0; i <= 5; i++) {
    Ascr[i] = tps(globval.CODvect[i]);
    for (j = 0; j <= 5; j++)
      Ascr[i] += globval.Ascr[i][j]*tps(0.0, j+1);
    Cell_Twiss(i0, i1, Ascr, false, false, dP);
  }

}


void FitTune(long qf, long qd, double nux, double nuy)
{
  long      i;
  iVector2  nq = {0,0};
  Vector2   nu = {0.0, 0.0};
  fitvect   qfbuf, qdbuf;

  /* Get elements for the first quadrupole family */
  nq[X_] = GetnKid(qf);
  for (i = 1; i <= nq[X_]; i++)
    qfbuf[i-1] = Elem_GetPos(qf, i);

  /* Get elements for the second quadrupole family */
  nq[Y_] = GetnKid(qd);
  for (i = 1; i <= nq[Y_]; i++)
    qdbuf[i - 1] = Elem_GetPos(qd, i);

  nu[X_] = nux; nu[Y_] = nuy;

  /* fit tunes */
  Ring_Fittune(nu, nueps, nq, qfbuf, qdbuf, nudkL, nuimax);
}


void FitChrom(long sf, long sd, double ksix, double ksiy)
{
  long      i;
  iVector2  ns = {0,0};
  fitvect   sfbuf, sdbuf;
  Vector2   ksi ={0.0, 0.0};

  /* Get elements for the first sextupole family */
  ns[X_] = GetnKid(sf);
  for (i = 1; i <= ns[X_]; i++)
    sfbuf[i-1] = Elem_GetPos(sf, i);

  /* Get elements for the second sextupole family */
  ns[Y_] = GetnKid(sd);
  for (i = 1; i <= ns[Y_]; i++)
    sdbuf[i-1] = Elem_GetPos(sd, i);

  ksi[X_] = ksix; ksi[Y_] = ksiy;

  /* Fit chromaticities */
  /*    Ring_Fitchrom(ksi, ksieps, ns, sfbuf, sdbuf, 1.0, 1);*/
  Ring_Fitchrom(ksi, ksieps, ns, sfbuf, sdbuf, ksidkpL, ksiimax);
}


void FitDisp(long q, long pos, double eta)
{
  long     i, nq;
  fitvect  qbuf;

  /* Get elements for the quadrupole family */
  nq = GetnKid(q);
  for (i = 1; i <= nq; i++)
    qbuf[i-1] = Elem_GetPos(q, i);

  Ring_FitDisp(pos, eta, dispeps, nq, qbuf, dispdkL, dispimax);
}


#define nfloq           4

void getfloqs(psVector &x)
{
  // Transform to Floquet space
  LinTrans(nfloq, globval.Ascrinv, x);
}

#undef nfloq


#define ntrack          4

// 4D tracking in normal or Floquet space over nmax turns

void track(const char *file_name,
	   double ic1, double ic2, double ic3, double ic4, double dp,
	   long int nmax, long int &lastn, long int &lastpos, int floqs,
	   double f_rf)
{
  /* Single particle tracking around closed orbit:

          Output                floqs

        Phase Space               0     [x, px, y, py, delta, ct]  
	Floquet Space             1     [x^, px^, y^, py^, delta, ct]
	Action-Angle Variables    2     [2Jx, phx, 2Jy, phiy, delta, ct]

  */
  long int   i;
  double     twoJx, twoJy, phix, phiy, scl_1 = 1.0, scl_2 = 1.0;
  psVector     x0, x1, x2, xf;
  FILE       *outf;

  bool  prt = false;

  if (floqs == 0) {
    scl_1 = 1e3; scl_2 = 1e3;
    x0[x_] = ic1; x0[px_] = ic2; x0[y_] = ic3; x0[py_] = ic4;
  } else if (floqs == 1) {
    scl_1 = 1.0; scl_2 = 1.0;
    x0[x_] = ic1; x0[px_] = ic2; x0[y_] = ic3; x0[py_] = ic4;
    LinTrans(4, globval.Ascr, x0);
  } else if (floqs == 2) {
    scl_1 = 1e6; scl_2 = 1.0;
    x0[x_] = sqrt(ic1)*cos(ic2); x0[px_] = -sqrt(ic1)*sin(ic2);
    x0[y_] = sqrt(ic3)*cos(ic4); x0[py_] = -sqrt(ic3)*sin(ic4);
    LinTrans(4, globval.Ascr, x0);
  }

  outf = file_write(file_name);
  fprintf(outf, "# Tracking with TRACY");
  getcod(dp, lastpos);
  if (floqs == 0)
    fprintf(outf, "\n");
  else if (floqs == 1) {
    Ring_GetTwiss(false, dp);
    fprintf(outf, "# (Floquet space)\n");
  } else if (floqs == 2) {
    Ring_GetTwiss(false, dp);
    fprintf(outf, "# (Action-Angle variables)\n");
  }
  fprintf(outf, "#\n");
  fprintf(outf, "#%3d%6ld% .1E% .1E% .1E% .1E% 7.5f% 7.5f\n",
	  1, nmax, 1e0, 1e0, 0e0, 0e0, globval.TotalTune[0],
	  globval.TotalTune[1]);
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
	    1e3*globval.CODvect[ct_]);
  } else {
    fprintf(outf, "         [%%]           [deg]\n");
    fprintf(outf, "#\n");
    fprintf(outf, "%4d %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e\n",
	    0, scl_1*ic1, scl_2*ic2, scl_1*ic3, scl_2*ic4, 1e2*dp, 
	    2.0*f_rf*180.0*globval.CODvect[ct_]/c0);
  }
  x2[x_] = x0[x_] + globval.CODvect[x_];
  x2[px_] = x0[px_] + globval.CODvect[px_];
  x2[y_] = x0[y_] + globval.CODvect[y_];
  x2[py_] = x0[py_] + globval.CODvect[py_];
  if (globval.Cavity_on) {
    x2[delta_] = dp + globval.CODvect[delta_]; x2[ct_] = globval.CODvect[ct_];
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

    Cell_Pass(0, globval.Cell_nLoc, x2, lastpos);

    for (i = x_; i <= py_; i++)
      xf[i] = x2[i] - globval.CODvect[i];

    for (i = delta_; i <= ct_; i++)
      if (globval.Cavity_on && (i != ct_)) 
	xf[i] = x2[i] - globval.CODvect[i];
      else
	xf[i] = x2[i];

    if (floqs == 1)
      getfloqs(xf);
    else if (floqs == 2) {
      getfloqs(xf);
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
	      1e2*xf[4], 1e3*xf[5]);
    else
      fprintf(outf,
	      "%4ld %23.16le %23.16le %23.16le %23.16le %23.16le %23.16le\n",
	      lastn, scl_1*xf[0], scl_2*xf[1], scl_1*xf[2], scl_2*xf[3],
	      1e2*xf[4], 2.0*f_rf*180.0*xf[5]/c0);
  } while ((lastn != nmax) && (lastpos == globval.Cell_nLoc));

  fclose(outf);
}

#undef ntrack


#define step            0.1
#define px              0.0
#define py              0.0
void track_(double r, struct LOC_getdynap *LINK)
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
    LinTrans(5, globval.Ascr, x);
  }
  for (i = 0; i <= 3; i++)
    x[i] += globval.CODvect[i];
  lastn = 0;
  do {
    lastn++;
    Cell_Pass(0, globval.Cell_nLoc, x, lastpos);
  } while (lastn != LINK->nturn && lastpos == globval.Cell_nLoc);
  LINK->lost = (lastn != LINK->nturn);
}
#undef step
#undef px
#undef py


/****************************************************************************/
/* void Trac(double x, double px, double y, double py, double dp, double ctau,
          long nmax, long pos, long &lastn, long &lastpos, FILE *outf1)

   Purpose:
      Single particle tracking w/ respect to the chamber centrum

   Input:
      x, px, y, py 4 transverses coordinates
      ctau         c*tau
      dp           energy offset
      nmax         number of turns
      pos          starting position for tracking
      aperture     global physical aperture

   Output:
      lastn       last n (should be nmax if  not lost)
      lastpos     last position in the ring

   Return:
       none

   Global variables:
       globval

   specific functions:
       Cell_Pass

   Comments:
       Absolute TRACKING with respect to the center of the vacuum vessel
       BUG: last printout is wrong because not at pos but at the end of
            the ring
       26/04/03 print output for phase space is for position pos now
       01/12/03 tracking from 0 to pos -1L instead of 0 to pos
       (wrong observation point)

****************************************************************************/
void Trac(double x, double px, double y, double py, double dp, double ctau,
          long nmax, long pos, long &lastn, long &lastpos, FILE *outf1)
{
  psVector x1;     /* tracking coordinates */

  /* Compute closed orbit : usefull if insertion devices */

  x1[0] = x; x1[1] = px;
  x1[2] = y; x1[3] = py;
  x1[4] =dp; x1[5] = ctau;

  lastn = 0L;

  (lastpos)=pos;
  if(trace) fprintf(outf1, "\n");
  fprintf(outf1, "%6ld %+10.5e %+10.5e %+10.5e %+10.5e %+10.5e %+10.5e \n", lastn,
                 x1[0], x1[1], x1[2], x1[3], x1[4], x1[5]);
  Cell_Pass(pos -1L, globval.Cell_nLoc, x1, lastpos);

  if (lastpos == globval.Cell_nLoc)
  do {
    (lastn)++;
    Cell_Pass(0L,pos-1L, x1, lastpos);
    if(!trace) {
      fprintf(outf1, "%6ld %+10.5e %+10.5e %+10.5e %+10.5e %+10.5e %+10.5e \n",
	      lastn, x1[0], x1[1], x1[2], x1[3], x1[4], x1[5]);
    }
    if (lastpos == pos-1L) Cell_Pass(pos-1L,globval.Cell_nLoc, x1, lastpos);
  }
  while (((lastn) < nmax) && ((lastpos) == globval.Cell_nLoc));

  if (lastpos == globval.Cell_nLoc) Cell_Pass(0L,pos, x1, lastpos);

  if (lastpos != pos) {
    printf("Trac: Particle lost \n");
    fprintf(stdout, "turn:%6ld plane: %1d"
	    " %+10.5g %+10.5g %+10.5g %+10.5g %+10.5g %+10.5g \n", 
	    lastn, status.lossplane, x1[0], x1[1], x1[2], x1[3], x1[4], x1[5]);
  }
}

/****************************************************************************/
/*bool chk_if_lost(double x0, double y0, double delta,
		 long int nturn, bool floqs)
		 
   Purpose:
       Binary search for dynamical aperture in Floquet space.

   Input:
       none

   Output:
       none

   Return:
       none

   Global variables:
       px_0, py_0

   Specific functions:
       chk_if_lost

   Comments:
       none

****************************************************************************/

#define nfloq     4
bool chk_if_lost(double x0, double y0, double delta,
		 long int nturn, bool floqs)
{
  long int  i, lastn, lastpos;
  psVector    x;

  bool  prt = false;

  x[x_] = x0; x[px_] = px_0; x[y_] = y0; x[py_] = py_0;
  x[delta_] = delta; x[ct_] = 0.0;
  if (floqs)
    // transform to phase space
    LinTrans(nfloq, globval.Ascr, x);  
  for (i = 0; i <= 3; i++)  
    x[i] += globval.CODvect[i];

  lastn = 0;
  if (prt) {
    printf("\n");
    printf("chk_if_lost:\n");
    printf("%4ld %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
	   lastn, 1e3*x[x_], 1e3*x[px_], 1e3*x[y_], 1e3*x[py_],
	   1e2*x[delta_], 1e3*x[ct_]);
  }
  do {
    lastn++;
    Cell_Pass(0, globval.Cell_nLoc, x, lastpos);
    if (prt)
      printf("%4ld %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
	     lastn, 1e3*x[x_], 1e3*x[px_], 1e3*x[y_], 1e3*x[py_],
	     1e2*x[delta_], 1e3*x[ct_]);
  } while ((lastn != nturn) && (lastpos == globval.Cell_nLoc));
  return(lastn != nturn);
}
#undef nfloq

/****************************************************************************/
/* void getdynap(double *r, double phi, double delta, double eps,
	      int nturn, bool floqs)

   Purpose:
       Binary search for dynamical aperture in Floquet space.


   Input:
       none

   Output:
       none

   Return:
       none

   Global variables:
       none

   Specific functions:
       chk_if_lost

   Comments:
       none

****************************************************************************/
void getdynap(double &r, double phi, double delta, double eps,
	      int nturn, bool floqs)
{
  /* Determine dynamical aperture by binary search. */
  double  rmin = 0.0, rmax = r;

  const bool    prt   = false;
  const double  r_reset = 1e-3, r0 = 10e-3;

  if (prt) printf("\n");

  while (!chk_if_lost(rmax*cos(phi), rmax*sin(phi), delta, nturn, floqs)) {
    if (rmax < r_reset) rmax = r0;
    rmax *= 2.0;
  }
  while (rmax-rmin >= eps) {
    r = rmin + (rmax-rmin)/2.0;
    if (prt) printf("getdynap: %6.3f %6.3f %6.3f\n",
		    1e3*rmin, 1e3*rmax, 1e3*r);
    if (! chk_if_lost(r*cos(phi), r*sin(phi), delta, nturn, floqs) )
      rmin = r;
    else
      rmax = r;
    if (rmin > rmax) {
      printf("getdynap: rmin > rmax\n");
      exit_(0);
    }

  }
  r = rmin;
}



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
void getcsAscr(void)
{
  long i, j;
  double phi;
  Matrix R;

  UnitMat(6, R);
  for (i = 1; i <= 2; i++) {
    phi = -atan2(globval.Ascr[i * 2 - 2][i * 2 - 1], globval.Ascr[i * 2 - 2]
		[i * 2 - 2]);
    R[i * 2 - 2][i * 2 - 2] = cos(phi);
    R[i * 2 - 1][i * 2 - 1] = R[i * 2 - 2][i * 2 - 2];
    R[i * 2 - 2][i * 2 - 1] = sin(phi);
    R[i * 2 - 1][i * 2 - 2] = -R[i * 2 - 2][i * 2 - 1];
  }
  MulRMat(6, globval.Ascr, R);
  for (i = 1; i <= 2; i++) {
    if (globval.Ascr[i * 2 - 2][i * 2 - 2] < 0.0) {
      for (j = 0; j <= 5; j++) {
	globval.Ascr[j][i * 2 - 2] = -globval.Ascr[j][i * 2 - 2];
	globval.Ascr[j][i * 2 - 1] = -globval.Ascr[j][i * 2 - 1];
      }
    }
  }
  if (!InvMat(6, globval.Ascrinv))
    printf("  *** Ascr is singular\n");
}


/****************************************************************************/
/* void dynap(double r, double delta, double eps, int npoint, int nturn,
	   double x[], double y[], bool floqs, bool print)

   Purpose:
       Determine the dynamical aperture by tracking using polar coordinates,
       and sampling in phase.
       Assumes mid-plane symmetry

   Input:
       r initial guess
       delta off momentum energy
       eps precision for binary search
       npoint sample number for phase coordinate
       nturn number of turn for computing da
       floqs true means Floquet space
       print true means Print out to the screen

   Output:
       x[] horizontal dynamics aperture
       y[] vertical dynamics aperture

   Return:
       none

   Global variables:
       none

   Specific functions:
       getdynap

   Comments:
       none

****************************************************************************/
void dynap(FILE *fp, double r, const double delta,
	   const double eps, const int npoint, const int nturn,
       double x[], double y[], const bool floqs, const bool cod, const bool print)

{
  /* Determine the dynamical aperture by tracking.
     Assumes mid-plane symmetry.                    */

  long int  i, lastpos;
  double    phi, x_min, x_max, y_min, y_max;

  if (cod)
    getcod(delta, lastpos);
  else
    globval.CODvect.zero();
  if (floqs) {
    Ring_GetTwiss(false, delta);
    if (print) {
      printf("\n");
      printf("Dynamical Aperture (Floquet space):\n");
      printf("     x^         y^\n");
      printf("\n");
    }
    fprintf(fp, "# Dynamical Aperture (Floquet space):\n");
    fprintf(fp, "#      x^         y^\n");
    fprintf(fp, "#\n");
  } else {
    if (print) {
      printf("\n");
      printf("Dynamical Aperture:\n");
      printf("     x      y\n");
      printf("    [mm]   [mm]\n");
      printf("\n");
    }
    fprintf(fp, "# Dynamical Aperture:\n");
    fprintf(fp, "#    x      y\n");
    fprintf(fp, "#   [mm]   [mm]\n");
    fprintf(fp, "#\n");
  }

  x_min = 0.0; x_max = 0.0; y_min = 0.0; y_max = 0.0;
  for (i = 0; i < npoint; i++) {
    phi = i*M_PI/(npoint-1);
    if (i == 0) 
      phi = 1e-3;
    else if (i == npoint-1)
      phi -= 1e-3;
    getdynap(r, phi, delta, eps, nturn, floqs);
    x[i] = r*cos(phi); y[i] = r*sin(phi);
    x_min = min(x[i], x_min); x_max = max(x[i], x_max);
    y_min = min(y[i], y_min); y_max = max(y[i], y_max);
    if (!floqs) {
      if (print)
        printf("  %6.2f %6.2f\n", 1e3*x[i], 1e3*y[i]);
      fprintf(fp, "  %6.2f %6.2f\n", 1e3*x[i], 1e3*y[i]);
    } else {
      if (print)
        printf("  %10.3e %10.3e\n", x[i], y[i]);
      fprintf(fp, "  %10.3e %10.3e\n", x[i], y[i]);
    }
    fflush(fp);
  }

  if (print) {
    printf("\n");
    printf("  x^: %6.2f - %5.2f y^: %6.2f - %5.2f mm\n",
	   1e3*x_min, 1e3*x_max, 1e3*y_min, 1e3*y_max);
  }
}

/****************************************************************************/
/* double get_aper(int n, double x[], double y[])

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
double get_aper(int n, double x[], double y[])
{
  int     i;
  double  A;

  A = 0.0;
  for (i = 2; i <= n; i++)
    A += x[i-2]*y[i-1] - x[i-1]*y[i-2];
  A += x[n-1]*y[0] - x[0]*y[n-1];
// x2 from mid-plane symmetry
  A = fabs(A);
//  printf("\n");
//  printf("  Dyn. Aper.: %5.1f mm^2\n", 1e6*A);
  return(A);
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
void GetPhi(long n, double *x, double *px, double *y, double *py)
{
  /* Calculates the linear phase */
  long i;

  for (i = 1; i <= n; i++) {
    x[i - 1] = GetArg(x[i - 1], px[i - 1], i * globval.TotalTune[0]);
    y[i - 1] = GetArg(y[i - 1], py[i - 1], i * globval.TotalTune[1]);
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

void SetTol(int Fnum, double dxrms, double dyrms, double drrms)
{
  int   i;
  long  k;

  for (i = 1; i <= GetnKid(Fnum); i++) {
    k = Elem_GetPos(Fnum, i);
    Lattice.Cell[k].Elem.M->PdSrms[X_] = dxrms;
    Lattice.Cell[k].Elem.M->PdSrnd[X_] = normranf();
    Lattice.Cell[k].Elem.M->PdSrms[Y_] = dyrms;
    Lattice.Cell[k].Elem.M->PdSrnd[Y_] = normranf();
    Lattice.Cell[k].Elem.M->PdTrms = drrms;
    Lattice.Cell[k].Elem.M->PdTrnd = normranf();
    Mpole_SetdS(Fnum, i); Mpole_SetdT(Fnum, i);
  }
}


void Scale_Tol(int Fnum, double dxrms, double dyrms, double drrms)
{
  int       Knum;
  long int  loc;

  for (Knum = 1; Knum <= GetnKid(Fnum); Knum++) {
    loc = Elem_GetPos(Fnum, Knum);
    Lattice.Cell[loc].Elem.M->PdSrms[X_] = dxrms; Lattice.Cell[loc].Elem.M->PdSrms[Y_] = dyrms;
    Lattice.Cell[loc].Elem.M->PdTrms    = drrms;
    Mpole_SetdS(Fnum, Knum); Mpole_SetdT(Fnum, Knum);
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
void SetaTol(int Fnum, int Knum, double dx, double dy, double dr)
{
  long int  loc;

  loc = Elem_GetPos(Fnum, Knum);
  Lattice.Cell[loc].Elem.M->PdSrms[0] = dx; Lattice.Cell[loc].Elem.M->PdSrnd[0] = 1e0;
  Lattice.Cell[loc].Elem.M->PdSrms[1] = dy; Lattice.Cell[loc].Elem.M->PdSrnd[1] = 1e0;
  Lattice.Cell[loc].Elem.M->PdTrms    = dr; Lattice.Cell[loc].Elem.M->PdTrnd    = 1e0;
  Mpole_SetdS(Fnum, Knum); Mpole_SetdT(Fnum, Knum);
}


void ini_aper(const double Dxmin, const double Dxmax, 
              const double Dymin, const double Dymax) 
{ 
  int  k; 
 
  for (k = 0; k <= globval.Cell_nLoc; k++) { 
    Lattice.Cell[k].maxampl[X_][0] = Dxmin; Lattice.Cell[k].maxampl[X_][1] = Dxmax; 
    Lattice.Cell[k].maxampl[Y_][0] = Dymin; Lattice.Cell[k].maxampl[Y_][1] = Dymax; 
  } 
} 
 
void set_aper(const int Fnum, const double Dxmin, const double Dxmax,
	      const double Dymin, const double Dymax)
{
  int       i;
  long int  loc;

  for (i = 1; i <= GetnKid(Fnum); i++) {
    loc = Elem_GetPos(Fnum, i);
    Lattice.Cell[loc].maxampl[X_][0] = Dxmin; Lattice.Cell[loc].maxampl[X_][1] = Dxmax;
    Lattice.Cell[loc].maxampl[Y_][0] = Dymin; Lattice.Cell[loc].maxampl[Y_][1] = Dymax;
  }
}


void LoadApertures(const char *ChamberFileName)
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
      if (Fnum > 0) set_aper(Fnum, Xmin, Xmax, Ymin, Ymax);
  } while (fgets(line, 128, ChamberFile ) != NULL);

  fclose(ChamberFile);
}


// Load tolerances from the file
void LoadTolerances(const char *TolFileName) 
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
	SetTol(Fnum, dx, dy, dr);
      } else {
	printf("LoadTolerances: undefined element %s\n", FamName);
	exit_(1);
      }
    }
  } while (fgets(line, 128, tolfile) != NULL);

  fclose(tolfile);
}


// Load tolerances from the file
void ScaleTolerances(const char *TolFileName, const double scl) 
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
	Scale_Tol(Fnum, scl*dx, scl*dy, scl*dr);
      } else {
	printf("ScaleTolerances: undefined element %s\n", FamName);
	exit_(1);
      }
    }
  } while (fgets(line, 128, tolfile) != NULL);
  fclose(tolfile);
}


void SetKpar(int Fnum, int Knum, int Order, double k)
{

  Lattice.Cell[Elem_GetPos(Fnum, Knum)].Elem.M->PBpar[Order+HOMmax] = k;
  Mpole_SetPB(Fnum, Knum, Order);
}


void SetL(int Fnum, int Knum, double L)
{

  Lattice.Cell[Elem_GetPos(Fnum, Knum)].Elem.PL = L;
}


void SetL(int Fnum, double L)
{
  int  i;

  for (i = 1; i <= GetnKid(Fnum); i++)
    Lattice.Cell[Elem_GetPos(Fnum, i)].Elem.PL = L;
}


void SetdKpar(int Fnum, int Knum, int Order, double dk)
{

  Lattice.Cell[Elem_GetPos(Fnum, Knum)].Elem.M->PBpar[Order+HOMmax] += dk;
  Mpole_SetPB(Fnum, Knum, Order);
}


void SetKLpar(int Fnum, int Knum, int Order, double kL)
{
  long int  loc;

  loc = Elem_GetPos(Fnum, Knum);
  if (Lattice.Cell[loc].Elem.PL != 0e0)
    Lattice.Cell[loc].Elem.M->PBpar[Order+HOMmax] = kL/Lattice.Cell[loc].Elem.PL;
  else
    Lattice.Cell[loc].Elem.M->PBpar[Order+HOMmax] = kL;
  Mpole_SetPB(Fnum, Knum, Order);
}


void SetdKLpar(int Fnum, int Knum, int Order, double dkL)
{
  long int  loc;

  loc = Elem_GetPos(Fnum, Knum);
  if (Lattice.Cell[loc].Elem.PL != 0e0)
    Lattice.Cell[loc].Elem.M->PBpar[Order + HOMmax] += dkL/Lattice.Cell[loc].Elem.PL;
  else
    Lattice.Cell[loc].Elem.M->PBpar[Order + HOMmax] += dkL;
  Mpole_SetPB(Fnum, Knum, Order);
}


void SetdKrpar(int Fnum, int Knum, int Order, double dkrel)
{
  long int  loc;

  loc = Elem_GetPos(Fnum, Knum);
  if (Order == Dip && Lattice.Cell[loc].Elem.M->Pthick == thick)
    Lattice.Cell[loc].Elem.M->PBpar[Dip+HOMmax] += dkrel*Lattice.Cell[loc].Elem.M->Pirho;
  else
    Lattice.Cell[loc].Elem.M->PBpar[Order+HOMmax]
      += dkrel*Lattice.Cell[loc].Elem.M->PBpar[Order+HOMmax];
  Mpole_SetPB(Fnum, Knum, Order);
}


void Setbn(int Fnum, int order, double bn)
{
  int i;

  for (i = 1; i <=  GetnKid(Fnum); i++)
    SetKpar(Fnum, i, order, bn);
}


void SetbnL(int Fnum, int order, double bnL)
{
  int i;

  for (i = 1; i <= GetnKid(Fnum); i++)
    SetKLpar(Fnum, i, order, bnL);
}


void Setdbn(int Fnum, int order, double dbn)
{
  int i;

  for (i = 1; i <= GetnKid(Fnum); i++)
    SetdKpar(Fnum, i, order, dbn);
}


void SetdbnL(int Fnum, int order, double dbnL)
{
  int i;

  for (i = 1; i <= GetnKid(Fnum); i++) {
    SetdKLpar(Fnum, i, order, dbnL);
  }
}


void Setbnr(int Fnum, long order, double bnr)
{
  int  i;

  for (i = 1; i <= GetnKid(Fnum); i++)
    SetdKrpar(Fnum, i, order, bnr);
}


void SetbnL_sys(int Fnum, int Order, double bnL_sys)
{
  int       Knum;
  long int  loc;

  for (Knum = 1; Knum <= GetnKid(Fnum); Knum++) {
    loc = Elem_GetPos(Fnum, Knum);
    if (Lattice.Cell[loc].Elem.PL != 0.0)
      Lattice.Cell[loc].Elem.M->PBsys[Order+HOMmax] = bnL_sys/Lattice.Cell[loc].Elem.PL;
    else
      Lattice.Cell[loc].Elem.M->PBsys[Order+HOMmax] = bnL_sys;
    Mpole_SetPB(Fnum, Knum, Order);
  }
}


void set_dbn_rel(const int type, const int n, const double dbn_rel)
{
  long int  j;
  double    dbn;

  printf("\n");
  printf("Setting Db_%d/b_%d = %6.1e for:\n", n, type, dbn_rel);
  printf("\n");
  for (j = 0; j <= globval.Cell_nLoc; j++)
    if ((Lattice.Cell[j].Elem.Pkind == Mpole) && (Lattice.Cell[j].Elem.M->n_design == type)) {
      printf("%s\n", Lattice.Cell[j].Elem.PName);
      dbn = dbn_rel*Lattice.Cell[j].Elem.M->PBpar[type+HOMmax];
      Lattice.Cell[j].Elem.M->PBrms[n+HOMmax] = dbn;
      Lattice.Cell[j].Elem.M->PBrnd[n+HOMmax] = normranf();
      Mpole_SetPB(Lattice.Cell[j].Fnum, Lattice.Cell[j].Knum, n);
    }
}


double GetKpar(int Fnum, int Knum, int Order)
{
  return (Lattice.Cell[Elem_GetPos(Fnum, Knum)].Elem.M->PBpar[Order+HOMmax]);
}


double GetL(int Fnum, int Knum)
{
  return (Lattice.Cell[Elem_GetPos(Fnum, Knum)].Elem.PL);
}


double GetKLpar(int Fnum, int Knum, int Order)
{
  long int  loc;

  loc = Elem_GetPos(Fnum, Knum);
  if (Lattice.Cell[loc].Elem.PL != 0e0)
    return (Lattice.Cell[loc].Elem.M->PBpar[Order+HOMmax]*Lattice.Cell[loc].Elem.PL);
  else
    return (Lattice.Cell[loc].Elem.M->PBpar[Order+HOMmax]);
}


void SetdKLrms(int Fnum, int Order, double dkLrms)
{
  long int  Knum, loc;

  for (Knum = 1; Knum <= GetnKid(Fnum); Knum++) {
    loc = Elem_GetPos(Fnum, Knum);
    if (Lattice.Cell[loc].Elem.PL != 0e0)
      Lattice.Cell[loc].Elem.M->PBrms[Order+HOMmax] = dkLrms/Lattice.Cell[loc].Elem.PL;
    else
      Lattice.Cell[loc].Elem.M->PBrms[Order+HOMmax] = dkLrms;
    Lattice.Cell[loc].Elem.M->PBrnd[Order+HOMmax] = normranf();
    Mpole_SetPB(Fnum, Knum, Order);
  }
}


void Setdkrrms(int Fnum, int Order, double dkrrms)
{
  long int  Knum, loc;

  for (Knum = 1; Knum <= GetnKid(Fnum); Knum++) {
    loc = Elem_GetPos(Fnum, Knum);
    if (Order == Dip && Lattice.Cell[loc].Elem.M->Pthick == thick)
      Lattice.Cell[loc].Elem.M->PBrms[Dip+HOMmax] = dkrrms*Lattice.Cell[loc].Elem.M->Pirho;
    else
      Lattice.Cell[loc].Elem.M->PBrms[Order+HOMmax]
	= dkrrms*Lattice.Cell[loc].Elem.M->PBpar[Order+HOMmax];
    Lattice.Cell[loc].Elem.M->PBrnd[Order+HOMmax] = normranf();
    Mpole_SetPB(Fnum, Knum, Order);
  }
}


void SetKL(int Fnum, int Order)
{
  long int  Knum;

  for (Knum = 1; Knum <= GetnKid(Fnum); Knum++)
    Mpole_SetPB(Fnum, Knum, Order);
}


void set_dx(const int type, const double sigma_x, const double sigma_y)
{
  long int  j;

  printf("\n");
  printf("Setting sigma_x,y = (%6.1e, %6.1e) for b_%d:\n",
	 sigma_x, sigma_y, type);
  printf("\n");
  for (j = 0; j <= globval.Cell_nLoc; j++)
    if ((Lattice.Cell[j].Elem.Pkind == Mpole) && (Lattice.Cell[j].Elem.M->n_design == type)) {
      printf("%s\n", Lattice.Cell[j].Elem.PName);
      Lattice.Cell[j].Elem.M->PdSrms[X_] = sigma_x;
      Lattice.Cell[j].Elem.M->PdSrms[Y_] = sigma_y;
      Lattice.Cell[j].Elem.M->PdSrnd[X_] = normranf();
      Lattice.Cell[j].Elem.M->PdSrnd[Y_] = normranf();
      Mpole_SetdS(Lattice.Cell[j].Fnum, Lattice.Cell[j].Knum);
    }
}


void SetBpmdS(int Fnum, double dxrms, double dyrms)
{
  long int  Knum, loc;

  for (Knum = 1; Knum <= GetnKid(Fnum); Knum++) {
    loc = Elem_GetPos(Fnum, Knum);
    Lattice.Cell[loc].dS[X_] = normranf()*dxrms; Lattice.Cell[loc].dS[Y_] = normranf()*dyrms;
  }
}


/****************************************/
/* Routines for closed orbit correction */
/****************************************/


/****************************************************************************/
void codstat(double *mean, double *sigma, double *xmax, long lastpos, bool all)
{
  long i, n, loc;
  int      j;
  Vector2  sum, sum2;

  for (j = 0; j < 2; j++) {
    sum[j] = 0e0; sum2[j] = 0e0; xmax[j] = 0e0;
  }

  n = 0;
  if (all) {
    for (i = 0; i < lastpos; i++) {
      n++;
      for (j = 0; j < 2; j++) {
	sum[j] += Lattice.Cell[i].BeamPos[j*2];
	sum2[j] += sqr(Lattice.Cell[i].BeamPos[j*2]);
	xmax[j] = max(xmax[j], fabs(Lattice.Cell[i].BeamPos[j*2]));
      }
    }
  } else {
    for (i = 1; i <= n_bpm_[X_]; i++) {
      n++;
      for (j = 0; j < 2; j++) {
	loc = bpms_[j][i];
	sum[j] += Lattice.Cell[loc].BeamPos[j*2];
	sum2[j] += sqr(Lattice.Cell[loc].BeamPos[j*2]);
	xmax[j] = max(xmax[j], fabs(Lattice.Cell[loc].BeamPos[j*2]));
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
void CodStatBpm(double *mean, double *sigma, double *xmax, long lastpos,
                long bpmdis[])
{
  long     i, j, m, n;
  Vector2  sum, sum2;
  double   TEMP;

  m= n= 0;
  for (j = 0; j <= 1; j++) {
    sum[j] = 0.0; sum2[j] = 0.0; xmax[j] = 0.0;
  }

  for (i = 0; i <= lastpos; i++) {
    if (Lattice.Cell[i].Fnum == globval.bpm) {
      if (! bpmdis[m]) {
	for (j = 1; j <= 2; j++) {
	  sum[j - 1] += Lattice.Cell[i].BeamPos[j * 2 - 2];
	  TEMP = Lattice.Cell[i].BeamPos[j * 2 - 2];
	  sum2[j - 1] += TEMP * TEMP;
	  xmax[j - 1] = max(xmax[j - 1], fabs(Lattice.Cell[i].BeamPos[j * 2 - 2]));
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


void ChamberOff(void)
{
  int i;

  for (i = 0; i <= globval.Cell_nLoc; i++) {
    Lattice.Cell[i].maxampl[X_][0] = -max_ampl; Lattice.Cell[i].maxampl[X_][1] = max_ampl;
    Lattice.Cell[i].maxampl[Y_][0] = -max_ampl; Lattice.Cell[i].maxampl[Y_][1] = max_ampl;
  }
  status.chambre = false;
}


void PrintCh(void)
{
  long       i = 0;
  struct tm  *newtime;
  FILE       *f;

  const  char  *fic    = "chambre.out";

  newtime = GetTime();

  f = file_write(fic);
  fprintf(f, "# TRACY II v.2.6 -- %s -- %s \n", fic, asctime2(newtime));
  fprintf(f, "#    name                s      -xch     +xch     zch\n");
  fprintf(f, "#                               [mm]     [mm]     [mm]\n");
  fprintf(f, "#\n");

  for (i = 0; i <= globval.Cell_nLoc; i++)
    fprintf(f, "%4ld %15s  %6.2f  %7.3f  %7.3f  %7.3f\n",
	    i, Lattice.Cell[i].Elem.PName, Lattice.Cell[i].S,
	    Lattice.Cell[i].maxampl[X_][0]*1E3, Lattice.Cell[i].maxampl[X_][1]*1E3,
	    Lattice.Cell[i].maxampl[Y_][1]*1E3);

  fclose(f);
}


/** function from soleilcommon.c **/

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
#define nterm  2
void GetChromTrac(long Nb, long Nbtour, double emax, double *xix, double *xiz)
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

  Trac_Simple(x, xp, z, zp, emax, 0.0, Nbtour, Tab, &status);
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

  Trac_Simple(x, xp, z, zp, -emax, 0.0, Nbtour, Tab, &status);
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
void GetTuneTrac(long Nbtour, double emax, double *nux, double *nuz)
{
  double Tab[6][NTURN], fx[nterm], fz[nterm];
  int nb_freq[2];
  bool status;

  double x = 1e-6, xp = 0.0, z = 1e-6, zp = 0.0;

  Trac_Simple(x, xp, z, zp, emax, 0.0, Nbtour, Tab, &status);
  Get_NAFF(nterm, Nbtour, Tab, fx, fz, nb_freq);

  *nux = (fabs (fx[0]) > 1e-8 ? fx[0] : fx[1]);
  *nuz = fz[0];
}
#undef nterm

/****************************************************************************/
/* void ttwiss(double *alpha, double *beta, double *eta, double *etap, double dP)

   Purpose:
      Calculate Twiss functions for transport line

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
       redundant with TransTwiss

****************************************************************************/
void ttwiss(const Vector2 &alpha, const Vector2 &beta,
            const Vector2 &eta, const Vector2 &etap, const double dP)
{
  TraceABN(0, globval.Cell_nLoc, alpha, beta, eta, etap, dP);
}

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
void findcodS(double dP)
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
    
  if (globval.Cavity_on){
      dim = 6;   /* 6D tracking*/
    fprintf(stdout,"Error looking for cod in 6D\n");
    exit_(1);
  }
    else{
      dim = 4; /* 4D tracking */
      vcod[1] = Lattice.Cell[0].Eta[0]*dP; vcod[2] = Lattice.Cell[0].Etap[0]*dP;
      vcod[3] = Lattice.Cell[0].Eta[1]*dP; vcod[4] = Lattice.Cell[0].Etap[1]*dP;
  }
  
  Newton_RaphsonS(ntrial, vcod, dim, tolx);

  if (status.codflag == false)
    fprintf(stdout, "Error No COD found\n");
  if (trace) {
    for (k = 1; k <= 6; k++)
      x0[k-1] = vcod[k];
    fprintf(stdout, "Before cod % .5e % .5e % .5e % .5e % .5e % .5e \n",
	    x0[0], x0[1], x0[2], x0[3], x0[4], x0[5]);
    Cell_Pass(0, globval.Cell_nLoc, x0, lastpos);
    fprintf(stdout, "After  cod % .5e % .5e % .5e % .5e % .5e % .5e \n",
	    x0[0], x0[1], x0[2], x0[3], x0[4], x0[5]);
    Cell_Pass(0, globval.Cell_nLoc, x0, lastpos);
  }
  free_dvector(vcod,1,6);
}

/****************************************************************************/
/* void findcod(double dP)

   Purpose: 
       Search for the closed orbit using a numerical method
       Algo: Newton_Raphson method
             Quadratic convergence
             May need a guess starting point
             Simple precision algorithm
      4D
      Starting point: linear closed orbit

      6D
      Starting point: zero
        if radiation on : x[5] is the synchroneous phase (equilibrium RF phase)
                     off: x[5] is zero

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
       Method introduced because of bad convergence of da for ID
       using RADIA maps

****************************************************************************/
void findcod(double dP)
{
  psVector        vcod;
  const int     ntrial = 40;  // maximum number of trials for closed orbit
  const double  tolx = 1e-10;  // numerical precision
  int           k, dim = 0;
  long          lastpos;

  // initializations
  for (k = 0; k <= 5; k++)
    vcod[k] = 0.0;  
    
  if (globval.Cavity_on){
    fprintf(stdout,"warning looking for cod in 6D\n");
    dim = 6;
  } else{ // starting point linear closed orbit
    dim = 4;
    vcod[0] = Lattice.Cell[0].Eta[0]*dP; vcod[1] = Lattice.Cell[0].Etap[0]*dP;
    vcod[2] = Lattice.Cell[0].Eta[1]*dP; vcod[3] = Lattice.Cell[0].Etap[1]*dP;
    vcod[4] = dP;  // energy offset 
  }
  
  Newton_Raphson(dim, vcod, ntrial, tolx);

  if (status.codflag == false)
    fprintf(stdout, "Error No COD found\n");
  
  CopyVec(6, vcod, globval.CODvect); // save closed orbit at the ring entrance

  if (trace)
  {
    fprintf(stdout,
       "Before cod2 % .5e % .5e % .5e % .5e % .5e % .5e \n",
       vcod[0], vcod[1], vcod[2], vcod[3], vcod[4], vcod[5]);
    Cell_Pass(0, globval.Cell_nLoc, vcod, lastpos);
    fprintf(stdout,
       "After  cod2 % .5e % .5e % .5e % .5e % .5e % .5e \n",
       vcod[0], vcod[1], vcod[2], vcod[3], vcod[4], vcod[5]);
  }
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

void computeFandJS(double *x, int n, double **fjac, double *fvect)
{
  int     i, k;
  long    lastpos = 0L;
  psVector  x0, fx, fx1, fx2;

  const double deps = 1e-8;  //stepsize for numerical differentiation

  for (i = 1; i <= 6; i++)
    x0[i - 1] = x[i];
  
  Cell_Pass(0, globval.Cell_nLoc, x0, lastpos);

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

    Cell_Pass(0, globval.Cell_nLoc, x0, lastpos);  // tracking along the ring
    for (i = 1; i <= 6; i++)
      fx1[i - 1] = x0[i - 1];

    for (i = 1; i <= 6; i++)
      x0[i - 1] = x[i];
    x0[5] = 0.0;
    x0[k] -= deps;  // differential step in coordinate k

    Cell_Pass(0, globval.Cell_nLoc, x0, lastpos);  // tracking along the ring
    for (i = 1; i <= 6; i++)
      fx2[i - 1] = x0[i - 1];

    for (i = 1; i <= n; i++)  // symmetric difference formula
      fjac[i][k + 1] = 0.5 * (fx1[i - 1] - fx2[i - 1]) / deps;
    //~ for (i = 1; i <= n; i++) // forward difference formula
    //~ fjac[i][k + 1] = (float) ((x0[i - 1] - fx[i - 1]) / deps);  
  }
}

/****************************************************************************/
/* void computeFand(int n, float *x, float **fjac, float *fvect)

   Purpose:       
       Tracks x over one turn. And computes the Jacobian matrix of the 
       transformation by numerical differentiation.
       using symmetric difference formula
       double precision algorithm

   Input:
       x vector for evaluation

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
void computeFandJ(int n, psVector &x, Matrix &fjac, psVector &fvect)
{
  int     i, k;
  long    lastpos = 0;
  psVector  x0, fx1, fx2;

  const double deps = 1e-8;  //stepsize for numerical differentiation

  CopyVec(6, x, x0);
  
  Cell_Pass(0, globval.Cell_nLoc, x0, lastpos);
  CopyVec(n, x0, fvect);
  
  // compute Jacobian matrix by numerical differentiation
  for (k = 0; k < n; k++) {
    CopyVec(6L, x, x0);
    x0[k] += deps;  // differential step in coordinate k

    Cell_Pass(0, globval.Cell_nLoc, x0, lastpos);  // tracking along the ring
      CopyVec(6L, x0, fx1);

      CopyVec(6L, x, x0);
    x0[k] -= deps;  // differential step in coordinate k

    Cell_Pass(0, globval.Cell_nLoc, x0, lastpos);  // tracking along the ring
      CopyVec(6L, x0, fx2);

    for (i = 0; i < n; i++)  // symmetric difference formula
      fjac[i][k] = 0.5 * (fx1[i] - fx2[i]) / deps;
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

void Newton_RaphsonS(int ntrial, double x[], int n, double tolx)
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
    computeFandJS(x, n, alpha, fvect);

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
      status.codflag = true;
      break;
    }
  }
  // check whever closed orbit found out
  if ((k >= ntrial) && (errx >= tolx * 100)) status.codflag = false;

  free_dmatrix(alpha,1,n,1,n); free_dvector(bet,1,n); free_dvector(fvect,1,n);
  free_ivector(indx,1,n);
}


/****************************************************************************/
/* int Newton_Raphson(int n, double x[], int ntrial, double tolx)
 
   Purpose:
       Newton_Rapson algorithm from Numerical Recipes
       double precision algorithm
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
       computeFandJ
     InvMat, LinTrans 

   Comments:
       none

****************************************************************************/
int Newton_Raphson (int n, psVector &x, int ntrial, double tolx)
{
  int k,  i;
  double  errx;
  psVector  bet, fvect;
  Matrix  alpha;

  errx = 0.0;

  for (k = 1; k <= ntrial; k++) {  // loop over number of iterations
    // supply function values at x in fvect and Jacobian matrix in fjac
    computeFandJ(n, x, alpha, fvect);

    // Jacobian - Id
    for (i = 0; i < n; i++)
      alpha[i][i] -= 1.0;
    for (i = 0; i < n; i++)
      bet[i] = x[i] - fvect[i];  // right side of linear equation
    // inverse matrix using gauss jordan method from Tracy (from NR)
    if (!InvMat((long) n,alpha))
      fprintf(stdout,"Matrix non inversible ...\n");    
    LinTrans((long) n, alpha, bet); // bet = alpha*bet
        errx = 0.0;  // check root convergence
    for (i = 0; i < n; i++)
    {    // update solution
      errx += fabs(bet[i]);
      x[i] += bet[i]; 
    }
    
    if (trace)
      fprintf(stdout,
         "%02d: cod2 % .5e % .5e % .5e % .5e % .5e % .5e  errx =% .5e\n",
         k, x[0], x[1], x[2], x[3], x[4], x[5], errx);
    if (errx <= tolx)
    {
      status.codflag = true;
        return 1;
    }
  }
  // check whever closed orbit found out
  if ((k >= ntrial) && (errx >= tolx))
  {
    status.codflag = false;
      return 1;
  }
  return 0;
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
