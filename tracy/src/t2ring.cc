/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

   t2ring.c -- Routines for closed beam lines

*/


void GetNu(Vector2 &nu, Matrix &M)
{
  /* Not assuming mid-plane symmetry, the charachteristic polynomial for a
     symplectic periodic matrix is given by

       P(lambda) = det(M-lambda*I)
                 = (lambda-lambda0)(lambda-1/lambda0)
		   (lambda-lambda1)(lambda-1/lambda1)

     It follows that

       P(1) = (2-x)(2-y),     P(-1) = (2+x)(2+y)

     where

       x = lambda0 + 1/lambda0 = 2 cos(2 pi nu_x)

     and similarly for y. Eliminating y

       x^2 + 4 b x + 4 c = 0

     where

       b = (P(1)-P(-1))/16,    c =(P(1)+P(-1))/8 - 1

     Solving for x

       x,y = -2 (b +/- sqrt(b^2-c))

     where the sign is given by

       trace(hor) > trace(ver)

     gives

       nu_x,y = arccos(x,y/2)/(2 pi)

     For mid-plane symmetry it simplies to

       nu_x = arccos((m11+m22)/2)/(2 pi)                                      */

  int    i;
  double sgn, detp, detm, b, c, tr[2], b2mc, x;
  Matrix M1;

  const int n = 4;

  CopyMat((long)n, M, M1);
  for (i = 0; i < n; i++)
    M1[i][i] -= 1e0;   
  detp = DetMat((long)n, M1);
  for (i = 0; i < n; i++)
    M1[i][i] += 2e0;
  detm = DetMat((long)n, M1);
  for (i = 0; i < n; i++)
    M1[i][i] -= 1e0;

  for (i = 0; i < 2; i++)
    tr[i] = (M1[2*i][2*i]+M1[2*i+1][2*i+1]);
  sgn = (tr[X_] > tr[Y_])? 1e0 : -1e0;

  b = (detp-detm)/16e0; c = (detp+detm)/8e0 - 1e0;
  b2mc = sqr(b) - c;

  if (b2mc < 0e0) {
    globval.stable = false; nu[0] = NAN; nu[1] = NAN;
    printf("\nGetNu: unstable\n");
    return;
  }

  for (i = 0; i < 2; i++) {
    if (i == 0)
      x = -2e0*(b-sgn*sqrt(b2mc));
    else
      x = -2e0*(b+sgn*sqrt(b2mc));
    if (fabs(x) <= 1e0) {
      nu[i] = acos(x/2e0)/(2e0*M_PI);
      if (M1[2*i][2*i+1] < 0e0) nu[i] = 1e0 - nu[i];
    } else {
      globval.stable = false; nu[i] = NAN;
      printf("\nGetNu: unstable plane %d\n", i);
      return;
    }
  }

  return;
}


void Cell_GetABGN(Matrix &M,
		  Vector2 &alpha, Vector2 &beta, Vector2 &gamma, Vector2 &nu)
{
  int    i;
  double c = 0e0, s = 0e0;

  globval.stable = true;
  for (i = 0; i < 2; i++) {
    c = (M[2*i][2*i]+M[2*i+1][2*i+1])/2e0;
    globval.stable = (fabs(c) < 1e0);
    if (globval.stable) {
      s = sqrt(1e0-sqr(c))*sgn(M[2*i][2*i+1]);
      alpha[i] = (M[2*i][2*i]-M[2*i+1][2*i+1])/(2e0*s);
      beta[i] = M[2*i][2*i+1]/s; gamma[i] = -M[2*i+1][2*i]/s;
    }
  }
  GetNu(nu, M);
  if (!globval.stable) printf("Cell_GetABGN: unstable\n");
}


void Cell_Geteta(long i0, long i1, bool ring, double dP)
{
  long int  i = 0, lastpos = 0;
  int       j = 0, k = 0;
  psVector    xref;
  psVector    codbuf[Cell_nLocMax+1];
  CellType  *cellp;

  const int n = 4;

  /* cod for the energy dP - globval.dPcommon / 2e0 */
  if (ring)
    GetCOD(globval.CODimax, globval.CODeps, dP-globval.dPcommon/2e0, lastpos);
  else { /* beam mode */
    CopyVec(n+2, globval.CODvect, xref); xref[4] = dP - globval.dPcommon/2e0;
    Cell_Pass(i0, i1, xref, lastpos);
  }

  /* Store chromatic orbit for elements i0 to i1 in codbuf */
  for (i = i0; i <= i1; i++)
    CopyVec(n+2, Cell[i].BeamPos, codbuf[i]);

  /* cod for the energy dP - globval.dPcommon / 2e0 */
  if (ring)
    GetCOD(globval.CODimax, globval.CODeps, dP+globval.dPcommon/2e0, lastpos);
  else { /* beam mode */
    CopyVec(n+2, globval.CODvect, xref); xref[4] = dP + globval.dPcommon/2e0;
    Cell_Pass(i0, i1, xref, lastpos);
  }

  for (i = i0; i <= i1; i++) {
    cellp = &Cell[i];
    for (j = 1; j <= 2; j++) {
      k = j*2 - 1;
      cellp->Eta[j-1] = (cellp->BeamPos[k-1]-codbuf[i][k-1])/globval.dPcommon;
      cellp->Etap[j-1] = (cellp->BeamPos[k]-codbuf[i][k])/globval.dPcommon;
    }
  }
}


void getprm(Matrix &Ascr, Vector2 &alpha, Vector2 &beta)
{
  int  i, j;

  for (i = 1; i <= 2; i++) {
    j = i*2 - 1;
    alpha[i-1] = -(Ascr[j-1][j-1]*Ascr[j][j-1] + Ascr[j-1][j]*Ascr[j][j]);
    beta[i-1] = sqr(Ascr[j-1][j-1]) + sqr(Ascr[j-1][j]);
  }
}


void dagetprm(ss_vect<tps> &Ascr, Vector2 &alpha, Vector2 &beta)
{
  int  i = 0, j = 0;

  for (i = 1; i <= 2; i++) {
    j = i*2 - 1;
    alpha[i-1] = -(getmat(Ascr, j, j)*getmat(Ascr, j+1, j)
                 + getmat(Ascr, j, j+1)*getmat(Ascr, j+1, j+1));
    beta[i-1] = sqr(getmat(Ascr, j, j)) + sqr(getmat(Ascr, j, j+1));
  }
}


void Cell_Twiss(long i0, long i1, ss_vect<tps> &Ascr, bool chroma, bool ring,
		double dP)
{
  long int      i = 0;
  int           j = 0, k = 0;
  Vector2       nu1, dnu;   /* absolute and relative phase advance */
  ss_vect<tps>  Ascr0, Ascr1;
  CellType      *cellp;

  const int n = 4;

  /* initialization */
  for (j = 0; j <= 1; j++) {
    nu1[j] = 0.0; dnu[j] = 0.0; 
  }

  if (globval.radiation) globval.dE = 0.0;

  cellp = &Cell[i0];
  dagetprm(Ascr, cellp->Alpha, cellp->Beta);
  memcpy(cellp->Nu, nu1, sizeof(Vector2));

  Ascr0 = Ascr;
  for (j = 0; j <= n+1; j++)
    Ascr0[j] += tps(globval.CODvect[j]);

  Ascr1 = Ascr0;
  for (i = i0; i <= i1; i++) {
    Elem_Pass(i, Ascr1); cellp = &Cell[i];
    dagetprm(Ascr1, cellp->Alpha, cellp->Beta);
    for (j = 0; j <= 1; j++) {
      k = (j+1)*2 - 1;
      dnu[j] =
	(GetAngle(getmat(Ascr1, k, k), getmat(Ascr1, k, k+1)) -
	 GetAngle(getmat(Ascr0, k, k), getmat(Ascr0, k, k+1)))/(2.0*M_PI);

      if ((cellp->Elem.PL >= 0.0) && (dnu[j] < -1e-16))
	dnu[j] += 1.0;
      else if ((cellp->Elem.PL < 0.0) && (dnu[j] > 1e-16))
	dnu[j] -= 1.0;

      nu1[j] += dnu[j];

      cellp->Nu[j] = nu1[j];
      // cellp->Eta[j] = getmat(Ascr1, k, 5)*getmat(Ascr1, 6, 6) -
      //                getmat(Ascr1, k, 6)*getmat(Ascr1, 6, 5);
      // cellp->Etap[j] = getmat(Ascr1, k+1, 5)*getmat(Ascr1, 6, 6) -
      //                 getmat(Ascr1, k+1, 6)*getmat(Ascr1, 6, 5);
      cellp->Eta[j] = getmat(Ascr1, k, 5);
      cellp->Etap[j] = getmat(Ascr1, k+1, 5);
    }
    Ascr0 = Ascr1;
  }

  if (chroma && !globval.Cavity_on) Cell_Geteta(i0, i1, ring, dP);
}


void Ring_Getchrom(double dP)
{
  long int  lastpos = 0;
  int       j;
  Vector2   alpha = {0.0, 0.0}, beta = {0.0, 0.0}, gamma = {0.0, 0.0};
  Vector2   nu = {0.0, 0.0}, nu0 = {0.0, 0.0};
  
  if (dP != 0e0)
    printf("\nRing_Getchrom: linear chromaticity for delta = %e\n", dP);
  
  GetCOD(globval.CODimax, globval.CODeps, dP-globval.dPcommon*0.5, lastpos);
  if (!status.codflag) {
    printf("\nRing_Getchrom: closed orbit finder failed\n");
    return;
  }
  Cell_GetABGN(globval.OneTurnMat, alpha, beta, gamma, nu0);
  if (!globval.stable) {
    printf("\nRing_Getchrom: unstable\n");
    return;
  }
  
  GetCOD(globval.CODimax, globval.CODeps, dP+globval.dPcommon*0.5, lastpos);
  if (!status.codflag) {
    printf("\nRing_Getchrom: closed orbit finder failed");
    return;
  }
  Cell_GetABGN(globval.OneTurnMat, alpha, beta, gamma, nu);
  if (!globval.stable) {
    printf("\n:Ring_Getchrom: unstable\n");
    return;
  }

  for (j = 0; j < 2; j++)
    globval.Chrom[j] = (nu[j]-nu0[j])/globval.dPcommon;
  
  status.chromflag = true;
}


void Ring_Twiss(bool chroma, double dP)
{
  long int      lastpos = 0;
  int           n = 0;
  Vector2       alpha={0.0, 0.0}, beta={0.0, 0.0};
  Vector2       gamma={0.0, 0.0}, nu={0.0, 0.0};
  Matrix        R;
  ss_vect<tps>  AScr;

  n = (globval.Cavity_on)? 6 : 4;

  GetCOD(globval.CODimax, globval.CODeps, dP, lastpos);

  if (!status.codflag) return;

  // Check if stable
  Cell_GetABGN(globval.OneTurnMat, alpha, beta, gamma, nu);
  if (!globval.stable) {
    printf("Ring_Twiss: unstable\n");
    return;
  }
  // Get eigenvalues and eigenvectors for the one turn transfer matrix
  GDiag(n, Cell[globval.Cell_nLoc].S, globval.Ascr, globval.Ascrinv, R,
        globval.OneTurnMat, globval.Omega, globval.Alphac);

  // putlinmat(n, globval.Ascr, AScr);
  putlinmat(6, globval.Ascr, AScr);
  if (!globval.Cavity_on) {
    // AScr[delta_] = 0.0; AScr[ct_] = 0.0;
    AScr[delta_] = tps(0e0, delta_+1); AScr[ct_] = 0e0;
  }

  Cell_Twiss(0, globval.Cell_nLoc, AScr, chroma, true, dP);

  memcpy(globval.TotalTune, Cell[globval.Cell_nLoc].Nu, sizeof(Vector2));
  status.tuneflag = true;

  if (chroma && !globval.Cavity_on) {
    Ring_Getchrom(dP); GetCOD(globval.CODimax, globval.CODeps, dP, lastpos);
  }
}


void Ring_GetTwiss(bool chroma, double dP)
{

  if (trace) printf("enter ring_gettwiss\n");
  Ring_Twiss(chroma, dP);
  globval.Alphac = globval.OneTurnMat[ct_][delta_]/Cell[globval.Cell_nLoc].S;
  if (trace) printf("exit ring_gettwiss\n");
}


/* Local variables for Ring_Fittune: */
struct LOC_Ring_Fittune
{
  jmp_buf _JL999;
};


void shiftk(long Elnum, double dk, struct LOC_Ring_Fittune *LINK)
{
  CellType   *cellp;
  elemtype   *elemp;
  MpoleType  *M;

  cellp = &Cell[Elnum]; elemp = &cellp->Elem; M = elemp->M;
  M->PBpar[Quad+HOMmax] += dk;
  Mpole_SetPB(cellp->Fnum, cellp->Knum, (long)Quad);
}


void checkifstable(struct LOC_Ring_Fittune *LINK)
{
  if (!globval.stable) {
    printf("\ncheckifstable: unstable\n");
    longjmp(LINK->_JL999, 1);
  }
}


void Ring_Fittune(Vector2 &nu, double eps, iVector2 &nq, long qf[], long qd[],
                  double dkL, long imax)
{
  struct LOC_Ring_Fittune V;

  int      i, j, k;
  Vector2  nu0, nu1;
  psVector  dkL1, dnu;
  Matrix A;

  if (setjmp(V._JL999)) return;

  if (trace)
    printf("  Tune fit, nux =%10.5f, nuy =%10.5f, eps =% .3E,"
	   " imax =%4ld, dkL = % .5E\n", nu[0], nu[1], eps, imax, dkL);
  Ring_GetTwiss(false, 0e0); checkifstable(&V);
  memcpy(nu0, globval.TotalTune, sizeof(Vector2));
  i = 0;
  do {
    i++;
    /* First vary kf then kd */
    for (j = 1; j <= 2L; j++) {
      for (k = 0; k < nq[j-1]; k++) {
        if (j == 1)
          shiftk(qf[k], dkL, &V); // new value for qf
        else
          shiftk(qd[k], dkL, &V); // new value for qd
      }
      Ring_GetTwiss(false, 0e0);
      nu1[0] = globval.TotalTune[0]; nu1[1] = globval.TotalTune[1];
//      GetCOD(globval.CODimax, globval.CODeps, 0e0, lastpos);
//      Cell_GetABGN(globval.OneTurnMat, alpha, beta, gamma, nu1);
      checkifstable(&V);
      for (k = 0; k <= 1; k++) {
        dnu[k] = nu1[k] - (long)nu1[k] - nu0[k] + (long)nu0[k];
        if (fabs(dnu[k]) > 0.5) dnu[k] = 1 - fabs(dnu[k]);
        A[k][j-1] = dnu[k]/dkL;
      }

      /* Restore strength */
      for (k = 0; k < nq[j-1]; k++) {
        if (j == 1)
          shiftk(qf[k], -dkL, &V);
        else
          shiftk(qd[k], -dkL, &V);
        }
      }

    if (!InvMat(2L, A)) {
      printf("  A is singular\n");
      return;
    }

    for (j = 0; j <= 1; j++)
      dkL1[j] = nu[j] - nu0[j];
    LinTrans(2L, A, dkL1);
    for (j = 1; j <= 2; j++) {
      for (k = 0; k < nq[j-1]; k++) {
	if (j == 1)
	  shiftk(qf[k], dkL1[j-1], &V);
	else
	  shiftk(qd[k], dkL1[j-1], &V);
      }
    }
    Ring_GetTwiss(false, 0e0); checkifstable(&V);
    memcpy(nu0, globval.TotalTune, sizeof(Vector2));
    if (trace)
      printf("  Nux = %10.6f%10.6f, Nuy = %10.6f%10.6f,"
	     " QF*L = % .5E, QD*L = % .5E @%3d\n",
	     nu0[0], nu1[0], nu0[1], nu1[1],
	     Elem_GetKval(Cell[qf[0]].Fnum, 1, (long)Quad),
	     Elem_GetKval(Cell[qd[0]].Fnum, 1, (long)Quad), i);
  } while (sqrt(sqr(nu[0]-nu0[0])+sqr(nu[1]-nu0[1])) >= eps && i != imax);
}


void shiftkp(long Elnum, double dkp)
{
  CellType  *cellp;
  elemtype  *elemp;
  MpoleType *M;

  cellp = &Cell[Elnum]; elemp = &cellp->Elem; M = elemp->M;
  M->PBpar[Sext+HOMmax] += dkp;
  Mpole_SetPB(cellp->Fnum, cellp->Knum, (long)Sext);
}


void Ring_Fitchrom(const double ksi[], const double eps, const long ns[],
		   const long sf[], const long sd[],
		   const double dkpL, const long imax)
{
  bool      rad;
  long int  lastpos;
  int       i, j, k;
  Vector2   ksi0;
  psVector    dkpL1, dksi;
  Matrix    A;

  if (trace)
    printf("  Chromaticity fit, ksix =%10.5f, ksiy =%10.5f, eps =% .3E"
	   ", imax =%4ld, dkpL =%10.5f\n", ksi[0], ksi[1], eps, imax, dkpL);

  /* Turn off radiation */
  rad = globval.radiation; globval.radiation = false;
  GetCOD(globval.CODimax, globval.CODeps, 0e0, lastpos); Ring_Getchrom(0e0);
  for (j = 0; j <= 1; j++)
    ksi0[j] = globval.Chrom[j];
  i = 0;
  do {
    i++;
    /* First vary sf then sd */
    for (j = 1; j <= 2; j++) {
      for (k = 0; k < ns[j-1]; k++) {
	if (j == 1)
	  shiftkp(sf[k], dkpL);
	else
	  shiftkp(sd[k], dkpL);
      }
      GetCOD(globval.CODimax, globval.CODeps, 0e0, lastpos); Ring_Getchrom(0e0);
      for (k = 0; k <= 1; k++) {
	dksi[k] = globval.Chrom[k] - ksi0[k];
	A[k][j-1] = dksi[k] / dkpL;
      }
      /* Restore strength */
      for (k = 0; k < ns[j-1]; k++) {
	if (j == 1)
	  shiftkp(sf[k], -dkpL);
	else
	  shiftkp(sd[k], -dkpL);
      }
    }
    if (!InvMat(2L, A)) {
      printf("  A is singular\n");
      goto _L999;
    }
    for (j = 0; j <= 1; j++)
      dkpL1[j] = ksi[j] - ksi0[j];
    LinTrans(2L, A, dkpL1);
    for (j = 1; j <= 2; j++) {
      for (k = 0; k < ns[j-1]; k++) {
	if (j == 1)
	  shiftkp(sf[k], dkpL1[j-1]);
	else
	  shiftkp(sd[k], dkpL1[j-1]);
      }
    }
    GetCOD(globval.CODimax, globval.CODeps, 0e0, lastpos); Ring_Getchrom(0e0);
    for (j = 0; j <= 1; j++)
      ksi0[j] = globval.Chrom[j];
    if (trace)
      printf("  ksix =%10.6f, ksiy =%10.6f, SF = % .5E, SD = % .5E @%3d\n",
	     ksi0[0], ksi0[1], Elem_GetKval(Cell[sf[0]].Fnum, 1, Sext),
	     Elem_GetKval(Cell[sd[0]].Fnum, 1, Sext), i);
  } while (sqrt(sqr(ksi[0]-ksi0[0])+sqr(ksi[1]-ksi0[1])) >= eps && i != imax);
_L999:
  /* Restore radiation */
  globval.radiation = rad;
}


/* Local variables for Ring_FitDisp: */
struct LOC_Ring_FitDisp
{
  jmp_buf _JL999;
};


static void shiftk_(long Elnum, double dk, struct LOC_Ring_FitDisp *LINK)
{
  CellType *cellp;
  elemtype *elemp;
  MpoleType *M;

  cellp = &Cell[Elnum];
  elemp = &cellp->Elem;
  M = elemp->M;
  M->PBpar[Quad+HOMmax] += dk;
  Mpole_SetPB(cellp->Fnum, cellp->Knum, (long)Quad);
}


void checkifstable_(struct LOC_Ring_FitDisp *LINK)
{
  if (!globval.stable) {
    printf("\ncheckifstable_: unstable\n");
    longjmp(LINK->_JL999, 1);
  }
}


void Ring_FitDisp(long pos, double eta, double eps, long nq, long q[],
                  double dkL, long imax)
{
  /*pos : integer; eta, eps : double;
                         nq : integer; var q : fitvect;
                         dkL : double; imax : integer*/
  struct LOC_Ring_FitDisp V;

  int     i, j;
  double  dkL1, Eta0, deta;
  bool    rad = false;

  if (setjmp(V._JL999)) goto _L999;

  if (trace)
    printf("  Dispersion fit, etax =%10.5f, eps =% .3E"
	   ", imax =%4ld, dkL =%10.5f\n",
	   eta, eps, imax, dkL);
  /* Turn off radiation */
  rad = globval.radiation; globval.radiation = false;
  Ring_GetTwiss(true, 0e0); checkifstable_(&V);
  Eta0 = Cell[pos].Eta[0];
  i = 0;
  while (fabs(eta-Eta0) > eps && i < imax) {
    i++;
    for (j = 0; j < nq; j++)
      shiftk_(q[j], dkL, &V);
    Ring_GetTwiss(true, 0e0); checkifstable_(&V);
    deta = Cell[pos].Eta[0] - Eta0;
    if (deta == 0.0) {
      printf("  deta is 0\n");
      goto _L999;
    }
    dkL1 = (eta-Eta0)*dkL/deta - dkL;
    for (j = 0; j < nq; j++)
      shiftk_(q[j], dkL1, &V);
    Ring_GetTwiss(true, 0e0); checkifstable_(&V);
    Eta0 = Cell[pos].Eta[0];
    if (trace)
      printf("  Dispersion = % .5E, kL =% .5E @%3d\n",
	     Eta0, Elem_GetKval(Cell[q[0]].Fnum, 1, (long)Quad), i);
  }
_L999:
  /* Restore radiation */
  globval.radiation = rad;
}
