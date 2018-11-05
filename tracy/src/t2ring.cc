/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

   t2ring.c: Functions for rings.                                             */


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

       x = (lambda0+1/lambda0)/2 = cos(2 pi nu_x)

     and similarly for y. Eliminating y

       x^2 + 4 b x + 4 c = 0

     where

       b = (P(1)-P(-1))/16,    c =(P(1)+P(-1))/8 - 1

     Solving for x

       x,y = -b +/- sqrt(b^2-c)

     where the sign is given by

       trace(hor) > trace(ver)

     gives

       nu_x,y = arccos(x,y)/(2 pi)

     For mid-plane symmetry it simplies to

       nu_x = arccos((m11+m22)/2)/(2 pi)                                      */

  int    i;
  double sgn, detp, detm, b, c, tr[2], b2mc, x;
  Matrix M1;

  const int n = 4;

  CopyMat(n, M, M1);
  for (i = 0; i < n; i++)
    M1[i][i] -= 1e0;   
  detp = DetMat(n, M1);
  for (i = 0; i < n; i++)
    M1[i][i] += 2e0;
  detm = DetMat(n, M1);
  for (i = 0; i < n; i++)
    M1[i][i] -= 1e0;

  for (i = 0; i < 2; i++)
    tr[i] = (M1[2*i][2*i]+M1[2*i+1][2*i+1]);
  sgn = (tr[X_] > tr[Y_])? 1e0 : -1e0;

  b = (detp-detm)/16e0; c = (detp+detm)/8e0 - 1e0;
  b2mc = sqr(b) - c;

  if (b2mc < 0e0) {
    globval.stable = false; nu[X_] = NAN; nu[Y_] = NAN;
    printf("\nGetNu: unstable\n");
    return;
  }

  for (i = 0; i < 2; i++) {
    if (i == 0)
      x = -b + sgn*sqrt(b2mc);
    else
      x = -b - sgn*sqrt(b2mc);
    if (fabs(x) <= 1e0) {
      nu[i] = acos(x)/(2e0*M_PI);
      if (M1[2*i][2*i+1] < 0e0) nu[i] = 1e0 - nu[i];
    } else {
      globval.stable = false; nu[i] = NAN;
      printf("\nGetNu: unstable plane %d %10.3e\n", i, x);
      return;
    }
  }

  return;
}


void Cell_GetABGN(Matrix &M,
		  Vector2 &alpha, Vector2 &beta, Vector2 &gamma, Vector2 &nu)
{
  int    k;
  double c = 0e0, s = 0e0;

  globval.stable = true;
  for (k = 0; k < 2; k++) {
    c = (M[2*k][2*k]+M[2*k+1][2*k+1])/2e0;
    globval.stable = (fabs(c) < 1e0);
    if (globval.stable) {
      s = sqrt(1e0-sqr(c))*sgn(M[2*k][2*k+1]);
      alpha[k] = (M[2*k][2*k]-M[2*k+1][2*k+1])/(2e0*s);
      beta[k] = M[2*k][2*k+1]/s; gamma[k] = -M[2*k+1][2*k]/s;
    }
  }
  GetNu(nu, M);
  if (!globval.stable) printf("Cell_GetABGN: unstable\n");
}


void Cell_Geteta(long i0, long i1, bool ring, double dP)
{
  long int i, lastpos;
  int      k;
  psVector xref;
  psVector codbuf[Cell_nLocMax+1];
  CellType *cellp;

  const int n = 4;

  if (ring)
    GetCOD(globval.CODimax, globval.CODeps, dP-globval.dPcommon/2e0, lastpos);
  else {
    CopyVec(n+2, globval.CODvect, xref); xref[4] = dP - globval.dPcommon/2e0;
    Cell_Pass(i0, i1, xref, lastpos);
  }

  for (i = i0; i <= i1; i++)
    CopyVec(n+2, Cell[i].BeamPos, codbuf[i]);

  if (ring)
    GetCOD(globval.CODimax, globval.CODeps, dP+globval.dPcommon/2e0, lastpos);
  else {
    CopyVec(n+2, globval.CODvect, xref); xref[4] = dP + globval.dPcommon/2e0;
    Cell_Pass(i0, i1, xref, lastpos);
  }

  for (i = i0; i <= i1; i++) {
    cellp = &Cell[i];
    for (k = 0; k < 2; k++) {
      cellp->Eta[k] = (cellp->BeamPos[2*k]-codbuf[i][2*k])/globval.dPcommon;
      cellp->Etap[k] =
	(cellp->BeamPos[2*k+1]-codbuf[i][2*k+1])/globval.dPcommon;
    }
  }
}


void getprm(Matrix &Ascr, Vector2 &alpha, Vector2 &beta)
{
  int k;

  for (k = 0; k < 2; k++) {
    alpha[k] =
      -(Ascr[2*k][2*k]*Ascr[2*k+1][2*k] + Ascr[2*k][2*k+1]*Ascr[2*k+1][2*k+1]);
    beta[k] = sqr(Ascr[2*k][2*k]) + sqr(Ascr[2*k][2*k+1]);
  }
}


void dagetprm(ss_vect<tps> &Ascr, Vector2 &alpha, Vector2 &beta)
{
  int k;

  for (k = 1; k <= 2; k++) {
    alpha[k-1] = -(getmat(Ascr, 2*k-1, 2*k-1)*getmat(Ascr, 2*k, 2*k-1)
                 + getmat(Ascr, 2*k-1, 2*k)*getmat(Ascr, 2*k, 2*k));
    beta[k-1] = sqr(getmat(Ascr, 2*k-1, 2*k-1)) + sqr(getmat(Ascr, 2*k-1, 2*k));
  }
}


void Cell_Twiss(long i0, long i1, ss_vect<tps> &Ascr, bool chroma, bool ring,
		double dP)
{
  long int     i;
  int          k;
  Vector2      nu1, dnu;
  ss_vect<tps> Ascr0, Ascr1;
  CellType     *cellp;

  const int n = 4;

  for (k = 0; k < 2; k++) {
    nu1[k] = 0e0; dnu[k] = 0e0; 
  }

  if (globval.radiation) globval.dE = 0e0;

  cellp = &Cell[i0];
  dagetprm(Ascr, cellp->Alpha, cellp->Beta);
  memcpy(cellp->Nu, nu1, sizeof(Vector2));

  Ascr0 = Ascr;
  for (k = 0; k < n+2; k++)
    Ascr0[k] += tps(globval.CODvect[k]);

  Ascr1 = Ascr0;
  for (i = i0; i <= i1; i++) {
    Elem_Pass(i, Ascr1); cellp = &Cell[i];
    dagetprm(Ascr1, cellp->Alpha, cellp->Beta);
    for (k = 1; k <= 2; k++) {
      dnu[k-1] =
	(GetAngle(getmat(Ascr1, 2*k-1, 2*k-1), getmat(Ascr1, 2*k-1, 2*k)) -
	 GetAngle(getmat(Ascr0, 2*k-1, 2*k-1), getmat(Ascr0, 2*k-1, 2*k)))
	/(2e0*M_PI);

      if ((cellp->Elem.PL >= 0e0) && (dnu[k-1] < -1e-16))
	dnu[k-1] += 1e0;
      else if ((cellp->Elem.PL < 0e0) && (dnu[k-1] > 1e-16))
	dnu[k-1] -= 1e0;

      nu1[k-1] += dnu[k-1];

      cellp->Nu[k-1] = nu1[k-1];
      // cellp->Eta[k-1] = getmat(Ascr1, k, 5)*getmat(Ascr1, 6, 6) -
      //                getmat(Ascr1, k, 6)*getmat(Ascr1, 6, 5);
      // cellp->Etap[k-1] = getmat(Ascr1, 2*k, 5)*getmat(Ascr1, 6, 6) -
      //                 getmat(Ascr1, 2*k, 6)*getmat(Ascr1, 6, 5);
      cellp->Eta[k-1] = getmat(Ascr1, 2*k-1, 5);
      cellp->Etap[k-1] = getmat(Ascr1, 2*k, 5);
    }
    Ascr0 = Ascr1;
  }

  if (chroma && !globval.Cavity_on) Cell_Geteta(i0, i1, ring, dP);
}


void Ring_Getchrom(double dP)
{
  long int lastpos;
  int      k;
  Vector2  alpha = {0.0, 0.0}, beta = {0.0, 0.0}, gamma = {0.0, 0.0};
  Vector2     nu = {0.0, 0.0},  nu0 = {0.0, 0.0};
  
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
    printf("\nRing_Getchrom: unstable\n");
    return;
  }

  for (k = 0; k < 2; k++)
    globval.Chrom[k] = (nu[k]-nu0[k])/globval.dPcommon;
  
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

  if (trace) printf("enter Ring_GetTwiss\n");
  Ring_Twiss(chroma, dP);
  globval.Alphac = globval.OneTurnMat[ct_][delta_]/Cell[globval.Cell_nLoc].S;
  if (trace) printf("exit Ring_GetTwiss\n");
}
