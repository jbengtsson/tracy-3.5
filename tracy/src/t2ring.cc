/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

   t2ring.c: Functions for rings.                                             */


void LatticeType::ChamberOff(void)
{
  int i;

  for (i = 0; i <= conf.Cell_nLoc; i++) {
    elems[i]->maxampl[X_][0] = -max_ampl;
    elems[i]->maxampl[X_][1] = max_ampl;
    elems[i]->maxampl[Y_][0] = -max_ampl;
    elems[i]->maxampl[Y_][1] = max_ampl;
  }
  conf.chambre = false;
}


void LatticeType::PrintCh(void)
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

  for (i = 0; i <= conf.Cell_nLoc; i++)
    fprintf(f, "%4ld %15s  %6.2f  %7.3f  %7.3f  %7.3f\n",
	    i, elems[i]->PName, elems[i]->S,
	    elems[i]->maxampl[X_][0]*1E3, elems[i]->maxampl[X_][1]*1E3,
	    elems[i]->maxampl[Y_][1]*1E3);

  fclose(f);
}


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
    stable = false; nu[X_] = NAN; nu[Y_] = NAN;
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
      stable = false; nu[i] = NAN;
      printf("\nGetNu: unstable %s plane %10.3e\n", (i == 0)? "hor" : "ver", x);
      return;
    }
  }

  return;
}


bool Cell_GetABGN(Matrix &M, Vector2 &alpha, Vector2 &beta, Vector2 &gamma,
		  Vector2 &nu)
{
  bool   stable;
  int    k;
  double c = 0e0, s = 0e0;

  stable = true;
  for (k = 0; k < 2; k++) {
    c = (M[2*k][2*k]+M[2*k+1][2*k+1])/2e0;
    stable = (fabs(c) < 1e0);
    if (stable) {
      s = sqrt(1e0-sqr(c))*sgn(M[2*k][2*k+1]);
      alpha[k] = (M[2*k][2*k]-M[2*k+1][2*k+1])/(2e0*s);
      beta[k] = M[2*k][2*k+1]/s;
      gamma[k] = -M[2*k+1][2*k]/s;
    }
  }
  GetNu(nu, M);
  if (!stable) printf("Cell_GetABGN: unstable\n");
  return stable;
}


void LatticeType::Cell_Geteta(long i0, long i1, bool ring, double dP)
{
  long int i, lastpos;
  int      k;
  psVector xref;
  psVector codbuf[Cell_nLocMax+1];
  ElemType *elemp;

  const int n = 4;

  if (trace) printf("\nCell_Geteta: ring = %d\n", ring);

  if (ring)
    GetCOD(conf.CODimax, conf.CODeps, dP-conf.dPcommon/2e0, lastpos);
  else {
    CopyVec(n+2, conf.CODvect, xref); xref[4] = dP - conf.dPcommon/2e0;
    Cell_Pass(i0, i1, xref, lastpos);
  }

  for (i = i0; i <= i1; i++)
    CopyVec(n+2, elems[i]->BeamPos, codbuf[i]);

  if (ring)
    GetCOD(conf.CODimax, conf.CODeps, dP+conf.dPcommon/2e0, lastpos);
  else {
    CopyVec(n+2, conf.CODvect, xref); xref[4] = dP + conf.dPcommon/2e0;
    Cell_Pass(i0, i1, xref, lastpos);
  }

  for (i = i0; i <= i1; i++) {
    elemp = elems[i];
    for (k = 0; k < 2; k++) {
      elemp->Eta[k] = (elemp->BeamPos[2*k]-codbuf[i][2*k])/conf.dPcommon;
      elemp->Etap[k] =
	(elemp->BeamPos[2*k+1]-codbuf[i][2*k+1])/conf.dPcommon;
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


void LatticeType::Cell_Twiss(long i0, long i1, ss_vect<tps> &Ascr, bool chroma,
			     bool ring, double dP)
{
  long int     i;
  int          k;
  Vector2      nu1, dnu;
  ss_vect<tps> Ascr0, Ascr1;
  ElemType     *elemp;

  const int n = 4;

  for (k = 0; k < 2; k++) {
    nu1[k] = 0e0; dnu[k] = 0e0; 
  }

  if (conf.radiation) conf.dE = 0e0;

  elemp = elems[i0];
  dagetprm(Ascr, elemp->Alpha, elemp->Beta);
  memcpy(elemp->Nu, nu1, sizeof(Vector2));

  Ascr0 = Ascr;
  for (k = 0; k < n+2; k++)
    Ascr0[k] += tps(conf.CODvect[k]);

  Ascr1 = Ascr0;
  for (i = i0; i <= i1; i++) {
    elems[i]->Elem_Pass(conf, Ascr1); elemp = elems[i];
    dagetprm(Ascr1, elemp->Alpha, elemp->Beta);
    for (k = 1; k <= 2; k++) {
      dnu[k-1] =
	(GetAngle(getmat(Ascr1, 2*k-1, 2*k-1), getmat(Ascr1, 2*k-1, 2*k)) -
	 GetAngle(getmat(Ascr0, 2*k-1, 2*k-1), getmat(Ascr0, 2*k-1, 2*k)))
	/(2e0*M_PI);

      if ((elemp->PL >= 0e0) && (dnu[k-1] < -1e-16))
	dnu[k-1] += 1e0;
      else if ((elemp->PL < 0e0) && (dnu[k-1] > 1e-16))
	dnu[k-1] -= 1e0;

      nu1[k-1] += dnu[k-1];

      elemp->Nu[k-1] = nu1[k-1];
#if 0
      // Approximate:
      //   A = A0*A1 => [a_16, a_26] = [eta_x, eta_px]*A_long^-1
      // i.e., [a_15, a_25] != [0, 0].
      elemp->Eta[k-1] =
	getmat(Ascr1, k, 5)*getmat(Ascr1, 6, 6) -
	getmat(Ascr1, k, 6)*getmat(Ascr1, 6, 5);
      elemp->Etap[k-1] =
	getmat(Ascr1, 2*k, 5)*getmat(Ascr1, 6, 6) -
	getmat(Ascr1, 2*k, 6)*getmat(Ascr1, 6, 5);
#else
      elemp->Eta[k-1] = getmat(Ascr1, 2*k-1, 5);
      elemp->Etap[k-1] = getmat(Ascr1, 2*k, 5);
#endif
    }
    Ascr0 = Ascr1;
  }

  if (chroma && !conf.Cavity_on) Cell_Geteta(i0, i1, ring, dP);
}


void LatticeType::Ring_Getchrom(double dP)
{
  long int lastpos;
  int      k;
  Vector2  alpha = {0.0, 0.0}, beta = {0.0, 0.0}, gamma = {0.0, 0.0};
  Vector2     nu = {0.0, 0.0},  nu0 = {0.0, 0.0};
  
  if (dP != 0e0)
    printf("\nRing_Getchrom: linear chromaticity for delta = %e\n", dP);
  
  GetCOD(conf.CODimax, conf.CODeps, dP-conf.dPcommon*0.5, lastpos);
  if (!conf.codflag) {
    printf("\nRing_Getchrom: closed orbit finder failed\n");
    return;
  }
  conf.stable = Cell_GetABGN(conf.OneTurnMat, alpha, beta, gamma, nu0);
  if (!stable) {
    printf("\nRing_Getchrom: unstable\n");
    return;
  }
  
  GetCOD(conf.CODimax, conf.CODeps, dP+conf.dPcommon*0.5, lastpos);
  if (!conf.codflag) {
    printf("\nRing_Getchrom: closed orbit finder failed");
    return;
  }
  Cell_GetABGN(conf.OneTurnMat, alpha, beta, gamma, nu);
  if (!conf.stable) {
    printf("\nRing_Getchrom: unstable\n");
    return;
  }

  for (k = 0; k < 2; k++)
    conf.Chrom[k] = (nu[k]-nu0[k])/conf.dPcommon;
  
  conf.chromflag = true;
}


void LatticeType::Ring_Twiss(bool chroma, double dP)
{
  long int     lastpos = 0;
  int          n = 0;
  Vector2      alpha={0.0, 0.0}, beta={0.0, 0.0};
  Vector2      gamma={0.0, 0.0}, nu={0.0, 0.0};
  Matrix       R;
  ss_vect<tps> AScr;

  n = (conf.Cavity_on)? 6 : 4;

  GetCOD(conf.CODimax, conf.CODeps, dP, lastpos);

  if (!conf.codflag) return;

  // Check if stable
  stable = Cell_GetABGN(conf.OneTurnMat, alpha, beta, gamma, nu);
  if (!stable) {
    printf("Ring_Twiss: unstable\n");
    return;
  }
  // Get eigenvalues and eigenvectors for the one turn transfer matrix
  GDiag(n, elems[conf.Cell_nLoc]->S, conf.Ascr, conf.Ascrinv, R,
	conf.OneTurnMat, conf.Omega, conf.Alphac);

  // AScr = putlinmat(n, conf.Ascr);
  AScr = putlinmat(6, conf.Ascr);
  if (!conf.Cavity_on) {
    // AScr[delta_] = 0.0; AScr[ct_] = 0.0;
    AScr[delta_] = tps(0e0, delta_+1); AScr[ct_] = 0e0;
  }

  Cell_Twiss(0, conf.Cell_nLoc, AScr, chroma, true, dP);

  memcpy(conf.TotalTune, elems[conf.Cell_nLoc]->Nu, sizeof(Vector2));
  conf.tuneflag = true;

  if (chroma && !conf.Cavity_on) {
    Ring_Getchrom(dP); GetCOD(conf.CODimax, conf.CODeps, dP, lastpos);
  }
}


void LatticeType::Ring_GetTwiss(bool chroma, double dP)
{

  if (trace) printf("enter Ring_GetTwiss\n");
  Ring_Twiss(chroma, dP);
  conf.Alphac = conf.OneTurnMat[ct_][delta_]/elems[conf.Cell_nLoc]->S;
  if (trace) printf("exit Ring_GetTwiss\n");
}


void LatticeType::TraceABN(long i0, long i1, const Vector2 &alpha,
			   const Vector2 &beta, const Vector2 &eta,
			   const Vector2 &etap, const double dP)
{
  long          i, j;
  double        sb;
  ss_vect<tps>  Ascr;

  UnitMat(6, conf.Ascr);
  for (i = 1; i <= 2; i++) {
    sb = sqrt(beta[i-1]); j = i*2 - 1;
    conf.Ascr[j-1][j-1] = sb;               conf.Ascr[j-1][j] = 0.0;
    conf.Ascr[j][j - 1] = -(alpha[i-1]/sb); conf.Ascr[j][j] = 1/sb;
  }
  conf.Ascr[0][4] = eta[0]; conf.Ascr[1][4] = etap[0];
  conf.Ascr[2][4] = eta[1]; conf.Ascr[3][4] = etap[1];

  for (i = 0; i < 6; i++)
    conf.CODvect[i] = 0.0;
  conf.CODvect[4] = dP;

  for (i = 0; i <= 5; i++) {
    Ascr[i] = tps(conf.CODvect[i]);
    for (j = 0; j <= 5; j++)
      Ascr[i] += conf.Ascr[i][j]*tps(0.0, j+1);
    Cell_Twiss(i0, i1, Ascr, false, false, dP);
  }

}


void LatticeType::ttwiss(const Vector2 &alpha, const Vector2 &beta,
			 const Vector2 &eta, const Vector2 &etap,
			 const double dP)
{
  TraceABN(0, conf.Cell_nLoc, alpha, beta, eta, etap, dP);
}


double get_curly_H(const double alpha_x, const double beta_x,
		   const double eta_x, const double etap_x)
{
  double curly_H, gamma_x;

  gamma_x = (1.0+sqr(alpha_x))/beta_x;

  curly_H = gamma_x*sqr(eta_x) + 2.0*alpha_x*eta_x*etap_x + beta_x*sqr(etap_x);

  return curly_H;
}


void get_dI_eta_5(const int k, std::vector<ElemType*> Elem)
{
  double       L, K, h, b2, alpha, beta, gamma, psi, eta, etap;
  ss_vect<tps> Id;
  MpoleType    *Mp;

  Id.identity();

  Mp = dynamic_cast<MpoleType*>(Elem[k]);

  L = Elem[k]->PL;
  h = Mp->Pirho;
  b2 = Mp->PBpar[Quad+HOMmax];
  K = b2 + sqr(Mp->Pirho);
  psi = sqrt(fabs(K))*L;
  alpha = Elem[k-1]->Alpha[X_]; beta = Elem[k-1]->Beta[X_];
  gamma = (1e0+sqr(alpha))/beta;
  eta = Elem[k-1]->Eta[X_]; etap = Elem[k-1]->Etap[X_];

  Elem[k]->dI[1] += L*eta*h;
  Elem[k]->dI[2] += L*sqr(h);
  Elem[k]->dI[3] += L*fabs(cube(h));

  if (K > 0e0) {
    Elem[k]->dI[4] +=
      h/K*(2e0*b2+sqr(h))
      *((eta*sqrt(K)*sin(psi)+etap*(1e0-cos(psi)))
	+ h/sqrt(K)*(psi-sin(psi)));

    Elem[k]->dI[5] +=
      L*fabs(cube(h))
      *(gamma*sqr(eta)+2e0*alpha*eta*etap+beta*sqr(etap))
      - 2e0*pow(h, 4)/(pow(K, 3e0/2e0))
      *(sqrt(K)*(alpha*eta+beta*etap)*(cos(psi)-1e0)
	+(gamma*eta+alpha*etap)*(psi-sin(psi)))
      + fabs(pow(h, 5))/(4e0*pow(K, 5e0/2e0))
      *(2e0*alpha*sqrt(K)*(4e0*cos(psi)-cos(2e0*psi)-3e0)
	+beta*K*(2e0*psi-sin(2e0*psi))
	+gamma*(6e0*psi-8e0*sin(psi)+sin(2e0*psi)));
  } else {
    K = fabs(K);

    Elem[k]->dI[4] +=
      h/K*(2e0*b2+sqr(h))
      *((eta*sqrt(K)*sinh(psi)-etap*(1e0-cosh(psi)))
	- h/sqrt(K)*(psi-sinh(psi)));

    Elem[k]->dI[5] +=
      L*fabs(cube(h))*
      (gamma*sqr(eta)+2e0*alpha*eta*etap+beta*sqr(etap))
      + 2e0*pow(h, 4)/(pow(K, 3e0/2e0))
      *(sqrt(K)*(alpha*eta+beta*etap)*(cosh(psi)-1e0)
	+(gamma*eta+alpha*etap)*(psi-sinh(psi)))
      + fabs(pow(h, 5))/(4e0*pow(K, 5e0/2e0))
      *(2e0*alpha*sqrt(K)*(4e0*cosh(psi)-cosh(2e0*psi)-3e0)
	-beta*K*(2e0*psi-sinh(2e0*psi))
	+gamma*(6e0*psi-8e0*sinh(psi)+sinh(2e0*psi)));
  }
}


void LatticeType::get_I(double I[], const bool prt)
{
  int j, k;

  for (k = 0; k <= 5; k++)
    I[k] = 0e0;

  if (prt) {
    printf("\nget_I:\n");
    printf("\n      name               s     curly_H      I_1        I_2"
	   "        I_3        I_4        I_5      alpha_x    beta_x"
	   "     eta_x      eta'_x     alpha_y    beta_y\n\n");
  }
  for (j = 0; j <= conf.Cell_nLoc; j++)
    if ((elems[j]->Pkind == drift) || (elems[j]->Pkind == Mpole) ||
	(elems[j]->Pkind == Wigl) ||
	(elems[j]->Pkind == marker)) {
      if (prt)
	printf("%5d %-10s %6.3f %10.3e",
	       j, elems[j]->PName, elems[j]->S,
	       elems[j]->curly_dH_x);
      for (k = 1; k <= 5; k++) {
	I[k] += elems[j]->dI[k];
	if (prt) printf(" %10.3e", elems[j]->dI[k]);
      }
      if (prt)
	printf(" %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	       elems[j]->Alpha[X_], elems[j]->Beta[X_],
	       elems[j]->Eta[X_], elems[j]->Etap[X_],
	       elems[j]->Alpha[Y_], elems[j]->Beta[Y_]);
    }
}


template<typename T>
void LatticeType::Elem_Pass_Lin(ss_vect<T> ps)
{
  long int  k;
  MpoleType *Mp;

  for (k = 0; k <= conf.Cell_nLoc; k++) {
    if (elems[k]->Pkind == Mpole) { 
      Mp = dynamic_cast<MpoleType*>(elems[k]);
      if ((Mp->Pthick == thick) && (Mp->Porder <= Quad)) {
	ps = is_double<ss_vect<T> >::ps(Mp->M_lin*ps);
	
	if (conf.emittance && !conf.Cavity_on
	    && (elems[k]->PL != 0e0) && (Mp->Pirho != 0e0))
	  get_dI_eta_5(k, elems);
      }
    } else
      elems[k]->Elem_Pass(conf, ps);
  }
}


void LatticeType::get_eps_x(double &eps_x, double &sigma_delta, double &U_0,
			    double J[], double tau[], double I[],
			    const bool prt)
{
  bool         cav, emit;
  int          k;
  ss_vect<tps> A;

  const double
    C_q_scl = 1e18*C_q/sqr(m_e),
    E_0     = 1e9*conf.Energy,
    C       = elems[conf.Cell_nLoc]->S,
    T_0     = C/c0;

  /* Note:

        T
       M  J M = J,

        -1       T           |  0  I |        T   | beta   -alpha |
       A   = -J A  J,    J = |       |,    A A  = |               |
                             | -I  0 |            | -alpha  gamma |

     Transform to Floquet Space:

        -1           T
       A   eta = -J A  J eta,

               -1      T  -1                T    T
       H~ = ( A   eta )  A   eta = ( J eta )  A A  ( J eta )

  */

  cav = conf.Cavity_on; emit = conf.emittance;

  conf.Cavity_on = false; conf.emittance = false;
  Ring_GetTwiss(false, 0.0);
  A = putlinmat(6, conf.Ascr); A += conf.CODvect;
  conf.emittance = true;
  Elem_Pass_Lin(A);
  get_I(I, false);

  U_0 = 1e9*C_gamma*pow(conf.Energy, 4)*I[2]/(2e0*M_PI);
  eps_x = C_q_scl*sqr(conf.Energy)*I[5]/(I[2]-I[4]);
  sigma_delta = sqrt(C_q_scl*sqr(conf.Energy)*I[3]/(2e0*I[2]+I[4]));
  J[X_] = 1e0 - I[4]/I[2]; J[Z_] = 2e0 + I[4]/I[2]; J[Y_] = 4e0 - J[X_] - J[Z_];

  for (k = 0; k < 3; k++)
    tau[k] = 4e0*M_PI*T_0/(C_gamma*cube(1e-9*E_0)*J[k]*I[2]);

  if (prt) {
    printf("\n  I[1..5]:");
    for (k = 1; k <= 5; k++)
      printf(" %10.3e", I[k]);
    printf("\n");

    printf("\n  U_0   [keV]    = %5.1f\n", 1e-3*U_0);
    printf("  eps_x [nm.rad] = %6.4f\n", 1e9*eps_x);
    printf("  sigma_delta    = %9.3e\n", sigma_delta);
    printf("  J              = [%5.3f, %5.3f, %5.3f]\n", J[X_], J[Y_], J[Z_]);
    printf("  tau   [msec]   = [%e, %e, %e]\n",
	   1e3*tau[X_], 1e3*tau[Y_], 1e3*tau[Z_]);
  }

  conf.Cavity_on = cav; conf.emittance = emit;
}


double get_code(const ConfigType &conf, ElemType &Cell)
{
  double    code;
  MpoleType *M;

  switch (Cell.Pkind) {
  case drift:
    code = 0.0;
    break;
  case Mpole:
    M = dynamic_cast<MpoleType*>(&Cell);
    if (M->Pirho != 0.0)
      code = sgn(M->Pirho)*0.5;
    else if (M->PBpar[Quad+HOMmax] != 0)
      code = sgn(M->PBpar[Quad+HOMmax]);
    else if (M->PBpar[Sext+HOMmax] != 0)
      code = 1.5*sgn(M->PBpar[Sext+HOMmax]);
    else if (M->PBpar[Oct+HOMmax] != 0)
      code = 1.75*sgn(M->PBpar[Oct+HOMmax]);
    else if (M->PBpar[Dec+HOMmax] != 0)
      code = 1.75*sgn(M->PBpar[Dec+HOMmax]);
    else if (M->PBpar[Dodec+HOMmax] != 0)
      code = 1.75*sgn(M->PBpar[Dodec+HOMmax]);
    else if (Cell.Fnum == conf.bpm)
      code = 2.0;
    else
      code = 0.0;
    break;
  default:
    code = 0.0;
    break;
  }

  return code;
}


void LatticeType::prt_lat(const int loc1, const int loc2, const char *fname,
			  const bool all)
{
  long int i = 0;
  double   I5 = 0e0;
  FILE     *outf;

  outf = file_write(fname);
  fprintf(outf,
	  "#        name             s     code"
	  "   alphax   betax     nux      etax    etapx");
  fprintf(outf,
	  "     alphay   betay     nuy      etay    etapy      I5\n");
  fprintf(outf,
	  "#                        [m]"
	  "                     [m]                [m]");
  fprintf(outf, "                        [m]                [m]\n");
  fprintf(outf, "#\n");

  for (i = loc1; i <= loc2; i++) {
    if (all || (elems[i]->Fnum == conf.bpm)) {
      fprintf(outf,
	      "%4ld %15s %9.5f %4.1f %9.5f %8.5f %8.5f %8.5f %8.5f"
	      " %9.5f %8.5f %8.5f %8.5f %8.5f  %8.2e\n",
	      i, elems[i]->PName, elems[i]->S, get_code(conf, *elems[i]),
	      elems[i]->Alpha[X_], elems[i]->Beta[X_], elems[i]->Nu[X_],
	      elems[i]->Eta[X_], elems[i]->Etap[X_], elems[i]->Alpha[Y_],
	      elems[i]->Beta[Y_], elems[i]->Nu[Y_], elems[i]->Eta[Y_],
	      elems[i]->Etap[Y_], I5);
    }
  }

//  fprintf(outf, "\n");
//  fprintf(outf, "# emittance: %5.3f nm.rad\n", get_eps_x());

  fclose(outf);
}


void LatticeType::prt_lat(const char *fname, const bool all)
{
  prt_lat(0, conf.Cell_nLoc, fname, all);
}


void LatticeType::Cell_Twiss(const long int i0, const long int i1)
{
  long int     i;
  int          k, nu_int[2];
  double       alpha[2], beta[2], dnu[2], eta[2], etap[2];
  ss_vect<tps> A;

  for (k = 0; k < 2; k++)
    nu_int[k] = 0;

  for (i = i0; i <= i1; i++) {
    A = putlinmat(6, elems[i]->A);
    get_ab(A, alpha, beta, dnu, eta, etap);

    for (k = 0; k < 2; k++) {
      elems[i]->Alpha[k] = alpha[k]; elems[i]->Beta[k] = beta[k];
      elems[i]->Nu[k] = nu_int[k] + dnu[k];

      if (i > i0) {
	if((elems[i]->Nu[k] < elems[i-1]->Nu[k])
	   && (elems[i]->PL >= 0e0)) {
	  elems[i]->Nu[k] += 1e0; nu_int[k] += 1;
	} else if((elems[i]->Nu[k] > elems[i-1]->Nu[k]) &&
		  (elems[i]->PL < 0e0))
	  nu_int[k] -= 1;
      }

      elems[i]->Eta[k] = eta[k]; elems[i]->Etap[k] = etap[k];
    }
  }
}


void LatticeType::prt_lat(const int loc1, const int loc2, const char *fname,
			  const bool all, const int n)
{
  long int        i = 0;
  int             j, k;
  double          s, h;
  double          alpha[2], beta[2], nu[2], dnu[2], eta[2], etap[2], dnu1[2];
  double          curly_H;
  MpoleType       *Mp;
  ss_vect<double> eta_Fl;
  ss_vect<tps>    A, A_CS;
  FILE            *outf;

  const double  c1 = 1e0/(2e0*(2e0-pow(2e0, 1e0/3e0))), c2 = 0.5e0-c1;
  const double  d1 = 2e0*c1, d2 = 1e0-2e0*d1;

  outf = file_write(fname);
  fprintf(outf, "#        name           s   code"
	        "    alphax   betax     nux       etax       etapx");
  fprintf(outf, "       alphay   betay     nuy      etay    etapy\n");
  fprintf(outf, "#                      [m]"
	        "                    [m]                 [m]");
  fprintf(outf, "                             [m]                [m]\n");
  fprintf(outf, "#\n");

  for (i = loc1; i <= loc2; i++) {
    if (all || (elems[i]->Fnum == conf.bpm)) {
      if ((i != 0) &&
	  ((elems[i]->Pkind == drift) ||
	   ((elems[i]->Pkind == Mpole) && (elems[i]->PL != 0e0)))) {
	Mp = dynamic_cast<MpoleType*>(elems[i]);

	for (k = 0; k < 2; k++) {
	  alpha[k] = elems[i-1]->Alpha[k];
	  beta[k] = elems[i-1]->Beta[k];
	  nu[k] = elems[i-1]->Nu[k];
	  eta[k] = elems[i-1]->Eta[k]; etap[k] = elems[i-1]->Etap[k];
	}

	A = get_A(alpha, beta, eta, etap);

	s = elems[i]->S - elems[i]->PL; h = elems[i]->PL/n;

	for (j = 1; j <= n; j++) {
	  s += h;

	  if (elems[i]->Pkind == drift)
	    Drift(conf, h, A);
	  else if (elems[i]->Pkind == Mpole) {
	    if ((j == 1) && (Mp->Pirho != 0e0))
	      EdgeFocus(conf, Mp->Pirho, Mp->PTx1, Mp->Pgap, A);

	    Drift(conf, c1*h, A);
	    thin_kick(conf, Quad, Mp->PBpar, d1*h, Mp->Pirho, Mp->Pirho, A);
	    Drift(conf, c2*h, A);
	    thin_kick(conf, Quad, Mp->PBpar, d2*h, Mp->Pirho, Mp->Pirho, A);
	    Drift(conf, c2*h, A);
	    thin_kick(conf, Quad, Mp->PBpar, d1*h, Mp->Pirho, Mp->Pirho, A);
	    Drift(conf, c1*h, A);

	    if ((j == n) && (Mp->Pirho != 0e0))
	      EdgeFocus(conf, Mp->Pirho, Mp->PTx2, Mp->Pgap, A);
	  }

	  get_ab(A, alpha, beta, dnu, eta, etap);

	  if(elems[i]->PL < 0e0)
	    for (k = 0; k < 2; k++)
	      dnu[k] -= 1e0;

	  A_CS = get_A_CS(2, A, dnu1);

	  eta_Fl.zero();
	  for (k = 0; k < 2; k++) {
	    eta_Fl[2*k] = eta[k]; eta_Fl[2*k+1] = etap[k];
	  }
	  eta_Fl = (Inv(A_CS)*eta_Fl).cst();
	  curly_H = sqr(eta_Fl[x_]) + sqr(eta_Fl[px_]);

	  fprintf(outf, "%4ld %15s %6.2f %4.1f"
		  " %9.5f %8.5f %8.5f %11.8f %11.8f"
		  " %9.5f %8.5f %8.5f %8.5f %8.5f %10.3e %10.3e %10.3e\n",
		  i, elems[i]->PName, s, get_code(conf, *elems[i]),
		  alpha[X_], beta[X_], nu[X_]+dnu[X_], eta[X_], etap[X_],
		  alpha[Y_], beta[Y_], nu[Y_]+dnu[Y_], eta[Y_], etap[Y_],
		  eta_Fl[x_], eta_Fl[px_], curly_H);
	}
      } else {
	A = get_A(elems[i]->Alpha, elems[i]->Beta, elems[i]->Eta,
		  elems[i]->Etap);

	eta_Fl.zero();
	for (k = 0; k < 2; k++) {
	  eta_Fl[2*k] = elems[i]->Eta[k];
	  eta_Fl[2*k+1] = elems[i]->Etap[k];
	}
	eta_Fl = (Inv(A)*eta_Fl).cst();
	curly_H = sqr(eta_Fl[x_]) + sqr(eta_Fl[px_]);

	fprintf(outf, "%4ld %15s %6.2f %4.1f"
		" %9.5f %8.5f %8.5f %11.8f %11.8f"
		" %9.5f %8.5f %8.5f %8.5f %8.5f %10.3e %10.3e %10.3e\n",
		i, elems[i]->PName, elems[i]->S, get_code(conf, *elems[i]),
		elems[i]->Alpha[X_], elems[i]->Beta[X_], elems[i]->Nu[X_],
		elems[i]->Eta[X_], elems[i]->Etap[X_],
		elems[i]->Alpha[Y_], elems[i]->Beta[Y_], elems[i]->Nu[Y_],
		elems[i]->Eta[Y_], elems[i]->Etap[Y_],
		eta_Fl[x_], eta_Fl[px_], curly_H);
      }
    }
  }

  fclose(outf);
}


void LatticeType::prt_lat(const char *fname, const bool all, const int n)
{
  prt_lat(0, conf.Cell_nLoc, fname, all, n);
}


void LatticeType::prt_chrom_lat(void)
{
  long int  i;
  double    dbeta_ddelta[Cell_nLocMax][2], detax_ddelta[Cell_nLocMax];
  double    ksi[Cell_nLocMax][2];
  MpoleType *M;
  FILE      *outf;

  printf("\nprt_chrom_lat:\n  calling Ring_GetTwiss with delta != 0\n");
  Ring_GetTwiss(true, conf.dPcommon);
  for (i = 0; i <= conf.Cell_nLoc; i++) {
    dbeta_ddelta[i][X_] = elems[i]->Beta[X_];
    dbeta_ddelta[i][Y_] = elems[i]->Beta[Y_];
    detax_ddelta[i] = elems[i]->Eta[X_];
  }
  printf("  calling Ring_GetTwiss with delta != 0\n");
  Ring_GetTwiss(true, -conf.dPcommon);
  ksi[0][X_] = 0.0; ksi[0][Y_] = 0.0;
  for (i = 0; i <= conf.Cell_nLoc; i++) {
    dbeta_ddelta[i][X_] -= elems[i]->Beta[X_];
    dbeta_ddelta[i][Y_] -= elems[i]->Beta[Y_];
    detax_ddelta[i] -= elems[i]->Eta[X_];
    dbeta_ddelta[i][X_] /= 2.0*conf.dPcommon;
    dbeta_ddelta[i][Y_] /= 2.0*conf.dPcommon;
    detax_ddelta[i] /= 2.0*conf.dPcommon;
    if (i != 0) {
      ksi[i][X_] = ksi[i-1][X_]; ksi[i][Y_] = ksi[i-1][Y_];
    }
    if (elems[i]->Pkind == Mpole) {
      M = dynamic_cast<MpoleType*>(elems[i]);
      ksi[i][X_] -=
	M->PBpar[Quad+HOMmax]*elems[i]->PL*elems[i]->Beta[X_]
	/(4.0*M_PI);
      ksi[i][Y_] +=
	M->PBpar[Quad+HOMmax]*elems[i]->PL*elems[i]->Beta[Y_]
	/(4.0*M_PI);
    }
  }

  outf = file_write("chromout");
  fprintf(outf, "#     name              s    code"
	        "  bx*ex  sqrt(bx*by)  dbx/dd*ex  bx*dex/dd"
	        "  by*ex  dby/dd*ex by*dex/dd  ksix  ksiy"
	        "  dbx/dd  bx/dd dex/dd\n");
  fprintf(outf, "#                      [m]          [m]"
	        "      [m]          [m]       [m]");
  fprintf(outf, "       [m]      [m]       [m]\n");
  fprintf(outf, "#\n");
  for (i = 0; i <= conf.Cell_nLoc; i++) {
    fprintf(outf,
	    "%4ld %15s %6.2f %4.1f  %6.3f  %8.3f    %8.3f   %8.3f"
	    "   %6.3f %8.3f   %8.3f  %5.2f  %5.2f  %6.3f  %6.3f  %6.3f\n",
	    i, elems[i]->PName, elems[i]->S, get_code(conf, *elems[i]),
	    elems[i]->Beta[X_]*elems[i]->Eta[X_],
	    sqrt(elems[i]->Beta[X_]*elems[i]->Beta[Y_]),
	    dbeta_ddelta[i][X_]*elems[i]->Eta[X_],
	    detax_ddelta[i]*elems[i]->Beta[X_],
	    elems[i]->Beta[Y_]*elems[i]->Eta[X_],
	    dbeta_ddelta[i][Y_]*elems[i]->Eta[X_],
	    detax_ddelta[i]*elems[i]->Beta[Y_],
	    ksi[i][X_], ksi[i][Y_],
	    dbeta_ddelta[i][X_], dbeta_ddelta[i][Y_], detax_ddelta[i]);
  }
  fclose(outf);
}

#if 0

/* Local variables for Ring_Fittune: */

struct LOC_Ring_Fittune
{
  jmp_buf _JL999;
};


void LatticeType::checkifstable_(struct LOC_Ring_FitDisp *LINK)
{
  if (!stable) {
    printf("  lattice is unstable\n");
    longjmp(LINK->_JL999, 1);
  }
}


void LatticeType::shiftk(long Elnum, double dk, struct LOC_Ring_Fittune *LINK)
{
  ElemType  *elemp;
  MpoleType *M;

  elemp = elems[Elnum];
  M = dynamic_cast<MpoleType*>(elemp);
  M->PBpar[Quad+HOMmax] += dk;
  M->Mpole_SetPB(lat, elemp->Fnum, elemp->Knum, (long)Quad);
}


void LatticeType::checkifstable(struct LOC_Ring_Fittune *LINK)
{
  if (!stable) {
    printf("  lattice is unstable\n");
    longjmp(LINK->_JL999, 1);
  }
}


void LatticeType::Ring_Fittune(Vector2 &nu, double eps, iVector2 &nq, long qf[],
			       long qd[], double dkL, long imax)
{
  struct LOC_Ring_Fittune V;

  int      i, j, k;
  Vector2  nu0, nu1;
  psVector  dkL1, dnu;
  Matrix A;

  const double dP = 0e0;

  if (setjmp(V._JL999)) return;

  if (trace)
    printf("  Tune fit, nux =%10.5f, nuy =%10.5f, eps =% .3E,"
	   " imax =%4ld, dkL = % .5E\n", nu[0], nu[1], eps, imax, dkL);
  Ring_GetTwiss(false, dP); checkifstable(&V);
  memcpy(nu0, conf.TotalTune, sizeof(Vector2));
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
      Ring_GetTwiss(false, dP);
      nu1[0] = conf.TotalTune[0]; nu1[1] = conf.TotalTune[1];
//      GetCOD(conf.CODimax, conf.CODeps, dP, lastpos);
//      Cell_GetABGN(conf.OneTurnMat, alpha, beta, gamma, nu1);
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
    Ring_GetTwiss(false, dP); checkifstable(&V);
    memcpy(nu0, conf.TotalTune, sizeof(Vector2));
    if (trace)
      printf("  Nux = %10.6f%10.6f, Nuy = %10.6f%10.6f,"
	     " QF*L = % .5E, QD*L = % .5E @%3d\n",
	     nu0[0], nu1[0], nu0[1], nu1[1],
	     Elem_GetKval(elems[qf[0]]->Fnum, 1, (long)Quad),
	     Elem_GetKval(elems[qd[0]]->Fnum, 1, (long)Quad), i);
  } while (sqrt(sqr(nu[0]-nu0[0])+sqr(nu[1]-nu0[1])) >= eps && i != imax);
}


void LatticeType::shiftkp(long Elnum, double dkp)
{
  ElemType  *elemp;
  MpoleType *M;

  elemp = elems[Elnum];
  M = dynamic_cast<MpoleType*>(elemp);
  M->PBpar[Sext+HOMmax] += dkp;
  M->Mpole_SetPB(lat, elemp->Fnum, elemp->Knum, (long)Sext);
}


void LatticeType::Ring_Fitchrom(Vector2 &ksi, double eps, iVector2 &ns,
				long sf[], long sd[], double dkpL, long imax)
{
  bool      rad;
  long int  lastpos;
  int       i, j, k;
  Vector2   ksi0;
  psVector    dkpL1, dksi;
  Matrix    A;

  const double dP = 0e0;

  if (trace)
    printf("  Chromaticity fit, ksix =%10.5f, ksiy =%10.5f, eps =% .3E"
	   ", imax =%4ld, dkpL =%10.5f\n", ksi[0], ksi[1], eps, imax, dkpL);

  /* Turn off radiation */
  rad = conf.radiation; conf.radiation = false;
  GetCOD(conf.CODimax, conf.CODeps, dP, lastpos); Ring_Getchrom(dP);
  for (j = 0; j <= 1; j++)
    ksi0[j] = conf.Chrom[j];
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
      GetCOD(conf.CODimax, conf.CODeps, dP, lastpos); Ring_Getchrom(dP);
      for (k = 0; k <= 1; k++) {
	dksi[k] = conf.Chrom[k] - ksi0[k];
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
    GetCOD(conf.CODimax, conf.CODeps, dP, lastpos); Ring_Getchrom(dP);
    for (j = 0; j <= 1; j++)
      ksi0[j] = conf.Chrom[j];
    if (trace)
      printf("  ksix =%10.6f, ksiy =%10.6f, SF = % .5E, SD = % .5E @%3d\n",
	     ksi0[0], ksi0[1], Elem_GetKval(elems[sf[0]]->Fnum, 1, (long)Sext),
	     Elem_GetKval(elems[sd[0]]->Fnum, 1, (long)Sext), i);
  } while (sqrt(sqr(ksi[0]-ksi0[0])+sqr(ksi[1]-ksi0[1])) >= eps && i != imax);
_L999:
  /* Restore radiation */
  conf.radiation = rad;
}


/* Local variables for Ring_FitDisp: */
struct LOC_Ring_FitDisp
{
  jmp_buf _JL999;
};


void LatticeType::shiftk_(long Elnum, double dk, struct LOC_Ring_FitDisp *LINK)
{
  ElemType  *elemp;
  MpoleType *M;

  elemp = elems[Elnum];
  M = dynamic_cast<MpoleType*>(elemp);
  M->PBpar[Quad+HOMmax] += dk;
  M->Mpole_SetPB(lat, elemp->Fnum, elemp->Knum, (long)Quad);
}


void LatticeType::Ring_FitDisp(long pos, double eta, double eps, long nq,
			       long q[], double dkL, long imax)
{
  /*pos : integer; eta, eps : double;
                         nq : integer; var q : fitvect;
                         dkL : double; imax : integer*/
  struct LOC_Ring_FitDisp V;

  int     i, j;
  double  dkL1, Eta0, deta;
  bool    rad = false;

  const double dP = 0e0;

  if (setjmp(V._JL999)) goto _L999;

  if (trace)
    printf("  Dispersion fit, etax =%10.5f, eps =% .3E"
	   ", imax =%4ld, dkL =%10.5f\n",
	   eta, eps, imax, dkL);
  /* Turn off radiation */
  rad = conf.radiation; conf.radiation = false;
  Ring_GetTwiss(true, dP); checkifstable_(&V);
  Eta0 = elems[pos]->Eta[0];
  i = 0;
  while (fabs(eta-Eta0) > eps && i < imax) {
    i++;
    for (j = 0; j < nq; j++)
      shiftk_(q[j], dkL, &V);
    Ring_GetTwiss(true, dP); checkifstable_(&V);
    deta = elems[pos]->Eta[0] - Eta0;
    if (deta == 0.0) {
      printf("  deta is 0\n");
      goto _L999;
    }
    dkL1 = (eta-Eta0)*dkL/deta - dkL;
    for (j = 0; j < nq; j++)
      shiftk_(q[j], dkL1, &V);
    Ring_GetTwiss(true, dP); checkifstable_(&V);
    Eta0 = elems[pos]->Eta[0];
    if (trace)
      printf("  Dispersion = % .5E, kL =% .5E @%3d\n",
	     Eta0, Elem_GetKval(elems[q[0]]->Fnum, 1, (long)Quad), i);
  }
_L999:
  /* Restore radiation */
  conf.radiation = rad;
}

#endif
