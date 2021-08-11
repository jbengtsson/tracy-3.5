
void LatticeType::ChamberOff(void)
{
  int i;

  for (i = 0; i <= conf.Cell_nLoc; i++) {
    elems[i]->maxampl[X_][0] = -max_ampl;
    elems[i]->maxampl[X_][1] =  max_ampl;
    elems[i]->maxampl[Y_][0] = -max_ampl;
    elems[i]->maxampl[Y_][1] =  max_ampl;
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
	    i, elems[i]->Name.c_str(), elems[i]->S,
	    elems[i]->maxampl[X_][0]*1E3, elems[i]->maxampl[X_][1]*1E3,
	    elems[i]->maxampl[Y_][1]*1E3);

  fclose(f);
}


void GetNu(std::vector<double> &nu, std::vector< std::vector<double> > &M)
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

  int       i;
  double    sgn, detp, detm, b, c, tr[2], b2mc, x;
  arma::mat M1(ps_dim, ps_dim);

  const int n = 4;

  M1 = stlmattomat(M);
  for (i = 0; i < n; i++)
    M1(i, i) -= 1e0;   
  detp = det(M1);
  for (i = 0; i < n; i++)
    M1(i, i) += 2e0;
  detm = det(M1);
  for (i = 0; i < n; i++)
    M1(i, i) -= 1e0;

  for (i = 0; i < 2; i++)
    tr[i] = (M1(2*i, 2*i)+M1(2*i+1, 2*i+1));
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
      if (M1(2*i, 2*i+1) < 0e0) nu[i] = 1e0 - nu[i];
    } else {
      stable = false; nu[i] = NAN;
      printf("\nGetNu: unstable %s plane %10.3e\n", (i == 0)? "hor" : "ver", x);
      return;
    }
  }

  return;
}


bool Cell_GetABGN(std::vector< std::vector<double> > &M,
		  std::vector<double> &alpha,
		  std::vector<double> &beta, std::vector<double> &gamma,
		  std::vector<double> &nu)
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
  long int        i, lastpos;
  int             k;
  ss_vect<double> xref;
  ss_vect<double> codbuf[Cell_nLocMax+1];
  ElemType        *elemp;

  if (conf.trace) printf("\nCell_Geteta: ring = %d\n", ring);

  if (ring)
    GetCOD(conf.CODimax, conf.CODeps, dP-conf.dPcommon/2e0, lastpos);
  else {
    xref = vectops(conf.CODvect); xref[4] = dP - conf.dPcommon/2e0;
    Cell_Pass(i0, i1, xref, lastpos);
  }

  for (i = i0; i <= i1; i++)
    codbuf[i] = elems[i]->BeamPos;

  if (ring)
    GetCOD(conf.CODimax, conf.CODeps, dP+conf.dPcommon/2e0, lastpos);
  else {
    xref = vectops(conf.CODvect); xref[4] = dP + conf.dPcommon/2e0;
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


void getprm(arma::mat &Ascr, std::vector<double> &alpha,
	    std::vector<double> &beta)
{
  int k;

  for (k = 0; k < 2; k++) {
    alpha[k] =
      -(Ascr(2*k, 2*k)*Ascr(2*k+1, 2*k) + Ascr(2*k, 2*k+1)*Ascr(2*k+1, 2*k+1));
    beta[k] = sqr(Ascr(2*k, 2*k)) + sqr(Ascr(2*k, 2*k+1));
  }
}


void dagetprm(ss_vect<tps> &Ascr, std::vector<double> &alpha,
	      std::vector<double> &beta)
{
  int k;

  for (k = 1; k <= 2; k++) {
    alpha[k-1] =
      -(get_m_ij(Ascr, 2*k-1, 2*k-1)*get_m_ij(Ascr, 2*k, 2*k-1)
	+ get_m_ij(Ascr, 2*k-1, 2*k)*get_m_ij(Ascr, 2*k, 2*k));
    beta[k-1] =
      sqr(get_m_ij(Ascr, 2*k-1, 2*k-1))
      + sqr(get_m_ij(Ascr, 2*k-1, 2*k));
  }
}


void LatticeType::Cell_Twiss(long i0, long i1, ss_vect<tps> &Ascr, bool chroma,
			     bool ring, double dP)
{
  long int            i;
  int                 k;
  std::vector<double> nu1, dnu;
  ss_vect<tps>        Ascr0, Ascr1;
  ElemType            *elemp;

  const int n = 4;

  for (k = 0; k < 2; k++) {
    nu1.push_back(0e0); dnu.push_back(0e0); 
  }

  if (conf.radiation) conf.dE = 0e0;

  elemp = elems[i0];
  dagetprm(Ascr, elemp->Alpha, elemp->Beta);
  elemp->Nu = nu1;

  Ascr0 = Ascr;
  for (k = 0; k < n+2; k++)
    Ascr0[k] += tps(conf.CODvect[k]);

  Ascr1 = Ascr0;
  for (i = i0; i <= i1; i++) {
    elems[i]->Elem_Pass(conf, Ascr1); elemp = elems[i];
    dagetprm(Ascr1, elemp->Alpha, elemp->Beta);
    for (k = 1; k <= 2; k++) {
      dnu[k-1] =
	(atan2(get_m_ij(Ascr1, 2*k-1, 2*k),
	       get_m_ij(Ascr1, 2*k-1, 2*k-1)) -
	 atan2(get_m_ij(Ascr0, 2*k-1, 2*k),
	       get_m_ij(Ascr0, 2*k-1, 2*k-1)))
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
	get_m_ij(Ascr1, k, 5)*get_m_ij(Ascr1, 6, 6) -
	get_m_ij(Ascr1, k, 6)*get_m_ij(Ascr1, 6, 5);
      elemp->Etap[k-1] =
	get_m_ij(Ascr1, 2*k, 5)*get_m_ij(Ascr1, 6, 6) -
	get_m_ij(Ascr1, 2*k, 6)*get_m_ij(Ascr1, 6, 5);
#else
      elemp->Eta[k-1] = get_m_ij(Ascr1, 2*k-1, 5);
      elemp->Etap[k-1] = get_m_ij(Ascr1, 2*k, 5);
#endif
    }
    Ascr0 = Ascr1;
  }

  if (chroma && !conf.Cavity_on) Cell_Geteta(i0, i1, ring, dP);
}


void LatticeType::Ring_Getchrom(double dP)
{
  long int
    lastpos;
  int
    k;
  std::vector<double>
    alpha = {0.0, 0.0},
    beta  = {0.0, 0.0},
    gamma = {0.0, 0.0},
    nu    = {0.0, 0.0},
    nu0   = {0.0, 0.0};
  
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
  long int
    lastpos = 0;
  int
    k,
    n = 0;
  std::vector<double>
    alpha = {0.0, 0.0},
    beta  = {0.0, 0.0},
    gamma = {0.0, 0.0},
    nu    = {0.0, 0.0};
  arma::mat
    M(tps_n, tps_n),
    R(tps_n, tps_n),
    A(tps_n, tps_n),
    Ainv(tps_n, tps_n);
  ss_vect<tps> Ascr;

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
  M = stlmattomat(conf.OneTurnMat);
  GDiag(n, elems[conf.Cell_nLoc]->S, A, Ainv, R, M, conf.Omega, conf.Alphac);
  conf.Ascr = mattostlmat(A);
  conf.Ascrinv = mattostlmat(Ainv);

  // Ascr = put_mat(n, conf.Ascr);
  Ascr = stlmattomap(conf.Ascr);
  if (!conf.Cavity_on) {
    // Ascr[delta_] = 0.0; Ascr[ct_] = 0.0;
    Ascr[delta_] = tps(0e0, delta_+1); Ascr[ct_] = 0e0;
  }

  Cell_Twiss(0, conf.Cell_nLoc, Ascr, chroma, true, dP);

  for (k = 0; k < 2; k++)
    conf.TotalTune[k] = elems[conf.Cell_nLoc]->Nu[k];
  conf.tuneflag = true;

  if (chroma && !conf.Cavity_on) {
    Ring_Getchrom(dP); GetCOD(conf.CODimax, conf.CODeps, dP, lastpos);
  }
}


void LatticeType::Ring_GetTwiss(bool chroma, double dP)
{

  if (conf.trace) printf("enter Ring_GetTwiss\n");
  Ring_Twiss(chroma, dP);
  conf.Alphac = conf.OneTurnMat[ct_][delta_]/elems[conf.Cell_nLoc]->S;
  if (conf.trace) printf("exit Ring_GetTwiss\n");
}


void LatticeType::TraceABN(long i0, long i1, const std::vector<double> &alpha,
			   const std::vector<double> &beta,
			   const std::vector<double> &eta,
			   const std::vector<double> &etap, const double dP)
{
  long         i;
  ss_vect<tps> Id, Ascr;

  Id.identity();
  Ascr.zero();
  for (i = 0; i < 2; i++) {
    Ascr[2*i] = sqrt(beta[i])*Id[2*i] + eta[i]*Id[delta_];
    Ascr[2*i+1] =
      -alpha[i]/sqrt(beta[i])*Id[2*i] + 1/sqrt(beta[i])*Id[2*i+1]
      + etap[i]*Id[delta_];
  }

  for (i = 0; i < ps_dim; i++)
    conf.CODvect[i] = 0e0;
  conf.CODvect[delta_] = dP;

  Ascr += stlvectops(conf.CODvect);
  Cell_Twiss(i0, i1, Ascr, false, false, dP);

  conf.Ascr = maptostlmat(Ascr);
}


void LatticeType::ttwiss(const std::vector<double> &alpha,
			 const std::vector<double> &beta,
			 const std::vector<double> &eta,
			 const std::vector<double> &etap,
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


void get_dI_eta_5( std::vector<ElemType*> elems, const int k)
{
  double       L, K, h, b2, alpha, beta, gamma, psi, eta, etap;
  ss_vect<tps> Id;
  MpoleType    *Mp;

  Id.identity();

  Mp = dynamic_cast<MpoleType*>(elems[k]);

  L = elems[k]->PL;
  h = Mp->Pirho;
  b2 = Mp->PBpar[Quad+HOMmax];
  K = b2 + sqr(Mp->Pirho);
  psi = sqrt(fabs(K))*L;
  alpha = elems[k-1]->Alpha[X_]; beta = elems[k-1]->Beta[X_];
  gamma = (1e0+sqr(alpha))/beta;
  eta = elems[k-1]->Eta[X_]; etap = elems[k-1]->Etap[X_];

  elems[k]->dI[1] += L*eta*h;
  elems[k]->dI[2] += L*sqr(h);
  elems[k]->dI[3] += L*fabs(cube(h));

  if (K > 0e0) {
    elems[k]->dI[4] +=
      h/K*(2e0*b2+sqr(h))
      *((eta*sqrt(K)*sin(psi)+etap*(1e0-cos(psi)))
	+ h/sqrt(K)*(psi-sin(psi)));

    elems[k]->dI[5] +=
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

    elems[k]->dI[4] +=
      h/K*(2e0*b2+sqr(h))
      *((eta*sqrt(K)*sinh(psi)-etap*(1e0-cosh(psi)))
	- h/sqrt(K)*(psi-sinh(psi)));

    elems[k]->dI[5] +=
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


void LatticeType::get_I(std::vector<double> &I, const bool prt)
{
  int j, k;

  for (k = 0; k <= 5; k++)
    I[k] = 0e0;

  if (prt)
    printf("\nget_I:\n"
	   "      name               s     curly_H            I_1"
	   "                     I_2                     I_3"
	   "                     I_4                     I_5"
	   "             alpha_x    beta_x     eta_x      eta'_x     alpha_y"
	   "    beta_y\n\n");
  for (j = 0; j <= conf.Cell_nLoc; j++)
    if ((elems[j]->Pkind == drift) || (elems[j]->Pkind == Mpole) ||
	(elems[j]->Pkind == Wigl) ||
	(elems[j]->Pkind == marker)) {
      if (prt)
	printf("%5d %-15s %6.3f %10.3e",
	       j, elems[j]->Name.c_str(), elems[j]->S,
	       elems[j]->curly_dH_x);
      for (k = 1; k <= 5; k++) {
	I[k] += elems[j]->dI[k];
	if (prt) printf(" %10.3e (%10.3e)", elems[j]->dI[k], I[k]);
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
      ps = mat_pass(Mp->M_elem, ps);

      if (conf.emittance && !conf.Cavity_on
	  && (elems[k]->PL != 0e0) && (Mp->Pirho != 0e0))
	get_dI_eta_5(elems, k);
    } else
      elems[k]->Elem_Pass(conf, ps);
  }
}


void LatticeType::get_eps_x(double &eps_x, double &sigma_delta, double &U_0,
			    std::vector<double> &J, std::vector<double> &tau,
			    std::vector<double> &I, const bool prt)
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
  A = stlmattomap(conf.Ascr); A += stlvectops(conf.CODvect);
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


double get_code(const ConfigType &conf, const ElemType &Cell)
{
  double    code;

  switch (Cell.Pkind) {
  case drift:
    code = 0.0;
    break;
  case Mpole:
    // Declare static to avoid compiler error.
    const static MpoleType *M = dynamic_cast<const MpoleType*>(&Cell);
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


void LatticeType::prt_lat1(const int loc1, const int loc2,
			   const std::string &fname, const bool all)
{
  long int i = 0;
  double   I5 = 0e0;
  FILE     *outf;

  outf = file_write(fname.c_str());
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
	      "%4ld %-10s %9.5f %4.1f %9.5f %8.5f %8.5f %8.5f %8.5f"
	      " %9.5f %8.5f %8.5f %8.5f %8.5f  %8.2e\n",
	      i, elems[i]->Name.c_str(), elems[i]->S, get_code(conf, *elems[i]),
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


void LatticeType::prt_lat2(const std::string &fname, const bool all)
{ prt_lat1(0, conf.Cell_nLoc, fname, all); }


void LatticeType::Cell_Twiss(const long int i0, const long int i1)
{
  long int            i;
  int                 k, nu_int[2];
  std::vector<double> alpha{0e0, 0e0}, beta{0e0, 0e0}, dnu{0e0, 0e0},
                      eta{0e0, 0e0}, etap{0e0, 0e0};
  ss_vect<tps>        A;

  for (k = 0; k < 2; k++)
    nu_int[k] = 0;

  for (i = i0; i <= i1; i++) {
    A = stlmattomap(elems[i]->A);
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


void LatticeType::prt_lat3(const int loc1, const int loc2,
			   const std::string &fname, const bool all,
			   const int n)
{
  long int            i = 0;
  int                 j, k;
  double              s, h;
  std::vector<double> alpha{0e0, 0e0}, beta{0e0, 0e0}, nu{0e0, 0e0},
		      dnu{0e0, 0e0}, eta{0e0, 0e0}, etap{0e0, 0e0},
		      dnu1{0e0, 0e0};
  double              curly_H;
  MpoleType           *Mp;
  ss_vect<double>     eta_Fl;
  ss_vect<tps>        A, A_CS;
  FILE                *outf;

  const double  c1 = 1e0/(2e0*(2e0-pow(2e0, 1e0/3e0))), c2 = 0.5e0-c1;
  const double  d1 = 2e0*c1, d2 = 1e0-2e0*d1;

  outf = file_write(fname.c_str());
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

	  fprintf(outf, "%4ld %-10s %6.2f %4.1f"
		  " %9.5f %8.5f %8.5f %11.8f %11.8f"
		  " %9.5f %8.5f %8.5f %8.5f %8.5f %10.3e %10.3e %10.3e\n",
		  i, elems[i]->Name.c_str(), s, get_code(conf, *elems[i]),
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

	fprintf(outf, "%4ld %-10s %6.2f %4.1f"
		" %9.5f %8.5f %8.5f %11.8f %11.8f"
		" %9.5f %8.5f %8.5f %8.5f %8.5f %10.3e %10.3e %10.3e\n",
		i, elems[i]->Name.c_str(), elems[i]->S,
		get_code(conf, *elems[i]), elems[i]->Alpha[X_],
		elems[i]->Beta[X_], elems[i]->Nu[X_], elems[i]->Eta[X_],
		elems[i]->Etap[X_], elems[i]->Alpha[Y_], elems[i]->Beta[Y_],
		elems[i]->Nu[Y_], elems[i]->Eta[Y_], elems[i]->Etap[Y_],
		eta_Fl[x_], eta_Fl[px_], curly_H);
      }
    }
  }

  fclose(outf);
}


void LatticeType::prt_lat4(const std::string &fname, const bool all,
			   const int n)
{ prt_lat3(0, conf.Cell_nLoc, fname, all, n); }


void LatticeType::prt_cod(const char *file_name, const bool all)
{
  long      i;
  double    b0, a0;
  FILE      *outf;
  struct tm *newtime;

  outf = file_write(file_name);

  /* Get time and date */
  newtime = GetTime();

  fprintf(outf,"# TRACY II v.2.6 -- %s -- %s \n",
	  file_name, asctime2(newtime));

  fprintf(outf, "#       name             s  code  betax   nux   betay   nuy"
	  "   xcod   ycod    dSx    dSy   dipx   dipy\n");
  fprintf(outf, "#                       [m]        [m]           [m]       "
	  "   [mm]   [mm]    [mm]   [mm] [mrad]  [mrad]\n");
  fprintf(outf, "#\n");

  for (i = 0; i <= conf.Cell_nLoc; i++) {
    if (all || (elems[i]->Fnum == conf.bpm)) {
      if (elems[i]->Pkind == Mpole) {
	b0 = Elem_GetKval(elems[i]->Fnum, elems[i]->Knum, Dip);
	a0 = Elem_GetKval(elems[i]->Fnum, elems[i]->Knum, -Dip);
      } else {
	b0 = 0e0;
	a0 = 0e0;
      }
      /* COD is in local coordinates */
      fprintf(outf,
	      "%4ld %-10s %6.2f %4.1f %6.3f %6.3f %6.3f %6.3f"
	      " %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",
	      i, elems[i]->Name.c_str(), elems[i]->S, get_code(conf, *elems[i]),
	      elems[i]->Beta[X_], elems[i]->Nu[X_],
	      elems[i]->Beta[Y_], elems[i]->Nu[Y_],
	      1e3*elems[i]->BeamPos[x_], 1e3*elems[i]->BeamPos[y_],
	      1e3*elems[i]->dS[X_], 1e3*elems[i]->dS[Y_], -1e3*b0, 1e3*a0);
    }
  }
  fclose(outf);
}


void LatticeType::prt_beampos(const char *file_name)
{
  long int k;
  FILE     *outf;

  outf = file_write(file_name);

  fprintf(outf, "#       name             s  code    xcod   ycod\n");
  fprintf(outf, "#                       [m]          [m]     [m]\n");
  fprintf(outf, "#\n");

  for (k = 0; k <= conf.Cell_nLoc; k++)
    fprintf(outf, "%4ld %-15s %6.2f %4.1f %12.5e %12.5e\n",
	    k, elems[k]->Name.c_str(), elems[k]->S,
	    get_code(conf, *elems[k]), elems[k]->BeamPos[x_],
	    elems[k]->BeamPos[y_]);

  fclose(outf);
}


void write_misalignments(const LatticeType &lat, const char* filename) {
  int       j;
  ElemType  *clp;
  MpoleType *M;
  FILE      *fp;

  fp = fopen(filename, "w");

  fprintf(fp, "# misalignment file");

  for (long int n = 0; n < lat.conf.Cell_nLoc; n++) {
    clp = lat.elems[n];
    fprintf(fp, "\n%li %i %.8e %.8e %.8e %.8e",
	    n, clp->Pkind, clp->dS[X_], clp->dS[Y_], clp->dT[X_],
	    clp->dT[Y_]); 
    
    switch(clp->Pkind) {
      case 2: // multipole
	M = dynamic_cast<MpoleType*>(clp);
         fprintf(fp, " %i", M->Porder);
        for (j = 0; j < 2*HOMmax+1; j++)
          fprintf(fp, " %.8e", M->PB[j]);
        break;
      default:
        break;  
    }
  }
  fclose(fp);
}


void LatticeType::prt_beamsizes(const int cnt)
{
  int   k;
  FILE  *fp;
  char fname [30];

  sprintf(fname,"%s_%d.out", beam_envelope_file, cnt);
  fp = file_write(fname);

  fprintf(fp,
	  "# k    name    s    s_xx    s_pxpx    s_xpx    s_yy    s_pypy"
	  "    s_ypy    theta_xy    s_xy\n");
  for(k = 0; k <= conf.Cell_nLoc; k++)
    fprintf(fp,"%4d %10s %e %e %e %e %e %e %e %e %e\n",
	    k, elems[k]->Name.c_str(), elems[k]->S,
	    elems[k]->sigma[x_][x_], elems[k]->sigma[px_][px_],
	    elems[k]->sigma[x_][px_],
	    elems[k]->sigma[y_][y_], elems[k]->sigma[py_][py_],
	    elems[k]->sigma[y_][py_],
	    radtodeg(atan2(2e0*elems[k]->sigma[x_][y_],
			   elems[k]->sigma[x_][x_]-elems[k]->sigma[y_][y_])
		     /2e0), elems[k]->sigma[x_][y_]);

  fclose(fp);

  sprintf(fname,"%s_%d.out","misalignments", cnt);
  write_misalignments(*this, fname);
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
	    i, elems[i]->Name.c_str(), elems[i]->S, get_code(conf, *elems[i]),
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


void LatticeType::Ring_Fittune(std::vector<double> &nu, double eps,
			       std::vector<int> &nq, long qf[],
			       long qd[], double dkL, long imax)
{
  struct LOC_Ring_Fittune V;

  int
    i,
    j,
    k;
  std::vector<double>
    nu0 = {0.0, 0.0},
    nu1 = {0.0, 0.0};
  ss_vect<double>
    dkL1,
    dnu;
  arma::mat
    A(tps_n, tps_n);

  const double dP = 0e0;

  if (setjmp(V._JL999)) return;

  if (trace)
    printf("  Tune fit, nux =%10.5f, nuy =%10.5f, eps =% .3E,"
	   " imax =%4ld, dkL = % .5E\n", nu[0], nu[1], eps, imax, dkL);
  Ring_GetTwiss(false, dP); checkifstable(&V);
  nu0 = conf.TotalTune;
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
    nu0 = conf.TotalTune;
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


void LatticeType::Ring_Fitchrom(std::vector<double> &ksi, double eps,
				std::vector<int> &ns, long sf[], long sd[],
				double dkpL, long imax)
{
  bool                rad;
  long int            lastpos;
  int                 i, j, k;
  std::vector<double> ksi0 = {0.0, 0.0};
  ss_vect<double>     dkpL1, dksi;
  arma::mat           A(tps_n, tps_n);

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


void findcod(LatticeType &lat, double dP)
{
  ss_vect<double> vcod;
  const int       ntrial = 40;  // maximum number of trials for closed orbit
  const double    tolx = 1e-10;  // numerical precision
  int             k, dim = 0;
  long            lastpos;

  // initializations
  for (k = 0; k <= 5; k++)
    vcod[k] = 0.0;  
    
  if (lat.conf.Cavity_on){
    fprintf(stdout,"warning looking for cod in 6D\n");
    dim = 6;
  } else{ // starting point linear closed orbit
    dim = 4;
    vcod[0] = lat.elems[0]->Eta[0]*dP; vcod[1] = lat.elems[0]->Etap[0]*dP;
    vcod[2] = lat.elems[0]->Eta[1]*dP; vcod[3] = lat.elems[0]->Etap[1]*dP;
    vcod[4] = dP;  // energy offset 
  }
  
  Newton_Raphson(lat, dim, vcod, ntrial, tolx);

  if (lat.conf.codflag == false)
    fprintf(stdout, "Error No COD found\n");
  
  lat.conf.CODvect = pstostlvec(vcod); // save closed orbit at the ring entrance

  if (lat.conf.trace) {
    fprintf(stdout,
       "Before cod2 % .5e % .5e % .5e % .5e % .5e % .5e \n",
       vcod[0], vcod[1], vcod[2], vcod[3], vcod[4], vcod[5]);
    lat.Cell_Pass(0, lat.conf.Cell_nLoc, vcod, lastpos);
    fprintf(stdout,
       "After  cod2 % .5e % .5e % .5e % .5e % .5e % .5e \n",
       vcod[0], vcod[1], vcod[2], vcod[3], vcod[4], vcod[5]);
  }
}


void prt_lin_map(const int n_DOF, const ss_vect<tps> &map)
{
  int i, j;

  std::cout << std::endl;
  for (i = 1; i <= 2*n_DOF; i++) {
    for (j = 1; j <= 2*n_DOF; j++)
      if (true)
	std::cout << std::scientific << std::setprecision(6)
	     << std::setw(14) << get_m_ij(map, i, j);
      else
	std::cout << std::scientific << std::setprecision(16)
	     << std::setw(24) << get_m_ij(map, i, j);
    std::cout << std::endl;
  }
}


ss_vect<tps> get_A(const std::vector<double> &alpha,
		   const std::vector<double> &beta,
		   const std::vector<double> &eta,
		   const std::vector<double> &etap)
{
  int          k;
  ss_vect<tps> A, Id;

  Id.identity();

  A.identity();
  for (k = 0; k < 2; k++) {
    A[2*k]  = sqrt(beta[k])*Id[2*k];
    A[2*k+1] = -alpha[k]/sqrt(beta[k])*Id[2*k] + 1.0/sqrt(beta[k])*Id[2*k+1];

    A[2*k] += eta[k]*Id[delta_]; A[2*k+1] += etap[k]*Id[delta_];
  }

  return A;
}


ss_vect<tps> get_A_CS(const int n, const ss_vect<tps> &A,
		      std::vector<double> &dnu)
{
  int          k;
  double       c, s;
  ss_vect<tps> Id, R;

  Id.identity(); R.identity(); get_dnu(n, A, dnu);
  for (k = 0; k < n; k++) {
    c = cos(2.0*M_PI*dnu[k]); s = sin(2.0*M_PI*dnu[k]);
    R[2*k] = c*Id[2*k] - s*Id[2*k+1]; R[2*k+1] = s*Id[2*k] + c*Id[2*k+1];
  }

  return A*R;
}


void get_ab(const ss_vect<tps> &A, std::vector<double> &alpha,
	    std::vector<double> &beta, std::vector<double> &dnu,
	    std::vector<double> &eta, std::vector<double> &etap)
{
  int          k;
  ss_vect<tps> A_Atp;

  A_Atp = A*tp_S(2, A);

  for (k = 0; k <= 1; k++) {
    eta[k] = A[2*k][delta_]; etap[k] = A[2*k+1][delta_];

    alpha[k] = -A_Atp[2*k][2*k+1]; beta[k] = A_Atp[2*k][2*k];
  }

  get_dnu(2, A, dnu);
}


void get_twoJ(const int n_DOF, const ss_vect<double> &ps, const ss_vect<tps> &A,
	      double twoJ[])
{
  int             j, no;
  long int        jj[ps_dim];
  ss_vect<double> z;

  no = no_tps;
  danot_(1);

  for (j = 0; j < ps_dim; j++)
    jj[j] = (j < 2*n_DOF)? 1 : 0;

  z = (PInv(A, jj)*ps).cst();

  for (j = 0; j < n_DOF; j++)
    twoJ[j] = sqr(z[2*j]) + sqr(z[2*j+1]);

  danot_(no);
}


void SetKLpar(LatticeType &lat, int Fnum, int Knum, int Order, double kL)
{
  long int  loc;
  MpoleType *M;

  loc = lat.Elem_GetPos(Fnum, Knum);
  M = dynamic_cast<MpoleType*>(lat.elems[loc]);
  if (lat.elems[loc]->PL != 0e0)
    M->PBpar[Order+HOMmax] = kL/lat.elems[loc]->PL;
  else
    M->PBpar[Order+HOMmax] = kL;
  lat.SetPB(Fnum, Knum, Order);
}


double GetKpar(LatticeType &lat, int Fnum, int Knum, int Order)
{
  long int  loc;
  MpoleType *M;


  loc = lat.Elem_GetPos(Fnum, Knum);
  M = dynamic_cast<MpoleType*>(lat.elems[loc]);
  return M->PBpar[Order+HOMmax];
}


void Trac(LatticeType &lat, double x, double px, double y, double py, double dp,
	  double ctau, long nmax, long pos, long &lastn, long &lastpos,
	  FILE *outf1)
{
  /* Compute closed orbit : usefull if insertion devices */

  ss_vect<double> x1;     /* tracking coordinates */

  x1[0] = x; x1[1] = px;
  x1[2] = y; x1[3] = py;
  x1[4] =dp; x1[5] = ctau;

  lastn = 0L;

  (lastpos)=pos;
  if(lat.conf.trace) fprintf(outf1, "\n");
  fprintf(outf1, "%6ld %+10.5e %+10.5e %+10.5e %+10.5e %+10.5e %+10.5e \n",
	  lastn, x1[0], x1[1], x1[2], x1[3], x1[4], x1[5]);
  lat.Cell_Pass(pos -1L, lat.conf.Cell_nLoc, x1, lastpos);

  if (lastpos == lat.conf.Cell_nLoc)
  do {
    (lastn)++;
    lat.Cell_Pass(0L, pos-1L, x1, lastpos);
    if(!lat.conf.trace) {
      fprintf(outf1, "%6ld %+10.5e %+10.5e %+10.5e %+10.5e %+10.5e %+10.5e \n",
	      lastn, x1[0], x1[1], x1[2], x1[3], x1[4], x1[5]);
    }
    if (lastpos == pos-1L)
      lat.Cell_Pass(pos-1L,lat.conf.Cell_nLoc, x1, lastpos);
  }
  while (((lastn) < nmax) && ((lastpos) == lat.conf.Cell_nLoc));

  if (lastpos == lat.conf.Cell_nLoc) lat.Cell_Pass(0L,pos, x1, lastpos);

  if (lastpos != pos) {
    printf("Trac: Particle lost \n");
    fprintf(stdout, "turn:%6ld plane: %1d"
	    " %+10.5g %+10.5g %+10.5g %+10.5g %+10.5g %+10.5g \n", 
	    lastn, lat.conf.lossplane, x1[0], x1[1], x1[2], x1[3], x1[4],
	    x1[5]);
  }
}


void LatticeType::print(const string &str)
{
  int j, k;

  printf("\n***************************************************************"
	 "***************\n");
  printf("\n");
  printf("  dPcommon     =  %9.3e  dPparticle   =  %9.3e"
	 "  Energy [GeV] = %.3f\n",
         conf.dPcommon, conf.dPparticle, conf.Energy);
  printf("  MaxAmplx [m] = %9.3e  MaxAmply [m] = %9.3e"
	 "  RFAccept [%%] = %4.2f\n",
         elems[0]->maxampl[X_][0], elems[0]->maxampl[Y_][0], conf.delta_RF*1e2);
  printf(" Cavity_On    =  %s    ", conf.Cavity_on ? "TRUE " : "FALSE");
  printf("  Radiation_On = %s     \n", conf.radiation ? "TRUE " : "FALSE");
  printf("  bpm          =  %3d        qt           = %3d        ",
	 conf.bpm, conf.qt);
  printf(" Chambre_On   = %s     \n", conf.chambre ? "TRUE " : "FALSE");
  printf("  hcorr        =  %3d        vcorr        = %3d\n\n",
	 conf.hcorr, conf.vcorr);
  printf("  alphac       =   %22.16e\n", conf.Alphac); 
  printf("  nux          =  %19.16f      nuz  =  %19.16f",
         conf.TotalTune[X_], conf.TotalTune[Y_]);
  if (conf.Cavity_on)
    printf("  omega  = %11.9f\n", conf.Omega);
  else {
    printf("\n");
    printf("  ksix         = %10.6f                ksiz = %10.6f\n",
	   conf.Chrom[X_], conf.Chrom[Y_]);
  }

  printf("\n  OneTurn matrix:\n");
  for (j = 0; j < ps_dim; j++) {
    for (k = 0; k < ps_dim; k++)
      printf("%14.6e", conf.OneTurnMat[j][k]);
    printf("\n");
  }
  fflush(stdout);
}


void LatticeType::GetEmittance(const int Fnum, const bool prt)
{
  // A. Chao "Evaluation of Beam Distribution Parameters in an Electron
  // Storage Ring" J. Appl. Phys 50 (2), 595-598.
  bool                emit, rad, cav, path;
  int                 i, j, h_RF;
  long int            lastpos, loc;
  double              C, theta, V_RF, phi0, gamma_z;
  double              sigma_s, sigma_delta;
  std::vector<double> nu = {0.0, 0.0, 0.0};
  arma::mat           Ascr(tps_n, tps_n), sigma(nv_tps, nv_tps);
  ss_vect<tps>        Ascr_map;
  CavityType          *Cp;

  // save state
  rad = conf.radiation; emit = conf.emittance;
  cav = conf.Cavity_on; path = conf.pathlength;

  C = elems[conf.Cell_nLoc]->S;

  // damped system
  conf.radiation = true; conf.emittance  = true;
  conf.Cavity_on = true; conf.pathlength = false;

  Ring_GetTwiss(false, 0.0);

  // radiation loss is computed in Cav_Pass

  loc = Elem_GetPos(Fnum, 1);
  Cp = dynamic_cast<CavityType*>(elems[loc]);

  conf.U0 = conf.dE*1e9*conf.Energy;
  V_RF = Cp->Pvolt;
  h_RF = Cp->Ph;
  phi0 = fabs(asin(conf.U0/V_RF));
  conf.delta_RF =
    sqrt(-V_RF*cos(M_PI-phi0)*(2.0-(M_PI-2.0*(M_PI-phi0))
    *tan(M_PI-phi0))/(fabs(conf.Alphac)*M_PI*h_RF*1e9*conf.Energy));

  // Compute diffusion coeffs. for eigenvectors [sigma_xx, sigma_yy, sigma_zz]
  Ascr_map = stlmattomap(conf.Ascr); Ascr_map += stlvectops(conf.CODvect);

  Cell_Pass(0, conf.Cell_nLoc, Ascr_map, lastpos);

  // K. Robinson "Radiation Effects in Circular Electron Accelerators"
  // Phys. Rev. 111 (2), 373-380.
  // Iu.F. Orlov, E.K. Tarasov "Damping of Oscillations in a Cyclic Electron
  // Accelerator" J. Exptl. Theoret. Phys. 34, 651-657 (1958).
  // To leading order:
  //   Sum_k(alpha_k) = 2*U_0/E
  // or
  //   J_x + J_y + J_z = 4.

  for (i = 0; i < ps_dim/2; i++) {
    // partition numbers
    conf.J[i] =
      2.0*(1.0+conf.CODvect[delta_])*conf.alpha_rad[i]/conf.dE;
    // damping times
    conf.tau[i] = -C/(c0*conf.alpha_rad[i]);
    // diffusion coeff. and emittance (alpha is for betatron amplitudes)
    conf.eps[i] = -conf.D_rad[i]/(2.0*conf.alpha_rad[i]);
    // fractional tunes
    nu[i]  = atan2(conf.wi[i*2], conf.wr[i*2])/(2.0*M_PI);
    if (nu[i] < 0.0) nu[i] = 1.0 + nu[i];
  }

  // undamped system
  conf.radiation = !false; conf.emittance = false;

  Ring_GetTwiss(false, 0.0);

  // Compute sigmas arround the lattice:
  //   Sigma = A diag[J_1, J_1, J_2, J_2, J_3, J_3] A^T
  for (i = 0; i < 6; i++) {
    Ascr_map[i] = tps(conf.CODvect[i]);
    for (j = 0; j < 6; j++)
      Ascr_map[i] += conf.Ascr[i][j]*sqrt(conf.eps[j/2])*tps(0.0, j+1);
  }
  // prt_lin_map(3, Ascr_map);
  for (loc = 0; loc <= conf.Cell_nLoc; loc++) {
    elems[loc]->Elem_Pass(conf, Ascr_map);
    // sigma = A x A^tp
    elems[loc]->sigma = mattostlmat(Ascr*trans(Ascr));
  }

  // A. W. Chao, M. J. Lee "Particle Distribution Parameters in an Electron
  // Storage Ring" J. Appl. Phys. 47 (10), 4453-4456 (1976).
  // observable tilt angle
  theta =
    atan2(2e0*elems[0]->sigma[x_][y_],
	  (elems[0]->sigma[x_][x_]-elems[0]->sigma[y_][y_]))/2e0;

  // longitudinal alpha and beta
  conf.alpha_z =
    -conf.Ascr[ct_][ct_]*conf.Ascr[delta_][ct_]
    - conf.Ascr[ct_][delta_]*conf.Ascr[delta_][delta_];
  conf.beta_z =
    sqr(conf.Ascr[ct_][ct_]) + sqr(conf.Ascr[ct_][delta_]);
  gamma_z = (1.0+sqr(conf.alpha_z))/conf.beta_z;

  // bunch size
  sigma_s = sqrt(conf.beta_z*conf.eps[Z_]);
  sigma_delta = sqrt(gamma_z*conf.eps[Z_]);

  if (prt) {
    printf("\nEmittance:\n");
    printf("\nBeam energy [GeV]:              "
	   "Eb          = %4.2f\n", conf.Energy);
    printf("Energy loss per turn [keV]:     "
	   "U0          = %3.1f\n",
	   1e-3*conf.U0);
    printf("Synchronous phase [deg]:        "
	   "phi0        = 180 - %4.2f\n",
	   radtodeg(phi0));
    printf("RF bucket height [%%]:           "
	   "delta_RF    = %4.2f\n", 1e2*conf.delta_RF);
    printf("\nEquilibrium emittance [m.rad]:  "
	   "eps         =  [%9.3e, %9.3e, %9.3e]\n",
            conf.eps[X_], conf.eps[Y_], conf.eps[Z_]);
    printf("Bunch length [mm]:              "
	   "sigma_s     =  %5.3f\n", 1e3*sigma_s);
    printf("Momentum spread:                "
	   "sigma_delta =  %9.3e\n", sigma_delta);
    printf("Longitudinal phase space:       "
	   "alpha_z     = %10.3e, beta_z = %9.3e\n",
	   conf.alpha_z, conf.beta_z);
    printf("Partition numbers:              "
	   "J           =  [%5.3f, %5.3f, %5.3f]\n",
            conf.J[X_], conf.J[Y_], conf.J[Z_]);
    printf("Damping times [msec]:           "
	   "tau         =  [%3.1f, %3.1f, %3.1f]\n",
	   1e3*conf.tau[X_], 1e3*conf.tau[Y_], 1e3*conf.tau[Z_]);
    printf("Diffusion coeffs:               "
	   "D           =  [%7.1e, %7.1e, %7.1e]\n",
	   conf.D_rad[X_], conf.D_rad[Y_], conf.D_rad[Z_]);
    printf("\nalphac:                         "
	   "alphac      = %11.3e\n", conf.Alphac);
    printf("\nFractional tunes:               "
	   "nu_x        =  [%7.5f, %7.5f, %12.5e]\n", nu[X_], nu[Y_], nu[Z_]);
    printf("                                "
	   "1-nu_x      =  [%7.5f, %7.5f, %12.5e]\n",
	   1e0-nu[X_], 1e0-nu[Y_], 1e0-nu[Z_]);
    printf("\nsigmas:                         "
	   "sigma_x     =  %5.1f  microns, sigma_px    = %5.1f urad\n",
	   1e6*sqrt(elems[0]->sigma[x_][x_]),
	   1e6*sqrt(elems[0]->sigma[px_][px_]));
    printf("                                "
	   "sigma_y     =  %5.1f  microns, sigma_py    = %5.1f urad\n",
	   1e6*sqrt(elems[0]->sigma[y_][y_]),
	   1e6*sqrt(elems[0]->sigma[py_][py_]));
    printf("                                "
	   "sigma_s     =  %6.2f mm,      sigma_delta = %8.2e\n",
	   1e3*sqrt(elems[0]->sigma[ct_][ct_]),
	   sqrt(elems[0]->sigma[delta_][delta_]));

    printf("\nBeam ellipse twist [rad]:       tw      = %5.3f\n", theta);
    printf("                   [deg]:       tw      = %5.3f\n",
	   radtodeg(theta));
  }

  // restore state
  conf.radiation = rad; conf.emittance  = emit;
  conf.Cavity_on = cav; conf.pathlength = path;
}


double get_dynap(LatticeType &lat, const double delta, const int n_aper,
		 const int n_track, const bool cod)
{
  char   str[max_str];
  int    i;
  double x_aper[n_aper], y_aper[n_aper], DA;
  FILE   *fp;

  const int prt = true;

  fp = file_write("dynap.out");
  dynap(fp, lat, 5e-3, 0.0, 0.1e-3, n_aper, n_track, x_aper, y_aper, false, cod,
	prt);
  fclose(fp);
  DA = get_aper(n_aper, x_aper, y_aper);

  if (true) {
    sprintf(str, "dynap_dp%3.1f.out", 1e2*delta);
    fp = file_write(str);
    dynap(fp, lat, 5e-3, delta, 0.1e-3, n_aper, n_track,
      x_aper, y_aper, false, cod, prt);
    fclose(fp);
    DA += get_aper(n_aper, x_aper, y_aper);

    for (i = 0; i < ps_dim; i++)
      lat.conf.CODvect[i] = 0.0;
    sprintf(str, "dynap_dp%3.1f.out", -1e2*delta);
    fp = file_write(str);
    dynap(fp, lat, 5e-3, -delta, 0.1e-3, n_aper,
      n_track, x_aper, y_aper, false, cod, prt);
    fclose(fp);
    DA += get_aper(n_aper, x_aper, y_aper);
  }

  return DA/3.0;
}


struct tm* GetTime()
{
  struct tm *whattime;
  /* Get time and date */
  time_t aclock;
  time(&aclock);                 /* Get time in seconds */
  whattime = localtime(&aclock);  /* Convert time to struct */
  return whattime;
}


void computeFandJ(LatticeType &lat, int n, ss_vect<double> &x, arma::mat &fjac,
		  ss_vect<double> &fvect)
{
  int             i, k;
  long            lastpos = 0;
  ss_vect<double> x0, fx1, fx2;

  //stepsize for numerical differentiation
  const double deps = 1e-8;

  x0 = x;
  
  lat.Cell_Pass(0, lat.conf.Cell_nLoc, x0, lastpos);
  fvect = x0;
  
  // compute Jacobian matrix by numerical differentiation
  for (k = 0; k < n; k++) {
    x0 = x;
    x0[k] += deps;  // differential step in coordinate k

    // tracking along the ring
    lat.Cell_Pass(0, lat.conf.Cell_nLoc, x0, lastpos);
    fx1 = x0;

    x0 = x;
    // differential step in coordinate k
    x0[k] -= deps;

    // tracking along the ring
    lat.Cell_Pass(0, lat.conf.Cell_nLoc, x0, lastpos);
    fx2 = x0;

    // symmetric difference formula
    for (i = 0; i < n; i++)
      fjac(i, k) = 0.5 * (fx1[i] - fx2[i]) / deps;
  }
}


int Newton_Raphson(LatticeType &lat, int n, ss_vect<double> &x, int ntrial,
		   double tolx)
{
  int             k, i;
  double          errx;
  ss_vect<double> bet, fvect;
  arma::mat       alpha   = arma::mat(tps_n, tps_n);
  arma::vec       bet_vec = arma::vec(tps_n);

  errx = 0.0;

  for (k = 1; k <= ntrial; k++) {    // loop over number of iterations
    // supply function values at x in fvect and Jacobian matrix in fjac
    computeFandJ(lat, n, x, alpha, fvect);

    // Jacobian - Id
    for (i = 0; i < n; i++)
      alpha(i, i) -= 1.0;
    for (i = 0; i < n; i++)
      bet[i] = x[i] - fvect[i];      // right side of linear equation
    // inverse matrix using gauss jordan method from Tracy (from NR)
    alpha = inv(alpha);
    bet_vec = alpha*pstovec(bet);
    bet = vectops(bet_vec);         // bet = alpha*bet
    errx = 0.0;                      // check root convergence
    for (i = 0; i < n; i++) {
      // update solution
      errx += fabs(bet[i]);
      x[i] += bet[i]; 
    }
    
    if (lat.conf.trace)
      fprintf(stdout,
         "%02d: cod2 % .5e % .5e % .5e % .5e % .5e % .5e  errx =% .5e\n",
         k, x[0], x[1], x[2], x[3], x[4], x[5], errx);
    if (errx <= tolx)
    {
      lat.conf.codflag = true;
        return 1;
    }
  }
  // check whever closed orbit found out
  if ((k >= ntrial) && (errx >= tolx))
  {
    lat.conf.codflag = false;
      return 1;
  }
  return 0;
}


void get_dnu(const int n, const ss_vect<tps> &A, std::vector<double> &dnu)
{
  int k;

  const double eps = 1e-15;

  for (k = 0; k < n; k++) {
    if (k < 2)
      dnu[k] = atan2(A[2*k][2*k+1], A[2*k][2*k])/(2e0*M_PI);
    else
      dnu[k] = -atan2(A[ct_][delta_], A[ct_][ct_])/(2e0*M_PI);
    if (dnu[k] < -eps) dnu[k] += 1e0;
  }
}


ss_vect<tps> tp_S(const int n_DOF, const ss_vect<tps> &A)
{
  int          j;
  long int     jj[ps_dim];
  ss_vect<tps> S;

  for (j = 1; j <= ps_dim; j++)
    jj[j-1] = (j <= 2*n_DOF)? 1 : 0;

  S = get_S(n_DOF);

  return -S*PInv(A, jj)*S;
}


void dynap(FILE *fp, LatticeType &lat, double r, const double delta,
	   const double eps, const int npoint, const int nturn, double x[],
	   double y[], const bool floqs, const bool cod, const bool print)

{
  /* Determine the dynamical aperture by tracking.
     Assumes mid-plane symmetry.                    */

  long int  i, lastpos;
  double    phi, x_min, x_max, y_min, y_max;

  if (cod)
    lat.getcod(delta, lastpos);
  else
    lat.conf.CODvect.assign(lat.conf.CODvect.size(), 0e0);
  if (floqs) {
    lat.Ring_GetTwiss(false, delta);
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
    getdynap(lat, r, phi, delta, eps, nturn, floqs);
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


ss_vect<tps> get_S(const int n_DOF)
{
  int          j;
  ss_vect<tps> S;

  S.zero();
  for (j = 0; j < n_DOF; j++) {
    S[2*j] = tps(0.0, 2*j+2); S[2*j+1] = -tps(0.0, 2*j+1);
  }

  return S;
}


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

    sprintf(result, "%.3s %.3s%3d %.2d:%.2d:%.2d %u",
	    wday_name[timeptr->tm_wday],
	    mon_name[timeptr->tm_mon],
	    timeptr->tm_mday, timeptr->tm_hour,
	    timeptr->tm_min, timeptr->tm_sec,
	    1900+timeptr->tm_year);
    return result;
}


#define nfloq     4
bool chk_if_lost(LatticeType &lat, double x0, double y0, double delta,
		 long int nturn, bool floqs)
{
  long int        i, lastn, lastpos;
  ss_vect<double> ps;
  arma::vec       ps_vec = arma::vec(tps_n);

  bool prt = false;

  ps[x_] = x0; ps[px_] = px_0; ps[y_] = y0; ps[py_] = py_0;
  ps[delta_] = delta; ps[ct_] = 0.0;
  if (floqs) {
    // Transform to phase space.
    ps_vec = stlmattomat(lat.conf.Ascr)*pstovec(ps);
    ps = vectops(ps_vec);
  }
  for (i = 0; i <= 3; i++)
    ps[i] += lat.conf.CODvect[i];

  lastn = 0;
  if (prt) {
    printf("\n");
    printf("chk_if_lost:\n");
    printf("%4ld %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
	   lastn, 1e3*ps[x_], 1e3*ps[px_], 1e3*ps[y_], 1e3*ps[py_],
	   1e2*ps[delta_], 1e3*ps[ct_]);
  }
  do {
    lastn++;
    lat.Cell_Pass(0, lat.conf.Cell_nLoc, ps, lastpos);
    if (prt)
      printf("%4ld %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
	     lastn, 1e3*ps[x_], 1e3*ps[px_], 1e3*ps[y_], 1e3*ps[py_],
	     1e2*ps[delta_], 1e3*ps[ct_]);
  } while ((lastn != nturn) && (lastpos == lat.conf.Cell_nLoc));
  return(lastn != nturn);
}
#undef nfloq


void getdynap(LatticeType &lat, double &r, double phi, double delta, double eps,
	      int nturn, bool floqs)
{
  /* Determine dynamical aperture by binary search. */
  double  rmin = 0.0, rmax = r;

  const bool    prt   = false;
  const double  r_reset = 1e-3, r0 = 10e-3;

  if (prt) printf("\n");

  while (!chk_if_lost(lat, rmax*cos(phi), rmax*sin(phi), delta, nturn, floqs)) {
    if (rmax < r_reset) rmax = r0;
    rmax *= 2.0;
  }
  while (rmax-rmin >= eps) {
    r = rmin + (rmax-rmin)/2.0;
    if (prt) printf("getdynap: %6.3f %6.3f %6.3f\n",
		    1e3*rmin, 1e3*rmax, 1e3*r);
    if (!chk_if_lost(lat, r*cos(phi), r*sin(phi), delta, nturn, floqs) )
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
