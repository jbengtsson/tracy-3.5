/* NSLS-II specific library

   J. Bengtsson  NSLS-II, BNL  2004 -

   T. Shaftan, I. Pinayev, Y. Luo, C. Montag, B. Nash

*/


// global params

const bool normal = true; // Normal or rectangular distribution

int    Fnum_Cart, n_iter_Cart;

double u_Touschek;  // argument for Touschek D(ksi)
double chi_m;       // argument for IBS D(ksi)

// IBS (Bjorken-Mtingwa)
double a_IBS, b_IBS, c_IBS, a_k_IBS, b_k_IBS;

// for IBS
int    i_, j_;
double **C_;

ss_vect<tps> map;
MNF_struct   MNF;


// conversion

void lwr_case(char str[])
{
  int k;

  for (k = 0; k < (int)strlen(str); k++)
    str[k] = tolower(str[k]);
}


void upr_case(char str[])
{
  int k;

  for (k = 0; k < (int)strlen(str); k++)
    str[k] = toupper(str[k]);
}


// only supported by Redhat
#if 0
// generate backtrace
void prt_trace (void)
{
  const int  max_entries = 20;

  void    *array[max_entries];
  size_t  size;
  char    **strings;
  size_t  i;

  size = backtrace(array, max_entries);
  strings = backtrace_symbols(array, size);

  printf("prt_trace: obtained %zd stack frames\n", size);

  for (i = 0; i < size; i++)
    printf ("%s\n", strings[i]);

  free (strings);
}
#endif


// file I/O

// C++

void file_rd(std::ifstream &inf, const string &file_name)
{

  inf.open(file_name.c_str(), std::ios::in);
  if (!inf.is_open()) {
    cout << "File not found: " << file_name << "\n";
    exit_(-1);
  }
}


void file_wr(std::ofstream &outf, const string &file_name)
{

  outf.open(file_name.c_str(), std::ios::out);
  if (!outf.is_open()) {
    cout << "Could not create file: " << file_name << "\n";
    exit_(-1);
  }
}


void file_rd(std::ifstream &inf, const char file_name[])
{

  inf.open(file_name, std::ios::in);
  if (!inf.is_open()) {
    printf("File not found: %s\n", file_name);
    exit_(-1);
  }
}


void file_wr(std::ofstream &outf, const char file_name[])
{

  outf.open(file_name, std::ios::out);
  if (!outf.is_open()) {
    printf("Could not create file: %s\n", file_name);
    exit_(-1);
  }
}


// C

FILE* file_read(const char file_name[])
{
  FILE      *fp;

  fp = fopen(file_name, "r");
  if (fp == NULL) {
    printf("File not found: %s\n", file_name);
    exit_(-1);
  }
  return(fp);
}


FILE* file_write(const char file_name[])
{
  FILE      *fp;

  fp = fopen(file_name, "w");
  if (fp == NULL) {
    printf("Could not create file: %s\n", file_name);
    exit_(-1);
  }
  return(fp);
}


void chk_cod(const bool cod, const char *proc_name)
{

  if (!cod) {
    printf("%s: closed orbit not found\n", proc_name);
//     exit_(1);
  }
}


void no_sxt(void)
{
  int       k;
  MpoleType *M;

  std::cout << std::endl;
  std::cout << "zeroing sextupoles" << std::endl;
  for (k = 0; k <= Lattice.param.Cell_nLoc; k++) {
    M = static_cast<MpoleType*>(&Lattice.Cell[k]);
    if ((Lattice.Cell[k].Kind == Mpole) && (M->order >= Sext))
      SetKpar(Lattice.Cell[k].Fnum, Lattice.Cell[k].Knum, Sext, 0.0);
  }
}


void get_map(const bool cod)
{
  long int  lastpos;

  map.identity();
  if (cod) {
    Lattice.getcod(0e0, lastpos);
    map += Lattice.param.CODvect;
  }
  Lattice.Cell_Pass(0, Lattice.param.Cell_nLoc, map, lastpos);
  if (cod) map -= Lattice.param.CODvect;
}


#if NO > 1

tps get_h(void)
{
  ss_vect<tps>  map1, R;

  // Parallel transport nonlinear kick to start of lattice,
  // assumes left to right evaluation.

  if (true)
    // Dragt-Finn factorization
    return LieFact_DF(Inv(MNF.A0*MNF.A1)*map*MNF.A0*MNF.A1, R)*R;
  else {
    // single Lie exponent
    danot_(1); map1 = map; danot_(no_tps);
    return LieFact(Inv(MNF.A0*MNF.A1)*map*Inv(map1)*MNF.A0*MNF.A1);
  }
}

#endif

void get_m2(const ss_vect<tps> &ps, tps m2[])
{
  int  i, j, k;

  k = 0;
  for (i = 0; i < 2*nd_tps; i++)
    for (j = i; j < 2*nd_tps; j++) {
      k++; m2[k-1] = ps[i]*ps[j];
    }
}


ss_vect<tps> get_S(const int n_DOF)
{
  int           j;
  ss_vect<tps>  S;

  S.zero();
  for (j = 0; j < n_DOF; j++) {
    S[2*j] = tps(0.0, 2*j+2); S[2*j+1] = -tps(0.0, 2*j+1);
  }

  return S;
}


ss_vect<tps> tp_S(const int n_DOF, const ss_vect<tps> &A)
{
  int           j, jj[ss_dim];
  ss_vect<tps>  S;

  for (j = 1; j <= ss_dim; j++)
    jj[j-1] = (j <= 2*n_DOF)? 1 : 0;

  S = get_S(n_DOF);

  return -S*PInv(A, jj)*S;
}


void get_dnu(const int n, const ss_vect<tps> &A, double dnu[])
{
  int  k;

  for (k = 0; k < n; k++) {
    dnu[k] = atan2(A[2*k][2*k+1], A[2*k][2*k])/(2.0*M_PI);
    if (dnu[k] < 0.0) dnu[k] += 1.0;
  }
}


void get_ab(const ss_vect<tps> &A,
	    double alpha[], double beta[], double dnu[],
	    double eta[], double etap[])
{
  int           k;
  ss_vect<tps>  A_Atp;

  A_Atp = A*tp_S(2, A);

  for (k = 0; k <= 1; k++) {
    eta[k] = A[2*k][delta_]; etap[k] = A[2*k+1][delta_];

    alpha[k] = -A_Atp[2*k][2*k+1]; beta[k] = A_Atp[2*k][2*k];
  }

  get_dnu(2, A, dnu);
}


ss_vect<tps> get_A(const double alpha[], const double beta[],
		   const double eta[], const double etap[])
{
  int           k;
  ss_vect<tps>  A, Id;

  Id.identity();

  A.identity();
  for (k = 0; k < 2; k++) {
    A[2*k]  = sqrt(beta[k])*Id[2*k];
    A[2*k+1] = -alpha[k]/sqrt(beta[k])*Id[2*k] + 1.0/sqrt(beta[k])*Id[2*k+1];

    A[2*k] += eta[k]*Id[delta_]; A[2*k+1] += etap[k]*Id[delta_];
  }

  return A;
}


ss_vect<tps> get_A_CS(const int n, const ss_vect<tps> &A, double dnu[])
{
  int           k;
  double        c, s;
  ss_vect<tps>  Id, R;

  Id.identity(); R.identity(); get_dnu(n, A, dnu);
  for (k = 0; k < n; k++) {
    c = cos(2.0*M_PI*dnu[k]); s = sin(2.0*M_PI*dnu[k]);
    R[2*k] = c*Id[2*k] - s*Id[2*k+1]; R[2*k+1] = s*Id[2*k] + c*Id[2*k+1];
  }

  return A*R;
}


void prt_lin_map(const int n_DOF, const ss_vect<tps> &map)
{
  int i, j;

  std::cout << std::endl;
  for (i = 1; i <= 2*n_DOF; i++) {
    for (j = 1; j <= 2*n_DOF; j++)
      if (true)
	std::cout << std::scientific << std::setprecision(6)
	     << std::setw(14) << getmat(map, i, j);
      else
	std::cout << std::scientific << std::setprecision(16)
	     << std::setw(24) << getmat(map, i, j);
    std::cout << std::endl;
  }
}


void get_twoJ(const int n_DOF, const ss_vect<double> &ps,
	      const ss_vect<tps> &A, double twoJ[])
{
  int              j;
  iVector          jj;
  ss_vect<double>  z;

  for (j = 0; j < nv_tps; j++)
    jj[j] = (j < 2*n_DOF)? 1 : 0;

  z = (PInv(A, jj)*ps).cst();

  for (j = 0; j < n_DOF; j++)
    twoJ[j] = sqr(z[2*j]) + sqr(z[2*j+1]);
}


double get_curly_H(const double alpha_x, const double beta_x,
		   const double eta_x, const double etap_x)
{
  double  curly_H, gamma_x;

  gamma_x = (1.0+sqr(alpha_x))/beta_x;

  curly_H = gamma_x*sqr(eta_x) + 2.0*alpha_x*eta_x*etap_x + beta_x*sqr(etap_x);

  return curly_H;
}


double LatticeType::get_eps_x(void)
{
  bool         cav, emit;
  long int     lastpos;
  double       eps_x;
  ss_vect<tps> A;

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

  cav = Lattice.param.Cavity_on; emit = Lattice.param.emittance;

  Lattice.param.Cavity_on = false; Lattice.param.emittance = false;

  Lattice.Ring_GetTwiss(false, 0.0);

  putlinmat(6, Lattice.param.Ascr, A); A += Lattice.param.CODvect;

  Lattice.param.emittance = true;

  Lattice.Cell_Pass(0, Lattice.param.Cell_nLoc, A, lastpos);

  eps_x = 1470.0*pow(Lattice.param.Energy, 2.0)*I5/(I2-I4);

  printf("\n");
  printf("eps_x = %5.3f nm.rad\n", eps_x);
  printf("J_x   = %5.3f, J_z = %5.3f\n", 1.0-I4/I2, 2.0+I4/I2);

  Lattice.param.Cavity_on = cav; Lattice.param.emittance = emit;

  return eps_x;
}


void LatticeType::GetEmittance(const int Fnum, const bool prt)
{
  // A. Chao "Evaluation of Beam Distribution Parameters in an Electron
  // Storage Ring" J. Appl. Phys 50 (2), 595-598.
  bool         emit, rad, cav, path;
  int          i, j, h_RF;
  long int     lastpos, loc;
  double       circ, theta, V_RF, phi0, gamma_z;
  double       sigma_s, sigma_delta;
  Vector3      nu;
  Matrix       Ascr;
  ss_vect<tps> Ascr_map;
  CavityType   *C;

  // save state
  rad = Lattice.param.radiation; emit = Lattice.param.emittance;
  cav = Lattice.param.Cavity_on; path = Lattice.param.pathlength;

  circ = Lattice.Cell[Lattice.param.Cell_nLoc].S;

  // damped system
  Lattice.param.radiation = true; Lattice.param.emittance  = true;
  Lattice.param.Cavity_on = true; Lattice.param.pathlength = false;

  Lattice.Ring_GetTwiss(false, 0.0);

  // radiation loss is computed in Cav_Pass

  Lattice.param.U0 = Lattice.param.dE*1e9*Lattice.param.Energy;
  C = static_cast<CavityType*>(&Lattice.Cell[Lattice.Elem_GetPos(Fnum, 1)]);
  V_RF = C->volt;
  h_RF = C->h;
  phi0 = fabs(asin(Lattice.param.U0/V_RF));
  Lattice.param.delta_RF =
    sqrt(-V_RF*cos(M_PI-phi0)*(2.0-(M_PI-2.0*(M_PI-phi0))
    *tan(M_PI-phi0))/(fabs(Lattice.param.Alphac)*M_PI*h_RF
		      *1e9*Lattice.param.Energy));

  // Compute diffusion coeffs. for eigenvectors [sigma_xx, sigma_yy, sigma_zz]
  putlinmat(6, Lattice.param.Ascr, Ascr_map); Ascr_map += Lattice.param.CODvect;

  Lattice.Cell_Pass(0, Lattice.param.Cell_nLoc, Ascr_map, lastpos);

  // K. Robinson "Radiation Effects in Circular Electron Accelerators"
  // Phys. Rev. 111 (2), 373-380.
  // Iu.F. Orlov, E.K. Tarasov "Damping of Oscillations in a Cyclic Electron
  // Accelerator" J. Exptl. Theoet. Phys. 34, 651-657 (1958).
  // To leading order:
  //   Sum_k(alpha_k) = 2*U_0/E
  // or
  //   J_x + J_y + J_z = 4.

  for (i = 0; i < DOF; i++) {
    // partition numbers
    Lattice.param.J[i] =
      2.0*(1.0+Lattice.param.CODvect[delta_])*Lattice.param.alpha_rad[i]
      /Lattice.param.dE;
    // damping times
    Lattice.param.tau[i] = -circ/(c0*Lattice.param.alpha_rad[i]);
    // diffusion coeff. and emittance (alpha is for betatron amplitudes)
    Lattice.param.eps[i] =
      -Lattice.param.D_rad[i]/(2.0*Lattice.param.alpha_rad[i]);
    // fractional tunes
    nu[i]  = atan2(Lattice.param.wi[i*2], Lattice.param.wr[i*2])/(2.0*M_PI);
    if (nu[i] < 0.0) nu[i] = 1.0 + nu[i];
  }

  // undamped system
  Lattice.param.radiation = false; Lattice.param.emittance = false;

  Lattice.Ring_GetTwiss(false, 0.0);

  // Compute sigmas arround the lattice:
  //   Sigma = A diag[J_1, J_1, J_2, J_2, J_3, J_3] A^T
  for (i = 0; i < 6; i++) {
    Ascr_map[i] = tps(Lattice.param.CODvect[i]);
    for (j = 0; j < 6; j++)
      Ascr_map[i] +=
	Lattice.param.Ascr[i][j]*sqrt(Lattice.param.eps[j/2])*tps(0.0, j+1);
  }
  for (loc = 0; loc <= Lattice.param.Cell_nLoc; loc++) {
    Elem_Pass(loc, Ascr_map);
    // sigma = A x A^tp
    getlinmat(6, Ascr_map, Lattice.Cell[loc].sigma);
    TpMat(6, Lattice.Cell[loc].sigma);
    getlinmat(6, Ascr_map, Ascr); MulLMat(6, Ascr, Lattice.Cell[loc].sigma);
  }

  // A. W. Chao, M. J. Lee "Particle Distribution Parameters in an Electron
  // Storage Ring" J. Appl. Phys. 47 (10), 4453-4456 (1976).
  // observable tilt angle
  theta = atan2(2e0*Lattice.Cell[0].sigma[x_][y_],
	  (Lattice.Cell[0].sigma[x_][x_]-Lattice.Cell[0].sigma[y_][y_]))/2e0;

  // longitudinal alpha and beta
  Lattice.param.alpha_z =
    -Lattice.param.Ascr[ct_][ct_]*Lattice.param.Ascr[delta_][ct_]
    - Lattice.param.Ascr[ct_][delta_]*Lattice.param.Ascr[delta_][delta_];
  Lattice.param.beta_z =
    sqr(Lattice.param.Ascr[ct_][ct_]) + sqr(Lattice.param.Ascr[ct_][delta_]);
  gamma_z = (1.0+sqr(Lattice.param.alpha_z))/Lattice.param.beta_z;

  // bunch size
  sigma_s = sqrt(Lattice.param.beta_z*Lattice.param.eps[Z_]);
  sigma_delta = sqrt(gamma_z*Lattice.param.eps[Z_]);

  if (prt) {
    printf("\nEmittance:\n");
    printf("\nBeam energy [GeV]:              "
	   "Eb          = %4.2f\n", Lattice.param.Energy);
    printf("Energy loss per turn [keV]:     "
	   "U0          = %3.1f\n",
	   1e-3*Lattice.param.U0);
    printf("Synchronous phase [deg]:        "
	   "phi0        = 180 - %4.2f\n",
	   phi0*180.0/M_PI);
    printf("RF bucket height [%%]:           "
	   "delta_RF    = %4.2f\n", 1e2*Lattice.param.delta_RF);
    printf("\n");
    printf("Equilibrium emittance [m.rad]:  "
	   "eps_x       = %9.3e, eps_y  = %9.3e, eps_z  = %9.3e\n",
	   Lattice.param.eps[X_], Lattice.param.eps[Y_],
	   Lattice.param.eps[Z_]);
    printf("Bunch length [mm]:              "
	   "sigma_s     = %5.3f\n", 1e3*sigma_s);
    printf("Momentum spread:                "
	   "sigma_delta = %9.3e\n", sigma_delta);
    printf("Partition numbers:              "
	   "J_x         = %5.3f,     J_y    = %5.3f,     J_z    = %5.3f\n",
            Lattice.param.J[X_], Lattice.param.J[Y_], Lattice.param.J[Z_]);
    printf("Damping times [msec]:           "
	   "tau_x       = %3.1f,      tau_y  = %3.1f,      tau_z  = %3.1f\n",
	   1e3*Lattice.param.tau[X_], 1e3*Lattice.param.tau[Y_],
	   1e3*Lattice.param.tau[Z_]);
    printf("\n");
    printf("alphac:                         "
	   "alphac      = %8.4e\n", Lattice.param.Alphac);
    printf("\n");
    printf("Fractional tunes:               "
	   "nu_x        = %7.5f, nu_y   = %7.5f, nu_z   = %7.5f\n",
	   nu[X_], nu[Y_], nu[Z_]);
    printf("                                "
	   "1-nu_x      = %7.5f, 1-nu_y = %7.5f, 1-nu_z = %7.5f\n",
	   1e0-nu[X_], 1e0-nu[Y_], 1e0-nu[Z_]);
    printf("\n");
    printf("sigmas:                         "
	   "sigma_x     = %5.1f microns, sigma_px    = %5.1f urad\n",
	   1e6*sqrt(Lattice.Cell[0].sigma[x_][x_]),
	   1e6*sqrt(Lattice.Cell[0].sigma[px_][px_]));
    printf("                                "
	   "sigma_y     = %5.1f microns, sigma_py    = %5.1f urad\n",
	   1e6*sqrt(Lattice.Cell[0].sigma[y_][y_]),
	   1e6*sqrt(Lattice.Cell[0].sigma[py_][py_]));
    printf("                                "
	   "sigma_s     = %5.2f mm,      sigma_delta = %8.2e\n",
	   1e3*sqrt(Lattice.Cell[0].sigma[ct_][ct_]),
	   sqrt(Lattice.Cell[0].sigma[delta_][delta_]));

    printf("\n");
    printf("Beam ellipse twist [rad]:       tw = %5.3f\n", theta);
    printf("                   [deg]:       tw = %5.3f\n", theta*180.0/M_PI);
  }

  // restore state
  Lattice.param.radiation = rad; Lattice.param.emittance  = emit;
  Lattice.param.Cavity_on = cav; Lattice.param.pathlength = path;
}


// output

double get_code(CellType &Cell)
{
  double    code;
  MpoleType *M;

  switch (Cell.Kind) {
  case drift:
    code = 0.0;
    break;
  case Mpole:
    M = static_cast<MpoleType*>(&Cell);
    if (M->irho != 0.0)
      code = sgn(M->irho)*0.5;
    else if (M->Bpar[Quad+HOMmax] != 0)
      code = sgn(M->Bpar[Quad+HOMmax]);
    else if (M->Bpar[Sext+HOMmax] != 0)
      code = 1.5*sgn(M->Bpar[Sext+HOMmax]);
    else if (M->Bpar[Oct+HOMmax] != 0)
      code = 1.75*sgn(M->Bpar[Oct+HOMmax]);
    else if (M->Bpar[Dec+HOMmax] != 0)
      code = 1.75*sgn(M->Bpar[Dec+HOMmax]);
    else if (M->Bpar[Dodec+HOMmax] != 0)
      code = 1.75*sgn(M->Bpar[Dodec+HOMmax]);
    else if (Cell.Fnum == Lattice.param.bpm)
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
			   const int Fnum, const bool all)
{
  long int i = 0;
  double   I5;
  FILE     *outf;

  outf = file_write(fname);
  fprintf(outf, "#        name             s     code"
	        "   alphax   betax     nux      etax    etapx");
  fprintf(outf, "     alphay   betay     nuy      etay    etapy      I5\n");
  fprintf(outf, "#                        [m]"
	        "                     [m]                [m]");
  fprintf(outf, "                        [m]                [m]\n");
  fprintf(outf, "#\n");

  I2 = 0.0; I5 = 0.0;
  for (i = loc1; i <= loc2; i++) {
    if (all || (Lattice.Cell[i].Fnum == Fnum)) {
      fprintf(outf, "%4ld %15s %9.5f %4.1f"
	      " %9.5f %8.5f %8.5f %8.5f %8.5f"
	      " %9.5f %8.5f %8.5f %8.5f %8.5f  %8.2e\n",
	      i, Lattice.Cell[i].Name, Lattice.Cell[i].S,
	      get_code(Lattice.Cell[i]),
	      Lattice.Cell[i].Alpha[X_], Lattice.Cell[i].Beta[X_],
	      Lattice.Cell[i].Nu[X_],
	      Lattice.Cell[i].Eta[X_], Lattice.Cell[i].Etap[X_],
	      Lattice.Cell[i].Alpha[Y_], Lattice.Cell[i].Beta[Y_],
	      Lattice.Cell[i].Nu[Y_],
	      Lattice.Cell[i].Eta[Y_], Lattice.Cell[i].Etap[Y_], I5);
    }
  }

//  fprintf(outf, "\n");
//  fprintf(outf, "# emittance: %5.3f nm.rad\n", get_eps_x());

  fclose(outf);
}


void LatticeType::prt_lat(const char *fname, const int Fnum, const bool all)
{
  prt_lat(0, Lattice.param.Cell_nLoc, fname, Fnum, all);
}


void LatticeType::Cell_Twiss(const long int i0, const long int i1) {
  long int     i;
  int          k, nu_int[2];
  double       alpha[2], beta[2], dnu[2], eta[2], etap[2];
  ss_vect<tps> A;

  for (k = 0; k < 2; k++)
    nu_int[k] = 0;

  for (i = i0; i <= i1; i++) {
    putlinmat(6, Lattice.Cell[i].A, A);
    get_ab(A, alpha, beta, dnu, eta, etap);

    for (k = 0; k < 2; k++) {
      Lattice.Cell[i].Alpha[k] = alpha[k]; Lattice.Cell[i].Beta[k] = beta[k];
      Lattice.Cell[i].Nu[k] = nu_int[k] + dnu[k];

      if (i > i0) {
	if((Lattice.Cell[i].Nu[k] < Lattice.Cell[i-1].Nu[k])
	   && (Lattice.Cell[i].L >= 0e0)) {
	  Lattice.Cell[i].Nu[k] += 1e0; nu_int[k] += 1;
	} else if((Lattice.Cell[i].Nu[k] > Lattice.Cell[i-1].Nu[k]) &&
		  (Lattice.Cell[i].L < 0e0))
	  nu_int[k] -= 1;
      }

      Lattice.Cell[i].Eta[k] = eta[k]; Lattice.Cell[i].Etap[k] = etap[k];
    }
  }
}


void LatticeType::prt_lat(const int loc1, const int loc2, const char *fname,
			   const int Fnum, const bool all, const int n)
{
  long int        i = 0;
  int             j, k;
  double          s, h;
  double          alpha[2], beta[2], nu[2], dnu[2], eta[2], etap[2], dnu1[2];
  double          curly_H;
  MpoleType       *M;
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
    if (all || (Lattice.Cell[i].Fnum == Fnum)) {
      if ((i != 0) &&
	  ((Lattice.Cell[i].Kind == drift) ||
	   ((Lattice.Cell[i].Kind == Mpole)
	    && (Lattice.Cell[i].L != 0e0)))) {
	M = static_cast<MpoleType*>(&Lattice.Cell[i]);

	for (k = 0; k < 2; k++) {
	  alpha[k] = Lattice.Cell[i-1].Alpha[k];
	  beta[k] = Lattice.Cell[i-1].Beta[k];
	  nu[k] = Lattice.Cell[i-1].Nu[k];
	  eta[k] = Lattice.Cell[i-1].Eta[k];
	  etap[k] = Lattice.Cell[i-1].Etap[k];
	}

	A = get_A(alpha, beta, eta, etap);

	s = Lattice.Cell[i].S - Lattice.Cell[i].L;
	h = Lattice.Cell[i].L/n;

	for (j = 1; j <= n; j++) {
	  s += h;

	  if (Lattice.Cell[i].Kind == drift)
	    Drift(h, A);
	  else if (Lattice.Cell[i].Kind == Mpole) {
	    if ((j == 1) && (M->irho != 0e0))
	      EdgeFocus(M->irho, M->Tx1, M->gap, A);

	    Drift(c1*h, A);
	    thin_kick(Quad, M->Bpar, d1*h, M->irho, M->irho, A);
	    Drift(c2*h, A);
	    thin_kick(Quad, M->Bpar, d2*h, M->irho, M->irho, A);
	    Drift(c2*h, A);
	    thin_kick(Quad, M->Bpar, d1*h, M->irho, M->irho, A);
	    Drift(c1*h, A);

	    if ((j == n) && (M->irho != 0e0))
	      EdgeFocus(M->irho, M->Tx2, M->gap, A);
	  }

	  get_ab(A, alpha, beta, dnu, eta, etap);

	  if(Lattice.Cell[i].L < 0e0)
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
		  i, Lattice.Cell[i].Name, s, get_code(Lattice.Cell[i]),
		  alpha[X_], beta[X_], nu[X_]+dnu[X_], eta[X_], etap[X_],
		  alpha[Y_], beta[Y_], nu[Y_]+dnu[Y_], eta[Y_], etap[Y_],
		  eta_Fl[x_], eta_Fl[px_], curly_H);
	}
      } else {
	A = get_A(Lattice.Cell[i].Alpha, Lattice.Cell[i].Beta,
		  Lattice.Cell[i].Eta, Lattice.Cell[i].Etap);

	eta_Fl.zero();
	for (k = 0; k < 2; k++) {
	  eta_Fl[2*k] = Lattice.Cell[i].Eta[k];
	  eta_Fl[2*k+1] = Lattice.Cell[i].Etap[k];
	}
	eta_Fl = (Inv(A)*eta_Fl).cst();
	curly_H = sqr(eta_Fl[x_]) + sqr(eta_Fl[px_]);

	fprintf(outf, "%4ld %15s %6.2f %4.1f"
		" %9.5f %8.5f %8.5f %11.8f %11.8f"
		" %9.5f %8.5f %8.5f %8.5f %8.5f %10.3e %10.3e %10.3e\n",
		i, Lattice.Cell[i].Name, Lattice.Cell[i].S,
		get_code(Lattice.Cell[i]),
		Lattice.Cell[i].Alpha[X_], Lattice.Cell[i].Beta[X_],
		Lattice.Cell[i].Nu[X_],
		Lattice.Cell[i].Eta[X_], Lattice.Cell[i].Etap[X_],
		Lattice.Cell[i].Alpha[Y_], Lattice.Cell[i].Beta[Y_],
		Lattice.Cell[i].Nu[Y_],
		Lattice.Cell[i].Eta[Y_], Lattice.Cell[i].Etap[Y_],
		eta_Fl[x_], eta_Fl[px_], curly_H);
      }
    }
  }

  fclose(outf);
}


void LatticeType::prt_lat(const char *fname, const int Fnum, const bool all,
			   const int n)
{
  prt_lat(0, Lattice.param.Cell_nLoc, fname, Fnum, all, n);
}


void LatticeType::prt_chrom_lat(void)
{
  long int  i;
  double    dbeta_ddelta[Cell_nLocMax][2], detax_ddelta[Cell_nLocMax];
  double    ksi[Cell_nLocMax][2];
  MpoleType *M;
  FILE      *outf;

  printf("\n");
  printf("prt_chrom_lat: calling Ring_GetTwiss with delta != 0\n");
  Lattice.Ring_GetTwiss(true, Lattice.param.dPcommon);
  for (i = 0; i <= Lattice.param.Cell_nLoc; i++) {
    dbeta_ddelta[i][X_] = Lattice.Cell[i].Beta[X_];
    dbeta_ddelta[i][Y_] = Lattice.Cell[i].Beta[Y_];
    detax_ddelta[i] = Lattice.Cell[i].Eta[X_];
  }
  printf("prt_chrom_lat: calling Ring_GetTwiss with delta != 0\n");
  Lattice.Ring_GetTwiss(true, -Lattice.param.dPcommon);
  ksi[0][X_] = 0.0; ksi[0][Y_] = 0.0;
  for (i = 0; i <= Lattice.param.Cell_nLoc; i++) {
    dbeta_ddelta[i][X_] -= Lattice.Cell[i].Beta[X_];
    dbeta_ddelta[i][Y_] -= Lattice.Cell[i].Beta[Y_];
    detax_ddelta[i] -= Lattice.Cell[i].Eta[X_];
    dbeta_ddelta[i][X_] /= 2.0*Lattice.param.dPcommon;
    dbeta_ddelta[i][Y_] /= 2.0*Lattice.param.dPcommon;
    detax_ddelta[i] /= 2.0*Lattice.param.dPcommon;
    if (i != 0) {
      ksi[i][X_] = ksi[i-1][X_]; ksi[i][Y_] = ksi[i-1][Y_];
    }
    if (Lattice.Cell[i].Kind == Mpole) {
      M = static_cast<MpoleType*>(&(Lattice.Cell[i]);
      ksi[i][X_] -=
	M->Bpar[Quad+HOMmax]
	*Lattice.Cell[i].L*Lattice.Cell[i].Beta[X_]/(4.0*M_PI);
      ksi[i][Y_] +=
	M->Bpar[Quad+HOMmax]
	*Lattice.Cell[i].L*Lattice.Cell[i].Beta[Y_]/(4.0*M_PI);
    }
  }

  outf = file_write("chromlat.out");
  fprintf(outf, "#     name              s    code"
	        "  bx*ex  sqrt(bx*by)  dbx/dd*ex  bx*dex/dd"
	        "  by*ex  dby/dd*ex by*dex/dd  ksix  ksiy"
	        "  dbx/dd  bx/dd dex/dd\n");
  fprintf(outf, "#                      [m]          [m]"
	        "      [m]          [m]       [m]");
  fprintf(outf, "       [m]      [m]       [m]\n");
  fprintf(outf, "#\n");
  for (i = 0; i <= Lattice.param.Cell_nLoc; i++) {
    fprintf(outf, "%4ld %15s %6.2f %4.1f"
	          "  %6.3f  %8.3f    %8.3f   %8.3f"
	          "   %6.3f %8.3f   %8.3f  %5.2f  %5.2f"
	          "  %6.3f  %6.3f  %6.3f\n",
	    i, Lattice.Cell[i].Name, Lattice.Cell[i].S,
	    get_code(Lattice.Cell[i]),
	    Lattice.Cell[i].Beta[X_]*Lattice.Cell[i].Eta[X_],
	    sqrt(Lattice.Cell[i].Beta[X_]*Lattice.Cell[i].Beta[Y_]),
	    dbeta_ddelta[i][X_]*Lattice.Cell[i].Eta[X_],
	    detax_ddelta[i]*Lattice.Cell[i].Beta[X_],
	    Lattice.Cell[i].Beta[Y_]*Lattice.Cell[i].Eta[X_],
	    dbeta_ddelta[i][Y_]*Lattice.Cell[i].Eta[X_],
	    detax_ddelta[i]*Lattice.Cell[i].Beta[Y_],
	    ksi[i][X_], ksi[i][Y_],
	    dbeta_ddelta[i][X_], dbeta_ddelta[i][Y_], detax_ddelta[i]);
  }
  fclose(outf);
}


void LatticeType::prt_cod(const char *file_name, const int Fnum,
			   const bool all)
{
  long   i;
  FILE   *outf;
  long   FORLIM;
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

  FORLIM = Lattice.param.Cell_nLoc;
  for (i = 0L; i <= FORLIM; i++) {
    if (all || (Lattice.Cell[i].Fnum == Fnum)) {
      /* COD is in local coordinates */
      fprintf(outf, "%4ld %.*s %6.2f %4.1f %6.3f %6.3f %6.3f %6.3f"
	      " %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",
	      i, SymbolLength, Lattice.Cell[i].Name, Lattice.Cell[i].S,
	      get_code(Lattice.Cell[i]),
	      Lattice.Cell[i].Beta[X_], Lattice.Cell[i].Nu[X_],
	      Lattice.Cell[i].Beta[Y_], Lattice.Cell[i].Nu[Y_],
	      1e3*Lattice.Cell[i].BeamPos[x_], 1e3*Lattice.Cell[i].BeamPos[y_],
	      1e3*Lattice.Cell[i].dS[X_], 1e3*Lattice.Cell[i].dS[Y_],
	      -1e3*Lattice.Elem_GetKval(Lattice.Cell[i].Fnum,
					Lattice.Cell[i].Knum,Dip),
	      1e3*Lattice.Elem_GetKval(Lattice.Cell[i].Fnum,
				       Lattice.Cell[i].Knum, -Dip));
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

  for (k = 0; k <= Lattice.param.Cell_nLoc; k++)
    fprintf(outf, "%4ld %.*s %6.2f %4.1f %12.5e %12.5e\n",
	    k, SymbolLength, Lattice.Cell[k].Name, Lattice.Cell[k].S,
	    get_code(Lattice.Cell[k]),
	    Lattice.Cell[k].BeamPos[x_], Lattice.Cell[k].BeamPos[y_]);

  fclose(outf);
}


// misalignments

void LatticeType::CheckAlignTol(const char *OutputFile)
  // check aligment errors of individual magnets on giders
  // the dT and roll angle are all printed out
{
  int  i, j;
  int  n_girders;
  int  gs_Fnum, ge_Fnum;
  int  gs_nKid, ge_nKid;
  int  dip_Fnum,dip_nKid;
  int  loc, loc_gs, loc_ge;
  char * name;
  double s;
  double dSsys[2], dSrms[2], dSrnd[2], dS[2], dT[2];
  std::fstream fout;

  gs_Fnum = Lattice.param.gs;   gs_nKid = Lattice.GetnKid(gs_Fnum);
  ge_Fnum = Lattice.param.ge;   ge_nKid = Lattice.GetnKid(ge_Fnum);
  if (gs_nKid == ge_nKid)
    n_girders= gs_nKid;
  else {
    std::cout << " The numbers of GS and GE not same. " << std::endl;
    exit (1);
  }

  fout.open(OutputFile,std::ios::out);
  if(!fout) {
    std::cout << "error in opening the file  " << std::endl;
    exit_(0);
  }

  fout << "Girders, Quads, Sexts:  " << std::endl;
  for (i = 1; i <= n_girders; i++){
    fout << i << ":" << std::endl;
    loc_gs = Lattice.Elem_GetPos(gs_Fnum, i);
    loc_ge = Lattice.Elem_GetPos(ge_Fnum, i);

    loc = loc_gs;
    dSsys[X_] = Lattice.Cell[loc].Elem.M->dSsys[X_];
    dSsys[Y_] = Lattice.Cell[loc].Elem.M->dSsys[Y_];
    dSrms[X_] = Lattice.Cell[loc].Elem.M->dSrms[X_];
    dSrms[Y_] = Lattice.Cell[loc].Elem.M->dSrms[Y_];
    dSrnd[X_] = Lattice.Cell[loc].Elem.M->dSrnd[X_];
    dSrnd[Y_] = Lattice.Cell[loc].Elem.M->dSrnd[Y_];
    dS[X_] = Lattice.Cell[loc].dS[X_]; dS[Y_] = Lattice.Cell[loc].dS[Y_];
    dT[0] = Lattice.Cell[loc].dT[0]; dT[1] = Lattice.Cell[loc].dT[1];
    s = Lattice.Cell[loc].S; name = Lattice.Cell[loc].Name;
    fout << "  " << name << "  " << loc << "   " << s
	 << "  " <<  dSsys[X_] << "  " <<  dSsys[Y_]
	 << "   " << dSrms[X_] << "  " <<  dSrms[Y_]
	 << "   " << dSrnd[X_] << "  " <<  dSrnd[Y_]
         << "   " << Lattice.Cell[loc].Elem.M->dTrms << "  "
	 << Lattice.Cell[loc].Elem.M->dTrnd << "   " << dS[X_]     << "  "
	 <<  dS[Y_] << "   " << atan2( dT[1], dT[0] )  << std::endl;

    for (j = loc_gs+1; j < loc_ge; j++) {
      if ((Lattice.Cell[j].Kind == Mpole) &&
	  (Lattice.Cell[j].Elem.M->n_design >= Quad ||
	   Lattice.Cell[j].Elem.M->n_design >= Sext)) {
        loc = j;
	dSsys[X_] = Lattice.Cell[loc].Elem.M->dSsys[X_];
	dSsys[Y_] = Lattice.Cell[loc].Elem.M->dSsys[Y_];
	dSrms[X_] = Lattice.Cell[loc].Elem.M->dSrms[X_];
	dSrms[Y_] = Lattice.Cell[loc].Elem.M->dSrms[Y_];
	dSrnd[X_] = Lattice.Cell[loc].Elem.M->dSrnd[X_];
	dSrnd[Y_] = Lattice.Cell[loc].Elem.M->dSrnd[Y_];
	dS[X_] = Lattice.Cell[loc].dS[X_]; dS[Y_] = Lattice.Cell[loc].dS[Y_];
	dT[0] = Lattice.Cell[loc].dT[0];   dT[1] = Lattice.Cell[loc].dT[1];
	s = Lattice.Cell[loc].S; name=Lattice.Cell[loc].Name;
	fout << "  " << name << "  " << loc << "   " << s
	     << "  " <<  dSsys[X_] << "  " <<  dSsys[Y_]
	     << "   " << dSrms[X_] << "  " <<  dSrms[Y_]
	     << "   " << dSrnd[X_] << "  " <<  dSrnd[Y_]
	     << "   " << Lattice.Cell[loc].Elem.M->dTrms << "  "
	     << Lattice.Cell[loc].Elem.M->dTrnd
	     << "   " << dS[X_] << "  " <<  dS[Y_]
	     << "   " << atan2( dT[1], dT[0] )  << std::endl;
      }
    }

    loc = loc_ge;
    dSsys[X_] = Lattice.Cell[loc].Elem.M->dSsys[X_];
    dSsys[Y_] = Lattice.Cell[loc].Elem.M->dSsys[Y_];
    dSrms[X_] = Lattice.Cell[loc].Elem.M->dSrms[X_];
    dSrms[Y_] = Lattice.Cell[loc].Elem.M->dSrms[Y_];
    dSrnd[X_] = Lattice.Cell[loc].Elem.M->dSrnd[X_];
    dSrnd[Y_] = Lattice.Cell[loc].Elem.M->dSrnd[Y_];
    dS[X_] = Lattice.Cell[loc].dS[X_]; dS[Y_] = Lattice.Cell[loc].dS[Y_];
    dT[0] = Lattice.Cell[loc].dT[0]; dT[1] = Lattice.Cell[loc].dT[1];
    s=Lattice.Cell[loc].S; name=Lattice.Cell[loc].Name;
    fout << "  " << name << "  " << loc << "   " << s
	 << "  " <<  dSsys[X_] << "  " <<  dSsys[Y_]
	 << "   " << dSrms[X_] << "  " <<  dSrms[Y_]
	 << "   " << dSrnd[X_] << "  " <<  dSrnd[Y_]
         << "   " << Lattice.Cell[loc].Elem.M->dTrms
	 << "  " << Lattice.Cell[loc].Elem.M->dTrnd
         << "   " << dS[X_]     << "  " <<  dS[Y_]
	 << "   " << atan2( dT[1], dT[0] )  << std::endl;

  }

  fout << "  " << std::endl;
  fout << "Dipoles:  " << std::endl;
  dip_Fnum = Lattice.Elem_Index("B1"); dip_nKid = Lattice.GetnKid(dip_Fnum);
  for (i = 1; i <= dip_nKid; i++){
    loc = Lattice.Elem_GetPos(dip_Fnum, i);
    dSsys[X_] = Lattice.Cell[loc].Elem.M->dSsys[X_];
    dSsys[Y_] = Lattice.Cell[loc].Elem.M->dSsys[Y_];
    dSrms[X_] = Lattice.Cell[loc].Elem.M->dSrms[X_];
    dSrms[Y_] = Lattice.Cell[loc].Elem.M->dSrms[Y_];
    dSrnd[X_] = Lattice.Cell[loc].Elem.M->dSrnd[X_];
    dSrnd[Y_] = Lattice.Cell[loc].Elem.M->dSrnd[Y_];
    dS[X_] = Lattice.Cell[loc].dS[X_]; dS[Y_] = Lattice.Cell[loc].dS[Y_];
    dT[0] = Lattice.Cell[loc].dT[0]; dT[1] = Lattice.Cell[loc].dT[1];
    s = Lattice.Cell[loc].S; name = Lattice.Cell[loc].Name;
    fout << "  " << name << "  " << loc << "   " << s
	 << "  " <<  dSsys[X_] << "  " <<  dSsys[Y_]
	 << "   " << dSrms[X_] << "  " <<  dSrms[Y_]
	 << "   " << dSrnd[X_] << "  " <<  dSrnd[Y_]
	 << "   " << Lattice.Cell[loc].Elem.M->dTrms
	 << "  " << Lattice.Cell[loc].Elem.M->dTrnd
	 << "   " << dS[X_]     << "  " <<  dS[Y_]
	 << "   " << atan2( dT[1], dT[0] )  << std::endl;
  }

  fout.close();
}


void misalign_rms_elem(const int Fnum, const int Knum,
		       const double dx_rms, const double dy_rms,
		       const double dr_rms, const bool new_rnd)
{
  long int   loc;
  MpoleType  *mp;

  loc = Lattice.Elem_GetPos(Fnum, Knum); mp = Lattice.Cell[loc].Elem.M;

  mp->dSrms[X_] = dx_rms; mp->dSrms[Y_] = dy_rms; mp->dTrms = dr_rms;
  if (new_rnd) {
    if (normal) {
      mp->dSrnd[X_] = normranf(); mp->dSrnd[Y_] = normranf();
      mp->dTrnd = normranf();
    } else {
      mp->dSrnd[X_] = ranf(); mp->dSrnd[Y_] = ranf();
      mp->dTrnd = ranf();
    }
  }

  Mpole_SetdS(Fnum, Knum); Mpole_SetdT(Fnum, Knum);
}

void misalign_sys_elem(const int Fnum, const int Knum,
		       const double dx_sys, const double dy_sys,
		       const double dr_sys)
{
  long int   loc;
  MpoleType  *mp;

  loc = Lattice.Elem_GetPos(Fnum, Knum); mp = Lattice.Cell[loc].Elem.M;

  mp->dSsys[X_] = dx_sys; mp->dSsys[Y_] = dy_sys; mp->dTsys = dr_sys;

  Mpole_SetdS(Fnum, Knum); Mpole_SetdT(Fnum, Knum);
}

void misalign_rms_fam(const int Fnum,
		      const double dx_rms, const double dy_rms,
		      const double dr_rms, const bool new_rnd)
{
  int  i;

  for (i = 1; i <= Lattice.GetnKid(Fnum); i++)
    misalign_rms_elem(Fnum, i, dx_rms, dy_rms, dr_rms, new_rnd);
}

void misalign_sys_fam(const int Fnum,
		      const double dx_sys, const double dy_sys,
		      const double dr_sys)
{
  int  i;

  for (i = 1; i <= Lattice.GetnKid(Fnum); i++)
    misalign_sys_elem(Fnum, i, dx_sys, dy_sys, dr_sys);
}

void misalign_rms_type(const int type,
		       const double dx_rms, const double dy_rms,
		       const double dr_rms, const bool new_rnd)
{
  long int   k;

  if ((type >= All) && (type <= HOMmax)) {
    for (k = 1; k <= Lattice.param.Cell_nLoc; k++) {
      if ((Lattice.Cell[k].Kind == Mpole) &&
	  ((type == Lattice.Cell[k].Elem.M->n_design) ||
	  ((type == All) &&
	   ((Lattice.Cell[k].Fnum != Lattice.param.gs) &&
	    (Lattice.Cell[k].Fnum != Lattice.param.ge))))) {
	// if all: skip girders
	misalign_rms_elem(Lattice.Cell[k].Fnum, Lattice.Cell[k].Knum,
			  dx_rms, dy_rms, dr_rms, new_rnd);
      }
    }
  } else {
    printf("misalign_rms_type: incorrect type %d\n", type); exit_(1);
  }
}

void misalign_sys_type(const int type,
		       const double dx_sys, const double dy_sys,
		       const double dr_sys)
{
  long int   k;

  if ((type >= All) && (type <= HOMmax)) {
    for (k = 1; k <= Lattice.param.Cell_nLoc; k++) {
      if ((Lattice.Cell[k].Kind == Mpole) &&
	  ((type == Lattice.Cell[k].Elem.M->n_design) ||
	  ((type == All) &&
	   ((Lattice.Cell[k].Fnum != Lattice.param.gs) &&
	    (Lattice.Cell[k].Fnum != Lattice.param.ge))))) {
	// if all: skip girders
	misalign_sys_elem(Lattice.Cell[k].Fnum, Lattice.Cell[k].Knum,
			  dx_sys, dy_sys, dr_sys);
      }
    }
  } else {
    printf("misalign_sys_type: incorrect type %d\n", type); exit_(1);
  }
}

void misalign_rms_girders(const int gs, const int ge,
			  const double dx_rms, const double dy_rms,
			  const double dr_rms, const bool new_rnd)
{
  int       i, k, n_girders, n_ge, n_gs;
  long int  loc_gs, loc_ge, j;
  double    s_gs, s_ge, dx_gs[2], dx_ge[2], s;

  n_gs = Lattice.GetnKid(gs); n_ge = Lattice.GetnKid(ge);

  if (n_gs == n_ge)
    n_girders = n_gs;
  else {
    std::cout << "set_girders: no of GS != no of GE" << std::endl;
    exit (1);
  }

  misalign_rms_fam(gs, dx_rms, dy_rms, dr_rms, new_rnd);
  misalign_rms_fam(ge, dx_rms, dy_rms, dr_rms, new_rnd);

  for (i = 1; i <= n_girders; i++) {
    loc_gs = Lattice.Elem_GetPos(gs, i); loc_ge = Lattice.Elem_GetPos(ge, i);
    s_gs = Lattice.Cell[loc_gs].S; s_ge = Lattice.Cell[loc_ge].S;

    // roll for a rigid boby
    // Note, girders needs to be introduced as gs->ge pairs
    Lattice.Cell[loc_ge].Elem.M->dTrnd = Lattice.Cell[loc_gs].Elem.M->dTrnd;
    Mpole_SetdT(ge, i);

    for (k = 0; k <= 1; k++) {
      dx_gs[k] = Lattice.Cell[loc_gs].dS[k];
      dx_ge[k] = Lattice.Cell[loc_ge].dS[k];
    }

    // move elements onto mis-aligned girder
    for (j = loc_gs+1; j < loc_ge; j++) {
      if ((Lattice.Cell[j].Kind == Mpole) ||
	  (Lattice.Cell[j].Fnum == Lattice.param.bpm)) {
        s = Lattice.Cell[j].S;
	for (k = 0; k <= 1; k++)
	  Lattice.Cell[j].Elem.M->dSsys[k]
	    = dx_gs[k] + (dx_ge[k]-dx_gs[k])*(s-s_gs)/(s_ge-s_gs);
	Lattice.Cell[j].Elem.M->dTsys =
	  Lattice.Cell[loc_gs].Elem.M->dTrms*Lattice.Cell[loc_gs].Elem.M->dTrnd;
      }
    }
  }
}


void misalign_sys_girders(const int gs, const int ge,
			  const double dx_sys, const double dy_sys,
			  const double dr_sys)
{
  int       i, k, n_girders, n_ge, n_gs;
  long int  loc_gs, loc_ge, j;
  double    s_gs, s_ge, dx_gs[2], dx_ge[2], s;

  n_gs = Lattice.GetnKid(gs); n_ge = Lattice.GetnKid(ge);

  if (n_gs == n_ge)
    n_girders = n_gs;
  else {
    std::cout << "set_girders: no of GS != no of GE" << std::endl;
    exit (1);
  }

  misalign_sys_fam(gs, dx_sys, dy_sys, dr_sys);
  misalign_sys_fam(ge, dx_sys, dy_sys, dr_sys);

  for (i = 1; i <= n_girders; i++) {
    loc_gs = Lattice.Elem_GetPos(gs, i); loc_ge = Lattice.Elem_GetPos(ge, i);
    s_gs = Lattice.Cell[loc_gs].S; s_ge = Lattice.Cell[loc_ge].S;

    // roll for a rigid boby
    // Note, girders needs to be introduced as gs->ge pairs
    Lattice.Cell[loc_ge].Elem.M->dTrnd = Lattice.Cell[loc_gs].Elem.M->dTrnd;
    Mpole_SetdT(ge, i);

    for (k = 0; k <= 1; k++) {
      dx_gs[k] = Lattice.Cell[loc_gs].dS[k];
      dx_ge[k] = Lattice.Cell[loc_ge].dS[k];
    }

    // move elements onto mis-aligned girder
    for (j = loc_gs+1; j < loc_ge; j++) {
      if ((Lattice.Cell[j].Kind == Mpole)
	  || (Lattice.Cell[j].Fnum == Lattice.param.bpm)) {
        s = Lattice.Cell[j].S;
	for (k = 0; k <= 1; k++)
	  Lattice.Cell[j].Elem.M->dSsys[k]
	    = dx_gs[k] + (dx_ge[k]-dx_gs[k])*(s-s_gs)/(s_ge-s_gs);
	Lattice.Cell[j].Elem.M->dTsys =
	  Lattice.Cell[loc_gs].Elem.M->dTrms
	  *Lattice.Cell[loc_gs].Elem.M->dTrnd;
      }
    }
  }
}


// apertures

void set_aper_elem(const int Fnum, const int Knum,
		   const double Dxmin, const double Dxmax,
		   const double Dymin, const double Dymax)
{
  int  k;

    k = Lattice.Elem_GetPos(Fnum, Knum);
    Lattice.Cell[k].maxampl[X_][0] = Dxmin;
    Lattice.Cell[k].maxampl[X_][1] = Dxmax;
    Lattice.Cell[k].maxampl[Y_][0] = Dymin;
    Lattice.Cell[k].maxampl[Y_][1] = Dymax;
 }

void set_aper_fam(const int Fnum,
		  const double Dxmin, const double Dxmax,
		  const double Dymin, const double Dymax)
{
  int k;

  for (k = 1; k <= Lattice.GetnKid(Fnum); k++)
    set_aper_elem(Fnum, k, Dxmin, Dxmax, Dymin, Dymax);
}

void set_aper_type(const int type, const double Dxmin, const double Dxmax,
		   const double Dymin, const double Dymax)
{
  long int   k;

  if (type >= All && type <= HOMmax) {
    for(k = 1; k <= Lattice.param.Cell_nLoc; k++)
      if (((Lattice.Cell[k].Kind == Mpole) &&
	   (Lattice.Cell[k].Elem.M->n_design == type)) || (type == All))
	set_aper_elem(Lattice.Cell[k].Fnum, Lattice.Cell[k].Knum,
		      Dxmin, Dxmax, Dymin, Dymax);
  } else
    printf("set_aper_type: bad design type %d\n", type);
}


double get_L(const int Fnum, const int Knum)
{
  return Lattice.Cell[Lattice.Elem_GetPos(Fnum, Knum)].L;
}


void set_L(const int Fnum, const int Knum, const double L)
{
  long int  loc;
  double    phi;
  CellType  *cellp;
  elemtype  *elemp;
  MpoleType *M;

  loc = Lattice.Elem_GetPos(Fnum, Knum);
  cellp = &Lattice.Cell[loc]; elemp = &cellp->Elem;
  if (cellp->Kind == Mpole) {
    M = elemp->M;
    if (M->irho != 0e0) {
      // Phi is constant.
      phi = cellp->L*M->irho; M->irho = phi/L;
      // M->Pc0 = sin(phi/2e0);
    }
  }
  cellp->L = L;
}


void set_L(const int Fnum, const double L)
{
  int  k;

  for (k = 1; k <= Lattice.GetnKid(Fnum); k++)
    set_L(Fnum, k, L);
}


void set_dL(const int Fnum, const int Knum, const double dL)
{

  Lattice.Cell[Lattice.Elem_GetPos(Fnum, Knum)].L += dL;
}


void set_dL(const int Fnum, const double dL)
{
  int  k;

  for (k = 1; k <= Lattice.GetnKid(Fnum); k++)
    Lattice.Cell[Lattice.Elem_GetPos(Fnum, k)].L += dL;
}


// Multipole components.

void set_bn(const string type, const int Fnum, const int Knum,
	    const int n, const double bn, const double an)
{
  elemtype elem;

  if (n < 1) {
    std::cout << "set_bn: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  // elem = Lattice.Cell[Lattice.Elem_GetPos(Fnum, Knum)].Elem;
  // if (string == "design") {
  //   elem.M->Bpar[HOMmax+n] = bn; elem.M->Bpar[HOMmax-n] = an;
  // } else if (string == "rms") {

  // } else if (string == "sys") {
  // }
}

void get_bn_design_elem(const int Fnum, const int Knum,
			const int n, double &bn, double &an)
{
  elemtype  elem;

  if (n < 1) {
    std::cout << "get_bn_design_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  elem = Lattice.Cell[Lattice.Elem_GetPos(Fnum, Knum)].Elem;

  bn = elem.M->Bpar[HOMmax+n]; an = elem.M->Bpar[HOMmax-n];
}


void get_bnL_design_elem(const int Fnum, const int Knum,
			 const int n, double &bnL, double &anL)
{
  CellType *cellp;
  elemtype *elemp;

  if (n < 1) {
    std::cout << "get_bnL_design_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  cellp = &Lattice.Cell[Lattice.Elem_GetPos(Fnum, Knum)]; elemp = &cellp->Elem;

  bnL = elemp->M->Bpar[HOMmax+n]; anL = elemp->M->Bpar[HOMmax-n];

  if (cellp->L != 0.0) {
    bnL *= cellp->L; anL *= cellp->L;
  }
}


void set_bn_design_elem(const int Fnum, const int Knum,
			const int n, const double bn, const double an)
{
  elemtype  elem;

  if (n < 1) {
    std::cout << "set_bn_design_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  elem = Lattice.Cell[Lattice.Elem_GetPos(Fnum, Knum)].Elem;

  elem.M->Bpar[HOMmax+n] = bn; elem.M->Bpar[HOMmax-n] = an;

  Mpole_SetB(Fnum, Knum, n); Mpole_SetB(Fnum, Knum, -n);
}


void set_dbn_design_elem(const int Fnum, const int Knum,
			 const int n, const double dbn, const double dan)
{
  elemtype  elem;

  if (n < 1) {
    std::cout << "set_dbn_design_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  elem = Lattice.Cell[Lattice.Elem_GetPos(Fnum, Knum)].Elem;

  elem.M->Bpar[HOMmax+n] += dbn; elem.M->Bpar[HOMmax-n] += dan;

  Mpole_SetB(Fnum, Knum, n); Mpole_SetB(Fnum, Knum, -n);
}


void set_bn_design_fam(const int Fnum,
		       const int n, const double bn, const double an)
{
  int k;

  if (n < 1) {
    std::cout << "set_bn_design_fam: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  for (k = 1; k <= Lattice.GetnKid(Fnum); k++)
    set_bn_design_elem(Fnum, k, n, bn, an);
}


void set_dbn_design_fam(const int Fnum,
			const int n, const double dbn, const double dan)
{
  int k;

  if (n < 1) {
    std::cout << "set_dbn_design_fam: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  for (k = 1; k <= Lattice.GetnKid(Fnum); k++)
    set_dbn_design_elem(Fnum, k, n, dbn, dan);
}


void set_bnL_design_elem(const int Fnum, const int Knum,
			 const int n, const double bnL, const double anL)
{
  CellType *cellp;
  elemtype *elemp;

  if (n < 1) {
    std::cout << "set_bnL_design_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  cellp = &Lattice.Cell[Lattice.Elem_GetPos(Fnum, Knum)]; elemp = &cellp->Elem;

  if (cellp->L != 0.0) {
    elemp->M->Bpar[HOMmax+n] = bnL/cellp->L;
    elemp->M->Bpar[HOMmax-n] = anL/cellp->L;
  } else {
    // thin kick
    elemp->M->Bpar[HOMmax+n] = bnL; elemp->M->Bpar[HOMmax-n] = anL;
  }

  Mpole_SetB(Fnum, Knum, n); Mpole_SetB(Fnum, Knum, -n);
}


void set_dbnL_design_elem(const int Fnum, const int Knum,
			  const int n, const double dbnL, const double danL)
{
  CellType *cellp;
  elemtype *elemp;

  if (n < 1) {
    std::cout << "set_dbnL_design_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  cellp = &Lattice.Cell[Lattice.Elem_GetPos(Fnum, Knum)]; elemp = &cellp->Elem;

  if (cellp->L != 0.0) {
    elemp->M->Bpar[HOMmax+n] += dbnL/cellp->L;
    elemp->M->Bpar[HOMmax-n] += danL/cellp->L;
  } else {
    // thin kick
    elemp->M->Bpar[HOMmax+n] += dbnL; elemp->M->Bpar[HOMmax-n] += danL;
  }

  Mpole_SetB(Fnum, Knum, n); Mpole_SetB(Fnum, Knum, -n);
}


void set_dbnL_design_fam(const int Fnum,
			 const int n, const double dbnL, const double danL)
{
  int k;

  if (n < 1) {
    std::cout << "set_dbnL_design_fam: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  for (k = 1; k <= Lattice.GetnKid(Fnum); k++)
    set_bnL_design_elem(Fnum, k, n, dbnL, danL);
}


void set_bnL_design_fam(const int Fnum,
			const int n, const double bnL, const double anL)
{
  int k;

  if (n < 1) {
    std::cout << "set_bnL_design_fam: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  for (k = 1; k <= Lattice.GetnKid(Fnum); k++)
    set_bnL_design_elem(Fnum, k, n, bnL, anL);
}


void set_bnL_design_type(const int type,
			 const int n, const double bnL, const double anL)
{
  long int  k;

  if (n < 1) {
    std::cout << "set_bnL_design_type: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  if ((type >= Dip) && (type <= HOMmax)) {
    for (k = 1; k <= Lattice.param.Cell_nLoc; k++)
      if ((Lattice.Cell[k].Kind == Mpole) &&
	  (Lattice.Cell[k].Elem.M->n_design == type))
	set_bnL_design_elem(Lattice.Cell[k].Fnum, Lattice.Cell[k].Knum,
			    n, bnL, anL);
  } else
    printf("Bad type argument to set_bnL_design_type()\n");
}


void set_bnL_sys_elem(const int Fnum, const int Knum,
		      const int n, const double bnL, const double anL)
{
  CellType *cellp;
  elemtype *elemp;

  if (n < 1) {
    std::cout << "set_bnL_sys_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  cellp = &Lattice.Cell[Lattice.Elem_GetPos(Fnum, Knum)]; elemp = &cellp->Elem;

  if (cellp->L != 0.0) {
    elemp->M->Bsys[HOMmax+n] = bnL/cellp->L;
    elemp->M->Bsys[HOMmax-n] = anL/cellp->L;
  } else {
    // thin kick
    elemp->M->Bsys[HOMmax+n] = bnL; elemp->M->Bsys[HOMmax-n] = anL;
  }

  Mpole_SetB(Fnum, Knum, n); Mpole_SetB(Fnum, Knum, -n);

  if (trace) {
    printf("set_bnL_sys_elem: %s %3d %e %e\n",
	   cellp->Name, n, elemp->M->Bsys[HOMmax+n], elemp->M->Bsys[HOMmax-n]);
    printf("                                      %e %e\n",
	   elemp->M->B[HOMmax+n], elemp->M->B[HOMmax-n]);
  }
}


void set_bnL_sys_fam(const int Fnum,
		     const int n, const double bnL, const double anL)
{
  int k;

  if (n < 1) {
    std::cout << "set_bnL_sys_fam: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  for (k = 1; k <= Lattice.GetnKid(Fnum); k++)
    set_bnL_sys_elem(Fnum, k, n, bnL, anL);
}


void set_bnL_sys_type(const int type,
		      const int n, const double bnL, const double anL)
{
  long int   k;

  if (n < 1) {
    std::cout << "set_bnL_sys_type: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  if (type >= Dip && type <= HOMmax) {
    for(k = 1; k <= Lattice.param.Cell_nLoc; k++)
      if ((Lattice.Cell[k].Kind == Mpole)
	  && (Lattice.Cell[k].Elem.M->n_design == type))
	set_bnL_sys_elem(Lattice.Cell[k].Fnum, Lattice.Cell[k].Knum,
			 n, bnL, anL);
  } else
    printf("Bad type argument to set_bnL_sys_type()\n");
}


void set_bnL_rms_elem(const int Fnum, const int Knum,
		      const int n, const double bnL, const double anL,
		      const bool new_rnd)
{
  CellType *cellp;
  elemtype *elemp;

  bool  prt = false;

  if (n < 1) {
    std::cout << "set_bnL_rms_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  cellp = &Lattice.Cell[Lattice.Elem_GetPos(Fnum, Knum)]; elemp = &cellp->Elem;

  if (cellp->L != 0.0) {
    elemp->M->Brms[HOMmax+n] = bnL/cellp->L;
    elemp->M->Brms[HOMmax-n] = anL/cellp->L;
  } else {
    // thin kick
    elemp->M->Brms[HOMmax+n] = bnL; elemp->M->Brms[HOMmax-n] = anL;
  }

  if(new_rnd){
    if (normal) {
      elemp->M->Brnd[HOMmax+n] = normranf();
      elemp->M->Brnd[HOMmax-n] = normranf();
    } else {
      elemp->M->Brnd[HOMmax+n] = ranf(); elemp->M->Brnd[HOMmax-n] = ranf();
    }
  }

  if (prt)
    printf("set_bnL_rms_elem:  Fnum = %d, Knum = %d"
	   ", bnL = %e, anL = %e %e %e\n",
	   Fnum, Knum, bnL, anL,
	   elemp->M->Brms[HOMmax+n], elemp->M->Brms[HOMmax-n]);

  Mpole_SetB(Fnum, Knum, n); Mpole_SetB(Fnum, Knum, -n);
}


void set_bnL_rms_fam(const int Fnum,
		     const int n, const double bnL, const double anL,
		     const bool new_rnd)
{
  int k;

  if (n < 1) {
    std::cout << "set_bnL_rms_fam: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  for (k = 1; k <= Lattice.GetnKid(Fnum); k++)
    set_bnL_rms_elem(Fnum, k, n, bnL, anL, new_rnd);
}


void set_bnL_rms_type(const int type,
		      const int n, const double bnL, const double anL,
		      const bool new_rnd)
{
  long int   k;

  if (n < 1) {
    std::cout << "get_bnL_rms_type: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  if (type >= Dip && type <= HOMmax) {
    for(k = 1; k <= Lattice.param.Cell_nLoc; k++)
      if ((Lattice.Cell[k].Kind == Mpole)
	  && (Lattice.Cell[k].Elem.M->n_design == type))
	set_bnL_rms_elem(Lattice.Cell[k].Fnum, Lattice.Cell[k].Knum,
			 n, bnL, anL, new_rnd);
  } else
    printf("Bad type argument to set_bnL_rms_type()\n");
}


void set_bnr_sys_elem(const int Fnum, const int Knum,
		      const int n, const double bnr, const double anr)
{
  int        nd;
  MpoleType  *mp;
  bool prt = false;

  if (n < 1) {
    std::cout << "set_bnr_sys_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  mp = Lattice.Cell[Lattice.Elem_GetPos(Fnum, Knum)].Elem.M; nd = mp->n_design;
  // errors are relative to design values for (Dip, Quad, Sext, ...)
  mp->Bsys[HOMmax+n] = bnr*mp->Bpar[HOMmax+nd];
  mp->Bsys[HOMmax-n] = anr*mp->Bpar[HOMmax+nd];

  Mpole_SetB(Fnum, Knum, n); Mpole_SetB(Fnum, Knum, -n);

  if (prt)
    printf("set the n=%d component of %s to %e %e %e\n",
	   n, Lattice.Cell[Lattice.Elem_GetPos(Fnum, Knum)].Name,
	   bnr, mp->Bpar[HOMmax+nd], mp->Bsys[HOMmax+n]);
}


void set_bnr_sys_fam(const int Fnum,
		     const int n, const double bnr, const double anr)
{
  int k;

  if (n < 1) {
    std::cout << "set_bnr_sys_fam: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  for (k = 1; k <= Lattice.GetnKid(Fnum); k++)
    set_bnr_sys_elem(Fnum, k, n, bnr, anr);
}


void set_bnr_sys_type(const int type,
		      const int n, const double bnr, const double anr)
{
  long int   k;

  if (n < 1) {
    std::cout << "set_bnr_sys_type: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  if (type >= Dip && type <= HOMmax) {
    for(k = 1; k <= Lattice.param.Cell_nLoc; k++)
      if ((Lattice.Cell[k].Kind == Mpole)
	  && (Lattice.Cell[k].Elem.M->n_design == type))
	set_bnr_sys_elem(Lattice.Cell[k].Fnum, Lattice.Cell[k].Knum,
			 n, bnr, anr);
  } else
    printf("Bad type argument to set_bnr_sys_type()\n");
}


void set_bnr_rms_elem(const int Fnum, const int Knum,
		      const int n, const double bnr, const double anr,
		      const bool new_rnd)
{
  int        nd;
  MpoleType  *mp;

  bool  prt = false;

  if (n < 1) {
    std::cout << "set_bnr_rms_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  mp = Lattice.Cell[Lattice.Elem_GetPos(Fnum, Knum)].Elem.M; nd = mp->n_design;
  // errors are relative to design values for (Dip, Quad, Sext, ...)
  if (nd == Dip) {
    mp->Brms[HOMmax+n] = bnr*mp->irho; mp->Brms[HOMmax-n] = anr*mp->irho;
  } else {
    mp->Brms[HOMmax+n] = bnr*mp->Bpar[HOMmax+nd];
    mp->Brms[HOMmax-n] = anr*mp->Bpar[HOMmax+nd];
  }

  if(new_rnd){
    if (normal) {
      mp->Brnd[HOMmax+n] = normranf(); mp->Brnd[HOMmax-n] = normranf();
    } else {
      mp->Brnd[HOMmax+n] = ranf(); mp->Brnd[HOMmax-n] = ranf();
    }
  }

  Mpole_SetB(Fnum, Knum, n); Mpole_SetB(Fnum, Knum, -n);

  if (prt) {
    printf("set_bnr_rms_elem:  Fnum = %d, Knum = %d, n = %d, n_design = %d"
	   ", new_rnd = %d, r_# = (%e, %e)\n",
	   Fnum, Knum, n, nd, new_rnd,
	   mp->Brnd[HOMmax+n], mp->Brnd[HOMmax-n]);
    printf("  (bnr, anr) = (%e, %e), Brms = (%e, %e), B = (%e, %e)\n",
	 bnr, anr, mp->Brms[HOMmax+n], mp->Brms[HOMmax-n],
	 mp->B[HOMmax+n], mp->B[HOMmax-n]);
  }
}


void set_bnr_rms_fam(const int Fnum,
		     const int n, const double bnr, const double anr,
		     const bool new_rnd)
{
  int k;

  if (n < 1) {
    std::cout << "set_bnr_rms_fam: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  for (k = 1; k <= Lattice.GetnKid(Fnum); k++)
    set_bnr_rms_elem(Fnum, k, n, bnr, anr, new_rnd);
}


void set_bnr_rms_type(const int type,
		      const int n, const double bnr, const double anr,
		      const bool new_rnd)
{
  long int   k;

  if (n < 1) {
    std::cout << "set_bnr_rms_type: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  if (type >= Dip && type <= HOMmax) {
    for(k = 1; k <= Lattice.param.Cell_nLoc; k++)
      if ((Lattice.Cell[k].Kind == Mpole)
	  && (Lattice.Cell[k].Elem.M->n_design == type))
	set_bnr_rms_elem(Lattice.Cell[k].Fnum, Lattice.Cell[k].Knum,
			 n, bnr, anr, new_rnd);
  } else
    printf("Bad type argument to set_bnr_rms_type()\n");
}


double get_Wiggler_BoBrho(const int Fnum, const int Knum)
{
  return Lattice.Cell[Lattice.Elem_GetPos(Fnum, Knum)].Elem.W->BoBrhoV[0];
}


void set_Wiggler_BoBrho(const int Fnum, const int Knum, const double BoBrhoV)
{
  Lattice.Cell[Lattice.Elem_GetPos(Fnum, Knum)].Elem.W->BoBrhoV[0] = BoBrhoV;
  Lattice.Cell[Lattice.Elem_GetPos(Fnum, Knum)].Elem.W->BW[HOMmax+Quad]
    = -sqr(BoBrhoV)/2.0;
  Wiggler_SetB(Fnum, Knum, Quad);
}


void set_Wiggler_BoBrho(const int Fnum, const double BoBrhoV)
{
  int  k;

  for (k = 1; k <= Lattice.GetnKid(Fnum); k++)
    set_Wiggler_BoBrho(Fnum, k, BoBrhoV);
}


void set_ID_scl(const int Fnum, const int Knum, const double scl)
{
  int           k;
  WigglerType*  W;

  switch (Lattice.Cell[Lattice.Elem_GetPos(Fnum, Knum)].Kind) {
  case Wigl:
    // scale the ID field
    W = Lattice.Cell[Lattice.Elem_GetPos(Fnum, Knum)].Elem.W;
    for (k = 0; k < W->n_harm; k++) {
      W->BoBrhoH[k] = scl*Lattice.ElemFam[Fnum-1].ElemF.W->BoBrhoH[k];
      W->BoBrhoV[k] = scl*Lattice.ElemFam[Fnum-1].ElemF.W->BoBrhoV[k];
    }
    break;
  case Insertion:
    Lattice.Cell[Lattice.Elem_GetPos(Fnum, Knum)].Elem.ID->scaling = scl;
    break;
  case FieldMap:
    Lattice.Cell[Lattice.Elem_GetPos(Fnum, Knum)].Elem.FM->scl = scl;
    break;
  default:
    std::cout << "set_ID_scl: unknown element type" << std::endl;
    exit_(1);
    break;
  }
}


void set_ID_scl(const int Fnum, const double scl)
{
  int  k;

  for (k = 1; k <= Lattice.GetnKid(Fnum); k++)
    set_ID_scl(Fnum, k, scl);
}


void LatticeType::SetFieldValues_fam(const int Fnum, const bool rms,
				     const double r0, const int n,
				     const double Bn, const double An,
				     const bool new_rnd)
{
  int     N;
  double  bnr, anr;

  N = Lattice.Cell[Lattice.Elem_GetPos(Fnum, 1)].Elem.M->n_design;
  if (r0 == 0.0) {
    // input is: (b_n L), (a_n L)
    if(rms)
      set_bnL_rms_fam(Fnum, n, Bn, An, new_rnd);
    else
      set_bnL_sys_fam(Fnum, n, Bn, An);
  } else {
    bnr = Bn/pow(r0, (double)(n-N)); anr = An/pow(r0, (double)(n-N));
    if(rms)
      set_bnr_rms_fam(Fnum, n, bnr, anr, new_rnd);
    else
      set_bnr_sys_fam(Fnum, n, bnr, anr);
  }
}


void LatticeType::SetFieldValues_type(const int N, const bool rms,
				      const double r0, const int n,
				      const double Bn, const double An,
				      const bool new_rnd)
{
  double  bnr, anr;

  if (r0 == 0.0) {
    // input is: (b_n L), (a_n L)
    if(rms)
      set_bnL_rms_type(N, n, Bn, An, new_rnd);
    else
      set_bnL_sys_type(N, n, Bn, An);
  } else {
    bnr = Bn/pow(r0, (double)(n-N)); anr = An/pow(r0, (double)(n-N));
    if(rms)
      set_bnr_rms_type(N, n, bnr, anr, new_rnd);
    else
      set_bnr_sys_type(N, n, bnr, anr);
  }
}


void LatticeType::SetFieldErrors(const char *name, const bool rms,
				 const double r0,
				 const int n, const double Bn, const double An,
				 const bool new_rnd)
{
  int     Fnum;

  if (strcmp("all", name) == 0) {
    printf("all: not yet implemented\n");
  } else if (strcmp("dip", name) == 0) {
    SetFieldValues_type(Dip, rms, r0, n, Bn, An, new_rnd);
  } else if (strcmp("quad", name) == 0) {
    SetFieldValues_type(Quad, rms, r0, n, Bn, An, new_rnd);
  } else if (strcmp("sext", name) == 0) {
    SetFieldValues_type(Sext, rms, r0, n, Bn, An, new_rnd);
  } else {
    Fnum = Lattice.Elem_Index(name);
    if(Fnum > 0)
      SetFieldValues_fam(Fnum, rms, r0, n, Bn, An, new_rnd);
    else
      printf("SetFieldErrors: undefined element %s\n", name);
  }
}


// closed orbit correction by n_orbit iterations
bool LatticeType::CorrectCOD(const int n_orbit, const double scl)
{
  bool            cod;
  int             i;
  long int        lastpos;
  Vector2         mean, sigma, max;
  ss_vect<double> ps;

  // ps.zero(); Lattice.Cell_Pass(0, Lattice.param.Cell_nLoc, ps, lastpos);
  // for (i = 1; i <= n_orbit; i++) {
  //   lstc(1, lastpos); lstc(2, lastpos);

  //   ps.zero(); Lattice.Cell_Pass(0, Lattice.param.Cell_nLoc, ps, lastpos);
  // }
  // if (false) Lattice.prt_cod("cod.out", Lattice.param.bpm, true);
 
  cod = Lattice.getcod(0e0, lastpos);
  if (cod) {
    codstat(mean, sigma, max, Lattice.param.Cell_nLoc, true);
    printf("\n");
    printf("Initial RMS orbit (all):    x = %7.1e mm, y = %7.1e mm\n",
	   1e3*sigma[X_], 1e3*sigma[Y_]);
    codstat(mean, sigma, max, Lattice.param.Cell_nLoc, false);
    printf("\n");
    printf("Initial RMS orbit (BPMs):   x = %7.1e mm, y = %7.1e mm\n",
	   1e3*sigma[X_], 1e3*sigma[Y_]);

    for (i = 1; i <= n_orbit; i++){
      lsoc(1, scl); lsoc(2, scl);
      cod = Lattice.getcod(0e0, lastpos);
      if (!cod) break;

      if (cod) {
	codstat(mean, sigma, max, Lattice.param.Cell_nLoc, false);
	printf("Corrected RMS orbit (BPMs): x = %7.1e mm, y = %7.1e mm\n",
	       1e3*sigma[X_], 1e3*sigma[Y_]);
	if (i == n_orbit) {
	  codstat(mean, sigma, max, Lattice.param.Cell_nLoc, true);
	  printf("\n");
	  printf("Corrected RMS orbit (all):  x = %7.1e mm, y = %7.1e mm\n",
		 1e3*sigma[X_], 1e3*sigma[Y_]);
	}
      }
    }
  }

  return cod;
}


void LatticeType::prt_beamsizes()
{
  int  k;
  FILE *fp;

  fp = file_write(beam_envelope_file);

  fprintf(fp,"# k    name    s    s_xx    s_pxpx    s_xpx    s_yy    s_pypy    "
	  "s_ypy    theta_xy\n");
  for(k = 0; k <= Lattice.param.Cell_nLoc; k++)
    fprintf(fp,"%4d %10s %e %e %e %e %e %e %e %e\n",
	    k, Lattice.Cell[k].Name, Lattice.Cell[k].S,
	    Lattice.Cell[k].sigma[x_][x_], Lattice.Cell[k].sigma[px_][px_],
	    Lattice.Cell[k].sigma[x_][px_],
	    Lattice.Cell[k].sigma[y_][y_], Lattice.Cell[k].sigma[py_][py_],
	    Lattice.Cell[k].sigma[y_][py_],
	    atan2(2e0*Lattice.Cell[k].sigma[x_][y_],
		  Lattice.Cell[k].sigma[x_][x_]-Lattice.Cell[k].sigma[y_][y_])
	    /2e0*180.0/M_PI);

  fclose(fp);
}


double f_int_Touschek(const double u)
{
  double  f;

  if (u > 0.0)
    f = (1.0/u-log(1.0/u)/2.0-1.0)*exp(-u_Touschek/u);
  else
    f = 0.0;

  return f;
}


double Touschek_loc(const long int i, const double gamma,
		    const double delta_RF,
		    const double eps_x, const double eps_y,
		    const double sigma_delta, const double sigma_s,
		    const bool ZAP_BS)
{
  int     k;
  double  sigma_x, sigma_y, sigma_xp, curly_H, dtau_inv;
  double  alpha[2], beta[2], eta[2], etap[2];

  if ((i < 0) || (i > Lattice.param.Cell_nLoc)) {
    std::cout << "Touschek_loc: undefined location " << i << std::endl;
    exit(1);
  }

  if (!ZAP_BS) {
    curly_H = get_curly_H(Lattice.Cell[i].Alpha[X_], Lattice.Cell[i].Beta[X_],
			  Lattice.Cell[i].Eta[X_], Lattice.Cell[i].Etap[X_]);

    // Compute beam sizes for given hor/ver emittance, sigma_s,
    // and sigma_delta (for x ~ 0): sigma_x0' = sqrt(eps_x/beta_x).
    sigma_x = sqrt(Lattice.Cell[i].Beta[X_]*eps_x+sqr(Lattice.Cell[i].Eta[X_]
						      *sigma_delta));
    sigma_y = sqrt(Lattice.Cell[i].Beta[Y_]*eps_y);
    sigma_xp = (eps_x/sigma_x)*sqrt(1e0+curly_H*sqr(sigma_delta)/eps_x);
  } else {
    // ZAP averages the optics functions over an element instead of the
    // integrand; incorrect.

    for (k = 0; k < 2; k++) {
      alpha[k] = (Lattice.Cell[i-1].Alpha[k]+Lattice.Cell[i].Alpha[k])/2e0;
      beta[k] = (Lattice.Cell[i-1].Beta[k]+Lattice.Cell[i].Beta[k])/2e0;
      eta[k] = (Lattice.Cell[i-1].Eta[k]+Lattice.Cell[i].Eta[k])/2e0;
      etap[k] = (Lattice.Cell[i-1].Etap[k]+Lattice.Cell[i].Etap[k])/2e0;
    }

    curly_H = get_curly_H(alpha[X_], beta[X_], eta[X_], etap[X_]);

    // Compute beam sizes for given hor/ver emittance, sigma_s,
    // and sigma_delta (for x ~ 0): sigma_x0' = sqrt(eps_x/beta_x).
    sigma_x = sqrt(beta[X_]*eps_x+sqr(eta[X_]*sigma_delta));
    sigma_y = sqrt(beta[Y_]*eps_y);
    sigma_xp = (eps_x/sigma_x)*sqrt(1e0+curly_H*sqr(sigma_delta)/eps_x);
  }

  u_Touschek = sqr(delta_RF/(gamma*sigma_xp));

  dtau_inv = dqromb(f_int_Touschek, 0e0, 1e0)/(sigma_x*sigma_y*sigma_xp);

  return dtau_inv;
}


double LatticeType::Touschek(const double Qb, const double delta_RF,
			      const double eps_x, const double eps_y,
			      const double sigma_delta, const double sigma_s)
{
  // Note, ZAP (LBL-21270) averages the optics functions over an element
  // instead of the integrand; incorrect.  Hence, the Touschek lifetime is
  // overestimated by ~20%.

  long int  i;
  double    p1, p2, dtau_inv, tau_inv;

  const bool    ZAP_BS = false;
  const double  gamma = 1e9*Lattice.param.Energy/m_e, N_e = Qb/q_e;

  printf("\n");
  printf("Qb = %4.2f nC, delta_RF = %4.2f%%"
	 ", eps_x = %9.3e m.rad, eps_y = %9.3e m.rad\n",
	 1e9*Qb, 1e2*delta_RF, eps_x, eps_y);
  printf("sigma_delta = %8.2e, sigma_s = %4.2f mm\n",
	 sigma_delta, 1e3*sigma_s);

  // Integrate around the lattice with Trapezoidal rule; for nonequal steps.

  p1 = Touschek_loc(0, gamma, delta_RF, eps_x, eps_y, sigma_delta, sigma_s,
		    ZAP_BS);

  tau_inv = 0e0;
  for(i = 1; i <= Lattice.param.Cell_nLoc; i++) {
    p2 = Touschek_loc(i, gamma, delta_RF, eps_x, eps_y, sigma_delta, sigma_s,
		      ZAP_BS);

    if (!ZAP_BS)
      dtau_inv = (p1+p2)/2e0;
    else
      dtau_inv = p2;

    tau_inv += dtau_inv*Lattice.Cell[i].L; p1 = p2;

    if (false) {
      dtau_inv *=
	N_e*sqr(r_e)*c0/(8.0*M_PI*cube(gamma)*sigma_s)/sqr(delta_RF);

      printf("%4ld %9.3e\n", i, dtau_inv);
    }
  }

  tau_inv *=
    N_e*sqr(r_e)*c0/(8.0*M_PI*cube(gamma)*sigma_s)
    /(sqr(delta_RF)*Lattice.Cell[Lattice.param.Cell_nLoc].S);

  printf("\n");
  printf("Touschek lifetime [hrs]: %10.3e\n", 1e0/(3600e0*tau_inv));

  return(1e0/(tau_inv));
}


void mom_aper(double &delta, double delta_RF, const long int k,
	      const int n_turn, const bool positive)
{
  // Binary search to determine momentum aperture at location k.
  int       j;
  long int  lastpos;
  double    delta_min, delta_max;
  psVector    x;

  const double  eps = 1e-4;

  delta_min = 0.0; delta_max = positive ? fabs(delta_RF) : -fabs(delta_RF);
  while (fabs(delta_max-delta_min) > eps) {
    delta = (delta_max+delta_min)/2.0;

    // propagate initial conditions
    CopyVec(6, Lattice.param.CODvect, x); Lattice.Cell_Pass(0, k, x, lastpos);
    // generate Touschek event
    x[delta_] += delta;

    // complete one turn
    Lattice.Cell_Pass(k+1, Lattice.param.Cell_nLoc, x, lastpos);
    if (lastpos < Lattice.param.Cell_nLoc)
      // particle lost
      delta_max = delta;
    else {
      // track
      for(j = 0; j < n_turn; j++) {
	Lattice.Cell_Pass(0, Lattice.param.Cell_nLoc, x, lastpos);

	if ((delta_max > delta_RF) || (lastpos < Lattice.param.Cell_nLoc)) {
	  // particle lost
	  delta_max = delta;
	  break;
	}
      }

      if ((delta_max <= delta_RF) && (lastpos == Lattice.param.Cell_nLoc))
	// particle not lost
	delta_min = delta;
    }
  }
}


double LatticeType::Touschek(const double Qb, const double delta_RF,
			      const bool consistent,
			      const double eps_x, const double eps_y,
			      const double sigma_delta, double sigma_s,
			      const int n_turn, const bool aper_on,
			      double sum_delta[][2], double sum2_delta[][2])
{
  bool     cav, aper;
  long int k;
  double   tau_inv, delta_p, delta_m, curly_H0, curly_H1, L;
  double   sigma_x, sigma_y, sigma_xp;

  const bool prt = false;

  //  const char  file_name[] = "Touschek.out";
  const double  eps = 1e-12, gamma = 1e9*Lattice.param.Energy/m_e, N_e = Qb/q_e;

  cav = Lattice.param.Cavity_on; aper = Lattice.param.Aperture_on;

  Lattice.param.Cavity_on = true;

  Lattice.Ring_GetTwiss(true, 0.0);

  Lattice.param.Aperture_on = aper_on;

  printf("\nQb = %4.2f nC, delta_RF = %4.2f%%"
	 ", eps_x = %9.3e m.rad, eps_y = %9.3e m.rad\n",
	 1e9*Qb, 1e2*delta_RF, eps_x, eps_y);
  printf("sigma_delta = %8.2e, sigma_s = %4.2f mm\n",
	 sigma_delta, 1e3*sigma_s);

  printf("\nMomentum aperture:\n");

  delta_p = delta_RF; mom_aper(delta_p, delta_RF, 0, n_turn, true);
  delta_m = -delta_RF; mom_aper(delta_m, delta_RF, 0, n_turn, false);
  delta_p = min(delta_RF, delta_p); delta_m = max(-delta_RF, delta_m);
  sum_delta[0][0] += delta_p; sum_delta[0][1] += delta_m;
  sum2_delta[0][0] += sqr(delta_p); sum2_delta[0][1] += sqr(delta_m);

  tau_inv = 0e0; curly_H0 = -1e30;
  for (k = 1; k <= Lattice.param.Cell_nLoc; k++) {
    L = Lattice.Cell[k].L;

    curly_H1 = get_curly_H(Lattice.Cell[k].Alpha[X_], Lattice.Cell[k].Beta[X_],
			   Lattice.Cell[k].Eta[X_], Lattice.Cell[k].Etap[X_]);

    if (fabs(curly_H0-curly_H1) > eps) {
      mom_aper(delta_p, delta_RF, k, n_turn, true);
      delta_m = -delta_p; mom_aper(delta_m, delta_RF, k, n_turn, false);
      delta_p = min(delta_RF, delta_p); delta_m = max(-delta_RF, delta_m);
      printf("%4ld %6.2f %3.2lf%% %3.2lf%%\n",
	     k, Lattice.Cell[k].S, 1e2*delta_p, 1e2*delta_m);
      curly_H0 = curly_H1;
    }

    sum_delta[k][0] += delta_p; sum_delta[k][1] += delta_m;
    sum2_delta[k][0] += sqr(delta_p); sum2_delta[k][1] += sqr(delta_m);
    if (prt)
      printf("%4ld %6.2f %3.2lf %3.2lf\n",
	     k, Lattice.Cell[k].S, 1e2*delta_p, 1e2*delta_m);

    if (!consistent) {
      // Compute beam sizes for given hor/ver emittance, sigma_s,
      // and sigma_delta (for x ~ 0): sigma_x0' = sqrt(eps_x/beta_x).
      sigma_x = sqrt(Lattice.Cell[k].Beta[X_]*eps_x
		     +sqr(sigma_delta*Lattice.Cell[k].Eta[X_]));
      sigma_y = sqrt(Lattice.Cell[k].Beta[Y_]*eps_y);
      sigma_xp = (eps_x/sigma_x)*sqrt(1e0+curly_H1*sqr(sigma_delta)/eps_x);
    } else {
      // use self-consistent beam sizes
      sigma_x = sqrt(Lattice.Cell[k].sigma[x_][x_]);
      sigma_y = sqrt(Lattice.Cell[k].sigma[y_][y_]);
      sigma_xp = sqrt(Lattice.Cell[k].sigma[px_][px_]);
      sigma_s = sqrt(Lattice.Cell[k].sigma[ct_][ct_]);
    }

    u_Touschek = sqr(delta_p/(gamma*sigma_xp));

    tau_inv +=
      dqromb(f_int_Touschek, 0e0, 1e0)
      /(sigma_x*sigma_xp*sigma_y*sqr(delta_p))*L;

    fflush(stdout);
  }

  tau_inv *=
    N_e*sqr(r_e)*c0/(8.0*M_PI*cube(gamma)*sigma_s)
    /Lattice.Cell[Lattice.param.Cell_nLoc].S;

  printf("\n");
  printf("Touschek lifetime [hrs]: %4.2f\n", 1e0/(3600e0*tau_inv));

  Lattice.param.Cavity_on = cav; Lattice.param.Aperture_on = aper;

  return 1/tau_inv;
}


double LatticeType::f_IBS(const double chi_m)
{
  // Interpolated integral (V. Litvinenko).

  double  f, ln_chi_m;

  const double A = 1.0, B = 0.579, C = 0.5;

  ln_chi_m = log(chi_m); f = A + B*ln_chi_m + C*sqr(ln_chi_m);

  return f;
}


double f_int_IBS(const double chi)
{
  // Integrand for numerical integration.

  return exp(-chi)*log(chi/chi_m)/chi;
}


double get_int_IBS(void)
{
  int     k;
  double  f;

  const int     n_decades = 30;
  const double  base      = 10e0;

  f = 0e0;
  for (k = 0; k <= n_decades; k++) {
    if (k == 0)
      f += dqromo(f_int_IBS, chi_m, pow(base, (double)k), dmidsql);
    else
      f += dqromb(f_int_IBS, pow(base, (double)(k-1)), pow(base, (double)k));
  }

  return f;
}


void LatticeType::IBS(const double Qb, const double eps_SR[], double eps[],
		       const bool prt1, const bool prt2)
{
  /* J. Le Duff "Single and Multiple Touschek Effects" (e.g. CERN 89-01)

     The equilibrium emittance is given by (tau is for amplitudes not
     invariants):

       sigma_delta^2 = sigma_delta_SR^2 + D_delta_IBS*tau_delta/2

     and

       D_x = <D_delta*curly_H> =>

       eps_x = eps_x_SR + eps_x_IBS = eps_x_SR + D_x_IBS*tau_x/2

     where

       D_x_IBS ~ 1/eps_x

     Multiplying with eps_x gives

       eps_x^2 = eps_x*eps_x_SR + (eps_x*D_x_IBS)*tau_x/2

     where the 2nd term is roughly a constant. The IBS limit is obtained by

       eps_x_SR -> 0 => eps_x_IBS = sqrt((eps_x*D_x_IBS)*tau_x/2)

     which is inserted into the original equation

       eps_x^2 = eps_x*eps_x_SR + eps_x_IBS^2

     and then solved for eps_x

       eps_x = eps_x_SR/2 + sqrt(eps_x_SR^2/4+eps_x_IBS^2)                   */

  long int  k;
  double    D_x, D_delta, b_max, L, gamma_z, a;
  double    sigma_x, sigma_xp, sigma_y, sigma_s, sigma_delta;
  double    incr, curly_H, eps_IBS[3];
  double    sigma_s_SR, sigma_delta_SR;

  const bool    integrate = false;
  const double  gamma = 1e9*Lattice.param.Energy/m_e, N_b = Qb/q_e;

  // bunch size
  gamma_z = (1.0+sqr(Lattice.param.alpha_z))/Lattice.param.beta_z;

  sigma_s_SR = sqrt(Lattice.param.beta_z*eps_SR[Z_]);
  sigma_delta_SR = sqrt(gamma_z*eps_SR[Z_]);

  sigma_s = sqrt(Lattice.param.beta_z*eps[Z_]);
  sigma_delta = sqrt(gamma_z*eps[Z_]);

  if (prt1) {
    printf("\nQb             = %4.2f nC,        Nb          = %9.3e\n",
	   1e9*Qb, N_b);
    printf("eps_x_SR       = %7.3f pm.rad, eps_x       = %7.3f pm.rad\n",
	   1e12*eps_SR[X_], 1e12*eps[X_]);
    printf("eps_y_SR       = %7.3f pm.rad, eps_y       = %7.3f pm.rad\n",
	   1e12*eps_SR[Y_], 1e12*eps[Y_]);
    printf("eps_z_SR       = %9.3e,      eps_z       = %9.3e\n",
	   eps_SR[Z_], eps[Z_]);
    printf("alpha_z        = %9.3e,      beta_z      = %9.3e\n",
	   Lattice.param.alpha_z, Lattice.param.beta_z);
    printf("sigma_s_SR     = %9.3e mm,   sigma_s     = %9.3e mm\n",
	   1e3*sigma_s_SR, 1e3*sigma_s);
    printf("sigma_delta_SR = %9.3e,      sigma_delta = %9.3e\n",
	   sigma_delta_SR, sigma_delta);
  }

  D_delta = 0.0; D_x = 0.0;
  for(k = 0; k <= Lattice.param.Cell_nLoc; k++) {
    L = Lattice.Cell[k].L;

    curly_H = get_curly_H(Lattice.Cell[k].Alpha[X_], Lattice.Cell[k].Beta[X_],
			  Lattice.Cell[k].Eta[X_], Lattice.Cell[k].Etap[X_]);

    // Compute beam sizes for given hor/ver emittance, sigma_s,
    // and sigma_delta (for x ~ 0): sigma_x0' = sqrt(eps_x/beta_x).
    sigma_x = sqrt(Lattice.Cell[k].Beta[X_]*eps[X_]
		   +sqr(Lattice.Cell[k].Eta[X_]*sigma_delta));
    sigma_xp = (eps[X_]/sigma_x)*sqrt(1.0+curly_H*sqr(sigma_delta)/eps[X_]);
    sigma_y = sqrt(Lattice.Cell[k].Beta[Y_]*eps[Y_]);

    b_max = 2.0*sqrt(M_PI)/pow(N_b/(sigma_x*sigma_y*sigma_s), 1.0/3.0);

    chi_m = r_e/(b_max*sqr(gamma*sigma_xp));

    if (!integrate)
      incr = f_IBS(chi_m)/sigma_y;
    else
      incr = get_int_IBS()/sigma_y;

    D_delta += incr*L; D_x += incr*curly_H*L;
  }

  a =
    N_b*sqr(r_e)*c0
    /(32.0*M_PI*cube(gamma)*sigma_s*Lattice.Cell[Lattice.param.Cell_nLoc].S);

  // eps_x*D_X
  D_x *= a;
  // Compute eps_IBS.
  eps_IBS[X_] = sqrt(D_x*Lattice.param.tau[X_]/2e0);
  // Solve for eps_x
  eps[X_] = eps_SR[X_]*(1.0+sqrt(1.0+4.0*sqr(eps_IBS[X_]/eps_SR[X_])))/2.0;

  // compute diffusion coeffs.
  D_x /= eps[X_]; D_delta *= a/eps[X_];

  eps[Y_] = eps_SR[Y_]/eps_SR[X_]*eps[X_];

  sigma_delta = sqrt(sqr(sigma_delta_SR)+D_delta*Lattice.param.tau[Z_]/2e0);
  eps[Z_] = sqr(sigma_delta)/gamma_z;

  sigma_s = sqrt(Lattice.param.beta_z*eps[Z_]);

  if (prt2) {
    printf("\nD_x         = %9.3e\n", D_x);
    printf("eps_x(IBS)  = %7.3f pm.rad\n", 1e12*eps_IBS[X_]);
    printf("eps_x       = %7.3f pm.rad\n", 1e12*eps[X_]);
    printf("eps_y       = %7.3f pm.rad\n", 1e12*eps[Y_]);
    printf("D_delta     = %9.3e\n", D_delta);
    printf("eps_z       = %9.3e\n", eps[Z_]);
    printf("sigma_s     = %9.3e mm\n", 1e3*sigma_s);
    printf("sigma_delta = %9.3e\n", sigma_delta);
  }
}


double f_int_IBS_BM(const double lambda)
{
  double  f;

  f =
    sqrt(lambda)*(a_k_IBS*lambda+b_k_IBS)
    /pow(cube(lambda)+a_IBS*sqr(lambda)+b_IBS*lambda+c_IBS, 3.0/2.0);

  return f;
}


double get_int_IBS_BM(void)
{
  int     k;
  double  f;

  const int     n    = 30;
  const double  decade = 10e0;

  f = 0e0;
  for (k = 0; k <= n; k++) {
    if (k == 0)
      f += dqromb(f_int_IBS_BM, 0e0, pow(decade, (double)k));
    else
      f += dqromo(f_int_IBS_BM, pow(decade, (double)(k-1)),
		  pow(decade, (double)k), dmidsql);
  }

  return f;
}


void LatticeType::IBS_BM(const double Qb, const double eps_SR[], double eps[],
			  const bool prt1, const bool prt2)
{
  // J. Bjorken, S. K. Mtingwa "Intrabeam Scattering" Part. Accel. 13, 115-143
  // (1983).
  // M. Conte, M. Martini "Intrabeam Scattering in the CERN Antiproton
  // Accumulator" Par. Accel. 17, 1-10 (1985).
  // F. Zimmermann "Intrabeam Scattering with Non-Ultrarelatvistic Corrections
  // and Vertical Dispersion" CERN-AB-2006-002.
  // Note, ZAP (LBL-21270) uses the Conte-Martini formalism.  However, it also
  // averages the optics functions over an element instead of the integrand;
  // incorrect.  Hence, the IBS effect is underestimated.

  long int  k;
  int       i;
  double    gamma_z, sigma_s_SR, sigma_delta_SR, sigma_s, sigma_delta;
  double    V, beta_m[2], sigma_m[2], alpha[2], beta[2], eta[2], etap[2];
  double    T_trans, rho, lambda_D, r_max;
  double    r_min, r_min_Cl, r_min_QM, log_Coulomb;
  double    L, curly_H[2], phi[2], dtau_inv[3], tau_inv[3], Gamma;
  double    a_BM, b_BM, c_BM, a2_CM, b2_CM, D_CM;
  double    D_x, D_delta, eps_IBS[3];

  const bool    ZAP_BS = false;
  const int     model = 2; // 1: Bjorken-Mtingwa, 2: Conte-Martini, 3: MAD-X
  const double  gamma = 1e9*Lattice.param.Energy/m_e;
  const double  beta_rel = sqrt(1e0-1e0/sqr(gamma));
  const double  N_b = Qb/q_e, q_i = 1e0;

  // bunch size
  gamma_z = (1e0+sqr(Lattice.param.alpha_z))/Lattice.param.beta_z;

  sigma_s_SR = sqrt(Lattice.param.beta_z*eps_SR[Z_]);
  sigma_delta_SR = sqrt(gamma_z*eps_SR[Z_]);

  sigma_s = sqrt(Lattice.param.beta_z*eps[Z_]); sigma_delta = sqrt(gamma_z*eps[Z_]);

  if (prt1) {
    printf("\nQb             = %4.2f nC,        Nb          = %9.3e\n",
	   1e9*Qb, N_b);
    printf("eps_x_SR       = %7.3f pm.rad, eps_x       = %7.3f pm.rad\n",
	   1e12*eps_SR[X_], 1e12*eps[X_]);
    printf("eps_y_SR       = %7.3f pm.rad, eps_y       = %7.3f pm.rad\n",
	   1e12*eps_SR[Y_], 1e12*eps[Y_]);
    printf("eps_z_SR       = %9.3e,      eps_z       = %9.3e\n",
	   eps_SR[Z_], eps[Z_]);
    printf("alpha_z        = %9.3e,      beta_z      = %9.3e\n",
	   Lattice.param.alpha_z, Lattice.param.beta_z);
    printf("sigma_s_SR     = %9.3e mm,   sigma_s     = %9.3e mm\n",
	   1e3*sigma_s_SR, 1e3*sigma_s);
    printf("sigma_delta_SR = %9.3e,      sigma_delta = %9.3e\n",
	   sigma_delta_SR, sigma_delta);
  }

  // Compute the Coulomb log

  for (i = 0; i < 2; i++)
    beta_m[i] = 0e0;

  for(k = 0; k <= Lattice.param.Cell_nLoc; k++)
    for (i = 0; i < 2; i++)
      beta_m[i] += Lattice.Cell[k].Beta[i]*Lattice.Cell[k].L;

  for (i = 0; i < 2; i++) {
    beta_m[i] /= Lattice.Cell[Lattice.param.Cell_nLoc].S;
    sigma_m[i] = sqrt(beta_m[i]*eps[i]);
  }

  V = 8e0*pow(M_PI, 3e0/2e0)*sigma_m[X_]*sigma_m[Y_]*sigma_s;
  rho = N_b/V;
  T_trans = (gamma*1e9*Lattice.param.Energy-m_e)*eps[X_]/beta_m[X_];
  lambda_D = 743.4e-2/q_i*sqrt(T_trans/(1e-6*rho));
  r_max = min(sigma_m[X_], lambda_D);

  r_min_Cl = 1.44e-9*sqr(q_i)/T_trans;

  if (!ZAP_BS)
    r_min_QM = 1.973e-13/(2e0*sqrt(T_trans*m_e));
  else
    // Bug in ZAP.
    r_min_QM = 1.973e-13/(2e0*sqrt(1e-12*T_trans*m_e));

  r_min = max(r_min_Cl, r_min_QM);

  log_Coulomb = log(r_max/r_min);

  for (i = 0; i < 3; i++)
    tau_inv[i] = 0e0;

  for(k = 1; k <= Lattice.param.Cell_nLoc; k++) {
    L = Lattice.Cell[k].L;

    for (i = 0; i < 2; i++){
      if (!ZAP_BS) {
	alpha[i] = Lattice.Cell[k].Alpha[i]; beta[i] = Lattice.Cell[k].Beta[i];
	eta[i] = Lattice.Cell[k].Eta[i]; etap[i] = Lattice.Cell[k].Etap[i];
      } else {
	// Note, ZAP averages the optics functions over an element instead of
	// the integrand; incorrect.
	alpha[i] = (Lattice.Cell[k-1].Alpha[i]+Lattice.Cell[k].Alpha[i])/2e0;
	beta[i] = (Lattice.Cell[k-1].Beta[i]+Lattice.Cell[k].Beta[i])/2e0;
	eta[i] = (Lattice.Cell[k-1].Eta[i]+Lattice.Cell[k].Eta[i])/2e0;
	etap[i] = (Lattice.Cell[k-1].Etap[i]+Lattice.Cell[k].Etap[i])/2e0;
      }

      curly_H[i] = get_curly_H(alpha[i], beta[i], eta[i], etap[i]);

      phi[i] = etap[i] + alpha[i]*eta[i]/beta[i];
    }

    Gamma = eps[X_]*eps[Y_]*sigma_delta*sigma_s;

    if ((model == 1) || (model == 2)) {
      a_BM = sqr(gamma)*(curly_H[X_]/eps[X_]+1e0/sqr(sigma_delta));

      b_BM =
	(beta[X_]/eps[X_]+beta[Y_]/eps[Y_])*sqr(gamma)
	*(sqr(eta[X_])/(eps[X_]*beta[X_])+1e0/sqr(sigma_delta))
	+ beta[X_]*beta[Y_]/(eps[X_]*eps[Y_])*sqr(gamma*phi[X_]);

      c_BM =
	beta[X_]*beta[Y_]/(eps[X_]*eps[Y_])*sqr(gamma)
	*(sqr(eta[X_])/(eps[X_]*beta[X_])+1e0/sqr(sigma_delta));
    }

    switch (model) {
    case 1:
      // Bjorken-Mtingwa
      a_IBS = a_BM; b_IBS = b_BM; c_IBS = c_BM;
      break;
    case 2:
      // Conte-Martini
      a_IBS = a_BM + beta[X_]/eps[X_] + beta[Y_]/eps[Y_];

      b_IBS = b_BM + beta[X_]*beta[Y_]/(eps[X_]*eps[Y_]);

      c_IBS = c_BM;
      break;
    case 3:
      // MAD-X
      a_IBS =
	sqr(gamma)
	*(curly_H[X_]/eps[X_]+curly_H[Y_]/eps[Y_]+1e0/sqr(sigma_delta))
	+ beta[X_]/eps[X_] + beta[Y_]/eps[Y_];

      b_IBS =
	(beta[X_]/eps[X_]+beta[Y_]/eps[Y_])*sqr(gamma)
	*(sqr(eta[X_])/(eps[X_]*beta[X_])
	  +sqr(eta[Y_])/(eps[Y_]*beta[Y_])
	  +1e0/sqr(sigma_delta))
	+ beta[X_]*beta[Y_]/(eps[X_]*eps[Y_])*sqr(gamma)
	*(sqr(phi[X_])+sqr(phi[Y_])+1e0/sqr(gamma));

      c_IBS =
	beta[X_]*beta[Y_]/(eps[X_]*eps[Y_])*sqr(gamma)
	*(sqr(eta[X_])/(eps[X_]*beta[X_])
	  +sqr(eta[Y_])/(eps[Y_]*beta[Y_])
	  +1e0/sqr(sigma_delta));
      break;
    }

    // Horizontal plane.
    switch (model) {
    case 1:
      a_k_IBS = 2e0*a_BM; b_k_IBS = b_BM;
      break;
    case 2:
      D_CM =
	sqr(gamma)
	*(sqr(eta[X_])/(beta[X_]*eps[X_])+beta[X_]/eps[X_]*sqr(phi[X_]));

      a2_CM =
	beta[X_]/eps[X_]
	*(6e0*beta[X_]/eps[X_]*sqr(gamma*phi[X_])-a_BM
	  +2e0*beta[X_]/eps[X_]-beta[Y_]/eps[Y_])/D_CM;

      b2_CM =
	beta[X_]/eps[X_]
	*(6e0*beta[X_]*beta[Y_]/(eps[X_]*eps[Y_])*sqr(gamma*phi[X_])+b_IBS
	  -3e0*beta[Y_]/eps[Y_]*a_BM)/D_CM;

      a_k_IBS =	2e0*a_BM - beta[X_]/eps[X_] - beta[Y_]/eps[Y_] + a2_CM;

      b_k_IBS = b_IBS - 3e0*beta[X_]*beta[Y_]/(eps[X_]*eps[Y_]) + b2_CM;
      break;
    case 3:
      a_k_IBS =
	2e0*sqr(gamma)
	*(curly_H[X_]/eps[X_]+curly_H[Y_]/eps[Y_]+1e0/sqr(sigma_delta))
	- beta[X_]*curly_H[Y_]/(curly_H[X_]*eps[Y_])
	+ beta[X_]/(curly_H[X_]*sqr(gamma))
	*(2e0*beta[X_]/eps[X_]-beta[Y_]/eps[Y_]-sqr(gamma/sigma_delta));

      b_k_IBS =
	(beta[X_]/eps[X_]+beta[Y_]/eps[Y_])*sqr(gamma)
	*(curly_H[X_]/eps[X_]+curly_H[Y_]/eps[Y_]+1e0/sqr(sigma_delta))
	- sqr(gamma)
	*(sqr(beta[X_]*phi[X_]/eps[X_])+sqr(beta[Y_]*phi[Y_]/eps[Y_]))
	+ (beta[X_]/eps[X_]-4e0*beta[Y_]/eps[Y_])*beta[X_]/eps[X_]
	+ beta[X_]/(sqr(gamma)*curly_H[X_])
	*(sqr(gamma/sigma_delta)
	  *(beta[X_]/eps[X_]-2e0*beta[Y_]/eps[Y_])
	  +beta[X_]*beta[Y_]/(eps[X_]*eps[Y_])
	  +sqr(gamma)*(2e0*sqr(beta[Y_]*phi[Y_]/eps[Y_])
		       -sqr(beta[X_]*phi[X_]/eps[X_])))
	+ beta[X_]*curly_H[Y_]/(eps[Y_]*curly_H[X_])
	*(beta[X_]/eps[X_]-2e0*beta[Y_]/eps[Y_]);
      break;
    }

    dtau_inv[X_] = sqr(gamma)*curly_H[X_]/eps[X_]*get_int_IBS_BM()/Gamma;
    tau_inv[X_] += dtau_inv[X_]*L;

    // Vertical plane.
    switch (model) {
    case 1:
      a_k_IBS = -a_BM; b_k_IBS = b_BM - 3e0*eps[Y_]/beta[Y_]*c_BM;
      break;
    case 2:
      a_k_IBS =	-(a_BM+beta[X_]/eps[X_]-2e0*beta[Y_]/eps[Y_]);

      b_k_IBS = b_IBS - 3e0*eps[Y_]/beta[Y_]*c_IBS;
      break;
    case 3:
      a_k_IBS =
	-sqr(gamma)*(curly_H[X_]/eps[X_]+2e0*curly_H[Y_]/eps[Y_]
		     +beta[X_]/beta[Y_]*curly_H[Y_]/eps[X_]
		     +1e0/sqr(sigma_delta))
	+ 2e0*pow(gamma, 4e0)*curly_H[Y_]/beta[Y_]
	*(curly_H[Y_]/eps[Y_]+curly_H[X_]/eps[X_])
	+ 2e0*pow(gamma, 4e0)*curly_H[Y_]/(beta[Y_]*sqr(sigma_delta))
	- (beta[X_]/eps[X_]-2e0*beta[Y_]/eps[Y_]);

      b_k_IBS =
	sqr(gamma)*(beta[Y_]/eps[Y_]-2e0*beta[X_]/eps[X_])
	*(curly_H[X_]/eps[X_]+1e0/sqr(sigma_delta))
	+ (beta[Y_]/eps[Y_]-4e0*beta[X_]/eps[X_])*sqr(gamma)
	*curly_H[Y_]/eps[Y_]
	+ beta[X_]*beta[Y_]/(eps[X_]*eps[Y_])
	+ sqr(gamma)
	*(2e0*sqr(beta[X_]*phi[X_]/eps[X_])-sqr(beta[Y_]*phi[Y_]/eps[Y_]))
	+ pow(gamma, 4e0)*curly_H[Y_]/beta[Y_]
	*(beta[X_]/eps[X_]+beta[Y_]/eps[Y_])
	*(curly_H[Y_]/eps[Y_]+1e0/sqr(sigma_delta))
	+ (beta[Y_]/eps[Y_]+beta[X_]/eps[X_])
	*pow(gamma, 4e0)*curly_H[X_]*curly_H[Y_]/(beta[Y_]*eps[X_])
	- pow(gamma, 4e0)*curly_H[Y_]/beta[Y_]
	*(sqr(beta[X_]*phi[X_]/eps[X_])+sqr(beta[Y_]*phi[Y_]/eps[Y_]));
      break;
    }

    dtau_inv[Y_] = beta[Y_]/eps[Y_]*get_int_IBS_BM()/Gamma;
    tau_inv[Y_] += dtau_inv[Y_]*L;

    // Longitudinal plane.
    switch (model) {
    case 1:
      a_k_IBS =	2e0*a_BM; b_k_IBS = b_BM;
      break;
    case 2:
      a_k_IBS =	2e0*a_BM - beta[X_]/eps[X_] - beta[Y_]/eps[Y_];

      b_k_IBS = b_BM - 2e0*beta[X_]*beta[Y_]/(eps[X_]*eps[Y_]);
      break;
    case 3:
      a_k_IBS =
	2e0*sqr(gamma)
	*(curly_H[X_]/eps[X_]+curly_H[Y_]/eps[Y_]+1e0/sqr(sigma_delta))
	- beta[X_]/eps[X_] - beta[Y_]/eps[Y_];

      b_k_IBS =
	(beta[X_]/eps[X_]+beta[Y_]/eps[Y_])*sqr(gamma)
	*(curly_H[X_]/eps[X_]+curly_H[Y_]/eps[Y_]+1e0/sqr(sigma_delta))
	- 2e0*beta[X_]*beta[Y_]/(eps[X_]*eps[Y_])
	- sqr(gamma)*(sqr(beta[X_]*phi[X_]/eps[X_])
		      +sqr(beta[Y_]*phi[Y_]/eps[Y_]));
      break;
    }

    dtau_inv[Z_] = sqr(gamma/sigma_delta)*get_int_IBS_BM()/Gamma;
    tau_inv[Z_] += dtau_inv[Z_]*L;

    if (false) {
      for (i = 0; i < 3; i++)
	dtau_inv[i] *=
	  sqr(r_e)*c0*N_b*log_Coulomb
	  /(M_PI*cube(2e0*beta_rel)*pow(gamma, 4e0));

      printf("%4ld %10.3e %10.3e %10.3e %5.3f\n",
      	     k, dtau_inv[Z_], dtau_inv[X_], dtau_inv[Y_], L);
    }
  }

  for (i = 0; i < 3; i++)
    tau_inv[i] *=
      sqr(r_e)*c0*N_b*log_Coulomb
      /(M_PI*cube(2e0*beta_rel)*pow(gamma, 4e0)*Lattice.Cell[Lattice.param.Cell_nLoc].S);

  D_x = eps[X_]*tau_inv[X_]; D_delta = eps[Z_]*tau_inv[Z_];

  // eps_x*D_x
  D_x *= eps[X_];
  // Compute eps_IBS.
  eps_IBS[X_] = sqrt(D_x*Lattice.param.tau[X_]/2e0);
  // Solve for eps_x
  eps[X_] = eps_SR[X_]*(1e0+sqrt(1e0+4e0*sqr(eps_IBS[X_]/eps_SR[X_])))/2e0;
  // Compute D_x.
  D_x /= eps[X_];

  eps[Y_] = eps_SR[Y_]/eps_SR[X_]*eps[X_];

  D_delta = tau_inv[Z_]*sqr(sigma_delta);
  sigma_delta =
    sqrt((D_delta+2e0/Lattice.param.tau[Z_]*sqr(sigma_delta_SR))
	 *Lattice.param.tau[Z_]/2e0);
  eps[Z_] = sqr(sigma_delta)/gamma_z;

  sigma_s = sqrt(Lattice.param.beta_z*eps[Z_]);

  if (prt2) {
    printf("\nCoulomb Log = %6.3f\n", log_Coulomb);
    printf("D_x         = %9.3e\n", D_x);
    printf("eps_x(IBS)  = %7.3f pm.rad\n", 1e12*eps_IBS[X_]);
    printf("eps_x       = %7.3f pm.rad\n", 1e12*eps[X_]);
    printf("eps_y       = %7.3f pm.rad\n", 1e12*eps[Y_]);
    printf("D_delta     = %9.3e\n", D_delta);
    printf("eps_z       = %9.3e\n", eps[Z_]);
    printf("sigma_s     = %9.3e mm\n", 1e3*sigma_s);
    printf("sigma_delta = %9.3e\n", sigma_delta);
  }
}


void rm_space(char *name)
{
  int  i, k;

  i = 0;
  while (name[i] == ' ')
    i++;

  for (k = i; k <= (int)strlen(name)+1; k++)
    name[k-i] = name[k];
}


void get_bn(const char file_name[], int n, const bool prt)
{
  char      line[max_str], str[max_str], str1[max_str], *token, *name, *p;
  int       n_prm, Fnum, Knum, order;
  double    bnL, bn, C, L;
  FILE      *inf, *fp_lat;

  inf = file_read(file_name); fp_lat = file_write("get_bn.lat");

  no_sxt();

  // if n = 0: go to last data set
  if (n == 0) {
    while (fgets(line, max_str, inf) != NULL )
      if (strstr(line, "n = ") != NULL)	sscanf(line, "n = %d", &n);

    fclose(inf); inf = file_read(file_name);
  }

  if (prt) {
    printf("\n");
    printf("reading values (n=%d): %s\n", n, file_name);
    printf("\n");
  }

  sprintf(str, "n = %d", n);
  do
    fgets(line, max_str, inf);
  while (strstr(line, str) == NULL);

  fprintf(fp_lat, "\n");
  n_prm = 0;
  while (fgets(line, max_str, inf) != NULL) {
    if (strcmp(line, "\n") == 0) break;
    n_prm++;
    name = strtok_r(line, "(", &p);
    rm_space(name);
    strcpy(str, name); Fnum = Lattice.Elem_Index(str);
    strcpy(str1, name); upr_case(str1);
    token = strtok_r(NULL, ")", &p); sscanf(token, "%d", &Knum);
    strtok_r(NULL, "=", &p); token = strtok_r(NULL, "\n", &p);
    sscanf(token, "%lf %d", &bnL, &order);
    if (prt) printf("%6s(%2d) = %10.6f %d\n", name, Knum, bnL, order);
   if (Fnum != 0) {
      if (order == 0)
        SetL(Fnum, bnL);
      else
        SetbnL(Fnum, order, bnL);

      L = GetL(Fnum, 1);
      if (Knum == 1) {
	if (order == 0)
	  fprintf(fp_lat, "%s: Drift, L = %8.6f;\n", str1, bnL);
	else {
	  bn = (L != 0.0)? bnL/L : bnL;
	  if (order == Quad)
	    fprintf(fp_lat, "%s: Quadrupole, L = %8.6f, K = %10.6f, N = Nquad"
		    ", Method = Meth;\n", str1, L, bn);
	  else if (order == Sext)
	    fprintf(fp_lat, "%s: Sextupole, L = %8.6f, K = %10.6f"
		    ", N = Nsext, Method = Meth;\n", str1, L, bn);
	  else {
	    fprintf(fp_lat, "%s: Multipole, L = %8.6f"
		    ", N = 1, Method = Meth,\n", str1, L);
	    fprintf(fp_lat, "     HOM = (%d, %13.6f, %3.1f);\n",
		    order, bn, 0.0);
	  }
	}
      }
    } else {
      printf("element %s not found\n", name);
      exit_(1);
    }
  }
  if (prt) printf("\n");

  C = Lattice.Cell[Lattice.param.Cell_nLoc].S; recalc_S();
  if (prt)
    printf("New Cell Length: %5.3f (%5.3f)\n",
	   Lattice.Cell[Lattice.param.Cell_nLoc].S, C);

  fclose(inf); fclose(fp_lat);
}


double LatticeType::get_dynap(const double delta, const int n_aper,
			       const int n_track, const bool cod)
{
  char      str[max_str];
  int       i;
  double    x_aper[n_aper], y_aper[n_aper], DA;
  FILE      *fp;

  const int  prt = true;

  fp = file_write("dynap.out");
  Lattice.dynap(fp, 5e-3, 0.0, 0.1e-3, n_aper, n_track, x_aper, y_aper, false,
		cod, prt);
  fclose(fp);
  DA = get_aper(n_aper, x_aper, y_aper);

  if (true) {
    sprintf(str, "dynap_dp%3.1f.out", 1e2*delta);
    fp = file_write(str);
    Lattice.dynap(fp, 5e-3, delta, 0.1e-3, n_aper, n_track,
		  x_aper, y_aper, false, cod, prt);
    fclose(fp);
    DA += get_aper(n_aper, x_aper, y_aper);

    for (i = 0; i < nv_; i++)
      Lattice.param.CODvect[i] = 0.0;
    sprintf(str, "dynap_dp%3.1f.out", -1e2*delta);
    fp = file_write(str);
    Lattice.dynap(fp, 5e-3, -delta, 0.1e-3, n_aper,
		  n_track, x_aper, y_aper, false, cod, prt);
    fclose(fp);
    DA += get_aper(n_aper, x_aper, y_aper);
  }

  return DA/3.0;
}


double get_chi2(long int n, double x[], double y[], long int m, psVector b)
{
  /* Compute chi2 for polynomial fit */

  int     i, j;
  double  sum, z;

  sum = 0.0;
  for (i = 0; i < n; i++) {
    z = 0.0;
    for (j = 0; j <= m; j++)
      z += b[j]*pow(x[i], (double)j);
    sum += pow(y[i]-z, 2.0);
  }
  return(sum/n);
}


void pol_fit(int n, double x[], double y[], int order, psVector &b,
	     double &sigma, const bool prt)
{
  /* Polynomial fit by linear chi-square */

  int     i, j, k;
  Matrix  T1;

  const	double sigma_k = 1.0, chi2 = 4.0;

  for (i = 0; i <= order; i++) {
    b[i] = 0.0;
    for (j = 0; j <= order; j++)
      T1[i][j] = 0.0;
  }

  for (i = 0; i < n; i++)
    for (j = 0; j <= order; j++) {
      b[j] += y[i]*pow(x[i], (double)j)/pow(sigma_k, 2.0);
      for (k = 0; k <= order; k++)
	T1[j][k] += pow(x[i], (double)(j+k))/pow(sigma_k, 2.0);
    }

  if (InvMat(order+1, T1)) {
    LinTrans(order+1, T1, b); sigma = get_chi2(n, x, y, order, b);
    if (prt) {
      printf("\n  n    Coeff.\n");
      for (i = 0; i <= order; i++)
	printf("%3d %10.3e+/-%8.2e\n", i, b[i], sigma*sqrt(chi2*T1[i][i]));
    }
  } else {
    printf("pol_fit: Matrix is singular\n");
    exit_(1);
  }
}


void LatticeType::get_ksi2(const double d_delta)
{
  const int     n_points = 20, order = 5;

  int       i, n;
  double    delta[2*n_points+1], nu[2][2*n_points+1], sigma;
  psVector  b;
  FILE      *fp;

  fp = file_write("chrom2.out");
  n = 0;
  for (i = -n_points; i <= n_points; i++) {
    n++; delta[n-1] = i*(double)d_delta/(double)n_points;
    Lattice.Ring_GetTwiss(false, delta[n-1]);
    nu[0][n-1] = Lattice.param.TotalTune[X_]; nu[1][n-1] = Lattice.param.TotalTune[Y_];
    fprintf(fp, "%5.2f %8.5f %8.5f\n", 1e2*delta[n-1], nu[0][n-1], nu[1][n-1]);
  }
  printf("\n");
  printf("Horizontal chromaticity:\n");
  pol_fit(n, delta, nu[0], order, b, sigma, true);
  printf("\n");
  printf("Vertical chromaticity:\n");
  pol_fit(n, delta, nu[1], order, b, sigma, true);
  printf("\n");
  fclose(fp);
}


bool find_nu(const int n, const double nus[], const double eps, double &nu)
{
  bool  lost;
  int   k;

  k = 0;
  while ((k < n) && (fabs(nus[k]-nu) > eps)) {
    if (trace) printf("nu = %7.5f(%7.5f) eps = %7.1e\n", nus[k], nu, eps);
    k++;
  }

  if (k < n) {
    if (trace) printf("nu = %7.5f(%7.5f)\n", nus[k], nu);
    nu = nus[k]; lost = false;
  } else {
    if (trace) printf("lost\n");
    lost = true;
  }

  return !lost;
}


bool get_nu(const double Ax, const double Ay, const double delta,
	    const double eps, double &nu_x, double &nu_y)
{
  const int  n_turn = 512, n_peaks = 5;

  bool      lost, ok_x, ok_y;
  char      str[max_str];
  int       i;
  long int  lastpos, lastn, n;
  double    x[n_turn], px[n_turn], y[n_turn], py[n_turn];
  double    nu[2][n_peaks], A[2][n_peaks];
  psVector    x0;
  FILE      *fp;

  const bool   prt = false;
  const char   file_name[] = "track.out";

  // complex FFT in Floquet space
  x0[x_] = Ax; x0[px_] = 0.0; x0[y_] = Ay; x0[py_] = 0.0;
  LinTrans(4, Lattice.param.Ascrinv, x0);
  Lattice.track(file_name, x0[x_], x0[px_], x0[y_], x0[py_], delta,
		n_turn, lastn, lastpos, 1, 0.0);
  if (lastn == n_turn) {
    GetTrack(file_name, &n, x, px, y, py);
    sin_FFT((int)n, x, px); sin_FFT((int)n, y, py);

    if (prt) {
      strcpy(str, file_name); strcat(str, ".fft");
      fp = file_write(str);
      for (i = 0; i < n; i++)
	fprintf(fp, "%5.3f %12.5e %8.5f %12.5e %8.5f\n",
		(double)i/(double)n, x[i], px[i], y[i], py[i]);
      fclose(fp);
    }

    GetPeaks1(n, x, n_peaks, nu[X_], A[X_]);
    GetPeaks1(n, y, n_peaks, nu[Y_], A[Y_]);

    ok_x = find_nu(n_peaks, nu[X_], eps, nu_x);
    ok_y = find_nu(n_peaks, nu[Y_], eps, nu_y);

    lost = !ok_x || !ok_y;
  } else
    lost = true;

  return !lost;
}


void LatticeType::dnu_dA(const double Ax_max, const double Ay_max,
			  const double delta, const int n_ampl)
{
  bool      ok;
  int       i;
  double    nu_x, nu_y, Ax, Ay, Jx, Jy;
  psVector    ps;
  FILE      *fp;

  const double  A_min  = 0.1e-3;
//  const double  eps0   = 0.04, eps   = 0.02;
//  const double  eps0   = 0.025, eps   = 0.02;
//   const double  eps0   = 0.04, eps   = 0.015;
  const double  eps = 0.01;

  Lattice.Ring_GetTwiss(false, 0.0);

  if (trace) printf("dnu_dAx\n");

  nu_x = fract(Lattice.param.TotalTune[X_]);
  nu_y = fract(Lattice.param.TotalTune[Y_]);

  fp = file_write("dnu_dAx.out");
  fprintf(fp, "#   A_x        A_y        J_x        J_y      nu_x    nu_y\n");
  fprintf(fp, "#\n");
  fprintf(fp, "%10.3e %10.3e %10.3e %10.3e %8.6f %8.6f\n",
	  0e0, 0e0, 0e0, 0e0, fract(nu_x), fract(nu_y));

  Ay = A_min;
  for (i = 1; i <= n_ampl; i++) {
    Ax = -i*Ax_max/n_ampl;
    ps[x_] = Ax; ps[px_] = 0e0; ps[y_] = Ay; ps[py_] = 0e0;
    Lattice.getfloqs(ps);
    Jx = (sqr(ps[x_])+sqr(ps[px_]))/2.0; Jy = (sqr(ps[y_])+sqr(ps[py_]))/2.0;
    ok = get_nu(Ax, Ay, delta, eps, nu_x, nu_y);
    if (ok)
      fprintf(fp, "%10.3e %10.3e %10.3e %10.3e %8.6f %8.6f\n",
	      1e3*Ax, 1e3*Ay, 1e6*Jx, 1e6*Jy, fract(nu_x), fract(nu_y));
    else
      fprintf(fp, "# %10.3e %10.3e particle lost\n", 1e3*Ax, 1e3*Ay);
  }

  if (trace) printf("\n");

  nu_x = fract(Lattice.param.TotalTune[X_]);
  nu_y = fract(Lattice.param.TotalTune[Y_]);

  fprintf(fp, "\n");
  fprintf(fp, "%10.3e %10.3e %10.3e %10.3e %8.6f %8.6f\n",
	  0e0, 0e0, 0e0, 0e0, fract(nu_x), fract(nu_y));

  Ay = A_min;
  for (i = 0; i <= n_ampl; i++) {
    Ax = i*Ax_max/n_ampl;
    ps[x_] = Ax; ps[px_] = 0e0; ps[y_] = Ay; ps[py_] = 0e0;
    Lattice.getfloqs(ps);
    Jx = (sqr(ps[x_])+sqr(ps[px_]))/2.0; Jy = (sqr(ps[y_])+sqr(ps[py_]))/2.0;
    ok = get_nu(Ax, Ay, delta, eps, nu_x, nu_y);
    if (ok)
      fprintf(fp, "%10.3e %10.3e %10.3e %10.3e %8.6f %8.6f\n",
	      1e3*Ax, 1e3*Ay, 1e6*Jx, 1e6*Jy, fract(nu_x), fract(nu_y));
    else
      fprintf(fp, "# %10.3e %10.3e particle lost\n", 1e3*Ax, 1e3*Ay);
  }

  fclose(fp);

  if (trace) printf("dnu_dAy\n");

  nu_x = fract(Lattice.param.TotalTune[X_]);
  nu_y = fract(Lattice.param.TotalTune[Y_]);

  fp = file_write("dnu_dAy.out");
  fprintf(fp, "#   A_x        A_y      nu_x    nu_y\n");
  fprintf(fp, "#\n");
  fprintf(fp, "%10.3e %10.3e %10.3e %10.3e %8.6f %8.6f\n",
	  0e0, 0e0, 0e0, 0e0, fract(nu_x), fract(nu_y));

  Ax = A_min;
  for (i = 1; i <= n_ampl; i++) {
    Ay = -i*Ay_max/n_ampl;
    Jx = pow(Ax, 2.0)/(2.0*Lattice.Cell[Lattice.param.Cell_nLoc].Beta[X_]);
    Jy = pow(Ay, 2.0)/(2.0*Lattice.Cell[Lattice.param.Cell_nLoc].Beta[Y_]);
    ok = get_nu(Ax, Ay, delta, eps, nu_x, nu_y);
    if (ok)
      fprintf(fp, "%10.3e %10.3e %10.3e %10.3e %8.6f %8.6f\n",
	      1e3*Ax, 1e3*Ay, 1e6*Jx, 1e6*Jy, fract(nu_x), fract(nu_y));
    else
      fprintf(fp, "# %10.3e %10.3e particle lost\n", 1e3*Ax, 1e3*Ay);
  }

  if (trace) printf("\n");

  nu_x = fract(Lattice.param.TotalTune[X_]);
  nu_y = fract(Lattice.param.TotalTune[Y_]);

  fprintf(fp, "\n");
  fprintf(fp, "%10.3e %10.3e %10.3e %10.3e %8.6f %8.6f\n",
	  0e0, 0e0, 0e0, 0e0, fract(nu_x), fract(nu_y));

  Ax = A_min;
  for (i = 0; i <= n_ampl; i++) {
    Ay = i*Ay_max/n_ampl;
    Jx = pow(Ax, 2.0)/(2.0*Lattice.Cell[Lattice.param.Cell_nLoc].Beta[X_]);
    Jy = pow(Ay, 2.0)/(2.0*Lattice.Cell[Lattice.param.Cell_nLoc].Beta[Y_]);
    ok = get_nu(Ax, Ay, delta, eps, nu_x, nu_y);
    if (ok)
      fprintf(fp, "%10.3e %10.3e %10.3e %10.3e %8.6f %8.6f\n",
	      1e3*Ax, 1e3*Ay, 1e6*Jx, 1e6*Jy, fract(nu_x), fract(nu_y));
    else
      fprintf(fp, "# %10.3e %10.3e particle lost\n", 1e3*Ax, 1e3*Ay);
  }

  fclose(fp);
}


bool LatticeType::orb_corr(const int n_orbit)
{
  bool      cod = false;
  int       i;
  long      lastpos;
  Vector2   xmean, xsigma, xmax;

  printf("\n");
  Lattice.param.CODvect.zero();
  for (i = 1; i <= n_orbit; i++) {
    cod = Lattice.getcod(0.0, lastpos);
    if (cod) {
      codstat(xmean, xsigma, xmax, Lattice.param.Cell_nLoc, false);
      printf("\n");
      printf("RMS orbit [mm]: (%8.1e+/-%7.1e, %8.1e+/-%7.1e)\n",
	     1e3*xmean[X_], 1e3*xsigma[X_], 1e3*xmean[Y_], 1e3*xsigma[Y_]);
      lsoc(1, 1e0); lsoc(2, 1e0);
      cod = Lattice.getcod(0.0, lastpos);
      if (cod) {
	codstat(xmean, xsigma, xmax, Lattice.param.Cell_nLoc, false);
	printf("RMS orbit [mm]: (%8.1e+/-%7.1e, %8.1e+/-%7.1e)\n",
	       1e3*xmean[X_], 1e3*xsigma[X_], 1e3*xmean[Y_], 1e3*xsigma[Y_]);
      } else
	printf("orb_corr: failed\n");
    } else
      printf("orb_corr: failed\n");
  }

  Lattice.prt_cod("orb_corr.out", Lattice.param.bpm, true);

  return cod;
}


void LatticeType::get_alphac(void)
{
  CellType  Cell;

  Lattice.getelem(Lattice.param.Cell_nLoc, &Cell);
  Lattice.param.Alphac = Lattice.param.OneTurnMat[ct_][delta_]/Cell.S;
}


void LatticeType::get_alphac2(void)
{
  /* Note, do not extract from M[5][4], i.e. around delta
     dependent fixed point.                                */

  const int     n_points = 5;
  const double  d_delta  = 2e-2;

  int       i, j, n;
  long int  lastpos;
  double    delta[2*n_points+1], alphac[2*n_points+1], sigma;
  psVector    x, b;
  CellType  Cell;

  Lattice.param.pathlength = false;
  Lattice.getelem(Lattice.param.Cell_nLoc, &Cell); n = 0;
  for (i = -n_points; i <= n_points; i++) {
    n++; delta[n-1] = i*(double)d_delta/(double)n_points;
    for (j = 0; j < nv_; j++)
      x[j] = 0.0;
    x[delta_] = delta[n-1];
    Lattice.Cell_Pass(0, Lattice.param.Cell_nLoc, x, lastpos);
    alphac[n-1] = x[ct_]/Cell.S;
  }
  pol_fit(n, delta, alphac, 3, b, sigma, true);
  printf("\n");
  printf("alphac = %10.3e + %10.3e*delta + %10.3e*delta^2\n",
	 b[1], b[2], b[3]);
}


double f_bend(double b0L[])
{
  long int lastpos;
  psVector   ps;

  const int   n_prt = 10;

  n_iter_Cart++;

  SetbnL_sys(Fnum_Cart, Dip, b0L[1]);

  ps.zero();
  Lattice.Cell_Pass(Lattice.Elem_GetPos(Fnum_Cart, 1)-1,
		    Lattice.Elem_GetPos(Fnum_Cart, 1), ps, lastpos);

  if (n_iter_Cart % n_prt == 0)
    std::cout << std::scientific << std::setprecision(3)
	      << std::setw(4) << n_iter_Cart
	      << std::setw(11) << ps[x_] << std::setw(11) << ps[px_]
	      << std::setw(11) << ps[ct_]
	      << std::setw(11) << b0L[1] << std::endl;

  return sqr(ps[x_]) + sqr(ps[px_]);
}


void bend_cal_Fam(const int Fnum)
{
  /* Adjusts b1L_sys to zero the orbit for a given gradient. */
  const int  n_prm = 1;

  int    iter;
  double *b0L, **xi, fret;
  psVector ps;

  const double ftol = 1e-15;

  b0L = dvector(1, n_prm); xi = dmatrix(1, n_prm, 1, n_prm);

  std::cout << std::endl;
  std::cout << "bend_cal: " << Lattice.ElemFam[Fnum-1].Name << ":"
	    << std::endl;

  Fnum_Cart = Fnum;  b0L[1] = 0.0; xi[1][1] = 1e-3;

  std::cout << std::endl;
  n_iter_Cart = 0;
  dpowell(b0L, xi, n_prm, ftol, &iter, &fret, f_bend);

  free_dvector(b0L, 1, n_prm); free_dmatrix(xi, 1, n_prm, 1, n_prm);
}


void bend_cal(void)
{
  long int  k;

  for (k = 1; k <= Lattice.param.Elem_nFam; k++)
    if ((Lattice.ElemFam[k-1].Kind == Mpole) &&
	(Lattice.ElemFam[k-1].ElemF.M->irho != 0.0) &&
	(Lattice.ElemFam[k-1].ElemF.M->Bpar[Quad+HOMmax] != 0.0))
      if (Lattice.ElemFam[k-1].nKid > 0) bend_cal_Fam(k);
}


double h_ijklm(const tps &h, const int i, const int j, const int k,
	       const int l, const int m)
{
  int      i1;
  iVector  jj;

  for (i1 = 0; i1 < nv_tps; i1++)
    jj[i1] = 0;
  jj[x_] = i; jj[px_] = j; jj[y_] = k; jj[py_] = l; jj[delta_] = m;
  return h[jj];
}


void set_tune(const char file_name1[], const char file_name2[], const int n)
{
  const int  n_b2 = 8;

  char      line[max_str], names[n_b2][max_str];
  int       j, k, Fnum;
  double    b2s[n_b2], nu[2];
  std::ifstream  inf1, inf2;
  FILE      *fp_lat;

  file_rd(inf1, file_name1); file_rd(inf2, file_name2);

  // skip 1st and 2nd line
  inf1.getline(line, max_str); inf2.getline(line, max_str);
  inf1.getline(line, max_str); inf2.getline(line, max_str);

  inf1.getline(line, max_str);
  sscanf(line, "%*s %*s %*s %s %s %s %s",
	 names[0], names[1], names[2], names[3]);
  inf2.getline(line, max_str);
  sscanf(line, "%*s %*s %*s %s %s %s %s",
	 names[4], names[5], names[6], names[7]);

  printf("\n");
  printf("quads:");
  for (k = 0; k <  n_b2; k++) {
    lwr_case(names[k]); rm_space(names[k]);
    printf(" %s", names[k]);
  }
  printf("\n");

  // skip 4th line
  inf1.getline(line, max_str); inf2.getline(line, max_str);

  fp_lat = file_write("set_tune.lat");
  while (inf1.getline(line, max_str)) {
    sscanf(line, "%d%lf %lf %lf %lf %lf %lf",
	   &j, &nu[X_], &nu[Y_], &b2s[0], &b2s[1], &b2s[2], &b2s[3]);
    inf2.getline(line, max_str);
    sscanf(line, "%d%lf %lf %lf %lf %lf %lf",
	   &j, &nu[X_], &nu[Y_], &b2s[4], &b2s[5], &b2s[6], &b2s[7]);

    if (j == n) {
      printf("\n");
      printf("%3d %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n",
	     n,
	     b2s[0], b2s[1], b2s[2], b2s[3], b2s[4], b2s[5], b2s[6], b2s[7]);

      for (k = 0; k <  n_b2; k++) {
	Fnum = Lattice.Elem_Index(names[k]);
	set_bn_design_fam(Fnum, Quad, b2s[k], 0.0);

	fprintf(fp_lat, "%s: Quadrupole, L = %8.6f, K = %10.6f, N = Nquad"
		", Method = Meth;\n",
		names[k], Lattice.ElemFam[Fnum-1].L, b2s[k]);
      }
      break;
    }
  }

  fclose(fp_lat);
}
