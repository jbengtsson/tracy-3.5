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
  const int max_entries = 20;

  void   *array[max_entries];
  size_t size;
  char   **strings;
  size_t i;

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
  FILE *fp;
  fp = fopen(file_name, "r");
  if (fp == NULL) {
    printf("File not found: %s\n", file_name);
    exit_(-1);
  }
  return(fp);
}


FILE* file_write(const char file_name[])
{
  FILE *fp;
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


void no_sxt(LatticeType &lat)
{
  int       k;
  MpoleType *M;

  // printf("\nzeroing sextupoles\n");
  for (k = 0; k <= lat.conf.Cell_nLoc; k++) {
    M = dynamic_cast<MpoleType*>(lat.elems[k]);
    if ((lat.elems[k]->Pkind == Mpole) && (M->Porder >= Sext))
      SetKpar(lat, lat.elems[k]->Fnum, lat.elems[k]->Knum, Sext, 0.0);
  }
}


void get_map(LatticeType &lat, const bool cod)
{
  long int lastpos;

  map.identity();
  if (cod) {
    lat.getcod(0e0, lastpos);
    map += lat.conf.CODvect;
  }
  lat.Cell_Pass(0, lat.conf.Cell_nLoc, map, lastpos);
  if (cod) map -= lat.conf.CODvect;
}


#if NO > 1

tps get_h(void)
{
  ss_vect<tps> map1, R;

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
  int i, j, k;

  k = 0;
  for (i = 0; i < 2*nd_tps; i++)
    for (j = i; j < 2*nd_tps; j++) {
      k++; m2[k-1] = ps[i]*ps[j];
    }
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


ss_vect<tps> tp_S(const int n_DOF, const ss_vect<tps> &A)
{
  int          j;
  long int     jj[ss_dim];
  ss_vect<tps> S;

  for (j = 1; j <= ss_dim; j++)
    jj[j-1] = (j <= 2*n_DOF)? 1 : 0;

  S = get_S(n_DOF);

  return -S*PInv(A, jj)*S;
}


void get_dnu(const int n, const ss_vect<tps> &A, double dnu[])
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


void get_ab(const ss_vect<tps> &A,
	    double alpha[], double beta[], double dnu[],
	    double eta[], double etap[])
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


ss_vect<tps> get_A(const double alpha[], const double beta[],
		   const double eta[], const double etap[])
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


ss_vect<tps> get_A_CS(const int n, const ss_vect<tps> &A, double dnu[])
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
  int             j, no;
  long int        jj[ss_dim];
  ss_vect<double> z;

  no = no_tps;
  danot_(1);

  for (j = 0; j < nv_tps; j++)
    jj[j] = (j < 2*n_DOF)? 1 : 0;

  z = (PInv(A, jj)*ps).cst();

  for (j = 0; j < n_DOF; j++)
    twoJ[j] = sqr(z[2*j]) + sqr(z[2*j+1]);

  danot_(no);
}


void GetEmittance(LatticeType &lat, const int Fnum, const bool prt)
{
  // A. Chao "Evaluation of Beam Distribution Parameters in an Electron
  // Storage Ring" J. Appl. Phys 50 (2), 595-598.
  bool         emit, rad, cav, path;
  int          i, j, h_RF;
  long int     lastpos, loc;
  double       C, theta, V_RF, phi0, gamma_z;
  double       sigma_s, sigma_delta;
  Vector3      nu;
  Matrix       Ascr;
  ss_vect<tps> Ascr_map;
  CavityType   *Cp;

  // save state
  rad = lat.conf.radiation; emit = lat.conf.emittance;
  cav = lat.conf.Cavity_on; path = lat.conf.pathlength;

  C = lat.elems[lat.conf.Cell_nLoc]->S;

  // damped system
  lat.conf.radiation = true; lat.conf.emittance  = true;
  lat.conf.Cavity_on = true; lat.conf.pathlength = false;

  lat.Ring_GetTwiss(false, 0.0);

  // radiation loss is computed in Cav_Pass

  loc = lat.Elem_GetPos(Fnum, 1);
  Cp = dynamic_cast<CavityType*>(lat.elems[loc]);

  lat.conf.U0 = lat.conf.dE*1e9*lat.conf.Energy;
  V_RF = Cp->Pvolt;
  h_RF = Cp->Ph;
  phi0 = fabs(asin(lat.conf.U0/V_RF));
  lat.conf.delta_RF =
    sqrt(-V_RF*cos(M_PI-phi0)*(2.0-(M_PI-2.0*(M_PI-phi0))
    *tan(M_PI-phi0))/(fabs(lat.conf.Alphac)*M_PI*h_RF*1e9*lat.conf.Energy));

  // Compute diffusion coeffs. for eigenvectors [sigma_xx, sigma_yy, sigma_zz]
  Ascr_map = putlinmat(6, lat.conf.Ascr); Ascr_map += lat.conf.CODvect;

  lat.Cell_Pass(0, lat.conf.Cell_nLoc, Ascr_map, lastpos);

  // K. Robinson "Radiation Effects in Circular Electron Accelerators"
  // Phys. Rev. 111 (2), 373-380.
  // Iu.F. Orlov, E.K. Tarasov "Damping of Oscillations in a Cyclic Electron
  // Accelerator" J. Exptl. Theoret. Phys. 34, 651-657 (1958).
  // To leading order:
  //   Sum_k(alpha_k) = 2*U_0/E
  // or
  //   J_x + J_y + J_z = 4.

  for (i = 0; i < DOF; i++) {
    // partition numbers
    lat.conf.J[i] =
      2.0*(1.0+lat.conf.CODvect[delta_])*lat.conf.alpha_rad[i]/lat.conf.dE;
    // damping times
    lat.conf.tau[i] = -C/(c0*lat.conf.alpha_rad[i]);
    // diffusion coeff. and emittance (alpha is for betatron amplitudes)
    lat.conf.eps[i] = -lat.conf.D_rad[i]/(2.0*lat.conf.alpha_rad[i]);
    // fractional tunes
    nu[i]  = atan2(lat.conf.wi[i*2], lat.conf.wr[i*2])/(2.0*M_PI);
    if (nu[i] < 0.0) nu[i] = 1.0 + nu[i];
  }

  // undamped system
  lat.conf.radiation = !false; lat.conf.emittance = false;

  lat.Ring_GetTwiss(false, 0.0);

  // Compute sigmas arround the lattice:
  //   Sigma = A diag[J_1, J_1, J_2, J_2, J_3, J_3] A^T
  for (i = 0; i < 6; i++) {
    Ascr_map[i] = tps(lat.conf.CODvect[i]);
    for (j = 0; j < 6; j++)
      Ascr_map[i] += lat.conf.Ascr[i][j]*sqrt(lat.conf.eps[j/2])*tps(0.0, j+1);
  }
  // prt_lin_map(3, Ascr_map);
  for (loc = 0; loc <= lat.conf.Cell_nLoc; loc++) {
    lat.elems[loc]->Elem_Pass(lat.conf, Ascr_map);
    // sigma = A x A^tp
    getlinmat(6, Ascr_map, lat.elems[loc]->sigma);
    TpMat(6, lat.elems[loc]->sigma);
    getlinmat(6, Ascr_map, Ascr); MulLMat(6, Ascr, lat.elems[loc]->sigma);
  }

  // A. W. Chao, M. J. Lee "Particle Distribution Parameters in an Electron
  // Storage Ring" J. Appl. Phys. 47 (10), 4453-4456 (1976).
  // observable tilt angle
  theta = atan2(2e0*lat.elems[0]->sigma[x_][y_],
	  (lat.elems[0]->sigma[x_][x_]-lat.elems[0]->sigma[y_][y_]))/2e0;

  // longitudinal alpha and beta
  lat.conf.alpha_z =
    -lat.conf.Ascr[ct_][ct_]*lat.conf.Ascr[delta_][ct_]
    - lat.conf.Ascr[ct_][delta_]*lat.conf.Ascr[delta_][delta_];
  lat.conf.beta_z = sqr(lat.conf.Ascr[ct_][ct_]) + sqr(lat.conf.Ascr[ct_][delta_]);
  gamma_z = (1.0+sqr(lat.conf.alpha_z))/lat.conf.beta_z;

  // bunch size
  sigma_s = sqrt(lat.conf.beta_z*lat.conf.eps[Z_]);
  sigma_delta = sqrt(gamma_z*lat.conf.eps[Z_]);

  if (prt) {
    printf("\nEmittance:\n");
    printf("\nBeam energy [GeV]:              "
	   "Eb          = %4.2f\n", lat.conf.Energy);
    printf("Energy loss per turn [keV]:     "
	   "U0          = %3.1f\n",
	   1e-3*lat.conf.U0);
    printf("Synchronous phase [deg]:        "
	   "phi0        = 180 - %4.2f\n",
	   phi0*180.0/M_PI);
    printf("RF bucket height [%%]:           "
	   "delta_RF    = %4.2f\n", 1e2*lat.conf.delta_RF);
    printf("\n");
    printf("Equilibrium emittance [m.rad]:  "
	   "eps         =  [%9.3e, %9.3e, %9.3e]\n",
            lat.conf.eps[X_], lat.conf.eps[Y_], lat.conf.eps[Z_]);
    printf("Bunch length [mm]:              "
	   "sigma_s     =  %5.3f\n", 1e3*sigma_s);
    printf("Momentum spread:                "
	   "sigma_delta =  %9.3e\n", sigma_delta);
    printf("Longitudinal phase space:       "
	   "alpha_z     = %10.3e, beta_z = %9.3e\n",
	   lat.conf.alpha_z, lat.conf.beta_z);
     printf("Partition numbers:              "
	   "J           =  [%5.3f, %5.3f, %5.3f]\n",
            lat.conf.J[X_], lat.conf.J[Y_], lat.conf.J[Z_]);
    printf("Damping times [msec]:           "
	   "tau         =  [%3.1f, %3.1f, %3.1f]\n",
	   1e3*lat.conf.tau[X_], 1e3*lat.conf.tau[Y_], 1e3*lat.conf.tau[Z_]);
    printf("Diffusion coeffs:               "
	   "D           =  [%7.1e, %7.1e, %7.1e]\n",
	   lat.conf.D_rad[X_], lat.conf.D_rad[Y_], lat.conf.D_rad[Z_]);
    printf("\n");
    printf("alphac:                         "
	   "alphac      = %11.3e\n", lat.conf.Alphac);
    printf("\n");
    printf("Fractional tunes:               "
	   "nu_x        =  [%7.5f, %7.5f, %12.5e]\n", nu[X_], nu[Y_], nu[Z_]);
    printf("                                "
	   "1-nu_x      =  [%7.5f, %7.5f, %12.5e]\n",
	   1e0-nu[X_], 1e0-nu[Y_], 1e0-nu[Z_]);
    printf("\n");
    printf("sigmas:                         "
	   "sigma_x     =  %5.1f  microns, sigma_px    = %5.1f urad\n",
	   1e6*sqrt(lat.elems[0]->sigma[x_][x_]),
	   1e6*sqrt(lat.elems[0]->sigma[px_][px_]));
    printf("                                "
	   "sigma_y     =  %5.1f  microns, sigma_py    = %5.1f urad\n",
	   1e6*sqrt(lat.elems[0]->sigma[y_][y_]),
	   1e6*sqrt(lat.elems[0]->sigma[py_][py_]));
    printf("                                "
	   "sigma_s     =  %6.2f mm,      sigma_delta = %8.2e\n",
	   1e3*sqrt(lat.elems[0]->sigma[ct_][ct_]),
	   sqrt(lat.elems[0]->sigma[delta_][delta_]));

    printf("\n");
    printf("Beam ellipse twist [rad]:       tw      = %5.3f\n", theta);
    printf("                   [deg]:       tw      = %5.3f\n",
	   theta*180.0/M_PI);
  }

  // restore state
  lat.conf.radiation = rad; lat.conf.emittance  = emit;
  lat.conf.Cavity_on = cav; lat.conf.pathlength = path;
}


// output

void printcod(LatticeType &lat, const char *fname)
{
  long     i;
  ElemType cell;
  FILE     *outf;
  long     FORLIM;

  outf = file_write(fname);

  fprintf(outf,
	  "#       name             s    betax   nux   betay   nuy      posx"
	  "          posy          dSx          dSy         dipx         dipy"
	  "      posx-dSx     posy-dSy\n");
  fprintf(outf,
	  "#                       [m]    [m]           [m]             [m]"
	  "           [m]           [m]          [m]         [rad]        [rad]"
	  "       [um]         [um]\n#\n");

  FORLIM = lat.conf.Cell_nLoc;
  for (i = 1; i <= FORLIM; i++) {
    lat.getelem(i, &cell);

    /* COD is in local coordinates */
    fprintf(outf, "%4ld:%15s ", i, cell.PName);
    fprintf(outf, "%6.2f%7.3f%7.3f%7.3f%7.3f % .5E % .5E % .5E % .5E % .5E"
	    " % .5E % .5E % .5E\n",
	    cell.S, cell.Beta[X_], cell.Nu[X_], cell.Beta[Y_], cell.Nu[Y_],
	    cell.BeamPos[x_], cell.BeamPos[y_], cell.dS[X_], cell.dS[Y_],
	    -lat.Elem_GetKval(cell.Fnum, cell.Knum, (long)Dip),
	     lat.Elem_GetKval(cell.Fnum, cell.Knum, (long)(-Dip)),
	    (cell.BeamPos[x_]-cell.dS[X_])*1.e6,
	    (cell.BeamPos[y_]-cell.dS[Y_])*1.e6);
  }
  if (outf != NULL)
    fclose(outf);
}


void prt_cod(LatticeType &lat, const char *file_name, const bool all)
{
  long      i;
  FILE      *outf;
  long      FORLIM;
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

  FORLIM = lat.conf.Cell_nLoc;
  for (i = 0L; i <= FORLIM; i++) {
    if (all || (lat.elems[i]->Fnum == lat.conf.bpm)) {
      /* COD is in local coordinates */
      fprintf(outf,
	      "%4ld %.*s %6.2f %4.1f %6.3f %6.3f %6.3f %6.3f"
	      " %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",
	      i, SymbolLength, lat.elems[i]->PName, lat.elems[i]->S,
	      get_code(lat.conf, *lat.elems[i]),
	      lat.elems[i]->Beta[X_], lat.elems[i]->Nu[X_],
	      lat.elems[i]->Beta[Y_], lat.elems[i]->Nu[Y_],
	      1e3*lat.elems[i]->BeamPos[x_], 1e3*lat.elems[i]->BeamPos[y_],
	      1e3*lat.elems[i]->dS[X_], 1e3*lat.elems[i]->dS[Y_],
	      -1e3*lat.Elem_GetKval(lat.elems[i]->Fnum,
				    lat.elems[i]->Knum, Dip),
	      1e3*lat.Elem_GetKval(lat.elems[i]->Fnum, lat.elems[i]->Knum,
				   -Dip));
    }
  }
  fclose(outf);
}


void prt_beampos(LatticeType &lat, const char *file_name)
{
  long int k;
  FILE     *outf;

  outf = file_write(file_name);

  fprintf(outf, "#       name             s  code    xcod   ycod\n");
  fprintf(outf, "#                       [m]          [m]     [m]\n");
  fprintf(outf, "#\n");

  for (k = 0; k <= lat.conf.Cell_nLoc; k++)
    fprintf(outf, "%4ld %.*s %6.2f %4.1f %12.5e %12.5e\n",
	    k, SymbolLength, lat.elems[k]->PName, lat.elems[k]->S,
	    get_code(lat.conf, *lat.elems[k]), lat.elems[k]->BeamPos[x_],
	    lat.elems[k]->BeamPos[y_]);

  fclose(outf);
}


// misalignments

void CheckAlignTol(LatticeType &lat, const char *OutputFile)
  // check aligment errors of individual magnets on giders
  // the dT and roll angle are all printed out
{
  int          i, j;
  int          n_girders;
  int          gs_Fnum, ge_Fnum;
  int          gs_nKid, ge_nKid;
  int          dip_Fnum,dip_nKid;
  int          loc, loc_gs, loc_ge;
  char         *name;
  double       s;
  double       PdSsys[2], PdSrms[2], PdSrnd[2], dS[2], dT[2];
  MpoleType    *M;
  std::fstream fout;

  gs_Fnum = lat.conf.gs;   gs_nKid = lat.GetnKid(gs_Fnum);
  ge_Fnum = lat.conf.ge;   ge_nKid = lat.GetnKid(ge_Fnum);
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
    loc_gs = lat.Elem_GetPos(gs_Fnum, i); loc_ge = lat.Elem_GetPos(ge_Fnum, i);

    loc = loc_gs;
    M = dynamic_cast<MpoleType*>(lat.elems[loc]);
    PdSsys[X_] = M->PdSsys[X_];
    PdSsys[Y_] = M->PdSsys[Y_];
    PdSrms[X_] = M->PdSrms[X_];
    PdSrms[Y_] = M->PdSrms[Y_];
    PdSrnd[X_] = M->PdSrnd[X_];
    PdSrnd[Y_] = M->PdSrnd[Y_];
    dS[X_] = lat.elems[loc]->dS[X_]; dS[Y_] = lat.elems[loc]->dS[Y_];
    dT[0] = lat.elems[loc]->dT[0]; dT[1] = lat.elems[loc]->dT[1];
    s = lat.elems[loc]->S; name = lat.elems[loc]->PName;
    fout << "  " << name << "  " << loc << "   " << s
	 << "  " <<  PdSsys[X_] << "  " <<  PdSsys[Y_]
	 << "   " << PdSrms[X_] << "  " <<  PdSrms[Y_]
	 << "   " << PdSrnd[X_] << "  " <<  PdSrnd[Y_]
         << "   " << M->PdTrms << "  "
	 << M->PdTrnd << "   " << dS[X_]     << "  " <<  dS[Y_]
	 << "   " << atan2( dT[1], dT[0] )  << std::endl;

    for (j = loc_gs+1; j < loc_ge; j++) {
      loc = j;
      M = dynamic_cast<MpoleType*>(lat.elems[loc]);
      if ((lat.elems[j]->Pkind == Mpole) &&
	  (M->n_design >= Quad || M->n_design >= Sext)) {
	PdSsys[X_] = M->PdSsys[X_];
	PdSsys[Y_] = M->PdSsys[Y_];
	PdSrms[X_] = M->PdSrms[X_];
	PdSrms[Y_] = M->PdSrms[Y_];
	PdSrnd[X_] = M->PdSrnd[X_];
	PdSrnd[Y_] = M->PdSrnd[Y_];
	dS[X_] = lat.elems[loc]->dS[X_]; dS[Y_] = lat.elems[loc]->dS[Y_];
	dT[0] = lat.elems[loc]->dT[0];   dT[1] = lat.elems[loc]->dT[1];
	s = lat.elems[loc]->S; name=lat.elems[loc]->PName;
	fout << "  " << name << "  " << loc << "   " << s
	     << "  " <<  PdSsys[X_] << "  " <<  PdSsys[Y_]
	     << "   " << PdSrms[X_] << "  " <<  PdSrms[Y_]
	     << "   " << PdSrnd[X_] << "  " <<  PdSrnd[Y_]
	     << "   " << M->PdTrms << "  "
	     << M->PdTrnd
	     << "   " << dS[X_] << "  " <<  dS[Y_]
	     << "   " << atan2( dT[1], dT[0] )  << std::endl;
      }
    }

    loc = loc_ge;
    M = dynamic_cast<MpoleType*>(lat.elems[loc]);
    PdSsys[X_] = M->PdSsys[X_];
    PdSsys[Y_] = M->PdSsys[Y_];
    PdSrms[X_] = M->PdSrms[X_];
    PdSrms[Y_] = M->PdSrms[Y_];
    PdSrnd[X_] = M->PdSrnd[X_];
    PdSrnd[Y_] = M->PdSrnd[Y_];
    dS[X_] = lat.elems[loc]->dS[X_]; dS[Y_] = lat.elems[loc]->dS[Y_];
    dT[0] = lat.elems[loc]->dT[0]; dT[1] = lat.elems[loc]->dT[1];
    s=lat.elems[loc]->S; name=lat.elems[loc]->PName;
    fout << "  " << name << "  " << loc << "   " << s
	 << "  " <<  PdSsys[X_] << "  " <<  PdSsys[Y_]
	 << "   " << PdSrms[X_] << "  " <<  PdSrms[Y_]
	 << "   " << PdSrnd[X_] << "  " <<  PdSrnd[Y_]
         << "   " << M->PdTrms
	 << "  " << M->PdTrnd
         << "   " << dS[X_]     << "  " <<  dS[Y_]
	 << "   " << atan2( dT[1], dT[0] )  << std::endl;

  }

  fout << "  " << std::endl;
  fout << "Dipoles:  " << std::endl;
  dip_Fnum = ElemIndex("B1"); dip_nKid = lat.GetnKid(dip_Fnum);
  for (i = 1; i <= dip_nKid; i++){
    loc = lat.Elem_GetPos(dip_Fnum, i);
    M = dynamic_cast<MpoleType*>(lat.elems[loc]);
    PdSsys[X_] = M->PdSsys[X_];
    PdSsys[Y_] = M->PdSsys[Y_];
    PdSrms[X_] = M->PdSrms[X_];
    PdSrms[Y_] = M->PdSrms[Y_];
    PdSrnd[X_] = M->PdSrnd[X_];
    PdSrnd[Y_] = M->PdSrnd[Y_];
    dS[X_] = lat.elems[loc]->dS[X_]; dS[Y_] = lat.elems[loc]->dS[Y_];
    dT[0] = lat.elems[loc]->dT[0]; dT[1] = lat.elems[loc]->dT[1];
    s = lat.elems[loc]->S; name = lat.elems[loc]->PName;
    fout << "  " << name << "  " << loc << "   " << s
	 << "  " <<  PdSsys[X_] << "  " <<  PdSsys[Y_]
	 << "   " << PdSrms[X_] << "  " <<  PdSrms[Y_]
	 << "   " << PdSrnd[X_] << "  " <<  PdSrnd[Y_]
	 << "   " << M->PdTrms
	 << "  " << M->PdTrnd
	 << "   " << dS[X_]     << "  " <<  dS[Y_]
	 << "   " << atan2( dT[1], dT[0] )  << std::endl;
  }

  fout.close();
}


void misalign_rms_elem(LatticeType &lat, const int Fnum, const int Knum,
		       const double dx_rms, const double dy_rms,
		       const double dr_rms, const bool new_rnd)
{
  long int  loc;
  MpoleType *mp;

  loc = lat.Elem_GetPos(Fnum, Knum);
  mp = dynamic_cast<MpoleType*>(lat.elems[loc]);

  mp->PdSrms[X_] = dx_rms; mp->PdSrms[Y_] = dy_rms; mp->PdTrms = dr_rms;
  if (new_rnd) {
    if (normal) {
      mp->PdSrnd[X_] = normranf(); mp->PdSrnd[Y_] = normranf();
      mp->PdTrnd = normranf();
    } else {
      mp->PdSrnd[X_] = ranf(); mp->PdSrnd[Y_] = ranf();
      mp->PdTrnd = ranf();
    }
  }

  lat.SetdS(Fnum, Knum); lat.SetdT(Fnum, Knum);
}

void misalign_sys_elem(LatticeType &lat, const int Fnum, const int Knum,
		       const double dx_sys, const double dy_sys,
		       const double dr_sys)
{
  long int  loc;
  MpoleType *mp;

  loc = lat.Elem_GetPos(Fnum, Knum);
  mp  = dynamic_cast<MpoleType*>(lat.elems[loc]);

  mp->PdSsys[X_] = dx_sys; mp->PdSsys[Y_] = dy_sys; mp->PdTsys = dr_sys;

  lat.SetdS(Fnum, Knum); lat.SetdT(Fnum, Knum);
}

void misalign_rms_fam(LatticeType &lat, const int Fnum,
		      const double dx_rms, const double dy_rms,
		      const double dr_rms, const bool new_rnd)
{
  int  i;

  for (i = 1; i <= lat.GetnKid(Fnum); i++)
    misalign_rms_elem(lat, Fnum, i, dx_rms, dy_rms, dr_rms, new_rnd);
}

void misalign_sys_fam(LatticeType &lat, const int Fnum,
		      const double dx_sys, const double dy_sys,
		      const double dr_sys)
{
  int  i;

  for (i = 1; i <= lat.GetnKid(Fnum); i++)
    misalign_sys_elem(lat, Fnum, i, dx_sys, dy_sys, dr_sys);
}

void misalign_rms_type(LatticeType &lat, const int type,
		       const double dx_rms, const double dy_rms,
		       const double dr_rms, const bool new_rnd)
{
  long int  k;
  MpoleType *M;

  if ((type >= All) && (type <= HOMmax)) {
    for (k = 1; k <= lat.conf.Cell_nLoc; k++) {
      M = dynamic_cast<MpoleType*>(lat.elems[k]);
      if ((lat.elems[k]->Pkind == Mpole) &&
	  ((type == M->n_design) || ((type == All) &&
	   ((lat.elems[k]->Fnum != lat.conf.gs)
	    && (lat.elems[k]->Fnum != lat.conf.ge))))) {
	// if all: skip girders
	misalign_rms_elem(lat, lat.elems[k]->Fnum, lat.elems[k]->Knum,
			  dx_rms, dy_rms, dr_rms, new_rnd);
      }
    }
  } else {
    printf("misalign_rms_type: incorrect type %d\n", type); exit_(1);
  }
}

void misalign_sys_type(LatticeType &lat, const int type,
		       const double dx_sys, const double dy_sys,
		       const double dr_sys)
{
  long int  k;
  MpoleType *M;

  if ((type >= All) && (type <= HOMmax)) {
    for (k = 1; k <= lat.conf.Cell_nLoc; k++) {
      M = dynamic_cast<MpoleType*>(lat.elems[k]);
      if ((lat.elems[k]->Pkind == Mpole) &&
	  ((type == M->n_design) || ((type == All) &&
	   ((lat.elems[k]->Fnum != lat.conf.gs)
	    && (lat.elems[k]->Fnum != lat.conf.ge))))) {
	// if all: skip girders
	misalign_sys_elem(lat, lat.elems[k]->Fnum, lat.elems[k]->Knum,
			  dx_sys, dy_sys, dr_sys);
      }
    }
  } else {
    printf("misalign_sys_type: incorrect type %d\n", type); exit_(1);
  }
}

void misalign_rms_girders(LatticeType &lat, const int gs, const int ge,
			  const double dx_rms, const double dy_rms,
			  const double dr_rms, const bool new_rnd)
{
  int       i, k, n_girders, n_ge, n_gs;
  long int  loc_gs, loc_ge, j;
  double    s_gs, s_ge, dx_gs[2], dx_ge[2], s;
  MpoleType *Mgs, *Mge, *Mj;

  n_gs = lat.GetnKid(gs); n_ge = lat.GetnKid(ge);

  if (n_gs == n_ge)
    n_girders = n_gs;
  else {
    std::cout << "set_girders: no of GS != no of GE" << std::endl;
    exit (1);
  }

  misalign_rms_fam(lat, gs, dx_rms, dy_rms, dr_rms, new_rnd);
  misalign_rms_fam(lat, ge, dx_rms, dy_rms, dr_rms, new_rnd);

  for (i = 1; i <= n_girders; i++) {
    loc_gs = lat.Elem_GetPos(gs, i);
    loc_ge = lat.Elem_GetPos(ge, i);
    s_gs = lat.elems[loc_gs]->S;
    s_ge = lat.elems[loc_ge]->S;

    // roll for a rigid boby
    // Note, girders needs to be introduced as gs->ge pairs
    Mgs = dynamic_cast<MpoleType*>(lat.elems[loc_gs]);
    Mge = dynamic_cast<MpoleType*>(lat.elems[loc_ge]);
    Mge->PdTrnd = Mgs->PdTrnd;
    lat.SetdT(ge, i);

    for (k = 0; k <= 1; k++) {
      dx_gs[k] = lat.elems[loc_gs]->dS[k]; dx_ge[k] = lat.elems[loc_ge]->dS[k];
    }

    // move elements onto mis-aligned girder
    for (j = loc_gs+1; j < loc_ge; j++) {
      if ((lat.elems[j]->Pkind == Mpole)
	  || (lat.elems[j]->Fnum == lat.conf.bpm)) {
	Mj = dynamic_cast<MpoleType*>(lat.elems[j]);
        s = lat.elems[j]->S;
	for (k = 0; k <= 1; k++)
	  Mj->PdSsys[k] = dx_gs[k] + (dx_ge[k]-dx_gs[k])*(s-s_gs)/(s_ge-s_gs);
	Mj->PdTsys = Mgs->PdTrms*Mgs->PdTrnd;
      }
    }
  }
}


void misalign_sys_girders(LatticeType &lat, const int gs, const int ge,
			  const double dx_sys, const double dy_sys,
			  const double dr_sys)
{
  int       i, k, n_girders, n_ge, n_gs;
  long int  loc_gs, loc_ge, j;
  double    s_gs, s_ge, dx_gs[2], dx_ge[2], s;
  MpoleType *Mgs, *Mge, *Mj;

  n_gs = lat.GetnKid(gs); n_ge = lat.GetnKid(ge);

  if (n_gs == n_ge)
    n_girders = n_gs;
  else {
    std::cout << "set_girders: no of GS != no of GE" << std::endl;
    exit (1);
  }

  misalign_sys_fam(lat, gs, dx_sys, dy_sys, dr_sys);
  misalign_sys_fam(lat, ge, dx_sys, dy_sys, dr_sys);

  for (i = 1; i <= n_girders; i++) {
    loc_gs = lat.Elem_GetPos(gs, i); loc_ge = lat.Elem_GetPos(ge, i);
    s_gs = lat.elems[loc_gs]->S; s_ge = lat.elems[loc_ge]->S;

    // roll for a rigid boby
    // Note, girders needs to be introduced as gs->ge pairs
    Mgs = dynamic_cast<MpoleType*>(lat.elems[loc_gs]);
    Mge = dynamic_cast<MpoleType*>(lat.elems[loc_ge]);
    Mge->PdTrnd = Mgs->PdTrnd;
    lat.SetdT(ge, i);

    for (k = 0; k <= 1; k++) {
      dx_gs[k] = lat.elems[loc_gs]->dS[k]; dx_ge[k] = lat.elems[loc_ge]->dS[k];
    }

    // move elements onto mis-aligned girder
    for (j = loc_gs+1; j < loc_ge; j++) {
      if ((lat.elems[j]->Pkind == Mpole)
	  || (lat.elems[j]->Fnum == lat.conf.bpm)) {
	Mj = dynamic_cast<MpoleType*>(lat.elems[j]);
        s = lat.elems[j]->S;
	for (k = 0; k <= 1; k++)
	  Mj->PdSsys[k] = dx_gs[k] + (dx_ge[k]-dx_gs[k])*(s-s_gs)/(s_ge-s_gs);
	Mj->PdTsys = Mgs->PdTrms*Mgs->PdTrnd;
      }
    }
  }
}


// apertures

void set_aper_elem(LatticeType &lat, const int Fnum, const int Knum,
		   const double Dxmin, const double Dxmax,
		   const double Dymin, const double Dymax)
{
  int k;

    k = lat.Elem_GetPos(Fnum, Knum);
    lat.elems[k]->maxampl[X_][0] = Dxmin; lat.elems[k]->maxampl[X_][1] = Dxmax;
    lat.elems[k]->maxampl[Y_][0] = Dymin; lat.elems[k]->maxampl[Y_][1] = Dymax;
 }

void set_aper_fam(LatticeType &lat, const int Fnum,
		  const double Dxmin, const double Dxmax,
		  const double Dymin, const double Dymax)
{
  int k;

  for (k = 1; k <= lat.GetnKid(Fnum); k++)
    set_aper_elem(lat, Fnum, k, Dxmin, Dxmax, Dymin, Dymax);
}

void set_aper_type(LatticeType &lat, const int type, const double Dxmin,
		   const double Dxmax, const double Dymin, const double Dymax)
{
  long int  k;
  MpoleType *M;

  if (type >= All && type <= HOMmax) {
    for(k = 1; k <= lat.conf.Cell_nLoc; k++) {
      M = dynamic_cast<MpoleType*>(lat.elems[k]);
      if (((lat.elems[k]->Pkind == Mpole) && (M->n_design == type))
	  || (type == All))
	set_aper_elem(lat, lat.elems[k]->Fnum, lat.elems[k]->Knum, Dxmin, Dxmax,
		      Dymin, Dymax);
    }
  } else
    printf("set_aper_type: bad design type %d\n", type);
}


double get_L(LatticeType &lat, const int Fnum, const int Knum)
{
  return lat.elems[lat.Elem_GetPos(Fnum, Knum)]->PL;
}


void set_L(LatticeType &lat, const int Fnum, const int Knum, const double L)
{
  long int  loc;
  double    phi;
  ElemType  *elemp;
  MpoleType *M;

  loc = lat.Elem_GetPos(Fnum, Knum);
  elemp = lat.elems[loc];
  if (elemp->Pkind == Mpole) {
    M = dynamic_cast<MpoleType*>(elemp);
    if (M->Pirho != 0e0) {
      // Phi is constant.
      phi = elemp->PL*M->Pirho; M->Pirho = phi/L;
      // M->Pc0 = sin(phi/2e0);
    }
  }
  elemp->PL = L;
}


void set_L(LatticeType &lat, const int Fnum, const double L)
{
  int k;

  for (k = 1; k <= lat.GetnKid(Fnum); k++)
    set_L(lat, Fnum, k, L);
}


void set_dL(LatticeType &lat, const int Fnum, const int Knum, const double dL)
{
  lat.elems[lat.Elem_GetPos(Fnum, Knum)]->PL += dL;
}


void set_dL(LatticeType &lat, const int Fnum, const double dL)
{
  int k;

  for (k = 1; k <= lat.GetnKid(Fnum); k++)
    set_dL(lat, Fnum, k, dL);
}


// Multipole Components.

void LatticeType::set_b_n(const set_mpole set_, const int Fnum, const int Knum,
			 const int n, const double b_n)
{
  //       / b_n, n > 0 
  // : set |
  //       \ a_n, n < 0

  ElemType  *elemp;
  MpoleType *M;

  if ((-HOMmax <= n) && (n <= HOMmax)) {
    elemp = elems[Elem_GetPos(Fnum, Knum)];
    M = dynamic_cast<MpoleType*>(elemp);

    switch (set_) {
    case set_mpole(b_n_):
      M->PBpar[HOMmax+n] = b_n;
      break;
    case set_mpole(db_n_):
      M->PBpar[HOMmax+n] += b_n;
      break;
    case set_mpole(b_nL_):
      M->PBpar[HOMmax+n] = (elemp->PL != 0e0)? b_n/elemp->PL : b_n;
      break;
    case set_mpole(db_nL_):
      M->PBpar[HOMmax+n] += (elemp->PL != 0e0)? b_n/elemp->PL : b_n;
      break;
    default:
      printf("\nset_b_n: unknown set_ %d\n", set_);
      exit(1);
      break;
    }

   SetPB(Fnum, Knum, n);
  } else {
    printf("set_b_n: n < 1 (%d)\n", n);
    exit(1);
  }
}

//------------------------------------------------------------------------------

void get_bn_design_elem(LatticeType &lat, const int Fnum, const int Knum,
			const int n, double &bn, double &an)
{
  ElemType  *elem;
  MpoleType *M;

  if (n < 1) {
    std::cout << "get_bn_design_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  elem = lat.elems[lat.Elem_GetPos(Fnum, Knum)];
  M = dynamic_cast<MpoleType*>(elem);

  bn = M->PBpar[HOMmax+n]; an = M->PBpar[HOMmax-n];
}


void get_bnL_design_elem(LatticeType &lat, const int Fnum, const int Knum,
			 const int n, double &bnL, double &anL)
{
  ElemType  *elem;
  MpoleType *M;

  if (n < 1) {
    std::cout << "get_bnL_design_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  elem = lat.elems[lat.Elem_GetPos(Fnum, Knum)];
  M = dynamic_cast<MpoleType*>(elem);

  bnL = M->PBpar[HOMmax+n]; anL = M->PBpar[HOMmax-n];

  if (elem->PL != 0.0) {
    bnL *= elem->PL; anL *= elem->PL;
  }
}


void set_bn_design_elem(LatticeType &lat, const int Fnum, const int Knum,
			const int n, const double bn, const double an)
{
  ElemType  *elem;
  MpoleType *M;

  if (n < 1) {
    std::cout << "set_bn_design_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  elem = lat.elems[lat.Elem_GetPos(Fnum, Knum)];
  M = dynamic_cast<MpoleType*>(elem);

  M->PBpar[HOMmax+n] = bn; M->PBpar[HOMmax-n] = an;

  lat.SetPB(Fnum, Knum, n); lat.SetPB(Fnum, Knum, -n);
}


void set_bn_design_fam(LatticeType &lat, const int Fnum, const int n,
		       const double bn, const double an)
{
  int k;

  if (n < 1) {
    std::cout << "set_bn_design_fam: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  for (k = 1; k <= lat.GetnKid(Fnum); k++)
    set_bn_design_elem(lat, Fnum, k, n, bn, an);
}


void set_dbn_design_elem(LatticeType &lat, const int Fnum, const int Knum,
			 const int n, const double dbn, const double dan)
{
  ElemType  *elem;
  MpoleType *M;

  if (n < 1) {
    std::cout << "set_dbn_design_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  elem = lat.elems[lat.Elem_GetPos(Fnum, Knum)];
  M = dynamic_cast<MpoleType*>(elem);

  M->PBpar[HOMmax+n] += dbn; M->PBpar[HOMmax-n] += dan;

  lat.SetPB(Fnum, Knum, n); lat.SetPB(Fnum, Knum, -n);
}


void set_dbn_design_fam(LatticeType &lat, const int Fnum, const int n,
			const double dbn, const double dan)
{
  int k;

  if (n < 1) {
    std::cout << "set_dbn_design_fam: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  for (k = 1; k <= lat.GetnKid(Fnum); k++)
    set_dbn_design_elem(lat, Fnum, k, n, dbn, dan);
}


void set_bnL_design_elem(LatticeType &lat, const int Fnum, const int Knum,
			 const int n, const double bnL, const double anL)
{
  ElemType  *elem;
  MpoleType *M;

  if (n < 1) {
    std::cout << "set_bnL_design_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  elem = lat.elems[lat.Elem_GetPos(Fnum, Knum)];
  M = dynamic_cast<MpoleType*>(elem);

  if (elem->PL != 0.0) {
    M->PBpar[HOMmax+n] = bnL/elem->PL;
    M->PBpar[HOMmax-n] = anL/elem->PL;
  } else {
    // thin kick
    M->PBpar[HOMmax+n] = bnL; M->PBpar[HOMmax-n] = anL;
  }

  lat.SetPB(Fnum, Knum, n); lat.SetPB(Fnum, Knum, -n);
}


void set_dbnL_design_elem(LatticeType &lat, const int Fnum, const int Knum,
			  const int n, const double dbnL, const double danL)
{
  ElemType  *elem;
  MpoleType *M;

  if (n < 1) {
    std::cout << "set_dbnL_design_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  elem = lat.elems[lat.Elem_GetPos(Fnum, Knum)];
  M = dynamic_cast<MpoleType*>(elem);

  if (elem->PL != 0.0) {
    M->PBpar[HOMmax+n] += dbnL/elem->PL;
    M->PBpar[HOMmax-n] += danL/elem->PL;
  } else {
    // thin kick
    M->PBpar[HOMmax+n] += dbnL; M->PBpar[HOMmax-n] += danL;
  }

  lat.SetPB(Fnum, Knum, n); lat.SetPB(Fnum, Knum, -n);
}


void set_dbnL_design_fam(LatticeType &lat, const int Fnum, const int n,
			 const double dbnL, const double danL)
{
  int k;

  if (n < 1) {
    std::cout << "set_dbnL_design_fam: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  for (k = 1; k <= lat.GetnKid(Fnum); k++)
    set_dbnL_design_elem(lat, Fnum, k, n, dbnL, danL);
}


void set_bnL_design_fam(LatticeType &lat, const int Fnum, const int n,
			const double bnL, const double anL)
{
  int k;

  if (n < 1) {
    std::cout << "set_bnL_design_fam: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  for (k = 1; k <= lat.GetnKid(Fnum); k++)
    set_bnL_design_elem(lat, Fnum, k, n, bnL, anL);
}

//------------------------------------------------------------------------------

void set_bnL_design_type(LatticeType &lat, const int type, const int n,
			 const double bnL, const double anL)
{
  long int  k;
  MpoleType *M;

  if (n < 1) {
    std::cout << "set_bnL_design_type: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  if ((type >= Dip) && (type <= HOMmax)) {
    for (k = 1; k <= lat.conf.Cell_nLoc; k++) {
      M = dynamic_cast<MpoleType*>(lat.elems[k]);
      if ((lat.elems[k]->Pkind == Mpole) && (M->n_design == type))
	set_bnL_design_elem(lat, lat.elems[k]->Fnum, lat.elems[k]->Knum, n, bnL,
			    anL);
    }
  } else
    printf("Bad type argument to set_bnL_design_type()\n");
}

//------------------------------------------------------------------------------

void set_bnL_sys_elem(LatticeType &lat, const int Fnum, const int Knum,
		      const int n, const double bnL, const double anL)
{
  ElemType  *elem;
  MpoleType *M;

  if (n < 1) {
    std::cout << "set_bnL_sys_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  elem = lat.elems[lat.Elem_GetPos(Fnum, Knum)];
  M = dynamic_cast<MpoleType*>(elem);

  if (elem->PL != 0.0) {
    M->PBsys[HOMmax+n] = bnL/elem->PL;
    M->PBsys[HOMmax-n] = anL/elem->PL;
  } else {
    // thin kick
    M->PBsys[HOMmax+n] = bnL; M->PBsys[HOMmax-n] = anL;
  }

  lat.SetPB(Fnum, Knum, n); lat.SetPB(Fnum, Knum, -n);

  if (trace) {
    printf("set_bnL_sys_elem: %s %3d %e %e\n",
	   elem->PName, n, M->PBsys[HOMmax+n], M->PBsys[HOMmax-n]);
    printf("                                      %e %e\n",
	   M->PB[HOMmax+n], M->PB[HOMmax-n]);
  }
}


void set_bnL_sys_fam(LatticeType &lat, const int Fnum, const int n,
		     const double bnL, const double anL)
{
  int k;

  if (n < 1) {
    std::cout << "set_bnL_sys_fam: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  for (k = 1; k <= lat.GetnKid(Fnum); k++)
    set_bnL_sys_elem(lat, Fnum, k, n, bnL, anL);
}


void set_bnL_sys_type(LatticeType &lat, const int type, const int n,
		      const double bnL, const double anL)
{
  long int  k;
  MpoleType *M;

  if (n < 1) {
    std::cout << "set_bnL_sys_type: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  if (type >= Dip && type <= HOMmax) {
    for(k = 1; k <= lat.conf.Cell_nLoc; k++) {
      M = dynamic_cast<MpoleType*>(lat.elems[k]);
      if ((lat.elems[k]->Pkind == Mpole) && (M->n_design == type))
	set_bnL_sys_elem(lat, lat.elems[k]->Fnum, lat.elems[k]->Knum, n, bnL,
			 anL);
    }
  } else
    printf("Bad type argument to set_bnL_sys_type()\n");
}


void set_bnL_rms_elem(LatticeType &lat, const int Fnum, const int Knum,
		      const int n, const double bnL, const double anL,
		      const bool new_rnd)
{
  ElemType  *elem;
  MpoleType *M;

  bool prt = false;

  if (n < 1) {
    std::cout << "set_bnL_rms_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  elem = lat.elems[lat.Elem_GetPos(Fnum, Knum)];
  M = dynamic_cast<MpoleType*>(elem);

  if (elem->PL != 0.0) {
    M->PBrms[HOMmax+n] = bnL/elem->PL;
    M->PBrms[HOMmax-n] = anL/elem->PL;
  } else {
    // thin kick
    M->PBrms[HOMmax+n] = bnL; M->PBrms[HOMmax-n] = anL;
  }

  if(new_rnd){
    if (normal) {
      M->PBrnd[HOMmax+n] = normranf();
      M->PBrnd[HOMmax-n] = normranf();
    } else {
      M->PBrnd[HOMmax+n] = ranf(); M->PBrnd[HOMmax-n] = ranf();
    }
  }

  if (prt)
    printf("set_bnL_rms_elem:  Fnum = %d, Knum = %d"
	   ", bnL = %e, anL = %e %e %e\n",
	   Fnum, Knum, bnL, anL,
	   M->PBrms[HOMmax+n], M->PBrms[HOMmax-n]);

  lat.SetPB(Fnum, Knum, n); lat.SetPB(Fnum, Knum, -n);
}


void set_bnL_rms_fam(LatticeType &lat, const int Fnum, const int n,
		     const double bnL, const double anL, const bool new_rnd)
{
  int k;

  if (n < 1) {
    std::cout << "set_bnL_rms_fam: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  for (k = 1; k <= lat.GetnKid(Fnum); k++)
    set_bnL_rms_elem(lat, Fnum, k, n, bnL, anL, new_rnd);
}


void set_bnL_rms_type(LatticeType &lat, const int type, const int n,
		      const double bnL, const double anL, const bool new_rnd)
{
  long int  k;
  MpoleType *M;

  if (n < 1) {
    std::cout << "get_bnL_rms_type: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  if (type >= Dip && type <= HOMmax) {
    for(k = 1; k <= lat.conf.Cell_nLoc; k++) {
      M = dynamic_cast<MpoleType*>(lat.elems[k]);
      if ((lat.elems[k]->Pkind == Mpole) && (M->n_design == type))
	set_bnL_rms_elem(lat, lat.elems[k]->Fnum, lat.elems[k]->Knum, n, bnL,
			 anL, new_rnd);
    }
  } else
    printf("Bad type argument to set_bnL_rms_type()\n");
}


void set_bnr_sys_elem(LatticeType &lat, const int Fnum, const int Knum,
		      const int n, const double bnr, const double anr)
{
  int       nd;
  MpoleType *M;
  bool      prt = false;

  if (n < 1) {
    std::cout << "set_bnr_sys_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  M = dynamic_cast<MpoleType*>(lat.elems[lat.Elem_GetPos(Fnum, Knum)]);
  nd = M->n_design;
  // errors are relative to design values for (Dip, Quad, Sext, ...)
  M->PBsys[HOMmax+n] = bnr*M->PBpar[HOMmax+nd];
  M->PBsys[HOMmax-n] = anr*M->PBpar[HOMmax+nd];

  lat.SetPB(Fnum, Knum, n); lat.SetPB(Fnum, Knum, -n);

  if (prt)
    printf("set the n=%d component of %s to %e %e %e\n",
	   n, lat.elems[lat.Elem_GetPos(Fnum, Knum)]->PName,
	   bnr, M->PBpar[HOMmax+nd], M->PBsys[HOMmax+n]);
}


void set_bnr_sys_fam(LatticeType &lat, const int Fnum, const int n,
		     const double bnr, const double anr)
{
  int k;

  if (n < 1) {
    std::cout << "set_bnr_sys_fam: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  for (k = 1; k <= lat.GetnKid(Fnum); k++)
    set_bnr_sys_elem(lat, Fnum, k, n, bnr, anr);
}


void set_bnr_sys_type(LatticeType &lat, const int type, const int n,
		      const double bnr, const double anr)
{
  long int  k;
  MpoleType *M;

  if (n < 1) {
    std::cout << "set_bnr_sys_type: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  if (type >= Dip && type <= HOMmax) {
    for(k = 1; k <= lat.conf.Cell_nLoc; k++) {
      M = dynamic_cast<MpoleType*>(lat.elems[k]);
      if ((lat.elems[k]->Pkind == Mpole) && (M->n_design == type))
	set_bnr_sys_elem(lat, lat.elems[k]->Fnum, lat.elems[k]->Knum, n, bnr,
			 anr);
    }
  } else
    printf("Bad type argument to set_bnr_sys_type()\n");
}


void set_bnr_rms_elem(LatticeType &lat, const int Fnum, const int Knum,
		      const int n, const double bnr, const double anr,
		      const bool new_rnd)
{
  int       nd;
  MpoleType *M;

  bool prt = false;

  if (n < 1) {
    std::cout << "set_bnr_rms_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  M = dynamic_cast<MpoleType*>(lat.elems[lat.Elem_GetPos(Fnum, Knum)]);
  nd = M->n_design;
  // errors are relative to design values for (Dip, Quad, Sext, ...)
  if (nd == Dip) {
    M->PBrms[HOMmax+n] = bnr*M->Pirho; M->PBrms[HOMmax-n] = anr*M->Pirho;
  } else {
    M->PBrms[HOMmax+n] = bnr*M->PBpar[HOMmax+nd];
    M->PBrms[HOMmax-n] = anr*M->PBpar[HOMmax+nd];
  }

  if(new_rnd){
    if (normal) {
      M->PBrnd[HOMmax+n] = normranf(); M->PBrnd[HOMmax-n] = normranf();
    } else {
      M->PBrnd[HOMmax+n] = ranf(); M->PBrnd[HOMmax-n] = ranf();
    }
  }

  lat.SetPB(Fnum, Knum, n); lat.SetPB(Fnum, Knum, -n);

  if (prt) {
    printf("set_bnr_rms_elem:  Fnum = %d, Knum = %d, n = %d, n_design = %d"
	   ", new_rnd = %d, r_# = (%e, %e)\n",
	   Fnum, Knum, n, nd, new_rnd,
	   M->PBrnd[HOMmax+n], M->PBrnd[HOMmax-n]);
    printf("  (bnr, anr) = (%e, %e), PBrms = (%e, %e), PB = (%e, %e)\n",
	   bnr, anr, M->PBrms[HOMmax+n], M->PBrms[HOMmax-n],
	   M->PB[HOMmax+n], M->PB[HOMmax-n]);
  }
}


void set_bnr_rms_fam(LatticeType &lat, const int Fnum, const int n,
		     const double bnr, const double anr, const bool new_rnd)
{
  int k;

  if (n < 1) {
    std::cout << "set_bnr_rms_fam: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  for (k = 1; k <= lat.GetnKid(Fnum); k++)
    set_bnr_rms_elem(lat, Fnum, k, n, bnr, anr, new_rnd);
}


void set_bnr_rms_type(LatticeType &lat, const int type, const int n,
		      const double bnr, const double anr, const bool new_rnd)
{
  long int  k;
  MpoleType *M;

  if (n < 1) {
    std::cout << "set_bnr_rms_type: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  if (type >= Dip && type <= HOMmax) {
    for(k = 1; k <= lat.conf.Cell_nLoc; k++) {
      M = dynamic_cast<MpoleType*>(lat.elems[k]);
      if ((lat.elems[k]->Pkind == Mpole) && (M->n_design == type))
	set_bnr_rms_elem(lat, lat.elems[k]->Fnum, lat.elems[k]->Knum, n, bnr,
			 anr, new_rnd);
    }
  } else
    printf("Bad type argument to set_bnr_rms_type()\n");
}

//------------------------------------------------------------------------------

double get_Wiggler_BoBrho(LatticeType &lat, const int Fnum, const int Knum)
{
  WigglerType *W;

  W = dynamic_cast<WigglerType*>(lat.elems[lat.Elem_GetPos(Fnum, Knum)]);
  return W->BoBrhoV[0];
}


void set_Wiggler_BoBrho(LatticeType &lat, const int Fnum, const int Knum,
			const double BoBrhoV)
{
  WigglerType *W;

  W = dynamic_cast<WigglerType*>(lat.elems[lat.Elem_GetPos(Fnum, Knum)]);
  W->BoBrhoV[0] = BoBrhoV;
  W->PBW[HOMmax+Quad] = -sqr(BoBrhoV)/2.0;
  lat.SetPB(Fnum, Knum, Quad);
}


void set_Wiggler_BoBrho(LatticeType &lat, const int Fnum, const double BoBrhoV)
{
  int k;

  for (k = 1; k <= lat.GetnKid(Fnum); k++)
    set_Wiggler_BoBrho(lat, Fnum, k, BoBrhoV);
}


void set_ID_scl(LatticeType &lat, const int Fnum, const int Knum,
		const double scl)
{
  int           k;
  WigglerType   *W, *Wp;
  InsertionType *ID;
  FieldMapType  *FM;

  switch (lat.elems[lat.Elem_GetPos(Fnum, Knum)]->Pkind) {
  case Wigl:
    // scale the ID field
    W = dynamic_cast<WigglerType*>(lat.elems[lat.Elem_GetPos(Fnum, Knum)]);
    Wp = dynamic_cast<WigglerType*>(lat.elemf[Fnum-1].ElemF);
    for (k = 0; k < W->n_harm; k++) {
      W->BoBrhoH[k] = scl*Wp->BoBrhoH[k];
      W->BoBrhoV[k] = scl*Wp->BoBrhoV[k];
    }
    break;
  case Insertion:
    ID = dynamic_cast<InsertionType*>(lat.elems[lat.Elem_GetPos(Fnum, Knum)]);
    ID->scaling = scl;
    break;
  case FieldMap:
    FM = dynamic_cast<FieldMapType*>(lat.elems[lat.Elem_GetPos(Fnum, Knum)]);
    FM->scl = scl;
    break;
  default:
    std::cout << "set_ID_scl: unknown element type" << std::endl;
    exit_(1);
    break;
  }
}


void set_ID_scl(LatticeType &lat, const int Fnum, const double scl)
{
  int  k;

  for (k = 1; k <= lat.GetnKid(Fnum); k++)
    set_ID_scl(lat, Fnum, k, scl);
}


void SetFieldValues_fam(LatticeType &lat, const int Fnum, const bool rms,
			const double r0, const int n, const double Bn,
			const double An, const bool new_rnd)
{
  int       N;
  double    bnr, anr;
  MpoleType *M;

  M = dynamic_cast<MpoleType*>(lat.elems[lat.Elem_GetPos(Fnum, 1)]);
  N = M->n_design;
  if (r0 == 0.0) {
    // input is: (b_n L), (a_n L)
    if(rms)
      set_bnL_rms_fam(lat, Fnum, n, Bn, An, new_rnd);
    else
      set_bnL_sys_fam(lat, Fnum, n, Bn, An);
  } else {
    bnr = Bn/pow(r0, (double)(n-N)); anr = An/pow(r0, (double)(n-N));
    if(rms)
      set_bnr_rms_fam(lat, Fnum, n, bnr, anr, new_rnd);
    else
      set_bnr_sys_fam(lat, Fnum, n, bnr, anr);
  }
}


void SetFieldValues_type(LatticeType &lat, const int N, const bool rms,
			 const double r0, const int n, const double Bn,
			 const double An, const bool new_rnd)
{
  double bnr, anr;

  if (r0 == 0.0) {
    // input is: (b_n L), (a_n L)
    if(rms)
      set_bnL_rms_type(lat, N, n, Bn, An, new_rnd);
    else
      set_bnL_sys_type(lat, N, n, Bn, An);
  } else {
    bnr = Bn/pow(r0, (double)(n-N)); anr = An/pow(r0, (double)(n-N));
    if(rms)
      set_bnr_rms_type(lat, N, n, bnr, anr, new_rnd);
    else
      set_bnr_sys_type(lat, N, n, bnr, anr);
  }
}


void SetFieldErrors(LatticeType &lat, const char *name, const bool rms,
		    const double r0, const int n, const double Bn,
		    const double An, const bool new_rnd)
{
  int Fnum;

  if (strcmp("all", name) == 0) {
    printf("all: not yet implemented\n");
  } else if (strcmp("dip", name) == 0) {
    SetFieldValues_type(lat, Dip, rms, r0, n, Bn, An, new_rnd);
  } else if (strcmp("quad", name) == 0) {
    SetFieldValues_type(lat, Quad, rms, r0, n, Bn, An, new_rnd);
  } else if (strcmp("sext", name) == 0) {
    SetFieldValues_type(lat, Sext, rms, r0, n, Bn, An, new_rnd);
  } else {
    Fnum = ElemIndex(name);
    if(Fnum > 0)
      SetFieldValues_fam(lat, Fnum, rms, r0, n, Bn, An, new_rnd);
    else
      printf("SetFieldErrors: undefined element %s\n", name);
  }
}


// closed orbit correction by n_orbit iterations
bool CorrectCOD(LatticeType &lat, const int n_orbit, const double scl)
{
  bool            cod;
  int             i;
  long int        lastpos;
  Vector2         mean, sigma, max;
  ss_vect<double> ps;

  // ps.zero(); lat.Cell_Pass(0, conf.Cell_nLoc, ps, lastpos);
  // for (i = 1; i <= n_orbit; i++) {
  //   lstc(1, lastpos); lstc(2, lastpos);

  //   ps.zero(); lat.Cell_Pass(0, conf.Cell_nLoc, ps, lastpos);
  // }
  // if (false) prt_cod("cod.out", conf.bpm, true);
 
  cod = lat.getcod(0e0, lastpos);
  if (cod) {
    codstat(lat, mean, sigma, max, lat.conf.Cell_nLoc, true);
    printf("\n");
    printf("Initial RMS orbit (all):    x = %7.1e mm, y = %7.1e mm\n",
	   1e3*sigma[X_], 1e3*sigma[Y_]);
    codstat(lat, mean, sigma, max, lat.conf.Cell_nLoc, false);
    printf("\n");
    printf("Initial RMS orbit (BPMs):   x = %7.1e mm, y = %7.1e mm\n",
	   1e3*sigma[X_], 1e3*sigma[Y_]);

    for (i = 1; i <= n_orbit; i++){
      lsoc(lat, 1, scl); lsoc(lat, 2, scl);
      cod = lat.getcod(0e0, lastpos);
      if (!cod) break;

      if (cod) {
	codstat(lat, mean, sigma, max, lat.conf.Cell_nLoc, false);
	printf("Corrected RMS orbit (BPMs): x = %7.1e mm, y = %7.1e mm\n",
	       1e3*sigma[X_], 1e3*sigma[Y_]);
	if (i == n_orbit) {
	  codstat(lat, mean, sigma, max, lat.conf.Cell_nLoc, true);
	  printf("\n");
	  printf("Corrected RMS orbit (all):  x = %7.1e mm, y = %7.1e mm\n",
		 1e3*sigma[X_], 1e3*sigma[Y_]);
	}
      }
    }
  }

  return cod;
}


void write_misalignments(LatticeType &lat, const char* filename) {
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

void prt_beamsizes(LatticeType &lat, const int cnt)
{
  int   k;
  FILE  *fp;
  char fname [30];

  sprintf(fname,"%s_%d.out", beam_envelope_file, cnt);
  fp = file_write(fname);

  fprintf(fp,
	  "# k    name    s    s_xx    s_pxpx    s_xpx    s_yy    s_pypy"
	  "    s_ypy    theta_xy    s_xy\n");
  for(k = 0; k <= lat.conf.Cell_nLoc; k++)
    fprintf(fp,"%4d %10s %e %e %e %e %e %e %e %e %e\n",
	    k, lat.elems[k]->PName, lat.elems[k]->S,
	    lat.elems[k]->sigma[x_][x_], lat.elems[k]->sigma[px_][px_],
	    lat.elems[k]->sigma[x_][px_],
	    lat.elems[k]->sigma[y_][y_], lat.elems[k]->sigma[py_][py_],
	    lat.elems[k]->sigma[y_][py_],
	    atan2(2e0*lat.elems[k]->sigma[x_][y_],
		  lat.elems[k]->sigma[x_][x_]-lat.elems[k]->sigma[y_][y_])
	    /2e0*180.0/M_PI,
	    lat.elems[k]->sigma[x_][y_]);

  fclose(fp);

  sprintf(fname,"%s_%d.out","misalignments", cnt);
  write_misalignments(lat, fname);
}


double f_int_Touschek(const double u)
{
  double f;

  if (u > 0.0)
    f = (1.0/u-log(1.0/u)/2.0-1.0)*exp(-u_Touschek/u);
  else
    f = 0.0;

  return f;
}


double Touschek_loc(LatticeType &lat, const long int i, const double gamma,
		    const double delta_RF,
		    const double eps_x, const double eps_y,
		    const double sigma_delta, const double sigma_s,
		    const bool ZAP_BS)
{
  int    k;
  double sigma_x, sigma_y, sigma_xp, curly_H, dtau_inv;
  double alpha[2], beta[2], eta[2], etap[2];

  if ((i < 0) || (i > lat.conf.Cell_nLoc)) {
    std::cout << "Touschek_loc: undefined location " << i << std::endl;
    exit(1);
  }

  if (!ZAP_BS) {
    curly_H = get_curly_H(lat.elems[i]->Alpha[X_], lat.elems[i]->Beta[X_],
			  lat.elems[i]->Eta[X_], lat.elems[i]->Etap[X_]);

    // Compute beam sizes for given hor/ver emittance, sigma_s,
    // and sigma_delta (for x ~ 0): sigma_x0' = sqrt(eps_x/beta_x).
    sigma_x =
      sqrt(lat.elems[i]->Beta[X_]*eps_x+sqr(lat.elems[i]->Eta[X_]*sigma_delta));
    sigma_y = sqrt(lat.elems[i]->Beta[Y_]*eps_y);
    sigma_xp = (eps_x/sigma_x)*sqrt(1e0+curly_H*sqr(sigma_delta)/eps_x);
  } else {
    // ZAP averages the optics functions over an element instead of the
    // integrand; incorrect.

    for (k = 0; k < 2; k++) {
      alpha[k] = (lat.elems[i-1]->Alpha[k]+lat.elems[i]->Alpha[k])/2e0;
      beta[k] = (lat.elems[i-1]->Beta[k]+lat.elems[i]->Beta[k])/2e0;
      eta[k] = (lat.elems[i-1]->Eta[k]+lat.elems[i]->Eta[k])/2e0;
      etap[k] = (lat.elems[i-1]->Etap[k]+lat.elems[i]->Etap[k])/2e0;
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


double Touschek(LatticeType &lat, const double Qb, const double delta_RF,
		const double eps_x, const double eps_y,
		const double sigma_delta, const double sigma_s)
{
  // Note, ZAP (LBL-21270) averages the optics functions over an element
  // instead of the integrand; incorrect.  Hence, the Touschek lifetime is
  // overestimated by ~20%.

  long int i;
  double   p1, p2, dtau_inv, tau_inv;

  const bool   ZAP_BS = false;
  const double gamma  = 1e9*lat.conf.Energy/m_e, N_e = Qb/q_e;

  printf("\n");
  printf("Qb = %4.2f nC, delta_RF = %4.2f%%"
	 ", eps_x = %9.3e m.rad, eps_y = %9.3e m.rad\n",
	 1e9*Qb, 1e2*delta_RF, eps_x, eps_y);
  printf("sigma_delta = %8.2e, sigma_s = %4.2f mm\n",
	 sigma_delta, 1e3*sigma_s);

  // Integrate around the lattice with Trapezoidal rule; for nonequal steps.

  p1 = Touschek_loc(lat, 0, gamma, delta_RF, eps_x, eps_y, sigma_delta, sigma_s,
		    ZAP_BS);

  tau_inv = 0e0;
  for(i = 1; i <= lat.conf.Cell_nLoc; i++) {
    p2 = Touschek_loc(lat, i, gamma, delta_RF, eps_x, eps_y, sigma_delta,
		      sigma_s,
		      ZAP_BS);

    if (!ZAP_BS)
      dtau_inv = (p1+p2)/2e0;
    else
      dtau_inv = p2;

    tau_inv += dtau_inv*lat.elems[i]->PL; p1 = p2;

    if (false) {
      dtau_inv *=
	N_e*sqr(r_e)*c0/(8.0*M_PI*cube(gamma)*sigma_s)/sqr(delta_RF);

      printf("%4ld %9.3e\n", i, dtau_inv);
    }
  }

  tau_inv *=
    N_e*sqr(r_e)*c0/(8.0*M_PI*cube(gamma)*sigma_s)
    /(sqr(delta_RF)*lat.elems[lat.conf.Cell_nLoc]->S);

  printf("\n");
  printf("Touschek lifetime [hrs]: %10.3e\n", 1e0/(3600e0*tau_inv));

  return(1e0/(tau_inv));
}


void mom_aper(LatticeType &lat, double &delta, double delta_RF,
	      const long int k, const int n_turn, const bool positive)
{
  // Binary search to determine momentum aperture at location k.
  int      j;
  long int lastpos;
  double   delta_min, delta_max;
  psVector x;

  const double eps = 1e-4;

  delta_min = 0.0; delta_max = positive ? fabs(delta_RF) : -fabs(delta_RF);
  while (fabs(delta_max-delta_min) > eps) {
    delta = (delta_max+delta_min)/2.0;

    // propagate initial conditions
    CopyVec(6, lat.conf.CODvect, x); lat.Cell_Pass(0, k, x, lastpos);
    // generate Touschek event
    x[delta_] += delta;

    // complete one turn
    lat.Cell_Pass(k+1, lat.conf.Cell_nLoc, x, lastpos);
    if (lastpos < lat.conf.Cell_nLoc)
      // particle lost
      delta_max = delta;
    else {
      // track
      for(j = 0; j < n_turn; j++) {
	lat.Cell_Pass(0, lat.conf.Cell_nLoc, x, lastpos);

	if ((delta_max > delta_RF) || (lastpos < lat.conf.Cell_nLoc)) {
	  // particle lost
	  delta_max = delta;
	  break;
	}
      }

      if ((delta_max <= delta_RF) && (lastpos == lat.conf.Cell_nLoc))
	// particle not lost
	delta_min = delta;
    }
  }
}


double Touschek(LatticeType &lat, const double Qb, const double delta_RF,
		const bool consistent,const double eps_x, const double eps_y,
		const double sigma_delta, double sigma_s, const int n_turn,
		const bool aper_on, double sum_delta[][2],
		double sum2_delta[][2])
{
  bool     cav, aper;
  long int k;
  double   tau_inv, delta_p, delta_m, curly_H0, curly_H1, L;
  double   sigma_x, sigma_y, sigma_xp;
  FILE     *outf;

  const bool   prt       = false;
  const string file_name = "touschek.out";
  const double
    eps   = 1e-12,
    gamma = 1e9*lat.conf.Energy/m_e,
    N_e   = Qb/q_e;

  cav = lat.conf.Cavity_on; aper = lat.conf.Aperture_on;

  lat.conf.Cavity_on = true;

  lat.Ring_GetTwiss(true, 0.0);

  lat.conf.Aperture_on = aper_on;

  printf("\nQb = %4.2f nC, delta_RF = %4.2f%%"
	 ", eps_x = %9.3e m.rad, eps_y = %9.3e m.rad\n",
	 1e9*Qb, 1e2*delta_RF, eps_x, eps_y);
  printf("sigma_delta = %8.2e, sigma_s = %4.2f mm\n",
	 sigma_delta, 1e3*sigma_s);

  printf("\nMomentum aperture:\n");

  delta_p = delta_RF; mom_aper(lat, delta_p, delta_RF, 0, n_turn, true);
  delta_m = -delta_RF; mom_aper(lat, delta_m, delta_RF, 0, n_turn, false);
  delta_p = min(delta_RF, delta_p); delta_m = max(-delta_RF, delta_m);
  sum_delta[0][0] += delta_p; sum_delta[0][1] += delta_m;
  sum2_delta[0][0] += sqr(delta_p); sum2_delta[0][1] += sqr(delta_m);

  outf = file_write(file_name.c_str());

  tau_inv = 0e0; curly_H0 = -1e30;
  for (k = 1; k <= lat.conf.Cell_nLoc; k++) {
    L = lat.elems[k]->PL;

    curly_H1 = get_curly_H(lat.elems[k]->Alpha[X_], lat.elems[k]->Beta[X_],
			   lat.elems[k]->Eta[X_], lat.elems[k]->Etap[X_]);

    if (fabs(curly_H0-curly_H1) > eps) {
      mom_aper(lat, delta_p, delta_RF, k, n_turn, true);
      delta_m = -delta_p; mom_aper(lat, delta_m, delta_RF, k, n_turn, false);
      delta_p = min(delta_RF, delta_p); delta_m = max(-delta_RF, delta_m);
      printf("%4ld %6.2f %3.2lf%% %3.2lf%%\n",
	     k, lat.elems[k]->S, 1e2*delta_p, 1e2*delta_m);
      curly_H0 = curly_H1;
    }

    sum_delta[k][X_] += delta_p; sum_delta[k][Y_] += delta_m;
    sum2_delta[k][X_] += sqr(delta_p); sum2_delta[k][Y_] += sqr(delta_m);
    fprintf(outf, "%4ld %7.2f %5.3f %6.3f\n",
	    k, lat.elems[k]->S, 1e2*sum_delta[k][X_], 1e2*sum_delta[k][Y_]);
    fflush(outf);
    if (prt)
      printf("%4ld %6.2f %3.2lf %3.2lf\n",
	     k, lat.elems[k]->S, 1e2*delta_p, 1e2*delta_m);

    if (!consistent) {
      // Compute beam sizes for given hor/ver emittance, sigma_s,
      // and sigma_delta (for x ~ 0): sigma_x0' = sqrt(eps_x/beta_x).
      sigma_x =
	sqrt(lat.elems[k]->Beta[X_]*eps_x
	     +sqr(sigma_delta*lat.elems[k]->Eta[X_]));
      sigma_y = sqrt(lat.elems[k]->Beta[Y_]*eps_y);
      sigma_xp = (eps_x/sigma_x)*sqrt(1e0+curly_H1*sqr(sigma_delta)/eps_x);
    } else {
      // use self-consistent beam sizes
      sigma_x = sqrt(lat.elems[k]->sigma[x_][x_]);
      sigma_y = sqrt(lat.elems[k]->sigma[y_][y_]);
      sigma_xp = sqrt(lat.elems[k]->sigma[px_][px_]);
      sigma_s = sqrt(lat.elems[k]->sigma[ct_][ct_]);
    }

    u_Touschek = sqr(delta_p/(gamma*sigma_xp));

    tau_inv +=
      dqromb(f_int_Touschek, 0e0, 1e0)
      /(sigma_x*sigma_xp*sigma_y*sqr(delta_p))*L;

    if (prt) fflush(stdout);
  }

  fclose(outf);

  tau_inv *=
    N_e*sqr(r_e)*c0/(8.0*M_PI*cube(gamma)*sigma_s)
    /lat.elems[lat.conf.Cell_nLoc]->S;

  printf("\n");
  printf("Touschek lifetime [hrs]: %4.2f\n", 1e0/(3600e0*tau_inv));

  lat.conf.Cavity_on = cav; lat.conf.Aperture_on = aper;

  return 1/tau_inv;
}


double f_IBS(const double chi_m)
{
  // Interpolated integral (V. Litvinenko).

  double f, ln_chi_m;

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
  int    k;
  double f;

  const int    n_decades = 30;
  const double base      = 10e0;

  f = 0e0;
  for (k = 0; k <= n_decades; k++) {
    if (k == 0)
      f += dqromo(f_int_IBS, chi_m, pow(base, (double)k), dmidsql);
    else
      f += dqromb(f_int_IBS, pow(base, (double)(k-1)), pow(base, (double)k));
  }

  return f;
}


void IBS(LatticeType &lat, const double Qb, const double eps_SR[], double eps[],
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

  long int k;
  double   D_x, D_delta, b_max, L, gamma_z, a;
  double   sigma_x, sigma_xp, sigma_y, sigma_s, sigma_delta;
  double   incr, curly_H, eps_IBS[3];
  double   sigma_s_SR, sigma_delta_SR;

  const bool   integrate = false;
  const double gamma = 1e9*lat.conf.Energy/m_e, N_b = Qb/q_e;

  // bunch size
  gamma_z = (1.0+sqr(lat.conf.alpha_z))/lat.conf.beta_z;

  sigma_s_SR = sqrt(lat.conf.beta_z*eps_SR[Z_]);
  sigma_delta_SR = sqrt(gamma_z*eps_SR[Z_]);

  sigma_s = sqrt(lat.conf.beta_z*eps[Z_]); sigma_delta = sqrt(gamma_z*eps[Z_]);

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
	   lat.conf.alpha_z, lat.conf.beta_z);
    printf("sigma_s_SR     = %9.3e mm,   sigma_s     = %9.3e mm\n",
	   1e3*sigma_s_SR, 1e3*sigma_s);
    printf("sigma_delta_SR = %9.3e,      sigma_delta = %9.3e\n",
	   sigma_delta_SR, sigma_delta);
  }

  D_delta = 0.0; D_x = 0.0;
  for(k = 0; k <= lat.conf.Cell_nLoc; k++) {
    L = lat.elems[k]->PL;

    curly_H = get_curly_H(lat.elems[k]->Alpha[X_], lat.elems[k]->Beta[X_],
			  lat.elems[k]->Eta[X_], lat.elems[k]->Etap[X_]);

    // Compute beam sizes for given hor/ver emittance, sigma_s,
    // and sigma_delta (for x ~ 0): sigma_x0' = sqrt(eps_x/beta_x).
    sigma_x =
      sqrt(lat.elems[k]->Beta[X_]*eps[X_]
	   +sqr(lat.elems[k]->Eta[X_]*sigma_delta));
    sigma_xp = (eps[X_]/sigma_x)*sqrt(1.0+curly_H*sqr(sigma_delta)/eps[X_]);
    sigma_y = sqrt(lat.elems[k]->Beta[Y_]*eps[Y_]);

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
    /(32.0*M_PI*cube(gamma)*sigma_s*lat.elems[lat.conf.Cell_nLoc]->S);

  // eps_x*D_X
  D_x *= a;
  // Compute eps_IBS.
  eps_IBS[X_] = sqrt(D_x*lat.conf.tau[X_]/2e0);
  // Solve for eps_x
  eps[X_] = eps_SR[X_]*(1.0+sqrt(1.0+4.0*sqr(eps_IBS[X_]/eps_SR[X_])))/2.0;

  // compute diffusion coeffs.
  D_x /= eps[X_]; D_delta *= a/eps[X_];

  eps[Y_] = eps_SR[Y_]/eps_SR[X_]*eps[X_];

  sigma_delta = sqrt(sqr(sigma_delta_SR)+D_delta*lat.conf.tau[Z_]/2e0);
  eps[Z_] = sqr(sigma_delta)/gamma_z;

  sigma_s = sqrt(lat.conf.beta_z*eps[Z_]);

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
  int    k;
  double f;

  const int    n    = 30;
  const double decade = 10e0;

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


void IBS_BM(LatticeType &lat, const double Qb, const double eps_SR[],
	    double eps[], const bool prt1, const bool prt2)
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

  long int k;
  int      i;
  double   gamma_z, sigma_s_SR, sigma_delta_SR, sigma_s, sigma_delta;
  double   V, beta_m[2], sigma_m[2], alpha[2], beta[2], eta[2], etap[2];
  double   T_trans, rho, lambda_D, r_max;
  double   r_min, r_min_Cl, r_min_QM, log_Coulomb;
  double   L, curly_H[2], phi[2], dtau_inv[3], tau_inv[3], Gamma;
  double   a_BM, b_BM, c_BM, a2_CM, b2_CM, D_CM;
  double   D_x, D_delta, eps_IBS[3];

  const bool   ZAP_BS = true;
  const int    model = 2; // 1: Bjorken-Mtingwa, 2: Conte-Martini, 3: MAD-X
  const double gamma = 1e9*lat.conf.Energy/m_e;
  const double beta_rel = sqrt(1e0-1e0/sqr(gamma));
  const double N_b = Qb/q_e, q_i = 1e0;

  // bunch size
  gamma_z = (1e0+sqr(lat.conf.alpha_z))/lat.conf.beta_z;

  sigma_s_SR = sqrt(lat.conf.beta_z*eps_SR[Z_]);
  sigma_delta_SR = sqrt(gamma_z*eps_SR[Z_]);

  sigma_s = sqrt(lat.conf.beta_z*eps[Z_]); sigma_delta = sqrt(gamma_z*eps[Z_]);

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
	   lat.conf.alpha_z, lat.conf.beta_z);
    printf("sigma_s_SR     = %9.3e mm,   sigma_s     = %9.3e mm\n",
	   1e3*sigma_s_SR, 1e3*sigma_s);
    printf("sigma_delta_SR = %9.3e,      sigma_delta = %9.3e\n",
	   sigma_delta_SR, sigma_delta);
  }

  // Compute the Coulomb log

  for (i = 0; i < 2; i++)
    beta_m[i] = 0e0;

  for(k = 0; k <= lat.conf.Cell_nLoc; k++)
    for (i = 0; i < 2; i++)
      beta_m[i] += lat.elems[k]->Beta[i]*lat.elems[k]->PL;

  for (i = 0; i < 2; i++) {
    beta_m[i] /= lat.elems[lat.conf.Cell_nLoc]->S;
    sigma_m[i] = sqrt(beta_m[i]*eps[i]);
  }

  V = 8e0*pow(M_PI, 3e0/2e0)*sigma_m[X_]*sigma_m[Y_]*sigma_s;
  rho = N_b/V;
  // Transverse temperature (p. 171, "ZAP User's Manual"):
  //   T_trans [eV] = 2*E_kin_trans [eV]->
  T_trans = (gamma*1e9*lat.conf.Energy-m_e)*eps[X_]/beta_m[X_];
  lambda_D = 743.4e-2/q_i*sqrt(T_trans/(1e-6*rho));
  r_max = min(sigma_m[X_], lambda_D);

  // Classical distance of closest approach:
  //   d = z*Z*e^2/(2*pi*eps_0*m*v^2) = z*Z*e/(4*pi*eps_0*T_trans [eV]).
  r_min_Cl = 1.44e-9*sqr(q_i)/T_trans;

  // Test particle de Broglie wavelength in the scattering center's rest system
  // p. 3, J.A. Krommes "An Introduction to the Physics of the Coulomb
  // Logarithm with Emphasis on Quantum-Mechanical Effects" arXiv:1806.04990
  // (2018):
  //   b_min = h/mu*u
  // where mu is the reduced mass. Similarly, from p. 34, "NRL Plasma
  // Formularly" (2013):
  //   b_min = h_bar/2*mu*u = h_bar*c0/sqrt(m*c0^2 [eV]*T_trans [eV])
  // i.e., the constant should be h_bar*c0 = 1.973e-7; not 1.973e-13.
  if (!ZAP_BS)
    r_min_QM = 1.973e-13/(2e0*sqrt(T_trans*m_e));
  else
    // ZAP implementation, a Bug (variable name is "COULOG"):
    //   [ETRANS] = eV, [E0] = Mev => *1e6, not *1e-6).
    r_min_QM = 1.973e-13/(2e0*sqrt(1e-12*T_trans*m_e));

  // Interpolation between classical and quantum-mechanical limits.
  r_min = max(r_min_Cl, r_min_QM);

  log_Coulomb = log(r_max/r_min);

  for (i = 0; i < 3; i++)
    tau_inv[i] = 0e0;

  for(k = 1; k <= lat.conf.Cell_nLoc; k++) {
    L = lat.elems[k]->PL;

    for (i = 0; i < 2; i++){
      if (!ZAP_BS) {
	alpha[i] = lat.elems[k]->Alpha[i]; beta[i] = lat.elems[k]->Beta[i];
	eta[i] = lat.elems[k]->Eta[i]; etap[i] = lat.elems[k]->Etap[i];
      } else {
	// Note, ZAP averages the optics functions over an element instead of
	// the integrand; incorrect.
	alpha[i] = (lat.elems[k-1]->Alpha[i]+lat.elems[k]->Alpha[i])/2e0;
	beta[i] = (lat.elems[k-1]->Beta[i]+lat.elems[k]->Beta[i])/2e0;
	eta[i] = (lat.elems[k-1]->Eta[i]+lat.elems[k]->Eta[i])/2e0;
	etap[i] = (lat.elems[k-1]->Etap[i]+lat.elems[k]->Etap[i])/2e0;
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
      /(M_PI*cube(2e0*beta_rel)*pow(gamma, 4e0)
	*lat.elems[lat.conf.Cell_nLoc]->S);

  D_x = eps[X_]*tau_inv[X_]; D_delta = eps[Z_]*tau_inv[Z_];

  // eps_x*D_x
  D_x *= eps[X_];
  // Compute eps_IBS.
  eps_IBS[X_] = sqrt(D_x*lat.conf.tau[X_]/2e0);
  // Solve for eps_x
  eps[X_] = eps_SR[X_]*(1e0+sqrt(1e0+4e0*sqr(eps_IBS[X_]/eps_SR[X_])))/2e0;
  // Compute D_x.
  D_x /= eps[X_];

  eps[Y_] = eps_SR[Y_]/eps_SR[X_]*eps[X_];

  D_delta = tau_inv[Z_]*sqr(sigma_delta);
  sigma_delta =
    sqrt((D_delta+2e0/lat.conf.tau[Z_]*sqr(sigma_delta_SR))
	 *lat.conf.tau[Z_]/2e0);
  eps[Z_] = sqr(sigma_delta)/gamma_z;

  sigma_s = sqrt(lat.conf.beta_z*eps[Z_]);

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
  int i, k;

  i = 0;
  while (name[i] == ' ')
    i++;

  for (k = i; k <= (int)strlen(name)+1; k++)
    name[k-i] = name[k];
}


void get_bn(LatticeType &lat, const char file_name[], int n, const bool prt)
{
  char   line[max_str], str[max_str], str1[max_str], *token, *name, *p;
  int    n_prm, Fnum, Knum, order;
  double bnL, bn, C, L;
  FILE   *inf, *fp_lat;

  inf = file_read(file_name); fp_lat = file_write("get_bn.lat");

  no_sxt(lat);

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
    strcpy(str, name); Fnum = ElemIndex(str);
    strcpy(str1, name); upr_case(str1);
    token = strtok_r(NULL, ")", &p); sscanf(token, "%d", &Knum);
    strtok_r(NULL, "=", &p); token = strtok_r(NULL, "\n", &p);
    sscanf(token, "%lf %d", &bnL, &order);
    if (prt) printf("%6s(%2d) = %10.6f %d\n", name, Knum, bnL, order);
   if (Fnum != 0) {
      if (order == 0)
        SetL(lat, Fnum, bnL);
      else
        SetbnL(lat, Fnum, order, bnL);

      L = GetL(lat, Fnum, 1);
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

  C = lat.elems[lat.conf.Cell_nLoc]->S; recalc_S(lat);
  if (prt)
    printf("New Cell Length: %5.3f (%5.3f)\n",
	   lat.elems[lat.conf.Cell_nLoc]->S, C);

  fclose(inf); fclose(fp_lat);
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

    for (i = 0; i < nv_; i++)
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


double get_chi2(long int n, double x[], double y[], long int m, psVector b)
{
  /* Compute chi2 for polynomial fit */

  int    i, j;
  double sum, z;

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

  int    i, j, k;
  Matrix T1;

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
	printf("%3d %10.3e +/-%8.2e\n", i, b[i], sigma*sqrt(chi2*T1[i][i]));
    }
  } else {
    printf("pol_fit: Matrix is singular\n");
    exit_(1);
  }
}


void get_ksi2(LatticeType &lat, const double d_delta)
{
  const int n_points = 20, order = 5;

  int      i, n;
  double   delta[2*n_points+1], nu[2][2*n_points+1], sigma;
  psVector b;
  FILE     *fp;

  fp = file_write("chrom2.out");
  n = 0;
  for (i = -n_points; i <= n_points; i++) {
    n++; delta[n-1] = i*(double)d_delta/(double)n_points;
    lat.Ring_GetTwiss(false, delta[n-1]);
    nu[0][n-1] = lat.conf.TotalTune[X_]; nu[1][n-1] = lat.conf.TotalTune[Y_];
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
  bool lost;
  int  k;

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


bool get_nu(LatticeType &lat, const double Ax, const double Ay,
	    const double delta,const double eps, double &nu_x, double &nu_y)
{
  const int n_turn = 512, n_peaks = 5;

  bool     lost, ok_x, ok_y;
  char     str[max_str];
  int      i;
  long int lastpos, lastn, n;
  double   x[n_turn], px[n_turn], y[n_turn], py[n_turn];
  double   nu[2][n_peaks], A[2][n_peaks];
  psVector x0;
  FILE     *fp;

  const bool prt = false;
  const char file_name[] = "track.out";

  // complex FFT in Floquet space
  x0[x_] = Ax; x0[px_] = 0.0; x0[y_] = Ay; x0[py_] = 0.0;
  LinTrans(4, lat.conf.Ascrinv, x0);
  track(file_name, lat, x0[x_], x0[px_], x0[y_], x0[py_], delta,
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


void dnu_dA(LatticeType &lat, const double Ax_max, const double Ay_max,
	    const double delta, const int n_ampl)
{
  bool     ok;
  int      i;
  double   nu_x, nu_y, Ax, Ay, Jx, Jy;
  psVector ps;
  FILE     *fp;

  const double A_min  = 0.1e-3;
//  const double  eps0   = 0.04, eps   = 0.02;
//  const double  eps0   = 0.025, eps   = 0.02;
//   const double  eps0   = 0.04, eps   = 0.015;
  const double eps = 0.01;

  lat.Ring_GetTwiss(false, 0.0);

  if (trace) printf("dnu_dAx\n");

  nu_x = fract(lat.conf.TotalTune[X_]); nu_y = fract(lat.conf.TotalTune[Y_]);

  fp = file_write("dnu_dAx.out");
  fprintf(fp, "#   A_x        A_y        J_x        J_y      nu_x    nu_y\n");
  fprintf(fp, "#\n");
  fprintf(fp, "%10.3e %10.3e %10.3e %10.3e %8.6f %8.6f\n",
	  0e0, 0e0, 0e0, 0e0, fract(nu_x), fract(nu_y));

  Ay = A_min;
  for (i = 1; i <= n_ampl; i++) {
    Ax = -i*Ax_max/n_ampl;
    ps[x_] = Ax; ps[px_] = 0e0; ps[y_] = Ay; ps[py_] = 0e0;
    getfloqs(lat, ps);
    Jx = (sqr(ps[x_])+sqr(ps[px_]))/2.0; Jy = (sqr(ps[y_])+sqr(ps[py_]))/2.0;
    ok = get_nu(lat, Ax, Ay, delta, eps, nu_x, nu_y);
    if (ok)
      fprintf(fp, "%10.3e %10.3e %10.3e %10.3e %8.6f %8.6f\n",
	      1e3*Ax, 1e3*Ay, 1e6*Jx, 1e6*Jy, fract(nu_x), fract(nu_y));
    else
      fprintf(fp, "# %10.3e %10.3e particle lost\n", 1e3*Ax, 1e3*Ay);
  }

  if (trace) printf("\n");

  nu_x = fract(lat.conf.TotalTune[X_]); nu_y = fract(lat.conf.TotalTune[Y_]);

  fprintf(fp, "\n");
  fprintf(fp, "%10.3e %10.3e %10.3e %10.3e %8.6f %8.6f\n",
	  0e0, 0e0, 0e0, 0e0, fract(nu_x), fract(nu_y));

  Ay = A_min;
  for (i = 0; i <= n_ampl; i++) {
    Ax = i*Ax_max/n_ampl;
    ps[x_] = Ax; ps[px_] = 0e0; ps[y_] = Ay; ps[py_] = 0e0;
    getfloqs(lat, ps);
    Jx = (sqr(ps[x_])+sqr(ps[px_]))/2.0; Jy = (sqr(ps[y_])+sqr(ps[py_]))/2.0;
    ok = get_nu(lat, Ax, Ay, delta, eps, nu_x, nu_y);
    if (ok)
      fprintf(fp, "%10.3e %10.3e %10.3e %10.3e %8.6f %8.6f\n",
	      1e3*Ax, 1e3*Ay, 1e6*Jx, 1e6*Jy, fract(nu_x), fract(nu_y));
    else
      fprintf(fp, "# %10.3e %10.3e particle lost\n", 1e3*Ax, 1e3*Ay);
  }

  fclose(fp);

  if (trace) printf("dnu_dAy\n");

  nu_x = fract(lat.conf.TotalTune[X_]); nu_y = fract(lat.conf.TotalTune[Y_]);

  fp = file_write("dnu_dAy.out");
  fprintf(fp, "#   A_x        A_y      nu_x    nu_y\n");
  fprintf(fp, "#\n");
  fprintf(fp, "%10.3e %10.3e %10.3e %10.3e %8.6f %8.6f\n",
	  0e0, 0e0, 0e0, 0e0, fract(nu_x), fract(nu_y));

  Ax = A_min;
  for (i = 1; i <= n_ampl; i++) {
    Ay = -i*Ay_max/n_ampl;
    Jx = pow(Ax, 2.0)/(2.0*lat.elems[lat.conf.Cell_nLoc]->Beta[X_]);
    Jy = pow(Ay, 2.0)/(2.0*lat.elems[lat.conf.Cell_nLoc]->Beta[Y_]);
    ok = get_nu(lat, Ax, Ay, delta, eps, nu_x, nu_y);
    if (ok)
      fprintf(fp, "%10.3e %10.3e %10.3e %10.3e %8.6f %8.6f\n",
	      1e3*Ax, 1e3*Ay, 1e6*Jx, 1e6*Jy, fract(nu_x), fract(nu_y));
    else
      fprintf(fp, "# %10.3e %10.3e particle lost\n", 1e3*Ax, 1e3*Ay);
  }

  if (trace) printf("\n");

  nu_x = fract(lat.conf.TotalTune[X_]); nu_y = fract(lat.conf.TotalTune[Y_]);

  fprintf(fp, "\n");
  fprintf(fp, "%10.3e %10.3e %10.3e %10.3e %8.6f %8.6f\n",
	  0e0, 0e0, 0e0, 0e0, fract(nu_x), fract(nu_y));

  Ax = A_min;
  for (i = 0; i <= n_ampl; i++) {
    Ay = i*Ay_max/n_ampl;
    Jx = pow(Ax, 2.0)/(2.0*lat.elems[lat.conf.Cell_nLoc]->Beta[X_]);
    Jy = pow(Ay, 2.0)/(2.0*lat.elems[lat.conf.Cell_nLoc]->Beta[Y_]);
    ok = get_nu(lat, Ax, Ay, delta, eps, nu_x, nu_y);
    if (ok)
      fprintf(fp, "%10.3e %10.3e %10.3e %10.3e %8.6f %8.6f\n",
	      1e3*Ax, 1e3*Ay, 1e6*Jx, 1e6*Jy, fract(nu_x), fract(nu_y));
    else
      fprintf(fp, "# %10.3e %10.3e particle lost\n", 1e3*Ax, 1e3*Ay);
  }

  fclose(fp);
}


bool orb_corr(LatticeType &lat, const int n_orbit)
{
  bool    cod = false;
  int     i;
  long    lastpos;
  Vector2 xmean, xsigma, xmax;

  printf("\n");
  lat.conf.CODvect.zero();
  for (i = 1; i <= n_orbit; i++) {
    cod = lat.getcod(0.0, lastpos);
    if (cod) {
      codstat(lat, xmean, xsigma, xmax, lat.conf.Cell_nLoc, false);
      printf("\n");
      printf("RMS orbit [mm]: (%8.1e+/-%7.1e, %8.1e+/-%7.1e)\n",
	     1e3*xmean[X_], 1e3*xsigma[X_], 1e3*xmean[Y_], 1e3*xsigma[Y_]);
      lsoc(lat, 1, 1e0); lsoc(lat, 2, 1e0);
      cod = lat.getcod(0.0, lastpos);
      if (cod) {
	codstat(lat, xmean, xsigma, xmax, lat.conf.Cell_nLoc, false);
	printf("RMS orbit [mm]: (%8.1e+/-%7.1e, %8.1e+/-%7.1e)\n",
	       1e3*xmean[X_], 1e3*xsigma[X_], 1e3*xmean[Y_], 1e3*xsigma[Y_]);
      } else
	printf("orb_corr: failed\n");
    } else
      printf("orb_corr: failed\n");
  }

  prt_cod(lat, "orb_corr.out", true);

  return cod;
}


void get_alphac(LatticeType &lat)
{
  ElemType Cell;

  lat.getelem(lat.conf.Cell_nLoc, &Cell);
  lat.conf.Alphac = lat.conf.OneTurnMat[ct_][delta_]/Cell.S;
}


void get_alphac2(LatticeType &lat)
{
  /* Note, do not extract from M[5][4], i.e. around delta
     dependent fixed point.                                */

  const int    n_points = 5;
  const double d_delta  = 2e-2;

  int      i, j, n;
  long int lastpos;
  double   delta[2*n_points+1], alphac[2*n_points+1], sigma;
  psVector x, b;
  ElemType Cell;

  lat.conf.pathlength = false;
  lat.getelem(lat.conf.Cell_nLoc, &Cell); n = 0;
  for (i = -n_points; i <= n_points; i++) {
    n++; delta[n-1] = i*(double)d_delta/(double)n_points;
    for (j = 0; j < nv_; j++)
      x[j] = 0.0;
    x[delta_] = delta[n-1];
    lat.Cell_Pass(0, lat.conf.Cell_nLoc, x, lastpos);
    alphac[n-1] = x[ct_]/Cell.S;
  }
  pol_fit(n, delta, alphac, 3, b, sigma, true);
  printf("\n");
  printf("alphac = %10.3e + %10.3e*delta + %10.3e*delta^2\n",
	 b[1], b[2], b[3]);
}


#if 0

double f_bend(double b0L[])
{
  long int lastpos;
  psVector ps;

  const int n_prt = 10;

  n_iter_Cart++;

  SetbnL_sys(lat, Fnum_Cart, Dip, b0L[1]);

  ps.zero();
  lat.Cell_Pass(lat.Elem_GetPos(Fnum_Cart, 1)-1, lat.Elem_GetPos(Fnum_Cart, 1),
  	    ps, lastpos);

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
  const int n_prm = 1;

  int      iter;
  double   *b0L, **xi, fret;
  psVector ps;

  const double ftol = 1e-15;

  b0L = dvector(1, n_prm); xi = dmatrix(1, n_prm, 1, n_prm);

  std::cout << std::endl;
  std::cout << "bend_cal: " << lat.elemf[Fnum-1].ElemF->PName << ":"
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
  MpoleType *M;

  for (k = 1; k <= lat.conf.Elem_nFam; k++) {
    M = dynamic_cast<MpoleType*>(lat.elemf[k-1].ElemF);
    if ((lat.elemf[k-1].ElemF->Pkind == Mpole) &&
	(M->Pirho != 0.0) && (M->PBpar[Quad+HOMmax] != 0.0))
      if (lat.elemf[k-1].nKid > 0) bend_cal_Fam(k);
  }
}

#endif

double h_ijklm(const tps &h, const int i, const int j, const int k,
	       const int l, const int m)
{
  int      i1;
  long int jj[ss_dim];

  for (i1 = 0; i1 < nv_tps; i1++)
    jj[i1] = 0;
  jj[x_] = i; jj[px_] = j; jj[y_] = k; jj[py_] = l; jj[delta_] = m;
  return h[jj];
}


void set_tune(LatticeType &lat, const char file_name1[],
	      const char file_name2[], const int n)
{
  const int n_b2 = 8;

  char          line[max_str], names[n_b2][max_str];
  int           j, k, Fnum;
  double        b2s[n_b2], nu[2];
  std::ifstream inf1, inf2;
  FILE          *fp_lat;

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
	Fnum = ElemIndex(names[k]);
	set_bn_design_fam(lat, Fnum, Quad, b2s[k], 0.0);

	fprintf(fp_lat, "%s: Quadrupole, L = %8.6f, K = %10.6f, N = Nquad"
		", Method = Meth;\n",
		names[k], lat.elemf[Fnum-1].ElemF->PL, b2s[k]);
      }
      break;
    }
  }

  fclose(fp_lat);
}


void get_map_twiss(const ss_vect<tps> &M,
		   double beta0[], double beta1[], double nu[], bool stable[])
{
  // Assumes that alpha_0 = alpha_1 = 0.
  int    k;
  double cosmu, sinmu;

  const bool prt = true;

  for (k = 0; k < 2; k++) {
    stable[k] = M[2*k][2*k]*M[2*k+1][2*k+1] <= 1e0;
    if (stable[k]) {
      cosmu = sqrt(M[2*k][2*k]*M[2*k+1][2*k+1]);
      if (M[2*k][2*k] < 0e0) cosmu = -cosmu;
      nu[k] = acos(cosmu)/(2e0*M_PI);
      if (M[2*k][2*k+1] < 0e0) nu[k] = 1e0 - nu[k];
      if (nu[k] != 0e0) {
	sinmu = sin(2e0*M_PI*nu[k]);
	beta0[k] = M[2*k][2*k+1]*M[2*k+1][2*k+1]/(cosmu*sinmu);
	beta1[k] = M[2*k][2*k]*M[2*k][2*k+1]/(cosmu*sinmu);
      } else {
	printf("\nget_map_twiss: dnu is zero for plane %d\n", k);
	exit(1);
      }
    } else {
      printf("\nget_map_twiss: M unstable for plane %d\n", k);
      exit(1);
    }
  }

  if (prt) {
    printf("\n  nu    = [%8.5f, %8.5f]\n", nu[X_], nu[Y_]);
    printf("  beta  = [%8.5f, %8.5f]\n", beta0[X_], beta0[Y_]);
    printf("  beta  = [%8.5f, %8.5f]\n", beta1[X_], beta1[Y_]);
  }
}


void set_map(MapType *Map)
{
  // Set phase-space rotation.
  int          k;
  double       dnu[2], alpha[2], beta[2], eta_x, etap_x, cosmu, sinmu;
  ss_vect<tps> Id, Id_eta;

  Id.identity();

  for (k = 0; k < 2; k++) {
    dnu[k] = Map->dnu[k]; alpha[k] = Map->alpha[k]; beta[k] = Map->beta[k];
  }
  eta_x = Map->eta_x; etap_x = Map->etap_x;

  Map->M.identity();
  for (k = 0; k < 2; k++) {
    cosmu = cos(2e0*M_PI*dnu[k]); sinmu = sin(2e0*M_PI*dnu[k]);
    Map->M[2*k]   = (cosmu+alpha[k]*sinmu)*Id[2*k] + beta[k]*sinmu*Id[2*k+1];
    Map->M[2*k+1] =
      -(1e0+sqr(alpha[k]))*sinmu/beta[k]*Id[2*k]
      + (cosmu-alpha[k]*sinmu)*Id[2*k+1];
  }

  // Zero linear dispersion contribution.
  Id_eta.identity();
  Id_eta[x_] = eta_x*Id[delta_]; Id_eta[px_] = etap_x*Id[delta_];
  Id_eta = Map->M*Id_eta;
  Map->M[x_]  -= (Id_eta[x_][delta_]-eta_x)*Id[delta_];
  Map->M[px_] -= (Id_eta[px_][delta_]-etap_x)*Id[delta_];
}


void set_map(LatticeType &lat, const int Fnum, const double dnu[])
{
  long int loc;
  int      j, k;
  MapType  *Map;

  for (j = 1; j <= lat.GetnKid(Fnum); j++) {
    loc = lat.Elem_GetPos(Fnum, j);
    Map = dynamic_cast<MapType*>(lat.elems[loc]);
    for (k = 0; k < 2; k++) {
      Map->dnu[k]   = dnu[k];
      Map->alpha[k] = lat.elems[loc]->Alpha[k];
      Map->beta[k]  = lat.elems[loc]->Beta[k];
    }
    Map->eta_x  = lat.elems[loc]->Eta[X_];
    Map->etap_x = lat.elems[loc]->Etap[X_];

    set_map(Map);
  }
}


void set_map_per(MapType *Map,
		 const double alpha0[], const double beta0[],
		 const double eta0[], const double etap0[])
{
  // Phase advance is set to zero.
  int          k;
  ss_vect<tps> Id, Id_eta;

  Id.identity(); Map->M.identity();
  for (k = 0; k < 2; k++)
    Map->M[2*k+1] += 2e0*alpha0[k]/beta0[k]*Id[2*k];

  Id_eta.identity();
  Id_eta[x_] = eta0[X_]*Id[delta_]; Id_eta[px_] = etap0[X_]*Id[delta_];
  Id_eta = Map->M*Id_eta;
  Map->M[px_] -= (Id_eta[px_][delta_]+etap0[X_])*Id[delta_];
}


void set_map_per(LatticeType &lat, const int Fnum, const double alpha0[],
		 const double beta0[], const double eta0[],
		 const double etap0[])
{
  int     j;
  MapType *Map;

  for (j = 1; j <= lat.GetnKid(Fnum); j++) {
    Map = dynamic_cast<MapType*>(lat.elems[lat.Elem_GetPos(Fnum, j)]);
    set_map_per(Map, alpha0, beta0, eta0, etap0);
  }
}


void set_map_reversal(LatticeType &lat, ElemType *Cell)
{
  long int loc, lastpos;
  MapType  *Map;

  danot_(1);
  Map = dynamic_cast<MapType*>(Cell);
  Map->M.identity();
  loc = lat.Elem_GetPos(Cell->Fnum, Cell->Knum);
  lat.Cell_Pass(0, loc-1, Map->M, lastpos);
  Map->M = Inv(Map->M);
}


void set_map_reversal(LatticeType &lat, const long int Fnum)
{
  int j;

  for (j = 1; j <= lat.GetnKid(Fnum); j++)
    set_map_reversal(lat, lat.elems[lat.Elem_GetPos(Fnum, j)]);
}


void setmp(long ilat, long m, long n, double rr, double bnoff, double cmn)
{
  // ilat index of element in lattice
  // m multipole to be set 
  // n base multipole for which the coefficient is defined
  //     ! n or m < 0 means skew!
  // rr reference radius for definition of coefficient
  // bnoff offset for base multipole 
  // cmn coefficient from magnet group in "units". normalization:
  // base multipole=10000
  //    cmn = (nAn/r) from magnet group
  // ==> bm = cmn*1e-4 / rr^(|m|-|n|) * (bn-bnoff)

  long      i;
  double    bn, dbn, mp;
  ElemType  elem;
  MpoleType *M;

  mp = cmn * 1e-4;
  for (i=abs(n); i<abs(m); i++) {mp/=rr;}

  //--> get multipole n of element at position ilat
  lat.getelem(ilat, &elem);
  M = dynamic_cast<MpoleType*>(&elem);
  bn = M->PBpar[n+HOMmax];

  if (bn>0) {dbn=bn-bnoff;} else {dbn=bn+bnoff;}

  mp=mp*dbn;
  //-->save mp as multipole m, skew if m<0
  M->PBpar[m+HOMmax] = mp;
  M->PB[m+HOMmax] =
    M->PBpar[m+HOMmax] + M->PBsys[m+HOMmax] +
    M->PBrms[m+HOMmax]*M->PBrnd[m+HOMmax];
  if (abs(m) > M->Porder && M->PB[m+HOMmax] != 0e0)
    M->Porder = abs(m);
}


void setmpall(LatticeType &lat, double rref)
{
  ElemType cell;
  long     i;
  bool     mset;
  
  if (true) {
    for (i = 0; i <= lat.conf.Cell_nLoc; i++) {
      lat.getelem(i, &cell);
      mset=false;
      if (cell.Pkind == Mpole) {
	if (strncmp(cell.PName,"bn",2) == 0) {
	  // BN magnet
	  setmp(i, 3, 1, rref , 0.0, 2.0940E1) ;
	  setmp(i, 5, 1, rref , 0.0, 1.1376E1) ;
          mset=true;       
	}
	if (strncmp(cell.PName,"vb",2) == 0) {
	  //VB magnet (mirrorplate version) 
	  setmp(i, 3, 1, rref , 0.0, -39.86) ;
	  setmp(i, 4, 2, rref , 0.0, -90.93) ;
	  setmp(i, 5, 1, rref , 0.0, -3.17) ;
	  setmp(i, 6, 2, rref , 0.0, -22.10) ;
          mset=true;       
	}
	if (strncmp(cell.PName,"qp",2) == 0) {
	  //QP main quads
	  setmp(i, 6, 2, rref , 6.75, -26) ; //dodeka suppressed at 55 T/m
	  setmp(i,10, 2, rref , 0.0 ,-1.65) ;
          mset=true;       
	}
	if (strncmp(cell.PName,"s",1) == 0) {
	  //S sextupole
	  setmp(i, 9, 3, rref , 0.0, -8.1);
	  setmp(i,15, 3, rref , 0.0, -24.9);
          mset=true;       
	}
	if (strncmp(cell.PName,"qa",2) == 0) {
	  //QA tuning quads and skew quads
	  setmp(i, 6, 2, rref , 0.0, 1435);
	  setmp(i,-6,-2, rref , 0.0,-1435);
          mset=true;       
	}
	if (strncmp(cell.PName,"cs",2) == 0) {
	  //CS skew quads
	  setmp(i, 6, 2, rref , 0.0, 1435);
	  setmp(i,-6,-2, rref , 0.0,-1435);
          mset=true;
	}
        if (mset) {
	  printf("%s systematic multipoles defined\n", cell.PName);
	}
      }
    }
  }
}
