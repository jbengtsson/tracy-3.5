#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


const bool
  mat_meth      = false,
  zero_b_3      = false,
  zero_b_4      = false,
  fit_nu        = false,
  fit_xi        = false,
  ps_rot        = false,
  chk_mpole_sym = false,  // Requires super period.
  chk_dnu       = false,  // Requires super period.
  comp_H_long   = false,
  Deta          = false,
  get_tol       = false,
  phiob_2xL     = false;  

const int
  n_aper  = 25,
  n_track = 1000;

const double
  dnu[] = {0.0, 0.0},
  nu[]  = {57.202/20.0+0.5/20.0, 20.7435/20.0-0.5/20.0};


void set_ps_rot(const string &fam_name, const double dnu_x, const double dnu_y)
{
  const int
    Fnum = ElemIndex(fam_name.c_str());
  const double
    dnu_0[] = {0e0, 0e0},
    dnu[]   = {dnu_x, dnu_y};

  set_map(Fnum, dnu_0);
  printf("\ntweak_nu:\n");
  printf("  dnu = [%8.5f, %8.5f]\n", dnu[X_], dnu[Y_]);
  set_map(Fnum, dnu);
  Ring_GetTwiss(true, 0e0);
  printglob();
 }


double* get_dnu_straight(const int loc)
{
  static double
    dnu[2] =  {
    Cell[globval.Cell_nLoc].Nu[X_]-Cell[loc].Nu[X_],
    Cell[globval.Cell_nLoc].Nu[Y_]-Cell[loc].Nu[Y_]
  };

  printf("\nget_dnu_straight:\n");
  printf(" nu               = [%7.5f, %7.5f]\n",
	 Cell[globval.Cell_nLoc].Nu[X_], Cell[globval.Cell_nLoc].Nu[Y_]);
  printf(" dnu              = [%7.5f, %7.5f]\n",
	 Cell[loc].Nu[X_], Cell[loc].Nu[Y_] );
  printf(" dnu_1/2_straight = [%7.5f, %7.5f]\n", dnu[X_], dnu[Y_]);

  return dnu;
}


void set_dnu_straight(const string &fam_name, const int loc)
{
  const double dnu_half_straight[] = {0.25, 0.125};
  // const double dnu_half_straight[] = {0.5, 0.25};

  double*       dnu;
  static double dnu_ps_rot[2];

  dnu = get_dnu_straight(loc);
  for (int k = 0; k < 2; k++)
    dnu_ps_rot[k] = dnu_half_straight[k] - dnu[k];

  printf("\nset_dnu_straight:\n");
  printf(" dnu_ps_rot = [%7.5f, %7.5f]\n", dnu_ps_rot[X_], dnu_ps_rot[Y_]);

  set_ps_rot(fam_name, dnu_ps_rot[X_], dnu_ps_rot[Y_]);
}


void fit_nu_jb
(const std::vector<int> &Fnum_b_2, const double dnu_x, const double dnu_y,
 const double db_2L)
{
  int    n_b_2, j, k;
  double **A, **U, **V, *w, *dnu, *db_2, b_2, a_2;

  const bool   prt = !false;
  const int    m   = 2;
  const double
    dnu_vec[]  = {dnu_x, dnu_y},
    svd_cut    = 1e-10;

  n_b_2 = Fnum_b_2.size();

  A = dmatrix(1, m, 1, n_b_2);
  U = dmatrix(1, m, 1, n_b_2);
  V = dmatrix(1, n_b_2, 1, n_b_2);
  w = dvector(1, n_b_2);
  dnu = dvector(1, m);
  db_2 = dvector(1, n_b_2);

  if (prt)
    printf("\nfit_nu_jb: dnu_sp = [%7.5f, %7.5f]\n", dnu_x, dnu_y);
  for (k = 1; k <= n_b_2; k++) {
    set_dbnL_design_fam(Fnum_b_2[k-1], Quad, db_2L, 0e0);
    Ring_GetTwiss(false, 0e0);
    if (prt)
      printf("\nfit_nu_jb: nu1+ = [%9.5f, %9.5f]\n",
	     globval.TotalTune[X_], globval.TotalTune[Y_]);

    for (j = 1; j <= m; j++)
      A[j][k] = globval.TotalTune[j-1];
    set_dbnL_design_fam(Fnum_b_2[k-1], Quad, -2e0*db_2L, 0e0);
    Ring_GetTwiss(false, 0e0);
    if (prt)
      printf("fit_nu_jb: nu1- = [%9.5f, %9.5f]\n",
	 globval.TotalTune[X_], globval.TotalTune[Y_]);
    for (j = 1; j <= m; j++) {
      A[j][k] -= globval.TotalTune[j-1];
      A[j][k] /= 2e0*db_2L;
    }

    set_dbnL_design_fam(Fnum_b_2[k-1], Quad, db_2L, 0e0);
  }

  Ring_GetTwiss(false, 0e0);
  if (prt)
    printf("\nfit_nu_jb: nu1  = [%9.5f, %9.5f]\n",
	   globval.TotalTune[X_], globval.TotalTune[Y_]);
  for (j = 1; j <= m; j++)
    dnu[j] = dnu_vec[j-1];

  dmcopy(A, m, n_b_2, U);
  dsvdcmp(U, m, n_b_2, w, V);

  if (prt) {
    printf("\nfit_nu_jb:\n  singular values:\n");
    for (j = 1; j <= n_b_2; j++) {
      printf("    %9.3e", w[j]);
      if (w[j] < svd_cut) {
	w[j] = 0e0;
	printf(" (zeroed)");
      }
      printf("\n");
    }
  }

  dsvbksb(U, w, V, m, n_b_2, dnu, db_2);

  if (prt) {
    dmdump(stdout, "\nA:", A, 2, 2, "%11.3e");
    dvdump(stdout, "\ndb_2:", db_2, 2, "%11.3e");
  }

  for (k = 1; k <= n_b_2; k++)
    set_dbnL_design_fam(Fnum_b_2[k-1], Quad, db_2[k], 0e0);

  if (prt) {
    printf("\n  b_2:\n");
    for (k = 0; k < n_b_2; k++) {
      get_bn_design_elem(Fnum_b_2[k], 1, Quad, b_2, a_2);
      printf("    %-8s %10.5f\n", ElemFam[Fnum_b_2[k]-1].ElemF.PName, b_2);
    }
    printf("\n");
  }

  free_dmatrix(A, 1, m, 1, n_b_2);
  free_dmatrix(U, 1, m, 1, n_b_2);
  free_dmatrix(V, 1, n_b_2, 1, n_b_2);
  free_dvector(w, 1, n_b_2);
  free_dvector(dnu, 1, m);
  free_dvector(db_2, 1, n_b_2);
}


void fit_xi_jb
(const std::vector<int> &Fnum_b_3, const double xi_x, const double xi_y,
 const double db_3L)
{
  int    n_b_3, j, k;
  double **A, **U, **V, *w, *dxi, *db_3, b_3, a_3;

  const bool   prt = !false;
  const int    m   = 2;
  const double
    xi[]  = {xi_x, xi_y},
    svd_cut = 1e-10;

  n_b_3 = Fnum_b_3.size();

  A = dmatrix(1, m, 1, n_b_3);
  U = dmatrix(1, m, 1, n_b_3);
  V = dmatrix(1, n_b_3, 1, n_b_3);
  w = dvector(1, n_b_3);
  dxi = dvector(1, m);
  db_3 = dvector(1, n_b_3);

  // Zero sextupoles to track linear chromaticity.
  if (false) no_sxt();

  for (k = 1; k <= n_b_3; k++) {
    set_dbnL_design_fam(Fnum_b_3[k-1], Sext, db_3L, 0e0);
    Ring_Getchrom(0e0);
    if (prt)
      printf("\nfit_xi_jb: xi1+ = [%9.5f, %9.5f]\n",
	     globval.Chrom[X_], globval.Chrom[Y_]);

    for (j = 1; j <= m; j++)
      A[j][k] = globval.Chrom[j-1];
    set_dbnL_design_fam(Fnum_b_3[k-1], Sext, -2e0*db_3L, 0e0);
    Ring_Getchrom(0e0);
    if (prt)
      printf("fit_xi_jb: xi1- = [%9.5f, %9.5f]\n",
	 globval.Chrom[X_], globval.Chrom[Y_]);
    for (j = 1; j <= m; j++) {
      A[j][k] -= globval.Chrom[j-1];
      A[j][k] /= 2e0*db_3L;
    }

    set_dbnL_design_fam(Fnum_b_3[k-1], Sext, db_3L, 0e0);
  }

  Ring_Getchrom(0e0);
  if (prt)
    printf("\nfit_xi_jb: xi1  = [%9.5f, %9.5f]\n",
	   globval.Chrom[X_], globval.Chrom[Y_]);
  for (j = 1; j <= m; j++)
    dxi[j] = -(globval.Chrom[j-1]-xi[j-1]);

  dmcopy(A, m, n_b_3, U);
  dsvdcmp(U, m, n_b_3, w, V);

  if (prt) {
    printf("\nfit_xi_jb:\n  singular values:\n");
    for (j = 1; j <= n_b_3; j++) {
      printf("    %9.3e", w[j]);
      if (w[j] < svd_cut) {
	w[j] = 0e0;
	printf(" (zeroed)");
      }
      printf("\n");
    }
  }

  dsvbksb(U, w, V, m, n_b_3, dxi, db_3);

  if (prt) {
    dmdump(stdout, "\nA:", A, 2, 2, "%11.3e");
    dvdump(stdout, "\ndb_3:", db_3, 2, "%11.3e");
  }

  for (k = 1; k <= n_b_3; k++)
    set_dbnL_design_fam(Fnum_b_3[k-1], Sext, db_3[k], 0e0);

  if (prt) {
    printf("\n  b_3:\n");
    for (k = 0; k < n_b_3; k++) {
      get_bn_design_elem(Fnum_b_3[k], 1, Sext, b_3, a_3);
      printf("    %-8s %10.5f\n", ElemFam[Fnum_b_3[k]-1].ElemF.PName, b_3);
    }
    printf("\n");
  }

  free_dmatrix(A, 1, m, 1, n_b_3);
  free_dmatrix(U, 1, m, 1, n_b_3);
  free_dmatrix(V, 1, n_b_3, 1, n_b_3);
  free_dvector(w, 1, n_b_3);
  free_dvector(dxi, 1, m);
  free_dvector(db_3, 1, n_b_3);
}


void track(const int n_turn, const double Ax, const double Ay)
{
  const string
    file_name = "track.dat";
  const double
    A[] = {Ax, Ay};

  long int        lastpos;
  ss_vect<double> ps;
  ofstream        outf;

  file_wr(outf, file_name.c_str());
  ps.zero();
  for (int k = 0; k < 2; k++)
    ps[2*k] = A[k];
  outf << "# turn          x                     p_x"
       << "                     y                     p_y"
       << "                   delta                   c*t\n"
       << "#              [m]                   [rad]"
       << "                   [m]                   [rad]"
       << "                                          [m] \n";
    outf << scientific << setprecision(14)
	 << setw(4) << 0 << setw(23) << ps << "\n";
  for (int k = 1; k <= n_turn; k++) {
    Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
    outf << scientific << setprecision(14)
	 << setw(4) << k << setw(23) << ps << "\n";
  }
  outf.close();
}


void chk_optics(const double alpha_x, const double beta_x,
		const double eta_x, const double etap_x,
		const double alpha_y, const double beta_y,
		const double eta_y, const double etap_y)
{
  Vector2 alpha, beta, eta, etap;

  alpha[X_] = alpha_x;
  alpha[Y_] = alpha_y;
  beta[X_]  = beta_x;
  beta[Y_]  = beta_y;
  eta[X_]   = eta_x;
  eta[Y_]   = eta_y;
  etap[X_]  = etap_x;
  etap[Y_]  = etap_y;

  ttwiss(alpha, beta, eta, etap, 0e0);
}


void chk_phi()
{
  int    k;
  double dphi, phi, mphi;

  printf("\n");
  phi = 0e0; mphi = 0e0;
  for (k = 0; k <= globval.Cell_nLoc; k++) {
    if ((Cell[k].Elem.Pkind == Mpole) &&
	(Cell[k].Elem.M->Pirho != 0e0)) {
      dphi = Cell[k].Elem.PL*Cell[k].Elem.M->Pirho*180e0/M_PI;
      phi += dphi;
      if (dphi < 0e0) mphi += dphi;
    }
  }
  printf("\nphi = %8.6f phi- = %8.6f phi+ = %8.6f\n", phi, mphi, phi-mphi);
}


void get_phiob_2xL_ratios(void)
{
  // Get 1/(rho*b_2) ratios.
  double b_2, a_2;
  
  printf("\n");
  for (int k = 0; k <= globval.Cell_nLoc; k++) {
    if (Cell[k].Elem.Pkind == Mpole) {
      get_bn_design_elem(Cell[k].Fnum, 1, Quad, b_2, a_2);
      if ((Cell[k].Elem.M->Pirho != 0e0) && (b_2 != 0e0)) {
	auto L = Cell[k].Elem.PL;
	auto irho = Cell[k].Elem.M->Pirho;
	auto phi = L*irho*180e0/M_PI;
	auto phiob_2xL = irho/b_2;
	printf("  %8s L = %7.5f phi = %8.5f b_2xL = %8.5f"
	       " phi/(b_2*L) = %12.5e\n",
	       Cell[k].Elem.PName, L, phi, b_2*L, phiob_2xL);
      }
    }
  }
}


void chk_mpole_Fam(const int Fnum)
{
  // Assumes that the multipoles are split in halfs - i.e., to obtain the
  // linear optics at the centre.
  int    loc;
  double dnu[2], dnu_0[2];

  printf("\n   name        s    beta_x*eta_x  beta_y* eta_x"
	 "   dnu_x    dnu_y\n");
  for (auto k = 0; k < 2; k++)
    dnu_0[k] = NAN;
  for (auto j = 1; j <= GetnKid(Fnum); j += 2) {
    loc = Elem_GetPos(Fnum, j);
    for (auto k = 0; k < 2; k++) {
      dnu[k] = Cell[loc].Nu[k] - dnu_0[k];
      dnu_0[k] = Cell[loc].Nu[k];
    }

    printf("  %.8s %7.3f   %8.5f       %8.5f    %8.5f %8.5f\n",
	   Cell[loc].Elem.PName, Cell[loc].S,
	   Cell[loc].Beta[X_]*Cell[loc].Eta[X_],
	   Cell[loc].Beta[Y_]*Cell[loc].Eta[X_], dnu[X_], dnu[Y_]);
  }
}


int get_ElemIndex(string elem_name)
{
  const int Fnum = ElemIndex(elem_name.c_str());
  const int n_Kid = GetnKid(Fnum);
  if (n_Kid != 0) {
    return Fnum;
  } else {
    printf("\nchk_mpole_Fam: *** no kids for %s\n", elem_name.c_str());
    exit(1);
  }
}


void chk_mpole(const int lat_case)
{
  int              k;
  std::vector<int> Fnum;

  switch (lat_case) {
  case 1:
    Fnum.push_back(get_ElemIndex("s1_h2"));
    Fnum.push_back(get_ElemIndex("s2_h2"));
    Fnum.push_back(get_ElemIndex("s3_h2"));
    Fnum.push_back(get_ElemIndex("s4_h2"));
    break;
  case 2:
    Fnum.push_back(get_ElemIndex("s1_n1"));
    Fnum.push_back(get_ElemIndex("s2_n1"));
    Fnum.push_back(get_ElemIndex("s3_n1"));
    Fnum.push_back(get_ElemIndex("s4_n1"));
    break;
  case 3:
    Fnum.push_back(get_ElemIndex("s1_n1"));
    Fnum.push_back(get_ElemIndex("s2_n1"));
    Fnum.push_back(get_ElemIndex("s3a_n1"));
    Fnum.push_back(get_ElemIndex("s3b_n1"));
    Fnum.push_back(get_ElemIndex("s3c_n1"));
    Fnum.push_back(get_ElemIndex("s4a_n1"));
    Fnum.push_back(get_ElemIndex("s4b_n1"));
    break;
  case 4:
    Fnum.push_back(get_ElemIndex("sfm"));
    Fnum.push_back(get_ElemIndex("sfi"));
    Fnum.push_back(get_ElemIndex("sdqd_1"));
    Fnum.push_back(get_ElemIndex("sdqd_2"));
    Fnum.push_back(get_ElemIndex("sdqd_3"));
    Fnum.push_back(get_ElemIndex("sdqd_4"));
    Fnum.push_back(get_ElemIndex("sdqd_5"));
    Fnum.push_back(get_ElemIndex("sdendq"));
    Fnum.push_back(get_ElemIndex("sfo"));
    break;
  default:
    printf("\nchk_mpole: unknown lattice type\n");
    exit(1);
    break;
  }

  Ring_GetTwiss(true, 0e0);
 
  printf("\nMultipole Scheme:\n");
  for (k = 0; k < (int)Fnum.size(); k++)
    chk_mpole_Fam(Fnum[k]);
}


void chk_dnu_straight(const string &fam_name)
{
  const int loc = Elem_GetPos(ElemIndex(fam_name.c_str()), 1);

  double dnu[2];
  
  for (auto k = 0; k < 2; k++)
    dnu[k] = 2*Cell[loc].Nu[k];
  printf("\n  dnu = [%7.5f, %7.5f]\n", dnu[X_], dnu[Y_]);
}


void prt_b_n(void)
{
  const string file_name = "lat_bn.out"; 

  FILE* outf;

  outf = file_write(file_name.c_str());

  fprintf(outf, "#        name           s   code   phi       b_2"
	  "          b_3          b_4\n");
  fprintf(outf, "#                      [m]        [deg]    [1/m^2]"
	  "      [1/m^3]      [1/m^4]\n");
  for (auto k = 0; k <= globval.Cell_nLoc; k++)
    if (Cell[k].Elem.Pkind != Mpole)
      fprintf(outf, "%4d %15s %6.2f %4.1f %6.3f %12.5e %12.5e %12.5e\n",
	      k, Cell[k].Elem.PName, Cell[k].S, get_code(Cell[k]),
	      0e0, 0e0, 0e0, 0e0);
    else {
      auto phi = Cell[k].Elem.M->Pirho*Cell[k].Elem.PL*180e0/M_PI;
      fprintf(outf, "%4d %15s %6.2f %4.1f %6.3f %12.5e %12.5e %12.5e\n",
	      k, Cell[k].Elem.PName, Cell[k].S, get_code(Cell[k]),
	      phi, Cell[k].Elem.M->PBpar[Quad+HOMmax],
	      Cell[k].Elem.M->PBpar[Sext+HOMmax],
	      Cell[k].Elem.M->PBpar[Oct+HOMmax]);
    }

  fclose(outf);
}


psVector compute_alpha_c(void)
{
  // Note, do not extract from M[5][4], i.e. around delta dependent fixed
  // point.

  const int    n_points = 5;
  const double d_delta  = 2e-2;

  int      i, j, n;
  long int lastpos;
  double   delta[2*n_points+1], alphac[2*n_points+1], sigma;
  psVector x, b;
  CellType Cell;

  globval.pathlength = false;
  getelem(globval.Cell_nLoc, &Cell); n = 0;
  for (i = -n_points; i <= n_points; i++) {
    n++; delta[n-1] = i*(double)d_delta/(double)n_points;
    for (j = 0; j < nv_; j++)
      x[j] = 0e0;
    x[delta_] = delta[n-1];
    Cell_Pass(0, globval.Cell_nLoc, x, lastpos);
    alphac[n-1] = x[ct_]/Cell.S;
  }
  pol_fit(n, delta, alphac, 3, b, sigma, true);

  return b;
}


void compute_alpha_bucket()
{
  psVector alpha_c;

  alpha_c = compute_alpha_c();
  printf("\n  alphac    = %10.3e %+10.3e*delta %+10.3e*delta^2\n",
	 alpha_c[1], alpha_c[2], alpha_c[3]);
  printf("  RF alpha bucket to 2nd order [%%]      = [%5.1f, %5.1f]\n",
	 -1e2*alpha_c[1]/alpha_c[2], 1e2*alpha_c[1]/(2e0*alpha_c[2]));

  printf("  Unstable fixet point to 3rd order [%%] =  %5.1f\n",
	 -1e2*alpha_c[2]/(2e0*alpha_c[3])*
	 (1e0-sqrt(1e0-4e0*alpha_c[1]*alpha_c[3]/sqr(alpha_c[2]))));

  printf("  Shift due to alpha^(3)_c [%%] = [%5.1f, %5.1f]\n",
	 1e2*cube(alpha_c[1]/alpha_c[2])*alpha_c[3],
	 -1e2*sqr(alpha_c[1])*alpha_c[3]/(48e0*cube(alpha_c[2])));
}


double H_long
(const double phi, const double delta, const int h_rf, const double V_rf,
 const double phi_0, const psVector &alpha_c, const int n_alpha_c)
{
  const double
    E_0 = 1e9*globval.Energy;

  double H;

  H = V_rf/E_0*(cos(phi+phi_0)+phi*sin(phi_0));
  for (auto i = 2; i <= n_alpha_c+1; i++)
    H += 2e0*pi*h_rf*alpha_c[i-1]*pow(delta, (double)i)/i;
  return H;
}


void prt_H_long
(const int Fnum_cav, const int n, const double phi_max,
 const double delta_max, const int n_alpha_c, const bool neg_alpha_c)
{
  const string
    file_name = "H_long.dat";
  const long int
    loc = Elem_GetPos(Fnum_cav, 1);
  const CavityType
    *C = Cell[loc].Elem.C;
  const int
    h_RF  = C->harm_num;
  const double
    E_0 = 1e9*globval.Energy,
    U_0 = globval.U0,
    V_RF = C->V_RF,
    phi_0 = -fabs(asin(globval.U0/V_RF));

  double   phi, delta, H, delta_RF;
  psVector alpha_c;
  FILE     *outf;

  outf = file_write(file_name.c_str());

  alpha_c = compute_alpha_c();

  delta_RF =
    sqrt(-V_RF*cos(M_PI+phi_0)*(2e0-(M_PI-2e0*(M_PI+phi_0))*tan(M_PI+phi_0))
	 /(alpha_c[1]*M_PI*h_RF*E_0));
  printf("\nU_0 [keV]        = %3.1f\n", 1e-3*U_0);
  printf("phi_0 [deg]      = 180 %4.2f\n", phi_0*180e0/M_PI);
  printf("RF bucket height = %4.2f\n", 1e2*delta_RF);

  for (auto i = -n; i <= n ; i++) {
    for (auto j = -n; j <= n ; j++) {
      phi = i*phi_max*M_PI/(n*180e0);
      delta = j*delta_max/n;
      H = H_long(phi, delta, h_RF, V_RF, M_PI+phi_0, alpha_c, n_alpha_c);
      fprintf(outf, "  %8.2f %10.5f, %13.5e\n", phi*180e0/M_PI, 1e2*delta, H);
    }
    fprintf(outf, "\n");
  }

  fclose(outf);
}


void compute_Deta(const double delta)
{
  // Evaluate derivative; to avoid effect of tune shift.
  double         h;
  vector<double> Deta[2];
  FILE           *outf;

  const double d_delta   = 1e-5;
  const string file_name = "Deta.out";

  outf = file_write(file_name.c_str());

  printf("\nOptics for delta = %10.3e\n", d_delta);
  Ring_GetTwiss(true, d_delta);
  printglob();
  for (auto k = 0; k <= globval.Cell_nLoc; k++) {
    Deta[x_].push_back(Cell[k].Eta[X_]);
    Deta[px_].push_back(Cell[k].Etap[X_]);
  }
  printf("\nOptics for delta = %10.3e\n", -d_delta);
  Ring_GetTwiss(true, -d_delta);
  printglob();
  fprintf(outf, "#  k name                 s   type     eta_x        eta'_x"
	  "   D_delta eta_x  D_delta eta'_x   eta_x/rho   eta_x/rho eta^(2)_x"
	  "  D_delta beta_x/rho\n"
	        "#                        [m]            [m]"
	  "                       [m]\n");
  for (auto k = 0; k <= globval.Cell_nLoc; k++) {
    if (Cell[k].Elem.Pkind == Mpole)
      h = Cell[k].Elem.M->Pirho;
    else
      h = 0e0;
    Deta[x_][k] -= Cell[k].Eta[X_];
    Deta[x_][k] /= (2e0*d_delta);
    Deta[px_][k] -= Cell[k].Etap[X_];
    Deta[px_][k] /= (2e0*d_delta);
    fprintf(outf, "%4d %10s %8.3f %4.1f %12.5e %12.5e %12.5e   %12.5e   %12.5e"
	    "     %12.5e         %12.5e\n",
	    k, Cell[k].Elem.PName, Cell[k].S, get_code(Cell[k]),
	    Cell[k].Eta[x_], Cell[k].Etap[x_], Deta[x_][k], Deta[px_][k],
	    h*Cell[k].Eta[x_],
	    pow(Cell[k].Etap[x_], 2)/2e0, h*Deta[x_][k]);
  }

  fclose(outf);
}


void prt_cod_1(const char *file_name)
{
  FILE *outf;

  outf = file_write(file_name);

  fprintf(outf, "#    name                 s    code        x_cod"
	  "                 p_x,cod                y_cod"
	  "                p_y,cod\n");
  fprintf(outf, "#                        [m]                [m]"
	  "                   [rad]                  [m]"
	  "                  [rad]\n");
  for (auto i = 0; i <= globval.Cell_nLoc; i++)
    fprintf(outf,
	    "%4d %-15s %9.5f %4.1f %21.14e %21.14e %21.14e %21.14e\n",
	    i, Cell[i].Elem.PName, Cell[i].S, get_code(Cell[i]),
	    Cell[i].BeamPos[x_], Cell[i].BeamPos[px_], Cell[i].BeamPos[y_],
	    Cell[i].BeamPos[py_]);
  fclose(outf);
}


double get_phi(const int Fnum)
{
  const auto  loc = Elem_GetPos(Fnum, 1);
  return Cell[loc].Elem.PL*Cell[loc].Elem.M->Pirho*180e0/M_PI;
}


void compute_phi(const int k)
{
  // max_4u_sp_jb_5.
  const string dip_1[] = {
    "d1_u6", "d1_u5", "d1_u4", "d1_u3", "d1_u2", "d1_u1", "d1_0", "d1_d1",
    "d1_d2", "d1_d3", "d1_d4", "d1_d5"};
  // m4U_250220_h02_09_01_01_tracy-2.
  const string dip_2[] = {
    "d1_h2_sl_dm1", "d1_h2_sl_dm2", "d1_h2_sl_dm3", "d1_h2_sl_dm4",
    "d1_h2_sl_dm5", "d1_h2_sl_ds0", "d1_h2_sl_ds1", "d1_h2_sl_ds2",
    "d1_h2_sl_ds3", "d1_h2_sl_ds4", "d1_h2_sl_ds5", "d1_h2_sl_ds6"};
  // m4U_250316_h03_01_01_01_tracy-2.
  const string dip_3[] = {
    "d1_h2_sl_dm5", "d1_h2_sl_dm4", "d1_h2_sl_dm3", "d1_h2_sl_dm2",
    "d1_h2_sl_dm1", "d1_h2_sl_d0", "d1_h2_sl_ds1", "d1_h2_sl_ds2",
    "d1_h2_sl_ds3", "d1_h2_sl_ds4", "d1_h2_sl_ds5"};
  // m4U_250505_h02_12_01_01_tracy-2.
  const string dip_4[] = {
    "d1_h2_sl_dm1", "d1_h2_sl_dm2", "d1_h2_sl_dm3", "d1_h2_sl_dm4",
    "d1_h2_sl_dm5", "d1_h2_sl_ds0", "d1_h2_sl_ds1", "d1_h2_sl_ds2",
    "d1_h2_sl_ds3", "d1_h2_sl_ds4", "d1_h2_sl_ds5", "d1_h2_sl_ds6"};
 
  double              phi = 0e0;
  std::vector<string> dip_names;

  switch (k) {
  case 1:
    dip_names.assign(dip_1, dip_1+12);
    break;
  case 2:
    dip_names.assign(dip_2, dip_2+12);
    break;
  case 3:
    dip_names.assign(dip_3, dip_3+11);
    break;
  case 4:
    dip_names.assign(dip_4, dip_4+12);
    break;
  default:
    printf("\n*** compute_phi: unknown case\n");
    exit(1);
    break;
  }

  for (auto k = 0; k < dip_names.size(); k++)
    phi += get_phi(ElemIndex(dip_names[k]));
  printf("\ncompute_phi: phi_d1 = %5.3f\n", phi);
}


void fit_nu_jb_2(const double nu_x, const double nu_y)
{
  const int lat = 0;

  std::vector<int> Fnum;

  switch (lat) {
  case 0:
    Fnum.push_back(ElemIndex("q1_h2"));
    Fnum.push_back(ElemIndex("q2_h2"));

    Fnum.push_back(ElemIndex("r1_h2"));
    Fnum.push_back(ElemIndex("r2_h2"));
    Fnum.push_back(ElemIndex("r3_h2"));

    Fnum.push_back(ElemIndex("d1_h2_sl_dm5"));
    Fnum.push_back(ElemIndex("d1_h2_sl_dm4"));
    Fnum.push_back(ElemIndex("d1_h2_sl_dm3"));
    Fnum.push_back(ElemIndex("d1_h2_sl_dm2"));
    Fnum.push_back(ElemIndex("d1_h2_sl_dm1"));
    Fnum.push_back(ElemIndex("d1_h2_sl_d0"));
    Fnum.push_back(ElemIndex("d1_h2_sl_ds1"));
    Fnum.push_back(ElemIndex("d1_h2_sl_ds2"));
    Fnum.push_back(ElemIndex("d1_h2_sl_ds3"));
    Fnum.push_back(ElemIndex("d1_h2_sl_ds4"));
    Fnum.push_back(ElemIndex("d1_h2_sl_ds5"));

    Fnum.push_back(ElemIndex("d2_h2_sl_df0"));
    Fnum.push_back(ElemIndex("d2_h2_sl_df1"));
    Fnum.push_back(ElemIndex("d2_h2_sl_df2"));
    Fnum.push_back(ElemIndex("d2_h2_sl_df3"));
    Fnum.push_back(ElemIndex("d2_h2_sl_df4"));
    Fnum.push_back(ElemIndex("d2_h2_sl_df5"));
    Fnum.push_back(ElemIndex("d2_h2_sl_df6"));

    Fnum.push_back(ElemIndex("d3_h2_sl_df0"));
    Fnum.push_back(ElemIndex("d3_h2_sl_df1"));
    Fnum.push_back(ElemIndex("d3_h2_sl_df2"));
    Fnum.push_back(ElemIndex("d3_h2_sl_df3"));
    Fnum.push_back(ElemIndex("d3_h2_sl_df4"));
    Fnum.push_back(ElemIndex("d3_h2_sl_df5"));
    Fnum.push_back(ElemIndex("d3_h2_sl_df6"));

  }

  fit_nu_jb(Fnum, nu_x, nu_y, 1e-3);
  Ring_GetTwiss(true, 0e0);
  printglob();
}


void fit_xi_jb_2(const double xi_x, const double xi_y)
{
  std::vector<int> Fnum;

#if 0
  Fnum.push_back(ElemIndex("s1_h2"));
  Fnum.push_back(ElemIndex("s2_h2"));
  Fnum.push_back(ElemIndex("s3_h2"));
  Fnum.push_back(ElemIndex("s4_h2"));
#else
  Fnum.push_back(ElemIndex("sf_f"));
  Fnum.push_back(ElemIndex("sd_d"));
#endif  

  fit_xi_jb(Fnum, xi_x, xi_y, 1e0);
  Ring_GetTwiss(true, 0e0);
  printglob();
}


void compute_beta_beat(const std::vector<std::vector<double>> &beta_ref)
{
  int    n = 0;
  double sum[2] = {0e0, 0e0}, sum_2[2] = {0e0, 0e0}, mean[2], sigma[2];

  Ring_GetTwiss(true, 0e0);

  for (auto j = 0; j <= globval.Cell_nLoc; j++) {
    n++;
    for (auto k = 0; k < 2; k++) {
      sum[k] += Cell[j].Beta[k] - beta_ref[j][k];
      sum_2[k] += sqr((Cell[j].Beta[k] - beta_ref[j][k])/Cell[j].Beta[k]);
    }
  }
  for (auto k = 0; k < 2; k++) {
    mean[k] = sum[k]/n;
    sigma[k] = (n*sum_2[k]-sqr(sum[k]) >= 0e0)?
      sqrt((n*sum_2[k]-sqr(sum[k]))/(n*(n-1e0))) : 0e0;
  }

  printf("beta-beat [%%] = [%6.3f +/- %5.3f, %6.3f +/- %5.3f]\n",
	 1e2*mean[X_], 1e2*sigma[X_], 1e2*mean[Y_], 1e2*sigma[Y_]);
}


void get_b_2_tol(const double db_2_rms, const int n_aper, const int n_track)
{
  const string
    file_name  = "dynap";
  const bool
    Floq_space = false,
    cod        = true,
    prt        = false;
  const double
    r_0   = 5e-3,
    dr    = 0.1e-3,
    delta = 0e0;

  stringstream                     str;
  double                           x_aper[n_aper], y_aper[n_aper], DA;
  std::vector<double>              beta = {0e0, 0e0};
  std::vector<std::vector<double>> beta_ref;
  FILE                             *fp;

  globval.Cavity_on = false;
  Ring_GetTwiss(true, 0e0);
  for (auto j = 0; j <= globval.Cell_nLoc; j++) {
    beta.assign({Cell[j].Beta[X_], Cell[j].Beta[Y_]});
    beta_ref.push_back(beta);
  }

  globval.Cavity_on = true;

  str << scientific << setprecision(2) << file_name << "_"
      << setw(8) << db_2_rms << ".out";
  fp = file_write(str.str().c_str());

  set_bnr_rms_type(Dip,  Quad, db_2_rms, 0e0, true);
  set_bnr_rms_type(Quad, Quad, db_2_rms, 0e0, true);

  if (false) {
    str.str("");
    str.clear();
    str << scientific << setprecision(2) << file_name << "_"
	<< setw(8) << db_2_rms << ".dat";
    prtmfile(str.str().c_str());
  }

  dynap(fp, r_0, delta, dr, n_aper, n_track, x_aper, y_aper, Floq_space, cod,
	prt);
  fclose(fp);
  DA = get_aper(n_aper, x_aper, y_aper);
  printf("db_2_rms = %8.2e DA [mm^2] = %8.2e ", db_2_rms, 1e6*DA);
  compute_beta_beat(beta_ref);

  globval.Cavity_on = false;
}


void set_state(void)
{
  globval.H_exact        = false;
  globval.quad_fringe    = false;
  globval.Cavity_on      = false;
  globval.radiation      = false;
  globval.emittance      = false;
  globval.IBS            = false;
  globval.pathlength     = false;
  globval.Aperture_on    = false;
  globval.Cart_Bend      = false;
  globval.dip_edge_fudge = true;
}


int main(int argc, char *argv[])
{
  const long seed = 1121;

  int              loc;
  double           I[6], eps_x, sigma_delta, U_0, J[3], tau[3];
  std::vector<int> bpm;

  iniranf(seed);
  setrancut(1.0);

  reverse_elem = false;

  globval.mat_meth = mat_meth;

  trace = false;

  FieldMap_filetype = 6;

  if (!true)
    Read_Lattice(argv[1]);
  else {
#if 0
    rdmfile(argv[1]);
#else
    rdmfile_at(argv[1]);
#endif
  }

  set_state();

#if 0
  long int        lastpos;
  ss_vect<double> ps;
  ss_vect<tps>    M;

  globval.radiation = true;

  ps[x_]     =  1e-6;
  ps[px_]    = -2e-6;
  ps[y_]     =  3e-6;
  ps[py_]    = -4e-6;
  ps[delta_] =  5e-6;
  ps[ct_]    = -6e-6;
  printf("\nx_0:\n");
  cout << scientific << setprecision(5) << setw(13) << ps << "\n";
  Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
  printf("\nx_1:\n");
  cout << scientific << setprecision(5) << setw(13) << ps << "\n";

  M.identity();
  printf("\nM:\n");
  prt_lin_map(3, M);
  Cell_Pass(0, globval.Cell_nLoc, M, lastpos);
  printf("\nM:\n");
  prt_lin_map(3, M);
  assert(false);

  if (false) {
    long lastpos;
    getcod(1e-3, lastpos);
    prt_cod_1("cod.out");
    exit(0);
  }
#endif

  chk_phi();

  if (phiob_2xL) get_phiob_2xL_ratios();

  if (zero_b_3)
    no_mult(Sext);
  if (zero_b_4)
    no_mult(Oct);

  Ring_GetTwiss(true, 0e0);
  printglob();

  if (!false)
    compute_alpha_bucket();

  if (false)
    compute_phi(4);

  if (!false) prt_b_n();

  if (false) {
    // A 1/2 ps_rot at the entrance & exit of the super period for a symmetric
    // approach.
    loc = Elem_GetPos(ElemIndex("lsborder"), 2);
    printf("\nloc = %d\n", loc);
    set_dnu_straight("ps_rot", loc);

    Ring_GetTwiss(true, 0e0);
    printglob();
  }

  if (fit_nu)
    fit_nu_jb_2(nu[X_]-globval.TotalTune[X_],
		nu[Y_]-globval.TotalTune[Y_]);

  if (fit_xi) {
    if (true)
      fit_xi_jb_2(0e0, 0e0);
    else
      fit_xi_jb_2(2e0/20e0, 2e0/20e0);
  }

  if (ps_rot) {
    // A 1/2 ps_rot at the entrance & exit of the super period for a symmetric
    // approach.
    set_ps_rot("ps_rot", dnu[X_]/2e0, dnu[Y_]/2e0);
  }

  prtmfile("flat_file.dat");
  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);
  prt_chrom_lat("chromlat.out");

  if (Deta)
    compute_Deta(2e-2);

  if (chk_mpole_sym)
    chk_mpole(2);

  if (chk_dnu)
    chk_dnu_straight("lsborder");

  if (!false) {
    if (!globval.mat_meth)
      GetEmittance(ElemIndex("cav"), false, true);
    else
      get_eps_x(eps_x, sigma_delta, U_0, J, tau, I, true);

    if (comp_H_long)
      prt_H_long(ElemIndex("cav"), 25, 180e0, 15e-2, 3, false);
  }

  if (false) {
    loc = Elem_GetPos(ElemIndex("lsborder"), 2);
    printf("\nloc = %d\n", loc);
    get_dnu_straight(loc);
  }

  if (false) {
    chk_optics(0.0, 13.24261, 0.0, 0.0, 0.0, 2.35728, 0.0, 0.0);
    prt_lat("linlat1.out", globval.bpm, true);
    prt_lat("linlat.out", globval.bpm, true, 10);
    prtmfile("flat_file.dat");
  }

  if (false) {
    globval.Cavity_on = false;
    track(100, -6e-3, 0e0);
  }

  if (get_tol) {
    if (true) {
      printf("\n");
      get_b_2_tol(0.00e-3, n_aper, n_track);
      get_b_2_tol(0.25e-3, n_aper, n_track);
      get_b_2_tol(0.50e-3, n_aper, n_track);
      get_b_2_tol(1.00e-3, n_aper, n_track);
      get_b_2_tol(2.50e-3, n_aper, n_track);
    }

    Ring_GetTwiss(true, 0e0);
    prt_lat("linlat1.out", globval.bpm, true);
    prt_lat("linlat.out", globval.bpm, true, 10);
  }
}
