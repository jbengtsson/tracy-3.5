#define NO 1

#include "tracy_lib.h"


int no_tps = NO;


const bool
  zero_b_3      = !false,
  zero_b_4      = !false,
  set_b_3       = false,
  ps_rot        = false,
  chk_mpole_sym = false,
  chk_dnu       = false, // Requires super period.
  comp_H_long   = false;

const double
  dnu[] = {0.0, 0.0};


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


void fit_xi_jb
(const std::vector<int> &Fnum_b3, const double ksi_x, const double ksi_y,
 const double db3L)
{
  int    n_b3, j, k;
  double **A, **U, **V, *w, *b, *x, b3, a3;

  const bool   prt = !false;
  const int    m   = 2;
  const double
    ksi0[]  = {ksi_x, ksi_y},
    svd_cut = 1e-10;

  n_b3 = Fnum_b3.size();

  A = dmatrix(1, m, 1, n_b3); U = dmatrix(1, m, 1, n_b3);
  V = dmatrix(1, n_b3, 1, n_b3);
  w = dvector(1, n_b3); b = dvector(1, m); x = dvector(1, n_b3);

  // Zero sextupoles to track linear chromaticity.
  if (false) no_sxt();

  for (k = 1; k <= n_b3; k++) {
    set_dbnL_design_fam(Fnum_b3[k-1], Sext, db3L, 0e0);
    Ring_Getchrom(0e0);
    if (prt)
      printf("\nfit_xi_jb: ksi1+ = [%9.5f, %9.5f]\n",
	     globval.Chrom[X_], globval.Chrom[Y_]);

    for (j = 1; j <= m; j++)
      A[j][k] = globval.Chrom[j-1];
    set_dbnL_design_fam(Fnum_b3[k-1], Sext, -2e0*db3L, 0e0);
    Ring_Getchrom(0e0);
    if (prt)
      printf("fit_xi_jb: ksi1- = [%9.5f, %9.5f]\n",
	 globval.Chrom[X_], globval.Chrom[Y_]);
    for (j = 1; j <= 2; j++) {
      A[j][k] -= globval.Chrom[j-1];
      A[j][k] /= 2e0*db3L;
    }

    set_dbnL_design_fam(Fnum_b3[k-1], Sext, db3L, 0e0);
  }

  Ring_Getchrom(0e0);
  if (prt)
    printf("\nfit_xi_jb: ksi1  = [%9.5f, %9.5f]\n",
	   globval.Chrom[X_], globval.Chrom[Y_]);
  for (j = 1; j <= 2; j++)
    b[j] = -(globval.Chrom[j-1]-ksi0[j-1]);

  dmcopy(A, m, n_b3, U);
  dsvdcmp(U, m, n_b3, w, V);

  if (prt) {
    printf("\nfit_xi_jb:\n  singular values:\n");
    for (j = 1; j <= n_b3; j++) {
      printf("    %9.3e", w[j]);
      if (w[j] < svd_cut) {
	w[j] = 0e0;
	printf(" (zeroed)");
      }
      printf("\n");
    }
  }

  dsvbksb(U, w, V, m, n_b3, b, x);

  if (prt) {
    dmdump(stdout, "\nA:", A, 2, 2, "%11.3e");
    dvdump(stdout, "\nx:", x, 2, "%11.3e");
  }

  for (k = 1; k <= n_b3; k++)
    set_dbnL_design_fam(Fnum_b3[k-1], Sext, x[k], 0e0);

  if (prt) {
    printf("\n  b3:\n");
    for (k = 0; k < n_b3; k++) {
      get_bn_design_elem(Fnum_b3[k], 1, Sext, b3, a3);
      printf("    %-8s %10.5f\n", ElemFam[Fnum_b3[k]-1].ElemF.PName, b3);
    }
    printf("\n");
  }

  free_dmatrix(A, 1, m, 1, n_b3); free_dmatrix(U, 1, m, 1, n_b3);
  free_dmatrix(V, 1, n_b3, 1, n_b3);
  free_dvector(w, 1, n_b3); free_dvector(b, 1, m); free_dvector(x, 1, n_b3);
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


void chk_mpole_Fam(const int Fnum)
{
  int n_Kids, k, loc[2];

  printf("\n   name        s     beta_x   beta_y   eta_x    dnu_x    dnu_y\n");
  n_Kids = GetnKid(Fnum);
  for (k = 1; k <= n_Kids; k++) {
    loc[0] = (k > 1)? Elem_GetPos(Fnum, k-1) : Elem_GetPos(Fnum, n_Kids);
    loc[1] = Elem_GetPos(Fnum, k);
    printf("  %.8s %7.3f %8.5f %8.5f %8.5f %8.5f %8.5f\n",
	   Cell[loc[1]].Elem.PName, Cell[loc[1]].S,
	   Cell[loc[1]].Beta[X_], Cell[loc[1]].Beta[Y_], Cell[loc[1]].Eta[X_],
	   (k > 1)? Cell[loc[1]].Nu[X_]-Cell[loc[0]].Nu[X_] : NAN,
	   (k > 1)?Cell[loc[1]].Nu[Y_]-Cell[loc[0]].Nu[Y_] : NAN);
  }
}


void chk_mpole(const int lat_case)
{
  int              k;
  std::vector<int> Fnum;

  switch (lat_case) {
  case 1:
    Fnum.push_back(ElemIndex("s1_h2"));
    Fnum.push_back(ElemIndex("s2_h2"));
    Fnum.push_back(ElemIndex("s3_h2"));
    Fnum.push_back(ElemIndex("s4_h2"));
    Fnum.push_back(ElemIndex("s5_h2"));
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
  printf("  RF alpha bucket [%%] = [%3.1f, %3.1f]\n",
	 1e2*alpha_c[1]/(2e0*alpha_c[2]), -1e2*alpha_c[1]/alpha_c[2]);
}


double H_long
(const double phi, const double delta, const int h_rf, const double V_rf,
 const double phi_0, const psVector &alpha_c)
{
  const int
    n_alpha_c = 3;
  const double
    E_0 = 1e9*globval.Energy;

  double H;

  H = V_rf/E_0*(cos(phi+phi_0)+phi*sin(phi_0));
  for (auto i = 2; i <= n_alpha_c+1; i++)
    H += 2e0*pi*h_rf*alpha_c[i-1]*pow(delta, (double)i)/i;
  return H;
}


void prt_H_long
(const string &cav_name, const int n, const double phi_max,
 const double delta_max, const bool neg_alpha_c)
{
  const string
    file_name = "H_long.dat";
  const long int
    loc = Elem_GetPos(ElemIndex(cav_name.c_str()), 1);
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
      H = H_long(phi, delta, h_RF, V_RF, M_PI+phi_0, alpha_c);
      fprintf(outf, "  %8.2f %10.5f, %13.5e\n", phi*180e0/M_PI, 1e2*delta, H);
    }
    fprintf(outf, "\n");
  }

  fclose(outf);
}


void compute_Deta_x(const double delta)
{
  // Evaluate derivative; to avoid effect of tune shift.
  double         h;
  vector<double> Deta_x;
  FILE           *outf;

  const double d_delta   = 1e-5;
  const string file_name = "Deta_x.out";

  outf = file_write(file_name.c_str());

  printf("\nOptics for delta = %10.3e\n", d_delta);
  Ring_GetTwiss(true, d_delta);
  printglob();
  for (auto k = 0; k <= globval.Cell_nLoc; k++)
    Deta_x.push_back(Cell[k].Eta[X_]);
  printf("\nOptics for delta = %10.3e\n", -d_delta);
  Ring_GetTwiss(true, -d_delta);
  printglob();
  fprintf(outf, "#  k     name             s    type    eta_x      eta'_x"
	        "    Ddeta_x/Ddelta\n"
	        "#                        [m]            [m]"
	        "                      [m]\n");
  for (auto k = 0; k <= globval.Cell_nLoc; k++) {
    if (Cell[k].Elem.Pkind == Mpole)
      h = Cell[k].Elem.M->Pirho;
    else
      h = 0e0;
    Deta_x[k] -= Cell[k].Eta[X_];
    Deta_x[k] /= (2e0*d_delta);
    fprintf(outf, "%4d %10s %8.3f %4.1f %12.5e %12.5e %12.5e\n",
	    k, Cell[k].Elem.PName, Cell[k].S, get_code(Cell[k]),
	    h*Cell[k].Eta[x_], pow(Cell[k].Etap[x_], 2)/2e0, h*Deta_x[k]);
  }

  fclose(outf);
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
  int              loc;
  double           I[6], eps_x, sigma_delta, U_0, J[3], tau[3];
  std::vector<int> bpm;

  reverse_elem = !false;

  globval.mat_meth = false;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  set_state();

  chk_phi();

  if (zero_b_3)
    no_mult(Sext);
  if (zero_b_4)
    no_mult(Oct);

  Ring_GetTwiss(true, 0e0);
  printglob();

  if (!false)
    compute_alpha_bucket();

  prt_b_n();

  if (false) {
    // A 1/2 ps_rot at the entrance & exit of the super period for a symmetric
    // approach.
    loc = Elem_GetPos(ElemIndex("lsborder"), 2);
    printf("\nloc = %d\n", loc);
    set_dnu_straight("ps_rot", loc);

    Ring_GetTwiss(true, 0e0);
    printglob();
  }

  if (set_b_3) {
    const int lat = 3;
    std::vector<int> Fnum;
    switch (lat) {
    case 1:
      Fnum.push_back(ElemIndex("s3"));
      Fnum.push_back(ElemIndex("s4"));
      break;
    case 2:
      Fnum.push_back(ElemIndex("s3_f1"));
      Fnum.push_back(ElemIndex("s4_f1"));
      break;
    case 3:
      Fnum.push_back(ElemIndex("s4_h2"));
      Fnum.push_back(ElemIndex("s5_h2"));
      break;
    }
    fit_xi_jb(Fnum, 0e0, 0e0, 1e0);

    Ring_GetTwiss(true, 0e0);
    printglob();
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

  if (!false)
    compute_Deta_x(2e-2);

  if (chk_mpole_sym)
    chk_mpole(1);

  if (chk_dnu)
    chk_dnu_straight("lsborder");

  if (!false) {
    if (!globval.mat_meth)
      GetEmittance(ElemIndex("cav"), false, true);
    else
      get_eps_x(eps_x, sigma_delta, U_0, J, tau, I, true);

    if (comp_H_long)
      prt_H_long("cav", 25, 180e0, 10e-2, false);
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
}
