#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


const bool
  set_dnu = false,
  mI_rot  = false,
  HOA_rot = false,
  prt_ms  = false,
  prt_dt  = false;

const int n_cell = 4;

const double
  // nu[]     = {62.74, 21.34},
  nu[]     = {0.0, -0.2},
  dnu_mI[] = {1.5-1.44129-0.0, 0.5-0.47593-0.0},
  nu_HOA[] = {19.0/8.0, 15.0/16.0};


double rad2deg(const double a) { return a*180e0/M_PI; }

double deg2rad(const double a) { return a*M_PI/180e0; }


void prt_name(FILE *outf, const char *name, const string &str, const int len)
{
  int j, k;

  j = 0;
  do {
    fprintf(outf, "%c", name[j]);
    j++;
  } while (name[j] != ' ');
  fprintf(outf, str.c_str());
  for (k = j; k < len; k++)
    fprintf(outf, " ");
}


void get_cod_rms(const double dx, const double dy,
		 const int n_seed, const bool all)
{
  bool                cod;
  int                 i, j, k, n, n_cod;
  std::vector<double> x1[6], x2[6], x_mean[6], x_sigma[6];
  FILE                *fp;

  const int n_cod_corr = 5;

  globval.Cavity_on = false;

  for (j = 0; j <= globval.Cell_nLoc; j++)
    for (k = 0; k < 6; k++) {
      x1[k].push_back(0e0); x2[k].push_back(0e0);
    }
  
  fp = file_write("cod_rms.out");
  
  n_cod = 0;
  for (i = 0; i < n_seed; i++) {
    printf("\norb_corr: seed no %d\n", i+1);

    misalign_rms_type(Dip,  dx, dy, 0e0, true);
    misalign_rms_type(Quad, dx, dy, 0e0, true);
    
    cod = orb_corr(n_cod_corr);

    if (cod) {
      n_cod++;

      n = 0;
      for (j = 0; j <= globval.Cell_nLoc; j++)
	if (all || ((Cell[j].Elem.Pkind == Mpole) &&
		    (Cell[j].Elem.M->n_design == Sext))) {
	  n++;
	  for (k = 0; k < 6; k++) {
	    x1[k][n-1] += Cell[j].BeamPos[k];
	    x2[k][n-1] += sqr(Cell[j].BeamPos[k]);
	  }
	}
    } else
      printf("orb_corr: failed\n");

    // Reset orbit trims.
    set_bn_design_fam(globval.hcorr, Dip, 0e0, 0e0);
    set_bn_design_fam(globval.vcorr, Dip, 0e0, 0e0);
  }

  printf("\nget_cod_rms: no of seeds %d, no of cods %d\n", n_seed, n_cod);

  n = 0;
  for (j = 0; j <= globval.Cell_nLoc; j++)
    if (all || ((Cell[j].Elem.Pkind == Mpole) &&
		(Cell[j].Elem.M->n_design == Sext))) {
      n++;
      for (k = 0; k < 6; k++) {
	x_mean[k].push_back(x1[k][n-1]/n_cod);
	x_sigma[k].push_back(sqrt((n_cod*x2[k][n-1]-sqr(x1[k][n-1]))
				  /(n_cod*(n_cod-1.0))));
      }
      fprintf(fp, "%8.3f %6.2f %10.3e +/- %10.3e %10.3e +/- %10.3e\n",
	      Cell[j].S, get_code(Cell[j]),
	      1e3*x_mean[x_][n-1], 1e3*x_sigma[x_][n-1],
	      1e3*x_mean[y_][n-1], 1e3*x_sigma[y_][n-1]);
    } else
      fprintf(fp, "%8.3f %6.2f\n", Cell[j].S, get_code(Cell[j]));
  
  fclose(fp);
}


void track(const double Ax, const double Ay)
{
  long int        lastpos;
  int             i;
  ss_vect<double> xt, xs;
  FILE            *fd;

  getcod(0e0, lastpos);

  fd = fopen("trackdat_oneturn.dat","w");
  fprintf(fd, "orbit %22.14e %22.14e %22.14e %22.14e %22.14e %22.14e\n",
	  Cell[0].BeamPos[0], Cell[0].BeamPos[1],
	  Cell[0].BeamPos[2], Cell[0].BeamPos[3],
	  Cell[0].BeamPos[4], Cell[0].BeamPos[5] );
  fprintf(fd, "orbit %22.14e %22.14e %22.14e %22.14e %22.14e %22.14e\n",
	  globval.CODvect[0], globval.CODvect[1],
	  globval.CODvect[2], globval.CODvect[3],
	  globval.CODvect[4], globval.CODvect[5]);

  xt.zero(); xt[x_] = Ax; xt[y_] = Ay; 

  fprintf(fd, "start %22.14e %22.14e %22.14e %22.14e %22.14e %22.14e\n",
	  xt[0], xt[1], xt[2], xt[3], xt[4], xt[5] );

  for (i = 0; i <= globval.Cell_nLoc; i++) {
    Cell_Pass(i, i, xt, lastpos);
    fprintf(fd, "%5d %22.14e %22.14e %22.14e %22.14e %22.14e %22.14e \n",
	    i, xt[0], xt[1], xt[2], xt[3], xt[4], xt[5]);
  }

  fclose(fd);
}


void fit_ksi1(const std::vector<int> &Fnum_b3,
	      const double ksi_x, const double ksi_y, const double db3L,
	      double svd[])
{
  int    n_b3, j, k, n_svd;
  double **A, **U, **V, *w, *b, *x, b3L, a3L;

  const bool   prt = false;
  const int    m   = 2;
  const double
    ksi0[]  = {ksi_x, ksi_y},
    svd_cut = 1e-10;

  n_b3 = Fnum_b3.size();

  A = dmatrix(1, m, 1, n_b3); U = dmatrix(1, m, 1, n_b3);
  V = dmatrix(1, n_b3, 1, n_b3);
  w = dvector(1, n_b3); b = dvector(1, m); x = dvector(1, n_b3);

  // Zero sextupoles to track linear chromaticity.
  no_sxt();

  for (k = 1; k <= n_b3; k++) {
    set_dbnL_design_fam(Fnum_b3[k-1], Sext, db3L, 0e0);
    Ring_Getchrom(0e0);
    if (prt)
      printf("\nfit_ksi1: ksi1+ = [%9.5f, %9.5f]\n",
	     globval.Chrom[X_], globval.Chrom[Y_]);

    for (j = 1; j <= m; j++)
      A[j][k] = globval.Chrom[j-1];
    set_dbnL_design_fam(Fnum_b3[k-1], Sext, -2e0*db3L, 0e0);
    Ring_Getchrom(0e0);
    if (prt)
      printf("fit_ksi1: ksi1- = [%9.5f, %9.5f]\n",
	 globval.Chrom[X_], globval.Chrom[Y_]);
    for (j = 1; j <= 2; j++) {
      A[j][k] -= globval.Chrom[j-1]; A[j][k] /= 2e0*db3L;
    }

    set_dbnL_design_fam(Fnum_b3[k-1], Sext, db3L, 0e0);
  }

  Ring_Getchrom(0e0);
  if (prt)
    printf("\nfit_ksi1: ksi1  = [%9.5f, %9.5f]\n",
	   globval.Chrom[X_], globval.Chrom[Y_]);
  for (j = 1; j <= 2; j++)
    b[j] = -(globval.Chrom[j-1]-ksi0[j-1]);

  dmcopy(A, m, n_b3, U); dsvdcmp(U, m, n_b3, w, V);

  printf("\nfit_ksi1:\n  singular values:\n");
  n_svd = 0;
  for (j = 1; j <= n_b3; j++) {
    printf("  %9.3e", w[j]);
    if (w[j] < svd_cut) {
      w[j] = 0e0;
      printf(" (zeroed)");
    } else {
      if (n_svd > 2) {
	printf("fit_ksi1: more than 2 non-zero singular values");
	exit(1);
      }
      svd[n_svd] = w[j];
      n_svd++;
    }
    printf("\n");
  }

  dsvbksb(U, w, V, m, n_b3, b, x);

  for (k = 1; k <= n_b3; k++)
    set_dbnL_design_fam(Fnum_b3[k-1], Sext, x[k], 0e0);

  if (prt) {
    printf("  b3:\n  ");
    for (k = 0; k < n_b3; k++) {
      get_bn_design_elem(Fnum_b3[k], 1, Sext, b3L, a3L);
      printf(" %9.5f", b3L);
    }
    printf("\n");
  }

  free_dmatrix(A, 1, m, 1, n_b3); free_dmatrix(U, 1, m, 1, n_b3);
  free_dmatrix(V, 1, n_b3, 1, n_b3);
  free_dvector(w, 1, n_b3); free_dvector(b, 1, m); free_dvector(x, 1, n_b3);
}


void fit_ksi1(const int lat_case, const double ksi_x, const double ksi_y)
{
  double           svd[2];
  std::vector<int> Fnum;

  // ESRF-U        1,
  // M-6HBAi-2-1-1 2,
  // M-6HBA-0-.-.  3.

  switch (lat_case) {
  case 1:
    Fnum.push_back(ElemIndex("sd1a"));
    Fnum.push_back(ElemIndex("sd1b"));
    Fnum.push_back(ElemIndex("sd1d"));
    Fnum.push_back(ElemIndex("sd1e"));
    Fnum.push_back(ElemIndex("sf2a"));
    Fnum.push_back(ElemIndex("sf2e"));
    Fnum.push_back(ElemIndex("sh1a"));
    Fnum.push_back(ElemIndex("sh2b"));
    Fnum.push_back(ElemIndex("sh3e"));
    break;
  case 2:
    Fnum.push_back(ElemIndex("sf1"));
    Fnum.push_back(ElemIndex("sd1"));
    Fnum.push_back(ElemIndex("sd2"));
    Fnum.push_back(ElemIndex("sh1"));
    Fnum.push_back(ElemIndex("sh2"));
    break;
  case 3:
    Fnum.push_back(ElemIndex("sf1"));
    Fnum.push_back(ElemIndex("sd1"));
    Fnum.push_back(ElemIndex("sd2"));
    Fnum.push_back(ElemIndex("sh1"));
    Fnum.push_back(ElemIndex("sh2"));
    break;
  default:
    printf("\nfit_ksi1: unknown lattice type\n");
    exit(1);
    break;
  }

  fit_ksi1(Fnum, 0e0, 0e0, 1e1, svd);
}


void chk_phi()
{
  int    k;
  double dphi, phi, mphi;

  printf("\n");
  phi = 0e0; mphi = 0e0;
  for (k = 0; k <= globval.Cell_nLoc; k++) {
    // if ((Cell[k].Elem.Pkind == Mpole) &&
    // 	(Cell[k].Elem.M->n_design == Dip)) {
    if ((Cell[k].Elem.Pkind == Mpole) &&
	(Cell[k].Elem.M->Pirho != 0e0)) {
      dphi = rad2deg(Cell[k].Elem.PL*Cell[k].Elem.M->Pirho);
      if (dphi != 0e0) {
	prt_name(stdout, Cell[k].Elem.PName, "", 8);
	printf(" %9.6f\n", dphi);
      }
      phi += dphi;
      if (dphi < 0e0) mphi += dphi;
    }
  }
  printf("\nphi = %8.6f phi- = %8.6f phi+ = %8.6f\n", phi, mphi, phi-mphi);
}


void dpath_length()
{
  int    k, loc;
  double phi, dL, L_tot, phi_tot, mphi, Lc1, phi2, rho2, L2, Lc2;

  printf("\n");
  L_tot = 0e0; phi_tot = 0e0; mphi = 0e0;
  for (k = 0; k <= globval.Cell_nLoc; k++) {
    if ((Cell[k].Elem.Pkind == Mpole) &&
	(Cell[k].Elem.M->Pirho != 0e0)) {
      phi = Cell[k].Elem.PL*Cell[k].Elem.M->Pirho;
      phi_tot += phi;
      if (phi < 0e0) {
	dL = Cell[k].Elem.PL - 2e0*sin(phi/2e0)/Cell[k].Elem.M->Pirho;
	L_tot += dL; mphi += phi;
	prt_name(stdout, Cell[k].Elem.PName, "", 8);
	printf(" phi [deg] = %9.6f L [m] = %9.6f rho [m] = %9.6f"
	       " dL [mm] = %9.6f L_tot [mm] = %9.6f\n",
	       rad2deg(phi), Cell[k].Elem.PL, 1e0/Cell[k].Elem.M->Pirho,
	       1e3*dL, 1e3*L_tot);
      }
    }
  }
  printf("\nphi = %8.6f phi- = %8.6f phi+ = %8.6f\n",
	 rad2deg(phi_tot), rad2deg(mphi), rad2deg(phi_tot-mphi));

  printf("\n");
  loc = Elem_GetPos(ElemIndex("dq1"), 1);
  phi = Cell[loc].Elem.PL*Cell[loc].Elem.M->Pirho;
  Lc1 = 2e0*sin(phi/2e0)/Cell[loc].Elem.M->Pirho;
  phi2 = phi + mphi/8e0; rho2 = Lc1/(2e0*sin(phi2/2e0)); L2 = rho2*phi2;
  Lc1 = 2e0*sin(phi/2e0)/Cell[loc].Elem.M->Pirho;
  Lc2 = 2e0*rho2*sin(phi2/2e0);
  prt_name(stdout, Cell[loc].Elem.PName, "", 8);
  printf(" phi = %9.6f L = %9.6f rho = %9.6f Lc1 = %9.6f\n",
	 rad2deg(phi), Cell[loc].Elem.PL, 1e0/Cell[loc].Elem.M->Pirho, Lc1);
  printf("         phi = %9.6f L = %9.6f rho = %9.6f Lc2 = %9.6f\n",
	 rad2deg(phi2), L2, rho2, Lc2);
  printf(" dL [mm] = %9.6f\n", 1e3*8e0*(L2-Cell[loc].Elem.PL));
}


void prt_symm(const std::vector<int> &Fam)
{
  long int loc;
  int      j, k;

  for (j = 0; j < (int)Fam.size(); j++) {
    printf("\n");
    for (k = 1; k <= GetnKid(Fam[j]); k++) {
      loc = Elem_GetPos(Fam[j], k);
      if (k % 2 == 0) loc -= 1;
      printf(" %5.1f %6.3f %6.3f %6.3f\n",
	     Cell[loc].S, Cell[loc].Beta[X_], Cell[loc].Beta[Y_],
	     Cell[loc].Eta[X_]);
    }
  }
}


void prt_quad(const std::vector<int> &Fam)
{
  long int loc;
  int      j;

  printf("\n");
  for (j = 0; j < (int)Fam.size(); j++) {
    loc = Elem_GetPos(Fam[j], 1);
    printf(" %4.1f %6.3f %6.3f %2d\n",
	   Cell[loc].S, Cell[loc].Beta[X_], Cell[loc].Beta[Y_],
	   GetnKid(Fam[j]));
  }
}


void prt_drift()
{
  int k;

  printf("\n");
  for (k = 0; k <= globval.Cell_nLoc; k++)
    if (Cell[k].Elem.Pkind == drift)
      printf("%3d %10s %13.10f\n", k, Cell[k].Elem.PName, Cell[k].Elem.PL);
}


void prt_dip()
{
  int                        j, k, loc;
  double                     phi, L, L1, phi1, phi_rel, phi_rel_tot;
  std::vector<int>           row;
  std::vector< vector<int> > Fnum;
  std::vector<double>        phi2;

  row.push_back(ElemIndex("bl1_1"));
  row.push_back(ElemIndex("bl1_2"));
  row.push_back(ElemIndex("bl1_3"));
  row.push_back(ElemIndex("bl1_4"));
  row.push_back(ElemIndex("bl1_5"));
  Fnum.push_back(row);
  row.clear();

  row.push_back(ElemIndex("bl2_1"));
  row.push_back(ElemIndex("bl2_2"));
  row.push_back(ElemIndex("bl2_3"));
  row.push_back(ElemIndex("bl2_4"));
  row.push_back(ElemIndex("bl2_5"));
  Fnum.push_back(row);
  row.clear();

  for (j = 0; j < (int)Fnum.size(); j++) {
    printf("\n");
    L1 = 0.0; phi1 = 0e0; phi2.push_back(0e0);
    for (k = 0; k < (int)Fnum[j].size(); k++) {
      loc = Elem_GetPos(Fnum[j][k], 1);
      L = Cell[loc].Elem.PL; phi = rad2deg(L*Cell[loc].Elem.M->Pirho);
      L1 += Cell[loc].Elem.PL; phi1 += phi;
      printf("%10s %13.10f %13.10f %13.10f %13.10f\n",
	     Cell[loc].Elem.PName, L,
	     phi, Cell[loc].Elem.M->PTx1, Cell[loc].Elem.M->PTx2);
    }
    printf("\nMagnet: L = %13.10f phi = %13.10f\n", L1, phi1);
    phi2[j] += phi1;
    printf("\nCell: phi = %13.10f\n", phi2[j]);
  }

  for (j = 0; j < (int)Fnum.size(); j++) {
    printf("\nphi ratios: \n");
    phi_rel_tot = 0e0;
    for (k = 0; k < (int)Fnum[j].size(); k++) {
      loc = Elem_GetPos(Fnum[j][k], 1);
      L = Cell[loc].Elem.PL; phi = rad2deg(L*Cell[loc].Elem.M->Pirho);
      phi_rel = phi/phi2[j];
      phi_rel_tot += phi_rel;
      printf(" %8.6f", phi_rel);
    }
    printf("\nTotal: %8.6f\n", phi_rel_tot);
  }
}


void chk_optics(const double alpha_x, const double beta_x,
		const double alpha_y, const double beta_y,
		const double eta_x, const double etap_x,
		const double eta_y, const double etap_y)
{
  Vector2 alpha, beta, eta, etap;

  alpha[X_] = alpha_x; alpha[Y_] = alpha_y;
  beta[X_]  = beta_x;  beta[Y_]  = beta_y;
  eta[X_]   = eta_x;   eta[Y_]   = eta_y;
  etap[X_]  = etap_x;  etap[Y_]  = etap_y;

  ttwiss(alpha, beta, eta, etap, 0e0);
}


void chk_mini_beta(const std::vector<int> &Fam)
{
  int    j, k, loc;
  double nu0[] = {0e0, 0e0};

  for (j = 0; j < (int)Fam.size(); j++) {
    printf("\n");
    for (k = 1; k <= GetnKid(Fam[j]); k++) {
      loc = Elem_GetPos(Fam[j], k);
      if (k % 2 == 1) {
	nu0[X_] = Cell[loc].Nu[X_]; nu0[Y_] = Cell[loc].Nu[Y_];
      } else {
	loc -= 1;
	printf(" %5.1f %6.3f %6.3f %6.3f %8.5f %8.5f\n",
	       Cell[loc].S, Cell[loc].Beta[X_], Cell[loc].Beta[Y_],
	       Cell[loc].Eta[X_],
	       Cell[loc].Nu[X_]-nu0[X_], Cell[loc].Nu[Y_]-nu0[Y_]);
      }
    }
  }
}


void chk_high_ord_achr(const int lat_case)
{
  int              k;
  double           dnu[2];
  std::vector<int> loc;

  // ESRF-U        1,
  // M-6HBAi-2-1-1 2,
  // M-6HBA-0-.-.  3.

  Ring_GetTwiss(true, 0e0);
 
  switch (lat_case) {
  case 1:
    dnu[X_] = 19.0/8.0; dnu[Y_] = 15.0/16.0;
    // loc.push_back(Elem_GetPos(ElemIndex("dc_1_01"),  1));
    loc.push_back(0);
    loc.push_back(Elem_GetPos(ElemIndex("idmarker"), 2));
    loc.push_back(Elem_GetPos(ElemIndex("idmarker"), 3));
    loc.push_back(Elem_GetPos(ElemIndex("idmarker"), 4));
    loc.push_back(Elem_GetPos(ElemIndex("idmarker"), 5));
   break;
  case 2:
    dnu[X_] = 11.0/8.0; dnu[Y_] = 15.0/16;
    loc.push_back(Elem_GetPos(ElemIndex("ls"), 1));
    loc.push_back(Elem_GetPos(ElemIndex("ss"), 1));
    loc.push_back(Elem_GetPos(ElemIndex("ss"), 2));
    loc.push_back(Elem_GetPos(ElemIndex("ss"), 3));
    loc.push_back(Elem_GetPos(ElemIndex("ls"), 2));
    break;
  default:
    printf("\nchk_high_ord_achr: unknown lattice type\n");
    exit(1);
    break;
  }

  printf("\nCell phase advance:\n");
  printf("Ideal:    [%7.5f, %7.5f]\n", dnu[X_], dnu[Y_]);
  for (k = 1; k < (int)loc.size(); k++)
    printf(" %9.5f %8.5f %8.5f %7.5f [%7.5f, %7.5f]\n",
	   Cell[loc[k]].S, Cell[loc[k]].Alpha[X_], Cell[loc[k]].Alpha[Y_],
	   Cell[loc[k]].S-Cell[loc[k-1]].S, 
	   Cell[loc[k]].Nu[X_]-Cell[loc[k-1]].Nu[X_], 
	   Cell[loc[k]].Nu[Y_]-Cell[loc[k-1]].Nu[Y_]);
}


void chk_mI_trans(const int lat_case)
{
  int Fnum, k, loc0, loc1;

  // ESRF-U        1,
  // M-6HBAi-2-1-1 2,
  // M-6HBA-0-.-.  3.

  Ring_GetTwiss(true, 0e0);
 
  switch (lat_case) {
  case 1:
    Fnum = ElemIndex("dispbumpcenter");
   break;
  case 2:
    Fnum = ElemIndex("sf1");
   break;
  case 3:
    Fnum = ElemIndex("sf1_ctr");
   break;
  default:
    printf("\nchk_mI_trans: unknown lattice type\n");
    exit(1);
    break;
  }

  printf("\nChromatic sextupole phase advance:\n");
  // for (k = 3; k <= GetnKid(Fnum); k += 4) {
  for (k = 2; k <= GetnKid(Fnum); k += 2) {
    loc0 = Elem_GetPos(Fnum, k-1); loc1 = Elem_GetPos(Fnum, k);
    printf(" %8s %7.3f [%7.5f, %7.5f]\n",
	   Cell[loc1].Elem.PName, Cell[loc1].S,
	   Cell[loc1].Nu[X_]-Cell[loc0].Nu[X_], 
	   Cell[loc1].Nu[Y_]-Cell[loc0].Nu[Y_]);
  }
}


void chk_lin_chrom(void)
{
  int Fnum, loc0, loc1;

  Ring_GetTwiss(true, 0e0);
 
  printf("\nchk_lin_chrom:\n");
  Fnum = ElemIndex("sf1");
  loc0 = Elem_GetPos(Fnum, 2); loc1 = Elem_GetPos(Fnum, 3);
  printf(" %8s [%7.5f, %7.5f]\n",
	 Cell[loc1].Elem.PName,
	 Cell[loc1].Nu[X_]-Cell[loc0].Nu[X_], 
	 Cell[loc1].Nu[Y_]-Cell[loc0].Nu[Y_]);
  loc1 = Elem_GetPos(Fnum, 1);
  printf(" %8s [%7.5f, %7.5f]\n",
	 Cell[loc1].Elem.PName, 2e0*Cell[loc1].Nu[X_], 2e0*Cell[loc1].Nu[Y_]);
}


void chk_drv_terms(void)
{
  int k;

  Ring_GetTwiss(true, 0e0);
 
  printf("\nh_10200 phase advance:\n");
  for (k = 0; k <= globval.Cell_nLoc; k++)
    if ((Cell[k].Elem.Pkind == Mpole) && (Cell[k].Elem.M->Porder == Sext))
      printf(" %8s %7.3f %7.5f\n",
	     Cell[k].Elem.PName, Cell[k].S, Cell[k].Nu[X_]+2.0*Cell[k].Nu[Y_]);
  k = globval.Cell_nLoc;
  printf(" %8s %7.3f %7.5f\n",
	 Cell[k].Elem.PName, Cell[k].S, Cell[k].Nu[X_]+2.0*Cell[k].Nu[Y_]);
}


void chk_mpole_Fam(const int Fnum, const bool exit)
{
  int k, loc;

  printf("\n");
  for (k = 1; k <= GetnKid(Fnum); k++) {
    loc = Elem_GetPos(Fnum, k);
    if (!exit && ((k-1) % 2 == 1)) loc -= 1;
    printf("%8s %7.3f %8.5f %8.5f %8.5f\n",
	   Cell[loc].Elem.PName, Cell[loc].S,
	   Cell[loc].Beta[X_], Cell[loc].Beta[Y_], Cell[loc].Eta[X_]);
  }
}


void chk_mpole(void)
{
  int              k;
  std::vector<int> Fnum;

  // M-6HBAi 12,
  // TBA-6x8 2,
  const int lat_case = 2;

  switch (lat_case) {
  case 1:
    // M-6HBAi.
    Fnum.push_back(ElemIndex("sextmark"));
    Fnum.push_back(ElemIndex("sd1"));
    Fnum.push_back(ElemIndex("sd2"));
    break;
  case 2:
    // H-6-BA.
    Fnum.push_back(ElemIndex("sfh"));
    Fnum.push_back(ElemIndex("sd1a"));
    Fnum.push_back(ElemIndex("sd1b"));
    break;
  default:
    printf("\nchk_mpole: unknown lattice type\n");
    exit(1);
    break;
  }

  Ring_GetTwiss(true, 0e0);
 
  printf("\nSextupole Scheme:\n");
  for (k = 0; k < (int)Fnum.size(); k++)
    chk_mpole_Fam(Fnum[k], false);
}


void chk_dip(void)
{
  int    k;
  double L, phi, L_sum, phi_sum;

  Ring_GetTwiss(true, 0e0);
 
  printf("\nLong grad dipole:\n");
  L_sum = 0e0; phi_sum = 0e0;
  for (k = 0; k <= globval.Cell_nLoc; k++) {
    if ((Cell[k].Elem.Pkind == Mpole) && (Cell[k].Elem.M->Pirho != 0e0)) {
      L = Cell[k].Elem.PL;
      phi = L*Cell[k].Elem.M->Pirho*180e0/M_PI;
      L_sum += L; phi_sum += phi;
      printf(" %6s %4.3f %7.3f %9.6f %9.6f %9.6f %9.6f %9.6f\n",
	     Cell[k].Elem.PName, L, 1e0/Cell[k].Elem.M->Pirho,
	     phi, Cell[k].Elem.M->PTx1, Cell[k].Elem.M->PTx2,
	     L_sum, phi_sum);
    }
  }
}


void dnu_mpole(void)
{
  long int         loc1, loc0 = 0;
  int              n;
  std::vector<int> Fnum;

  const int lat_case = 1;

  switch (lat_case) {
  case 1:
    // H-6-BA.
    Fnum.push_back(ElemIndex("sf"));
    // Fnum.push_back(ElemIndex("sda"));
    // Fnum.push_back(ElemIndex("sdb"));
    break;
  default:
    printf("\ndnu_mpole: unknown lattice type\n");
    exit(1);
    break;
  }

  Ring_GetTwiss(true, 0e0);
 
  printf("\nMultipole Phase Advances:\n");
  for (n = 1; n <= GetnKid(Fnum[0]); n++) {
    loc1 = Elem_GetPos(Fnum[0], n);
    if (n == 1)
      printf("%10s %7.5f %7.5f\n",
	     Cell[loc1].Elem.PName, Cell[loc1].Nu[X_], Cell[loc1].Nu[Y_]);
    else
      printf("%10s %7.5f %7.5f\n",
	     Cell[loc1].Elem.PName,
	     Cell[loc1].Nu[X_]-Cell[loc0].Nu[X_],
	     Cell[loc1].Nu[Y_]-Cell[loc0].Nu[Y_]);
    loc0 = loc1;
  }
}


double get_pole_tip_field(const double Brho, const double R_ref,
			  const int n, const double bn)
{
  return Brho*bn*pow(R_ref, n-1);
}


void pole_tip_field(const double R_ref)
{
  int    k, n;
  double phi, b1, bn, an;

  const double Brho = globval.Energy*1e9/c0;

  for (k = 0; k <= globval.Cell_nLoc; k++)
    if (Cell[k].Elem.Pkind == Mpole) {
      switch (Cell[k].Elem.M->n_design) {
      case Dip:
	phi = Cell[k].Elem.PL*Cell[k].Elem.M->Pirho;
	b1 = Cell[k].Elem.M->Pirho;
	n = Quad;
	get_bn_design_elem(Cell[k].Fnum, 1, n, bn, an);
	printf("%10s L = %5.3f b_%1d = %7.3f     B^ = %6.3f phi = %6.3f\n",
	       Cell[k].Elem.PName, Cell[k].Elem.PL,
	       1, b1, Brho*b1, phi*180e0/M_PI);
	printf("                          b_%1d = %7.3f     B^ = %6.3f\n",
	       n, bn, get_pole_tip_field(Brho, R_ref, n, bn));
	break;
      case Quad:
	n = Cell[k].Elem.M->n_design;
	get_bn_design_elem(Cell[k].Fnum, 1, n, bn, an);
	printf("%10s L = %5.3f b_%1d = %7.3f     B^ = %6.3f\n",
	       Cell[k].Elem.PName, Cell[k].Elem.PL, n, bn,
	       get_pole_tip_field(Brho, R_ref, n, bn));
	break;
      case Sext:
	n = Cell[k].Elem.M->n_design;
	get_bn_design_elem(Cell[k].Fnum, 1, n, bn, an);
	printf("%10s L = %5.3f b_%1d = %11.3e B^ = %6.3f\n",
	       Cell[k].Elem.PName, Cell[k].Elem.PL, n, bn,
	       get_pole_tip_field(Brho, R_ref, n, bn));
	break;
      case Dodec:
	n = Cell[k].Elem.M->n_design - 2;
	get_bn_design_elem(Cell[k].Fnum, 1, n, bn, an);
	printf("%10s L = %5.3f b_%1d = %11.3e B^ = %6.3f\n",
	       Cell[k].Elem.PName, Cell[k].Elem.PL, n, bn,
	       get_pole_tip_field(Brho, R_ref, n, bn));
	n = Cell[k].Elem.M->n_design;
	get_bn_design_elem(Cell[k].Fnum, 1, n, bn, an);
	printf("%10s L = %5.3f b_%1d = %11.3e B^ = %6.3f\n",
	       Cell[k].Elem.PName, Cell[k].Elem.PL, n, bn,
	       get_pole_tip_field(Brho, R_ref, n, bn));
	break;
      default:
	printf("\npole_tip_field: unknown lattice type\n");
	exit(1);
	break;
      }
    }
}


void get_dbeta_deta(const double delta)
{
  // Evaluate derivative; to avoid effect of tune shift.
  int            j, k;
  vector<double> dbeta[2], deta_x;
  FILE           *outf;

  const double d_delta = 1e-5;

  const string file_name = "dbeta_deta.out";

  outf = file_write(file_name.c_str());

  printf("\nOptics for delta = %10.3e\n", d_delta);
  Ring_GetTwiss(true, d_delta); printglob();
  for (k = 0; k <= globval.Cell_nLoc; k++) {
    for (j = 0; j < 2; j++)
      dbeta[j].push_back(Cell[k].Beta[j]);
    deta_x.push_back(Cell[k].Eta[X_]);
  }
  printf("\nOptics for delta = %10.3e\n", -d_delta);
  Ring_GetTwiss(true, -d_delta); printglob();
  for (k = 0; k <= globval.Cell_nLoc; k++) {
    for (j = 0; j < 2; j++) {
      dbeta[j][k] -= Cell[k].Beta[j]; dbeta[j][k] /= (2e0*d_delta);
    }
    deta_x[k] -= Cell[k].Eta[X_]; deta_x[k] /= (2e0*d_delta);
    fprintf(outf, "%4d %10s %8.3f %4.1f %12.5e %12.5e %12.5e\n",
	    k, Cell[k].Elem.PName, Cell[k].S, get_code(Cell[k]),
	    dbeta[X_][k], dbeta[Y_][k], deta_x[k]);
  }

  fclose(outf);
}


ss_vect<tps> get_sympl_form(const int dof)
{
  int          k;
  ss_vect<tps> Id, omega;

  Id.identity(); omega.zero();
  for (k = 0; k < dof; k++) {
    omega[2*k] = Id[2*k+1]; omega[2*k+1] = -Id[2*k];
  }
  return omega;
}


void A_At_pass(void)
{
  long int     lastpos;
  int          i;
  ss_vect<tps> A, A_Atp;

  A.identity(); putlinmat(4, globval.Ascr, A);
  A_Atp = A*tp_S(2, A);
  printf("\n    alpha_x  beta_x    alpha_y  beta_y:\n"
	 "  %9.5f %8.5f %9.5f %8.5f\n",
	 -A_Atp[x_][px_], A_Atp[x_][x_], -A_Atp[y_][py_], A_Atp[y_][y_]);
  for (i = 0; i <= globval.Cell_nLoc; i++) {
    Cell_Pass(i, i, A_Atp, lastpos);
    A_Atp = tp_S(2, A_Atp);
    Cell_Pass(i, i, A_Atp, lastpos);
    A_Atp = tp_S(2, A_Atp);
    printf("  %9.5f %8.5f %9.5f %8.5f\n",
	   -A_Atp[x_][px_], A_Atp[x_][x_], -A_Atp[y_][py_], A_Atp[y_][y_]);
  }
}


void curly_H_s(void)
{
  long int        lastpos;
  int             i;
  double          dnu[2];
  ss_vect<double> eta, eta_Fl;
  ss_vect<tps>    A;
  FILE            *outf;

  outf = file_write("curly_H_s.out");
  printf("\n");
  eta.zero(); eta[0] = Cell[0].Eta[X_]; eta[1] = Cell[0].Etap[X_];
  A.identity(); putlinmat(2, globval.Ascr, A);
  for (i = 0; i <= globval.Cell_nLoc; i++) {
    eta.zero(); eta[0] = Cell[i].Eta[X_]; eta[1] = Cell[i].Etap[X_];
    Cell_Pass(i, i, A, lastpos); A = get_A_CS(2, A, dnu);

    eta_Fl = (Inv(A)*eta).cst();

    fprintf(outf, "  %6.3f %10.3e %10.3e\n",
	    Cell[i].S, eta_Fl[x_], eta_Fl[px_]);

  }
  fclose(outf);
}


void prt_eta_Fl(void)
{
  long int        lastpos;
  int             i, k;
  double          dphi, s, mu_x, alpha1_x, beta1_x, curly_H;
  ss_vect<double> eta0, eta_Fl0, eta_Fl, omega_M_eta, A_Atp_omega_M_eta;
  ss_vect<tps>    Id, A, M, R, A_Atp0, A_Atp, Omega;
  FILE            *outf;

  const int    n_step = 25;
  const double
    L        = 0.75,
    phi      = 5.0,
    rho      = L/(5.0*M_PI/180e0),
    alpha0_x = 0.0,
    beta0_x  = 0.19177,
    gamma0_x = (1e0+sqr(alpha0_x))/beta0_x;

  Id.identity();

  outf = file_write("eta_Fl.out");

  A.identity(); putlinmat(2, globval.Ascr, A);

  eta0.zero(); eta0[x_] = Cell[0].Eta[X_]; eta0[px_] = Cell[0].Etap[X_];
  eta_Fl0 = (Inv(A)*eta0).cst();

  A_Atp0.identity();
  A_Atp0[x_] = Cell[0].Beta[X_]*Id[x_] - Cell[0].Alpha[X_]*Id[px_];
  A_Atp0[px_] =
    -Cell[0].Alpha[X_]*Id[x_]
    + (1e0+sqr(Cell[0].Alpha[X_]))*Id[px_]/Cell[0].Beta[X_];

  Omega.identity(); Omega[x_] = Id[px_]; Omega[px_] = -Id[x_];

  R.identity(); M.identity();

  for (i = 0; i <= n_step; i++) {
    s = i*L/n_step;

    mu_x = atan(s/beta0_x);

    M[x_] =
      cos(s/rho)*Id[x_] + rho*sin(s/rho)*Id[px_] + rho*(1e0-cos(s/rho));
    M[px_] = -sin(s/rho)/rho*Id[x_] + cos(s/rho)*Id[px_] + sin(s/rho);

    A_Atp = (M-M.cst())*A_Atp0*tp_S(1, M-M.cst());
    beta1_x = A_Atp[x_][x_]; alpha1_x = -A_Atp[px_][x_];

    if (true) {
      R[x_] =
	cos(mu_x)*Id[x_] + sin(mu_x)*Id[px_]
	+ rho*(1e0-cos(s/rho))/sqrt(beta1_x);
      R[px_] =
	-sin(mu_x)*Id[x_] + cos(mu_x)*Id[px_]
	+ (beta1_x*sin(s/rho)+alpha1_x*rho*(1e0-cos(s/rho)))/sqrt(beta1_x);

      eta_Fl = (R*eta_Fl0).cst();

      curly_H = sqr(eta_Fl[x_]) + sqr(eta_Fl[px_]);
    } else {
      omega_M_eta = (Omega*M*eta0).cst();
      A_Atp_omega_M_eta = (A_Atp*omega_M_eta).cst();
      curly_H = 0e0;
      for (k = 0; k < 2; k++)
	curly_H += omega_M_eta[k]*A_Atp_omega_M_eta[k];
    }


    fprintf(outf, "  %6.3f %6.3f  %6.3f %6.3f %10.3e %10.3e %10.3e\n",
	    s, mu_x/(2e0*M_PI), alpha1_x, beta1_x, eta_Fl[x_], eta_Fl[px_],
	    curly_H);

  }

  fclose(outf);
}


void track(const string fname, const int n, const double x, const double p_x,
	   const double y, const double p_y, const double delta)
{
  long int        lastpos;
  int             k;
  ss_vect<double> ps;
  ofstream        outf;

  file_wr(outf, fname.c_str());

  ps.zero();
  ps[x_] = x; ps[px_] = p_x; ps[y_] = y; ps[py_] = p_y; ps[delta_] = delta;

  outf << std::scientific << std::setprecision(6)
       << "\n" << std::setw(14) << ps << "\n"; 
  for (k = 1; k <= n; k++) {
    Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);  
    outf << std::scientific << std::setprecision(6)
	 << std::setw(14) << ps << "\n";
  }

  outf.close();
}


void prt_mat(const int n, const Matrix &A)
{
  int i, j;

  printf("matrix:\n");
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++)
      printf(" %18.15f", A[i][j]);
    printf("\n");
  }
}


void get_disp(void)
{
  int             k;
  long int        jj[ss_dim], lastn, lastpos;
  double          twoJ[2], curly_H[2], ds[2], ds0, ds_hat, delta_mean[2],
                  delta_hat, phi_x, f_rf, m_56, m_65, m_66, alpha_s, beta_s,
                  gamma_s, nu_s, alpha_c, C;
  ss_vect<double> eta, A, ps, D;
  ss_vect<tps>    Ascr, Id, M;
  ofstream        outf;

  f_rf = Cell[Elem_GetPos(ElemIndex("cav"), 1)].Elem.C->Pfreq;
  printf("\nf_rf = %10.3e\n", f_rf);

  Id.identity();

  globval.Cavity_on = !false; globval.radiation = false;
  Ring_GetTwiss(true, 0e0); printglob();
  nu_s = -globval.TotalTune[Z_];
  alpha_c = globval.Alphac;
  C = Cell[globval.Cell_nLoc].S;

#if 0
  // Expanded.
  printf("\n  m_66, m_65 = %21.14e %21.14e\n",
	 1e0, -sqr(2e0*M_PI*nu_s)/(alpha_c*C));
  printf("  m_56, m_55 = %21.14e %21.14e\n", -alpha_c*C, 1e0);
#else
  // "Exact".
  alpha_s =
    -globval.Ascr[ct_][ct_]*globval.Ascr[delta_][ct_]
    - globval.Ascr[ct_][delta_]*globval.Ascr[delta_][delta_];
  beta_s = sqr(globval.Ascr[ct_][ct_]) + sqr(globval.Ascr[ct_][delta_]);
  gamma_s = (1e0+sqr(alpha_s))/beta_s;

  printf("\n  nu_s       = %15.8e\n", nu_s);
  printf("  alpha_s    = %10.3e\n", alpha_s);
  printf("  beta_s     = %10.3e\n", beta_s);
  printf("  gamma_s    = %10.3e\n", gamma_s);

  printf("\n  m_66, m_65 = %21.14e %21.14e\n",
	 cos(2e0*M_PI*nu_s)-alpha_s*sin(2e0*M_PI*nu_s),
	 -gamma_s*sin(2e0*M_PI*nu_s));
  printf("  m_56, m_55 = %21.14e %21.14e\n",
	 beta_s*sin(2e0*M_PI*nu_s),
	 cos(2e0*M_PI*nu_s)+alpha_s*sin(2e0*M_PI*nu_s));
#endif

  globval.Cavity_on = false; globval.radiation = false;
  Ring_GetTwiss(true, 0e0); printglob();

  putlinmat(6, globval.OneTurnMat, M);
  D[x_] = M[x_][x_]*M[px_][delta_] - M[px_][x_]*M[x_][delta_];
  D[px_] = M[x_][px_]*M[px_][delta_] - M[px_][px_]*M[x_][delta_];
  printf("\n  m_51, m_52 = %13.6e %13.6e\n", D[x_], D[px_]);
  printf("  m_61, m_62 = %13.6e %13.6e\n",
	 -gamma_s*sin(2e0*M_PI*nu_s)*D[x_], -gamma_s*sin(2e0*M_PI*nu_s)*D[px_]);

  A.zero();
  A[x_] = 10e-6; A[px_] = 0e-6; A[y_] = 0e-6; A[py_] = 0e-6; A[delta_] = 0e-3;
  Ascr.zero();
  putlinmat(4, globval.Ascr, Ascr);
  get_twoJ(2, A, Ascr, twoJ);

  eta.zero();
  eta[x_] = Cell[globval.Cell_nLoc].Eta[X_];
  eta[px_] = Cell[globval.Cell_nLoc].Etap[X_];
  get_twoJ(1, eta, Ascr, curly_H);

  ds0 =
    Cell[globval.Cell_nLoc].Etap[X_]*A[x_]
    - Cell[globval.Cell_nLoc].Eta[X_]*A[px_];
  for (k = 0; k < 2; k++)
    ds[k] = M_PI*globval.Chrom[k]*twoJ[k];
  ds_hat = sqrt(twoJ[X_]*curly_H[X_]);

  for (k = 0; k < 2; k++)
    delta_mean[k] = ds[k]/(alpha_c*C);
  delta_hat =
    sqr(2e0*M_PI*nu_s)*ds_hat
    /(alpha_c*C*sin(M_PI*globval.TotalTune[X_]));

  alpha_s = (1e0-cos(2e0*M_PI*nu_s))/sin(2e0*M_PI*nu_s);
  beta_s  = alpha_c*C/sin(2e0*M_PI*nu_s);
  M.identity();
  M[ct_] =
    (cos(2e0*M_PI*nu_s)+alpha_s*sin(2e0*M_PI*nu_s))*Id[ct_]
    + beta_s*sin(2e0*M_PI*nu_s)*Id[delta_];
  M[delta_] =
    -(1e0+sqr(alpha_s))*sin(2e0*M_PI*nu_s)/beta_s*Id[ct_]
    + (cos(2e0*M_PI*nu_s)-alpha_s*sin(2e0*M_PI*nu_s))*Id[delta_];
  prt_lin_map(3, M);

  printf("\n  alpha_c                   = %9.3e\n", alpha_c);
  printf("  A_x                       = %9.3e [micron]\n", 1e6*A[X_]);
  printf("  2*J                       = %9.3e %9.3e\n", twoJ[X_], twoJ[Y_]);
  printf("  curly_H                   = %9.3e\n", curly_H[X_]);
  printf("\n  ds0                       = %7.5f [micron]\n", 1e6*ds0);
  printf("  ds = 2*pi*ksi*J           = %9.3e %9.3e [micron]\n",
	 1e6*ds[X_], 1e6*ds[Y_]);
  printf("  ds^ = sqrt(2*J_x*curly_H) = %7.5f [micron]\n", 1e6*ds_hat);
  printf("\n  nu_s                      = %10.5e\n", nu_s);
  printf("  delta_mean                = %9.3e %9.3e\n",
	 delta_mean[X_], delta_mean[Y_]);
  printf("  delta_hat                 = %9.3e\n", delta_hat);
  printf("\n  nu_s                      = %10.3e\n", nu_s);
  printf("  alpha_s                   = %10.3e\n", alpha_s);
  printf("  beta_s                    = %10.3e\n", beta_s);

  if (!false) {
    globval.Cavity_on = false; globval.radiation = false;
    Ring_GetTwiss(true, 0e0); printglob();

    printf("\ndet{M}-1 = %12.5e\n", DetMat(6, globval.OneTurnMat)-1e0);

    globval.Cavity_on = true; globval.radiation = false;
    Ring_GetTwiss(true, 0e0); printglob();

    printf("\ndet{M}-1 = %12.5e\n", DetMat(6, globval.OneTurnMat)-1e0);

    globval.alpha_z =
      -globval.Ascr[ct_][ct_]*globval.Ascr[delta_][ct_]
      - globval.Ascr[ct_][delta_]*globval.Ascr[delta_][delta_];
    globval.beta_z =
      sqr(globval.Ascr[ct_][ct_]) + sqr(globval.Ascr[ct_][delta_]);
    globval.TotalTune[Z_] = fabs(globval.TotalTune[Z_]);

    printf("\nnu_z    = %12.5e\n", globval.TotalTune[Z_]);
    printf("beta_z  = %12.5e %12.5e\n",
	   globval.beta_z, alpha_c*C/sin(2e0*M_PI*globval.TotalTune[Z_]));
    printf("alpha_z = %12.5e %12.5e %12.5e+O(nu_s)^2\n",
	   globval.alpha_z,
	   (1e0-cos(2e0*M_PI*globval.TotalTune[Z_]))
	   /sin(2e0*M_PI*globval.TotalTune[Z_]),
	   M_PI*globval.TotalTune[Z_]);

    M.zero();
    M[ct_] =
      (cos(2e0*M_PI*globval.TotalTune[Z_])
       +globval.alpha_z*sin(2e0*M_PI*globval.TotalTune[Z_]))*Id[ct_]
      + globval.beta_z*sin(2e0*M_PI*globval.TotalTune[Z_])*Id[delta_];
    M[delta_] =
      -(1e0+sqr(globval.alpha_z))/globval.beta_z
      *sin(2e0*M_PI*globval.TotalTune[Z_])*Id[ct_]
      + (cos(2e0*M_PI*globval.TotalTune[Z_])-globval.alpha_z
	 *sin(2e0*M_PI*globval.TotalTune[Z_]))*Id[delta_];

    phi_x =
      atan2(
	    Cell[globval.Cell_nLoc].Alpha[X_]*Cell[globval.Cell_nLoc].Eta[X_]
	    +Cell[globval.Cell_nLoc].Beta[X_]*Cell[globval.Cell_nLoc].Etap[X_],
	    Cell[globval.Cell_nLoc].Eta[X_]);

    prt_lin_map(3, M);
    printf("\ndet{M}-1 = %12.5e\n",
	   M[ct_][ct_]*M[delta_][delta_]-M[ct_][delta_]*M[delta_][ct_]-1e0);
    printf("phi_x  = %12.5e\n", phi_x*180e0/M_PI);
  }

  globval.Cavity_on = !false; globval.radiation = false;
  track("track.out", A[x_], A[px_], A[y_], A[py_], A[delta_], 2000,
	lastn, lastpos,	0, 0*f_rf);

  if (!false) {
    // Standard Map.
    file_wr(outf, "std_map.out");

    printf("\n");
    prtmat(6, globval.OneTurnMat);
    m_56 = globval.OneTurnMat[ct_][delta_];
    m_65 = globval.OneTurnMat[delta_][ct_];

    ps.zero();
    for (k = 1; k <= 2000; k++) {
      ps[delta_] +=
	-sqr(2e0*M_PI*nu_s)/(alpha_c*C)*ps[ct_];
      ps[ct_] +=
	alpha_c*C*(ps[delta_]+ds_hat*sin(k*2e0*M_PI*globval.TotalTune[X_]));
      outf << scientific << setprecision(5)
	   << setw(5) << k << setw(13) << ps << "\n";
    }
    outf.close();
  }
}


void get_Poincare_Map(void)
{
  int          k;
  double       C, mu[3], alpha[3], beta[3], gamma[3], tau[3];
  ss_vect<tps> Id, M, M_rad;

  if (false) no_sxt();

  Id.identity();

  globval.Cavity_on = true; globval.radiation = true;
  Ring_GetTwiss(true, 0e0); printglob();

  C = Cell[globval.Cell_nLoc].S;
  mu[Z_] = -2e0*M_PI*globval.Omega;
  alpha[Z_] =
    -globval.Ascr[ct_][ct_]*globval.Ascr[delta_][ct_]
    - globval.Ascr[ct_][delta_]*globval.Ascr[delta_][delta_];
  beta[Z_] = sqr(globval.Ascr[ct_][ct_]) + sqr(globval.Ascr[ct_][delta_]);
  gamma[Z_] = (1e0+sqr(alpha[Z_]))/beta[Z_];

  printf("\n  %12.5e %12.5e %12.5e\n", -globval.Omega, alpha[Z_], beta[Z_]);
  printf("  %12.5e %12.5e\n",
	 -globval.Ascr[ct_][ct_]*globval.Ascr[delta_][ct_],
	 -globval.Ascr[ct_][delta_]*globval.Ascr[delta_][delta_]);

  for (k = 0; k < 2; k++) {
    mu[k] = 2e0*M_PI*globval.TotalTune[k];
    alpha[k] = Cell[0].Alpha[k]; beta[k] = Cell[0].Beta[k];
    gamma[k] = (1e0+sqr(alpha[k]))/beta[k];
  }

  M.zero();
  for (k = 0; k < 3; k++) {
    if (k < 2) {
      M[2*k] =
	(cos(mu[k])+alpha[k]*sin(mu[k]))*Id[2*k] + beta[k]*sin(mu[k])*Id[2*k+1];
      M[2*k+1] =
	-gamma[k]*sin(mu[k])*Id[2*k]
	+ (cos(mu[k])-alpha[k]*sin(mu[k]))*Id[2*k+1];
    } else {
#if 1
      M[2*k+1] =
	(cos(mu[k])+alpha[k]*sin(mu[k]))*Id[2*k+1] + beta[k]*sin(mu[k])*Id[2*k];
      M[2*k] =
	-gamma[k]*sin(mu[k])*Id[2*k+1]
	+ (cos(mu[k])-alpha[k]*sin(mu[k]))*Id[2*k];
#else
      M[2*k+1] = Id[2*k+1] + globval.OneTurnMat[ct_][delta_]*Id[2*k];
      M[2*k] =
	-sqr(mu[Z_])/globval.OneTurnMat[ct_][delta_]*Id[2*k+1]
	+ (1e0-sqr(mu[Z_]))*Id[2*k];
#endif
    }
  }
  M[x_] += globval.OneTurnMat[x_][delta_]* Id[delta_];
  M[px_] += globval.OneTurnMat[px_][delta_]*Id[delta_];
  M[ct_] +=
    (M[x_][x_]*M[px_][delta_]-M[px_][x_]*M[x_][delta_])*Id[x_]
    +(M[x_][px_]*M[px_][delta_]-M[px_][px_]*M[x_][delta_])*Id[px_];
  M[delta_] +=
    -sqr(mu[Z_])*M[ct_][x_]/M[ct_][delta_]*Id[x_]
    -sqr(mu[Z_])*M[ct_][px_]/M[ct_][delta_]*Id[px_];

  M_rad.zero();
  for (k = 0; k < 3; k++) {
    tau[k] = -C/(c0*globval.alpha_rad[k]);
    M_rad[2*k] = exp(-C/(c0*tau[k]))*Id[2*k];
    M_rad[2*k+1] = exp(-C/(c0*tau[k]))*Id[2*k+1];
  }
  M = M_rad*M;
  printf("\n  %12.5e %12.5e %12.5e\n", tau[X_], tau[Y_], tau[Z_]);

  prt_lin_map(3, M);
}


void get_matrix(const string &name, const double delta)
{
  int          k;
  double       L, rho, b2, K[2], psi[2];
  elemtype     Elem;
  ss_vect<tps> Id, map;

  // prt_mat(6, globval.OneTurnMat);

  Id.identity();

  Elem = Cell[Elem_GetPos(ElemIndex(name.c_str()), 1)].Elem;
  L = Elem.PL;
  rho = 1e0/Elem.M->Pirho;
  b2  = Elem.M->PBpar[Quad+HOMmax];
  printf("\n  L = %7.5f rho = %7.5f b_2 = %7.5f \n", L, rho, b2);
  K[X_] = b2 + 1e0/sqr(rho); K[Y_] = b2;
  for (k = 0; k < 2; k++)
    psi[k] = sqrt(fabs(K[k])/(1e0+delta))*L;

  map.identity();
  if (K[X_] >= 0e0) {
    map[x_] =
      cos(psi[X_])*Id[x_] + sin(psi[X_])/(sqrt(K[X_]*(1e0+delta)))*Id[px_]
      + (1e0-cos(psi[X_]))/(rho*K[X_])*Id[delta_];
    map[px_] =
      -sqrt(K[X_]*(1e0+delta))*sin(psi[X_])*Id[x_] + cos(psi[X_])*Id[px_]
      + sin(psi[X_])*sqrt(1e0+delta)/(rho*sqrt(K[X_]))*Id[delta_];
    map[y_] =
      (psi[Y_] != 0e0)?
      cosh(psi[Y_])*Id[y_] + sinh(psi[Y_])/(sqrt(K[Y_]*(1e0+delta)))*Id[py_]
      :
      cosh(psi[Y_])*Id[y_] + L*(1e0+delta)*Id[py_];
    map[py_] =
      sqrt(K[Y_]*(1e0+delta))*sinh(psi[Y_])*Id[y_] + cosh(psi[Y_])*Id[py_];
    map[ct_] +=
      sin(psi[X_])*sqrt(1e0+delta)/(rho*sqrt(K[X_]))*Id[x_]
      + (1e0-cos(psi[X_]))/(rho*K[X_])*Id[px_]
      + (psi[X_]-sin(psi[X_]))*sqrt(1e0+delta)
      /(sqr(rho)*pow(K[X_], 3e0/2e0))*Id[delta_];
  } else {
    printf("\nK_x < 0\n");
    K[X_] = -K[X_]; K[Y_] = -K[Y_];
    map[x_] =
      cosh(psi[X_])*Id[x_] + sinh(psi[X_])/(sqrt(K[X_]*(1e0+delta)))*Id[px_]
      - (1e0-cosh(psi[X_]))/(rho*K[X_])*Id[delta_];
    map[px_] =
      sqrt(K[X_]*(1e0+delta))*sinh(psi[X_])*Id[x_] + cosh(psi[X_])*Id[px_]
      + sinh(psi[X_])*sqrt(1e0+delta)/(rho*sqrt(K[X_]))*Id[delta_];
    map[y_] =
      cos(psi[Y_])*Id[y_] + sin(psi[Y_])/(sqrt(K[Y_]*(1e0+delta)))*Id[py_];
    map[py_] =
      -sqrt(K[Y_]*(1e0+delta))*sin(psi[Y_])*Id[y_] + cos(psi[Y_])*Id[py_];
    map[ct_] =
      sinh(psi[X_])*sqrt(1e0+delta)/(rho*sqrt(K[X_]))*Id[x_]
      - (1e0-cosh(psi[X_]))/(rho*K[X_])*Id[px_]
      - (psi[X_]-sinh(psi[X_]))*sqrt(1e0+delta)
      /(sqr(rho)*pow(K[X_], 3e0/2e0))*Id[delta_];
  }

  prt_lin_map(3, map);
}


void get_eta(void)
{
  int             k;
  long int        jj[ss_dim];
  ss_vect<double> D, eta;
  ss_vect<tps>    Id, M;

  Id.identity();

  Ring_GetTwiss(true, 0e0); printglob();
  prt_lat("linlat1.out", globval.bpm, true);

  for (k = 0; k < ss_dim; k++)
    jj[k] = 0;
  jj[x_] = 1; jj[px_] = 1;
  putlinmat(2, globval.OneTurnMat, M);
  prt_lin_map(1, PInv(Id-M, jj));

  M[x_] =
    (1e0+Cell[globval.Cell_nLoc].Alpha[X_]/tan(M_PI*globval.TotalTune[X_]))
    /2e0*Id[x_]
    + Cell[globval.Cell_nLoc].Beta[X_]/(2e0*tan(M_PI*globval.TotalTune[X_]))
    *Id[px_];
  M[px_] =
    -(1e0+sqr(Cell[globval.Cell_nLoc].Alpha[X_]))
    /(2e0*Cell[globval.Cell_nLoc].Beta[X_]*tan(M_PI*globval.TotalTune[X_]))
    *Id[x_]
    + (1e0-Cell[globval.Cell_nLoc].Alpha[X_]/tan(M_PI*globval.TotalTune[X_]))
    /2e0*Id[px_];

  for (k = 0; k < 2; k++)
    D[k] = globval.OneTurnMat[k][delta_];

  eta = (M*D).cst();

  prt_lin_map(1, M);
  prt_lin_map(1, D);
  prt_lin_map(1, eta);

  printf("eta = %13.6e %13.6e\n", eta[x_], eta[px_]);
}


void orm(const string &bpm, const int i,
	 const string &corr, const int j)
{
  long int loc_bpm, loc_corr;
  double   nu, spiq, betai, betaj, nui, nuj, A_ij;

  nu = globval.TotalTune[X_]; spiq = sin(M_PI*nu);
  loc_bpm = Elem_GetPos(ElemIndex(bpm.c_str()), i);
  betai = Cell[loc_bpm].Beta[X_]; nui = Cell[loc_bpm].Nu[X_];
  loc_corr = Elem_GetPos(ElemIndex(corr.c_str()), j);
  betaj = Cell[loc_corr].Beta[X_]; nuj = Cell[loc_corr].Nu[X_];
  A_ij = sqrt(betai*betaj)/(2e0*spiq)*cos(nu*M_PI-fabs(2e0*M_PI*(nui-nuj)));

  printf("\norm:     A_ij = %12.5e\n", A_ij);
}


void orm_num(const string &bpm, const int i,
	     const string &corr, const int j, const double eps)
{
  long int        lastpos, loc_bpm, loc_corr;
  double          A_ij, x0, x1;

  loc_bpm = Elem_GetPos(ElemIndex(bpm.c_str()), i);
  loc_corr = Elem_GetPos(ElemIndex(corr.c_str()), j);
  set_dbnL_design_elem(Cell[loc_corr].Fnum, Cell[loc_corr].Knum, Dip, eps, 0e0);
  getcod(0.0, lastpos);
  x1 = Cell[loc_bpm].BeamPos[x_];
  set_dbnL_design_elem(Cell[loc_corr].Fnum, Cell[loc_corr].Knum, Dip, -2e0*eps,
		       0e0);
  getcod(0.0, lastpos);
  x0 = Cell[loc_bpm].BeamPos[x_];
  set_dbnL_design_elem(Cell[loc_corr].Fnum, Cell[loc_corr].Knum, Dip, eps, 0e0);
  A_ij = (x1-x0)/(2e0*eps);

  printf("orm_num: A_ij = %12.5e\n", A_ij);
}


void trm(const string &bpm, const int i,
	 const string &corr, const int j)
{
  long int loc_bpm, loc_corr;
  double   betai, betaj, nui, nuj, A_ij;

  loc_bpm = Elem_GetPos(ElemIndex(bpm.c_str()), i);
  betai = Cell[loc_bpm].Beta[X_]; nui = Cell[loc_bpm].Nu[X_];
  loc_corr = Elem_GetPos(ElemIndex(corr.c_str()), j);
  betaj = Cell[loc_corr].Beta[X_]; nuj = Cell[loc_corr].Nu[X_];
  A_ij =
    (loc_bpm > loc_corr)?
    sqrt(betai*betaj)*sin(2e0*M_PI*(nui-nuj)) : 0e0;

  printf("\ntrm:     A_ij = %12.5e\n", A_ij);
}


void trm_num(const string &bpm, const int i,
	     const string &corr, const int j, const double eps)
{
  long int        lastpos, loc_bpm, loc_corr;
  double          A_ij;
  ss_vect<double> ps0, ps1;

  loc_bpm = Elem_GetPos(ElemIndex(bpm.c_str()), i);
  loc_corr = Elem_GetPos(ElemIndex(corr.c_str()), j);
  set_dbnL_design_elem(Cell[loc_corr].Fnum, Cell[loc_corr].Knum, Dip, eps, 0e0);
  ps1.zero();
  Cell_Pass(0, loc_bpm, ps1, lastpos);
  set_dbnL_design_elem(Cell[loc_corr].Fnum, Cell[loc_corr].Knum, Dip, -2e0*eps,
		       0e0);
  ps0.zero();
  Cell_Pass(0, loc_bpm, ps0, lastpos);
  set_dbnL_design_elem(Cell[loc_corr].Fnum, Cell[loc_corr].Knum, Dip, eps, 0e0);
  A_ij = (ps1[x_]-ps0[x_])/(2e0*eps);

  printf("trm_num: A_ij = %12.5e\n", A_ij);
}


void wtf(void)
{
  globval.Cavity_on = true; globval.radiation = !false;
  Ring_GetTwiss(true, 0e0); printglob();
  globval.alpha_z =
    -globval.Ascr[ct_][ct_]*globval.Ascr[delta_][ct_]
    - globval.Ascr[ct_][delta_]*globval.Ascr[delta_][delta_];
  globval.beta_z = sqr(globval.Ascr[ct_][ct_]) + sqr(globval.Ascr[ct_][delta_]);
 
  printf("\nLattice Parameters:\n  alpha = [%9.5f, %9.5f, %9.5f]\n",
	 Cell[0].Alpha[X_], Cell[0].Alpha[Y_], globval.alpha_z);
  printf("  beta  = [%9.5f, %9.5f, %9.5f]\n",
	 Cell[0].Beta[X_], Cell[0].Beta[Y_], globval.beta_z);
}


int main(int argc, char *argv[])
{
  bool             tweak;
  long int         lastn, lastpos, loc, loc2;
  int              k, b2_fam[2], b3_fam[2], lat_case;
  double           b2[2], a2, b3[2], b3L[2], a3, a3L, f_rf, dx, dnu[3], I[6];
  tps              a;
  Matrix           M;
  std::vector<int> Fam;
  ss_vect<tps>     Ascr, A_Atp, Id, Ms;
  ostringstream    str;

  const long   seed   = 1121;
  const int    n_turn = 2064;
  const double delta  = 1.5e-2;
  //                   nu[]    = { 102.18/20.0, 68.30/20.0 };
  // const std::string q_fam[] = { "qfe", "qde" }, s_fam[] = { "sfh", "sd" };
  //                   nu[]    = { 39.1/12.0, 15.25/12.0 };
  //                   // nu[]    = { 3.266+0.01, 1.275 };
  // const std::string q_fam[] = { "qm2b", "qm3" }, s_fam[] = { "sfh", "sd" };
  const std::string
    q_fam[] = { "qf03", "qd04" },
    s_fam[] = { "sfh", "sd" };

  const double R_ref = 5e-3;

  // 1: DIAMOND, 2: NSLS-II, 3: Oleg I, 4: Oleg II.
  FieldMap_filetype = 4; sympl = !false;

  reverse_elem = !false;

  trace = !true;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  globval.H_exact    = false; globval.quad_fringe    = false;
  globval.Cavity_on  = false; globval.radiation      = false;
  globval.emittance  = false; globval.IBS            = false;
  globval.pathlength = false; globval.bpm            = 0;
  globval.Cart_Bend  = false; globval.dip_edge_fudge = true;

  if (false) no_sxt();

  globval.Cavity_on = false; globval.radiation = false;
  Ring_GetTwiss(true, 0e0); printglob();

  if (!false) {
    if (!false) {
      ss_vect<double> ps;
      ss_vect<tps>    A;

      trace = true;

      if (!true) {
	ps.zero(); ps[delta_] = 1e-3;
	printf("\n");
	Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
      } else {
	globval.emittance = true;

	putlinmat(6, globval.Ascr, A);
	printf("\n");
	Cell_Pass(0, globval.Cell_nLoc, A, lastpos);
      }
      exit(0);
    }

    get_eps_x();
    // get_I(I, true);
    exit(0);
  }

  if (false) {
    wtf();
    exit(0);
  }

  if (false) {
    get_Poincare_Map();
    exit(0);
  }

  if (false) {
    long int     lastpos;
    ss_vect<tps> map;

    globval.Cavity_on = !false; globval.radiation = false;

    map.identity();
    Cell_Pass(0, globval.Cell_nLoc, map, lastpos);
    prt_lin_map(3, map);
    exit(0);
  }


  if (false) {
    int             k;
    ss_vect<double> ps;
    ofstream        outf;

    file_wr(outf, "track.out");
    globval.Cavity_on = !false; globval.radiation = !false;
    ps.zero();
    ps[x_] = 0e-3; ps[px_] = 0e-3; ps[y_] = 0e-3; ps[py_] = 0e-3;
    ps[ct_] = 0e-3;
    for (k = 0; k <= globval.Cell_nLoc; k++) {
      Cell_Pass(k, k, ps, lastpos);
      outf << setw(4) << k
	   << fixed << setprecision(5) << setw(9) << Cell[k].S
	   << " " << setw(10) << Cell[k].Elem.PName
	   << scientific << setprecision(14)
	   << setw(22) << ps << "\n";
    }
    outf.close();
    exit(0);
  }

  if (false) {
    orm("bpm_11", 3, "ch_11", 1);
    orm_num("bpm_11", 3, "ch_11", 1, 1e-8);
    printf("\n");
    trm("bpm_11", 3, "ch_11", 1);
    trm_num("bpm_11", 3, "ch_11", 1, 1e-8);
    exit(0);
  }

  if (false) {
    get_matrix("dq1", 0e0);
    exit(0);
  }

  if (false) {
    get_eta();
    exit(0);
  }

  if (false) {
    globval.Cavity_on  = !false; globval.radiation = !false;
    track("track.out", 10, 0e0, 0e0, 0e0, 0e0, 0e0);
    exit(0);
  }

  if (false) {
    A_At_pass();
    exit(0);
  }

  if (false) {
    curly_H_s();
    exit(0);
  }

  if (false) {
    prt_eta_Fl();
    exit(0);
  }

  if (false) {
    chk_lin_chrom();
    // exit(0);
  }
  
  if (prt_dt) {
    printf("Lattice Case (1..3)? ");
    scanf("%d", &lat_case);

    fit_ksi1(lat_case, 0e0, 0e0);

    chk_mI_trans(lat_case);

    no_sxt();
    Ring_GetTwiss(true, 0e0); printglob();

    prt_lat("linlat1.out", globval.bpm, true);
    prt_lat("linlat.out", globval.bpm, true, 10);
    prt_chrom_lat();

    exit(0);
  }

  if (false) {
    loc = Elem_GetPos(ElemIndex("bb"), 1);
    map.identity();
    // Tweak to remain within field map range at entrance.
    tweak = true;
    if (tweak) {
      dx = -1.4e-3; map[x_] += dx;
    }
    Cell_Pass(loc, loc, map, lastpos);
    if (tweak) map[x_] -= dx;
    prt_lin_map(3, map);
    getlinmat(6, map, M);
    printf("\n1-Det: %9.3e\n", 1e0-DetMat(6, M));
    exit(0);
  }

  if (false) {
    chk_phi();
    // exit(0);
  }

  if (false) {
    dpath_length();
    exit(0);
  }

  if (false) {
    prt_drift();
    exit(0);
  }

  if (false) {
    prt_dip();
    exit(0);
  }

  if (false) {
    const double
      alpha0[] = { 0.91959,  3.19580},
      beta0[]  = { 1.27047, 17.77859},
      eta0[]   = { 0.0,      0.0},
      etap0[]  = { 0.0,      0.0},
      alpha1[] = { 0.81112,  0.45975},
      beta1[]  = { 4.56189,  2.62113},
      eta1[]   = { 0.07297,  0.0},
      etap1[]  = {-0.01063,  0.0};

    set_map_per(ElemIndex("ps_per"), alpha1, beta1, eta1, etap1);
    if (!true) {
      Ring_GetTwiss(true, 0e0); printglob();
    } else
      ttwiss(alpha0, beta0, eta0, etap0, 0e0);

    prt_lat("linlat1.out", globval.bpm, true);
    prt_lat("linlat.out", globval.bpm, true, 10);

    exit(0);
  }

  if (mI_rot) {
    set_map(ElemIndex("mI_rot"), dnu_mI);
    Ring_GetTwiss(true, 0e0); printglob();
  }

  if (HOA_rot) {
    loc = Elem_GetPos(ElemIndex("HOA_1_rot"), 1);
    for (k = 0; k < 2; k++)
      dnu[k] = fract(nu_HOA[k]) - fract(Cell[loc].Nu[k]);
    printf("\ndnu HOA_1: [%7.5f %7.5f]\n", dnu[X_], dnu[Y_]);
    set_map(ElemIndex("HOA_1_rot"), dnu);

    loc2 = Elem_GetPos(ElemIndex("HOA_2_rot"), 1);
    for (k = 0; k < 2; k++)
      dnu[k] = fract(nu_HOA[k]) - fract(Cell[loc2].Nu[k]-Cell[loc].Nu[k]);
    printf("\ndnu HOA_2: [%7.5f %7.5f]\n", dnu[X_], dnu[Y_]);
    set_map(ElemIndex("HOA_2_rot"), dnu);

    Ring_GetTwiss(true, 0e0); printglob();
  }

  if (set_dnu) {
    dnu[X_] = 0.0; dnu[Y_] = 0.0;
    set_map(ElemIndex("ps_rot"), dnu);
    Ring_GetTwiss(true, 0e0); printglob();
    for (k = 0; k < 2; k++)
      if (!true)
	dnu[k] = nu[k]/n_cell - globval.TotalTune[k];
      else
	dnu[k] = nu[k];
    printf("\ntune set to:\n  dnu     = [%8.5f, %8.5f]\n", dnu[X_], dnu[Y_]);
    printf("  nu_cell = [%8.5f, %8.5f]\n", nu[X_]/n_cell, nu[Y_]/n_cell);
    printf("  nu      = [%8.5f, %8.5f]\n", nu[X_], nu[Y_]);
    set_map(ElemIndex("ps_rot"), dnu);
    Ring_GetTwiss(true, 0e0); printglob();
  }

  prtmfile("flat_file.dat");
  // globval.bpm = ElemIndex("bpm");
  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);
  prt_chrom_lat();

  if (prt_ms) {
    loc = Elem_GetPos(ElemIndex("ms"), 1);
    printf("\n%10s:\n  {{%12.10f, %12.10f}, {%12.10f, %12.10f},"
	   " {%12.10f, %12.10f}, {%3.1f, %3.1f}};\n",
	   Cell[loc].Elem.PName,
	   Cell[loc].Alpha[X_], Cell[loc].Alpha[Y_],
	   Cell[loc].Beta[X_], Cell[loc].Beta[Y_],
	   Cell[loc].Eta[X_], Cell[loc].Eta[Y_],
	   Cell[loc].Etap[X_], Cell[loc].Etap[Y_]);
    printf("\n%10s:\n  (%12.10f, %12.10f, %12.10f, %12.10f,"
	   "\n  %12.10f, %12.10f, %3.1f, %3.1f)\n",
	   Cell[loc].Elem.PName,
	   Cell[loc].Alpha[X_], Cell[loc].Beta[X_],
	   Cell[loc].Alpha[Y_], Cell[loc].Beta[Y_],
	   Cell[loc].Eta[X_], Cell[loc].Etap[X_], 
	   Cell[loc].Eta[Y_], Cell[loc].Etap[Y_]);
    exit(0);
  }

  if (false) {
    loc = Elem_GetPos(ElemIndex("lgb00"), 1);
    printf("\n%10s:  \n  %13.10f %13.10f %13.10f %13.10f %13.10f %13.10f\n",
	   Cell[loc].Elem.PName,
	   Cell[loc].Alpha[X_], Cell[loc].Beta[X_],
	   Cell[loc].Eta[X_], Cell[loc].Etap[X_],
	   Cell[loc].Alpha[Y_], Cell[loc].Beta[Y_]);
    exit(0);
  }

  if (false) {
    printf("Lattice Case (1..3)? ");
    scanf("%d", &lat_case);

    chk_high_ord_achr(lat_case);
    // exit(0);
  }

  if (false) {
    chk_drv_terms();
    exit(0);
  }

  if (false) {
    chk_mpole();
    prt_lat("linlat1.out", globval.bpm, true);
    prt_lat("linlat.out", globval.bpm, true, 10);
    exit(0);
  }

  if (false) {
    pole_tip_field(R_ref);
    exit(0);
  }

  if (false) {
    chk_dip();
    exit(0);
  }

  if (false) {
    chk_optics(0.0, 9.634816, 0.0, 5.954390, 0.0, 0.0, 0.0, 0.0);
    prt_lat("linlat1.out", globval.bpm, true);
    prt_lat("linlat.out", globval.bpm, true, 10);
    exit(0);
  }

  if (false) {
    globval.Cavity_on = false; globval.radiation = false;
    track(0e-3, 0e-3);
    exit(0);
  }

  if (false) {
    dnu_mpole();
    exit(0);
  }

  if (false) {
    get_dbeta_deta(1e-4);
    // exit(0);
  }

  globval.Cavity_on = false; globval.radiation = false;
  Ring_GetTwiss(true, 0e0); printglob();

  if (false) get_alphac2();

  if (false) {
    iniranf(seed); setrancut(1e0);
    // globval.bpm = ElemIndex("mon");
    globval.bpm = ElemIndex("bpm");
    globval.hcorr = ElemIndex("ch"); globval.vcorr = ElemIndex("cv");
    // ALS-U.
    // globval.bpm = ElemIndex("bpm_m");
    // globval.hcorr = ElemIndex("corr_h"); globval.vcorr = ElemIndex("corr_v");

    gcmat(globval.bpm, globval.hcorr, 1);
    gcmat(globval.bpm, globval.vcorr, 2);

    get_cod_rms(50e-6, 50e-6, 100, true);

    exit(0);
  }

  if (false) {
    Fam.push_back(ElemIndex("ts1b"));
    // Fam.push_back(ElemIndex("ts1d"));
    chk_mini_beta(Fam);
    exit(0);
  }

  if (false) {
    // Fam.push_back(ElemIndex("s1b"));
    // Fam.push_back(ElemIndex("s1d"));
    // Fam.push_back(ElemIndex("s2b"));
    // Fam.push_back(ElemIndex("s2d"));
    // Fam.push_back(ElemIndex("sx1"));
    // Fam.push_back(ElemIndex("sy1"));

  switch (2) {
  case 1:
    // DIAMOND.
    Fam.push_back(ElemIndex("ts1a"));
    Fam.push_back(ElemIndex("ts1ab"));
    Fam.push_back(ElemIndex("ts2a"));
    Fam.push_back(ElemIndex("ts2ab"));
    Fam.push_back(ElemIndex("ts1b"));
    Fam.push_back(ElemIndex("ts2b"));
    Fam.push_back(ElemIndex("ts1c"));
    Fam.push_back(ElemIndex("ts2c"));
    Fam.push_back(ElemIndex("ts1d"));
    Fam.push_back(ElemIndex("ts2d"));
    Fam.push_back(ElemIndex("ts1e"));
    Fam.push_back(ElemIndex("ts2e"));

    Fam.push_back(ElemIndex("s1"));
    Fam.push_back(ElemIndex("s2"));
    Fam.push_back(ElemIndex("s3"));
    Fam.push_back(ElemIndex("s4"));
    Fam.push_back(ElemIndex("s5"));
    break;
  case 2:
    // DIAMOND-II, 6-BA.
    Fam.push_back(ElemIndex("sd1"));
    Fam.push_back(ElemIndex("sd2"));
    Fam.push_back(ElemIndex("sd3"));
    // Fam.push_back(ElemIndex("sd4"));
    Fam.push_back(ElemIndex("sf21"));
    Fam.push_back(ElemIndex("sd31"));
    // Fam.push_back(ElemIndex("sd41"));
    Fam.push_back(ElemIndex("sf1"));
    // Fam.push_back(ElemIndex("sf2"));
    Fam.push_back(ElemIndex("sh1a"));
    Fam.push_back(ElemIndex("sh1e"));
    break;
  default:
    printf("\nmain: unknown case\n");
    exit(1);
    break;
  }

    prt_symm(Fam);
  }

  if (false) {
    Fam.push_back(ElemIndex("q1_2"));
    Fam.push_back(ElemIndex("q1_2m"));
    Fam.push_back(ElemIndex("q1b"));
    Fam.push_back(ElemIndex("q2b"));
    Fam.push_back(ElemIndex("q1d"));
    Fam.push_back(ElemIndex("q2d"));
    Fam.push_back(ElemIndex("q3d"));

    prt_quad(Fam);
  }

  if (true) GetEmittance(ElemIndex("cav"), true);

  if (false) {
    Id.identity();
    Ms[x_] =
      (cos(-2e0*M_PI*globval.TotalTune[Z_])
       +globval.alpha_z*sin(-2e0*M_PI*globval.TotalTune[Z_]))*Id[x_]
      + globval.beta_z*sin(-2e0*M_PI*globval.TotalTune[Z_])*Id[px_];
    Ms[px_] =
      -(1e0+sqr(globval.alpha_z))/globval.beta_z
      *sin(-2e0*M_PI*globval.TotalTune[Z_])*Id[x_]
      +(cos(-2e0*M_PI*globval.TotalTune[Z_])
	-globval.alpha_z*sin(-2e0*M_PI*globval.TotalTune[Z_]))*Id[px_];
    prt_lin_map(1, Ms);
    Ms = exp(-Cell[globval.Cell_nLoc].S/(c0*globval.tau[Z_]))*Ms;
    prt_lin_map(1, Ms);
    printf("\nDet = %17.10e",
	   Ms[x_][x_]*Ms[px_][px_]-Ms[x_][px_]*Ms[px_][x_]);

    // Variables needs to be changed too.
    // putlinmat(6, globval.Ascr, Ascr);
    // a = Ascr[delta_]; Ascr[delta_] = Ascr[ct_]; Ascr[ct_] = a;
    // Ascr = get_A_CS(3, Ascr, dnu);
    // A_Atp = Ascr*tp_S(3, Ascr);
    // printf("\nA_Atp\n");
    // prt_lin_map(3, A_Atp);

    exit(0);
  }

  if (false) {
    // b2_fam[0] = ElemIndex(q_fam[0].c_str());
    // b2_fam[1] = ElemIndex(q_fam[1].c_str());
    // FitTune(b2_fam[0], b2_fam[1], nu[X_], nu[Y_]);
    // get_bn_design_elem(b2_fam[0], 1, Quad, b2[0], a2);
    // get_bn_design_elem(b2_fam[1], 1, Quad, b2[1], a2);

    // printf("\nnu_x = %8.5f nu_y = %8.5f\n",
    // 	   globval.TotalTune[X_], globval.TotalTune[Y_]);
    // printf("  %s = %8.5f  %s = %8.5f\n",
    // 	   q_fam[0].c_str(), b2[0], q_fam[1].c_str(), b2[1]);

    // Ring_GetTwiss(true, 0e0); printglob();
  }

  if (false) {
    // b3_fam[0] = ElemIndex(s_fam[0].c_str());
    // b3_fam[1] = ElemIndex(s_fam[1].c_str());
    // FitChrom(b3[0], b3[1], 0e0, 0e0);
    // get_bn_design_elem(b3_fam[0], 1, Sext, b3[0], a3);
    // get_bn_design_elem(b3_fam[1], 1, Sext, b3[1], a3);
    // get_bnL_design_elem(b3_fam[0], 1, Sext, b3L[0], a3L);
    // get_bnL_design_elem(b3_fam[1], 1, Sext, b3L[1], a3L);

    // printf("\n%s = %10.5f (%10.5f), %s = %10.5f (%10.5f)\n",
    // 	   s_fam[0].c_str(), b3[0], b3L[0], s_fam[1].c_str(), b3[1], b3L[1]);

    // Ring_GetTwiss(true, 0e0); printglob();
  }

  if (false) {
    get_disp();
    exit(0);
  }

  if (true) {
    globval.Cavity_on = true;
    get_dynap(delta, 25, n_turn, false);
  }

}
