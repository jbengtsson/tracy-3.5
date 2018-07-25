#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


const double dnu[] = {0.03, 0.02};


void prt_name(FILE * outf, const char *name)
{
    int j, k, len;

    len = strlen(name);

    j = 0;
    do {
	fprintf(outf, "%c", name[j]);
	j++;
    } while ((j < len) && (name[j] != ' '));
    fprintf(outf, ",");
    for (k = j; k < len; k++)
	fprintf(outf, "%c", name[k]);
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


void chk_high_ord_achr(void)
{
  int              k;
  double           dnu[2];
  std::vector<int> loc;

  // M-6HBAi 1,
  // H-6BA   2,
  // H-8BA   3,
  // RB-6BA  4.
  const int lat_case = 1;

  Ring_GetTwiss(true, 0e0);
 
  switch (lat_case) {
  case 1:
    dnu[X_] = 19.0/8.0; dnu[Y_] = 15.0/16.0;
    loc.push_back(Elem_GetPos(ElemIndex("dc_1_01"),  1));
    loc.push_back(Elem_GetPos(ElemIndex("idmarker"), 2));
    loc.push_back(Elem_GetPos(ElemIndex("idmarker"), 3));
    loc.push_back(Elem_GetPos(ElemIndex("idmarker"), 4));
   break;
  case 2:
    dnu[X_] = 19.0/8.0; dnu[Y_] = 15.0/16.0;
    loc.push_back(Elem_GetPos(ElemIndex("ss1"), 1));
    loc.push_back(Elem_GetPos(ElemIndex("ss1"), 3));
    loc.push_back(Elem_GetPos(ElemIndex("ss1"), 5));
    loc.push_back(globval.Cell_nLoc);
    break;
  case 3:
    dnu[X_] = 19.0/8.0; dnu[Y_] = 15.0/16.0;
    loc.push_back(Elem_GetPos(ElemIndex("du1"), 1));
    loc.push_back(Elem_GetPos(ElemIndex("du1"), 3));
    loc.push_back(Elem_GetPos(ElemIndex("du1"), 5));
    loc.push_back(globval.Cell_nLoc);
    break;
  case 4:
    dnu[X_] = 23.0/8.0; dnu[Y_] = 19.0/16.0;
    loc.push_back(Elem_GetPos(ElemIndex("dss1"), 1));
    loc.push_back(Elem_GetPos(ElemIndex("dss1"), 3));
    loc.push_back(Elem_GetPos(ElemIndex("dss1"), 5));
    loc.push_back(globval.Cell_nLoc);
    break;
  }

  printf("\nCell phase advance:\n");
  printf("Ideal:    [%7.5f, %7.5f]\n", dnu[X_], dnu[Y_]);
  for (k = 0; k < (int)loc.size(); k++)
    if (k == 0)
      printf("\n %9.5f %8.5f %8.5f %7.5f [%7.5f, %7.5f]\n",
	     Cell[loc[k]].S, Cell[loc[k]].Alpha[X_], Cell[loc[k]].Alpha[Y_],
	     Cell[loc[0]].S,  Cell[loc[0]].Nu[X_], Cell[loc[0]].Nu[Y_]);
    else
      printf(" %9.5f %8.5f %8.5f %7.5f [%7.5f, %7.5f]\n",
	     Cell[loc[k]].S, Cell[loc[k]].Alpha[X_], Cell[loc[k]].Alpha[Y_],
	     Cell[loc[k]].S-Cell[loc[k-1]].S, 
	     Cell[loc[k]].Nu[X_]-Cell[loc[k-1]].Nu[X_], 
	     Cell[loc[k]].Nu[Y_]-Cell[loc[k-1]].Nu[Y_]);
}


void chk_mI_trans(void)
{
  int Fnum, k, loc0, loc1;

  // M-6HBAi 1,
  const int lat_case = 1;

  Ring_GetTwiss(true, 0e0);
 
  switch (lat_case) {
  case 1:
    Fnum = ElemIndex("sextmark");
   break;
  case 2:
    Fnum = ElemIndex("sf");
   break;
  }

  printf("\nChromatic sextupole phase advance:\n");
  for (k = 2; k <= GetnKid(Fnum); k += 2) {
    loc0 = Elem_GetPos(Fnum, k-1); loc1 = Elem_GetPos(Fnum, k);
    printf(" %8s %7.3f [%7.5f, %7.5f]\n",
	   Cell[loc1].Elem.PName, Cell[loc1].S,
	   Cell[loc1].Nu[X_]-Cell[loc0].Nu[X_], 
	   Cell[loc1].Nu[Y_]-Cell[loc0].Nu[Y_]);
  }
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

  // SLS-2   1,
  // M-6HBAi 2,
  const int lat_case = 2;

  switch (lat_case) {
  case 1:
    // SLS-2.
    Fnum.push_back(ElemIndex("sd"));
    Fnum.push_back(ElemIndex("sf"));
    // Fnum.push_back(ElemIndex("sdm"));
    // Fnum.push_back(ElemIndex("sfm"));

    Fnum.push_back(ElemIndex("sxx"));
    Fnum.push_back(ElemIndex("sxy"));
    Fnum.push_back(ElemIndex("syy"));

    Fnum.push_back(ElemIndex("ocxm"));
    // Fnum.push_back(ElemIndex("ocx1"));
    // Fnum.push_back(ElemIndex("ocx2"));

    Fnum.push_back(ElemIndex("ocym"));
    // Fnum.push_back(ElemIndex("ocy1"));
    // Fnum.push_back(ElemIndex("ocy2"));

    Fnum.push_back(ElemIndex("oxx"));
    Fnum.push_back(ElemIndex("oxy"));
    Fnum.push_back(ElemIndex("oyy"));
    break;
  case 2:
    // M-6HBAi.
    Fnum.push_back(ElemIndex("sextmark"));

    Fnum.push_back(ElemIndex("sd1"));
    Fnum.push_back(ElemIndex("sd2"));
    // Fnum.push_back(ElemIndex("sd1a"));
    // Fnum.push_back(ElemIndex("sd2a"));
    // Fnum.push_back(ElemIndex("sd1b"));
    // Fnum.push_back(ElemIndex("sd2b"));
    break;
  case 3:
    // H-6-BA.
    Fnum.push_back(ElemIndex("sf"));
    Fnum.push_back(ElemIndex("sda"));
    Fnum.push_back(ElemIndex("sdb"));

    // Fnum.push_back(ElemIndex("s1"));
    // Fnum.push_back(ElemIndex("s2"));
    // Fnum.push_back(ElemIndex("s3"));
    // Fnum.push_back(ElemIndex("s4"));
    // Fnum.push_back(ElemIndex("s5"));

    Fnum.push_back(ElemIndex("o1a"));
    Fnum.push_back(ElemIndex("o2a"));
    Fnum.push_back(ElemIndex("o1b"));
    Fnum.push_back(ElemIndex("o2b"));
    Fnum.push_back(ElemIndex("o3"));

    Fnum.push_back(ElemIndex("o4"));
    Fnum.push_back(ElemIndex("o5"));
    Fnum.push_back(ElemIndex("o6"));
    break;
  case 4:
    // H-8-BA_II.
    Fnum.push_back(ElemIndex("sf")); Fnum.push_back(ElemIndex("sd"));
    Fnum.push_back(ElemIndex("s1")); Fnum.push_back(ElemIndex("s2"));
    Fnum.push_back(ElemIndex("s3")); Fnum.push_back(ElemIndex("s4"));
    Fnum.push_back(ElemIndex("s5")); Fnum.push_back(ElemIndex("s6"));
    break;
  case 5:
    // RB-6-BA.
    Fnum.push_back(ElemIndex("sd"));  Fnum.push_back(ElemIndex("sfm"));
    Fnum.push_back(ElemIndex("sdm")); Fnum.push_back(ElemIndex("sxx"));
    Fnum.push_back(ElemIndex("sxy1"));
    Fnum.push_back(ElemIndex("sxy2"));
    Fnum.push_back(ElemIndex("sxy3"));
    Fnum.push_back(ElemIndex("syy1"));
    Fnum.push_back(ElemIndex("syy2"));
    Fnum.push_back(ElemIndex("syy3"));
    break;
  case 6:
    // ALS-U.
    Fnum.push_back(ElemIndex("sf1"));
    Fnum.push_back(ElemIndex("sf2"));
    Fnum.push_back(ElemIndex("sf3"));
    Fnum.push_back(ElemIndex("sd1h"));
    Fnum.push_back(ElemIndex("sd2h"));
    Fnum.push_back(ElemIndex("sd3h"));
    Fnum.push_back(ElemIndex("sd4h"));
    Fnum.push_back(ElemIndex("sh1"));
    Fnum.push_back(ElemIndex("sh2"));
    Fnum.push_back(ElemIndex("sh3"));
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
      }
    }
}


void get_dbeta_deta(const double delta)
{
  int            j, k;
  vector<double> dbeta[2], deta[2];
  FILE           *outf;

  const string file_name = "dbeta_deta.out";

  outf = file_write(file_name.c_str());

  printf("\nOptics for delta = %4.2f%%\n", 1e2*delta);
  Ring_GetTwiss(true, delta); printglob();
  for (k = 0; k <= globval.Cell_nLoc; k++) {
    for (j = 0; j < 2; j++)
      dbeta[j].push_back(Cell[k].Beta[j]);
    deta[0].push_back(Cell[k].Eta[X_]); deta[1].push_back(Cell[k].Etap[X_]);
  }
  printf("\nOptics for delta = %4.2f%%\n", 0e0);
  Ring_GetTwiss(true, 0e0); printglob();
  for (k = 0; k <= globval.Cell_nLoc; k++) {
    fprintf(outf, "%4d %10s %8.3f %4.1f %10.5f %10.5f %8.5f %8.5f"
	    " %10.5f %10.5f %8.5f %8.5f\n",
	    k, Cell[k].Elem.PName, Cell[k].S, get_code(Cell[k]),
	    Cell[k].Beta[X_], Cell[k].Beta[Y_],
	    Cell[k].Eta[X_], Cell[k].Etap[X_],
	    dbeta[X_][k], dbeta[Y_][k], deta[0][k], deta[1][k]);
  }

  fclose(outf);
}


int main(int argc, char *argv[])
{
  bool             tweak;
  long int         lastn, lastpos, loc;
  int              b2_fam[2], b3_fam[2];
  double           b2[2], a2, b3[2], b3L[2], a3, a3L, f_rf, dx;
  Matrix           M;
  std::vector<int> Fam;
  ostringstream    str;

  const long        seed   = 1121;
  const int         n_turn = 2064;
  const double      delta  = 3e-2,
  //                   nu[]    = { 102.18/20.0, 68.30/20.0 };
  // const std::string q_fam[] = { "qfe", "qde" }, s_fam[] = { "sfh", "sd" };
  //                   nu[]    = { 39.1/12.0, 15.25/12.0 };
  //                   // nu[]    = { 3.266+0.01, 1.275 };
  // const std::string q_fam[] = { "qm2b", "qm3" }, s_fam[] = { "sfh", "sd" };
                    nu[]    = { 9.2, 3.64 };
  const std::string q_fam[] = { "qf03", "qd04" }, s_fam[] = { "sfh", "sd" };

  const double R_ref = 5e-3;

  // 1: DIAMOND, 2: NSLS-II, 3: Oleg I, 4: Oleg II.
  FieldMap_filetype = 4; sympl = !false;

  reverse_elem = false;

  trace = !true;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  globval.H_exact        = false;  globval.quad_fringe = false;
  globval.Cavity_on      = false;  globval.radiation   = false;
  globval.emittance      = false;  globval.IBS         = false;
  globval.pathlength     = false;  globval.bpm         = 0;
  globval.dip_edge_fudge = true;

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
    prt_drift();
    exit(0);
  }

  if (false) {
    chk_high_ord_achr();
    exit(0);
  }

  if (false) {
    chk_mI_trans();
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
    chk_optics(-0.3026596977, 11.1525671, -0.1377651491, 0.3779633018,
	       0.7946624468, -0.02521684799, 0.0, 0.0);
    prt_lat("linlat1.out", globval.bpm, true);
    prt_lat("linlat.out", globval.bpm, true, 10);
    exit(0);
  }

  if (false) {
    globval.Cavity_on = true;
    track(6e-3, 0.1e-3);
    exit(0);
  }

  if (false) {
    dnu_mpole();
    exit(0);
  }

  if (!false) {
    Ring_GetTwiss(true, 0e0); printglob();
    set_map("ps_rot", dnu[X_], dnu[Y_]);
    Ring_GetTwiss(true, 0e0); printglob();
  }

  if (false) {
    get_dbeta_deta(1e-2);
    exit(0);
  }

  if (false) no_sxt();

  globval.Cavity_on = false; globval.radiation = false;
  Ring_GetTwiss(true, 0e0); printglob();

  if (false) get_alphac2();

  prtmfile("flat_file.dat");

  // globval.bpm = ElemIndex("bpm");
  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);
  prt_chrom_lat();

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

  if (true)  GetEmittance(ElemIndex("cav"), true);

  if (false) {
    b2_fam[0] = ElemIndex(q_fam[0].c_str());
    b2_fam[1] = ElemIndex(q_fam[1].c_str());
    FitTune(b2_fam[0], b2_fam[1], nu[X_], nu[Y_]);
    get_bn_design_elem(b2_fam[0], 1, Quad, b2[0], a2);
    get_bn_design_elem(b2_fam[1], 1, Quad, b2[1], a2);

    printf("\nnu_x = %8.5f nu_y = %8.5f\n",
	   globval.TotalTune[X_], globval.TotalTune[Y_]);
    printf("  %s = %8.5f  %s = %8.5f\n",
	   q_fam[0].c_str(), b2[0], q_fam[1].c_str(), b2[1]);

    Ring_GetTwiss(true, 0e0); printglob();
  }

  if (false) {
    b3_fam[0] = ElemIndex(s_fam[0].c_str());
    b3_fam[1] = ElemIndex(s_fam[1].c_str());
    FitChrom(b3[0], b3[1], 0e0, 0e0);
    get_bn_design_elem(b3_fam[0], 1, Sext, b3[0], a3);
    get_bn_design_elem(b3_fam[1], 1, Sext, b3[1], a3);
    get_bnL_design_elem(b3_fam[0], 1, Sext, b3L[0], a3L);
    get_bnL_design_elem(b3_fam[1], 1, Sext, b3L[1], a3L);

    printf("\n%s = %10.5f (%10.5f), %s = %10.5f (%10.5f)\n",
	   s_fam[0].c_str(), b3[0], b3L[0], s_fam[1].c_str(), b3[1], b3L[1]);

    Ring_GetTwiss(true, 0e0); printglob();
  }

  prtmfile("flat_file.fit");

  if (false) {
    globval.Cavity_on = false; globval.radiation = false;

    f_rf = Cell[Elem_GetPos(ElemIndex("cav"), 1)].Elem.C->Pfreq;
    printf("\nf_rf = %10.3e\n", f_rf);

    globval.Cavity_on = true;
    // Synchro-betatron resonance for "101pm_above_coupres_tracy.lat".
    // track("track.out", 2.6e-3, 0e0, 1e-6, 0e0, 0e0, n_turn, lastn, lastpos,
    // 	  0, 0*f_rf);
    // track("track.out", 1e-6, 0e0, 1.9e-3, 0e0, 0e0, n_turn, lastn, lastpos,
    // 	  0, 0*f_rf);
    
    // track("track.out", 1e-3, 0e0, 1e-3, 0e0, 0e0, 10*n_turn, lastn, lastpos,
    // 	  0, f_rf);

    // lattice/101pm_s7o7_a_tracy.lat.
    track("track.out", 5.5e-3, 0e0, 1e-6, 0e0, 0e0, n_turn, lastn, lastpos,
    	  0, 0*f_rf);
  }

  if (true) {
    globval.Cavity_on = true;
    get_dynap(delta, 25, n_turn, false);
  }

}
