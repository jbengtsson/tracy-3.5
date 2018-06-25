#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


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
  MpoleType           *M;
  FILE                *fp;

  const int n_cod_corr = 5;

  Lattice.param.Cavity_on = false;

  for (j = 0; j <= Lattice.param.Cell_nLoc; j++)
    for (k = 0; k < 6; k++) {
      x1[k].push_back(0e0); x2[k].push_back(0e0);
    }
  
  fp = file_write("cod_rms.out");
  
  n_cod = 0;
  for (i = 0; i < n_seed; i++) {
    printf("\norb_corr: seed no %d\n", i+1);

    misalign_rms_type(Dip,  dx, dy, 0e0, true);
    misalign_rms_type(Quad, dx, dy, 0e0, true);
    
    cod = Lattice.orb_corr(n_cod_corr);

    if (cod) {
      n_cod++;

      n = 0;
      for (j = 0; j <= Lattice.param.Cell_nLoc; j++) {
	if (all || (Lattice.Cell[j]->Elem.Kind == Mpole)) {
	  M = static_cast<MpoleType*>(Lattice.Cell[j]);
	  if (M->n_design == Sext) {
	    n++;
	    for (k = 0; k < 6; k++) {
	      x1[k][n-1] += Lattice.Cell[j]->BeamPos[k];
	      x2[k][n-1] += sqr(Lattice.Cell[j]->BeamPos[k]);
	    }
	  }
	}
      }
    } else
      printf("orb_corr: failed\n");

    // Reset orbit trims.
    set_bn_design_fam(Lattice.param.hcorr, Dip, 0e0, 0e0);
    set_bn_design_fam(Lattice.param.vcorr, Dip, 0e0, 0e0);
  }

  printf("\nget_cod_rms: no of seeds %d, no of cods %d\n", n_seed, n_cod);

  n = 0;
  for (j = 0; j <= Lattice.param.Cell_nLoc; j++)
    if (all || (Lattice.Cell[j]->Elem.Kind == Mpole)) {
      M = static_cast<MpoleType*>(Lattice.Cell[j]);
      if (M->n_design == Sext) {
	n++;
	for (k = 0; k < 6; k++) {
	  x_mean[k].push_back(x1[k][n-1]/n_cod);
	  x_sigma[k].push_back(sqrt((n_cod*x2[k][n-1]-sqr(x1[k][n-1]))
				    /(n_cod*(n_cod-1.0))));
	}
	fprintf(fp, "%8.3f %6.2f %10.3e +/- %10.3e %10.3e +/- %10.3e\n",
		Lattice.Cell[j]->S, get_code(Lattice.Cell[j]),
		1e3*x_mean[x_][n-1], 1e3*x_sigma[x_][n-1],
		1e3*x_mean[y_][n-1], 1e3*x_sigma[y_][n-1]);
      }
    } else
      fprintf(fp, "%8.3f %6.2f\n", Lattice.Cell[j]->S,
	      get_code(Lattice.Cell[j]));
  
  fclose(fp);
}


void track(const double Ax, const double Ay)
{
  long int        lastpos;
  int             i;
  ss_vect<double> xt, xs;
  FILE            *fd;

  Lattice.getcod(0e0, lastpos);

  fd = fopen("trackdat_oneturn.dat","w");
  fprintf(fd, "orbit %22.14e %22.14e %22.14e %22.14e %22.14e %22.14e\n",
	  Lattice.Cell[0]->BeamPos[0], Lattice.Cell[0]->BeamPos[1],
	  Lattice.Cell[0]->BeamPos[2], Lattice.Cell[0]->BeamPos[3],
	  Lattice.Cell[0]->BeamPos[4], Lattice.Cell[0]->BeamPos[5] );
  fprintf(fd, "orbit %22.14e %22.14e %22.14e %22.14e %22.14e %22.14e\n",
	  Lattice.param.CODvect[0], Lattice.param.CODvect[1],
	  Lattice.param.CODvect[2], Lattice.param.CODvect[3],
	  Lattice.param.CODvect[4], Lattice.param.CODvect[5]);

  xt.zero(); xt[x_] = Ax; xt[y_] = Ay; 

  fprintf(fd, "start %22.14e %22.14e %22.14e %22.14e %22.14e %22.14e\n",
	  xt[0], xt[1], xt[2], xt[3], xt[4], xt[5] );

  for (i = 0; i <= Lattice.param.Cell_nLoc; i++) {
    Lattice.Cell_Pass(i, i, xt, lastpos);
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
    for (k = 1; k <= Lattice.GetnKid(Fam[j]); k++) {
      loc = Lattice.Elem_GetPos(Fam[j], k);
      if (k % 2 == 0) loc -= 1;
      printf(" %5.1f %6.3f %6.3f %6.3f\n",
	     Lattice.Cell[loc]->S, Lattice.Cell[loc]->Beta[X_],
	     Lattice.Cell[loc]->Beta[Y_], Lattice.Cell[loc]->Eta[X_]);
    }
  }
}


void prt_quad(const std::vector<int> &Fam)
{
  long int loc;
  int      j;

  printf("\n");
  for (j = 0; j < (int)Fam.size(); j++) {
    loc = Lattice.Elem_GetPos(Fam[j], 1);
    printf(" %4.1f %6.3f %6.3f %2d\n",
	   Lattice.Cell[loc]->S, Lattice.Cell[loc]->Beta[X_],
	   Lattice.Cell[loc]->Beta[Y_], Lattice.GetnKid(Fam[j]));
  }
}


void chk_optics(const double alpha_x, const double alpha_y,
		const double beta_x, const double beta_y,
		const double eta_x, const double etap_x,
		const double eta_y, const double etap_y)
{
  Vector2 alpha, beta, eta, etap;

  alpha[X_] = alpha_x; alpha[Y_] = alpha_y;
  beta[X_]  = beta_x;  beta[Y_]  = beta_y;
  eta[X_]   = eta_x;   eta[Y_]   = eta_y;
  etap[X_]  = etap_x;  etap[Y_]  = etap_y;

  Lattice.ttwiss(alpha, beta, eta, etap, 0e0);
}


void chk_mini_beta(const std::vector<int> &Fam)
{
  int    j, k, loc;
  double nu0[] = {0e0, 0e0};

  for (j = 0; j < (int)Fam.size(); j++) {
    printf("\n");
    for (k = 1; k <= Lattice.GetnKid(Fam[j]); k++) {
      loc = Lattice.Elem_GetPos(Fam[j], k);
      if (k % 2 == 1) {
	nu0[X_] = Lattice.Cell[loc]->Nu[X_];
	nu0[Y_] = Lattice.Cell[loc]->Nu[Y_];
      } else {
	loc -= 1;
	printf(" %5.1f %6.3f %6.3f %6.3f %8.5f %8.5f\n",
	       Lattice.Cell[loc]->S, Lattice.Cell[loc]->Beta[X_],
	       Lattice.Cell[loc]->Beta[Y_], Lattice.Cell[loc]->Eta[X_],
	       Lattice.Cell[loc]->Nu[X_]-nu0[X_],
	       Lattice.Cell[loc]->Nu[Y_]-nu0[Y_]);
      }
    }
  }
}


void chk_high_ord_achr(void)
{
  int              k;
  double           dnu[2];
  std::vector<int> loc;

  // D-TBA  0,
  // H-6BA  1,
  // H-8BA  2,
  // RB-6BA 3.
  const int lat_case = 1;

  Lattice.Ring_GetTwiss(true, 0e0);
 
  switch (lat_case) {
  case 0:
    dnu[X_] = 0.0; dnu[Y_] = 0.0;
    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("dr_01"),     1));
    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("idmarker"),  2));
    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("idmarker"),  3));
    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("idmarker"),  4));

    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("dr_01"),     5));
    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("idmarker"),  6));
    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("idmarker"),  7));
    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("idmarker"),  8));

    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("dr_01"),     9));
    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("idmarker"), 10));
    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("idmarker"), 11));
    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("idmarker"), 12));

    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("dr_01"),    13));
    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("idmarker"), 14));
    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("idmarker"), 15));
    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("idmarker"), 16));

    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("dr_01"),    17));
    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("idmarker"), 18));
    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("idmarker"), 19));
    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("idmarker"), 20));

    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("dr_01"),    21));
    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("idmarker"), 22));
    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("idmarker"), 23));
    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("idmarker"), 24));
   break;
  case 1:
    dnu[X_] = 19.0/8.0; dnu[Y_] = 15.0/16.0;
    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("ss1"), 1));
    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("ss1"), 3));
    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("ss1"), 5));
    loc.push_back(Lattice.param.Cell_nLoc);
    break;
  case 2:
    dnu[X_] = 19.0/8.0; dnu[Y_] = 15.0/16.0;
    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("du1"), 1));
    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("du1"), 3));
    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("du1"), 5));
    loc.push_back(Lattice.param.Cell_nLoc);
    break;
  case 3:
    dnu[X_] = 23.0/8.0; dnu[Y_] = 19.0/16.0;
    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("dss1"), 1));
    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("dss1"), 3));
    loc.push_back(Lattice.Elem_GetPos(Lattice.Elem_Index("dss1"), 5));
    loc.push_back(Lattice.param.Cell_nLoc);
    break;
  }

  printf("\nCell phase advance:\n");
  printf("Ideal:    [%7.5f, %7.5f]\n", dnu[X_], dnu[Y_]);
  for (k = 0; k < (int)loc.size(); k++)
    if (k == 0)
      printf("\n %9.5f %8.5f %8.5f %7.5f [%7.5f, %7.5f]\n",
	     Lattice.Cell[loc[k]]->S,
	     Lattice.Cell[loc[k]]->Alpha[X_], Lattice.Cell[loc[k]]->Alpha[Y_],
	     Lattice.Cell[loc[0]]->S,
	     Lattice.Cell[loc[0]]->Nu[X_], Lattice.Cell[loc[0]]->Nu[Y_]);
    else
      printf(" %9.5f %8.5f %8.5f %7.5f [%7.5f, %7.5f]\n",
	     Lattice.Cell[loc[k]]->S,
	     Lattice.Cell[loc[k]]->Alpha[X_], Lattice.Cell[loc[k]]->Alpha[Y_],
	     Lattice.Cell[loc[k]]->S-Lattice.Cell[loc[k-1]]->S, 
	     Lattice.Cell[loc[k]]->Nu[X_]-Lattice.Cell[loc[k-1]]->Nu[X_], 
	     Lattice.Cell[loc[k]]->Nu[Y_]-Lattice.Cell[loc[k-1]]->Nu[Y_]);
}


void chk_mpole_Fam(const int Fnum, const bool exit)
{
  int k, loc;

  printf("\n");
  for (k = 1; k <= Lattice.GetnKid(Fnum); k++) {
    loc = Lattice.Elem_GetPos(Fnum, k);
    if (!exit && ((k-1) % 2 == 1)) loc -= 1;
    printf("%8s %7.3f %8.5f %8.5f\n",
	   Lattice.Cell[loc]->Name, Lattice.Cell[loc]->S,
	   Lattice.Cell[loc]->Beta[X_], Lattice.Cell[loc]->Beta[Y_]);
  }
}


void chk_mpole(void)
{
  int              k;
  std::vector<int> Fnum;

  const int lat_case = 1;

  switch (lat_case) {
  case 1:
    // H-6-BA.
    Fnum.push_back(Lattice.Elem_Index("sf"));
    Fnum.push_back(Lattice.Elem_Index("sda"));
    Fnum.push_back(Lattice.Elem_Index("sdb"));

    Fnum.push_back(Lattice.Elem_Index("s1"));
    Fnum.push_back(Lattice.Elem_Index("s2"));
    Fnum.push_back(Lattice.Elem_Index("s3"));
    Fnum.push_back(Lattice.Elem_Index("s4"));
    Fnum.push_back(Lattice.Elem_Index("s5"));

    Fnum.push_back(Lattice.Elem_Index("o1a"));
    Fnum.push_back(Lattice.Elem_Index("o2a"));
    Fnum.push_back(Lattice.Elem_Index("o3"));
    Fnum.push_back(Lattice.Elem_Index("o1b"));
    Fnum.push_back(Lattice.Elem_Index("o2b"));

    Fnum.push_back(Lattice.Elem_Index("o4"));
    Fnum.push_back(Lattice.Elem_Index("o5"));
    Fnum.push_back(Lattice.Elem_Index("o6"));
    break;
  case 2:
    // H-8-BA_II.
    Fnum.push_back(Lattice.Elem_Index("sf"));
    Fnum.push_back(Lattice.Elem_Index("sd"));
    Fnum.push_back(Lattice.Elem_Index("s1"));
    Fnum.push_back(Lattice.Elem_Index("s2"));
    Fnum.push_back(Lattice.Elem_Index("s3"));
    Fnum.push_back(Lattice.Elem_Index("s4"));
    Fnum.push_back(Lattice.Elem_Index("s5"));
    Fnum.push_back(Lattice.Elem_Index("s6"));
    break;
  case 3:
    // RB-6-BA.
    Fnum.push_back(Lattice.Elem_Index("sd"));
    Fnum.push_back(Lattice.Elem_Index("sfm"));
    Fnum.push_back(Lattice.Elem_Index("sdm"));
    Fnum.push_back(Lattice.Elem_Index("sxx"));
    Fnum.push_back(Lattice.Elem_Index("sxy1"));
    Fnum.push_back(Lattice.Elem_Index("sxy2"));
    Fnum.push_back(Lattice.Elem_Index("sxy3"));
    Fnum.push_back(Lattice.Elem_Index("syy1"));
    Fnum.push_back(Lattice.Elem_Index("syy2"));
    Fnum.push_back(Lattice.Elem_Index("syy3"));
    break;
  }

  Lattice.Ring_GetTwiss(true, 0e0);
 
  printf("\nSextupole Scheme:\n");
  for (k = 0; k < (int)Fnum.size(); k++)
    chk_mpole_Fam(Fnum[k], false);
}


void chk_dip(void)
{
  int       k;
  double    L, phi, L_sum, phi_sum;
  MpoleType *M;

  Lattice.Ring_GetTwiss(true, 0e0);
 
  printf("\nLong grad dipole:\n");
  L_sum = 0e0; phi_sum = 0e0;
  for (k = 0; k <= Lattice.param.Cell_nLoc; k++) {
    if (Lattice.Cell[k]->Elem.Kind == Mpole) {
      M = static_cast<MpoleType*>(Lattice.Cell[k]);
      if (M->irho != 0e0) {
	L = Lattice.Cell[k]->L;
	phi = L*M->irho*180e0/M_PI;
	L_sum += L; phi_sum += phi;
	printf(" %6s %4.3f %7.3f %9.6f %9.6f %9.6f %9.6f %9.6f\n",
	       Lattice.Cell[k]->Name, L, 1e0/M->irho, phi, M->Tx1, M->Tx2,
	       L_sum, phi_sum);
      }
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
    Fnum.push_back(Lattice.Elem_Index("sf"));
    // Fnum.push_back(Lattice.Elem_Index("sda"));
    // Fnum.push_back(Lattice.Elem_Index("sdb"));
    break;
  }

  Lattice.Ring_GetTwiss(true, 0e0);
 
  printf("\nMultipole Phase Advances:\n");
  for (n = 1; n <= Lattice.GetnKid(Fnum[0]); n++) {
    loc1 = Lattice.Elem_GetPos(Fnum[0], n);
    if (n == 1)
      printf("%10s %7.5f %7.5f\n",
	     Lattice.Cell[loc1]->Name,
	     Lattice.Cell[loc1]->Nu[X_], Lattice.Cell[loc1]->Nu[Y_]);
    else
      printf("%10s %7.5f %7.5f\n",
	     Lattice.Cell[loc1]->Name,
	     Lattice.Cell[loc1]->Nu[X_]-Lattice.Cell[loc0]->Nu[X_],
	     Lattice.Cell[loc1]->Nu[Y_]-Lattice.Cell[loc0]->Nu[Y_]);
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
  int        k, n;
  double     phi, b1, bn, an;
  MpoleType  *M;

  const double Brho = Lattice.param.Energy*1e9/c0;

  for (k = 0; k <= Lattice.param.Cell_nLoc; k++)
    if (Lattice.Cell[k]->Elem.Kind == Mpole) {
      M = static_cast<MpoleType*>(Lattice.Cell[k]);
      switch (M->n_design) {
      case Dip:
	phi = Lattice.Cell[k]->L*M->irho; b1 = M->irho;	n = Quad;
	get_bn_design_elem(Lattice.Cell[k]->Fnum, 1, n, bn, an);
	printf("%10s L = %5.3f b_%1d = %7.3f     B^ = %6.3f phi = %6.3f\n",
	       Lattice.Cell[k]->Name, Lattice.Cell[k]->L,
	       1, b1, Brho*b1, phi*180e0/M_PI);
	printf("                          b_%1d = %7.3f     B^ = %6.3f\n",
	       n, bn, get_pole_tip_field(Brho, R_ref, n, bn));
	break;
      case Quad:
	n = M->n_design;
	get_bn_design_elem(Lattice.Cell[k]->Fnum, 1, n, bn, an);
	printf("%10s L = %5.3f b_%1d = %7.3f     B^ = %6.3f\n",
	       Lattice.Cell[k]->Name, Lattice.Cell[k]->L, n, bn,
	       get_pole_tip_field(Brho, R_ref, n, bn));
	break;
      case Sext:
	n = M->n_design;
	get_bn_design_elem(Lattice.Cell[k]->Fnum, 1, n, bn, an);
	printf("%10s L = %5.3f b_%1d = %11.3e B^ = %6.3f\n",
	       Lattice.Cell[k]->Name, Lattice.Cell[k]->L, n, bn,
	       get_pole_tip_field(Brho, R_ref, n, bn));
	break;
      case Dodec:
	n = M->n_design - 2;
	get_bn_design_elem(Lattice.Cell[k]->Fnum, 1, n, bn, an);
	printf("%10s L = %5.3f b_%1d = %11.3e B^ = %6.3f\n",
	       Lattice.Cell[k]->Name, Lattice.Cell[k]->L, n, bn,
	       get_pole_tip_field(Brho, R_ref, n, bn));
	n = M->n_design;
	get_bn_design_elem(Lattice.Cell[k]->Fnum, 1, n, bn, an);
	printf("%10s L = %5.3f b_%1d = %11.3e B^ = %6.3f\n",
	       Lattice.Cell[k]->Name, Lattice.Cell[k]->L, n, bn,
	       get_pole_tip_field(Brho, R_ref, n, bn));
	break;
      }
    }
}


int main(int argc, char *argv[])
{
  bool             tweak;
  long int         lastn, lastpos, loc;
  int              b2_fam[2], b3_fam[2];
  double           b2[2], a2, b3[2], b3L[2], a3, a3L, f_rf, dx;
  Matrix           M;
  std::vector<int> Fam;
  CavityType       *C;
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

  Lattice.param.H_exact        = false; Lattice.param.quad_fringe  = false;
  Lattice.param.Cavity_on      = false; Lattice.param.radiation    = false;
  Lattice.param.emittance      = false; Lattice.param.IBS          = false;
  Lattice.param.pathlength     = false; Lattice.param.bpm          = 0;
  Lattice.param.dip_edge_fudge = true;  Lattice.param.reverse_elem = !false;

  // 1: DIAMOND, 3: Oleg I, 4: Oleg II.
  FieldMap_filetype = 1; sympl = false;

  trace = !true;

  if (true)
     Lattice.Read_Lattice(argv[1]);
  else
     Lattice.rdmfile(argv[1]);

  if (false) {
    loc = Lattice.Elem_GetPos(Lattice.Elem_Index("bb"), 1);
    map.identity();
    // Tweak to remain within field map range at entrance.
    tweak = true;
    if (tweak) {
      dx = -1.4e-3; map[x_] += dx;
    }
    Lattice.Cell_Pass(loc, loc, map, lastpos);
    if (tweak) map[x_] -= dx;
    prt_lin_map(3, map);
    getlinmat(6, map, M);
    printf("\n1-Det: %9.3e\n", 1e0-DetMat(6, M));
    exit(0);
  }

  if (false) {
    chk_high_ord_achr();
    exit(0);
  }

  if (false) {
    chk_mpole();
    Lattice.prt_lat("linlat1.out", Lattice.param.bpm, true);
    Lattice.prt_lat("linlat.out", Lattice.param.bpm, true, 10);
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
    chk_optics(0.0, 0.0, 9.21924, 3.17296, 0.0, 0.0, 0.0, 0.0);
    Lattice.prt_lat("linlat1.out", Lattice.param.bpm, true);
    Lattice.prt_lat("linlat.out", Lattice.param.bpm, true, 10);
    exit(0);
  }

  if (false) {
    Lattice.param.Cavity_on = true;
    track(6e-3, 0.1e-3);
    exit(0);
  }

  if (false) {
    Lattice.Ring_GetTwiss(true, 0e0); printglob();
    Lattice.set_map("ps_rot", 0.55/6.0, -0.07/6.0);
    Lattice.Ring_GetTwiss(true, 0e0); printglob();
    exit(0);
  }

  if (false) {
    dnu_mpole();
    exit(0);
  }

  if (false) no_sxt();

  Lattice.param.Cavity_on = false; Lattice.param.radiation = false;
  Lattice.Ring_GetTwiss(true, 0e0); printglob();

  if (false) Lattice.get_alphac2();

  Lattice.prtmfile("flat_file.dat");

  // Lattice.param.bpm = Lattice.Elem_Index("bpm");
  Lattice.prt_lat("linlat1.out", Lattice.param.bpm, true);
  Lattice.prt_lat("linlat.out", Lattice.param.bpm, true, 10);
  Lattice.prt_chrom_lat();

  if (false) {
    iniranf(seed); setrancut(1e0);
    // Lattice.param.bpm = Lattice.Elem_Index("mon");
    Lattice.param.bpm = Lattice.Elem_Index("bpm");
    Lattice.param.hcorr = Lattice.Elem_Index("ch");
    Lattice.param.vcorr = Lattice.Elem_Index("cv");
    // ALS-U.
    // Lattice.param.bpm = Lattice.Elem_Index("bpm_m");
    // Lattice.param.hcorr = Lattice.Elem_Index("corr_h");
    // Lattice.param.vcorr = Lattice.Elem_Index("corr_v");

    gcmat(Lattice.param.bpm, Lattice.param.hcorr, 1);
    gcmat(Lattice.param.bpm, Lattice.param.vcorr, 2);

    get_cod_rms(50e-6, 50e-6, 100, true);

    exit(0);
  }

  if (false) {
    Fam.push_back(Lattice.Elem_Index("ts1b"));
    // Fam.push_back(Lattice.Elem_Index("ts1d"));
    chk_mini_beta(Fam);
    exit(0);
  }

  if (false) {
    // Fam.push_back(Lattice.Elem_Index("s1b"));
    // Fam.push_back(Lattice.Elem_Index("s1d"));
    // Fam.push_back(Lattice.Elem_Index("s2b"));
    // Fam.push_back(Lattice.Elem_Index("s2d"));
    // Fam.push_back(Lattice.Elem_Index("sx1"));
    // Fam.push_back(Lattice.Elem_Index("sy1"));

  switch (2) {
  case 1:
    // DIAMOND.
    Fam.push_back(Lattice.Elem_Index("ts1a"));
    Fam.push_back(Lattice.Elem_Index("ts1ab"));
    Fam.push_back(Lattice.Elem_Index("ts2a"));
    Fam.push_back(Lattice.Elem_Index("ts2ab"));
    Fam.push_back(Lattice.Elem_Index("ts1b"));
    Fam.push_back(Lattice.Elem_Index("ts2b"));
    Fam.push_back(Lattice.Elem_Index("ts1c"));
    Fam.push_back(Lattice.Elem_Index("ts2c"));
    Fam.push_back(Lattice.Elem_Index("ts1d"));
    Fam.push_back(Lattice.Elem_Index("ts2d"));
    Fam.push_back(Lattice.Elem_Index("ts1e"));
    Fam.push_back(Lattice.Elem_Index("ts2e"));

    Fam.push_back(Lattice.Elem_Index("s1"));
    Fam.push_back(Lattice.Elem_Index("s2"));
    Fam.push_back(Lattice.Elem_Index("s3"));
    Fam.push_back(Lattice.Elem_Index("s4"));
    Fam.push_back(Lattice.Elem_Index("s5"));
    break;
  case 2:
    // DIAMOND-II, 6-BA.
    Fam.push_back(Lattice.Elem_Index("sd1"));
    Fam.push_back(Lattice.Elem_Index("sd2"));
    Fam.push_back(Lattice.Elem_Index("sd3"));
    // Fam.push_back(Lattice.Elem_Index("sd4"));
    Fam.push_back(Lattice.Elem_Index("sf21"));
    Fam.push_back(Lattice.Elem_Index("sd31"));
    // Fam.push_back(Lattice.Elem_Index("sd41"));
    Fam.push_back(Lattice.Elem_Index("sf1"));
    // Fam.push_back(Lattice.Elem_Index("sf2"));
    Fam.push_back(Lattice.Elem_Index("sh1a"));
    Fam.push_back(Lattice.Elem_Index("sh1e"));
    break;
  }

    prt_symm(Fam);
  }

  if (false) {
    Fam.push_back(Lattice.Elem_Index("q1_2"));
    Fam.push_back(Lattice.Elem_Index("q1_2m"));
    Fam.push_back(Lattice.Elem_Index("q1b"));
    Fam.push_back(Lattice.Elem_Index("q2b"));
    Fam.push_back(Lattice.Elem_Index("q1d"));
    Fam.push_back(Lattice.Elem_Index("q2d"));
    Fam.push_back(Lattice.Elem_Index("q3d"));

    prt_quad(Fam);
  }

  if (true) Lattice.GetEmittance(Lattice.Elem_Index("cav"), true);

  if (false) {
    b2_fam[0] = Lattice.Elem_Index(q_fam[0].c_str());
    b2_fam[1] = Lattice.Elem_Index(q_fam[1].c_str());
    Lattice.FitTune(b2_fam[0], b2_fam[1], nu[X_], nu[Y_]);
    get_bn_design_elem(b2_fam[0], 1, Quad, b2[0], a2);
    get_bn_design_elem(b2_fam[1], 1, Quad, b2[1], a2);

    printf("\nnu_x = %8.5f nu_y = %8.5f\n",
	   Lattice.param.TotalTune[X_], Lattice.param.TotalTune[Y_]);
    printf("  %s = %8.5f  %s = %8.5f\n",
	   q_fam[0].c_str(), b2[0], q_fam[1].c_str(), b2[1]);

    Lattice.Ring_GetTwiss(true, 0e0); printglob();
  }

  if (false) {
    b3_fam[0] = Lattice.Elem_Index(s_fam[0].c_str());
    b3_fam[1] = Lattice.Elem_Index(s_fam[1].c_str());
    Lattice.FitChrom(b3[0], b3[1], 0e0, 0e0);
    get_bn_design_elem(b3_fam[0], 1, Sext, b3[0], a3);
    get_bn_design_elem(b3_fam[1], 1, Sext, b3[1], a3);
    get_bnL_design_elem(b3_fam[0], 1, Sext, b3L[0], a3L);
    get_bnL_design_elem(b3_fam[1], 1, Sext, b3L[1], a3L);

    printf("\n%s = %10.5f (%10.5f), %s = %10.5f (%10.5f)\n",
	   s_fam[0].c_str(), b3[0], b3L[0], s_fam[1].c_str(), b3[1], b3L[1]);

    Lattice.Ring_GetTwiss(true, 0e0); printglob();
  }

  Lattice.prtmfile("flat_file.fit");

  if (false) {
    Lattice.param.Cavity_on = false; Lattice.param.radiation = false;

    C = static_cast<CavityType*>
      (Lattice.Cell[Lattice.Elem_GetPos(Lattice.Elem_Index("cav"), 1)]);

    f_rf = C->freq;
    printf("\nf_rf = %10.3e\n", f_rf);

    Lattice.param.Cavity_on = true;
    // Synchro-betatron resonance for "101pm_above_coupres_tracy.lat".
    // track("track.out", 2.6e-3, 0e0, 1e-6, 0e0, 0e0, n_turn, lastn, lastpos,
    // 	  0, 0*f_rf);
    // track("track.out", 1e-6, 0e0, 1.9e-3, 0e0, 0e0, n_turn, lastn, lastpos,
    // 	  0, 0*f_rf);
    
    // track("track.out", 1e-3, 0e0, 1e-3, 0e0, 0e0, 10*n_turn, lastn, lastpos,
    // 	  0, f_rf);

    // lattice/101pm_s7o7_a_tracy.lat.
    Lattice.track("track.out", 5.5e-3, 0e0, 1e-6, 0e0, 0e0, n_turn, lastn,
		  lastpos, 0, 0*f_rf);
  }

  if (true) {
    // set_map("M", 0.0, 0.0);

    Lattice.param.Cavity_on = true;
    Lattice.get_dynap(delta, 25, n_turn, false);
  }

}
