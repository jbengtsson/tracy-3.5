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

  ttwiss(alpha, beta, eta, etap, 0e0);
}


void chk_mini_beta(const std::vector<int> &Fam)
{
  int    j, k, loc;
  double nu0[2];

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


void chk_high_oord_achr(void )
{
  int loc[4];

  Ring_GetTwiss(true, 0e0);
 
  loc[0] = Elem_GetPos(ElemIndex("idmarker_end"), 1);
  loc[1] = Elem_GetPos(ElemIndex("ss1"), 1);
  loc[2] = Elem_GetPos(ElemIndex("idmarker"), 2);
  loc[3] = globval.Cell_nLoc;

  printf("\nCell phase advance:\n");
  printf("Ideal:  [%7.5f, %7.5f]\n", 19.0/8.0, 15.0/16.0);
  printf("\n        [%7.5f, %7.5f]\n",
	 Cell[loc[0]].Nu[X_], Cell[loc[0]].Nu[Y_]);
  printf("        [%7.5f, %7.5f]\n",
	 Cell[loc[1]].Nu[X_]-Cell[loc[0]].Nu[X_], 
	 Cell[loc[1]].Nu[Y_]-Cell[loc[0]].Nu[Y_]);
  printf("        [%7.5f, %7.5f]\n",
	 Cell[loc[2]].Nu[X_]-Cell[loc[1]].Nu[X_], 
	 Cell[loc[2]].Nu[Y_]-Cell[loc[1]].Nu[Y_]);
  printf("        [%7.5f, %7.5f]\n",
	 Cell[loc[3]].Nu[X_]-Cell[loc[2]].Nu[X_], 
	 Cell[loc[3]].Nu[Y_]-Cell[loc[2]].Nu[Y_]);
}


void chk_b3_Fam(const int Fnum, const bool exit)
{
  int k, loc;

  printf("\n");
  for (k = 1; k <= GetnKid(Fnum); k++) {
    loc = Elem_GetPos(Fnum, k);
    if (!exit && ((k-1) % 2 == 1)) loc -= 1;
    printf("%8s %7.5f %7.5f\n",
	   Cell[loc].Elem.PName, Cell[loc].Beta[X_], Cell[loc].Beta[Y_]);
  }
}


void chk_b3(void )
{
  int k;

  const int n_b3 = 9,
    Fnum[] =
      {ElemIndex("sfa"), ElemIndex("sfb"),
       ElemIndex("sda"), ElemIndex("sdb"),
       ElemIndex("s1"),
       ElemIndex("s2a"), ElemIndex("s2b"),
       ElemIndex("s3"),  ElemIndex("s4"), ElemIndex("s5"), ElemIndex("s6")};

  Ring_GetTwiss(true, 0e0);
 
  printf("\nSextupole Scheme:\n");
  for (k = 0; k < n_b3; k++)
    if (k == 2)
      chk_b3_Fam(Fnum[k], false);
    else
      chk_b3_Fam(Fnum[k], true);
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
  const double      delta  = 3.0e-2,
  //                   nu[]    = { 102.18/20.0, 68.30/20.0 };
  // const std::string q_fam[] = { "qfe", "qde" }, s_fam[] = { "sfh",  "sd" };
  //                   nu[]    = { 39.1/12.0, 15.25/12.0 };
  //                   // nu[]    = { 3.266+0.01, 1.275 };
  // const std::string q_fam[] = { "qm2b", "qm3" }, s_fam[] = { "sfh",  "sd" };
                    nu[]    = { 9.2, 3.64 };
  const std::string q_fam[] = { "qf03", "qd04" }, s_fam[] = { "sfh",  "sd" };

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  // 1: DIAMOND, 3: Oleg I, 4: Oleg II.
  FieldMap_filetype = 1; sympl = false;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

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
    chk_high_oord_achr();
    exit(0);
  }

  if (false) {
    chk_b3();
    exit(0);
  }

  if (false) {
    chk_optics(0.0, 0.0, 13.04171, 8.795924, 1.397388e-03, 0.0, 0.0, 0.0);
    prt_lat("linlat1.out", globval.bpm, true);
    prt_lat("linlat.out", globval.bpm, true, 10);
    exit(0);
  }

  if (false) {
    globval.Cavity_on = true;
    track(6e-3, 0.1e-3);
    exit(0);
  }

  if (true) no_sxt();

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

  if (true) GetEmittance(ElemIndex("cav"), true);

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
