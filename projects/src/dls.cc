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
      for (j = 0; j <= Lattice.param.Cell_nLoc; j++)
	if (all || ((Lattice.Cell[j].Kind == Mpole) &&
		    (Lattice.Cell[j].Elem.M->n_design == Sext))) {
	  n++;
	  for (k = 0; k < 6; k++) {
	    x1[k][n-1] += Lattice.Cell[j].BeamPos[k];
	    x2[k][n-1] += sqr(Lattice.Cell[j].BeamPos[k]);
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
    if (all || ((Lattice.Cell[j].Kind == Mpole) &&
		(Lattice.Cell[j].Elem.M->n_design == Sext))) {
      n++;
      for (k = 0; k < 6; k++) {
	x_mean[k].push_back(x1[k][n-1]/n_cod);
	x_sigma[k].push_back(sqrt((n_cod*x2[k][n-1]-sqr(x1[k][n-1]))
				  /(n_cod*(n_cod-1.0))));
      }
      fprintf(fp, "%8.3f %6.2f %10.3e +/- %10.3e %10.3e +/- %10.3e\n",
	      Lattice.Cell[j].S, get_code(Lattice.Cell[j]),
	      1e3*x_mean[x_][n-1], 1e3*x_sigma[x_][n-1],
	      1e3*x_mean[y_][n-1], 1e3*x_sigma[y_][n-1]);
    } else
      fprintf(fp, "%8.3f %6.2f\n",
	      Lattice.Cell[j].S, get_code(Lattice.Cell[j]));
  
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
	  Lattice.Cell[0].BeamPos[0], Lattice.Cell[0].BeamPos[1],
	  Lattice.Cell[0].BeamPos[2], Lattice.Cell[0].BeamPos[3],
	  Lattice.Cell[0].BeamPos[4], Lattice.Cell[0].BeamPos[5] );
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
	     Lattice.Cell[loc].S, Lattice.Cell[loc].Beta[X_],
	     Lattice.Cell[loc].Beta[Y_],
	     Lattice.Cell[loc].Eta[X_]);
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
	   Lattice.Cell[loc].S, Lattice.Cell[loc].Beta[X_],
	   Lattice.Cell[loc].Beta[Y_],
	   Lattice.GetnKid(Fam[j]));
  }
}


int main(int argc, char *argv[])
{
  long int         lastn, lastpos;
  int              b2_fam[2], b3_fam[2];
  double           b2[2], a2, b3[2], b3L[2], a3, a3L, f_rf;
  double           alpha[2], beta[2], eta[2], etap[2];
  std::vector<int> Fam;
  ostringstream    str;

  const long        seed   = 1121;
  const int         n_turn = 2064;
  const double      delta  = 1.5e-2,
  //                   nu[]    = { 102.18/20.0, 68.30/20.0 };
  // const std::string q_fam[] = { "qfe", "qde" }, s_fam[] = { "sfh",  "sd" };
  //                   nu[]    = { 39.1/12.0, 15.25/12.0 };
  //                   // nu[]    = { 3.266+0.01, 1.275 };
  // const std::string q_fam[] = { "qm2b", "qm3" }, s_fam[] = { "sfh",  "sd" };
                    nu[]    = { 51.21/6.0, 17.35/6.0 };
  const std::string q_fam[] = { "q1b", "q1d" }, s_fam[] = { "sfh",  "sd" };

  Lattice.param.H_exact    = false; Lattice.param.quad_fringe = false;
  Lattice.param.Cavity_on  = false; Lattice.param.radiation   = false;
  Lattice.param.emittance  = false; Lattice.param.IBS         = false;
  Lattice.param.pathlength = false; Lattice.param.bpm         = 0;

  if (true) {
    Lattice.Read_Lattice(argv[1]);
  } else
    Lattice.rdmfile(argv[1]);

  if (false) {
    Lattice.param.Cavity_on = true;
    track(6e-3, 0.1e-3);
    exit(0);
  }

  if (false) no_sxt();

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

    gcmat(Lattice.param.bpm, Lattice.param.hcorr, 1);
    gcmat(Lattice.param.bpm, Lattice.param.vcorr, 2);

    get_cod_rms(50e-6, 50e-6, 100, true);

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

    f_rf =
      Lattice.Cell[Lattice.Elem_GetPos(Lattice.Elem_Index("cav"), 1)]
      .Elem.C->freq;
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
    Lattice.track("track.out", 5.5e-3, 0e0, 1e-6, 0e0, 0e0, n_turn,
		  lastn, lastpos, 0, 0*f_rf);
  }

  if (true) {
    Lattice.param.Cavity_on = true;
    Lattice.get_dynap(delta, 25, n_turn, false);
  }

}
