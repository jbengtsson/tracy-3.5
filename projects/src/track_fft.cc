#define NO 1

#include "tracy_lib.h"

int no_tps = NO;

const char home_dir[] = "/home/simon";




// copied here from old 2011 nsls-ii_lib.cc because no longer included in 2017
// version
void LoadAlignTol(const char *AlignFile, const bool Scale_it,
		  const double Scale, const bool new_rnd, const int seed)
{
  char      line[max_str], Name[max_str],  type[max_str];
  int       j, k, Fnum, seed_val;
  long int  loc;
  double    dx, dy, dr;  // x and y misalignments [m] and roll error [rad]
  double    dr_deg;
  bool      rms = false, set_rnd;
  FILE      *fp;

  const bool  prt = true;

  if (prt) {
    printf("\n");
    printf("reading in %s\n", AlignFile);
  }

  fp = file_read(AlignFile);

  printf("\n");
  if (new_rnd)
    printf("set alignment errors\n");
  else
    printf("scale alignment errors: %4.2f\n", Scale);

  set_rnd = false;
  while (fgets(line, max_str, fp) != NULL) {
    if (prt) printf("%s", line);

    if ((strstr(line, "#") == NULL) && (strcmp(line, "\r\n") != 0)) {
      sscanf(line, "%s", Name);
      //check for whether to set seed
      if (strcmp("seed", Name) == 0) {
	set_rnd = true;
	sscanf(line, "%*s %d", &seed_val);
	seed_val += 2*seed;
	printf("setting random seed to %d\n", seed_val);
	iniranf(seed_val); 
      } else {
	sscanf(line,"%*s %s %lf %lf %lf", type, &dx, &dy, &dr);
	dr_deg = dr*180.0/M_PI;

	if (strcmp(type, "rms") == 0){
	  rms = true;
	  printf("<rms> ");
	}
	else if (strcmp(type, "sys") == 0){
	  rms = false;
	  printf("<sys> ");
	}
	else {
	  printf("LoadAlignTol: element %s:  need to specify rms or sys\n",
		 Name);
	  exit_(1);
	}

	if (rms && !set_rnd) {
	  printf("LoadAlignTol: seed not defined\n");
	  exit_(1);
	}

	if (Scale_it) {
	  dx *= Scale; dy *= Scale; dr *= Scale;
	} 

	if (strcmp("all", Name) == 0) {
	  printf("misaligning all:         dx = %e, dy = %e, dr = %e\n",
		 dx, dy, dr);
	  if(rms)
	    misalign_rms_type(All, dx, dy, dr_deg, new_rnd);
	  else
	    misalign_sys_type(All, dx, dy, dr_deg);
	} else if (strcmp("girder", Name) == 0) {
	  printf("misaligning girders:     dx = %e, dy = %e, dr = %e\n",
		 dx, dy, dr);
	  if (rms)
	    misalign_rms_girders(globval.gs, globval.ge, dx, dy, dr_deg,
				 new_rnd);
	  else
	    misalign_sys_girders(globval.gs, globval.ge, dx, dy, dr_deg);
	} else if (strcmp("dipole", Name) == 0) {
	  printf("misaligning dipoles:     dx = %e, dy = %e, dr = %e\n",
		 dx, dy, dr);
	  if (rms)
	    misalign_rms_type(Dip, dx, dy, dr_deg, new_rnd);
	  else
	    misalign_sys_type(Dip, dx, dy, dr_deg);
	} else if (strcmp("quad", Name) == 0) {
	  printf("misaligning quadrupoles: dx = %e, dy = %e, dr = %e\n",
		 dx, dy, dr);
	  if (rms)
	    misalign_rms_type(Quad, dx, dy, dr_deg, new_rnd);
	  else
	    misalign_sys_type(Quad, dx, dy, dr_deg);
	} else if (strcmp("sext", Name) == 0) {
	  printf("misaligning sextupoles:  dx = %e, dy = %e, dr = %e\n",
		 dx, dy, dr);
	  if (rms)
	    misalign_rms_type(Sext, dx, dy, dr_deg, new_rnd);
	  else
	    misalign_sys_type(Sext, dx, dy, dr_deg);
	} else if (strcmp("bpm", Name) == 0) {
	  printf("misaligning bpms:        dx = %e, dy = %e, dr = %e\n",
		 dx, dy, dr);
	  for (k = 0; k < 2; k++)
	    for (j = 1; j <= n_bpm_[k]; j++) {
	      loc = bpms_[k][j];
	      if (rms)
		misalign_rms_elem(Lattice.Cell[loc].Fnum, Lattice.Cell[loc].Knum,
				  dx, dy, dr_deg, new_rnd);
	      else
		misalign_sys_elem(Lattice.Cell[loc].Fnum, Lattice.Cell[loc].Knum,
				  dx, dy, dr_deg);
	    }
	} else {
	  Fnum = Lattice.Elem_Index(Name);
	  if(Fnum > 0) {
	    printf("misaligning all %s:  dx = %e, dy = %e, dr = %e\n",
		   Name, dx, dy, dr);
	    if (rms)
	      misalign_rms_fam(Fnum, dx, dy, dr_deg, new_rnd);
	    else
	      misalign_sys_fam(Fnum, dx, dy, dr_deg);
	  } else 
	    printf("LoadAlignTol: undefined element %s\n", Name);
	}
      }
    } else
      printf("%s", line);
  }

  fclose(fp);
}





//copied here form nsls-ii_lib.cc; needed for LoadFieldErr_scl
char* get_prm_scl(void)
{
  char  *prm;

  prm = strtok(NULL, " \t");
  if ((prm == NULL) || (strcmp(prm, "\r\n") == 0)) {
    printf("get_prm: incorrect format\n");
    exit_(1);
  }

  return prm;
}





//copied here form nsls-ii_lib.cc to add incrementing seed value
void LoadFieldErr_scl(const char *FieldErrorFile, const bool Scale_it,
		      const double Scale, const bool new_rnd, const int m) 
{  
  bool    rms, set_rnd;
  char    line[max_str], name[max_str], type[max_str], *prm;
  int     k, n, seed_val;
  double  Bn, An, r0;
  FILE    *inf;

  const bool  prt = true;

  inf = file_read(FieldErrorFile);

  set_rnd = false; 
  printf("\n");
  while (fgets(line, max_str, inf) != NULL) {
    if (strstr(line, "#") == NULL) {
      // check for whether to set new seed
      sscanf(line, "%s", name); 
      if (strcmp("seed", name) == 0) {
	set_rnd = true;
	sscanf(line, "%*s %d", &seed_val); 
	seed_val += 2*m;
	printf("setting random seed to %d\n", seed_val);
	iniranf(seed_val); 
      } else {
	sscanf(line, "%*s %s %lf", type, &r0);
	printf("%-4s %3s %7.1le", name, type, r0);
	rms = (strcmp("rms", type) == 0)? true : false;
	if (rms && !set_rnd) {
	  printf("LoadFieldErr: seed not defined\n");
	  exit_(1);
	}
	// skip first three parameters
	strtok(line, " \t");
	for (k = 1; k <= 2; k++)
	  strtok(NULL, " \t");
	while (((prm = strtok(NULL, " \t")) != NULL) &&
	       (strcmp(prm, "\r\n") != 0)) {
	  sscanf(prm, "%d", &n);
	  prm = get_prm_scl();
	  sscanf(prm, "%lf", &Bn);
	  prm = get_prm_scl(); 
	  sscanf(prm, "%lf", &An);
	  if (Scale_it)
	    {Bn *= Scale; An *= Scale;}
	  if (prt)
	    printf(" %2d %9.1e %9.1e\n", n, Bn, An);
	  // convert to normalized multipole components
	  SetFieldErrors(name, rms, r0, n, Bn, An, true);
	}
      }
    } else
      printf("%s", line);
  }

  fclose(inf);
}




//adapted from orb_corr() in nsls-ii_lib.cc (was broken after generalization to N families of BPMs/CORRs)
//
//for n_orbit>=1 get error orbit, correct orbit n_orbit times, print results
//for n_orbit==0 get only error orbit w/o correction -> use this to get COD before OCO
//(i.e. see effect of errors without any correction -> amplification)
//
bool orb_corr_scl(const int n_orbit)
{
  bool      cod = false;
  int       i;
  long      lastpos;
  Vector2   xmean, xsigma, xmax;

  double scl = 1e0;

  printf("\n");

//  FitTune(qf, qd, nu_x, nu_y);
//  printf("\n");
//  printf("  qf = %8.5f qd = %8.5f\n",
//	   GetKpar(qf, 1, Quad), GetKpar(qd, 1, Quad));

  int n_orbit2;
  if (n_orbit == 0) {
    n_orbit2 = 1;
  } else
    n_orbit2 = n_orbit;
  
  globval.CODvect.zero();
  for (i = 1; i <= n_orbit2; i++) {
    cod = Lattice.getcod(0.0, lastpos);
    if (cod) {
      codstat(xmean, xsigma, xmax, globval.Cell_nLoc, false); //false = take values only at BPM positions
      printf("\n");
      printf("RMS orbit [mm]: %8.1e +/- %7.1e, %8.1e +/- %7.1e\n", 
	     1e3*xmean[X_], 1e3*xsigma[X_], 1e3*xmean[Y_], 1e3*xsigma[Y_]);
      if (n_orbit != 0) {
	// J.B. 08/24/17 ->
	// The call to:
	//   gcmat(Lattice.Elem_Index("bpm_m"), Lattice.Elem_Index("corr_h"), 1); gcmat(Lattice.Elem_Index("bpm_m"), Lattice.Elem_Index("corr_v"), 2);
   	// configures the bpm system.
	// lsoc(1, Lattice.Elem_Index("bpm_m"), Lattice.Elem_Index("corr_h"), 1);  //updated from older T3 version
	// lsoc(1, Lattice.Elem_Index("bpm_m"), Lattice.Elem_Index("corr_v"), 2);  //updated from older T3 version
	lsoc(1, scl); lsoc(2, scl);
	// -> J.B. 08/24/17:
	cod = Lattice.getcod(0.0, lastpos);
	if (cod) {
	  codstat(xmean, xsigma, xmax, globval.Cell_nLoc, false); //false = take values only at BPM positions
	  printf("RMS orbit [mm]: %8.1e +/- %7.1e, %8.1e +/- %7.1e\n", 
		 1e3*xmean[X_], 1e3*xsigma[X_], 1e3*xmean[Y_], 1e3*xsigma[Y_]);
	} else
	  printf("orb_corr: failed\n");
      }
    } else
      printf("orb_corr: failed\n");
  }
  
  Lattice.prt_cod("orb_corr.out", Lattice.Elem_Index("bpm_m"), true);  //updated from older T3 version

  return cod;
}





void track_fft(const int n_turn,
	       const double x, const double y, const double delta)
{
  const int  n_peaks = 5;
  // some peculiar problem using n_turn to declare array sizes
  // -> phase out FFT.cc with NumRec
  const int  n_max = 2048;

  char            str[max_str];
  long int        n, lastn, lastpos;
  int             i, n_x, n_y;
  double          twoJx[n_max], phix[n_max], twoJy[n_max], phiy[n_max];
  double          nu[2][n_peaks], A[2][n_peaks], f, nu0[2][2], A0[2];
  ss_vect<double> ps;
  FILE            *fp;

  const int   max_no_peak = 12, n_peak = 8;

  const char  file_name[] = "track_fft.out";

  ps[x_] = x; ps[px_] = 0.0; ps[y_] = y; ps[py_] = 0.0; Lattice.getfloqs(ps);
  twoJx[0] = sqr(ps[x_]) + sqr(ps[px_]);
  twoJy[0] = sqr(ps[y_]) + sqr(ps[py_]);
  phix[0] = atan2(ps[px_], ps[x_]); phiy[0] = atan2(ps[py_], ps[y_]);

  Lattice.track(file_name, twoJx[0], phix[0], twoJy[0], phiy[0], delta,
		n_turn, lastn, lastpos,
		2,
		Lattice.Cell[Lattice.Elem_GetPos(Lattice.Elem_Index("cav"), 1)]
		.Elem.C->Pfreq);

  if (lastn == n_turn) {
    GetTrack(file_name, &n, twoJx, phix, twoJy, phiy);
    printf("\n");
    printf("Read %ld turns\n", n);

    sin_FFT(n, twoJx); sin_FFT(n, twoJy); sin_FFT(n, phix); sin_FFT(n, phiy);

    GetPeaks(n, twoJx, 1, nu[X_], A[X_]); GetPeaks(n, twoJy, 1, nu[Y_], A[Y_]);
    A0[X_] = A[X_][0]; A0[Y_] = A[Y_][0];
    GetPeaks(n, phix, 1, nu[X_], A[X_]); GetPeaks(n, phiy, 1, nu[Y_], A[Y_]);
    nu0[X_][1] = nu[X_][0]; nu0[Y_][1] = nu[Y_][0];

    GetTrack(file_name, &n, twoJx, phix, twoJy, phiy);
    rm_mean(n, twoJx); rm_mean(n, twoJy);
    sin_FFT(n, twoJx); sin_FFT(n, twoJy); sin_FFT(n, phix); sin_FFT(n, phiy);

    GetPeaks(n, twoJx, n_peaks, nu[X_], A[X_]);
    GetPeaks(n, twoJy, n_peaks, nu[Y_], A[Y_]);

//  strtok(file_name, "c");
//  printf("%s\n", name);
//  exit(0);

    strcpy(str, file_name); strcat(str, ".fft");
    fp = file_write(str);
    for (i = 0; i < n/2+1; i++)
      fprintf(fp, "%7.5f %12.5e %8.5f %12.5e %8.5f\n",
	      (double)i/(double)n, twoJx[i], phix[i], twoJy[i], phiy[i]);
    fclose(fp);

    fp = file_write("FFT.out");

    fprintf(fp, "nu_x = %7.5f, 1-nu_x = %7.5f, nu_y = %7.5f, 1-nu_y = %7.5f\n",
	    nu0[X_][1], 1.0-nu0[X_][1], nu0[Y_][1], 1.0-nu0[Y_][1]);

    fprintf(fp, "\n");
    fprintf(fp, "Horizontal plane:\n");
    fprintf(fp, "\n");
    fprintf(fp, "%1d %8.6f %9.3e\n", 0, nu0[X_][1], A0[X_]);
    for (i = 0; i < n_peak; i++) {
      FindRes(max_no_peak, nu0[X_][1], nu0[Y_][1], nu[X_][i], &n_x, &n_y);
      f = fabs(n_x*nu0[X_][1]+n_y*nu0[Y_][1]); f -= (int)f;
      if (f > 0.5) f = 1.0 - f;
      fprintf(fp, "%1d %8.5f %9.3e %2d %2d %8.5f\n",
	     i+1, nu[X_][i], A[X_][i], n_x, n_y, fabs(nu[X_][i])-f);
    }
    
    fprintf(fp, "\n");
    fprintf(fp, "Vertical plane:\n");
    fprintf(fp, "\n");
    fprintf(fp, "%1d %8.5f %9.3e\n", 0, nu0[Y_][1], A0[Y_]);
    for (i = 0; i < n_peak; i++) {
      FindRes(max_no_peak, nu0[X_][1], nu0[Y_][1], nu[Y_][i], &n_x, &n_y);
      f = fabs(n_x*nu0[X_][1]+n_y*nu0[Y_][1]); f -= (int)f;
      if (f > 0.5) f = 1.0 - f;
      fprintf(fp, "%1d %8.5f %9.3e %2d %2d %8.5f\n",
	      i+1, nu[Y_][i], A[Y_][i], n_x, n_y, fabs(nu[Y_][i])-f);
    }

    fclose(fp);

  } else
    printf("particle lost on turn no %ld", lastn);
}




/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////




int main(int argc, char *argv[])
{
  
  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false;
  globval.pathlength = false; globval.bpm         = 0;

  Lattice.Read_Lattice(argv[1]);

  // to turn off sextupoles:
  //no_sxt();

  Lattice.Ring_GetTwiss(true, 0.0); printglob();

  int n_turn;
  double x0, y0, px0, py0, delta0, f_rf;

  // M4 parameters
  if(true) {
    x0        = -4.665e-3;
    px0       = 1.067e-3 - 1.167e-3;
    y0        = 0.00e-3;
    py0       = 0.0e-3;
  } else {
    // M5 parameters
    x0        = 0.001e-3;
    px0       = 0.0;
    y0        = 0.001e-3;
    py0       = 0.0e-3;
  }

  delta0    = -0.5e-2;
  n_turn    = 1000;
  f_rf      =
    Lattice.Cell[Lattice.Elem_GetPos(Lattice.Elem_Index("cav"), 1)]
    .Elem.C->Pfreq;
  
  
  
  if (false) {
    const long  seed = 1121;
    iniranf(seed); setrancut(2.0);
    Lattice.Ring_GetTwiss(true, 0e-2); printglob(); //gettwiss computes one-turn matrix arg=(w or w/o chromat, dp/p)
    globval.gs = Lattice.Elem_Index("GS");
    globval.ge = Lattice.Elem_Index("GE");
    // compute response matrix (needed for OCO)
    gcmat(Lattice.Elem_Index("bpm_m"), Lattice.Elem_Index("corr_h"), 1);
    gcmat(Lattice.Elem_Index("bpm_m"), Lattice.Elem_Index("corr_v"), 2);
    // reset orbit trims
    zero_trims();
    LoadFieldErr_scl("/home/simon/projects/in/lattice/FieldErr.5e-4.dat", false,
		     1.0, true, 20); //last number is seed no.
    LoadAlignTol("/home/simon/projects/in/lattice/AlignErr.required+.dat",
		 false, 1.0, true, 20); //last number is seed no.
    bool       cod;
    cod = orb_corr_scl(3);  // use orb_corr_scl(0) to show orbit deviations BEFORE correction -> ampl. factor
    // use orb_corr_scl(3) to correct orbit in 3 iterations
    Lattice.GetEmittance(Lattice.Elem_Index("cav"), true);
  }
  


  // FOR TxT TRACKING AND FFT (results given in action-angle coordinates)
  if (false) {
    track_fft(n_turn, x0, y0, delta0);  //routine defined above
  }


  // FOR TxT TRACKING AND PLOTTING VS TURN NO.
  if (!true) {
    long int lastn, lastpos;
    globval.Cavity_on  = true; globval.radiation   = true;
    Lattice.track("track_fft.out", x0, px0, y0, py0, delta0, n_turn,
		  lastn, lastpos, 0, f_rf);  //track is in physlib.cc
    //                                                                       ^ floqs
    // floqs: 0 = phase space            (x, px, y, py, delta, ct)
    //        1 = floquet space          (x^, px^, y^, py^, delta, ct)
    //        2 = action-angle variables (2Jx, phix, 2Jy, phiy, delta, ct)
  }
  
  
  // FOR ELEMENTxELEMENT TRACKING
  if (!false) {

    long int  j, lastn, lastpos;
    struct    tm *newtime;
    FILE      *fp;
    
    globval.Cavity_on  = true; globval.radiation   = true;
    Lattice.track("track_fft.out", x0, px0, y0, py0, delta0, 1,
		  lastn, lastpos, 0, f_rf); // track_fft.out makes no sense here
    
    fp = file_write("track_ExE.out");
    
    newtime = GetTime();
    
    fprintf(fp,"# TRACY III v.3.5 -- %s -- %s \n",
	    "track_ExE.out", asctime2(newtime));
    fprintf(fp, "#         s   name              code  x          xp         "
	    "y          yp         delta      ctau\n");
    fprintf(fp,"#        [m]                          [mm]       [mrad]     "
	    "[mm]       [mrad]     [%%]        [?]\n");
    fprintf(fp, "#\n");
    
    for (j = 0; j <= globval.Cell_nLoc; j++){
      fprintf(fp, "%4li %8.3f %s %6.2f %10.3e %10.3e %10.3e %10.3e %10.3e "
	      "%10.3e\n",
	      j, Lattice.Cell[j].S, Lattice.Cell[j].Elem.PName,
	      get_code(Lattice.Cell[j]),
	      Lattice.Cell[j].BeamPos[0]*1e3, Lattice.Cell[j].BeamPos[1]*1e3,
	      Lattice.Cell[j].BeamPos[2]*1e3, Lattice.Cell[j].BeamPos[3]*1e3,
	      Lattice.Cell[j].BeamPos[4]*100, Lattice.Cell[j].BeamPos[5]);
	}
    
    fclose(fp);
    
  }
  
}
