#define NO 1

#include "tracy_lib.h"

int no_tps = NO; // arbitrary TPSA order is defined locally

extern bool  freq_map;

const double  delta  = 4.5e-2; // delta for off-momentum DA



/////////////////////////////////////////////////////////////////////////////////////////////////



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




void get_cod_rms_data(const int n_seed, const int nfam, const int fnums[], const double dx, const double dy, const double dr,
		      double x_mean[][6], double x_sigma[][6], double theta_mean[][2], double theta_sigma[][2])
{
  const int  n_elem = globval.Cell_nLoc+1;

  bool       cod;
  long int   j;
  int        i, k, ncorr[2];
  double     x1[n_elem][6], x2[n_elem][6], theta1[n_elem][6], theta2[n_elem][6], a1L, b1L;;

  for (j = 0; j <= globval.Cell_nLoc; j++){
    for (k = 0; k < 6; k++) {
      x1[j][k] = 0.0; x2[j][k] = 0.0;
    }
    for (k = 0; k < 2; k++){
      theta1[j][k] = 0;
      theta2[j][k] = 0;
    }
  }  
   
  for (i = 0; i < n_seed; i++) {
    // reset orbit trims
    zero_trims();

    for (j = 0; j < nfam; j++)
      misalign_rms_fam(fnums[j], dx, dy, dr, true);
    
    cod = orb_corr(3);
    
    if (cod) {
      
      for (k = 0; k < 2; k++)
	ncorr[k] = 0;
      
      for (j = 0; j <= globval.Cell_nLoc; j++){ // read back beam pos at bpm

	if ( (Cell[j].Elem.Pkind == Mpole) || (Cell[j].Elem.Pkind == drift) ) {
	  for (k = 0; k < 6; k++) {
	    x1[j][k] += Cell[j].BeamPos[k];
	    x2[j][k] += sqr(Cell[j].BeamPos[k]);
	  }
	  
	  if ( (Cell[j].Fnum == ElemIndex("corr_h")) || (Cell[j].Fnum == ElemIndex("corr_v")) ){ // read back corr strength at corr
	    get_bnL_design_elem(Cell[j].Fnum, Cell[j].Knum, Dip, b1L, a1L);
	    if ( Cell[j].Fnum == ElemIndex("corr_h") ){
	      ncorr[X_]++;
	      theta1[ncorr[X_]-1][X_] += b1L;
	      theta2[ncorr[X_]-1][X_] += sqr(b1L);
	    } else if ( Cell[j].Fnum == ElemIndex("corr_v") ){
	      ncorr[Y_]++;
	      theta1[ncorr[Y_]-1][Y_] += a1L;
	      theta2[ncorr[Y_]-1][Y_] += sqr(a1L);
	    }
	  }
	}
      }
      
    } else
      printf("orb_corr: failed\n");    
  }
  
  for (k = 0; k < 2; k++)
    ncorr[k] = 0;
  
  for (j = 0; j <= globval.Cell_nLoc; j++){
        
    if ( (Cell[j].Elem.Pkind == Mpole) || (Cell[j].Elem.Pkind == drift) ) {
      for (k = 0; k < 6; k++) {
	x_mean[j][k] = x1[j][k]/n_seed;
	x_sigma[j][k] = sqrt((n_seed*x2[j][k]-sqr(x1[j][k]))
			     /(n_seed*(n_seed-1.0)));
      }
      
      if ( Cell[j].Fnum == ElemIndex("corr_h") ){
	ncorr[X_]++;
	theta_mean[ncorr[X_]-1][X_] = theta1[ncorr[X_]-1][X_]/n_seed;
	theta_sigma[ncorr[X_]-1][X_] = sqrt((n_seed*theta2[ncorr[X_]-1][X_]-sqr(theta1[ncorr[X_]-1][X_]))
					    /(n_seed*(n_seed-1.0)));
	
      } else if ( Cell[j].Fnum == ElemIndex("corr_v") ){
	ncorr[Y_]++;
	theta_mean[ncorr[Y_]-1][Y_] = theta1[ncorr[Y_]-1][Y_]/n_seed;
	theta_sigma[ncorr[Y_]-1][Y_] = sqrt((n_seed*theta2[ncorr[Y_]-1][Y_]-sqr(theta1[ncorr[Y_]-1][Y_]))
					    /(n_seed*(n_seed-1.0)));
      }
    }
  }
}



void prt_cod_rms_data(const char name[], double x_mean[][6], double x_sigma[][6], double theta_mean[][2], double theta_sigma[][2])
{
  long     j;
  int      k, ncorr[2];
  FILE     *fp;
  struct    tm *newtime;

  fp = file_write(name);


  /* Get time and date */
  newtime = GetTime();

  fprintf(fp,"# TRACY III v.3.5 -- %s -- %s \n",
	  name, asctime2(newtime));
  fprintf(fp,"#         s   name              code  xcod_mean +/-  xcod_rms   ycod_mean +/-  ycod_rms"
	     "   dipx_mean +/-  dipx_sigma dipy_mean +/-  dipy_sigma\n");
  fprintf(fp,"#        [m]                             [mm]           [mm]       [mm]           [mm]"
	     "      [mrad]         [mrad]     [mrad]         [mrad]\n");
  fprintf(fp, "#\n");


  for (k = 0; k < 2; k++)
    ncorr[k] = 0;
  
  for (j = 0; j <= globval.Cell_nLoc; j++){
   
    fprintf(fp, "%4li %8.3f %s %6.2f %10.3e +/- %10.3e %10.3e +/- %10.3e",
	    j, Cell[j].S, Cell[j].Elem.PName, get_code(Cell[j]),
	    1e3*x_mean[j][x_], 1e3*x_sigma[j][x_],
	    1e3*x_mean[j][y_], 1e3*x_sigma[j][y_]);
   
    if ( Cell[j].Fnum == ElemIndex("corr_h") ){
      ncorr[X_]++;
      fprintf(fp, " %10.3e +/- %10.3e %10.3e +/- %10.3e\n",
	      1e3*theta_mean[ncorr[X_]-1][X_], 1e3*theta_sigma[ncorr[X_]-1][X_], 0.0, 0.0);

    } else if ( Cell[j].Fnum == ElemIndex("corr_v") ){
      ncorr[Y_]++;
      fprintf(fp, " %10.3e +/- %10.3e %10.3e +/- %10.3e\n",
	      0.0, 0.0, 1e3*theta_mean[ncorr[Y_]-1][Y_], 1e3*theta_sigma[ncorr[Y_]-1][Y_]);
    } else 
      fprintf(fp, " %10.3e +/- %10.3e %10.3e +/- %10.3e\n",
	      0.0, 0.0, 0.0, 0.0);
  }

  fclose(fp);
}


void add_family( const char *name, int &nfam, int fnums[] )
{
  int fnum;
  fnum = ElemIndex(name); // if the family doesn't exist code bails here
  if (fnum != 0)
    fnums[nfam++] = fnum;
  else {                  // hence this is useless
    printf("Family not defined. %s %d\n", name, fnum);  
    exit(1);
  }
}


//copied here form nsls-ii_lib.cc to make changes
void get_cod_rms_scl(const double dx, const double dy, const double dr,
		     const int n_seed)
{
  const int  n_elem = globval.Cell_nLoc+1;

  int      fnums[25], nfam; 
  double   x_mean[n_elem][6], x_sigma[n_elem][6], theta_mean[n_elem][2], theta_sigma[n_elem][2];

  nfam = 0;
  add_family("qf", nfam, fnums);
  add_family("qfm", nfam, fnums);
  add_family("qfend", nfam, fnums);
  add_family("qdend", nfam, fnums);
  add_family("sfi", nfam, fnums);
  add_family("sfo", nfam, fnums);
  add_family("sfm", nfam, fnums);
  add_family("sd", nfam, fnums);
  add_family("sdend", nfam, fnums);

  get_cod_rms_data(n_seed, nfam, fnums, dx, dy, dr, x_mean, x_sigma, theta_mean, theta_sigma);
  prt_cod_rms_data("cod_rms.out", x_mean, x_sigma, theta_mean, theta_sigma);

}



/////////////////////////////////////////////////////////////////////////////////////////////////




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
    cod = getcod(0.0, lastpos);
    if (cod) {
      codstat(xmean, xsigma, xmax, globval.Cell_nLoc, false); //false = take values only at BPM positions
      printf("\n");
      printf("RMS orbit [mm]: %8.1e +/- %7.1e, %8.1e +/- %7.1e\n", 
	     1e3*xmean[X_], 1e3*xsigma[X_], 1e3*xmean[Y_], 1e3*xsigma[Y_]);
      if (n_orbit != 0) {
	// lsoc(1, ElemIndex("bpm_m"), ElemIndex("corr_h"), 1);  //updated from older T3 version
	// lsoc(1, ElemIndex("bpm_m"), ElemIndex("corr_v"), 2);  //updated from older T3 version
	cod = getcod(0.0, lastpos);
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
  
  prt_cod("orb_corr.out", ElemIndex("bpm_m"), true);  //updated from older T3 version

  return cod;
}






/////////////////////////////////////////////////////////////////////////////////////////////////




//this is a mashup of get_cod_rms_scl
//instead of assigning errors to families we load and apply AlignErr.dat and FieldErr.dat
//a correction is then applied and the results saved
//this is repeated for n_seed seeds and overall statistics saved in cod_rms.out
void get_cod_rms_scl_new(const int n_seed)
{
  
  const int  n_elem = globval.Cell_nLoc+1;
  
  double     x_mean[n_elem][6], x_sigma[n_elem][6], theta_mean[n_elem][2], theta_sigma[n_elem][2];
  bool       cod;
  long int   j;
  int        i, k, ncorr[2], l;
  double     x1[n_elem][6], x2[n_elem][6], theta1[n_elem][6], theta2[n_elem][6], a1L, b1L;


  /////////////////////////////////////////////////////////////////////
  // this is required for output to beam_envelope_stats.out
  FILE  *fp;
  fp = file_write("beam_envelope_stats.out");
  double  sumyoverx[n_elem], sumtwist[n_elem];
  double  sumyoverxsq[n_elem], sumtwistsq[n_elem];
  double  maxyoverx[n_elem], maxtwist[n_elem];
  double  avgyoverx[n_elem], avgtwist[n_elem];
  double  sigyoverx[n_elem], sigtwist[n_elem];
  /////////////////////////////////////////////////////////////////////
  

  for (j = 0; j <= globval.Cell_nLoc; j++){
    for (k = 0; k < 6; k++) {
      x1[j][k] = 0.0; x2[j][k] = 0.0;
    }
    for (k = 0; k < 2; k++){
      theta1[j][k] = 0;
      theta2[j][k] = 0;
    }
  }  
  
  for (i = 0; i < n_seed; i++) {
    
    cout << endl << "*** COD ITERATION " << i+1 << " of " << n_seed << " ***" << endl;
    
    // reset orbit trims
    zero_trims();
    
    LoadFieldErr_scl("/home/simon/projects/in/lattice/FieldErr.m4.20140514.dat", false, 1.0, true, i);
    LoadAlignTol("/home/simon/projects/in/lattice/AlignErr.m4.20140130.dat", false, 1.0, true, i);
    
    // if COD fails because of large misalignments, introduce intermediate steps to approach solution
    if (false) {
      int  steps = 3;
      for (l = 1; l <=steps; l++) {
	cout << endl << "*** COD STEPPING " << l << " of " << steps << " (in iteration " << i+1 << "/" << n_seed << " ) ***" << endl;
	LoadAlignTol("/home/simon/projects/in/lattice/AlignErr.m4.20140130.dat", true, (double)l/double(steps), false, i);
	cod = orb_corr_scl(3);
      }
    } else
      cod = orb_corr_scl(3);  // use orb_corr_scl(0) to show orbit deviations BEFORE correction -> ampl. factor
                              // use orb_corr_scl(3) to correct orbit in 3 iterations
    // get coupling, ver. disp., etc. BEFORE/after correction
    GetEmittance(ElemIndex("cav"), true);
    prt_beamsizes(); // writes beam_envelope.out; contains sigmamatrix(s), theta(s)


    /////////////////////////////////////////////////////////////////////
    // this should get theta(s) and kappa(s) statistics over many seeds
    // (derived from prt_beamsizes)
    if (true) {
      // acquire data for each seed -> sum of seed values, sum of squares of seed values, max values
      for (k = 0; k <= globval.Cell_nLoc; k++) {
	sumyoverx[k] += sqrt(Cell[k].sigma[y_][y_])/sqrt(Cell[k].sigma[x_][x_]); // yoverx = sigma_y/sigma_x = kappa(s)
	sumtwist[k] += atan2(2e0*Cell[k].sigma[x_][y_],
			 Cell[k].sigma[x_][x_]-Cell[k].sigma[y_][y_])/2e0*180.0/M_PI;
	sumyoverxsq[k] += Cell[k].sigma[y_][y_]/Cell[k].sigma[x_][x_];
	sumtwistsq[k] += sqr(atan2(2e0*Cell[k].sigma[x_][y_],
			 Cell[k].sigma[x_][x_]-Cell[k].sigma[y_][y_])/2e0*180.0/M_PI);
	if ( (sqrt(Cell[k].sigma[y_][y_])/sqrt(Cell[k].sigma[x_][x_])) > maxyoverx[k] )
	  maxyoverx[k] = sqrt(Cell[k].sigma[y_][y_])/sqrt(Cell[k].sigma[x_][x_]);
	if ( abs(atan2(2e0*Cell[k].sigma[x_][y_],
			Cell[k].sigma[x_][x_]-Cell[k].sigma[y_][y_])/2e0*180.0/M_PI) > abs(maxtwist[k]) )
	  maxtwist[k] = abs(atan2(2e0*Cell[k].sigma[x_][y_],
			      Cell[k].sigma[x_][x_]-Cell[k].sigma[y_][y_])/2e0*180.0/M_PI);
      }    
      // then (after final seed) calculate statistics
      if ( i == (n_seed-1) ) {
	// get averages and rms from sums and sums of squares, respectively
	for (k = 0; k <= globval.Cell_nLoc; k++) {
	  avgyoverx[k] = sumyoverx[k]/n_seed;
	  avgtwist[k] = sumtwist[k]/n_seed;
	  sigyoverx[k] = sqrt( (sumyoverxsq[k]/n_seed) - sqr(avgyoverx[k]) );
	  sigtwist[k] = sqrt( (sumtwistsq[k]/n_seed) - sqr(avgtwist[k]) );
	}
	// finally, dump results into file
	fprintf(fp,"#  k name             s"
		"     mean(sigy/sigx) rms(sigy/sigx) max(sigy/sigx)"
		" mean(twist) rms(twist) max(abs(twist))\n");
	fprintf(fp,"#                     [m]"
		"   []           []           []"
		"           [deg]        [deg]        [deg]\n");
	fprintf(fp,"#\n");
	for(k = 0; k <= globval.Cell_nLoc; k++){
	  fprintf(fp,"%4d %10s %6.3f %e %e %e %e %e %e\n",
		  k, Cell[k].Elem.PName, Cell[k].S,
		  avgyoverx[k], sigyoverx[k], maxyoverx[k],
		  avgtwist[k], sigtwist[k], maxtwist[k]);
	}
	fclose(fp);
      }
    }
    /////////////////////////////////////////////////////////////////////
    
    
    if (cod) {
      
      for (k = 0; k < 2; k++)
	ncorr[k] = 0;
      
      for (j = 0; j <= globval.Cell_nLoc; j++){ // read back beam pos at bpm
	
	if ( (Cell[j].Elem.Pkind == Mpole) || (Cell[j].Elem.Pkind == drift) ) {
	  for (k = 0; k < 6; k++) {
	    x1[j][k] += Cell[j].BeamPos[k];
	    x2[j][k] += sqr(Cell[j].BeamPos[k]);
	  }
	  
	  if ( (Cell[j].Fnum == ElemIndex("corr_h") ) || (Cell[j].Fnum == ElemIndex("corr_v")) ){ // read back corr strength at corr
	    get_bnL_design_elem(Cell[j].Fnum, Cell[j].Knum, Dip, b1L, a1L);
	    if ( Cell[j].Fnum == ElemIndex("corr_h") ){
	      ncorr[X_]++;
	      theta1[ncorr[X_]-1][X_] += b1L;
	      theta2[ncorr[X_]-1][X_] += sqr(b1L);
	    } else if ( Cell[j].Fnum == ElemIndex("corr_v") ){
	      ncorr[Y_]++;
	      theta1[ncorr[Y_]-1][Y_] += a1L;
	      theta2[ncorr[Y_]-1][Y_] += sqr(a1L);
	    }
	  }
	}
      }
      
    } else
      printf("error: orb_corr: failed\n");    
  }
  
  for (k = 0; k < 2; k++)
    ncorr[k] = 0;
  
  for (j = 0; j <= globval.Cell_nLoc; j++){
    
    if ( (Cell[j].Elem.Pkind == Mpole) || (Cell[j].Elem.Pkind == drift) ) {
      for (k = 0; k < 6; k++) {
	x_mean[j][k] = x1[j][k]/n_seed;
	x_sigma[j][k] = sqrt((n_seed*x2[j][k]-sqr(x1[j][k]))
			     /(n_seed*(n_seed-1.0)));
      }
      
      if ( Cell[j].Fnum == ElemIndex("corr_h") ){
	ncorr[X_]++;
	theta_mean[ncorr[X_]-1][X_] = theta1[ncorr[X_]-1][X_]/n_seed;
	theta_sigma[ncorr[X_]-1][X_] = sqrt((n_seed*theta2[ncorr[X_]-1][X_]-sqr(theta1[ncorr[X_]-1][X_]))
					    /(n_seed*(n_seed-1.0)));
	
      } else if ( Cell[j].Fnum == ElemIndex("corr_v") ){
	ncorr[Y_]++;
	theta_mean[ncorr[Y_]-1][Y_] = theta1[ncorr[Y_]-1][Y_]/n_seed;
	theta_sigma[ncorr[Y_]-1][Y_] = sqrt((n_seed*theta2[ncorr[Y_]-1][Y_]-sqr(theta1[ncorr[Y_]-1][Y_]))
					    /(n_seed*(n_seed-1.0)));
      }
    }
  }
  
  prt_cod_rms_data("cod_rms.out", x_mean, x_sigma, theta_mean, theta_sigma);
  // this writes the statictics over N seeds to file
  
}




/////////////////////////////////////////////////////////////////////////////////////////////////




double get_dynap_scl(const double delta, const int n_track2)
{
  char      str[max_str];
  int       i;
  double    x_aper[n_aper], y_aper[n_aper], DA;
  FILE      *fp;

  const int  prt = true;

  fp = file_write("dynap.out");
  dynap(fp, 5e-3, 0.0, 0.1e-3, n_aper, n_track2, x_aper, y_aper, false, prt,
	true);
  fclose(fp);
  DA = get_aper(n_aper, x_aper, y_aper);

  if (true) {
    sprintf(str, "dynap_dp%3.1f.out", 1e2*delta);
    fp = file_write(str);
    dynap(fp, 5e-3, delta, 0.1e-3, n_aper, n_track2,
	  x_aper, y_aper, false, prt, true);
    fclose(fp);
    DA += get_aper(n_aper, x_aper, y_aper);

    for (i = 0; i < nv_; i++)
      globval.CODvect[i] = 0.0;
    sprintf(str, "dynap_dp%3.1f.out", -1e2*delta);
    fp = file_write(str);
    dynap(fp, 5e-3, -delta, 0.1e-3, n_aper,
	  n_track2, x_aper, y_aper, false, prt, true);
    fclose(fp);
    DA += get_aper(n_aper, x_aper, y_aper);
  }

  return DA/3.0;
}




/////////////////////////////////////////////////////////////////////////////////////////////////




void get_matching_params_scl()
{
  double nux, nuy, betax, betay;
  
  nux = globval.TotalTune[X_];
  nuy = globval.TotalTune[Y_];
  betax = globval.OneTurnMat[0][1]/sin(2*pi*nux);
  betay = globval.OneTurnMat[2][3]/sin(2*pi*nuy);
  
  printf("\n");
  printf("beta_x* = %10.9e\n", betax);
  printf("beta_y* = %10.9e\n", betay);
  printf("\n");
}




/////////////////////////////////////////////////////////////////////////////////////////////////




void prt_ZAP()
{
  long int  k;
  FILE      *outf;

  outf = file_write("ZAPLAT.DAT");

  fprintf(outf, "%ld %7.5f\n", globval.Cell_nLoc+1, Cell[globval.Cell_nLoc].S);
  fprintf(outf, "One super period\n");
  for (k = 0; k <= globval.Cell_nLoc; k++) {
    fprintf(outf, "%10.5f %6.3f %7.3f %6.3f %7.3f %6.3f %6.3f %5.3f\n",
	    Cell[k].S,
	    Cell[k].Beta[X_], Cell[k].Alpha[X_],
	    Cell[k].Beta[Y_], Cell[k].Alpha[Y_],
	    Cell[k].Eta[X_], Cell[k].Etap[X_], Cell[k].maxampl[X_][1]);
  }

  fprintf(outf, "0\n");

  fclose(outf);
}



/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////




int main(int argc, char *argv[])
{
  long int  lastpos;

  const long  seed = 1121;

  iniranf(seed); setrancut(2.0);

  // turn on globval.Cavity_on and globval.radiation to get proper synchr radiation damping
  // IDs accounted too if: wiggler model and symplectic integrator (method = 1)
  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; 
  globval.pathlength = false; //globval.bpm         = 0;
  
  
  // overview, on energy: 25-12
  //const double  x_max_FMA = 25e-3, y_max_FMA = 12e-3;
  //const int     n_x = 40, n_y = 40, n_tr = 2046;  // M5: ~24h
  // overview, on energy: 50-20 (to show island)
  //const double  x_max_FMA = 50e-3, y_max_FMA = 20e-3;
  //const int     n_x = 40, n_y = 40, n_tr = 2046;  // M5: ~5-6h
  // overview, off energy: 7-3
  //const double  x_max_FMA = 3e-3, delta_FMA = 7e-2;
  //const int     n_x = 28, n_dp = 56, n_tr = 2046;  // M5: ~16h
  // zoom, off energy: 4.5% - 1.5mm
  //const double  x_max_FMA = 1.5e-3, delta_FMA = 4.5e-2;
  //const int     n_x = 28, n_dp = 56, n_tr = 2046;
 
  // zoom, on energy: 8-2.5
  //const double  x_max_FMA = 8e-3, y_max_FMA = 2.5e-3;
  //const int     n_x = 64, n_y = 15, n_tr = 2046;
  // zoom, off energy: 7-3
  //const double  x_max_FMA = 3e-3, delta_FMA = 7e-2;
  //const double  x_max_FMA = 1.5e-3, delta_FMA = 5e-2;
  //const int     n_x = 28/2, n_dp = 56/2, n_tr = 2046;

  //const double  x_max_FMA = 20e-3, y_max_FMA = 8e-3;
  //const int     n_x = 40, n_y = 40, n_tr = 2046;  // M4: 7h
  //const int     n_x = 20, n_y = 20, n_tr = 2046;  // ~6h ???
  //const int     n_x = 10, n_y = 10, n_tr = 2046;  // ~1h
  //const double  x_max_FMA = 3e-3, delta_FMA = 7e-2;
  //const int     n_x = 28, n_dp = 56, n_tr = 2046;  // M4: 43h
  //const int     n_x = 14, n_dp = 28, n_tr = 2046;  // M4: 12h
  const double  x_max_FMA = 11e-3, y_max_FMA = 4e-3;
  const int     n_x = 40, n_y = 40, n_tr = 2046;  // M4: 65h
  //const int     n_x = 20, n_y = 20, n_tr = 2046;  // M4: 14h
  //const double  x_max_FMA = 1.5e-3, delta_FMA = 5e-2;
  //const int     n_x = 28, n_dp = 56, n_tr = 2046;  // M4: 43h
  //const int     n_x = 14, n_dp = 28, n_tr = 2046;  // M4: 11.5h



  if (true)
    Read_Lattice(argv[1]); //sets some globval params
  else
    rdmfile("flat_file.dat"); //instead of reading lattice file, get data from flat file

  //no_sxt(); //turns off sextupoles

  Ring_GetTwiss(true, 0e-2); printglob(); //gettwiss computes one-turn matrix arg=(w or w/o chromat, dp/p)
  //get_matching_params_scl();
  //get_alphac2();
  //GetEmittance(ElemIndex("cav"), true);
  //prt_beamsizes(); // writes beam_envelope.out; contains sigmamatrix(s), theta(s)

  prt_lat("linlat.out", ElemIndex("bpm_m"), true);  //updated from older T3 version
  //prt_lat("linlat.out", 0, true);  //just testing...

  //prtmfile("flat_file.dat"); // writes flat file
  //prt_chrom_lat(); //writes chromatic functions into chromlat.out
    
  //globval.Aperture_on = true;
  //LoadApers("/home/simon/projects/in/lattice/Apertures_wSeptum.dat", 1, 1);
  //prt_ZAP(); //writes input file for ZAP

  // to check that there are no extra kicks anywhere (ID models!)
  //getcod(0.0, lastpos);
  //prt_cod("cod.out", ElemIndex("bpm_m"), true);  //updated from older T3 version


  if (false) {
    //globval.bpm = ElemIndex("bpm_m");  // broken in new T3 version
    //globval.hcorr = ElemIndex("corr_h"); globval.vcorr = ElemIndex("corr_v");  // broken in new T3 version
    globval.gs = ElemIndex("GS"); globval.ge = ElemIndex("GE");

    //prints a specific closed orbit with corrector strengths
    //getcod(0.0, lastpos);
    //prt_cod("cod.out", ElemIndex("bpm_m"), true);  //updated from older T3 version
   

    // compute response matrix (needed for OCO)
    gcmat(ElemIndex("bpm_m"), ElemIndex("corr_h"), 1); gcmat(ElemIndex("bpm_m"), ElemIndex("corr_v"), 2);
    // print response matrix (routine in lsoc.cc)
    prt_gcmat(1);  prt_gcmat(2);
    // gets response matrix, does svd, evaluates correction for N seed orbits
    //get_cod_rms_scl(dx, dy, dr, n_seed)
    //get_cod_rms_scl(100e-6, 100e-6, 1.0e-3, 100); //trim coils aren't reset when finished

    // for alignments errors check LoadAlignTol (in nsls_ii_lib.cc) and AlignErr.dat
    //CorrectCOD_N("/home/simon/projects/in/lattice/AlignErr.dat", 3, 1, 1);
    //CorrectCOD_N("/home/simon/projects/in/lattice/AlignErr.m4.20140130.dat", 3, 1, 1);
    // 3rd parameter is no. of error steps (use if orb_corr fails)

    // for field errors check LoadFieldErr(in nsls_ii_lib.cc) and FieldErr.dat
    //LoadFieldErr("/home/simon/projects/in/lattice/FieldErr.unshunted.dat", true, 1.0, true);
    //LoadFieldErr("/home/simon/projects/in/lattice/FieldErr.m4.20140514.dat", true, 1.0, true);
    // for alignments errors check LoadAlignTol (in nsls_ii_lib.cc) and AlignErr.dat
    //LoadAlignTol("/home/simon/projects/in/lattice/AlignErr.required.dat", true, 1.0, true, 1);
    //LoadAlignTol("/home/simon/projects/in/lattice/AlignErr.m4.20140130.dat", true, 1.0, true, 1);
    //finds a specific closed orbit with corrector strengths and prints all to file
    //getcod(0.0, lastpos);
    //prt_cod("cod_err.out", ElemIndex("bpm_m"), true);  //updated from older T3 version
     
    // load alignment errors and field errors, correct orbit, repeat N times, and get statistics
    get_cod_rms_scl_new(50); //trim coils aren't reset when finished
    
    // for aperture limitations use LoadApers (in nsls_ii_lib.cc) and Apertures.dat
    //globval.Aperture_on = true;
    //LoadApers("/home/simon/projects/in/lattice/Apertures.dat", 1, 1);
    
  }
  
  if (false) {
    cout << endl;
    cout << "computing tune shifts" << endl;
    dnu_dA(12e-3, 5e-3, 0.0, 25);  // the final argument 25 defines the number of amplitude steps
    get_ksi2(delta); // this gets the chromas and writes them into chrom2.out
  }
  
  if (false) {
    fmap(n_x, n_y, n_tr, x_max_FMA, y_max_FMA, 0.0, true, false);
    //fmapdp(n_x, n_dp, n_tr, x_max_FMA, -delta_FMA, 1e-3, true, false); // always use -delta_FMA (+delta_FMA appears broken)
  } else {
    globval.Cavity_on = true; // this gives longitudinal motion
    globval.radiation = false; // this adds ripple around long. ellipse (needs many turns to resolve damp.)
    //globval.Aperture_on = true;
    //LoadApers("/home/simon/projects/in/lattice/Apertures_wSeptum.dat", 1, 1);
    //get_dynap_scl(delta, 512);
  }
  
  
  //
  // IBS & TOUSCHEK
  //
  int     k, n_turns;
  double  sigma_s, sigma_delta, tau, eps[3];
  FILE    *outf;
  const double  Qb   = 5e-9;
  
  if (true) {
    double  sum_delta[globval.Cell_nLoc+1][2];
    double  sum2_delta[globval.Cell_nLoc+1][2];
    
    GetEmittance(ElemIndex("cav"), true);
    
    // initialize momentum aperture arrays
    for(k = 0; k <= globval.Cell_nLoc; k++){
      sum_delta[k][0] = 0.0; sum_delta[k][1] = 0.0;
      sum2_delta[k][0] = 0.0; sum2_delta[k][1] = 0.0;
    }
    
    globval.eps[X_] = 0.320e-9;
    globval.eps[Y_] = 8e-12;
    sigma_delta     = 0.7727e-03;
    sigma_s         = 8.760e-3;

    // approx. (alpha_z << 1)
    globval.eps[Z_] = sigma_delta*sigma_s;
    globval.alpha_z = 0.0;
    globval.beta_z = sqr(sigma_s)/globval.eps[Z_];


    // INCLUDE LC (LC changes sigma_s and eps_z, but has no influence on sigma_delta)
    if (false) {
      double  newLength, bunchLengthening;
      newLength = 50.0e-3;
      bunchLengthening = newLength/sigma_s;
      sigma_s = newLength;
      globval.eps[Z_] *= bunchLengthening;
      globval.beta_z *= bunchLengthening;  // alpha_z and gamma_z do not change
    }
    
    globval.delta_RF = 7.062e-2; //globval.delta_RF given by cav voltage in lattice file

    //Touschek(Qb, globval.delta_RF, globval.eps[X_], globval.eps[Y_],
    //     sigma_delta, sigma_s);
          

    // IBS
    if (true) {       
      // initialize eps_IBS with eps_SR
      for(k = 0; k < 3; k++)
	eps[k] = globval.eps[k];
      for(k = 0; k < 16; k++){ //prototype (looping because IBS routine doesn't check convergence)
	cout << endl << "*** IBS iteration step " << k << " ***";
	IBS_BM(Qb, globval.eps, eps, true, true);  // use IBS_BM instead of old IBS
      }
    }
    
    
    // TOUSCHEK TRACKING
    if (false) {       
      globval.eps[X_] = 0.320e-9 + 0.0e-12;
      globval.eps[Y_] = 8e-12;
      sigma_delta     = 0.7727e-3;
      sigma_s         = 8.760e-3;

      Touschek(Qb, globval.delta_RF, globval.eps[X_], globval.eps[Y_],
	       sigma_delta, sigma_s);
      
      n_turns = 446; // track for one synchr.osc. -> 1/nu_s (M4 bare @ 1.8MV -> 446, 20130515)
                     //                                     (M5 bare @ 560 kV -> 419)
      
      globval.Aperture_on = true;
      // LoadApers("/home/simon/projects/in/lattice/Apertures.dat", 1, 1);
      
      //globval.delta_RF = 15e-2; //set globval.delta_RF very high to get lattice MA only
      
      tau = Touschek(Qb, globval.delta_RF, false,
		     globval.eps[X_], globval.eps[Y_],
		     sigma_delta, sigma_s,
		     n_turns, true, sum_delta, sum2_delta); //the TRUE flag requires apertures loaded
      
      printf("Touschek lifetime = %10.3e hrs\n", tau/3600.0);
      
      outf = file_write("mom_aper.out");
      for(k = 0; k <= globval.Cell_nLoc; k++)
	fprintf(outf, "%4d %7.2f %5.3f %6.3f\n",
		k, Cell[k].S, 1e2*sum_delta[k][0], 1e2*sum_delta[k][1]);
      fclose(outf);
    }
  
  }
}
