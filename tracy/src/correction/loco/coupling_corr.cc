// Coupling / vertical-dispersion correction. See correction/loco/coupling_corr.h.


corr::coupling_corr::~coupling_corr()
{
  if (SkewRespMat != 0) free_dmatrix(SkewRespMat, 1, N_COUPLE, 1, N_SKEW);
  if (VertCouple != 0) free_dvector(VertCouple, 1, N_COUPLE);
  if (SkewStrengthCorr != 0) free_dvector(SkewStrengthCorr, 1, N_SKEW);
  if (b != 0) free_dvector(b, 1, N_COUPLE);
  if (w != 0) free_dvector(w, 1, N_SKEW);
  if (V != 0) free_dmatrix(V, 1, N_SKEW, 1, N_SKEW);
  if (U != 0) free_dmatrix(U, 1, N_COUPLE, 1, N_SKEW);
  if (eta_y != 0) free_dvector(eta_y, 1, N_BPM);
}


// Read eta values from the file
void corr::coupling_corr::read_eta(const char *TolFileName)
{
  char    line[128], Name[32];
  int     j;
  double  dx, dr;
  FILE    *tolfile;

  tolfile = file_read(TolFileName);

  do
    fgets(line, 128, tolfile);
  while (strstr(line, "#") != NULL);

  printf("\nReading target eta values from file %s:\n",TolFileName);
  j=1;
  do {
    if (strstr(line, "#") == NULL) {
      sscanf(line,"%s %lf %lf", Name, &dx, &dr);
      if (j <= N_BPM) {
	eta_y[j]=dr/1e3;
	printf("%d %s %f %f\n", j, Name, dx, dr);
      } else {
	printf("GetEta: number of BPMs exceeded %d %d\n",N_BPM,j);
	exit_(1);
      }
      j++;
    }
  } while (fgets(line, 128, tolfile) != NULL);
  if (j <= N_BPM) {
    printf("GetEta: number of BPMs too small %d %d\n",N_BPM,j);
    exit_(1);
  }
  fclose(tolfile);
  printf("\n");
}


// "LOCO" for off-diagonal.
void corr::coupling_corr::find_model_matrix(const coupling_cfg &cfg,
					    const double deta_y_max,
					    const double deta_y_offset)
{
  //  Ring_GetTwiss(true, 0.0) should be called in advance
  int      i, j, k;
  long int loc;
  double   nuX, nuY, alpha, eta_y_min, eta_y_max;
  double   *etaSQ;
  double   **betaSQ, **nuSQ, **betaBPM, **nuBPM;
  double   **betaHC, **nuHC, **betaVC, **nuVC;
  FILE     *SkewMatFile, *fp;

  const int    Xi = 1, Yi = 2;
  const double pi = M_PI, twopi = 2.0*M_PI;


  etaSQ = dvector(1, N_SKEW); betaSQ = dmatrix(1, N_SKEW, 1, 2);
  nuSQ = dmatrix(1, N_SKEW, 1, 2);
  betaBPM = dmatrix(1, N_BPM, 1, 2); nuBPM = dmatrix(1, N_BPM, 1, 2);
  betaHC = dmatrix(1, N_HCOR, 1, 2); nuHC = dmatrix(1, N_HCOR, 1, 2);
  betaVC = dmatrix(1, N_VCOR, 1, 2); nuVC = dmatrix(1, N_VCOR, 1, 2);

  nuX = globval.TotalTune[X_]; nuY = globval.TotalTune[Y_];

  for (i = 1; i <= N_SKEW; i++) {
    loc = Elem_GetPos(globval.qt, i);
    etaSQ[i] = Cell[loc].Eta[X_];
    betaSQ[i][Xi] = Cell[loc].Beta[X_]; betaSQ[i][Yi] = Cell[loc].Beta[Y_];
    nuSQ[i][Xi] = Cell[loc].Nu[X_]; nuSQ[i][Yi] = Cell[loc].Nu[Y_];
  } // for i=1..N_SKEW

  for (i = 1; i <= N_BPM; i++) {
    betaBPM[i][Xi] = Cell[bpm_loc[i-1]].Beta[X_];
    betaBPM[i][Yi] = Cell[bpm_loc[i-1]].Beta[Y_];
    nuBPM[i][Xi] = Cell[bpm_loc[i-1]].Nu[X_];
    nuBPM[i][Yi] = Cell[bpm_loc[i-1]].Nu[Y_];
  } // for i=1..N_BPM

  for (i = 1; i <= N_HCOR; i++) {
    betaHC[i][Xi] = Cell[h_corr[i-1]].Beta[X_];
    betaHC[i][Yi] = Cell[h_corr[i-1]].Beta[Y_];
    nuHC[i][Xi] = Cell[h_corr[i-1]].Nu[X_];
    nuHC[i][Yi] = Cell[h_corr[i-1]].Nu[Y_];
  } // for i=1..N_HCOR

  for (i = 1; i <= N_VCOR; i++) {
    betaVC[i][Xi] = Cell[v_corr[i-1]].Beta[X_];
    betaVC[i][Yi] = Cell[v_corr[i-1]].Beta[Y_];
    nuVC[i][Xi] = Cell[v_corr[i-1]].Nu[X_];
    nuVC[i][Yi] = Cell[v_corr[i-1]].Nu[Y_];
  } // for i=1..N_VCOR


  for (i = 1; i <= N_SKEW; i++) {
    // looking for term for vertical dispersion
    alpha = etaSQ[i];
    // printf("For skew quad %3d kick is %9.2e\n",i,alpha);
    for (j = 1; j <= N_BPM; j++) {
      SkewRespMat[j][i] =
	cfg.VDweight*0.5*alpha*sqrt(betaSQ[i][Yi]*betaBPM[j][Yi])
	*cos(twopi*fabs(nuSQ[i][Yi]-nuBPM[j][Yi])-pi*nuY)/sin(pi*nuY);
    } // for (j=1; j<=N_BPM; j++)

    // looking for coupling of horizontal trim to vertical BPM
    for (k = 1; k <= N_HCOR; k++) {
      // find v-kick by i-th skew quad due to the k-th h-trim
      alpha = 0.5*sqrt(betaSQ[i][Xi]*betaHC[k][Xi])*
	cos(twopi*fabs(nuSQ[i][Xi]-nuHC[k][Xi])-pi*nuX)/sin(pi*nuX);
      // find vertical orbit due to the kick
      for (j = 1; j <= N_BPM; j++)
	// Block stride is N_BPM, not N_HCOR: each h-trim contributes a full
	// N_BPM-long orbit vector, and skew_stat's readers assume that layout.
	SkewRespMat[N_BPM+(k-1)*N_BPM+j][i] =
          cfg.HVweight*0.5*alpha*sqrt(betaSQ[i][Yi]*betaBPM[j][Yi])*
	  cos(twopi*fabs(nuSQ[i][Yi]-nuBPM[j][Yi])-pi*nuY)/sin(pi*nuY);
    } //for (k=1; k<=N_HCOR; k++)

   //loking for coupling of vertical trim to horizontal BPM
    for (k = 1; k <= N_VCOR; k++) {
      // find h-kick by i-th skew quad due to the k-th v-trim
      alpha = 0.5*sqrt(betaSQ[i][Yi]*betaVC[k][Yi])*
	cos(twopi*fabs(nuSQ[i][Yi]-nuVC[k][Yi])-pi*nuY)/sin(pi*nuY);
      // find horizontal orbit due to the kick
      for (j = 1; j <= N_BPM; j++)
	// Block stride is N_BPM, not N_VCOR (see the h-trim block above).
	SkewRespMat[N_BPM+N_BPM*N_HCOR+(k-1)*N_BPM+j][i] =
          cfg.VHweight*0.5*alpha*sqrt(betaSQ[i][Xi]*betaBPM[j][Xi])*
	  cos(twopi*fabs(nuSQ[i][Xi]-nuBPM[j][Xi])-pi*nuX)/sin(pi*nuX);
    } //for (k=1; k<=N_VCOR; k++)
  } // for i=1..N_SKEW

  SkewMatFile = file_write(SkewMatFileName);
  for (i = 1; i <= N_SKEW; i++) {
    for (j = 1; j <= N_COUPLE; j++)
      fprintf(SkewMatFile, "%9.2e ", SkewRespMat[j][i]);
    fprintf(SkewMatFile, "\n");
  }
  fclose(SkewMatFile);

  fp = file_write(deta_y_FileName);

  if (deta_y_max < 0.) {
    read_eta("eta_file.dat");
  }
  eta_y_max = -1e8;
  eta_y_min =  1e8;
  for (j = 1; j <= N_BPM; j++) {
    if (deta_y_max > 0.) {
      eta_y[j] = 0.0;
      for (i = 1; i <= N_SKEW; i++)
	if (i % cfg.SQ_per_scell == 0) {
	  eta_y[j] += 0.5*etaSQ[i]*sqrt(betaSQ[i][Yi]*betaBPM[j][Yi])
	      *cos(twopi*fabs(nuSQ[i][Yi]-nuBPM[j][Yi])-pi*nuY)/sin(pi*nuY);
	}
    }
    eta_y_max = max(eta_y[j], eta_y_max);
    eta_y_min = min(eta_y[j], eta_y_min);
  }
  for (j = 1; j <= N_BPM; j++) {
    eta_y[j] = (eta_y[j] - eta_y_min)/(eta_y_max - eta_y_min);
    eta_y[j]+=deta_y_offset;
  }

  fprintf(fp, "# nbpm %d SQ_per_scell %d etaymin %10.3e mm etaymax %10.3e mm"
	  " detaymax %10.3e mm detayoffset %10.3e mm\n",
	  N_BPM, cfg.SQ_per_scell, 1e3*eta_y_min, 1e3*eta_y_max,
	  1e3*deta_y_max, 1e3*deta_y_max*deta_y_offset);
  for (j = 1; j <= N_BPM; j++) {
    eta_y[j] *= fabs(deta_y_max);
    fprintf(fp, "%6.3f %10.3e\n", Cell[bpm_loc[j-1]].S, 1e3*eta_y[j]);
  }
  fclose(fp);

  free_dvector(etaSQ, 1, N_SKEW); free_dmatrix(betaSQ, 1, N_SKEW, 1, 2);
  free_dmatrix(nuSQ, 1, N_SKEW, 1, 2);
  free_dmatrix(betaBPM, 1, N_BPM, 1, 2); free_dmatrix(nuBPM, 1, N_BPM, 1, 2);
  free_dmatrix(betaHC, 1, N_HCOR, 1, 2); free_dmatrix(nuHC, 1, N_HCOR, 1, 2);
  free_dmatrix(betaVC, 1, N_VCOR, 1, 2); free_dmatrix(nuVC, 1, N_VCOR, 1, 2);
} // find_model_matrix


void corr::coupling_corr::ini_skew_cor(const coupling_cfg &cfg,
				       const double deta_y_max,
				       const double deta_y_offset)
{
  // Collect skew trims from param file "qt"
  int k;

  // No of skew quads, BPMs, and correctors
  N_SKEW = GetnKid(globval.qt);

  N_BPM = 0;
  for (k = 1; k <= GetnKid(globval.bpm); k++) {
    N_BPM++;

    if (N_BPM > max_bpm) {
      printf("ini_skew_cor: max no of BPMs exceeded %d (%d)\n",
	     N_BPM, max_bpm);
      exit_(1);
    }

    bpm_loc[N_BPM-1] = Elem_GetPos(globval.bpm, k);
  }

  N_HCOR = 0;
  h_corr[N_HCOR++] = Elem_GetPos(globval.hcorr, 1);
  h_corr[N_HCOR++] = Elem_GetPos(globval.hcorr, GetnKid(globval.hcorr)/3);
  h_corr[N_HCOR++] = Elem_GetPos(globval.hcorr, 2*GetnKid(globval.hcorr)/3);

  N_VCOR = 0;
  v_corr[N_VCOR++] = Elem_GetPos(globval.vcorr, 1);
  v_corr[N_VCOR++] = Elem_GetPos(globval.vcorr, GetnKid(globval.vcorr)/3);
  v_corr[N_VCOR++] = Elem_GetPos(globval.vcorr, 2*GetnKid(globval.vcorr)/3);

  N_COUPLE = N_BPM*(1+N_HCOR+N_VCOR);

  SkewRespMat = dmatrix(1, N_COUPLE, 1, N_SKEW);
  VertCouple = dvector(1, N_COUPLE);
  SkewStrengthCorr = dvector(1, N_SKEW);
  b = dvector(1, N_COUPLE); w = dvector(1, N_SKEW);
  V = dmatrix(1, N_SKEW, 1, N_SKEW); U = dmatrix(1, N_COUPLE, 1, N_SKEW);
  eta_y = dvector(1, N_BPM);

  printf("\n");
  printf("Number of trims:                   horizontal = %d, vertical = %d\n",
	 N_HCOR, N_VCOR);
  printf("Number of BPMs:                    %6d\n", N_BPM);
  printf("Number of skew quads:              %6d\n", N_SKEW);
  printf("Number of elements in skew vector: %6d\n", N_COUPLE);

  // find matrix
  Ring_GetTwiss(true, 0.0);

  printf("\n");
  printf("Looking for response matrix\n");
  find_model_matrix(cfg, deta_y_max, deta_y_offset);

  printf("Looking for SVD matrices\n");
  corr::svd_decomp_cut(SkewRespMat, N_COUPLE, N_SKEW, U, w, V, cfg.qt_s_cut,
		       true);
}


void corr::coupling_corr::find_coup_vector(const coupling_cfg &cfg,
					   double *VertCouple)
{
  int    i, j;
  double *resp;

  resp = dvector(1, N_BPM);

  // Find vertical dispersion
  Cell_Geteta(0, globval.Cell_nLoc, true, 0e0);

  for (i = 1; i <= N_BPM; i++)
    VertCouple[i] = cfg.VDweight*Cell[bpm_loc[i-1]].Eta[Y_];
  // Finished finding vertical dispersion

  // Off diagonal terms for horizontal trims: kick with "+Dip", read y_.
  for (j = 1; j <= N_HCOR; j++) {
    corr::measure_orm_column(Cell[h_corr[j-1]].Fnum, Cell[h_corr[j-1]].Knum,
			     +Dip, cfg.kick, bpm_loc, N_BPM, y_, resp);

    for (i = 1; i <= N_BPM; i++)
      VertCouple[N_BPM+(j-1)*N_BPM+i] = -cfg.HVweight*resp[i]; // sign reversal
  } // hcorr cycle


  // Off diagonal terms for vertical trims: kick with "-Dip", read x_.
  for (j = 1; j <= N_VCOR; j++) {
    corr::measure_orm_column(Cell[v_corr[j-1]].Fnum, Cell[v_corr[j-1]].Knum,
			     -Dip, cfg.kick, bpm_loc, N_BPM, x_, resp);

    for (i = 1; i <= N_BPM; i++)
      VertCouple[N_BPM+N_BPM*N_HCOR+(j-1)*N_BPM+i] = cfg.VHweight*resp[i];
  } // vcorr cycle

  free_dvector(resp, 1, N_BPM);
} // find_coup_vector


void corr::coupling_corr::skew_stat(const coupling_cfg &cfg,
				    double VertCouple[], const int cnt)
{
  int    i;
  FILE *outf = NULL;
  char fname[30];

  double max, mean, rms, sk;

  if (cnt>=0) {
    snprintf(fname, sizeof(fname), "%s_%d.out",skew_FileName,cnt);
    outf = file_write(fname);
    fprintf(outf, "# qt s [m] etax [m] name kl [1/m]\n");
  }

  // statistics for skew quadrupoles
  max = 0.0; rms = mean = 0.0;
  for(i = 1; i <= N_SKEW; i++) {
    sk = GetKLpar(globval.qt, i, -Quad);
    if (cnt>=0)
      fprintf(outf, "%4d %7.3f %12.5e %s %12.5e\n",
	       i,Cell[Elem_GetPos(globval.qt,i)].S,
	       Cell[Elem_GetPos(globval.qt,i)].Eta[X_],
	       Cell[Elem_GetPos(globval.qt,i)].Elem.PName,sk);
    if (fabs(sk) > max) max = fabs(sk);
    rms += sqr(sk);
    mean += sk;
  }
  mean = mean/N_SKEW;
  rms = sqrt(-mean*mean+rms/N_SKEW);

  if (cnt>=0)
    fprintf(outf,"# Max Mean Rms skew strength: %8.2e/%8.2e+/-%8.2e 1/m\n",
	    max, mean, rms);
  else
    printf("Max Mean Rms skew strength: %8.2e/%8.2e+/-%8.2e 1/m\n",
	   max, mean, rms);

  // statistics for vertical dispersion function
  max = 0.0; rms = mean = 0.0;
  for(i = 1; i <= N_BPM; i++) {
    if (fabs(VertCouple[i]/cfg.VDweight) > max)
      max = fabs(VertCouple[i]/cfg.VDweight);
    rms += sqr(VertCouple[i]/cfg.VDweight);
    mean += VertCouple[i]/cfg.VDweight;
  }
  mean = mean/N_BPM;
  rms = sqrt(-mean*mean+rms/N_BPM);
  if (cnt>=0)
    fprintf(outf,
	    "# Max Mean Rms vertical dispersion: %8.2e/%8.2e+/-%8.2e mm\n",
	    1e3*max, 1e3*mean, 1e3*rms);
  else
    printf("Max Mean Rms vertical dispersion: %8.2e/%8.2e+/-%8.2e mm\n",
	   1e3*max, 1e3*mean,1e3*rms);

  // statistics for off diagonal terms of response matrix (trims->bpms)
  max = 0.0; rms = mean = 0.0;
  for(i = N_BPM+1; i <= N_BPM*(1+N_HCOR); i++) {
    if (fabs(VertCouple[i]/cfg.HVweight) > max)
      max = fabs(VertCouple[i]/cfg.HVweight);
    rms += sqr(VertCouple[i]/cfg.HVweight);
    mean += VertCouple[i]/cfg.HVweight;
  }
  mean = mean/(N_HCOR*N_BPM);
  rms = sqrt(-mean*mean+rms/(N_HCOR*N_BPM));
  if (cnt>=0)
    fprintf(outf,
	    "# Max Mean Rms horizontal coupling: %8.2e/%8.2e+/-%8.2e mm/mrad\n",
	    max, mean, rms);
  else
    printf("Max Mean Rms horizontal coupling: %8.2e/%8.2e+/-%8.2e mm/mrad\n",
	   max, mean, rms);

  max = 0.0; rms = mean = 0.0;
  for(i = N_BPM*(1+N_HCOR)+1; i <= N_COUPLE; i++) {
    if (fabs(VertCouple[i]/cfg.VHweight) > max)
      max = fabs(VertCouple[i]/cfg.VHweight);
    rms += sqr(VertCouple[i]/cfg.VHweight);
    mean += VertCouple[i]/cfg.VHweight;
  }
  mean = mean/(N_VCOR*N_BPM);
  rms = sqrt(-mean*mean+rms/(N_VCOR*N_BPM));
  if (cnt>=0) {
    fprintf(outf,
	    "# Max Mean Rms vertical coupling: %8.2e/%8.2e+/-%8.2e mm/mrad\n",
	    max, mean, rms);
    fclose(outf);
  } else
    printf("Max Mean Rms vertical coupling: %8.2e/%8.2e+/-%8.2e mm/mrad\n",
	   max, mean, rms);
}


void corr::coupling_corr::corr_eps_y(const coupling_cfg &cfg, const int cnt)
{
  int  i, j;
  FILE *outf;
  char fname[30];
  int qtnr;
  double qtpos, qtkl, qteta;
  char qtnam[20];
  FILE *cinf;
  int loc;

  // Clear skew quad setpoints
  set_bnL_design_fam(globval.qt, Quad, 0.0, 0.0);

  // Find coupling vector
  printf("\n");
  printf("Looking for coupling error\n");
  find_coup_vector(cfg, VertCouple);

  //Find and print coupling statistics
  printf("\n");
  printf("Before correction\n");
  skew_stat(cfg, VertCouple, -1);

  // Coupling Correction
  printf("\n");
  for (i = 1; i <= cfg.n_lin; i++) {
    printf("Looking for correction\n");

    //Find Correcting Settings to skew quadrupoles
    for (j = 1; j <= N_BPM; j++)
      b[j] = cfg.VDweight*eta_y[j] - VertCouple[j];

    for (j = N_BPM+1; j <= N_COUPLE; j++)
      b[j] = -VertCouple[j];

    corr::svd_backsub(U, w, V, N_COUPLE, N_SKEW, b, SkewStrengthCorr);

    printf("Applying correction\n");
    // Add correction
    for (j = 1; j <= N_SKEW; j++)
      SetdKLpar(globval.qt, j, -Quad, SkewStrengthCorr[j]);

    printf("\n");
    printf("Looking for coupling error\n");
    // Find coupling vector
    find_coup_vector(cfg, VertCouple);

    printf("\n");
    printf("After run %d of correction\n", i);
    // Find and print coupling statistics
    skew_stat(cfg, VertCouple, -1);

  } // End of coupling Correction

  if (cfg.qt_from_file) {
    printf("\n");
    printf("Reading skew quad values from file 'qt_file.dat':\n");
    printf("\n");
    cinf = fopen("qt_file.dat" , "r");
    for(j = 1; j <= N_SKEW; j++) {
      fscanf(cinf, "%d %lg %lg %s %lg", &qtnr, &qtpos, &qteta, qtnam, &qtkl);
      printf("%d %d %s %f\n", j, qtnr, qtnam, qtkl);
      loc = Elem_GetPos(globval.qt, j);
      printf("%d %e %e %e %e\n",
	     loc, Cell[loc].Elem.PL,
	     qtkl,Cell[loc].Elem.M->PBpar[-Quad + HOMmax],
	     Cell[loc].Elem.M->PBpar[Quad + HOMmax]);
      SetdKLpar(globval.qt, j, -Quad, qtkl);
      loc = Elem_GetPos(globval.qt, j);
      printf("%d %e %e %e %e\n",
	     loc, Cell[loc].Elem.PL, qtkl,
	     Cell[loc].Elem.M->PBpar[-Quad + HOMmax],
	     Cell[loc].Elem.M->PBpar[Quad + HOMmax]);
    }
    printf("\n");
    Ring_GetTwiss(true, 0.0); printglob();
    printf("\n");
    printf("Looking for coupling error\n");
    // Find coupling vector
    find_coup_vector(cfg, VertCouple);
    printf("\n");
    printf("After application of skew values from file 'qt_file.dat'\n");
    skew_stat(cfg, VertCouple, -1);
  }

  skew_stat(cfg, VertCouple, cnt);

  snprintf(fname, sizeof(fname), "%s_%d.out",eta_y_FileName,cnt);
  outf = file_write(fname);

  fprintf(outf, "# nr s [m] name nuy etay [mm] etapy [mrad]\n");
  for (i = 0; i <= globval.Cell_nLoc; i++)
    fprintf(outf, "%4d %7.3f %s %6.3f %10.3e %10.3e\n",
	    i, Cell[i].S, Cell[i].Elem.PName,
	    Cell[i].Nu[Y_], 1e3*Cell[i].Eta[Y_], 1e3*Cell[i].Etap[Y_]);
  fclose(outf);

  find_coup_vector(cfg, VertCouple);
}
