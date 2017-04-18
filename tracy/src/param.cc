
bool        param_data_type::DA_bare      = false;
bool        param_data_type::freq_map     = false;
int         param_data_type::n_orbit      = 5;
int         param_data_type::n_scale      = 1;
std::string param_data_type::loc_Fam_name = "";
int         param_data_type::n_cell       = -1;

int param_data_type::n_lin         =  3;
int param_data_type::SQ_per_scell  =  2;
int param_data_type::BPM_per_scell = 12;
int param_data_type::HCM_per_scell = 12;
int param_data_type::VCM_per_scell = 12;

double param_data_type::kick       = 0.01e-3;
int    param_data_type::n_stat     = 1;

double param_data_type::VDweight   = 1e3,
       param_data_type::HVweight   = 1e0,
       param_data_type::VHweight   = 1e0;

int    param_data_type::n_track_DA = 512,
       param_data_type::n_aper_DA  = 15,
       param_data_type::n_delta_DA = 12;
double param_data_type::delta_DA   = 3e-2;

int    param_data_type::n_x        = 50,
       param_data_type::n_y        = 30,
       param_data_type::n_dp       = 25,
       param_data_type::n_tr       = 2064;
double param_data_type::x_max_FMA  = 20e-3,
       param_data_type::y_max_FMA  = 6e-3,
       param_data_type::delta_FMA  = 3e-2;

bool   param_data_type::bba        = false;


void param_data_type::get_param(const string &param_file)
{
  char              *s, name[max_str], line[max_str], str[max_str], *p;
  string            lat_file, flat_file;
  double            f_prm;
  std::ifstream     inf;
  std::stringstream sstr;

  const bool  prt = true;

  if (prt) {
    std::cout << std::endl;
    std::cout << "get_param: " << param_file << std::endl;
  }

  file_rd(inf, param_file.c_str());

  // read param file
  ae_file = ""; fe_file = ""; ap_file = "";

  if (prt) std::cout << std::endl;

  while (!inf.getline(line, max_str).eof()) {
    if (prt) std::cout << line << std::endl;

    if (strstr(line, "#") == NULL) {
      sscanf(line, "%s", name);

      if (strcmp("energy", name) == 0) {
	sscanf(line, "%*s %lf", &globval.Energy);
      } else if (strcmp("in_dir", name) == 0){
        sscanf(line, "%*s %s", str);
	in_dir = str;
      } else if (strcmp("ae_file", name) == 0){
        sscanf(line, "%*s %s", str);
	sstr.clear(); sstr.str("");
	sstr << in_dir << str; ae_file = sstr.str();
      } else if (strcmp("fe_file", name) == 0) {
        sscanf(line, "%*s %s", str);
	sstr.clear(); sstr.str("");
	sstr << in_dir << str; fe_file = sstr.str();
      } else if (strcmp("ap_file", name) == 0) {
        sscanf(line, "%*s %s", str);
	sstr.clear(); sstr.str("");
	sstr << in_dir << str; ap_file = sstr.str();
      } else if (strcmp("lat_file", name) == 0) {
        sscanf(line, "%*s %s", str);
        sstr.clear(); sstr.str("");
	sstr << in_dir << str; lat_FileName = sstr.str();
        Read_Lattice(lat_FileName.c_str());
      } else if (strcmp("flat_file", name) == 0) {
	sscanf(line, "%*s %s", str);
	sstr.clear(); sstr.str("");
	sstr << str << "flat_file.dat"; flat_file = sstr.str();
	rdmfile(flat_file.c_str());
      } else if (strcmp("s_cut", name) == 0) {
	sscanf(line, "%*s %lf", &f_prm);
	setrancut(f_prm);
      } else if (strcmp("n_stat", name) == 0)
	sscanf(line, "%*s %d", &n_stat);
      else if (strcmp("n_aper", name) == 0)
	sscanf(line, "%*s %d", &n_aper_DA);
      else if (strcmp("loc_name", name) == 0) {
        sscanf(line, "%*s %s", str);
	sstr.clear(); sstr.str(""); sstr << str; loc_Fam_name = sstr.str();
      } else if (strcmp("n_cell", name) == 0)
	sscanf(line, "%*s %d", &n_cell);
      else if (strcmp("n_scale", name) == 0)
	sscanf(line, "%*s %d", &n_scale);
      else if (strcmp("n_orbit", name) == 0)
	sscanf(line, "%*s %d", &n_orbit);
      else if (strcmp("bpm_names", name) == 0) {
	strtok_r(line, " \r", &p); s = strtok_r(NULL, " \r", &p);
	while (s != NULL) {
	  bpm_Fam_names.push_back(s); s = strtok_r(NULL, " \r", &p);
	}
      } else if (strcmp("h_corrs", name) == 0) {
	strtok_r(line, " \r", &p); s = strtok_r(NULL, " \r", &p);
	while (s != NULL) {
	    corr_Fam_names[X_].push_back(s); s = strtok_r(NULL, " \r", &p);
	}
      } else if (strcmp("v_corrs", name) == 0) {
	strtok_r(line, " \r", &p); s = strtok_r(NULL, " \r", &p);
	while (s != NULL) {
	  corr_Fam_names[Y_].push_back(s); s = strtok_r(NULL, " \r", &p);
	}
      } else if (strcmp("gs", name) == 0) {
	sscanf(line, "%*s %s", str);
	globval.gs = ElemIndex(str);
      } else if (strcmp("ge", name) == 0) {
	sscanf(line, "%*s %s", str);
	globval.ge = ElemIndex(str);
      } else if (strcmp("DA_bare", name) == 0) {
	sscanf(line, "%*s %s", str);
	DA_bare = (strcmp(str, "true") == 0)? true : false;
      } else if (strcmp("n_track", name) == 0)
	sscanf(line, "%*s %d", &n_track_DA);
      else if (strcmp("n_delta", name) == 0)
	sscanf(line, "%*s %d", &n_delta_DA);
      else if (strcmp("delta", name) == 0)
	sscanf(line, "%*s %lf", &delta_DA);
      else if (strcmp("freq_map", name) == 0) {
	sscanf(line, "%*s %s %d %d %d %d %lf %lf %lf",
	       str, &n_x, &n_y, &n_dp, &n_tr,
	       &x_max_FMA, &y_max_FMA, &delta_FMA);
	freq_map = (strcmp(str, "true") == 0)? true : false;
      } else if (strcmp("qt", name) == 0) {
	sscanf(line, "%*s %s", str);
	globval.qt = ElemIndex(str);
      } else if (strcmp("disp_wave_y", name) == 0)
	sscanf(line, "%*s %lf", &disp_wave_y);
      else if (strcmp("n_lin", name) == 0)
	sscanf(line, "%*s %d", &n_lin);
      else if (strcmp("VDweight", name) == 0)
	sscanf(line, "%*s %lf", &VDweight);
      else if (strcmp("HVweight", name) == 0)
	sscanf(line, "%*s %lf", &HVweight);
      else if (strcmp("VHweight", name) == 0)
	sscanf(line, "%*s %lf", &VHweight);
      else if (strcmp("N_calls", name) == 0) // ID correction parameters
	sscanf(line, "%*s %d", &N_calls);
      else if (strcmp("N_steps", name) == 0)
	sscanf(line, "%*s %d", &N_steps);
      else if (strcmp("ID_quads", name) == 0) {
	strtok_r(line, " \r", &p); s = strtok_r(NULL, " \r", &p); N_Fam = 0;
	while (s != NULL) {
	  N_Fam++;
	  if (N_Fam <= N_Fam_max) {
	    Q_Fam[N_Fam-1] = ElemIndex(s); s = strtok_r(NULL, " \r", &p);
	  } else {
	    printf("get_param: N_Fam_max exceeded (%d)\n", N_Fam_max);
	    exit(1);
	  }
	}
      } else {
	std::cout << "bad line in " << param_file << ": " << line << std::endl;
        exit_(1);
      }
    }
  }

  inf.close();
}


void param_data_type::get_bare(void)
{
  // Store optics function values at the sextupoles.
  long int j, k;

  n_sext = 0;
  for (j = 0; j <= globval.Cell_nLoc; j++) {
    if ((Cell[j].Elem.Pkind == Mpole) && (Cell[j].Elem.M->n_design >= Sext)) {
      n_sext++; sexts[n_sext-1] = j;
      for (k = 0; k < 2; k++) {
	betas0_[n_sext-1][k] = Cell[j].Beta[k];
	nus0_[n_sext-1][k] = Cell[j].Nu[k];
      }
    }
  }

  nu0_[X_] = globval.TotalTune[X_]; nu0_[Y_] = globval.TotalTune[Y_];
}


void param_data_type::get_dbeta_dnu(double m_dbeta[], double s_dbeta[],
				    double m_dnu[], double s_dnu[])
{
  int       k;
  long int  j, ind;
  double    dbeta, dnu;

  Ring_GetTwiss(false, 0.0);

  for (k = 0; k <= 1; k++) {
    m_dbeta[k] = 0.0; s_dbeta[k] = 0.0; m_dnu[k] = 0.0; s_dnu[k] = 0.0;
  }

  for (j = 0; j < n_sext; j++) {
    ind = sexts[j];
    for (k = 0; k <= 1; k++) {
      dbeta = (Cell[ind].Beta[k]-betas0_[j][k])/betas0_[j][k];
      m_dbeta[k] += dbeta; s_dbeta[k] += sqr(dbeta);
      dnu = Cell[ind].Nu[k] - nus0_[j][k];
      m_dnu[k] += dnu; s_dnu[k] += sqr(dnu);
    }
  }

  for (k = 0; k <= 1; k++) {
    m_dbeta[k] /= n_sext; m_dnu[k] /= n_sext;
    s_dbeta[k] = sqrt((s_dbeta[k]-n_sext*sqr(m_dbeta[k]))/(n_sext-1));
    s_dnu[k] = sqrt((s_dnu[k]-n_sext*sqr(m_dnu[k]))/(n_sext-1));
  }
}


void param_data_type::FindSQ_SVDmat(double **SkewRespMat, double **U,
				    double **V, double *w, int N_COUPLE,
				    int N_SKEW)
{
  int i, j;

  const double  cut = 1e-10; // cut value for SVD singular values

  for (i = 1; i <= N_COUPLE; i++)
    for (j = 1; j <= N_SKEW; j++)
      U[i][j] = SkewRespMat[i][j];

  // prepare matrices for SVD
  dsvdcmp(U, N_COUPLE, N_SKEW, w, V);

  // zero singular values
  printf("\n");
  printf("singular values:\n");
  printf("\n");
  for (i = 1; i <= N_SKEW; i++) {
    printf("%11.3e", w[i]);
    if (w[i] < cut ) {
      w[i] = 0.0;
      printf(" (zeroed)");
      if (i % 8 == 0) printf("\n");
    }
  }
  if (i % 8 != 0) printf("\n");
}


void param_data_type::FindMatrix(double **SkewRespMat, const double deta_y_max)
{
  //  Ring_GetTwiss(true, 0.0) should be called in advance
  int      i, j, k;
  long int loc;
  double   nuX, nuY, alpha, eta_y_max;
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
      SkewRespMat[j][i] = VDweight*0.5*alpha*sqrt(betaSQ[i][Yi]*betaBPM[j][Yi])
	*cos(twopi*fabs(nuSQ[i][Yi]-nuBPM[j][Yi])-pi*nuY)/sin(pi*nuY);
    } // for (j=1; j<=N_BPM; j++)

    // loking for coupling of horizontal trim to vertical BPM
    for (k = 1; k <= N_HCOR; k++) {
      // find v-kick by i-th skew quad due to the k-th h-trim
      alpha = 0.5*sqrt(betaSQ[i][Xi]*betaHC[k][Xi])*
	cos(twopi*fabs(nuSQ[i][Xi]-nuHC[k][Xi])-pi*nuX)/sin(pi*nuX);
      // find vertical orbit due to the kick
      for (j = 1; j <= N_BPM; j++)
	SkewRespMat[N_BPM+(k-1)*N_HCOR+j][i] =
          HVweight*0.5*alpha*sqrt(betaSQ[i][Yi]*betaBPM[j][Yi])*
	  cos(twopi*fabs(nuSQ[i][Yi]-nuBPM[j][Yi])-pi*nuY)/sin(pi*nuY);
    } //for (k=1; k<=N_HCOR; k++)

   //loking for coupling of vertical trim to horizontal BPM
    for (k = 1; k <= N_VCOR; k++) {
      // find h-kick by i-th skew quad due to the k-th v-trim
      alpha = 0.5*sqrt(betaSQ[i][Yi]*betaVC[k][Yi])*
	cos(twopi*fabs(nuSQ[i][Yi]-nuVC[k][Yi])-pi*nuY)/sin(pi*nuY);
      // find horizontal orbit due to the kick
      for (j = 1; j <= N_BPM; j++)
	SkewRespMat[N_BPM+N_BPM*N_HCOR+(k-1)*N_VCOR+j][i] =
          VHweight*0.5*alpha*sqrt(betaSQ[i][Xi]*betaBPM[j][Xi])*
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
  eta_y_max = 0.0;
  for (j = 1; j <= N_BPM; j++) {
    eta_y[j] = 0.0;
    for (i = 1; i <= N_SKEW; i++)
      if (i % SQ_per_scell == 0) {
	eta_y[j] += 0.5*etaSQ[i]*sqrt(betaSQ[i][Yi]*betaBPM[j][Yi])
	    *cos(twopi*fabs(nuSQ[i][Yi]-nuBPM[j][Yi])-pi*nuY)/sin(pi*nuY);
      }
    eta_y_max = max(fabs(eta_y[j]), eta_y_max);
  }
  for (j = 1; j <= N_BPM; j++)
    eta_y[j] /= eta_y_max;

  for (j = 1; j <= N_BPM; j++) {
    eta_y[j] *= deta_y_max;
    fprintf(fp, "%6.3f %10.3e\n", Cell[bpm_loc[j-1]].S, 1e3*eta_y[j]);
  }
  fclose(fp);

  free_dvector(etaSQ, 1, N_SKEW); free_dmatrix(betaSQ, 1, N_SKEW, 1, 2);
  free_dmatrix(nuSQ, 1, N_SKEW, 1, 2);
  free_dmatrix(betaBPM, 1, N_BPM, 1, 2); free_dmatrix(nuBPM, 1, N_BPM, 1, 2);
  free_dmatrix(betaHC, 1, N_HCOR, 1, 2); free_dmatrix(nuHC, 1, N_HCOR, 1, 2);
  free_dmatrix(betaVC, 1, N_VCOR, 1, 2); free_dmatrix(nuVC, 1, N_VCOR, 1, 2);
} // FindMatrix


void param_data_type::ini_skew_cor(const double deta_y_max)
{
  int k;

  std::cout << "ini_skew_cor: out-of-date (globval.hcorr...)" << std::endl;
  exit(1);

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
  FindMatrix(SkewRespMat, deta_y_max);

  printf("Looking for SVD matrices\n");
  FindSQ_SVDmat(SkewRespMat, U, V, w, N_COUPLE, N_SKEW);
}


void param_data_type::FindCoupVector(double *VertCouple)
{
  bool   cod;
  long   i, j;
  long   lastpos;
  double *orbitP, *orbitN;

  orbitP = dvector(1, N_BPM); orbitN = dvector(1, N_BPM);

  // Find vertical dispersion
  Cell_Geteta(0, globval.Cell_nLoc, true, 0e0);

  for (i = 1; i <= N_BPM; i++)
    VertCouple[i] = VDweight*Cell[bpm_loc[i-1]].Eta[Y_];
  // Finished finding vertical dispersion

  // Find off diagonal terms for horizontal trims
  for (j = 1; j <= N_HCOR; j++) {
    // positive kick: "+Dip" for horizontal
    SetdKLpar(Cell[h_corr[j-1]].Fnum, Cell[h_corr[j-1]].Knum, +Dip, kick);
    cod = getcod(0.0, lastpos); chk_cod(cod, "FindCoupVector");
    for (i = 1; i <= N_BPM; i++)
      orbitP[i] = Cell[bpm_loc[i-1]].BeamPos[y_];

    //negative kick: "+Dip" for horizontal
    SetdKLpar(Cell[h_corr[j-1]].Fnum, Cell[h_corr[j-1]].Knum, +Dip, -2*kick);
    cod = getcod(0.0, lastpos); chk_cod(cod, "FindCoupVector");
    for (i = 1; i <= N_BPM; i++)
      orbitN[i] = Cell[bpm_loc[i-1]].BeamPos[y_];

    // restore trim valueL: "+Dip" for horizontal
    SetdKLpar(Cell[h_corr[j-1]].Fnum, Cell[h_corr[j-1]].Knum, +Dip, kick);

    for (i = 1; i <= N_BPM; i++)
      VertCouple[N_BPM+(j-1)*N_HCOR+i] =
	HVweight*(orbitN[i]-orbitP[i])*0.5/kick; // sign reversal
  } // hcorr cycle


  // Find off diagonal terms for vertical trims
  for (j = 1; j <= N_VCOR; j++){
    // positive kick: "-Dip" for vertical
    SetdKLpar(Cell[v_corr[j-1]].Fnum, Cell[v_corr[j-1]].Knum, -Dip, kick);
    cod = getcod(0.0, lastpos); chk_cod(cod, "FindCoupVector");
    for (i = 1;  i <= N_BPM; i++)
      orbitP[i] = Cell[bpm_loc[i-1]].BeamPos[x_];

    // negative kick: "-Dip" for vertical
    SetdKLpar(Cell[v_corr[j-1]].Fnum, Cell[v_corr[j-1]].Knum, -Dip, -2*kick);
    cod = getcod(0.0, lastpos); chk_cod(cod, "FindCoupVector");
    for (i = 1; i <= N_BPM; i++)
      orbitN[i] = Cell[bpm_loc[i-1]].BeamPos[x_];

    // restore corrector: "-Dip" for vertical
    SetdKLpar(Cell[v_corr[j-1]].Fnum, Cell[v_corr[j-1]].Knum, -Dip, kick);

    for (i = 1; i <= N_BPM; i++)
      VertCouple[N_BPM+N_BPM*N_HCOR+(j-1)*N_VCOR+i] =
	VHweight*(orbitP[i]-orbitN[i])*0.5/kick;
  } // vcorr cycle

  free_dvector(orbitP, 1, N_BPM); free_dvector(orbitN, 1, N_BPM);
} // FindCoupVector


void param_data_type::SkewStat(double VertCouple[])
{
  int    i;
  double max, rms, sk;

  // statistics for skew quadrupoles
  max = 0.0; rms = 0.0;
  for(i = 1; i <= N_SKEW; i++) {
    sk = GetKLpar(globval.qt, i, -Quad);
    if (fabs(sk) > max) max = fabs(sk);
    rms += sqr(sk);
  }
  rms = sqrt(rms/N_SKEW);
  printf("Rms skew strength:       %8.2e+/-%8.2e\n", max, rms);

  // statistics for vertical dispersion function
  max = 0.0; rms = 0.0;
  for(i = 1; i <= N_BPM; i++) {
    if (fabs(VertCouple[i]) > max) max = fabs(VertCouple[i]/VDweight);
    rms += sqr(VertCouple[i]/VDweight);
  }
  rms = sqrt(rms/N_BPM);
  printf("Max vertical dispersion: %8.2e+/-%8.2e mm\n", 1e3*max, 1e3*rms);

  // statistics for off diagonal terms of response matrix (trims->bpms)
  max = 0.0; rms = 0.0;
  for(i = N_BPM+1; i <= N_BPM*(1+N_HCOR); i++) {
    if (fabs(VertCouple[i]) > max) max = fabs(VertCouple[i]/HVweight);
    rms += sqr(VertCouple[i]/HVweight);
  }
  rms = sqrt(rms/(N_HCOR*N_BPM));
  printf("Max horizontal coupling: %8.2e+/-%8.2e mm/mrad\n", max, rms);

  max = 0.0; rms = 0.0;
  for(i = N_BPM*(1+N_HCOR)+1; i <= N_COUPLE; i++) {
    if (fabs(VertCouple[i]) > max) max = fabs(VertCouple[i]/VHweight);
    rms += sqr(VertCouple[i]/VHweight);
  }
  rms = sqrt(rms/(N_VCOR*N_BPM));
  printf("Max vertical coupling:   %8.2e+/-%8.2e mm/mrad\n", max, rms);
}


void param_data_type::corr_eps_y(void)
{
  int  i, j;
  FILE *outf;

  // Clear skew quad setpoints
  set_bnL_design_fam(globval.qt, Quad, 0.0, 0.0);

  // Find coupling vector
  printf("\n");
  printf("Looking for coupling error\n");
  FindCoupVector(VertCouple);

  //Find and print coupling statistics
  printf("\n");
  printf("Before correction\n");
  SkewStat(VertCouple);

  // Coupling Correction
  printf("\n");
  for (i = 1; i <= n_lin; i++) {
    printf("Looking for correction\n");

    //Find Correcting Settings to skew quadrupoles
    for (j = 1; j <= N_BPM; j++)
      b[j] = VDweight*eta_y[j] - VertCouple[j];

    for (j = N_BPM+1; j <= N_COUPLE; j++)
      b[j] = -VertCouple[j];

    dsvbksb(U, w, V, N_COUPLE, N_SKEW, b, SkewStrengthCorr);

    printf("Applying correction\n");
    // Add correction
    for (j = 1; j <= N_SKEW; j++)
      SetdKLpar(globval.qt, j, -Quad, SkewStrengthCorr[j]);

    printf("\n");
    printf("Looking for coupling error\n");
    // Find coupling vector
    FindCoupVector(VertCouple);

    printf("\n");
    printf("After run %d of correction\n", i);
    // Find and print coupling statistics
    SkewStat(VertCouple);
  } // End of coupling Correction

  outf = file_write(eta_y_FileName);
  for (i = 0; i <= globval.Cell_nLoc; i++)
    fprintf(outf, "%4d %7.3f %s %6.3f %10.3e %10.3e\n",
	    i, Cell[i].S, Cell[i].Elem.PName,
	    Cell[i].Nu[Y_], 1e3*Cell[i].Eta[Y_], 1e3*Cell[i].Etap[Y_]);
  fclose(outf);

  FindCoupVector(VertCouple);
}


void param_data_type::get_IDs(void)
{
  int k;

  printf("\n");
  n_ID_Fams = 0;
  for (k = 0; k < globval.Elem_nFam; k++)
    switch (ElemFam[k].ElemF.Pkind) {
    case Wigl:
      printf("found ID family:   %s %12.5e\n",
	     ElemFam[k].ElemF.PName, ElemFam[k].ElemF.W->BoBrhoV[0]);
      n_ID_Fams++; ID_Fams[n_ID_Fams-1] = k + 1;
      break;
    case Insertion:
      printf("found ID family:   %s %12.5e",
	     ElemFam[k].ElemF.PName, ElemFam[k].ElemF.ID->scaling);
      if (ElemFam[k].ElemF.ID->scaling != 0e0) {
	printf("\n");
	n_ID_Fams++; ID_Fams[n_ID_Fams-1] = k + 1;
      } else
	printf("  not included\n");
      break;
    case FieldMap:
      printf("found ID family:   %s %12.5e\n",
	     ElemFam[k].ElemF.PName, ElemFam[k].ElemF.FM->scl);
      n_ID_Fams++; ID_Fams[n_ID_Fams-1] = k + 1;
      break;
    default:
      break;
    }
}


void set_ID_scl(const int Fnum, const double scl);


void param_data_type::set_IDs(const double scl)
{
  int k;

  printf("\n");
  for (k = 0; k < n_ID_Fams; k++) {
    switch (ElemFam[ID_Fams[k]-1].ElemF.Pkind) {
    case Wigl:
      printf("setting ID family: %s %12.5e\n",
	     ElemFam[ID_Fams[k]-1].ElemF.PName,
	     scl*ElemFam[ID_Fams[k]-1].ElemF.W->BoBrhoV[0]);

      set_Wiggler_BoBrho(ID_Fams[k],
			 scl*ElemFam[ID_Fams[k]-1].ElemF.W->BoBrhoV[0]);
      break;
    case Insertion:
      printf("setting ID family: %s %12.5e\n",
	     ElemFam[ID_Fams[k]-1].ElemF.PName, scl);

      set_ID_scl(ID_Fams[k], scl);
      break;
    case FieldMap:
      printf("setting ID family: %s %12.5e\n",
	     ElemFam[ID_Fams[k]-1].ElemF.PName, scl);

      set_ID_scl(ID_Fams[k], scl);
      break;
    default:
      std::cout << "set_IDs: unknown element type" << std::endl;
      exit_(1);
      break;
    }
  }
}


void param_data_type::reset_quads(void)
{
  int k;

  if (N_Fam > N_Fam_max) {
    printf("reset_quads: N_Fam > N_Fam_max: %d (%d)\n", N_Fam, N_Fam_max);
    exit_(0);
  }

  for (k = 0; k < N_Fam; k++) {
    // Note, actual values can differ from the original values
/*    printf("setting quad family: %s %12.5e\n",
	   ElemFam[Q_Fam[k]-1].ElemF.PName,
	   ElemFam[Q_Fam[k]-1].ElemF.M->PBpar[HOMmax+Quad]);

    set_bn_design_fam(Q_Fam[k], Quad,
		       ElemFam[Q_Fam[k]-1].ElemF.M->PBpar[HOMmax+Quad], 0.0);*/

    printf("setting quad family: %s %12.5e\n",
	   ElemFam[Q_Fam[k]-1].ElemF.PName, b2[k]);

    set_bn_design_fam(Q_Fam[k], Quad, b2[k], 0.0);
  }
}


void param_data_type::SVD(const int m, const int n, double **M,
			  double beta_nu[], double b2Ls_[], const double s_cut,
			  const bool first)
{
  int i, j;

  if (trace) {
    printf("\n");
    printf("SVD: first = %1d, m = %1d n = %1d\n", first, m, n);
  }

  if (first) {
    for (i = 1; i <= m; i++)
      for (j = 1; j <= n; j++)
	U1[i][j] = M[i][j];

    dsvdcmp(U1, m, n, w1, V1);

    if (first) {
      printf("\n");
      printf("singular values:\n");
      printf("\n");
    }

    for (i = 1; i <= n; i++) {
      if (first) printf("%11.3e", w1[i]);
      if (w1[i] < s_cut) {
	w1[i] = 0.0;
	if (first) printf(" (zeroed)");
      }
      if (first) if (i % 8 == 0) printf("\n");
    }
    if (first) if (n % 8 != 0) printf("\n");
  }

  dsvbksb(U1, w1, V1, m, n, beta_nu, b2Ls_);
}


void param_data_type::quad_config()
{
  int    i, j;
  double an;

  if (N_Fam > N_Fam_max) {
    printf("quad_config: N_Fam > N_Fam_max: %d (%d)\n", N_Fam, N_Fam_max);
    exit_(0);
  }

  Nquad = 0;
  for (i = 0; i < N_Fam; i++) {
    for (j = 1; j <= GetnKid(Q_Fam[i]); j++) {
      Nquad++;

      if (Nquad > n_b2_max) {
        printf("quad_config: max no of quadrupoles exceeded %d (%d)\n",
               Nquad, n_b2_max);
        exit_(1);
      }

      quad_prms[Nquad-1] = Elem_GetPos(Q_Fam[i], j);

      if (j == 1) get_bn_design_elem(Q_Fam[i], j, Quad, b2[i], an);
    }
  }

  printf("\n");
  printf("quad_config: Nquad = %d\n", Nquad);
}


bool param_data_type::get_SQ(void)
{
  int  j, k;
//  Vector2  alpha3[3], beta3[3], nu3[3], eta3[3], etap3[3];
  FILE *outf = NULL;

  /* Note, IDs are split for evaluation of the driving terms at the center:
       id1  1, 2
       id2  1, 2
       ...                                                                  */

  // Get Twiss params, no dispersion
  Ring_GetTwiss(false, 0e0);

  if (!status.codflag || !globval.stable) return false;

  // Get global tunes
  Nu_X = globval.TotalTune[X_]; Nu_Y = globval.TotalTune[Y_];

  if (trace) {
    printf("\n");
    printf("nu_x = %8.12f, nu_y = %8.12f\n", Nu_X, Nu_Y);

    // Get Twiss params in sext
    printf("\n");
    printf("Lattice functions at sextupoles:\n");

    outf = file_write("latfunS.out");

    fprintf(outf, "s betax nux betay nuy\n");
  }

  Nsext = 0;
  for (k = 0; k < globval.Cell_nLoc; k++) {
    if ((Cell[k].Elem.Pkind == Mpole) && (Cell[k].Elem.M->n_design == Sext)) {
      Nsext++;

      if (Nsext > n_b3_max) {
        printf("get_SQ: max no of sextupoles exceeded %d (%d)\n",
               Nsext, n_b3_max);
        exit_(1);
      }

      Ss[Nsext-1] = Cell[k].S; S_locs[Nsext-1] = k;

      for (j = 0; j <= 1; j++) {
	sb[j][Nsext-1] = Cell[k].Beta[j];
	sNu[j][Nsext-1] = Cell[k].Nu[j] - nu_0[j];
      }

      if (trace) {
	printf("%-8s %7.3f %8.5f %8.5f %8.5f %8.5f\n",
	       Cell[k].Elem.PName, Ss[Nsext-1],
	       sb[X_][Nsext-1], sNu[X_][Nsext-1]-nu_0[X_],
	       sb[Y_][Nsext-1], sNu[Y_][Nsext-1]-nu_0[Y_]);
	fprintf(outf, "%-8s %7.3f %8.5f %8.5f %8.5f %8.5f\n",
		Cell[k].Elem.PName, Ss[Nsext-1],
		sb[X_][Nsext-1], sNu[X_][Nsext-1]-nu_0[X_],
		sb[Y_][Nsext-1], sNu[Y_][Nsext-1]-nu_0[Y_]);
      }
    }
  }

  if (trace) fclose(outf);

  // Number of sexts in the ring
  printf("No of sextupoles = %d\n", Nsext);

  if (trace) {
    // Get Twiss params in quads
    printf("\n");
    printf("Lattice functions at quadrupoles:\n");

    outf = file_write("latfunQ.out");

    fprintf(outf, "s name betax nux betay nuy\n");
  }

  for (k = 0; k < Nquad; k++) {
    Sq[k] = Cell[quad_prms[k]].S;
    for (j = 0; j <= 1; j++) {
      // does not work for machine file (get_twiss_3)...
//       if (Cell[quad_prms[k]].Elem.M->Pthick == thick) {
// 	get_twiss3(quad_prms[k], alpha3, beta3, nu3, eta3, etap3);
// 	qb[j][k] = beta3[Y_][j]; qNu[j][k] = nu3[Y_][j] - nu_0[j];
//       } else {
	qb[j][k] = Cell[quad_prms[k]].Beta[j];
	qNu[j][k] = Cell[quad_prms[k]].Nu[j] - nu_0[j];
//       }
    }

    if (trace) {
      printf("%-8s %7.3f %8.5f %8.5f %8.5f %8.5f\n",
	     Cell[quad_prms[k]].Elem.PName, Sq[k], qb[X_][k],
	     qNu[X_][k], qb[Y_][k], qNu[Y_][k]);

      fprintf(outf, "%-8s %7.3f %8.5f %8.5f %8.5f %8.5f\n",
	      Cell[quad_prms[k]].Elem.PName, Sq[k], qb[X_][k],
	      qNu[X_][k], qb[Y_][k], qNu[Y_][k]);
    }
  }

  if (trace) fclose(outf);

  // Number of quads in the ring
  printf("No of quads      = %d\n", Nquad);

  return true;
}


double param_data_type::Bet(double bq, double nus, double nuq, double NuQ)
{
  return bq*cos(2.0*M_PI*(2.0*fabs(nus-nuq)-NuQ))/(2.0*sin(2.0*M_PI*NuQ));
}


double param_data_type::Nus(double bq, double nus, double nuq, double NuQ)
{
  double Nu, sgn;

  sgn = ((nus-nuq) <= 0)? -1: 1;

  Nu = -bq*sgn*(sin(2.0*M_PI*NuQ)+sin(2.0*M_PI*(2.0*fabs(nus-nuq)-NuQ)))
       /(8.0*M_PI*sin(2.0*M_PI*NuQ));

  return Nu;
}


void param_data_type::A_matrix(void)
{
  int    k, j;
  double BtX, BtY, NuX, NuY;

  // Defining Twiss in undisturbed quads
  for (k = 0; k < Nquad; k++)
    for (j = 0; j <= 1; j++) {
      qb0[j][k] = qb[j][k]; qNu0[j][k] = qNu[j][k];
    }

  // Defining Twiss in undisturbed sexts
  for (k = 0; k < Nsext; k++)
    for (j = 0; j <= 1; j++)
      sNu0[j][k] = sNu[j][k];

  // Now creating matrix A in X=A*B2L
  for (k = 1; k <= Nsext; k++) {
    for (j = 1; j <= Nquad; j++) {
      BtX = Bet(qb0[X_][j-1], sNu0[X_][k-1], qNu0[X_][j-1], Nu_X0);
      NuX = -Nus(qb0[X_][j-1], sNu0[X_][k-1], qNu0[X_][j-1], Nu_X0);
      BtY = -Bet(qb0[Y_][j-1], sNu0[Y_][k-1], qNu0[Y_][j-1], Nu_Y0);
      NuY = Nus(qb0[Y_][j-1], sNu0[Y_][k-1], qNu0[Y_][j-1], Nu_Y0);
      A1[k][j] = scl_dbeta*BtX;
      A1[k+Nsext][j] = scl_dbeta*BtY;
      A1[k+2*Nsext][j] = scl_dnu*NuX;
      A1[k+3*Nsext][j] = scl_dnu*NuY;
    }
  }
  // Now adding 2 more constraints for global tunes
  for (j = 1; j <= Nquad; j++) {
    A1[4*Nsext+1][j] = -scl_nu*qb0[X_][j-1]/(4.0*M_PI);
    A1[4*Nsext+2][j] =  scl_nu*qb0[Y_][j-1]/(4.0*M_PI);
  }

  if (trace) {
    printf("\n");
    printf("AA:\n");
    printf("\n");
    for (k = 1; k <= Nconstr; k++) {
      for (j = 1; j <= Nquad; j++)
	printf(" %10.3e", A1[k][j]);
      printf("\n");
    }
  }
}


void param_data_type::X_vector(const bool first)
{
  int k;

  dnu0[X_] = globval.TotalTune[X_] - Nu_X0;
  dnu0[Y_] = globval.TotalTune[Y_] - Nu_Y0;

  if (first) {
    // Initial fill of X
    for (k = 1; k <= Nsext; k++) {
      Xsext0[k]         = sb[X_][k-1];  Xsext0[k+Nsext]   = sb[Y_][k-1];
      Xsext0[k+2*Nsext] = sNu[X_][k-1]; Xsext0[k+3*Nsext] = sNu[Y_][k-1];
    }
    Xsext0[4*Nsext+1] = 0.0; Xsext0[4*Nsext+2] = 0.0;
  } else {
    // Now substracting from X in X=A*B2L
    for (k = 1; k <= Nsext; k++) {
      Xsext[k]         = scl_dbeta*(Xsext0[k]-sb[X_][k-1])/sb[X_][k-1];
      Xsext[k+Nsext]   = scl_dbeta*(Xsext0[k+Nsext]-sb[Y_][k-1])/sb[Y_][k-1];
      Xsext[k+2*Nsext] = scl_dnu*(Xsext0[k+2*Nsext]-sNu[X_][k-1]+dnu0[X_]/2.0);
      Xsext[k+3*Nsext] = scl_dnu*(Xsext0[k+3*Nsext]-sNu[Y_][k-1]+dnu0[Y_]/2.0);
    }
    Xsext[4*Nsext+1] = scl_nu*(Nu_X0-globval.TotalTune[X_]);
    Xsext[4*Nsext+2] = scl_nu*(Nu_Y0-globval.TotalTune[Y_]);
  }

  if (trace) {
    printf("\n");
    printf("X:\n");
    printf("\n");
    if (first) {
      for (k = 1; k <= Nconstr; k++) {
	printf(" %10.3e", Xsext0[k]);
	if (k % 10 == 0)  printf("\n");
      }
      if (Nconstr % 10 != 0) printf("\n");
    } else {
      for (k = 1; k <= Nconstr; k++) {
	printf(" %10.3e", Xsext[k]);
	if (k % 10 == 0)  printf("\n");
      }
      if (Nconstr % 10 != 0) printf("\n");
    }
  }
}


void param_data_type::ini_ID_corr(const bool IDs)
{
  int k;

  if (IDs) {
    // store ID families
    get_IDs();

    // zero ID's
    set_IDs(0.0);
  }

  // Configuring quads (1 --> C means thin quads located in the middle of 1s)
  quad_config();

  // Configuring quads (1 --> C means thin quads located in the middle of 1s)
  // Read Betas and Nus
  get_SQ(); Nconstr = 4*Nsext + 2;

  // Note, allocated vectors and matrices are deallocated in ID_corr
  Xsext = dvector(1, Nconstr); Xsext0 = dvector(1, Nconstr);
  b2Ls_ = dvector(1, Nquad); A1 = dmatrix(1, Nconstr, 1, Nquad);
  U1 = dmatrix(1, Nconstr, 1, Nquad); w1 = dvector(1, Nquad);
  V1 = dmatrix(1, Nquad, 1, Nquad);

  for (k = 1; k <= Nquad; k++)
    b2Ls_[k] = 0.0;

  // shift zero point to center of ID
//  nu_0[X_] = Cell[id_loc].Nu[X_]; nu_0[Y_] = Cell[id_loc].Nu[Y_];
  nu_0[X_] = 0.0; nu_0[Y_] = 0.0;

  // Defining undisturbed tunes
  Nu_X0 = globval.TotalTune[X_]; Nu_Y0 = globval.TotalTune[Y_];

  // Set-up matrix A in X=A*b2Ls_
  A_matrix();

  // Now fill the X in X=A*b2Ls_
  X_vector(true);
}


void param_data_type::W_diag(void)
{
  double bxf, byf, nxf, nyf, b2Lsum;
  int    k;

  bxf = 0.0; byf = 0.0; nxf = 0.0; nyf = 0.0;
  for (k = 1; k <= Nsext; k++) {
    bxf += sqr(Xsext[k]);
    byf += sqr(Xsext[k+Nsext]);
    nxf += sqr(Xsext[k+2*Nsext]);
    nyf += sqr(Xsext[k+3*Nsext]);
  }

  dnu0[X_] = globval.TotalTune[X_] - Nu_X0;
  dnu0[Y_] = globval.TotalTune[Y_] - Nu_Y0;

  b2Lsum = 0.0;
  for (k = 1; k <= Nquad; k++)
    b2Lsum += sqr(b2Ls_[k]);

  printf("\n");
  printf("Residuals: beta [%%], dnu : \n");
  printf("dbeta_x: %6.2f dbeta_y: %6.2f nu_x: %12.6e nu_y: %12.6e\n",
	 sqrt(bxf/Nsext)*1e2, sqrt(byf/Nsext)*1e2,
	 sqrt(nxf/Nsext), sqrt(nyf/Nsext));
  printf("tune shift: dnu_x = %8.5f, dnu_y = %8.5f\n", dnu0[X_], dnu0[Y_]);
  printf("Sum b2Ls_: %12.6e\n", sqrt(b2Lsum/Nquad));
}


bool param_data_type::ID_corr(const int N_calls, const int N_steps,
			      const bool IDs)
{
  int    i, j, k, Fnum;
  double b2L, a2L, L;
  FILE   *outf;

  printf("\n");
  printf("ID matching begins!\n");

  outf = file_write("ID_corr.out");
  for (i = 1; i <= N_steps; i++) { //This brings ID strength in steps
    if (IDs) set_IDs((double)i/(double)N_steps);

    get_SQ();                               // Read Betas and Nus
    X_vector(false);                        // Fill in dX in dX=A*db2Ls_
    W_diag();                               // Get statistics
    for (j = 1; j <= N_calls; j++) {
      SVD(Nconstr, Nquad, A1, Xsext, b2Ls_, 1e0, j == 1);

      if ((i == N_steps) && (j == N_calls)) fprintf(outf, "#b_2:\n");

      // add quad strengths (db2Ls_)
      for (k = 1; k <= Nquad; k++) {
	set_dbnL_design_elem(Cell[quad_prms[k-1]].Fnum,
			     Cell[quad_prms[k-1]].Knum, Quad,
			     -ID_step*b2Ls_[k], 0.0);

	if ((i == N_steps) && (j == N_calls)) {
	  Fnum = Cell[quad_prms[k-1]].Fnum; L = Cell[quad_prms[k-1]].Elem.PL;
	  get_bnL_design_elem(Fnum, Cell[quad_prms[k-1]].Knum, Quad, b2L, a2L);
	  // ElemFam not defined for flat file.
	  // fprintf(outf, "%10s %6.2f %3d %8.5f\n",
	  // 	  Cell[quad_prms[k-1]].Elem.PName, Cell[quad_prms[k-1]].S, k,
	  // 	  b2L-ElemFam[Fnum-1].ElemF.M->PBpar[HOMmax+Quad]*L);
	  fprintf(outf, "%10s %6.2f %3d %8.5f\n",
		  Cell[quad_prms[k-1]].Elem.PName, Cell[quad_prms[k-1]].S, k,
		  b2L);
	}
      }

      printf("\n");
      printf("Iteration: %2d\n", j);
      if (get_SQ()) {
	X_vector(false);                    // Fill in dX in dX=A*db2Ls_
	W_diag();                           // Get statistics

	printglob();
      } else {
	printf("ID_corr: correction failed\n");
	// restore lattice
	if (IDs) set_IDs(0.0);
	reset_quads();
	return false;
      }
    }
  }
  fclose(outf);

  outf = file_write("ID_corr_res.out");
  fprintf(outf, "# dbeta_x/beta_x  dbeta_y/beta_y  dnu_x  dnu_y\n");
  fprintf(outf, "#      [%%]             [%%]\n");
  fprintf(outf, "#\n");
  for (k = 1; k <= Nsext; k++)
    fprintf(outf, "%6.1f %6.2f %6.2f %10.3e %10.3e\n",
	    Ss[k-1], 1e2*Xsext[k], 1e2*Xsext[k+Nsext],
	    Xsext[k+2*Nsext], Xsext[k+3*Nsext]);
  fclose(outf);

  // Allow for repeated calls to ID_corr, allocation is done in ini_ID_corr.
  if (false) {
    free_dvector(Xsext, 1, Nconstr); free_dvector(Xsext0, 1, Nconstr);
    free_dvector(b2Ls_, 1, Nquad); free_dmatrix(A1, 1, Nconstr, 1, Nquad);
    free_dmatrix(U1, 1, Nconstr, 1, Nquad); free_dvector(w1, 1, Nquad);
    free_dmatrix(V1, 1, Nquad, 1, Nquad);
  }

  printf("\n");
  printf("ID matching ends!\n");

  return true;
}


char* get_prm(char **p)
{
  char *prm;

  prm = strtok_r(NULL, " \t\r", p);
  if (prm == NULL) {
    std::cout << std::endl;
    std::cout << "get_prm: incorrect format" << std::endl;
    exit_(1);
  }

  return prm;
}


void param_data_type::LoadAlignTol(const bool Scale_it, const double Scale,
				   const bool new_rnd, const int seed) const
{
  char     line[max_str], Name[max_str],  type[max_str];
  int      j, k, Fnum, seed_val;
  long int loc;
  double   dx, dy, dr;  // x and y misalignments [m] and roll error [rad]
  double   dr_deg;
  bool     rms = false, set_rnd;
  FILE     *fp;

  const bool prt = true;

  if (prt) printf("\nreading in %s\n", ae_file.c_str());

  fp = file_read(ae_file.c_str());

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
		misalign_rms_elem(Cell[loc].Fnum, Cell[loc].Knum,
				  dx, dy, dr_deg, new_rnd);
	      else
		misalign_sys_elem(Cell[loc].Fnum, Cell[loc].Knum,
				  dx, dy, dr_deg);
	    }
	} else {
	  Fnum = ElemIndex(Name);
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


void param_data_type::LoadFieldErr(const bool Scale_it, const double Scale,
				   const bool new_rnd) const
{
  bool          rms, set_rnd;
  char          line[max_str], name[max_str], type[max_str], *prm, *p;
  int           k, n, seed_val;
  double        Bn, An, r0;
  std::ifstream inf;

  file_rd(inf, fe_file.c_str());

  set_rnd = false;
  std::cout << std::endl;
  while (!inf.getline(line, max_str).eof()) {
    if (strstr(line, "#") == NULL) {
      // New seed?
      sscanf(line, "%s", name);
      if (strcmp("seed", name) == 0) {
	set_rnd = true;
	sscanf(line, "%*s %d", &seed_val);
	std::cout << "setting random seed to " << seed_val << std::endl;
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
	prm = strtok_r(line, " \t", &p);
	for (k = 1; k <= 2; k++)
	  prm = strtok_r(NULL, " \t", &p);
	while (((prm = strtok_r(NULL, " \t", &p)) != NULL) &&
	       (strcmp(prm, "\r\n") != 0)) {
	  sscanf(prm, "%d", &n);
	  prm = get_prm(&p); sscanf(prm, "%lf", &Bn);
	  prm = get_prm(&p); sscanf(prm, "%lf", &An);
	  if (Scale_it) {
	    Bn *= Scale; An *= Scale;
	  }
	  printf(" %2d %9.1e %9.1e\n", n, Bn, An);
	  // convert to normalized multipole components
	  SetFieldErrors(name, rms, r0, n, Bn, An, true);
	}
      }
    } else
    std::cout << line << std::endl;
  }

  inf.close();
}


void param_data_type::LoadApers(const double scl_x, const double scl_y) const
{
  char   line[max_str], Name[max_str];
  int    Fnum;
  double dxmin, dxmax, dymin, dymax;  // min and max x and apertures
  FILE   *fp;

  bool prt = true;

  fp = file_read(ap_file.c_str());

  printf("\n");
  printf("...Load and Set Apertures.\n");

  while (fgets(line, max_str, fp) != NULL) {
    if (strstr(line, "#") == NULL) {
      sscanf(line,"%s %lf %lf %lf %lf",
	     Name, &dxmin, &dxmax, &dymin, &dymax);
      dxmin *= scl_x; dxmax *= scl_x; dymin *= scl_y; dymax *= scl_y;
      if (strcmp("all", Name)==0) {
	if(prt)
	  printf("setting all apertures to"
		 " dxmin = %e, dxmax = %e, dymin = %e, dymax = %e\n",
		 dxmin, dxmax, dymin, dymax);
	set_aper_type(All, dxmin, dxmax, dymin, dymax);
	//	ini_aper(dxmin, dxmax, dymin, dymax);
      } else if (strcmp("quad", Name)==0) {
	if(prt)
	  printf("setting apertures at all quads to"
		 " dxmin = %e, dxmax = %e, dymin = %e, dymax = %e\n",
		 dxmin, dxmax, dymin, dymax);
	set_aper_type(Quad, dxmin, dxmax, dymin, dymax);
      } else if (strcmp("sext", Name) == 0) {
	if(prt)
	  printf("setting apertures at all sextupoles to"
		 " dxmin = %e, dxmax = %e, dymin = %e, dymax = %e\n",
		 dxmin, dxmax, dymin, dymax);
	set_aper_type(Sext, dxmin, dxmax, dymin, dymax);
      } else {
	Fnum = ElemIndex(Name);
	if(Fnum > 0) {
	  if(prt)
	    printf("setting apertures at all %s to"
		   " dxmin = %e, dxmax = %e, dymin = %e, dymax = %e\n",
		   Name, dxmin, dxmax, dymin, dymax);
	  set_aper_fam(Fnum, dxmin, dxmax, dymin, dymax);
	} else
	  printf("LoadApers: lattice does not contain element %s\n", Name);
      }
    } else
      printf("%s", line);
  }

  fclose(fp);
}


void param_data_type::Align_BPMs(const int n) const
{
  // Align BPMs to adjacent multipoles.

  bool     aligned;
  int      i, j, k;
  long int loc;

  const int n_step = 5;

  printf("\n");
  for (i = 1; i <= GetnKid(globval.bpm); i++) {
    loc = Elem_GetPos(globval.bpm, i);

    if ((loc == 1) || (loc == globval.Cell_nLoc)) {
      printf("Align_BPMs: BPM at entrance or exit of lattice: %ld\n", loc);
      exit_(1);
    }

    j = 1; aligned = false;
    do {
      if ((Cell[loc-j].Elem.Pkind == Mpole) &&
	  (Cell[loc-j].Elem.M->n_design == n)) {
	for (k = 0; k <= 1; k++)
	  Cell[loc].Elem.M->PdSsys[k] = Cell[loc-j].dS[k];
	printf("aligned BPM no %1d to %s\n", i, Cell[loc-j].Elem.PName);
	aligned = true; break;
      } else if ((Cell[loc+j].Elem.Pkind == Mpole) &&
		 (Cell[loc+j].Elem.M->n_design == n)) {
	for (k = 0; k <= 1; k++)
	  Cell[loc].Elem.M->PdSsys[k] = Cell[loc+j].dS[k];
	printf("aligned BPM no %1d to %s\n", i, Cell[loc+j].Elem.PName);
	aligned = true; break;
      }

      j++;
    } while (j <= n_step);

    if (aligned)
      Mpole_SetdS(globval.bpm, i);
    else
      printf("Align_BPMs: no multipole adjacent to BPM no %d\n", i);
  }
}


bool param_data_type::CorrectCOD_N(const int n_orbit, const int k)
{
  bool     cod = false;
  int      i, j;
  long int loc;
  double   m_dbeta[2], s_dbeta[2], m_dnu[2], s_dnu[2];

  // Clear trim setpoints
  for (j = 0; j < 2; j++)
    for (i = 1; i <= n_corr_[j]; i++) {
      loc = corrs_[j][i];
      set_bnL_design_elem(Cell[loc].Fnum, Cell[loc].Knum, Dip, 0.0, 0.0);
    }

  // load misalignments
  LoadAlignTol(true, 1.0, true, k);
  for (i = 1; i <= n_scale; i++) {
    // Scale the rms values
    LoadAlignTol(true, (double)i/(double)n_scale, false, k);

    if (bba) {
      // Beam based alignment
      Align_BPMs(Quad);
    }

    // get_traject();
    
    cod = CorrectCOD(n_orbit, 1e0);

    if (!cod) break;

    get_dbeta_dnu(m_dbeta, s_dbeta, m_dnu, s_dnu);
    printf("\n");
    printf("RMS dbeta_x/beta_x = %4.2f%%,   dbeta_y/beta_y = %4.2f%%\n",
	   1e2*s_dbeta[X_], 1e2*s_dbeta[Y_]);
    printf("RMS dnu_x          = %7.5f, dnu_y          = %7.5f\n",
	   s_dnu[X_], s_dnu[Y_]);
  }

  return cod;
}


void param_data_type::ini_COD_corr(const int n_bpm_Fam,
				   const std::string bpm_names[],
				   const int n_hcorr_Fam,
				   const std::string hcorr_names[],
				   const int n_vcorr_Fam,
				   const std::string vcorr_names[],
				   const bool svd)
{
  int i, j, Fnum, n_bpm, n_hcorr, n_vcorr;

  n_bpm = 0;
  for (i = 0; i < n_bpm_Fam; i++)
    n_bpm += GetnKid(ElemIndex(bpm_names[i]));

  n_hcorr = 0;
  for (i = 0; i < n_hcorr_Fam; i++)
    n_hcorr += GetnKid(ElemIndex(hcorr_names[i]));

  n_vcorr = 0;
  for (i = 0; i < n_vcorr_Fam; i++)
    n_vcorr += GetnKid(ElemIndex(vcorr_names[i]));

  long int  bpms[n_bpm], hcorrs[n_hcorr], vcorrs[n_vcorr];

  n_bpm = 0;
  for (i = 0; i < n_bpm_Fam; i++) {
    Fnum = ElemIndex(bpm_names[i]);
    for (j = 1; j <= GetnKid(Fnum); j++)
      bpms[n_bpm++] = Elem_GetPos(Fnum, j);
  }

  n_hcorr = 0;
  for (i = 0; i < n_hcorr_Fam; i++) {
    Fnum = ElemIndex(hcorr_names[i]);
    for (j = 1; j <= GetnKid(Fnum); j++)
      hcorrs[n_hcorr++] = Elem_GetPos(Fnum, j);
  }

  n_vcorr = 0;
  for (i = 0; i < n_vcorr_Fam; i++) {
    Fnum = ElemIndex(vcorr_names[i]);
    for (j = 1; j <= GetnKid(Fnum); j++)
      vcorrs[n_vcorr++] = Elem_GetPos(Fnum, j);
  }

  std::cout << std::endl;
  std::cout << "ini_COD_corr: n_bpm = " << n_bpm << ", n_hcorr = " << n_hcorr
       << ", n_vcorr = " << n_vcorr << std::endl;

  gcmat(n_bpm, bpms, n_hcorr, hcorrs, 1, svd);
  gcmat(n_bpm, bpms, n_vcorr, vcorrs, 2, svd);

  if (true) {
    gtcmat(n_bpm, bpms, n_hcorr, hcorrs, 1, svd);
    gtcmat(n_bpm, bpms, n_vcorr, vcorrs, 2, svd);
  }
}


bool param_data_type::cod_corr(const int n_cell, const double scl,
			       orb_corr_type orb_corr[])
{
  bool            cod = false;
  long int        lastpos;
  double          m_dbeta[2], s_dbeta[2], m_dnu[2], s_dnu[2];
  ss_vect<double> ps;

  orb_corr[X_].clr_trims(); orb_corr[Y_].clr_trims();
  
  cod = getcod(0e0, lastpos);

  if (!cod) {
    printf("\ncould not find closed orbit; threading beam\n");
    orb_corr[X_].clr_trims(); orb_corr[Y_].clr_trims();
    thread_beam(n_cell, loc_Fam_name, bpm_Fam_names, corr_Fam_names,
		n_orbit, scl);

    // prt_cod("cod.out", globval.bpm, true);    
  }

  cod = cod_correct(n_orbit, scl, orb_corr);

  get_dbeta_dnu(m_dbeta, s_dbeta, m_dnu, s_dnu);
  printf("\ncod_corr: rms dbeta_x/beta_x = %4.2f%%"
	 ",   dbeta_y/beta_y = %4.2f%%\n",
	 1e2*s_dbeta[X_], 1e2*s_dbeta[Y_]);
  printf("          rms dnu_x          = %7.5f, dnu_y          = %7.5f\n",
	 s_dnu[X_], s_dnu[Y_]);

  prt_cod("cod.out", globval.bpm, true);    

  return cod;
}


void param_data_type::Orb_and_Trim_Stat(void)
{
  int     i, j, N;
  int     SextCounter = 0;
  int     bins[5] = { 0, 0, 0, 0, 0 };
  double  bin = 40.0e-6; // bin size for stat
  double  tr; // trim strength
  Vector2 Sext_max, Sext_sigma, TrimMax, orb;

  Sext_max[X_] = 0.0; Sext_max[Y_] = 0.0;
  Sext_sigma[X_] = 0.0; Sext_sigma[Y_] = 0.0;
  TrimMax[X_] = 0.0; TrimMax[Y_] = 0.0;
  N = globval.Cell_nLoc; SextCounter = 0;
  for (i = 0; i <= N; i++) {
    if ((Cell[i].Elem.Pkind == Mpole) && (Cell[i].Elem.M->n_design == Sext)) {
      SextCounter++;
      orb[X_] = Cell[i].BeamPos[x_]; orb[Y_] = Cell[i].BeamPos[y_];
      Sext_sigma[X_] += orb[X_]*orb[X_]; Sext_sigma[Y_] += orb[Y_]*orb[Y_];
      if (fabs(orb[X_]) > Sext_max[X_]) Sext_max[X_] = fabs(orb[X_]);
      if (fabs(orb[Y_]) > Sext_max[Y_]) Sext_max[Y_] = fabs(orb[Y_]);
      j = (int) (sqrt(orb[X_]*orb[X_]+orb[Y_]*orb[Y_])/bin);
      if (j > 4) j = 4;
      if (j >= 0)
	bins[j]++;
      else
	printf("\nOrb_and_Trim_Stat: negative bin %d\n", j);
    } // sextupole handling

    if (Cell[i].Fnum == globval.hcorr) {
      tr = Cell[i].Elem.M->PBpar[HOMmax+Dip];
      if (fabs(tr) > TrimMax[X_]) TrimMax[X_] = fabs(tr);
    } // horizontal trim handling
    if (Cell[i].Fnum == globval.vcorr) {
      tr = Cell[i].Elem.M->PBpar[HOMmax-Dip];
      if (fabs(tr) > TrimMax[Y_]) TrimMax[Y_] = fabs(tr);
    } // vertical trim handling
  } // looking throught the cells

  Sext_sigma[X_] = sqrt(Sext_sigma[X_]/SextCounter);
  Sext_sigma[Y_] = sqrt(Sext_sigma[Y_]/SextCounter);

  printf("In sextupoles maximal horizontal orbit is:"
	 " %5.3f mm with sigma %5.3f mm\n",
          1e3*Sext_max[X_], 1e3*Sext_sigma[X_]);
  printf("and maximal vertical orbit is:            "
	 " %5.3f mm with sigma %5.3f mm.\n",
	 1e3*Sext_max[Y_], 1e3*Sext_sigma[Y_]);

  for (i = 0; i < 4;  i++) {
    printf("There are %3d sextupoles with offset between ", bins[i]);
    printf(" %5.3f mm and %5.3f mm\n", i*bin*1e3, (i+1)*bin*1e3);
  }
  printf("There are %3d sextupoles with offset ", bins[4]);
  printf("more than %5.3f mm \n", 4e3*bin);
  printf("Maximal hcorr is %6.3f mrad, maximal vcorr is %6.3f mrad\n",
	 1e3*TrimMax[X_], 1e3*TrimMax[Y_]);
}


void param_data_type::prt_cod_corr_lat(void)
{
  int  i;
  FILE *CodCorLatFile;

  CodCorLatFile = file_write(CodCorLatFileName);

  fprintf(CodCorLatFile, "#    name     s   sqrt(BxBy) betaX    nuX"
	  "    betaY    nuY     etaX etaX*betaY nuX-nuY \n");
  fprintf(CodCorLatFile, "#            [m]     [m]      [m]             [m]"
	  "              [m]     [m*m] \n");

  for (i = 0; i <= globval.Cell_nLoc; i++){
    fprintf(CodCorLatFile, "%4d:", i);

    if (i == 0)
      fprintf(CodCorLatFile, "%.*s", 6, "begin ");
    else
      fprintf(CodCorLatFile, "%.*s", 6, Cell[i].Elem.PName);

    fprintf(CodCorLatFile, "%7.3f  %5.2f    %5.2f  %7.4f  %5.2f  %7.4f"
	    "  %6.3f  %6.3f  %6.3f\n",
	    Cell[i].S, sqrt(Cell[i].Beta[X_]*Cell[i].Beta[Y_]),
            Cell[i].Beta[X_], Cell[i].Nu[X_], Cell[i].Beta[Y_], Cell[i].Nu[Y_],
	    Cell[i].Eta[X_], Cell[i].Eta[X_]*Cell[i].Beta[Y_],
            Cell[i].Nu[X_]-Cell[i].Nu[Y_]);
  }
  fclose(CodCorLatFile);
}


void param_data_type::err_and_corr_init(const string &param_file,
					orb_corr_type orb_corr[])
{
  globval.Cavity_on   = false; globval.radiation = false;
  globval.Aperture_on = false;

  get_param(param_file);

  Ring_GetTwiss(true, 0.0); printglob();

  get_bare();

  if (ae_file != "") {
    if (bba) {
      Align_BPMs(Quad);
    }

    cod_ini(bpm_Fam_names, corr_Fam_names, orb_corr);
  }

  if (N_calls > 0) ini_ID_corr(false);

  if (n_lin > 0) ini_skew_cor(disp_wave_y);
}


void param_data_type::err_and_corr_exit(orb_corr_type orb_corr[])
{
  int j;

  if (ae_file != "") {
    for (j = 0; j < 2; j++)
      orb_corr[j].dealloc();
  }
}


void get_bn2(const string file_name1, const string file_name2, int n,
	     const bool prt)
{
  char   line[max_str], str[max_str], str1[max_str], *token, *name, *p;
  int    n_prm, Fnum, Knum, order;
  double bnL, bn, C, L;
  FILE   *inf, *fp_lat;

  inf = file_read(file_name1.c_str()); fp_lat = file_write(file_name2.c_str());

  // if n = 0: go to last data set
  if (n == 0) {
    while (fgets(line, max_str, inf) != NULL )
      if (strstr(line, "n = ") != NULL)	sscanf(line, "n = %d", &n);

    fclose(inf); inf = file_read(file_name1.c_str());
  }

  if (prt) {
    printf("\n");
    printf("reading values (n=%d): %s\n", n, file_name1.c_str());
    printf("\n");
  }

  sprintf(str, "n = %d", n);
  do
    fgets(line, max_str, inf);
  while (strstr(line, str) == NULL);

  fprintf(fp_lat, "\n");
  n_prm = 0;
  while (fgets(line, max_str, inf) != NULL) {
    if (strcmp(line, "\n") == 0) break;
    n_prm++;
    name = strtok_r(line, "(", &p);
    rm_space(name);
    strcpy(str, name); Fnum = ElemIndex(str);
    strcpy(str1, name); upr_case(str1);
    token = strtok_r(NULL, ")", &p); sscanf(token, "%d", &Knum);
    strtok_r(NULL, "=", &p); token = strtok_r(NULL, "\n", &p);
    sscanf(token, "%lf %d", &bnL, &order);
    if (prt) printf("%6s(%2d) = %10.6f %d\n", name, Knum, bnL, order);

    if (Fnum != 0) {
      if (order == 0)
        SetL(Fnum, bnL);
      else
        SetbnL(Fnum, order, bnL);

      L = GetL(Fnum, 1);
      if (Knum == 1) {
	if (order == 0)
	  fprintf(fp_lat, "%s: Drift, L = %8.6f;\n", str1, bnL);
	else {
	  bn = (L != 0.0)? bnL/L : bnL;
	  if (order == Quad)
	    fprintf(fp_lat, "%s: Quadrupole, L = %8.6f, K = %19.16f,"
		    " N = Nquad, Method = Meth;\n", str1, L, bn);
	  else if (order == Sext)
	    fprintf(fp_lat, "%s: Sextupole, L = %8.6f, K = %19.16f"
		    ", N = Nsext, Method = Meth;\n", str1, L, bn);
	  else {
	    fprintf(fp_lat, "%s: Multipole, L = %8.6f"
		    ", N = 1, Method = Meth,\n", str1, L);
	    fprintf(fp_lat, "     HOM = (%d, %19.16f, %3.1f);\n",
		    order, bn, 0.0);
	  }
	}
      }
    } else {
      printf("element %s not found\n", name);
      exit_(1);
    }
  }
  if (prt) printf("\n");

  C = Cell[globval.Cell_nLoc].S; recalc_S();
  if (prt)
    printf("New Cell Length: %5.3f (%5.3f)\n", Cell[globval.Cell_nLoc].S, C);

  fclose(inf); fclose(fp_lat);
}
