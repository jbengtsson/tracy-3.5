// TODO: naming conventions and refactor - functionally, to set engeniring tolerances and correct including coping with IDs (id_corr is NOT LOCO)

// Define static variables.

bool        param_data_type::DA_bare      = false;
bool        param_data_type::freq_map     = false;
int         param_data_type::n_orbit      = 5;
int         param_data_type::n_scale      = 1;
std::string param_data_type::loc_Fam_name = "";
int         param_data_type::n_cell       = -1;
int         param_data_type::n_thread     = -1;

int param_data_type::n_lin         =  3;
int param_data_type::SQ_per_scell  =  1;
int param_data_type::BPM_per_scell = 10;
int param_data_type::HCM_per_scell = 10;
int param_data_type::VCM_per_scell = 10;

double param_data_type::kick       = 0.01e-3;
double param_data_type::v_maxkick  = 1.0e-3;
double param_data_type::h_maxkick  = 1.0e-3;
double param_data_type::h_cut      = 1.0e-4;
double param_data_type::v_cut      = 1.0e-4;
int    param_data_type::n_stat     = 1;
int    param_data_type::n_meth     = 0;

std::vector<double> bn_an[HOMmax+HOMmax+1];

double param_data_type::ID_s_cut    = 1e1;
  
double param_data_type::VDweight    = 1e3,
       param_data_type::HVweight    = 1e0,
       param_data_type::VHweight    = 1e0,
       param_data_type::qt_s_cut    = 1e0,
       param_data_type::disp_wave_y = 0e0,
       param_data_type::disp_wave_o = 0e0;
int    param_data_type::qt_from_file = 0;

double param_data_type::TuneX       = 0e0,
       param_data_type::TuneY       = 0e0,
       param_data_type::ChromX      = 1e6,
       param_data_type::ChromY      = 1e6;

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

//>>>> string copy functions
void TracyStrcpy(char *elem, char *pname) {
  long i;
  strncpy(elem, pname, SymbolLength);
  elem[SymbolLength] = '\0';

  i = SymbolLength - 1; // remove trailing spaces
  while (i >= 0 && elem[i] == ' ') {
    elem[i] = '\0';
    i--;
  }
}

void MyStrcpy (char *elem, char *pname, long leng) {
  long i;

  strncpy(elem, pname, leng); elem[leng]='\0';
  i = leng-1; // remove trailing spaces
  while ( elem[i] == ' ' ) {
     elem[i] = '\0';
     i--;
  }
}

#define seps 1E-6

void param_data_type::GirderSetup() {
  // Extracted to correction/girder_model; kept as a delegator during the refactor.
  girders.GirderSetup();
}

void param_data_type::SetCorMis(double gxrms, double gyrms, double gtrms,
				double jxrms, double jyrms, double exrms,
				double eyrms, double etrms, double rancutx,
				double rancuty, double rancutt, long iseed)
{
  // Extracted to correction/girder_model; kept as a delegator during the refactor.
  girders.SetCorMis(gxrms, gyrms, gtrms, jxrms, jyrms, exrms, eyrms, etrms,
		    rancutx, rancuty, rancutt, iseed);
}

void param_data_type::CorMis_in(double *gdxrms, double *gdzrms, double *gdarms, double *jdxrms, double *jdzrms, double *edxrms, double *edzrms, double *edarms, double *bdxrms, double *bdzrms, double *bdarms, double *rancutx, double *rancuty, double *rancutt, long *iseed, long *iseednr)
{
  // Extracted to correction/girder_model; kept as a delegator during the refactor.
  corr::CorMis_in(gdxrms, gdzrms, gdarms, jdxrms, jdzrms, edxrms, edzrms, edarms,
		  bdxrms, bdzrms, bdarms, rancutx, rancuty, rancutt, iseed,
		  iseednr);
}

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
      } else if (strcmp("at_flat_file", name) == 0) {
	sscanf(line, "%*s %s", str);
	sstr.clear(); sstr.str("");
	sstr << str << "flat_file.dat"; flat_file = sstr.str();
	rdmfile_at(flat_file.c_str());
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
        else if (strcmp("n_meth", name) == 0)
	sscanf(line, "%*s %d", &n_meth);
        else if (strcmp("h_maxkick", name) == 0)
	sscanf(line, "%*s %lf", &h_maxkick);
        else if (strcmp("v_maxkick", name) == 0)
	sscanf(line, "%*s %lf", &v_maxkick);
        else if (strcmp("h_cut", name) == 0)
	sscanf(line, "%*s %lf", &h_cut);
        else if (strcmp("v_cut", name) == 0)
	sscanf(line, "%*s %lf", &v_cut);
      else if (strcmp("n_aper", name) == 0)
	sscanf(line, "%*s %d", &n_aper_DA);
      else if (strcmp("loc_name", name) == 0) {
        sscanf(line, "%*s %s", str);
	sstr.clear(); sstr.str(""); sstr << str; loc_Fam_name = sstr.str();
      } else if (strcmp("n_cell", name) == 0)
	sscanf(line, "%*s %d", &n_cell);
      else if (strcmp("n_thread", name) == 0)
	sscanf(line, "%*s %d", &n_thread);
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
      } else if (strcmp("bpm", name) == 0) {
	sscanf(line, "%*s %s", str);
	globval.bpm = ElemIndex(str);
      } else if (strcmp("hcorr", name) == 0) {
	sscanf(line, "%*s %s", str);
	globval.hcorr = ElemIndex(str);
      } else if (strcmp("vcorr", name) == 0) {
	sscanf(line, "%*s %s", str);
	globval.vcorr = ElemIndex(str);
      } else if (strcmp("qt", name) == 0) {
	sscanf(line, "%*s %s", str);
	globval.qt = ElemIndex(str);
      } else if (strcmp("nux", name) == 0)
	sscanf(line, "%*s %le", &TuneX);
      else if (strcmp("nuy", name) == 0)
	sscanf(line, "%*s %le", &TuneY);
      else if (strcmp("six", name) == 0)
	sscanf(line, "%*s %le", &ChromX);
      else if (strcmp("siy", name) == 0)
	sscanf(line, "%*s %le", &ChromY);
      else if (strcmp("qt_s_cut", name) == 0)
	sscanf(line, "%*s %le", &qt_s_cut);
      else if (strcmp("disp_wave_y", name) == 0)
	sscanf(line, "%*s %lf", &disp_wave_y);
      else if (strcmp("disp_wave_o", name) == 0)
	sscanf(line, "%*s %lf", &disp_wave_o);
      else if (strcmp("qt_from_file", name) == 0)
	sscanf(line, "%*s %d", &qt_from_file);
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
	if (trace) {
	  printf("\nID_quads:\n");
	  for (auto k = 0; k < N_Fam; k++)
	    printf("  %10d\n", Q_Fam[k]);
	}
      } else if (strcmp("ID_s_cut", name) == 0)
	sscanf(line, "%*s %le", &ID_s_cut);
      else {
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

  nu0_[X_] = globval.TotalTune[X_];
  nu0_[Y_] = globval.TotalTune[Y_];
}


void param_data_type::get_dbeta_dnu(double m_dbeta[], double s_dbeta[],
				    double m_dnu[], double s_dnu[])
{
  // Extracted to correction/corr_utils; kept as a delegator during the refactor.
  corr::get_dbeta_dnu(m_dbeta, s_dbeta, m_dnu, s_dnu, n_sext, sexts, betas0_,
		      nus0_);
}


void param_data_type::FindSQ_SVDmat(double **SkewRespMat, double **U,
				    double **V, double *w, int N_COUPLE,
				    int N_SKEW)
{
  int i, j;

  for (i = 1; i <= N_COUPLE; i++)
    for (j = 1; j <= N_SKEW; j++)
      U[i][j] = SkewRespMat[i][j];

  // prepare matrices for SVD
  dsvdcmp(U, N_COUPLE, N_SKEW, w, V);

  printf("\n");
  printf("singular values: s_cut = %10.3e\n", qt_s_cut);
  printf("\n");

  // zero singular values
  printf("\n");
  printf("singular values:\n");
  printf("\n");
  for (i = 1; i <= N_SKEW; i++) {
    printf("%11.3e", w[i]);
    if (w[i] < qt_s_cut) {
      w[i] = 0.0;
      printf(" (zeroed)");
      if (i % 8 == 0) printf("\n");
    }
  }
  if (i % 8 != 0) printf("\n");
}

// Read eta values from the file
void param_data_type::ReadEta(const char *TolFileName) 
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
void param_data_type::FindMatrix(double **SkewRespMat, const double deta_y_max,
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
      SkewRespMat[j][i] = VDweight*0.5*alpha*sqrt(betaSQ[i][Yi]*betaBPM[j][Yi])
	*cos(twopi*fabs(nuSQ[i][Yi]-nuBPM[j][Yi])-pi*nuY)/sin(pi*nuY);
    } // for (j=1; j<=N_BPM; j++)

    // looking for coupling of horizontal trim to vertical BPM
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

  if (deta_y_max < 0.) {
    ReadEta("eta_file.dat");
  }
  eta_y_max = -1e8;
  eta_y_min =  1e8;
  for (j = 1; j <= N_BPM; j++) {
    if (deta_y_max > 0.) {
      eta_y[j] = 0.0;
      for (i = 1; i <= N_SKEW; i++)
	if (i % SQ_per_scell == 0) {
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
	  N_BPM, SQ_per_scell, 1e3*eta_y_min, 1e3*eta_y_max, 1e3*deta_y_max,
	  1e3*deta_y_max*deta_y_offset);
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
} // FindMatrix


void param_data_type::ini_skew_cor(const double deta_y_max,
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
  FindMatrix(SkewRespMat, deta_y_max, deta_y_offset);

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


void param_data_type::SkewStat(double VertCouple[], const int cnt)
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
    if (fabs(VertCouple[i]/VDweight) > max) max = fabs(VertCouple[i]/VDweight);
    rms += sqr(VertCouple[i]/VDweight);
    mean += VertCouple[i]/VDweight;
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
    if (fabs(VertCouple[i]/HVweight) > max) max = fabs(VertCouple[i]/HVweight);
    rms += sqr(VertCouple[i]/HVweight);
    mean += VertCouple[i]/HVweight;
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
    if (fabs(VertCouple[i]/VHweight) > max) max = fabs(VertCouple[i]/VHweight);
    rms += sqr(VertCouple[i]/VHweight);
    mean += VertCouple[i]/VHweight;
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


void param_data_type::corr_eps_y(const int cnt)
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
  FindCoupVector(VertCouple);

  //Find and print coupling statistics
  printf("\n");
  printf("Before correction\n");
  SkewStat(VertCouple, -1);

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
    SkewStat(VertCouple, -1);

  } // End of coupling Correction

  if (qt_from_file) {
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
    FindCoupVector(VertCouple);
    printf("\n");
    printf("After application of skew values from file 'qt_file.dat'\n");
    SkewStat(VertCouple, -1);
  }
      
  SkewStat(VertCouple, cnt);

  snprintf(fname, sizeof(fname), "%s_%d.out",eta_y_FileName,cnt);
  outf = file_write(fname);

  fprintf(outf, "# nr s [m] name nuy etay [mm] etapy [mrad]\n");
  for (i = 0; i <= globval.Cell_nLoc; i++)
    fprintf(outf, "%4d %7.3f %s %6.3f %10.3e %10.3e\n",
	    i, Cell[i].S, Cell[i].Elem.PName,
	    Cell[i].Nu[Y_], 1e3*Cell[i].Eta[Y_], 1e3*Cell[i].Etap[Y_]);
  fclose(outf);

  FindCoupVector(VertCouple);
}


void param_data_type::reset_quads(void)
{
  // Extracted to correction/id_corr; kept as a delegator during the refactor.
  id.reset_quads(N_Fam, Q_Fam);
}


// Initializing ID correction (NOT LOCO).
void param_data_type::ini_ID_corr(const bool IDs)
{
  // Extracted to correction/id_corr; kept as a delegator during the refactor.
  id.ini_ID_corr(IDs, N_Fam, Q_Fam);
}


bool param_data_type::ID_corr(const int N_calls, const int N_steps,
			      const bool IDs, const int cnt)
{
  // Extracted to correction/id_corr; kept as a delegator during the refactor.
  return id.ID_corr(N_calls, N_steps, IDs, cnt, N_Fam, Q_Fam, ID_s_cut);
}


void param_data_type::ReadCorMis(const bool Scale_it, const double Scale) const
{
  // Extracted to correction/error_model; kept as a delegator during the refactor.
  corr::ReadCorMis(Scale_it, Scale);
}

void param_data_type::LoadAlignTol(const bool Scale_it, const double Scale,
				   const bool new_rnd, const int seed) const
{
  // Extracted to correction/error_model; kept as a delegator during the refactor.
  corr::LoadAlignTol(ae_file, Scale_it, Scale, new_rnd, seed);
}


void param_data_type::LoadFieldErr(const bool Scale_it, const double Scale,
				   const bool new_rnd) const
{
  // Extracted to correction/error_model; kept as a delegator during the refactor.
  corr::LoadFieldErr(fe_file, Scale_it, Scale, new_rnd);
}


void param_data_type::LoadApers(const double scl_x, const double scl_y) const
{
  // Extracted to correction/error_model; kept as a delegator during the refactor.
  corr::LoadApers(ap_file, scl_x, scl_y);
}


void param_data_type::Align_BPMs(const int n, const double bdxrms,
				 const double bdzrms, const double bdarms) const
{
  // Align BPMs to adjacent multipoles.

  bool     aligned;
  int      i, j, k;
  long int loc;

  const int n_step = 25;

  // printf("Align_BPMs entered %d\n", GetnKid(globval.bpm));
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
	if (bdxrms >=0.) {
	  Cell[loc].Elem.M->PdSrms[0] = bdxrms;
	  Cell[loc].Elem.M->PdSrnd[0] = normranf();
	} 
	if (bdzrms >=0.) {
	  Cell[loc].Elem.M->PdSrms[1] = bdzrms;
	  Cell[loc].Elem.M->PdSrnd[1] = normranf();
	}
	if (bdarms >=0.) {
	  Cell[loc].Elem.M->PdTrms = bdarms;
	  Cell[loc].Elem.M->PdTrnd = normranf();
	}
	printf("aligned BPM no %1d to %s with BBA"
	       " error x= %f um z= %f um dt= %f urad\n",
	       i, Cell[loc-j].Elem.PName,
	       Cell[loc].Elem.M->PdSrms[0]*Cell[loc].Elem.M->PdSrnd[0]*1e6,
	       Cell[loc].Elem.M->PdSrms[1]*Cell[loc].Elem.M->PdSrnd[0]*1e6,
	       dtor(Cell[loc].Elem.M->PdTrms*Cell[loc].Elem.M->PdTrnd*1e6));
	aligned = true; break;
      } else if ((Cell[loc+j].Elem.Pkind == Mpole) &&
		 (Cell[loc+j].Elem.M->n_design == n)) {
	for (k = 0; k <= 1; k++)
	  Cell[loc].Elem.M->PdSsys[k] = Cell[loc+j].dS[k];
	if (bdxrms >=0.) {
	  Cell[loc].Elem.M->PdSrms[0] = bdxrms;
	  Cell[loc].Elem.M->PdSrnd[0] = normranf();
	} 
	if (bdzrms >=0.) {
	  Cell[loc].Elem.M->PdSrms[1] = bdzrms;
	  Cell[loc].Elem.M->PdSrnd[1] = normranf();
	}
	if (bdarms >=0.) {
	  Cell[loc].Elem.M->PdTrms = bdarms;
	  Cell[loc].Elem.M->PdTrnd = normranf();
	}
	printf("aligned BPM no %1d to %s with BBA"
	       " error x= %f um z= %f um dt= %f urad\n",
	       i,
	       Cell[loc+j].Elem.PName,Cell[loc].Elem.M->PdSrms[0]
	       *Cell[loc].Elem.M->PdSrnd[0]*1e6,
	       Cell[loc].Elem.M->PdSrms[1]*Cell[loc].Elem.M->PdSrnd[0]*1e6,
	       dtor(Cell[loc].Elem.M->PdTrms*Cell[loc].Elem.M->PdTrnd*1e6));
	aligned = true;
	break;
      }

      j++;
    } while (j <= n_step);

    if (aligned) {
      Mpole_SetdS(globval.bpm, i);
      Mpole_SetdT(globval.bpm, i);
    } else
      printf("Align_BPMs: no multipole adjacent to BPM no %d\n", i);
  }
}


void param_data_type::zero_mult(void)
{
  // Extracted to correction/corr_utils; kept as a delegator during the refactor.
  corr::zero_mult(bn_an);
}


void param_data_type::restore_mult(void)
{
  // Extracted to correction/corr_utils; kept as a delegator during the refactor.
  corr::restore_mult(bn_an);
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
      Align_BPMs(Quad,-1.,-1.,-1.);
    }

    // get_traject();
    
    zero_mult();

    cod = CorrectCOD(n_orbit, 1e0);

    restore_mult();

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


void param_data_type::ini_COD_corr
(const int n_bpm_Fam, const std::string bpm_names[],const int n_hcorr_Fam,
 const std::string hcorr_names[], const int n_vcorr_Fam,
 const std::string vcorr_names[], const bool svd)
{
  // Extracted to correction/orbit_corr; kept as a delegator during the refactor.
  corr::ini_COD_corr(n_bpm_Fam, bpm_names, n_hcorr_Fam, hcorr_names,
		     n_vcorr_Fam, vcorr_names, svd);
}


bool param_data_type::cod_corr
(const int n_cell, const double scl,
 const double h_maxkick, const double v_maxkick,
 orb_corr_type orb_corr[])
{
  // Extracted to correction/orbit_corr; kept as a delegator during the refactor.
  return corr::cod_corr(*this, n_cell, scl, h_maxkick, v_maxkick, orb_corr);
}


void param_data_type::Orb_and_Trim_Stat(orb_corr_type orb_corr[])
{
  // Extracted to correction/orbit_corr; kept as a delegator during the refactor.
  corr::Orb_and_Trim_Stat(orb_corr);
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
  double TotalTuneX,TotalTuneY;
  double dk;
  iVector2 nq;
  Vector2 nu;
  fitvect qfbuf, qdbuf;

  double ChromaX,ChromaY;
  double dks;
  iVector2 ns;
  Vector2 si;
  fitvect  sfbuf, sdbuf;

  long i;
  
  globval.Cavity_on   = false;
  globval.radiation   = false;
  globval.Aperture_on = false;

  get_param(param_file);

  Ring_GetTwiss(true, 0.0);
  printglob();

  // Fit tunes to TuneX and TuneY
  if (TuneX*TuneY > 0) {
    printf("\nparam_data_type::err_and_corr_init: fitting nu.\n");
    dk=1e-3;
    nq[0]=nq[1]=0;
    nu[0]=TuneX;
    nu[1]=TuneY;
    for (i = 0; i <= globval.Cell_nLoc; i++) {
      if ( Cell[i].Elem.Pkind == Mpole ) {
	if (strncmp(Cell[i].Elem.PName,"qax",3) == 0){
	  qfbuf[nq[0]]=i;
	  nq[0]++;
	}
	if (strncmp(Cell[i].Elem.PName,"qay",3) == 0){
	  qdbuf[nq[1]]=i;
	  nq[1]++;
	}
      }
    }

    printf("Fittune: nq[0]=%ld nq[1]=%ld\n",nq[0],nq[1]);
    TotalTuneX=globval.TotalTune[0];
    TotalTuneY=globval.TotalTune[1];
    Ring_Fittune(nu, (double)1e-4, nq, qfbuf, qdbuf, dk, 50L);
    printf("Fittune: nux= %f dnux= %f nuy= %f dnuy= %f\n",
	   globval.TotalTune[0], globval.TotalTune[0]-TotalTuneX,
	   globval.TotalTune[1], globval.TotalTune[1]-TotalTuneY);

    Ring_GetTwiss(true, 0.0); printglob();
  }

  // Fit chromaticities to ChromX and ChromY
  if (ChromX*ChromY < 1e6) {
    printf("\nparam_data_type::err_and_corr_init: fitting chi^(1)\n");
    dks=1e-3;
    ns[0]=ns[1]=0;
    si[0]=ChromX;
    si[1]=ChromY;
    for (i = 0; i <= globval.Cell_nLoc; i++) {
      if ( Cell[i].Elem.Pkind == Mpole ) {
	if (strncmp(Cell[i].Elem.PName,"sf",2) == 0){
	  sfbuf[ns[0]]=i;
	  ns[0]++;
	}
	if (strncmp(Cell[i].Elem.PName,"sd",2) == 0){
	  sdbuf[ns[1]]=i;
	  ns[1]++;
	}
      }
    }

    printf("Fitchrom: ns[0]=%ld ns[1]=%ld\n",ns[0],ns[1]);
    ChromaX=globval.Chrom[0];
    ChromaY=globval.Chrom[1];
    Ring_Fitchrom(si, 1e-4, ns, sfbuf, sdbuf, dks, 50L);
    printf("Fitchrom: six= %f dsix= %f siy= %f dsiy= %f\n",
	   globval.Chrom[0], globval.Chrom[0]-ChromaX, globval.Chrom[1],
	   globval.Chrom[1]-ChromaY);

    Ring_GetTwiss(true, 0.0); printglob();
  }

  get_bare();

  cod_ini(bpm_Fam_names, corr_Fam_names, orb_corr);

  if ((ae_file != "") && bba) Align_BPMs(Sext,-1.,-1.,-1.);

  if (N_calls > 0) ini_ID_corr(false);

  if (n_lin > 0) ini_skew_cor(disp_wave_y, disp_wave_o);
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

  snprintf(str, sizeof(str), "n = %d", n);
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
