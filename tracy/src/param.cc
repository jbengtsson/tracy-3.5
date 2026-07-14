// TODO: naming conventions and refactor - functionally, to set engeniring tolerances and correct including coping with IDs (id_corr is NOT LOCO)

// The param.dat knobs and their defaults moved to correction/corr_config.h as
// corr::config_data members (they were param_data_type statics); get_param moved
// to correction/corr_config.cc. param_data_type derives from config_data, so its
// callers still reach the knobs as params.<knob>. See correction_refactor.md.

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

void param_data_type::get_bare(void)
{
  // Extracted to correction/config; kept as a delegator during the refactor.
  bare.capture();
}


void param_data_type::get_dbeta_dnu(double m_dbeta[], double s_dbeta[],
				    double m_dnu[], double s_dnu[])
{
  // Extracted to correction/corr_utils; kept as a delegator during the refactor.
  corr::get_dbeta_dnu(m_dbeta, s_dbeta, m_dnu, s_dnu, bare);
}


// Delegates to correction/orbit_corr.


void param_data_type::ini_skew_cor(const double deta_y_max,
				   const double deta_y_offset)
{
  skew.ini_skew_cor(coupling_config(), deta_y_max, deta_y_offset);
}


void param_data_type::corr_eps_y(const int cnt)
{
  skew.corr_eps_y(coupling_config(), cnt);
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
  return corr::cod_corr(orbit_config(), bare, n_cell, scl, h_maxkick, v_maxkick,
			orb_corr);
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
