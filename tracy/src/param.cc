// TODO: naming conventions and refactor - functionally, to set engeniring tolerances and correct including coping with IDs (id_corr is NOT LOCO)

// param_data_type owns get_param and every param.dat knob (declared in param.h):
// reading the full-machine study config is a param-module job, not a correction
// one. The per-corrector slices it projects (corr::orbit_cfg / corr::coupling_cfg
// / corr::bare_optics) live under correction/.

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
  corr::LoadAlignTol(ae_file, true, 1.0, true, k);
  for (i = 1; i <= n_scale; i++) {
    // Scale the rms values
    corr::LoadAlignTol(ae_file, true, (double)i/(double)n_scale, false, k);

    if (bba) {
      // Beam based alignment
      Align_BPMs(Quad,-1.,-1.,-1.);
    }

    // get_traject();
    
    corr::zero_mult(bn_an);

    cod = CorrectCOD(n_orbit, 1e0);

    corr::restore_mult(bn_an);

    if (!cod) break;

    corr::get_dbeta_dnu(m_dbeta, s_dbeta, m_dnu, s_dnu, bare);
    printf("\n");
    printf("RMS dbeta_x/beta_x = %4.2f%%,   dbeta_y/beta_y = %4.2f%%\n",
	   1e2*s_dbeta[X_], 1e2*s_dbeta[Y_]);
    printf("RMS dnu_x          = %7.5f, dnu_y          = %7.5f\n",
	   s_dnu[X_], s_dnu[Y_]);
  }

  return cod;
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


corr::orbit_cfg param_data_type::orbit_config(void) const
{
  corr::orbit_cfg cfg;

  cfg.loc_Fam_name       = loc_Fam_name;
  cfg.bpm_Fam_names      = bpm_Fam_names;
  cfg.corr_Fam_names[X_] = corr_Fam_names[X_];
  cfg.corr_Fam_names[Y_] = corr_Fam_names[Y_];
  cfg.n_thread           = n_thread;
  cfg.n_orbit            = n_orbit;

  return cfg;
}


corr::coupling_cfg param_data_type::coupling_config(void) const
{
  corr::coupling_cfg cfg;

  cfg.VDweight     = VDweight;
  cfg.HVweight     = HVweight;
  cfg.VHweight     = VHweight;
  cfg.qt_s_cut     = qt_s_cut;
  cfg.kick         = kick;
  cfg.n_lin        = n_lin;
  cfg.SQ_per_scell = SQ_per_scell;
  cfg.qt_from_file = qt_from_file;

  return cfg;
}


void param_data_type::get_param(const std::string &param_file)
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
      else if (strcmp("tune_fams", name) == 0) {
	strtok_r(line, " \r", &p);
	s = strtok_r(NULL, " \r", &p); if (s) tune_fam[0] = s;
	s = strtok_r(NULL, " \r", &p); if (s) tune_fam[1] = s;
      } else if (strcmp("chrom_fams", name) == 0) {
	strtok_r(line, " \r", &p);
	s = strtok_r(NULL, " \r", &p); if (s) chrom_fam[0] = s;
	s = strtok_r(NULL, " \r", &p); if (s) chrom_fam[1] = s;
      }
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


void param_data_type::err_and_corr_init(const string &param_file)
{
  globval.Cavity_on   = false;
  globval.radiation   = false;
  globval.Aperture_on = false;

  get_param(param_file);

  Ring_GetTwiss(true, 0.0);
  printglob();

  // Fit tunes to TuneX and TuneY with the tune-quad families (tune_fams keyword,
  // defaulting to the historical SLS-2 qax/qay).
  if (TuneX*TuneY > 0)
    corr::fit_tune(tune_fam[0], tune_fam[1], TuneX, TuneY);

  // Fit chromaticities to ChromX and ChromY with the chroma-sextupole families
  // (chrom_fams keyword, defaulting to the historical SLS-2 sf/sd; on m4U set it
  // to s2_n1/s4_n1 — note s3_n1 overflows fitvect[200], see the parked fit bugs).
  if (ChromX*ChromY < 1e6)
    corr::fit_chrom(chrom_fam[0], chrom_fam[1], ChromX, ChromY);

  bare.capture();

  orbits.alloc(bpm_Fam_names, corr_Fam_names);

  if ((ae_file != "") && bba) Align_BPMs(Sext,-1.,-1.,-1.);

  if (N_calls > 0) id.ini_ID_corr(false, N_Fam, Q_Fam);

  if (n_lin > 0)
    skew.ini_skew_cor(coupling_config(), disp_wave_y, disp_wave_o);
}


void param_data_type::err_and_corr_exit(void)
{
  // The ae_file guard is inherited: with no alignment errors the correctors are
  // still allocated by err_and_corr_init, just never freed. Kept as-is — this
  // step only moves ownership.
  if (ae_file != "") orbits.dealloc();
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
