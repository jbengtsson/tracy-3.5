// Correction configuration and reference state — see correction/corr_config.h.

namespace corr {

void bare_optics::capture(void)
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
}


orbit_cfg config_data::orbit_config(void) const
{
  orbit_cfg cfg;

  cfg.loc_Fam_name       = loc_Fam_name;
  cfg.bpm_Fam_names      = bpm_Fam_names;
  cfg.corr_Fam_names[X_] = corr_Fam_names[X_];
  cfg.corr_Fam_names[Y_] = corr_Fam_names[Y_];
  cfg.n_thread           = n_thread;
  cfg.n_orbit            = n_orbit;

  return cfg;
}


coupling_cfg config_data::coupling_config(void) const
{
  coupling_cfg cfg;

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


void config_data::get_param(const std::string &param_file)
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

}  // namespace corr
