// Machine error model — see correction/error_model.h.

namespace corr {

// Next whitespace-delimited token from a strtok_r stream; used by LoadFieldErr.
static char* get_prm(char **p)
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


void ReadCorMis(const bool Scale_it, const double Scale)
{
  FILE *fi,*fo;
  const char cormisin[] = "cormis.txt";
  const char cormisout[] = "cormis.out";
  long i,ii;
  CellType Cell;
  double s1, s2, dx, dy, dt;
  char elem[8];
  double dxbn06, dybn06, dtbn06;

  dxbn06=dybn06=dtbn06=0.;

  /* Opening file */
  if ((fo = fopen(cormisout, "w")) == NULL) {
    fprintf(stdout, "cormisout: error while opening file %s\n", cormisout);
    exit_(1);
  }
  /* Opening file */
  if ((fi = fopen(cormisin, "r")) == NULL) {
    fprintf(stdout, "cormisin: error while opening file %s\n", cormisin);
    exit_(1);
  }

  for (i = 0; i <= globval.Cell_nLoc; i++)
  {
    getelem(i, &Cell);
    if (Cell.Elem.Pkind == Mpole)
    {
      fscanf(fi, "%ld %lf %lf %lf %lf %lf %s \n",
	     &ii, &s1, &s2, &dx, &dy, &dt, elem);
      dx/=1e6; dy/=1e6; dt/=1e6;

      if (i == ii) {
	if (Scale_it) {
	  dx *= Scale; dy *= Scale; dt *= Scale;
	}

	if (strcmp("bn06",elem) == 0) {
	  dxbn06=dx; dybn06=dy; dtbn06=dt;
	}
	if ((strcmp("vb",elem) == 0) || (strcmp("vbm",elem)) ==0 ) {
	  dx=dxbn06; dy=dybn06; dt=dtbn06;
	}

        Cell.Elem.M->PdSsys[0] = dx;
        Cell.Elem.M->PdSsys[1] = dy;
        Cell.Elem.M->PdTsys    = dt;

        putelem(ii, &Cell);
	Mpole_SetdS(Cell.Fnum, Cell.Knum);
	Mpole_SetdT(Cell.Fnum, Cell.Knum);

        fprintf(fo, "%ld %lf %lf %lf %lf %lf %s \n",
		ii,  s1, s2, dx*1e6, dy*1e6, dt*1e6, Cell.Elem.PName);
      }
    }
  }
  fclose(fo);
}

void LoadAlignTol(const std::string &ae_file, const bool Scale_it,
		  const double Scale, const bool new_rnd, const int seed)
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


void LoadFieldErr(const std::string &fe_file, const bool Scale_it,
		  const double Scale, const bool new_rnd)
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


void LoadApers(const std::string &ap_file, const double scl_x,
	       const double scl_y)
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

}  // namespace corr
