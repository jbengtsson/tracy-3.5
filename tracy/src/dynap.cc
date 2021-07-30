

bool DA_data_type::track(LatticeType &lat, const param_data_type &params,
			 const double x, const double px, const double y,
			 const double py, const double delta, const double f_rf,
			 const bool prt)
{
  long int        i, lastpos;
  ss_vect<double> ps;
  std::ofstream   os;

  ps[x_] = x; ps[px_] = px; ps[y_] = y; ps[py_] = py;
  ps[delta_] = delta; ps[ct_] = 0.0;

  if (prt) {
    os.open("track.out", std::ios::out);
    os << "# Tracking with Thor" << std::endl;
    os << "#" << std::endl;
    os << "#  n       x           p_x          y            p_y  "
       << "       delta         cdt" << std::endl;
    if (f_rf == 0.0) {
      os << "#         [mm]        [mrad]       [mm]         [mrad]"
	 << "                    [mm]" << std::endl;
      os << std::scientific << std::setprecision(16)
	 << std::setw(4) << 0
	 << std::setw(24) << 1e3*ps[x_] << std::setw(24) << 1e3*ps[px_]
	 << std::setw(24) << 1e3*ps[y_] << std::setw(24) << 1e3*ps[py_]
	 << std::setw(24) << 1e2*ps[delta_] 
	 << std::setw(24) << 1e3*ps[ct_] << std::endl;
    } else {
      os << "#         [mm]        [mrad]       [mm]         [mrad]"
	 << "                    [deg]" << std::endl;
      os << std::scientific << std::setprecision(16)
	 << std::setw(4) << 0
	 << std::setw(24) << 1e3*ps[x_] << std::setw(24) << 1e3*ps[px_]
	 << std::setw(24) << 1e3*ps[y_] << std::setw(24) << 1e3*ps[py_]
	 << std::setw(24) << 1e2*ps[delta_] 
	 << std::setw(24) << 2.0*f_rf*180.0*ps[ct_]/c0 << std::endl;
    }
    os << "#" << std::endl;
  }

  for (i = 1; i <= params.n_track_DA; i++) {
    lat.Cell_Pass(0, lat.conf.Cell_nLoc, ps, lastpos);
    if (lastpos == lat.conf.Cell_nLoc) {
      if (prt) {
	if (f_rf == 0.0)
	  os << std::scientific << std::setprecision(16)
	     << std::setw(4) << i
	     << std::setw(24) << 1e3*ps[x_] << std::setw(24) << 1e3*ps[px_]
	     << std::setw(24) << 1e3*ps[y_] << std::setw(24) << 1e3*ps[py_]
	     << std::setw(24) << 1e2*ps[delta_] 
	     << std::setw(24) << 1e3*ps[ct_] << std::endl;
	else
	  os <<std:: scientific << std::setprecision(16)
	     << std::setw(4) << i
	     << std::setw(24) << 1e3*ps[x_] << std::setw(24) << 1e3*ps[px_]
	     << std::setw(24) << 1e3*ps[y_] << std::setw(24) << 1e3*ps[py_]
	     << std::setw(24) << 1e2*ps[delta_] 
	     << std::setw(24) << 2.0*f_rf*180.0*ps[ct_]/c0 << std::endl;
      }
    } else
      return false;
  }
  if (prt) os.close();

  return true;
}


void DA_data_type::get_r_stable(LatticeType &lat, const param_data_type &params,
				double &r, const double phi,
				const double delta, const double eps)
{
  /* Binary search for dynamic aperture. */
  bool   lost = false;
  double r_min = 0.0, r_max = r;

  while (!lost ) {
    lost = ! track(lat, params, r_max*cos(phi), 0.0, r_max*sin(phi), 0.0, delta,
		   0, false);
    r_max *= 2.0;
  }
  while (r_max-r_min >= eps) {
    r = r_min + (r_max-r_min)/2.0;
    lost = !track(lat, params, r*cos(phi), 0.0, r*sin(phi), 0.0, delta, 0,
		  false);
    if (!lost)
      r_min = r;
    else
      r_max = r;
  }
  r = r_min + (r_max-r_min)/2.0;
}


double DA_data_type::get_dynap(LatticeType &lat, param_data_type &params,
			       FILE *fp, const double r, const double delta,
			       const double eps, double x_min[], double x_max[])
{
  /* Determine the dynamic aperture by tracking.
     Assumes mid-plane symmetry.                                              */

  int    i, j;
  double r1, phi, x0[2] = {0e0, 0e0}, x1[2] = {0e0, 0e0}, x2[2] = {0e0, 0e0};
  double DA;

  fprintf(fp, "\n");
  fprintf(fp, "# Dynamic Aperture:\n");
  fprintf(fp, "#    x      y\n");
  fprintf(fp, "#   [mm]   [mm]\n");
  fprintf(fp, "#\n");

  for (i = 0; i < 2; i++) {
    x_min[i] = 0.0; x_max[i] = 0.0;
  }

  DA = 0.0; r1 = r;
  for (i = 0; i < params.n_aper_DA; i++) {
    phi = i*pi/(params.n_aper_DA-1);
    if (i == 0)
      phi = 1e-3;
    else if (i == params.n_aper_DA-1)
      phi -= 1e-3;
    get_r_stable(lat, params, r1, phi, delta, eps);
    x2[X_] = r1*cos(phi); x2[Y_] = r1*sin(phi);
    for (j = 0; j <= 1; j++) {
      x_min[j] = min(x2[j], x_min[j]); x_max[j] = max(x2[j], x_max[j]);
    }
    if (i == 0) {
      x0[X_] = x2[X_]; x0[Y_] = x2[Y_];
    } else
      DA += x1[X_]*x2[Y_] - x2[X_]*x1[Y_];

    fprintf(fp, "  %6.2f %6.2f\n", 1e3*x2[X_], 1e3*x2[Y_]);

    x1[X_] = x2[X_]; x1[Y_] = x2[Y_];
  }
  DA += x2[X_]*x0[Y_] - x0[X_]*x2[Y_];
  // x2 from mid-plane symmetry
  DA = fabs(DA)/sqrt(lat.elems[lat.conf.Cell_nLoc]->Beta[X_]
       *lat.elems[lat.conf.Cell_nLoc]->Beta[Y_]);

  fprintf(fp, "\n");
  fprintf(fp, "# DA^ = %6.1f mm^2"
	  ", x^ = %6.2f - %5.2f mm, y^ = %6.2f - %5.2f mm\n",
	  1e6*DA, 1e3*x_min[X_], 1e3*x_max[X_], 1e3*x_min[Y_], 1e3*x_max[Y_]);

  fflush(fp);

  return DA;
} 


void DA_data_type::get_DA_bare(LatticeType &lat, param_data_type &params)
{
  char   str[max_str];
  int    j;
  double DA, x_min[2], x_max[2], x_hat[2], d;
  FILE   *DA_bare, *fp;

  lat.conf.Cavity_on = true;

  DA_bare = file_write("DA_bare.out");

  fprintf(DA_bare, "# beta_x = %4.2f, beta_y = %5.2f\n",
	  lat.elems[lat.conf.Cell_nLoc]->Beta[X_], lat.elems[lat.conf.Cell_nLoc]->Beta[Y_]);
  fprintf(DA_bare, "#\n");
  fprintf(DA_bare, "# Ideal lattice\n");
  fprintf(DA_bare, "#\n");
  fprintf(DA_bare, "# delta   DA      Ax        Ay      x^    y^\n");
  fprintf(DA_bare, "#  [%%]  [mm^2] [mm.mrad] [mm.mrad] [mm]  [mm]\n");

  for (j = 0; j <= params.n_delta_DA; j++) {
    d = (params.n_delta_DA > 0)?
      (double)j/(double)params.n_delta_DA*params.delta_DA : 0.0;

    sprintf(str, "DA_bare_%4.2f.out", 1e2*d); fp = file_write(str);

    DA = get_dynap(lat, params, fp, 10e-3, d, 0.1e-3, x_min, x_max); 

    fclose(fp);

    x_hat[X_] = (x_max[X_]-x_min[X_])/2.0; x_hat[Y_] = x_max[Y_];

    fprintf(DA_bare, "  %5.2f %6.1f   %4.1f      %4.1f   %4.1f  %4.1f\n", 
	    1e2*d, 1e6*DA,
	    1e6*sqr(x_hat[X_])/lat.elems[lat.conf.Cell_nLoc]->Beta[X_],
	    1e6*sqr(x_hat[Y_])/lat.elems[lat.conf.Cell_nLoc]->Beta[Y_],
	    1e3*x_hat[X_], 1e3*x_hat[Y_]);
  
    fflush(DA_bare);
  }

  fclose(DA_bare);
}


void DA_data_type::get_DA_real(LatticeType &lat, param_data_type &params,
			       orb_corr_type orb_corr[])
{
  bool     cod = false;
  char     str[max_str];
  long int lastpos;
  int      j, k;
  double   DA, x_min[2], x_max[2], d[params.n_delta_DA+1], x_hat;
  double   DA_m[params.n_delta_DA+1], DA_s[params.n_delta_DA+1];
  double   x_hat_m[params.n_delta_DA+1][2], x_hat_s[params.n_delta_DA+1][2];
  long     i;
  double   gdxrms, gdzrms, gdarms, jdxrms, jdzrms, edxrms, edzrms, edarms;
  double   bdxrms, bdzrms, bdarms;
  double   rancutx, rancuty, rancutt;
  long     iseednr, iseed[iseednrmax]={0L};
  char     fname[30];

  fitvect  qfbuf, qdbuf;
  std::vector<int>
    nq = {0, 0},
    ns = {0, 0};
  std::vector<double>
    nu = {0.0, 0.0},
    si = {0.0, 0.0};
  double   dk;
  double   TotalTuneX,TotalTuneY;

  fitvect  sfbuf, sdbuf;
  double   dks;
  double   ChromaX, ChromaY;

  ElemType cell, *WITH;
  FILE     *DA_real = NULL, *fp[params.n_delta_DA+1];

  const int n_cell = 20;

  gdxrms = gdzrms = gdarms = jdxrms = jdzrms = edxrms = edzrms = edarms
         = bdxrms = bdzrms = bdarms = 0e0;
  
  for (j = 0; j <= params.n_delta_DA; j++) {
    d[j] = (params.n_delta_DA > 0)?
      (double)j/(double)params.n_delta_DA*params.delta_DA : 0.0;
    sprintf(str, "DA_real_%4.2f.out", 1e2*d[j]); fp[j] = file_write(str);

    DA_m[j] = 0e0; DA_s[j] = 0e0;

    for (k = 0; k <= 1; k++) {
      x_hat_m[j][k] = 0e0; x_hat_s[j][k] = 0e0;
    }
  }

  DA_real = file_write("DA_real.out");
  fprintf(DA_real, "# beta_x = %4.2f, beta_y = %5.2f\n",
	  lat.elems[lat.conf.Cell_nLoc]->Beta[X_],
	  lat.elems[lat.conf.Cell_nLoc]->Beta[Y_]);
  fprintf(DA_real, "#\n");
  fprintf(DA_real, "# Real lattice\n");
  fprintf(DA_real, "#\n");
  fprintf(DA_real, "# deltaP        DA               Ax              Ay"
	  "            x^              y^\n");
  fprintf(DA_real, "#   %%          mm^2              mm              mm"
	    "            mm              mm\n");

  if (params.n_meth == 1) {
    printf("Entering GirderSetup\n");
    params.GirderSetup(lat);
  }

  for (j = 1; j <= params.n_stat; j++) {
    lat.conf.Cavity_on = false;

    if (params.fe_file != "") params.LoadFieldErr(lat, false, 1e0, true);
    if (params.ae_file != "") {
      // Load misalignments; set seed, no scaling of rms errors.
      if (trace) printf("get_DA_real: n_meth = %d", params.n_meth);
      if (params.n_meth == 0) {
        printf("entering LoadAlignTol\n");
        params.LoadAlignTol(lat, false, 1e0, true, j);
	bdxrms = bdzrms = bdarms = -1e0;
      } else if (params.n_meth == 1) {
	// printf("entering ReadCormis\n");
	// params.ReadCorMis(false,1e0);
	printf("Entering CorMis_in\n");
	params.CorMis_in(&gdxrms, &gdzrms, &gdarms, &jdxrms, &jdzrms, &edxrms,
			 &edzrms, &edarms, &bdxrms, &bdzrms, &bdarms, &rancutx,
			 &rancuty, &rancutt, iseed, &iseednr);
	if (params.n_stat > iseednr) {
	  printf("n_stat %d exceeds iseednr %ld\n", params.n_stat, iseednr);
	  exit(1);
	}
	printf("Entering SetCorMis\n");
	params.SetCorMis(lat, gdxrms, gdzrms, gdarms, jdxrms, jdzrms, edxrms,
			 edzrms, edarms, rancutx, rancuty, rancutt, iseed[j-1]);
      }

      // Beam based alignment with respect to sextupoles with errors bdxrms,
      // bdzrms, bdarms
      if (params.bba) {
        params.Align_BPMs(lat, Sext, bdxrms, bdzrms, bdarms);
      }
      cod = params.cod_corr(lat, n_cell, 1e0, params.h_maxkick,
			    params.v_maxkick, params.n_bits, orb_corr);
    } else
      cod = lat.getcod(0e0, lastpos);

    params.Orb_and_Trim_Stat(lat, orb_corr);

    if (params.N_calls > 0) {
      params.ID_corr(lat, params.N_calls, params.N_steps, false, j);
      cod = params.cod_corr(lat, n_cell, 1e0, params.h_maxkick,
			    params.v_maxkick, params.n_bits, orb_corr);
    }

    params.Orb_and_Trim_Stat(lat, orb_corr);

    printf("\n");
    if (cod) {
      printf("err_and_corr: orbit correction completed\n");

      sprintf(fname,"linlat_%d.out",j);
      lat.prt_lat(fname, true);
      sprintf(fname,"cod_%d.out",j);
      prt_cod(lat, fname, true);
      sprintf(fname,"cod_%d.dat",j);
      printcod(lat, fname);
      if (trace && (j == 1)) {
	orb_corr[X_].prt_svdmat(lat);
	orb_corr[Y_].prt_svdmat(lat);
      }
 
      lat.Ring_GetTwiss(true, 0.0); printglob(lat);

      GetEmittance(lat, ElemIndex("cav"), true);

      if (params.n_lin > 0) {
	params.corr_eps_y(lat, j);
	if (params.N_calls > 0) {
	  params.ID_corr(lat, params.N_calls, params.N_steps, false, j);
	  params.cod_corr(lat, n_cell, 1e0, params.h_maxkick, params.v_maxkick,
			  params.n_bits, orb_corr);
	}
 	lat.Ring_GetTwiss(true, 0.0); printglob(lat);
	GetEmittance(lat, ElemIndex("cav"), true);
      }

      ///////////////////////////////
      // Fit tunes to TuneX and TuneY
      
      if (params.TuneX*params.TuneY > 0) {
	dk=1e-3;
	nq[0]=nq[1]=0;
	nu[0]=params.TuneX;
	nu[1]=params.TuneY;
	for (i = 0; i <= lat.conf.Cell_nLoc; i++) {
	  WITH = lat.elems[i];
	  if ( WITH->Pkind == Mpole ) {
	    if (strncmp(lat.elems[i]->PName,"qax",3) == 0){
	      qfbuf[nq[0]]=i;
	      nq[0]++;
	    }
	    if (strncmp(lat.elems[i]->PName,"qay",3) == 0){
	      qdbuf[nq[1]]=i;
	      nq[1]++;
	    }
	  }
	}

	printf("Fittune: nq[0]=%ld nq[1]=%ld\n",nq[0],nq[1]);
	TotalTuneX=lat.conf.TotalTune[0];
	TotalTuneY=lat.conf.TotalTune[1];
	lat.Ring_Fittune(nu, (double)1e-4, nq, qfbuf, qdbuf, dk, 50L);
	printf("Fittune: nux= %f dnux= %f nuy= %f dnuy= %f\n",
	       lat.conf.TotalTune[0], lat.conf.TotalTune[0]-TotalTuneX,
	       lat.conf.TotalTune[1], lat.conf.TotalTune[1]-TotalTuneY);

	lat.Ring_GetTwiss(true, 0.0); printglob(lat);
	GetEmittance(lat, ElemIndex("cav"), true);
      }

      // Fit chromaticities to ChromX and ChromY

      if (params.ChromX*params.ChromY < 1e6) {
	dks=1e-3;
	ns[0]=ns[1]=0;
	si[0]=params.ChromX;
	si[1]=params.ChromY;
	for (i = 0; i <= lat.conf.Cell_nLoc; i++) {
	  WITH = lat.elems[i];
	  if ( WITH->Pkind == Mpole ) {
	    if (strncmp(lat.elems[i]->PName,"sf",2) == 0){
	      sfbuf[ns[0]]=i;
	      ns[0]++;
	    }
	    if (strncmp(lat.elems[i]->PName,"sd",2) == 0){
	      sdbuf[ns[1]]=i;
	      ns[1]++;
	    }
	  }
	}

	printf("Fitchrom: ns[0]=%ld ns[1]=%ld\n",ns[0],ns[1]);
	ChromaX=lat.conf.Chrom[0];
	ChromaY=lat.conf.Chrom[1];
	lat.Ring_Fitchrom(si, (double)1e-4, ns, sfbuf, sdbuf, dks, 50L);
	printf("Fitchrom: six= %f dsix= %f siy= %f dsiy= %f\n",
	       lat.conf.Chrom[0], lat.conf.Chrom[0]-ChromaX, lat.conf.Chrom[1],
	       lat.conf.Chrom[1]-ChromaY);

	lat.Ring_GetTwiss(true, 0.0); printglob(lat);
	GetEmittance(lat, ElemIndex("cav"), true);
      }
      
      // End of tune and chromaticity fit
      ///////////////////////////////////
      
      prt_beamsizes(lat, j);

      if (params.ap_file != "") params.LoadApers(lat, 1.0, 1.0);

      lat.conf.Cavity_on = true;

      // Define multipoles.
      // setmpall(0.01);

      for (i = 0; i <= lat.conf.Cell_nLoc; i++)
	lat.getelem(i, &cell);
      if (cell.Pkind == Mpole)
	printf("%ld %lf %lf %lf %s \n",
	       i, cell.dS[0]*1e6, cell.dS[1]*1e6, cell.dT[1]*1e6,
	       cell.PName);

      for (k = 0; k <= params.n_delta_DA; k++) {
	DA = get_dynap(lat, params, fp[k], 10e-3, d[k], 0.1e-3, x_min, x_max);
	DA_m[k] += DA; DA_s[k] += sqr(DA);
	x_hat = (x_max[X_]-x_min[X_])/2.0;
	x_hat_m[k][X_] += x_hat; x_hat_s[k][X_] += sqr(x_hat);
	x_hat = x_max[Y_];
	x_hat_m[k][Y_] += x_hat; x_hat_s[k][Y_] += sqr(x_hat);
      }

      if (params.n_lin > 0) {
	// reset skew quads
	printf("resetting skew quad family: %s\n",
	       lat.elems[lat.Elem_GetPos(lat.conf.qt,1)]->PName);
        set_bnL_design_fam(lat, lat.conf.qt, Quad, 0.0, 0.0);
      }
      if (params.N_calls > 0) params.reset_quads(lat);  
    } else
      chk_cod(cod, "err_and_corr");
  }

  for (j = 0; j <= params.n_delta_DA; j++) {
    fclose(fp[j]);

    get_mean_sigma(params.n_stat, DA_m[j], DA_s[j]);
    for (k = 0; k <= 1; k++)
      get_mean_sigma(params.n_stat, x_hat_m[j][k], x_hat_s[j][k]);

    fprintf(DA_real, "  %6.3f  %7.2f \xB1 %7.2f"
	    "   %6.2f \xB1 %6.3f    %6.2f \xB1 %6.3f"
	    "   %5.2f \xB1 %6.3f     %5.2f \xB1 %6.3f\n", 
	    d[j]*1e2,
	    1e6*DA_m[j], 1e6*DA_s[j],
	    1e6*sqr(x_hat_m[j][X_])/lat.elems[lat.conf.Cell_nLoc]->Beta[X_],
	    1e6*sqr(x_hat_s[j][X_])/lat.elems[lat.conf.Cell_nLoc]->Beta[X_],
	    1e6*sqr(x_hat_m[j][Y_])/lat.elems[lat.conf.Cell_nLoc]->Beta[Y_],
	    1e6*sqr(x_hat_s[j][Y_])/lat.elems[lat.conf.Cell_nLoc]->Beta[Y_],
	    1e3*x_hat_m[j][X_], 1e3*x_hat_s[j][X_],
	    1e3*x_hat_m[j][Y_], 1e3*x_hat_s[j][Y_]);
  }

  fclose(DA_real);
}


void DA_data_type::get_mean_sigma(const int n, double &m, double &s)
{
  m = (n > 0)? m/n : 0e0; 
  s = (n > 1)? sqrt((s-n*sqr(m))/(n-1)) : 0e0;
}
