

bool DA_data_type::track(const param_data_type &params,
			 const double x, const double px,
			 const double y, const double py,
			 const double delta,
			 const double f_rf, const bool prt)
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

  for (i = 1; i <= n_track; i++) {
    Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
    if (lastpos == globval.Cell_nLoc) {
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


void DA_data_type::get_r_stable(const param_data_type &params,
				double &r, const double phi,
				const double delta, const double eps)
{
  /* Binary search for dynamic aperture. */
  bool   lost = false;
  double r_min = 0.0, r_max = r;

  while (!lost ) {
    lost = ! track(params, r_max*cos(phi), 0.0, r_max*sin(phi), 0.0, delta,
		   0, false);
    r_max *= 2.0;
  }
  while (r_max-r_min >= eps) {
    r = r_min + (r_max-r_min)/2.0;
    lost = !track(params, r*cos(phi), 0.0, r*sin(phi), 0.0, delta, 0, false);
    if (!lost)
      r_min = r;
    else
      r_max = r;
  }
  r = r_min + (r_max-r_min)/2.0;
}


double DA_data_type::get_dynap(param_data_type &params,
			       FILE *fp, const double r, const double delta,
			       const double eps, double x_min[], double x_max[])
{
  /* Determine the dynamic aperture by tracking.
     Assumes mid-plane symmetry.                                              */

  int    i, j;
  double r1, phi, x0[2] = {0e0, 0e0}, x1[2], x2[2], DA;

  fprintf(fp, "\n");
  fprintf(fp, "# Dynamic Aperture:\n");
  fprintf(fp, "#    x      y\n");
  fprintf(fp, "#   [mm]   [mm]\n");
  fprintf(fp, "#\n");

  for (i = 0; i < 2; i++) {
    x_min[i] = 0.0; x_max[i] = 0.0;
  }

  DA = 0.0; r1 = r;
  for (i = 0; i < n_aper; i++) {
    phi = i*pi/(n_aper-1);
    if (i == 0)
      phi = 1e-3;
    else if (i == n_aper-1)
      phi -= 1e-3;
    get_r_stable(params, r1, phi, delta, eps);
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
  DA = fabs(DA)/sqrt(Cell[globval.Cell_nLoc].Beta[X_]
       *Cell[globval.Cell_nLoc].Beta[Y_]);

  fprintf(fp, "\n");
  fprintf(fp, "# DA^ = %6.1f mm^2"
	  ", x^ = %6.2f - %5.2f mm, y^ = %6.2f - %5.2f mm\n",
	  1e6*DA, 1e3*x_min[X_], 1e3*x_max[X_], 1e3*x_min[Y_], 1e3*x_max[Y_]);

  fflush(fp);

  return DA;
} 


void DA_data_type::get_DA_bare(param_data_type &params)
{
  char   str[max_str];
  int    j;
  double DA, x_min[2], x_max[2], x_hat[2], d;
  FILE   *DA_bare, *fp;

  globval.Cavity_on = true;

  DA_bare = file_write("DA_bare.out");

  fprintf(DA_bare, "# beta_x = %4.2f, beta_y = %5.2f\n",
	  Cell[globval.Cell_nLoc].Beta[X_], Cell[globval.Cell_nLoc].Beta[Y_]);
  fprintf(DA_bare, "#\n");
  fprintf(DA_bare, "# Ideal lattice\n");
  fprintf(DA_bare, "#\n");
  fprintf(DA_bare, "# delta   DA      Ax        Ay      x^    y^\n");
  fprintf(DA_bare, "#  [%%]  [mm^2] [mm.mrad] [mm.mrad] [mm]  [mm]\n");

  for (j = 0; j <= params.n_delta_DA; j++) {
    d = (params.n_delta_DA > 0)?
      (double)j/(double)params.n_delta_DA*params.delta_DA : 0.0;

    sprintf(str, "DA_bare_%4.2f.out", 1e2*d); fp = file_write(str);

    DA = get_dynap(params, fp, 10e-3, d, 0.1e-3, x_min, x_max); 

    fclose(fp);

    x_hat[X_] = (x_max[X_]-x_min[X_])/2.0; x_hat[Y_] = x_max[Y_];

    fprintf(DA_bare, "  %5.2f %6.1f   %4.1f      %4.1f   %4.1f  %4.1f\n", 
	    1e2*d, 1e6*DA,
	    1e6*sqr(x_hat[X_])/Cell[globval.Cell_nLoc].Beta[X_],
	    1e6*sqr(x_hat[Y_])/Cell[globval.Cell_nLoc].Beta[Y_],
	    1e3*x_hat[X_], 1e3*x_hat[Y_]);
  
    fflush(DA_bare);
  }

  fclose(DA_bare);
}


void DA_data_type::get_DA_real(param_data_type &params,
			       orb_corr_type orb_corr[])
{
  bool   cod = false;
  char   str[max_str];
  int    j, k;
  double DA, x_min[2], x_max[2], d[params.n_delta_DA+1], x_hat;
  double DA_m[params.n_delta_DA+1], DA_s[params.n_delta_DA+1];
  double x_hat_m[params.n_delta_DA+1][2], x_hat_s[params.n_delta_DA+1][2];

  const int    n_cell = 20;
  const double scl    = 0.5;
  
  FILE *DA_real = NULL, *fp[params.n_delta_DA+1];

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
	  Cell[globval.Cell_nLoc].Beta[X_],
	  Cell[globval.Cell_nLoc].Beta[Y_]);
  fprintf(DA_real, "#\n");
  fprintf(DA_real, "# Real lattice\n");
  fprintf(DA_real, "#\n");
  fprintf(DA_real, "# deltaP        DA               Ax              Ay"
	  "            x^              y^\n");
  fprintf(DA_real, "#   %%          mm^2              mm              mm"
	    "            mm              mm\n");

  for (j = 1; j <= params.n_stat; j++) {
    globval.Cavity_on = false;

    if (params.fe_file != "") params.LoadFieldErr(false, 1e0, true);

    if (params.ae_file != "") cod = params.cod_corr(n_cell, scl, j, orb_corr);
      
    // cod = CorrectCOD(n_orbit, 0.3);

    printf("\n");
    if (cod) {
      printf("err_and_corr: orbit correction completed\n");
      if (j == 1) prt_cod("cod.out", globval.bpm, true);
 
      Ring_GetTwiss(true, 0.0); printglob();

      GetEmittance(ElemIndex("cav"), true);

      if (params.n_lin > 0) params.corr_eps_y();

      if (params.ap_file != "") params.LoadApers(1.0, 1.0);

      Ring_GetTwiss(true, 0.0); printglob();

      GetEmittance(ElemIndex("cav"), true);

      prt_beamsizes();

      globval.Cavity_on = true;
      for (k = 0; k <= params.n_delta_DA; k++) {
	DA = get_dynap(params, fp[k], 10e-3, d[k], 0.1e-3, x_min, x_max);
	DA_m[k] += DA; DA_s[k] += sqr(DA);
	x_hat = (x_max[X_]-x_min[X_])/2.0;
	x_hat_m[k][X_] += x_hat; x_hat_s[k][X_] += sqr(x_hat);
	x_hat = x_max[Y_];
	x_hat_m[k][Y_] += x_hat; x_hat_s[k][Y_] += sqr(x_hat);
      }

     if (params.n_lin > 0)
	// reset skew quads
	// set_bn_design_fam(globval.qt, Quad, 0.0, 0.0);

      if (params.N_calls > 0) params.reset_quads();  
    } else
      chk_cod(cod, "error_and_correction");
  }

  for (j = 0; j <= params.n_delta_DA; j++) {
    fclose(fp[j]);

    get_mean_sigma(params.n_stat, DA_m[j], DA_s[j]);
    for (k = 0; k <= 1; k++)
      get_mean_sigma(params.n_stat, x_hat_m[j][k], x_hat_s[j][k]);

    fprintf(DA_real, "  %5.2f  %6.1f \xB1 %6.1f"
	    "   %5.1f \xB1 %5.2f    %5.1f \xB1 %5.2f"
	    "   %4.1f \xB1 %5.2f     %4.1f \xB1 %5.2f\n", 
	    d[j]*1e2,
	    1e6*DA_m[j], 1e6*DA_s[j],
	    1e6*sqr(x_hat_m[j][X_])/Cell[globval.Cell_nLoc].Beta[X_],
	    1e6*sqr(x_hat_s[j][X_])/Cell[globval.Cell_nLoc].Beta[X_],
	    1e6*sqr(x_hat_m[j][Y_])/Cell[globval.Cell_nLoc].Beta[Y_],
	    1e6*sqr(x_hat_s[j][Y_])/Cell[globval.Cell_nLoc].Beta[Y_],
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
