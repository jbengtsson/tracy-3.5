

bool track(const double x, const double px, const double y, const double py,
	   const double delta, const long int n, const double f_rf,
	   const bool prt)
{
  long int         i, lastpos;
  ss_vect<double>  ps;
  std::ofstream         os;

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

  for (i = 1; i <= n; i++) {
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


void get_r_stable(double &r, const double phi, const double delta,
		  const long int n, const double eps)
{
  /* Binary search for dynamic aperture. */
  bool    lost = false;
  double  r_min = 0.0, r_max = r;

  while (!lost ) {
    lost = ! track(r_max*cos(phi), 0.0, r_max*sin(phi), 0.0, delta,
		   n, 0, false);
    r_max *= 2.0;
  }
  while (r_max-r_min >= eps) {
    r = r_min + (r_max-r_min)/2.0;
    lost = !track(r*cos(phi), 0.0, r*sin(phi), 0.0, delta, n, 0, false);
    if (!lost)
      r_min = r;
    else
      r_max = r;
  }
  r = r_min + (r_max-r_min)/2.0;
}


double get_dynap(FILE *fp, const double r, const double delta, const int n,
		 const double eps, const int n_pts,
		 double x_min[], double x_max[])
{
  /* Determine the dynamic aperture by tracking.
     Assumes mid-plane symmetry.                                              */

  int           i, j;
  double        r1, phi, x0[2] = {0e0, 0e0}, x1[2], x2[2], DA;

  fprintf(fp, "\n");
  fprintf(fp, "# Dynamic Aperture:\n");
  fprintf(fp, "#    x      y\n");
  fprintf(fp, "#   [mm]   [mm]\n");
  fprintf(fp, "#\n");

  for (i = 0; i < 2; i++) {
    x_min[i] = 0.0; x_max[i] = 0.0;
  }

  DA = 0.0; r1 = r;
  for (i = 0; i < n_pts; i++) {
    phi = i*pi/(n_pts-1);
    if (i == 0)
      phi = 1e-3;
    else if (i == n_pts-1)
      phi -= 1e-3;
    get_r_stable(r1, phi, delta, n, eps);
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

void get_DA_bare(const double delta, const int n_delta)
{
  char    str[max_str];
  int     j;
  double  DA, x_min[2], x_max[2], x_hat[2], d;
  FILE    *DA_bare, *fp;

  globval.Cavity_on = true;

  DA_bare = file_write("DA_bare.out");

  fprintf(DA_bare, "# beta_x = %4.2f, beta_y = %5.2f\n",
	  Cell[globval.Cell_nLoc].Beta[X_], Cell[globval.Cell_nLoc].Beta[Y_]);
  fprintf(DA_bare, "#\n");
  fprintf(DA_bare, "# Ideal lattice\n");
  fprintf(DA_bare, "#\n");
  fprintf(DA_bare, "# delta   DA      Ax        Ay      x^    y^\n");
  fprintf(DA_bare, "#  [%%]  [mm^2] [mm.mrad] [mm.mrad] [mm]  [mm]\n");

  for (j = 0; j <= n_delta; j++) {
    d = (n_delta > 0)? (double)j/(double)n_delta*delta : 0.0;

    sprintf(str, "DA_bare_%4.2f.out", 1e2*d); fp = file_write(str);

    DA = get_dynap(fp, 10e-3, d, n_track, 0.1e-3, n_aper, x_min, x_max); 

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


void get_mean_sigma(const int n, double &m, double &s)
{

  m = (n > 0)? m/n : 0e0; 
  s = (n > 1)? sqrt((s-n*sqr(m))/(n-1)) : 0e0;
}
