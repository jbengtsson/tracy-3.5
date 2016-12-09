
bool track(const double x, const double px, const double y, const double py,
	   const double delta, const long int n, const double f_rf,
	   const bool prt);

void get_r_stable(double &r, const double phi, const double delta,
		  const long int n, const double eps);

double get_dynap(FILE *fp, const double r, const double delta, const int n,
		 const double eps, const int n_pts,
		 double x_min[], double x_max[]);

void get_DA_bare(const double delta, const int n_delta);

void get_mean_sigma(const int n, double &m, double &s);
