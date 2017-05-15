#define NO 1

#include "tracy_lib.h"

int no_tps = NO;

const int  n_bpm_max = 100;

char    bpm_names[n_bpm_max][max_str];
int     n_bpm, n_turn_, n_stats, jj[ss_dim];
double  data[n_bpm_max][2][2048];
double  alpha_mean[2], alpha_sigma[2];
double  betas_sum[n_bpm_max][2], betas_sum2[n_bpm_max][2];
double  betas_mean[n_bpm_max][2], betas_sigma[n_bpm_max][2];
double  dnus_sum[n_bpm_max][2], dnus_sum2[n_bpm_max][2];
double  dnus_mean[n_bpm_max][2], dnus_sigma[n_bpm_max][2];
double  twoJ[n_bpm_max][2][2048], phi[n_bpm_max][2][2048];
double  phi0[n_bpm_max][2][2048];
double  tune_mean[2], tune_sigma[2];

// Kalman filter.
ss_vect<tps>  Id, A, A_tp, H, H_tp, R, Q, P, K_;

ofstream  outf_optics;


void get_bpm_name(char *name)
{
  int  k;

  k = 0;
  do {
    if (name[k] == '-')
      name[k] = '_';
    else
      name[k] = tolower(name[k]);

    k++;
  } while (name[k] != '\0');

}


void rd_tbt(const char *file_name)
{
  const int  str_len = 1500;

  char      line[str_len];
  int       j, k;
  ifstream  inf;

  const bool  prt = false;
  const int   n_print = 8;

  inf.open(file_name);

  inf.getline(line, str_len);
  inf.getline(line, str_len); sscanf(line, "%d %d", &n_bpm, &n_turn_);

  cout << endl;
  cout << "no of BPMs = " << n_bpm << ", no of turns = " << n_turn_ << endl;

  inf.getline(line, str_len);

  if (prt) cout << endl;
  for (j = 0; j < n_bpm; j++) {
    inf.getline(line, str_len); sscanf(line, "%s", bpm_names[j]);
    get_bpm_name(bpm_names[j]);
    if (prt) {
      cout << " " << bpm_names[j];
      if ((j+1) % n_print == 0) cout << endl;
    }
  }
  if (prt) if (n_bpm % n_print != 0) cout << endl;

  inf.getline(line, str_len);

  if (prt) cout << endl;
  for (k = 0; k < n_turn_; k++) {
    if (prt) cout << endl;
    inf.getline(line, str_len);
    for (j = 0; j < n_bpm-1; j++) {
      if (j == 0)
	data[j][X_][k]= 1e-3*atof(strtok(line, " "));
      else
	data[j][X_][k] = 1e-3*atof(strtok(NULL, " "));

      if (prt) {
	cout << fixed << setprecision(6) << setw(10) << 1e3*data[j][X_][k];
	if ((j+1) % n_print == 0) cout << endl;
      }
    }
    j = n_bpm - 1;
    inf.getline(line, str_len);
    data[j][X_][k]= 1e-3*atof(strtok(line, " "));
    if (prt) {
      cout << fixed << setprecision(6) << setw(10) << 1e3*data[j][X_][k];
      if (k % n_print != 0) cout << endl;
    }
  }


  inf.getline(line, str_len);

  if (prt) cout << endl;
  for (k = 0; k < n_turn_; k++) {
    if (prt) cout << endl;
    inf.getline(line, str_len);
    for (j = 0; j < n_bpm-1; j++) {
      if (j == 0)
	data[j][Y_][k]= 1e-3*atof(strtok(line, " "));
      else
	data[j][Y_][k] = 1e-3*atof(strtok(NULL, " "));

      if (prt) {
	cout << fixed << setprecision(6) << setw(10) << 1e3*data[j][Y_][k];
	if ((j+1) % n_print == 0) cout << endl;
      }
    }
    j = n_bpm - 1;
    inf.getline(line, str_len);
    data[j][Y_][k]= 1e-3*atof(strtok(line, " "));
    if (prt) {
      cout << fixed << setprecision(6) << setw(10) << 1e3*data[j][Y_][k];
      if (k % n_print != 0) cout << endl;
    }
  }

  inf.close();
}


void DFT(double *x, const int n, const int sgn)
{
  int              j, k;
  complex<double>  X[n];

  const complex<double>  I = complex<double>(0e0, 1e0);

  for (j = 0; j <= n/2; j++) {
    X[j] = 0e0;
    for (k = 0; k < n; k++)
      X[j] += x[2*k+1]*exp((double)sgn*I*2e0*M_PI*(double)(k*j)/(double)n);
  }

  for (j = 0; j < n/2; j++) {
    x[2*j+1] = real(X[j]); x[2*(j+1)] = imag(X[j]);
  }
}


void FFT(const int n, const double x[], double A[], double phi[],
	 const int window)
{
  int     i;
  double  *xi;

  xi = dvector(1, 2*n);

  for (i = 0; i < n; i++) {
    switch (window) {
    case 1:
      // Rectangular.
      xi[2*i+1] = x[i];
      break;
    case 2:
      // Sine.
      xi[2*i+1] = sin((double)i/(double)(n-1)*M_PI)*x[i];
      break;
    case 3:
      // Sine^2.
      xi[2*i+1] = sqr(sin((double)i/(double)(n-1)*M_PI))*x[i];
      break;
    default:
      cout << "FFT: not implemented" << endl;
      exit(1);
      break;
    }

    xi[2*(i+1)] = 0e0;
  }

  dfour1(xi, (unsigned long)n, 1);

  for (i = 0; i < n; i++) {
    A[i] = sqrt(sqr(xi[2*i+1])+sqr(xi[2*(i+1)]))*2e0/n;
    phi[i]= -atan2(xi[2*(i+1)], xi[2*i+1]);
  }

  free_dvector(xi, 1, 2*n);
}


void FFT(const int n, const double x[], complex<double> X[], const int window)
{
  int     i;
  double  *xi;

  xi = dvector(1, 2*n);

  for (i = 0; i < n; i++) {
    switch (window) {
    case 1:
      // Rectangular.
      xi[2*i+1] = x[i];
      break;
    case 2:
      // Sine.
      xi[2*i+1] = sin((double)i/(double)(n-1)*M_PI)*x[i];
      break;
    case 3:
      // Sine^2.
      xi[2*i+1] = sqr(sin((double)i/(double)(n-1)*M_PI))*x[i];
      break;
    default:
      cout << "FFT: not implemented" << endl;
      exit(1);
      break;
    }

    xi[2*(i+1)] = 0e0;
  }

  dfour1(xi, (unsigned long)n, 1);

  for (i = 0; i < n; i++)
    X[i] = complex<double>(xi[2*(i+1)], xi[2*i+1]);

  free_dvector(xi, 1, 2*n);
}


void get_ind(const int n, const int k, int &ind1, int &ind3)
{
  // Spectrum for real signal is irror symmetric at k = (0, n/2).
  if (k == 0) {
    ind1 = 1; ind3 = 1;
  } else if (k == n/2) {
    ind1 = n/2-1; ind3 = n/2-1;
  } else {
    ind1 = k - 1; ind3 = k + 1;
  }
}


double get_nu(const int n, const double A[], const int k, const int window)
{
  int     ind, ind1, ind3;
  double  A1, A2, nu;

  get_ind(n, k, ind1, ind3);
  if (A[ind3] > A[ind1]) {
    A1 = A[k]; A2 = A[ind3]; ind = k;
  } else {
    A1 = A[ind1]; A2 = A[k];
    // Special case for 0 frequency.
    ind = (k != 0)? ind1 : -1;
  }
  // Avoid division by zero.
  if (A1+A2 != 0e0)
    switch (window) {
    case 1:
      nu = (ind+A2/(A1+A2))/n;
      break;
    case 2:
      nu = (ind-0.5e0+2e0*A2/(A1+A2))/n;
      break;
    case 3:
      nu = (ind-1e0+3e0*A2/(A1+A2))/n;
      break;
    default:
      cout << "get_nu: not defined" << endl;
      break;
    }
  else
    nu = 0e0;

  return nu;
}


double get_A(const int n, const double A[], const double nu, int k,
	     const int window)
{
  double corr;

  switch (window) {
  case 1:
    corr = Sinc(M_PI*(k-nu*n));
    break;
  case 2:
    corr = (Sinc(M_PI*(k+0.5e0-nu*n))+Sinc(M_PI*(k-0.5e0-nu*n)))/2e0;
    break;
  case 3:
    cout << "get_A: not implemented" << endl;
    exit(1);
    break;
  default:
    cout << "get_A: not defined" << endl;
    break;
  }

  return A[k]/corr;
}


void get_alpha(const int n, const complex<double> X[], const double nu,
		 const int k, double &delta, double &alpha)
{
  // M. Bertocco, C. Offelli, D. Petri "Analysis of Damped Sinusoidal Signals
  // via a Frequency-Domain Interpolation Algorithm" IEEE 43 (2), 245-250
  // (1994).
  int              ind1, ind3, d;
  complex<double>  rho, z;

  const complex<double>  I = complex<double>(0e0, 1e0);

  get_ind(n, k, ind1, ind3);

  if (abs(X[ind3]) > abs(X[ind1])) {
    d = 1; rho = X[ind3]/X[k];
  } else {
    d = -1; rho = X[ind1]/X[k];
  }
  z = (1e0-rho)/(1e0-rho*exp(-I*2e0*M_PI*(double)d/(double)n));
  delta = n*arg(z)/(2e0*M_PI); alpha = n*log(abs(z))/(2e0*M_PI);
}


int get_peak(const int n, const double A[])
{
  int     k, ind1, ind2, ind3;
  double  peak;

  k = 0; peak = 0e0;
  for (ind2 = 0; ind2 <= n/2; ind2++) {
    get_ind(n, ind2, ind1, ind3);
    if ((A[ind2] > peak) && (A[ind1] < A[ind2]) && (A[ind2] > A[ind3])) {
      peak = A[ind2]; k = ind2;
    }
  }

  return k;
}


double get_phi(const int n, const int k, const double nu, const double phi[])
{
  double  phi_nu;

  phi_nu = phi[k] - (n*nu-k)*M_PI;
  if (phi_nu > M_PI)
    phi_nu -= 2.0*M_PI;
  else if (phi_nu < -M_PI)
    phi_nu += 2.0*M_PI;

  return phi_nu;
}


void get_nu(const int n, const double x[], double &nu, double &A_nu,
	    double &phi_nu, double &delta, double &alpha, const int window)
{
  int              k;
  double           A[n], phi[n];
  complex<double>  X[n];

  FFT(n, x, A, phi, window);

  k = get_peak(n, A);
  nu = get_nu(n, A, k, window); A_nu = get_A(n, A, nu, k, window);
  phi_nu = get_phi(n, k, nu, phi);

  // Rectangular window.
  FFT(n, x, X, 1); get_alpha(n, X, nu, k, delta, alpha);
}


void get_nus(const int n_bpm, const int cut, const int n,
	     const int window)
{
  long int  loc;
  int       i, j, k;
  double    x[n], tunes[n_bpm][2], alpha[2], delta[2];
  double    twoJ, beta, As[n_bpm][2], phis[n_bpm][2], nus[n_bpm][2], dnu[2];
  double    tune_sum[2], tune_sum2[2];
  double    alpha_sum[2], alpha_sum2[2];
  double    twoJ_sum[2], twoJ_sum2[2], twoJ_mean[2], twoJ_sigma[2];
  double    phi0[2], phi0_sum[2], phi0_sum2[2], phi0_mean[2], phi0_sigma[2];

  const bool    prt = false;
  const int     sgn[] = {1, -1};
  const double  beta_pinger[] = {6.92e0, 6.76e0};

  for (j = 0; j < 2; j++) {
    tune_sum[j] = 0e0; tune_sum2[j] = 0e0;
    alpha_sum[j] = 0e0; alpha_sum2[j] = 0e0;
    twoJ_sum[j] = 0e0; twoJ_sum2[j] = 0e0;
    phi0_sum[j] = 0e0; phi0_sum2[j] = 0e0;
  }

  printf("\n");
  for (i = 0; i < n_bpm; i++) {
    loc = Elem_GetPos(ElemIndex(bpm_names[i]), 1);

    for (j = 0; j < 2; j++) {
      for (k = cut; k < n+cut; k++)
	x[k-cut] = data[i][j][k];

      rm_mean(n, x);

      get_nu(n, x, tunes[i][j], As[i][j], phis[i][j], delta[j], alpha[j],
	     window);

      if (sgn[j] < 0) phis[i][j] = -phis[i][j];
      if (phis[i][j] < 0e0) phis[i][j] += 2e0*M_PI;
      nus[i][j] = phis[i][j]/(2e0*M_PI);

      tune_sum[j] += tunes[i][j]; tune_sum2[j] += sqr(tunes[i][j]);
      alpha_sum[j] += alpha[j]; alpha_sum2[j] += sqr(alpha[j]);

      twoJ = sqr(As[i][j])/Cell[loc].Beta[j];
      twoJ_sum[j] += twoJ; twoJ_sum2[j] += sqr(twoJ);

      phi0[j] = (nus[i][j]-(Cell[loc].Nu[j]-(int)Cell[loc].Nu[j]))*2e0*M_PI;
      if (phi0[j] < 0e0) phi0[j] += 2e0*M_PI;
      phi0_sum[j] += phi0[j]; phi0_sum2[j] += sqr(phi0[j]);
    }

    // if (prt) printf("[%8.6f, %8.6f]\n", tunes[i][X_], tunes[i][Y_]);
  }

  for (j = 0; j < 2; j++) {
    twoJ_mean[j] = twoJ_sum[j]/n_bpm;
    twoJ_sigma[j] =
      sqrt((n_bpm*twoJ_sum2[j]-sqr(twoJ_sum[j]))/(n_bpm*(n_bpm-1e0)));

    phi0_mean[j] = phi0_sum[j]/n_bpm;
    phi0_sigma[j] =
      sqrt((n_bpm*phi0_sum2[j]-sqr(phi0_sum[j]))/(n_bpm*(n_bpm-1e0)));
  }

  cout << endl;
  cout << scientific << setprecision(3)
       << "twoJ = [" << twoJ_mean[X_] << "+/-" << twoJ_sigma[X_]
       << ", " << twoJ_mean[Y_] << "+/-" << twoJ_sigma[Y_] << "]"
       << fixed 
       << ", phi0 = [" << phi0_mean[X_] << "+/-" << phi0_sigma[X_]
       << ", " << phi0_mean[Y_] << "+/-" << phi0_sigma[Y_] << "]" << endl;
  cout << fixed << setprecision(3)
       << "A0   = [" << 1e3*sqrt(twoJ_mean[X_]*beta_pinger[X_]) << ", "
       << 1e3*sqrt(twoJ_mean[Y_]*beta_pinger[Y_]) << "] mm" << endl;

  // Normalize.
  if (prt) {
    cout << endl;
    cout << " bpm        A               nu            nu (model)" << endl;
  }
  for (i = 0; i < n_bpm; i++) {
    loc = Elem_GetPos(ElemIndex(bpm_names[i]), 1);

    for (j = 0; j < 2; j++) {
      beta = sqr(As[i][j])/twoJ_mean[j];

      nus[i][j] -= phi0_mean[j]/(2e0*M_PI);
      if (nus[i][j] < 0e0) nus[i][j] += 1e0;

      dnu[j] = nus[i][j] - (Cell[loc].Nu[j]-(int)Cell[loc].Nu[j]);
      if (dnu[j] < -0.5e0) dnu[j] += 1e0;
      if (dnu[j] > 0.5e0) dnu[j] -= 1e0;

      betas_sum[i][j] += beta; betas_sum2[i][j] += sqr(beta);
      dnus_sum[i][j] += dnu[j]; dnus_sum2[i][j] += sqr(dnu[j]);
    }

    outf_optics << fixed
		<< setw(4) << i+1
		<< setprecision(3) << setw(8) << Cell[loc].S
		<< setprecision(5) << setw(9) << dnu[X_]
		<< setw(9) << dnu[Y_] << endl;

    if (prt) {
      cout << fixed << setprecision(3)
	   << setw(3) << i+1
	   << "  ["
	   << setw(6) << 1e3*As[i][X_] << ", "
	   << setw(5) << 1e3*As[i][Y_] << "]  [" 
	   << setw(6) << nus[i][X_] << ", "
	   << setw(5) << nus[i][Y_] << "]  ["
	   << setw(6) << Cell[loc].Nu[X_]-(int)Cell[loc].Nu[X_] << ", "
	   << setw(5) << Cell[loc].Nu[Y_]-(int)Cell[loc].Nu[Y_] << "]"
	   << endl;
    }
  }

  for (j = 0; j < 2; j++) {
    tune_mean[j] = tune_sum[j]/n_bpm;
    if (sgn[j] < 0) tune_mean[j] = 1e0 - tune_mean[j];
    tune_sigma[j] =
      sqrt((n_bpm*tune_sum2[j]-sqr(tune_sum[j]))/(n_bpm*(n_bpm-1e0)));

    alpha_mean[j] = alpha_sum[j]/n_bpm;
    alpha_sigma[j] =
      sqrt((n_bpm*alpha_sum2[j]-sqr(alpha_sum[j]))/(n_bpm*(n_bpm-1e0)));
  }

  cout << endl;
  cout << fixed << setprecision(6)
       << "nu    = [" << tune_mean[X_] << "+/-" << tune_sigma[X_]
       << ", " << tune_mean[Y_] << "+/-" << tune_sigma[Y_] << "]" << endl;
  cout << fixed << setprecision(6)
       << "alpha = [" << alpha_mean[X_] << "+/-" << alpha_sigma[X_]
       << ", " << alpha_mean[Y_] << "+/-" << alpha_sigma[Y_] << "]" << endl;

  cout << fixed << setprecision(5)
       << setw(8) << nus[6][X_]-nus[5][X_]
       << setw(8) << nus[6][Y_]-nus[5][Y_] << endl;
}


void get_stats(const int n_bpm)
{
  long int  loc;
  int       j, k;
  double    dbeta[2], dnu[2];
  ofstream  outf;

  const bool    prt = false;
  const double  dbeta_max = 5.0, dnu_max = 0.05;

  if (prt) {
    cout << endl;
    cout << " bpm                A                               "
	 << "nu                              dnu" << endl;
  }
  for (j = 0; j < n_bpm; j++) {
    for (k = 0; k < 2; k++) {
      betas_mean[j][k] = betas_sum[j][k]/n_stats;
      betas_sigma[j][k] = 
	sqrt((n_stats*betas_sum2[j][k]-sqr(betas_sum[j][k]))
	     /(n_stats*(n_stats-1e0)));
      dnus_mean[j][k] = dnus_sum[j][k]/n_stats;
      dnus_sigma[j][k] = 
	sqrt((n_stats*dnus_sum2[j][k]-sqr(dnus_sum[j][k]))
	     /(n_stats*(n_stats-1e0)));
    }

    if (prt)
      cout << fixed << setprecision(3)
	   << setw(3) << j+1 << "  ["
	   << setw(5) << betas_mean[j][X_] << "+/-"
	   << setw(5) << betas_sigma[j][X_] << ", "
	   << setw(4) << betas_mean[j][Y_] << "+/-"
	   << setw(5) << betas_sigma[j][Y_] << "]  [" 
	   << setw(5) << dnus_mean[j][X_] << "+/-"
	   << setw(5) << dnus_sigma[j][X_] << ", "
	   << setw(4) << dnus_mean[j][Y_] << "+/-"
	   << setw(4) << dnus_sigma[j][Y_] << "]" << endl;
  }

  outf.open("tbt.out");

  outf << endl;
  outf << "# bpm  s [m]                 beta [m]                           nu"
       << endl;
  for (j = 0; j < n_bpm; j++) {
    loc = Elem_GetPos(ElemIndex(bpm_names[j]), 1);
    for (k = 0; k < 2; k++) {
      // dbeta[k] = betas_mean[j][k] - Cell[loc].Beta[k];
      dbeta[k] = betas_mean[j][k];
      if (betas_sigma[j][k] > dbeta_max) {
	dbeta[k] = 0e0; betas_sigma[j][k] = 0e0;
      }

      dnu[k] = dnus_mean[j][k] - (Cell[loc].Nu[k]-(int)Cell[loc].Nu[k]);
      if (dnus_sigma[j][k] > dnu_max) {
	cout << endl;
	cout << "BPM # " << j << " excluded, plane = " << k << endl;
	dnu[k] = 0e0; dnus_sigma[j][k] = 0e0;
      }
    }

    outf << fixed << setprecision(3)
	 << setw(4) << j+1
	 << setw(8) << Cell[loc].S
	 << setw(8) << dbeta[X_] << " +/- "
	 << setw(5) << betas_sigma[j][X_]
	 << setw(8) << dbeta[Y_] << " +/- "
	 << setw(5) << betas_sigma[j][Y_]
	 << setw(7) << dnus_mean[j][X_] << " +/- "
	 << setw(5) << dnus_sigma[j][X_]
	 << setw(7) << dnus_mean[j][Y_] << " +/- "
	 << setw(5) << dnus_sigma[j][Y_]
	 << setw(8) << Cell[loc].Beta[X_]
	 << setw(8) << Cell[loc].Beta[Y_] << endl;
  }

  outf.close();
}


void prt_FFT(const int n, const int cut, const double x[], const double y[],
	     const int window)
{
  int       j, k;
  double    A[2][n], phi[2][n], x1[2][n];
  ofstream  outf;

  outf.open("sls.out");

  for (j = cut; j < n+cut; j++) {
    x1[X_][j-cut] = x[j]; x1[Y_][j-cut] = y[j];

    outf << scientific << setprecision(3)
	 << setw(5) << j+1
	 << setw(11) << x[j] << setw(11) << y[j]
	 << endl;
  }

  outf.close();

  for (k = 0; k < 2; k++)
    FFT(n, x1[k], A[k], phi[k], window);

  outf.open("sls_fft.out");

  for (k = 0; k <= n/2; k++)
    outf << scientific << setprecision(3)
	 << setw(5) << k+1
	 << setw(10) << (double)k/(double)n
	 << setw(10) << A[X_][k] << setw(10) << A[Y_][k]
	 << endl;

  outf.close();
}


void Kalman_Filter1(void)
{
  // Kalman filter in Floquet space

  ss_vect<tps>  Pm, Sm;

  // predictor
  Pm = A*P*A_tp + Q;

  //corrector
  Sm = H*Pm*H_tp + R;
  K_ = Pm*H_tp*PInv(Sm, jj);
  P = (Id-K_*H)*Pm;
}


void Kalman_Filter2(ss_vect<double> &x, const ss_vect<double> &z)
{
  // Kalman filter in Floquet space

  ss_vect<tps>  xm;

  // predictor
  xm = (A*x).cst();

  //corrector
  x = (xm+K_*(z-H*xm)).cst();
}


void ss_est(const int n_bpm, const int cut, const int n)
{
  long int         loc1, loc2;
  int              i, j, k;
  double           twoJ_sum[2], twoJ_sum2[2], twoJ_mean[2], twoJ_sigma[2];
  double           dnu[n_bpm_max][2];
  double           phi_sum[2], phi_sum2[2], phi_mean[2], phi_sigma[2], dphi;
  ss_vect<double>  ps[n_bpm], z[n_bpm], dps, ps0;
  ss_vect<tps>     map, Ascr;
  ofstream         outf1, outf2;

  const int     n_DOF = 2, prt_bpm = 6;
  const double  sigma_bpm = 100e-6, sigma_model = 10e-6;


  Id.identity();

  for (k = 0; k < ss_dim; k++)
    jj[k] = (k < 4)? 1 : 0;

  loc1 = Elem_GetPos(ElemIndex(bpm_names[prt_bpm-1]), 1);
	
  cout << endl;
  cout << "Logging at: " << Cell[loc1].Elem.PName << endl;

  for (j = 0; j < n_bpm; j++)
    for (k = 0; k < 2; k++)
      rm_mean(n, data[j][k]);

  map.identity(); putlinmat(2*n_DOF, globval.OneTurnMat, map);
  Ascr.identity(); putlinmat(2*n_DOF, globval.Ascr, Ascr);
  A = PInv(Ascr, jj)*map*Ascr;

  cout << endl;
  cout << "A:" << endl;
  prt_lin_map(n_DOF, A);

  A_tp = tp_S(n_DOF, A);
  H.zero(); H[x_] = Id[x_]; H[y_] = Id[y_];
  // H is singular.
  H_tp = H;

  Q.zero(); R.zero();
  for (j = 0; j < 2*n_DOF; j++) {
    R[j] = sqr(sigma_bpm)*Id[j]; Q[j] = sqr(sigma_model)*Id[j];
  }

  cout << endl;
  cout << "Q:" << endl;
  prt_lin_map(n_DOF, Q);
  cout << endl;
  cout << "R:" << endl;
  prt_lin_map(n_DOF, R);

  P.identity();

  outf1.open("KF.out"); outf2.open("optics.out");

  for (i = 0; i < n_bpm; i++) {
    ps[i].zero();
    for (k = 0; k < 2; k++) {
      ps[i][2*k] = data[i][k][cut];
      ps[i][2*k+1] =
	(data[i][k][cut+1]-A[2*k][2*k]*data[i][k][cut])/A[2*k][2*k+1];
    }
  }

  for (j = cut; j < n+cut; j++) {
    for (k = 0; k < n_DOF; k++) {
      twoJ_sum[k] = 0e0; twoJ_sum2[k] = 0e0;
      phi_sum[k] = 0e0; phi_sum2[k] = 0e0;
    }

    for (i = 0; i < n_bpm; i++) {
      z[i].zero();
      for (k = 0; k < 2; k++) {
	z[i][2*k] = data[i][k][j];
	if (true) z[i][2*k] *= exp(-(j-cut)*2e0*M_PI*alpha_mean[k]/n);
      }

      if (i == prt_bpm-1) Kalman_Filter1();

      Kalman_Filter2(ps[i], z[i]);

      loc1 = Elem_GetPos(ElemIndex(bpm_names[i]), 1);
	
      for (k = 0; k < 2; k++)
	if (true) data[i][k][j] = ps[i][2*k];

      for (k = 0; k < n_DOF; k++) {
	twoJ[i][k][j] = (sqr(ps[i][2*k])+sqr(ps[i][2*k+1]));
	phi[i][k][j] = -atan2(ps[i][2*k+1], ps[i][2*k]);
	if (phi[i][k][j] < 0e0) phi[i][k][j] += 2e0*M_PI;

	phi0[i][k][j] = phi[i][k][j] - (j-cut)*tune_mean[k]*2e0*M_PI;
	phi0[i][k][j] += -(int)(phi0[i][k][j]/(2e0*M_PI))*2e0*M_PI;
	if (phi0[i][k][j] < -M_PI) phi0[i][k][j] += 2e0*M_PI;

	if (i > 0) {
	  loc2 = Elem_GetPos(ElemIndex(bpm_names[i-1]), 1);

	  dphi = phi[i][k][j] - phi[i-1][k][j];
	  if (dphi < 0e0) dphi += 2e0*M_PI;
	  dnu[i][k] = dphi/(2e0*M_PI) - (Cell[loc1].Nu[k]-Cell[loc2].Nu[k]);
	}
      }

      if (i == prt_bpm-1) {
	dps = ps[i] - z[i];
	
	outf1 << setw(4) << j+1;
	for (k = 0; k < 2*n_DOF; k++)
	  outf1 << scientific << setprecision(5) << setw(13) << ps[i][k];
	outf1 << scientific << setprecision(5)
	      << setw(13) << dps[x_] << setw(13) << dps[y_];
	for (k = 0; k < 2*n_DOF; k++)
	  outf1 << scientific << setprecision(5) << setw(12) << P[k][k];
	for (k = 0; k < n_DOF; k++)
	  outf1 << scientific << setprecision(5)
		<< setw(12) << twoJ[i][k][j]
		<< fixed << setw(9) << phi[i][k][j]
		<< setw(9) << phi0[i][k][j];
	outf1 << endl;
      }
    }

    for (k = 0; k < n_DOF; k++) {
      twoJ_mean[k] = twoJ_sum[k]/n_bpm;
      twoJ_sigma[k] = 
	sqrt((n_bpm*twoJ_sum2[k]-sqr(twoJ_sum[k]))/(n_bpm*(n_bpm-1e0)));
      phi_mean[k] = phi_sum[k]/n_bpm;
      phi_sigma[k] = 
	sqrt((n_bpm*phi_sum2[k]-sqr(phi_sum[k]))/(n_bpm*(n_bpm-1e0)));
    }

    outf2 << setw(3) << j+1;
    for (i = 1; i < 10; i++)
      outf2 << scientific << setprecision(5)
	    << setw(12) << twoJ[i][X_][j] << setw(12) << twoJ[i][Y_][j]
	    << fixed << setw(9) << dnu[i][X_] << setw(9) << dnu[i][Y_];
    outf2  << endl;
  }

  outf1.close(); outf2.close();

  loc1 = Elem_GetPos(ElemIndex(bpm_names[5]), 1);
  loc2 = Elem_GetPos(ElemIndex(bpm_names[6]), 1);

  cout << endl;
  cout << fixed << setprecision(5)
       << setw(8) << Cell[loc2].Nu[X_]-Cell[loc1].Nu[X_]
       << setw(8) << Cell[loc2].Nu[Y_]-Cell[loc1].Nu[Y_]
       << endl;

  cout << endl;
  cout << "K:" << endl;
  prt_lin_map(n_DOF, K_);

  cout << endl;
  cout << "P:" << endl;
  prt_lin_map(n_DOF, P);

  cout << endl;
  cout << "Sigmas:";
  for (k = 0; k < 2*n_DOF; k++)
    cout << scientific << setprecision(2) << setw(10) << sqrt(P[k][k]);
  cout << endl;

  cout << endl;
  cout << "Check:" << endl;
//  K_.identity();
  P = A*P*A_tp + Q;
//  P = (Id-K_*H)*P*tp_lin(4, Id-K_*H) + K_*R*tp_lin(4, K_);
  P = (Id-K_*H)*P;
  prt_lin_map(n_DOF, P);
}


void get_b1ob2_dnu(const int n, const ss_vect<double> ps1[],
		   const ss_vect<double> ps2[],
		   double b1ob2[], double dnu[])
{
  // Estimate beta_1/beta_2 and Dnu by tracking data from two adjacent BPMs.

  int     j, k;
  double  x1_sqr, x2_sqr, x1x2;

  cout << endl;
  for (j = 0; j < 2; j++) {
    x1_sqr = 0.0; x2_sqr = 0.0; x1x2 = 0.0;
    for (k = 0; k < n; k++) {
      x1_sqr += sqr(ps1[k][2*j]); x2_sqr += sqr(ps2[k][2*j]);
      x1x2 += ps1[k][2*j]*ps2[k][2*j];
    }

    x1_sqr /= n; x2_sqr /= n; x1x2 /= n;

    b1ob2[j] = x1_sqr/x2_sqr;
    dnu[j] = acos(x1x2/sqrt(x1_sqr*x2_sqr))/(2.0*M_PI);

    cout << scientific << setprecision(3)
	 << "b1ob2 = " << b1ob2[j] << ", dnu = " << dnu[j] << endl;
  }
}


void ss_est2(const int n, const int bpm1, const int bpm2)
{
  long int         loc1, loc2;
  int              j, k;
  double           b1ob2[2], dnu[2];
  ss_vect<double>  ps1[n], ps2[n];

  for (j = 0; j < n; j++)
    for (k = 0; k < 2; k++) {
      ps1[k][j] = data[bpm1-1][k][j]; ps2[k][j] = data[bpm2-1][k][j];
    }

  get_b1ob2_dnu(n, ps1, ps2, b1ob2, dnu);

  loc1 = Elem_GetPos(ElemIndex(bpm_names[bpm1-1]), 1);
  loc2 = Elem_GetPos(ElemIndex(bpm_names[bpm2-1]), 1);

  cout << endl;
  cout << scientific << setprecision(3)
       << "b1ob2 = " << Cell[loc1].Beta[X_]/Cell[loc2].Beta[X_]
       << ", dnu = " << Cell[loc2].Nu[X_]-Cell[loc1].Nu[X_] << endl;
  cout << scientific << setprecision(3)
       << "b1ob2 = " << Cell[loc1].Beta[Y_]/Cell[loc2].Beta[Y_]
       << ", dnu = " << Cell[loc2].Nu[Y_]-Cell[loc1].Nu[Y_] << endl;
}


int main(int argc, char *argv[])
{
  int  j, k, window, cut, n_turn, bpm1, bpm2;


  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  iniranf(1111); setrancut(5.0);

  // sls_ri_f6cwo_20.435_8.737_gset7

  Read_Lattice(argv[1]);

  // globval.Cavity_on = true;

  Ring_GetTwiss(true, 0.0); printglob();

  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);

  switch (2) {
  case 1:
    window = 2; cut = 0*5; n_turn = 2*1024; bpm1 = 6;

    rd_tbt("sls_tbt/tbt_090513_215959.log");

    prt_FFT(n_turn, cut, data[bpm1-1][X_], data[bpm1-1][Y_], window);
    break;
  case 2:
    window = 2; cut = 0; n_turn = 2048;

    for (j = 0; j < n_bpm; j++)
      for (k = 0; k < 2; k++) {
	betas_sum[j][k] = 0e0; betas_sum2[j][k] = 0e0;
	dnus_sum[j][k] = 0e0; dnus_sum2[j][k] = 0e0;
      }

    outf_optics.open("tbt_optics.out");

    if (true) {
      n_stats = 1;
      rd_tbt("sls_tbt/tbt_090513_215619.log");
      get_nus(n_bpm, cut, n_turn, window);

      n_stats += 1;
      rd_tbt("sls_tbt/tbt_090513_215631.log");
      get_nus(n_bpm, cut, n_turn, window);

      n_stats += 1;
      rd_tbt("sls_tbt/tbt_090513_215652.log");
      get_nus(n_bpm, cut, n_turn, window);
    } else {
      n_stats = 1;
      rd_tbt("sls_tbt/tbt_090513_215959.log");
      get_nus(n_bpm, cut, n_turn, window);

      n_stats += 1;
      rd_tbt("sls_tbt/tbt_090513_220010.log");
      get_nus(n_bpm, cut, n_turn, window);

      n_stats += 1;
      rd_tbt("sls_tbt/tbt_090513_220021.log");
      get_nus(n_bpm, cut, n_turn, window);
    }

    outf_optics.close();

    get_stats(n_bpm);
    break;
  case 3:
    window = 2; cut = 0; n_turn = 2048; bpm1 = 3;

    if (true) {
      switch (1) {
      case 1:
	  rd_tbt("sls_tbt/tbt_090513_215619.log");
	break;
      case 2:
	  rd_tbt("sls_tbt/tbt_090513_215631.log");
	break;
      case 3:
	  rd_tbt("sls_tbt/tbt_090513_215652.log");
	break;
      }
    } else {
      switch (1) {
      case 1:
	  rd_tbt("sls_tbt/tbt_090513_215959.log");
	break;
      case 2:
	  rd_tbt("sls_tbt/tbt_090513_220010.log");
	break;
      case 3:
	  rd_tbt("sls_tbt/tbt_090513_220021.log");
	break;
      }
    }

    get_nus(n_bpm, cut, n_turn, window);
    // ss_est(n_bpm, cut, n_turn);

    if (true) {
      for (k = 0; k < 2; k++) {
	rm_mean(n_turn, data[bpm1-1][k]);
	rm_mean(n_turn, twoJ[bpm1-1][k]);
      }

      prt_FFT(n_turn, cut, data[bpm1-1][X_], data[bpm1-1][Y_], window);
//       prt_FFT(n_turn, cut, twoJ[bpm1-1][X_], twoJ[bpm1-1][Y_], window);
    }
    break;
  case 4:
    window = 2; cut = 10; n_turn = 1024;

    for (j = 0; j < n_bpm; j++)
      for (k = 0; k < 2; k++) {
	betas_sum[j][k] = 0e0; betas_sum2[j][k] = 0e0;
	dnus_sum[j][k] = 0e0; dnus_sum2[j][k] = 0e0;
      }

    outf_optics.open("tbt_optics.out");

    n_stats = 1;
    rd_tbt("sls_tbt/tbt_090513_215619.log");
    ss_est(n_bpm, cut, n_turn);
    get_nus(n_bpm, cut, n_turn, window);
    
    n_stats += 1;
    rd_tbt("sls_tbt/tbt_090513_215631.log");
    ss_est(n_bpm, cut, n_turn);
    get_nus(n_bpm, cut, n_turn, window);
    
    n_stats += 1;
    rd_tbt("sls_tbt/tbt_090513_215652.log");
    ss_est(n_bpm, cut, n_turn);
    get_nus(n_bpm, cut, n_turn, window);
    
    outf_optics.close();

    get_stats(n_bpm);
    break;
  case 5:
    window = 2; cut = 0; n_turn = 2048; bpm1 = 3;

    rd_tbt("sls_tbt/tbt_090513_215619.log");
//     rd_tbt("sls_tbt/tbt_090513_215959.log");
    get_nus(n_bpm, cut, n_turn, window);
    // ss_est(n_bpm, cut, n_turn);

    for (k = 0; k < 2; k++)
      rm_mean(n_turn, twoJ[bpm1-1][k]);

//     prt_FFT(n_turn, cut, twoJ[bpm1-1][X_], twoJ[bpm1-1][Y_], window);
    prt_FFT(n_turn, cut, data[bpm1-1][X_], data[bpm1-1][Y_], window);
    break;
  case 6:
    n_turn = 2048; bpm1 = 5; bpm2 = 6;

    rd_tbt("sls_tbt/tbt_090513_215619.log");

    ss_est2(n_turn, bpm1, bpm2);
    break;
  }
}
