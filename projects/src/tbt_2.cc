#define NO 1

#include "tracy_lib.h"

int no_tps = NO;

const int  n_bpm_max = 100;

char    bpm_names[n_bpm_max][max_str];
int     n_bpm, n_turn_, n_stats, jj[ss_dim];
double  data[n_bpm_max][2][2048];
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


void prt_lin_map(const int n_DOF, const ss_vect<tps> &map)
{
  int  i, j;

  cout << endl;
  for (i = 1; i <= 2*n_DOF; i++) {
    for (j = 1; j <= 2*n_DOF; j++)
      if (true) 
	cout << scientific << setprecision(5)
	     << setw(13) << getmat(map, i, j);
      else
	cout << scientific << setprecision(16)
	     << setw(24) << getmat(map, i, j);
    cout << endl;
  }
}


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
	    double &phi_nu, const int window)
{
  int              k;
  double           A[n], phi[n];
  complex<double>  X[n];

  FFT(n, x, A, phi, window);

  k = get_peak(n, A);
  nu = get_nu(n, A, k, window); A_nu = get_A(n, A, nu, k, window);
  phi_nu = get_phi(n, k, nu, phi);
}


void get_nus(const int n_bpm, const int cut, const int n,
	     const int window)
{
  long int  loc;
  int       i, j, k;
  double    x[n], tunes[n_bpm][2];
  double    twoJ, beta, As[n_bpm][2], phis[n_bpm][2], nus[n_bpm][2], dnu[2];
  double    tune_sum[2], tune_sum2[2];
  double    twoJ_sum[2], twoJ_sum2[2], twoJ_mean[2], twoJ_sigma[2];
  double    phi0[2], phi0_sum[2], phi0_sum2[2], phi0_mean[2], phi0_sigma[2];

  const bool    prt = false;
  const int     sgn[] = {1, -1};
  const double  beta_pinger[] = {6.92e0, 6.76e0};

  for (j = 0; j < 2; j++) {
    tune_sum[j] = 0e0; tune_sum2[j] = 0e0;
    twoJ_sum[j] = 0e0; twoJ_sum2[j] = 0e0;
    phi0_sum[j] = 0e0; phi0_sum2[j] = 0e0;
  }

  for (i = 0; i < n_bpm; i++) {
    loc = Elem_GetPos(ElemIndex(bpm_names[i]), 1);

    for (j = 0; j < 2; j++) {
      for (k = cut; k < n+cut; k++)
	x[k-cut] = data[i][j][k];

      rm_mean(n, x);

      get_nu(n, x, tunes[i][j], As[i][j], phis[i][j], window);

      if (sgn[j] < 0) phis[i][j] = -phis[i][j];
      if (phis[i][j] < 0e0) phis[i][j] += 2e0*M_PI;
      nus[i][j] = phis[i][j]/(2e0*M_PI);

      tune_sum[j] += tunes[i][j]; tune_sum2[j] += sqr(tunes[i][j]);

      twoJ = sqr(As[i][j])/Cell[loc].Beta[j];
      twoJ_sum[j] += twoJ; twoJ_sum2[j] += sqr(twoJ);

      phi0[j] = (nus[i][j]-(Cell[loc].Nu[j]-(int)Cell[loc].Nu[j]))*2e0*M_PI;
      if (phi0[j] < 0e0) phi0[j] += 2e0*M_PI;
      phi0_sum[j] += phi0[j]; phi0_sum2[j] += sqr(phi0[j]);
    }
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
    cout << " bpm       A              nu" << endl;
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
	   << setw(6) << As[i][X_] << ", "
	   << setw(5) << As[i][Y_] << "]  [" 
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
  }

  cout << endl;
  cout << fixed << setprecision(6)
       << "nu    = [" << tune_mean[X_] << "+/-" << tune_sigma[X_]
       << ", " << tune_mean[Y_] << "+/-" << tune_sigma[Y_] << "]" << endl;

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
      dbeta[k] = betas_mean[j][k] - Cell[loc].Beta[k];
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
	 << setw(5) << dnus_sigma[j][Y_] << endl;
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


int main(int argc, char *argv[])
{
  int  j, k, window, cut, n_turn;


  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  iniranf(1111); setrancut(5.0);

  Read_Lattice(argv[1]);

  Ring_GetTwiss(true, 0.0); printglob();

  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);

  window = 2; cut = 0; n_turn = 2048/2;

  for (j = 0; j < n_bpm; j++)
    for (k = 0; k < 2; k++) {
      betas_sum[j][k] = 0e0; betas_sum2[j][k] = 0e0;
      dnus_sum[j][k] = 0e0; dnus_sum2[j][k] = 0e0;
    }

  outf_optics.open("tbt_optics.out");

  n_stats = 1;
  rd_tbt("SLS_TBT/tbt_090513_215619.log");
  get_nus(n_bpm, cut, n_turn, window);

  n_stats += 1;
  rd_tbt("SLS_TBT/tbt_090513_215631.log");
  get_nus(n_bpm, cut, n_turn, window);

  n_stats += 1;
  rd_tbt("SLS_TBT/tbt_090513_215652.log");
  get_nus(n_bpm, cut, n_turn, window);

  outf_optics.close();

  get_stats(n_bpm);
}
