#include <stdio.h>
#include <stddef.h>
#include <string.h>

#include <cstdlib>
#include <cfloat>
#include <cctype>
#include <cmath>
#include <complex>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <istream>
#include <string>
#include <vector>

using namespace std;

extern "C" {
    double *dvector(long nl, long nh);
    void free_dvector(double *v, long nl, long nh);
    void dfour1(double data[], unsigned long nn, int isign);
}

#define sqr(x) ((x)*(x))

enum spatial_index { X_ = 0, Y_ = 1, Z_ = 2 };

const int ss_dim = 6;

typedef double ss_vect[ss_dim];


struct lin_opt_type {
  public:
    std::vector < std::string > name;
    std::vector < int >loc;
     std::vector < double >s, alpha[2], beta[2], nu[2], eta[2], etap[2];

    void rd_data(const string & file_name);
};

struct bpm_data_type {
  public:
    int n_bpm, n_turn;
     std::vector < std::string > name;
     std::vector < int >loc;
     std::vector < std::vector < double > >data[2];

    void rd_bpm_names(ifstream & inf, const lin_opt_type & lin_opt);
    void rd_bpm_data(const int plane, ifstream & inf);
    void rd_tbt(const char *file_name, lin_opt_type & lin_opt);
};

struct est_lin_opt_type {
  public:
    int n_stats;
    double alpha_mean[2], alpha_sigma[2], tune_mean[2], tune_sigma[2];
     std::vector < double >
	beta[2], beta_sum[2], beta_sum2[2], beta_mean[2], beta_sigma[2],
	nu[2], dnu_sum[2], dnu_sum2[2], dnu_mean[2], dnu_sigma[2];

     std::vector < std::vector < double > >twoJ[2], phi[2], phi0[2];

    void zero(const int n);
    void get_stats(const bpm_data_type & bpm_data,
		   const lin_opt_type & lin_opt);
};


void get_line(ifstream & inf, stringstream & str)
{
    string line;

    getline(inf, line);
    str.clear();
    str.str("");
    str << line;
}


int get_loc(const string & name, const lin_opt_type & lin_opt)
{
    int k;

    k = -1;
    do {
	k++;
    } while ((k < (int) lin_opt.loc.size())
	     && (name != lin_opt.name[k]));
    return k;
}


void lin_opt_type::rd_data(const string & file_name)
{
    char comma;
    string name, line;
    int n, k;
    double s, type, alpha[2], beta[2], nu[2], eta[2], etap[2];
    stringstream str;
    ifstream inf;

    const bool prt = false;

    inf.open(file_name.c_str());

    if (prt)
	printf("\n");
    while (getline(inf, line)) {
	// Skip comments; lines starting with "#".
	if (line.compare(0, 1, "#") != 0) {
	    str.clear();
	    str.str("");
	    str << line;
	    if (prt)
		cout << str.str() << "\n";
	    str >> n >> comma >> name >> s >> comma >> type >>
		comma >> alpha[X_] >> comma >> beta[X_] >> comma >> nu[X_]
		>> comma >> eta[X_] >> comma >> etap[X_] >> comma >>
		alpha[Y_] >> comma >> beta[Y_] >> comma >> nu[Y_] >> comma
		>> eta[Y_] >> comma >> etap[Y_];
	    name.erase(name.end() - 1);
	    if (prt)
		printf("%4d, %-15s, %9.5f, %4.1f,"
		       " %9.5f, %8.5f, %8.5f, %8.5f, %8.5f,"
		       " %9.5f, %8.5f, %8.5f, %8.5f, %8.5f\n",
		       n, name.c_str(), s, type,
		       alpha[X_], beta[X_], nu[X_], eta[X_], etap[X_],
		       alpha[Y_], beta[Y_], nu[Y_], eta[Y_], etap[Y_]);
	    this->loc.push_back(n);
	    this->s.push_back(s);
	    this->name.push_back(name);
	    for (k = 0; k < 2; k++) {
		this->alpha[k].push_back(alpha[k]);
		this->beta[k].push_back(beta[k]);
		this->nu[k].push_back(nu[k]);
		this->eta[k].push_back(eta[k]);
		this->etap[k].push_back(etap[k]);
	    }
	}
    }

    inf.close();
}


void get_bpm_name(string & name)
{
    int k;

    k = 0;
    do {
	name[k] = (name[k] != '-') ? tolower(name[k]) : '_';
	k++;
    } while (k < (int) name.size());
}


void bpm_data_type::rd_bpm_names(ifstream & inf, const lin_opt_type & lin_opt)
{
    string line, name1;
    int j, k;
    stringstream str;

    const bool prt = false;
    const int n_print = 8;

    // Skip 1st line.
    getline(inf, line);
    // Read no of BPMs and no of turns.
    get_line(inf, str);
    str >> n_bpm >> n_turn;
    // Skip 3rd line.
    getline(inf, line);
    cout << "\nno of BPMs = " << n_bpm
	<< ", no of turns = " << n_turn << "\n";

    for (k = 0; k < 2; k++) {
	data[k].resize(n_bpm);
	for (j = 0; j < n_bpm; j++)
	    data[k][j].resize(n_turn);
    }

    if (prt)
	cout << "\n";
    for (j = 0; j < n_bpm; j++) {
	get_line(inf, str);
	str >> name1;
	get_bpm_name(name1);
	name.push_back(name1);
	loc.push_back(get_loc(name1, lin_opt));
	if (prt) {
	    cout << " " << name1[j];
	    if ((j + 1) % n_print == 0)
		cout << "\n";
	}
    }
    if (prt && (n_bpm % n_print != 0))
	cout << "\n";
}


void bpm_data_type::rd_bpm_data(const int plane, ifstream & inf)
{
    string line;
    int j, k;
    stringstream sstr;

    const bool prt = false;
    const int n_print = 8;

    if (prt)
	cout << "\n";
    // Skip line.
    get_line(inf, sstr);
    for (k = 0; k < n_turn; k++) {
	if (prt)
	    cout << "\n";
	get_line(inf, sstr);
	for (j = 0; j < n_bpm - 1; j++) {
	    sstr >> data[plane][j][k];
	    data[plane][j][k] *= 1e-3;
	    if (prt) {
		cout << fixed << setprecision(6)
		    << setw(10) << 1e3 * data[plane][j][k];
		if ((j + 1) % n_print == 0)
		    cout << "\n";
	    }
	}
	// Read data for last BPM.
	j = n_bpm - 1;
	get_line(inf, sstr);
	sstr >> data[plane][j][k];
	data[plane][j][k] *= 1e-3;
	if (prt) {
	    cout << fixed << setprecision(6)
		<< setw(10) << 1e3 * data[plane][j][k];
	    if (k % n_print == 0)
		cout << "\n";
	}
	if (prt && k % n_print != 0)
	  cout << "\n";
    }

}


void bpm_data_type::rd_tbt(const char *file_name, lin_opt_type & lin_opt)
{
    string line;
    stringstream sstr;
    ifstream inf;

    inf.open(file_name);

    this->rd_bpm_names(inf, lin_opt);
    this->rd_bpm_data(0, inf);
    this->rd_bpm_data(1, inf);

    inf.close();
}


void est_lin_opt_type::zero(const int n)
{
    int j, k;

    for (k = 0; k < 2; k++) {
	beta[k].resize(n);
	beta_sum[k].resize(n);
	beta_sum2[k].resize(n);
	beta_mean[k].resize(n);
	beta_sigma[k].resize(n);
	nu[k].resize(n);
	dnu_sum[k].resize(n);
	dnu_sum2[k].resize(n);
	dnu_mean[k].resize(n);
	dnu_sigma[k].resize(n);
    }

    for (j = 0; j < (int) beta_sum[X_].size(); j++)
	for (k = 0; k < 2; k++) {
	    beta_sum[k][j] = 0e0;
	    beta_sum2[k][j] = 0e0;
	    dnu_sum[k][j] = 0e0;
	    dnu_sum2[k][j] = 0e0;
	}
}


void DFT(double *x, const int n, const int sgn)
{
    int j, k;
    complex < double >X[n];

    const complex < double >I = complex < double >(0e0, 1e0);

    for (j = 0; j <= n / 2; j++) {
	X[j] = 0e0;
	for (k = 0; k < n; k++)
	    X[j] +=
	      x[2 * k + 1]* exp((double) sgn * I * 2e0 * M_PI
				* (double) (k * j) / (double) n);
    }
    for (j = 0; j < n / 2; j++) {
	x[2 * j + 1] = real(X[j]);
	x[2 * (j + 1)] = imag(X[j]);
    }
}


void
FFT(const int n, const double x[], double A[], double phi[],
    const int window)
{
    int i;
    double *xi;

    xi = dvector(1, 2 * n);

    for (i = 0; i < n; i++) {
	switch (window) {
	case 1:
	    // Rectangular.
	    xi[2 * i + 1] = x[i];
	    break;
	case 2:
	    // Sine.
	    xi[2 * i + 1] =
		sin((double) i / (double) (n - 1) * M_PI) * x[i];
	    break;
	case 3:
	    // Sine^2.
	    xi[2 * i + 1] =
		sqr(sin((double) i / (double) (n - 1) * M_PI)) * x[i];
	    break;
	default:
	    cout << "FFT: not implemented\n";
	    exit(1);
	    break;
	}

	xi[2 * (i + 1)] = 0e0;
    }

    dfour1(xi, (unsigned long) n, 1);

    for (i = 0; i < n; i++) {
	A[i] = sqrt(sqr(xi[2 * i + 1]) + sqr(xi[2 * (i + 1)])) * 2e0 / n;
	phi[i] = -atan2(xi[2 * (i + 1)], xi[2 * i + 1]);
    }

    free_dvector(xi, 1, 2 * n);
}


void
FFT(const int n, const double x[], complex < double >X[], const int window)
{
    int i;
    double *xi;

    xi = dvector(1, 2 * n);

    for (i = 0; i < n; i++) {
	switch (window) {
	case 1:
	    // Rectangular.
	    xi[2 * i + 1] = x[i];
	    break;
	case 2:
	    // Sine.
	    xi[2 * i + 1] =
		sin((double) i / (double) (n - 1) * M_PI) * x[i];
	    break;
	case 3:
	    // Sine^2.
	    xi[2 * i + 1] =
		sqr(sin((double) i / (double) (n - 1) * M_PI)) * x[i];
	    break;
	default:
	    cout << "FFT: not implemented" << "\n";
	    exit(1);
	    break;
	}

	xi[2 * (i + 1)] = 0e0;
    }

    dfour1(xi, (unsigned long) n, 1);

    for (i = 0; i < n; i++)
	X[i] = complex < double >(xi[2 * (i + 1)], xi[2 * i + 1]);

    free_dvector(xi, 1, 2 * n);
}


void get_ind(const int n, const int k, int &ind1, int &ind3)
{
    // Spectrum for real signal is mirror symmetric at k = (0, n/2).
    if (k == 0) {
	ind1 = 1;
	ind3 = 1;
    } else if (k == n / 2) {
	ind1 = n / 2 - 1;
	ind3 = n / 2 - 1;
    } else {
	ind1 = k - 1;
	ind3 = k + 1;
    }
}


double get_nu(const int n, const double A[], const int k, const int window)
{
    int ind, ind1, ind3;
    double A1, A2, nu = 0e0;

    get_ind(n, k, ind1, ind3);
    if (A[ind3] > A[ind1]) {
	A1 = A[k];
	A2 = A[ind3];
	ind = k;
    } else {
	A1 = A[ind1];
	A2 = A[k];
	// Special case for 0 frequency.
	ind = (k != 0) ? ind1 : -1;
    }
    // Avoid division by zero.
    if (A1 + A2 != 0e0)
	switch (window) {
	case 1:
	    nu = (ind + A2 / (A1 + A2)) / n;
	    break;
	case 2:
	    nu = (ind - 0.5e0 + 2e0 * A2 / (A1 + A2)) / n;
	    break;
	case 3:
	    nu = (ind - 1e0 + 3e0 * A2 / (A1 + A2)) / n;
	    break;
	default:
	    cout << "get_nu: not defined\n";
	    break;
    } else
	nu = 0e0;

    return nu;
}


double sinc(double omega)
{
    return (omega != 0e0) ? sin(omega) / omega : 1e0;
}


double
get_A(const int n, const double A[], const double nu, int k,
      const int window)
{
    double corr = 0e0;

    switch (window) {
    case 1:
	corr = sinc(M_PI * (k - nu * n));
	break;
    case 2:
	corr =
	    (sinc(M_PI * (k + 0.5e0 - nu * n)) +
	     sinc(M_PI * (k - 0.5e0 - nu * n))) / 2e0;
	break;
    case 3:
	cout << "get_A: not implemented\n";
	exit(1);
	break;
    default:
	cout << "get_A: not defined\n";
	break;
    }

    return A[k] / corr;
}


void
get_alpha(const int n, const complex < double >X[], const double nu,
	  const int k, double &delta, double &alpha)
{
    // M. Bertocco, C. Offelli, D. Petri "Analysis of Damped Sinusoidal
    // Signals
    // via a Frequency-Domain Interpolation Algorithm" IEEE 43 (2),
    // 245-250
    // (1994).
    int ind1, ind3, d;
    complex < double >rho, z;

    const complex < double >I = complex < double >(0e0, 1e0);

    get_ind(n, k, ind1, ind3);

    if (abs(X[ind3]) > abs(X[ind1])) {
	d = 1;
	rho = X[ind3] / X[k];
    } else {
	d = -1;
	rho = X[ind1] / X[k];
    }
    z = (1e0 - rho) / (1e0 -
		       rho * exp(-I * 2e0 * M_PI * (double) d /
				 (double) n));
    delta = n * arg(z) / (2e0 * M_PI);
    alpha = n * log(abs(z)) / (2e0 * M_PI);
}


int get_peak(const int n, const double A[])
{
    int k, ind1, ind2, ind3;
    double peak;

    k = 0;
    peak = 0e0;
    for (ind2 = 0; ind2 <= n / 2; ind2++) {
	get_ind(n, ind2, ind1, ind3);
	if ((A[ind2] > peak) && (A[ind1] < A[ind2]) && (A[ind2] > A[ind3])) {
	    peak = A[ind2];
	    k = ind2;
	}
   }

    return k;
}


double
get_phi(const int n, const int k, const double nu, const double phi[])
{
    double phi_nu;

    phi_nu = phi[k] - (n * nu - k) * M_PI;
    if (phi_nu > M_PI)
	phi_nu -= 2.0 * M_PI;
    else if (phi_nu < -M_PI)
	phi_nu += 2.0 * M_PI;

    return phi_nu;
}


void
get_nu(const int n, const double x[], double &nu, double &A_nu,
       double &phi_nu, double &delta, double &alpha, const int window)
{
    int k;
    double A[n], phi[n];
    complex < double >X[n];

    FFT(n, x, A, phi, window);

    k = get_peak(n, A);
    nu = get_nu(n, A, k, window);
    A_nu = get_A(n, A, nu, k, window);
    phi_nu = get_phi(n, k, nu, phi);

    // Rectangular window.
    FFT(n, x, X, 1);
    get_alpha(n, X, nu, k, delta, alpha);
}


void rm_mean(long int n, double x[])
{
    long int i;
    double mean;

    mean = 0e0;
    for (i = 0; i < n; i++)
	mean += x[i];
    mean /= n;
    for (i = 0; i < n; i++)
	x[i] -= mean;
}


void get_nus(ofstream & outf, const int cut, const int n, const int window,
	     const bpm_data_type & bpm_data, const lin_opt_type & lin_opt,
	     est_lin_opt_type & est_lin_opt)
{
    long int loc;
    int i, j, k;
    double x[n], tunes[bpm_data.n_bpm][2], alpha[2], delta[2];
    double twoJ, beta, As[bpm_data.n_bpm][2], phis[bpm_data.n_bpm][2];
    double nus[bpm_data.n_bpm][2], dnu[2];
    double tune_sum[2], tune_sum2[2];
    double alpha_sum[2], alpha_sum2[2];
    double twoJ_sum[2], twoJ_sum2[2], twoJ_mean[2], twoJ_sigma[2];
    double phi0[2], phi0_sum[2], phi0_sum2[2], phi0_mean[2], phi0_sigma[2];

    const bool prt = false;
    const int sgn[] = { 1, -1 };
    const double beta_pinger[] = { 6.92e0, 6.76e0 };

    for (j = 0; j < 2; j++) {
	tune_sum[j] = 0e0;
	tune_sum2[j] = 0e0;
	alpha_sum[j] = 0e0;
	alpha_sum2[j] = 0e0;
	twoJ_sum[j] = 0e0;
	twoJ_sum2[j] = 0e0;
	phi0_sum[j] = 0e0;
	phi0_sum2[j] = 0e0;
    }

    printf("\n");
    for (i = 0; i < bpm_data.n_bpm; i++) {
	loc = bpm_data.loc[i];

	for (j = 0; j < 2; j++) {
	    for (k = cut; k < n + cut; k++)
		x[k - cut] = bpm_data.data[j][i][k];

	    rm_mean(n, x);

	    get_nu(n, x, tunes[i][j], As[i][j], phis[i][j], delta[j],
		   alpha[j], window);

	    if (sgn[j] < 0)
		phis[i][j] = -phis[i][j];
	    if (phis[i][j] < 0e0)
		phis[i][j] += 2e0 * M_PI;
	    nus[i][j] = phis[i][j] / (2e0 * M_PI);

	    tune_sum[j] += tunes[i][j];
	    tune_sum2[j] += sqr(tunes[i][j]);
	    alpha_sum[j] += alpha[j];
	    alpha_sum2[j] += sqr(alpha[j]);

	    twoJ = sqr(As[i][j]) / lin_opt.beta[j][loc];
	    twoJ_sum[j] += twoJ;
	    twoJ_sum2[j] += sqr(twoJ);

	    phi0[j] =
		(nus[i][j] -
		 (lin_opt.nu[j][loc] -
		  (int) lin_opt.nu[j][loc])) * 2e0 * M_PI;
	    if (phi0[j] < 0e0)
		phi0[j] += 2e0 * M_PI;
	    phi0_sum[j] += phi0[j];
	    phi0_sum2[j] += sqr(phi0[j]);
	}

	// if (prt) printf("[%8.6f, %8.6f]\n", tunes[i][X_],
	// tunes[i][Y_]);
    }

    for (j = 0; j < 2; j++) {
	twoJ_mean[j] = twoJ_sum[j] / bpm_data.n_bpm;
	twoJ_sigma[j] =
	    sqrt((bpm_data.n_bpm * twoJ_sum2[j] - sqr(twoJ_sum[j]))
		 / (bpm_data.n_bpm * (bpm_data.n_bpm - 1e0)));

	phi0_mean[j] = phi0_sum[j] / bpm_data.n_bpm;
	phi0_sigma[j] =
	    sqrt((bpm_data.n_bpm * phi0_sum2[j] - sqr(phi0_sum[j]))
		 / (bpm_data.n_bpm * (bpm_data.n_bpm - 1e0)));
    }

    cout << scientific << setprecision(3)
	<< "\ntwoJ  = [" << twoJ_mean[X_] << "+/-" << twoJ_sigma[X_]
	<< ", " << twoJ_mean[Y_] << "+/-" << twoJ_sigma[Y_] << "]"
	<< fixed
	<< ", phi0 = [" << phi0_mean[X_] << "+/-" << phi0_sigma[X_]
	<< ", " << phi0_mean[Y_] << "+/-" << phi0_sigma[Y_] << "]\n";
    cout << fixed << setprecision(3)
	<< "A0    = [" << 1e3 * sqrt(twoJ_mean[X_] *
				    beta_pinger[X_]) << ", " << 1e3 *
	sqrt(twoJ_mean[Y_] * beta_pinger[Y_]) << "] mm\n";

    // Normalize.
    if (prt)
	cout << "\n bpm        A               nu            nu (model)\n";
    for (i = 0; i < bpm_data.n_bpm; i++) {
	loc = bpm_data.loc[i];

	for (j = 0; j < 2; j++) {
	    beta = sqr(As[i][j]) / twoJ_mean[j];

	    nus[i][j] -= phi0_mean[j] / (2e0 * M_PI);
	    if (nus[i][j] < 0e0)
		nus[i][j] += 1e0;

	    dnu[j] =
		nus[i][j] - (lin_opt.nu[j][loc] -
			     (int) lin_opt.nu[j][loc]);
	    if (dnu[j] < -0.5e0)
		dnu[j] += 1e0;
	    if (dnu[j] > 0.5e0)
		dnu[j] -= 1e0;

	    est_lin_opt.beta_sum[j][i] += beta;
	    est_lin_opt.beta_sum2[j][i] += sqr(beta);
	    est_lin_opt.dnu_sum[j][i] += dnu[j];
	    est_lin_opt.dnu_sum2[j][i] += sqr(dnu[j]);
	}

	outf << fixed << setw(4) << i + 1
	    << setprecision(3) << setw(8) << lin_opt.s[loc]
	    << setprecision(5) << setw(9) << dnu[X_]
	    << setw(9) << dnu[Y_] << "\n";

	if (prt) {
	    cout << fixed << setprecision(3)
		<< setw(3) << i + 1
		<< "  ["
		<< setw(6) << 1e3 * As[i][X_] << ", "
		<< setw(5) << 1e3 * As[i][Y_] << "]  ["
		<< setw(6) << nus[i][X_] << ", "
		<< setw(5) << nus[i][Y_] << "]  [" << setw(6)
		<< lin_opt.nu[X_][loc] - (int) lin_opt.nu[X_][loc] << ", "
		<< setw(5)
		<< lin_opt.nu[Y_][loc] -
		(int) lin_opt.nu[Y_][loc] << "]\n";
	}
    }

    for (j = 0; j < 2; j++) {
	est_lin_opt.tune_mean[j] = tune_sum[j] / bpm_data.n_bpm;
	if (sgn[j] < 0)
	    est_lin_opt.tune_mean[j] = 1e0 - est_lin_opt.tune_mean[j];
	est_lin_opt.tune_sigma[j] =
	    sqrt((bpm_data.n_bpm * tune_sum2[j] - sqr(tune_sum[j]))
		 / (bpm_data.n_bpm * (bpm_data.n_bpm - 1e0)));

	est_lin_opt.alpha_mean[j] = alpha_sum[j] / bpm_data.n_bpm;
	est_lin_opt.alpha_sigma[j] =
	    sqrt((bpm_data.n_bpm * alpha_sum2[j] - sqr(alpha_sum[j]))
		 / (bpm_data.n_bpm * (bpm_data.n_bpm - 1e0)));
    }

    cout << fixed << setprecision(6)
	<< "\nnu    = [" << est_lin_opt.tune_mean[X_] << "+/-"
	<< est_lin_opt.tune_sigma[X_]
	<< ", " << est_lin_opt.tune_mean[Y_] << "+/-"
	<< est_lin_opt.tune_sigma[Y_] << "]\n";
    cout << fixed << setprecision(6)
	<< "alpha = [" << est_lin_opt.alpha_mean[X_] << "+/-"
	<< est_lin_opt.alpha_sigma[X_]
	<< ", " << est_lin_opt.alpha_mean[Y_] << "+/-"
	<< est_lin_opt.alpha_sigma[Y_] << "]\n";

    cout << fixed << setprecision(5)
	<< setw(8) << nus[6][X_] - nus[5][X_]
	<< setw(8) << nus[6][Y_] - nus[5][Y_] << "\n";
}


void get_m_s(const int n, const double sum, const double sum2,
	     double &mean, double &sigma)
{
    mean = sum / n;
    sigma = sqrt((n * sum2 - sqr(sum)) / (n * (n - 1e0)));
}


void est_lin_opt_type::get_stats(const bpm_data_type & bpm_data,
				 const lin_opt_type & lin_opt)
{
    long int loc;
    int j, k;
    double dbeta[2], dnu[2];
    ofstream outf;

    const bool prt = false;
    const double dbeta_max = 5.0, dnu_max = 0.05;

    if (prt)
	cout << "\n bpm                A                               "
	    << "    dnu\n";
    for (j = 0; j < bpm_data.n_bpm; j++) {
	for (k = 0; k < 2; k++) {
	    get_m_s(n_stats, beta_sum[k][j], beta_sum2[k][j],
		    beta_mean[k][j], beta_sigma[k][j]);
	    get_m_s(n_stats, dnu_sum[k][j], dnu_sum2[k][j], dnu_mean[k][j],
		    dnu_sigma[k][j]);
	}

	if (prt)
	    cout << fixed << setprecision(3)
		<< setw(3) << j + 1 << "  ["
		<< setw(6) << beta_mean[X_][j] << "+/-"
		<< setw(6) << beta_sigma[X_][j] << ", "
		<< setw(6) << beta_mean[Y_][j] << "+/-"
		<< setw(6) << beta_sigma[Y_][j] << "]  ["
		<< setw(6) << dnu_mean[X_][j] << "+/-"
		<< setw(5) << dnu_sigma[X_][j] << ", "
		<< setw(6) << dnu_mean[Y_][j] << "+/-"
		<< setw(5) << dnu_sigma[Y_][j] << "]\n";
    }

    outf.open("tbt.out");

    outf << "\n# bpm  s [m]                 beta [m]"
	"                           nu\n";
    for (j = 0; j < bpm_data.n_bpm; j++) {
	loc = bpm_data.loc[j];
	for (k = 0; k < 2; k++) {
	    // dbeta[k] = beta_mean[k][j] - lin_opt.beta[k][loc];
	    dbeta[k] = beta_mean[k][j];
	    if (beta_sigma[k][j] > dbeta_max) {
		dbeta[k] = 0e0;
		beta_sigma[k][j] = 0e0;
	    }

	    dnu[k] =
		dnu_mean[k][j] - (lin_opt.nu[k][loc] -
				  (int) lin_opt.nu[k][loc]);
	    if (dnu_sigma[k][j] > dnu_max) {
		cout << "\nBPM # " << j << " excluded, plane = " << k
		    << "\n";
		dnu[k] = 0e0;
		dnu_sigma[k][j] = 0e0;
	    }
	}

	outf << fixed << setprecision(3)
	    << setw(4) << j + 1 << setw(8) << lin_opt.s[loc]
	    << setw(8) << dbeta[X_] << " +/- "
	    << setw(5) << beta_sigma[X_][j]
	    << setw(8) << dbeta[Y_] << " +/- "
	    << setw(5) << beta_sigma[Y_][j]
	    << setw(7) << dnu_mean[X_][j] << " +/- "
	    << setw(5) << dnu_sigma[X_][j]
	    << setw(7) << dnu_mean[Y_][j] << " +/- "
	    << setw(5) << dnu_sigma[Y_][j]
	    << setw(8) << lin_opt.beta[X_][loc]
	    << setw(8) << lin_opt.beta[Y_][loc] << "\n";
    }

    outf.close();
}


void
prt_FFT(const int n, const int cut, const std::vector < double >&x,
	const std::vector < double >&y, const int window)
{
    int j, k;
    double A[2][n], phi[2][n], x1[2][n];
    ofstream outf;

    outf.open("sls.out");

    for (j = cut; j < n + cut; j++) {
	x1[X_][j - cut] = x[j];
	x1[Y_][j - cut] = y[j];

	outf << scientific << setprecision(3)
	    << setw(5) << j +
	    1 << setw(11) << x[j] << setw(11) << y[j] << "\n";
    }

    outf.close();

    for (k = 0; k < 2; k++)
	FFT(n, x1[k], A[k], phi[k], window);

    outf.open("sls_fft.out");

    for (k = 0; k <= n / 2; k++)
	outf << scientific << setprecision(3)
	    << setw(5) << k + 1
	    << setw(10) << (double) k / (double) n
	    << setw(10) << A[X_][k] << setw(10) << A[Y_][k] << "\n";

    outf.close();
}


void
get_b1ob2_dnu(const int n, const ss_vect ps1[],
	      const ss_vect ps2[], double b1ob2[], double dnu[])
{
    // Estimate beta_1/beta_2 and Dnu by tracking data from two adjacent
    // BPMs.

    int j, k;
    double x1_sqr, x2_sqr, x1x2;

    cout << "\n";
    for (j = 0; j < 2; j++) {
	x1_sqr = 0.0;
	x2_sqr = 0.0;
	x1x2 = 0.0;
	for (k = 0; k < n; k++) {
	    x1_sqr += sqr(ps1[k][j]);
	    x2_sqr += sqr(ps2[k][j]);
	    x1x2 += ps1[k][j] * ps2[k][j];
	}

	x1_sqr /= n;
	x2_sqr /= n;
	x1x2 /= n;

	b1ob2[j] = x1_sqr / x2_sqr;
	dnu[j] = acos(x1x2 / sqrt(x1_sqr * x2_sqr)) / (2.0 * M_PI);

	cout << scientific << setprecision(3)
	    << "b1ob2 = " << b1ob2[j] << ", dnu = " << dnu[j] << "\n";
    }
}


void ss_est(const int n, const int bpm1, const int bpm2,
	    const bpm_data_type & bpm_data, const lin_opt_type & lin_opt)
{
    long int loc1, loc2;
    int j, k;
    double b1ob2[2], dnu[2];
    ss_vect ps1[n], ps2[n];

    for (j = 0; j < n; j++)
	for (k = 0; k < 2; k++) {
	    ps1[j][k] = bpm_data.data[k][bpm1 - 1][j];
	    ps2[j][k] = bpm_data.data[k][bpm2 - 1][j];
	}

    get_b1ob2_dnu(n, ps1, ps2, b1ob2, dnu);

    loc1 = bpm_data.loc[bpm1 - 1];
    loc2 = bpm_data.loc[bpm2 - 1];

    cout << "\n";
    cout << scientific << setprecision(3)
	<< "b1ob2 = " << lin_opt.beta[X_][loc1] / lin_opt.beta[X_][loc2]
	<< ", dnu = " << lin_opt.nu[X_][loc2] -
	lin_opt.nu[X_][loc1] << "\n";
    cout << scientific << setprecision(3)
	<< "b1ob2 = " << lin_opt.beta[Y_][loc1] / lin_opt.beta[Y_][loc2]
	<< ", dnu = " << lin_opt.nu[Y_][loc2] -
	lin_opt.nu[Y_][loc1] << "\n";
}


void prt_name(FILE * outf, const char *name)
{
    int j, k, len;

    len = strlen(name);

    j = 0;
    do {
	fprintf(outf, "%c", name[j]);
	j++;
    } while ((j < len) && (name[j] != ' '));
    fprintf(outf, ",");
    for (k = j; k < len; k++)
	fprintf(outf, "%c", name[k]);
}


int main(int argc, char *argv[])
{
    int window, cut, n_turn, bpm1, bpm2;
    bpm_data_type bpm_data;
    lin_opt_type lin_opt;
    est_lin_opt_type est_lin_opt;
    ofstream outf;

    const string file_name = "linlat_maxlab.out";

    // sls_ri_f6cwo_20.435_8.737_gset7
    lin_opt.rd_data(file_name);

    // T-b-T data.
    window = 2;
    cut = 0 * 5;
    n_turn = 2 * 1024;
    bpm1 = 6;

    bpm_data.rd_tbt("sls_tbt/tbt_090513_215959.log", lin_opt);

    prt_FFT(n_turn, cut, bpm_data.data[X_][bpm1 - 1],
	    bpm_data.data[Y_][bpm1 - 1], window);

    // Linear optics.
    window = 2;
    cut = 0;
    n_turn = 2048;

    est_lin_opt.zero(lin_opt.beta[X_].size());

    outf.open("tbt_optics.out");

    est_lin_opt.n_stats = 1;
    bpm_data.rd_tbt("sls_tbt/tbt_090513_215619.log", lin_opt);
    get_nus(outf, cut, n_turn, window, bpm_data, lin_opt, est_lin_opt);

    est_lin_opt.n_stats += 1;
    bpm_data.rd_tbt("sls_tbt/tbt_090513_215631.log", lin_opt);
    get_nus(outf, cut, n_turn, window, bpm_data, lin_opt, est_lin_opt);

    est_lin_opt.n_stats += 1;
    bpm_data.rd_tbt("sls_tbt/tbt_090513_215652.log", lin_opt);
    get_nus(outf, cut, n_turn, window, bpm_data, lin_opt, est_lin_opt);

    outf.close();

    est_lin_opt.get_stats(bpm_data, lin_opt);

    // Phase space.
    n_turn = 2048;
    bpm1 = 5;
    bpm2 = 6;

    bpm_data.rd_tbt("sls_tbt/tbt_090513_215619.log", lin_opt);
    ss_est(n_turn, bpm1, bpm2, bpm_data, lin_opt);
}
