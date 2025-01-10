import sys, re
import math, cmath
import numpy as np

'''

Numpy:

  Arrays   a[i, j, k]

  Slices   a[i:j:k], i, i+k, i+2k,..., j

  Negative i, j means n+i, n+j; e.g. for n = 10, a[-2:10], 8, 9.

'''

# Global constants.
X_ = 0; Y_ = 1; Z_ = 2
x_ = 0; px_ = 1; y_ = 2; py_ = 3


def sqr(x): return x**2

def printf(format, *args): sys.stdout.write(format % args)

def fprintf(outf, format, *args): outf.write(format % args)


class lin_opt_type (object):
    def __init__(self):
        self.loc   = []
        self.name  = []
        self.s     = np.zeros(0)
        self.alpha = np.zeros((2, 0))
        self.beta  = np.zeros((2, 0))
        self.nu    = np.zeros((2, 0))
        self.eta   = np.zeros((2, 0))
        self.etap  = np.zeros((2, 0))

    def rd_data(self, file_name):
        alpha = np.zeros(2); beta = np.zeros(2); nu = np.zeros(2);
        eta   = np.zeros(2); etap = np.zeros(2);
        prt = False

        inf = open(file_name, 'r')
        if prt: printf('\n')
        for line in inf:
            # Skip comments; lines starting with '#'.
            if line[0] != '#':
                [n, name, s, type,
                 alpha[X_], beta[X_], nu[X_], eta[X_], etap[X_],
                 alpha[Y_], beta[Y_], nu[Y_], eta[Y_], etap[Y_]] \
                     = line.strip('\n').split(',')
                self.loc.append(int(n))
                self.name  = np.append(self.name, name.lstrip())
                self.s     = np.append(self.s, float(s))
                self.alpha = np.append(self.alpha, [[alpha[X_]], [alpha[Y_]]],
                                       axis=1)
                self.beta  = np.append(self.beta, [[beta[X_]], [beta[Y_]]],
                                       axis=1)
                self.nu    = np.append(self.nu, [[nu[X_]], [nu[Y_]]], axis=1)
                self.eta   = np.append(self.eta, [[eta[X_]], [eta[Y_]]], axis=1)
                self.etap  = np.append(self.etap, [[etap[X_]], [etap[Y_]]],
                                       axis=1)
                if prt:
                    printf('%4d, %-15s, %9.5f,'
                           ' %9.5f, %8.5f, %8.5f, %8.5f, %8.5f,'
                           ' %9.5f, %8.5f, %8.5f, %8.5f, %8.5f\n',
                           n, name, s,
                           alpha[X_], beta[X_], nu[X_], eta[X_], etap[X_],
                           alpha[Y_], beta[Y_], nu[Y_], eta[Y_], etap[Y_])
        inf.close()


class bpm_data_type (object):
    def __init__(self):
        self.n_bpm  = 0
        self.n_turn = 0
        self.name   = []
        self.loc    = []
        self.data   = np.zeros((2, 0, 0))

    def get_loc(self, name, lin_opt):
        k = 0
        while k < len(lin_opt.loc) and name != lin_opt.name[k]:
            k += 1
        return k

    def get_bpm_name(self, name):
        name = re.sub ('-', '_', name).lower()
        return name

    def rd_bpm_names(self, inf, lin_opt):
        prt = False
        n_print = 8
        # Skip 1st line.
        line = inf.readline()
        # Read no of BPMs and no of turns.
        [self.n_bpm, self.n_turn] = inf.readline().strip('\n').split()
        self.n_bpm  = int(self.n_bpm)
        self.n_turn = int(self.n_turn)
        # Skip 3rd line.
        line = inf.readline()
        printf('\nno of BPMs = %d, no of turns = %d \n',
               self.n_bpm, self.n_turn)

        self.data.resize((2, self.n_bpm, self.n_turn))

        if prt: printf('\n')
        for j in range(self.n_bpm):
            name = (inf.readline().strip('\n').split())[0]
            name = self.get_bpm_name(name)
            self.name.append(name)
            self.loc.append(self.get_loc(name, lin_opt))
            if prt:
                printf(' %s', self.name[j])
                if (j+1) % n_print == 0: printf('\n')
        if prt and (self.n_bpm % n_print != 0): printf('\n')

    def rd_bpm_data(self, plane, inf):
        prt = False
        n_print = 8
        if prt: printf('\n')
        # Skip line.
        line = inf.readline()
        for k in range(self.n_turn):
            if prt: printf('\n')
            data = inf.readline().strip('\n').split()
	    # Read last BPM; on separate line.
            data.append(inf.readline().strip('\n').split()[0])
            self.data[plane, :, k] = 1e-3*np.array(data, dtype=float)
            if prt:
                for j in range(self.n_bpm):
                    printf('%10.6f', 1e3*self.data[plane, j, k])
                    if (j+1) % n_print == 0: printf('\n')
                if n_bpm % n_print != 0: printf('\n')

    def rd_tbt(self, file_name, lin_opt):
        inf = open(file_name, 'r')
        self.rd_bpm_names(inf, lin_opt)
        self.rd_bpm_data(0, inf)
        self.rd_bpm_data(1, inf)
        inf.close()


class est_lin_opt_type (object):
    def __init__(self):
        self.beta_pinger = np.zeros(2)
        self.n_stats     = 0e0
        self.alpha_mean  = np.zeros(2);      self.alpha_sigma = np.zeros(2)
        self.tune_mean   = np.zeros(2);      self.tune_sigma  = np.zeros(2)
        self.beta        = np.zeros((2, 0))
        self.beta_sum    = np.zeros((2, 0)); self.beta_sum2   = np.zeros((2, 0))
        self.beta_mean   = np.zeros((2, 0)); self.beta_sigma  = np.zeros((2, 0))
        self.nu          = np.zeros((2, 0))
        self.dnu_sum     = np.zeros((2, 0)); self.dnu_sum2    = np.zeros((2, 0))
        self.dnu_mean    = np.zeros((2, 0)); self.dnu_sigma   = np.zeros((2, 0))
        self.twoJ        = np.zeros((2, 0))
        self.phi         = np.zeros((2, 0))
        self.phi0        = np.zeros((2, 0))

    def zero(self, n):
        self.beta.resize((2, n))
        self.beta_sum.resize((2, n));  self.beta_sum2.resize((2, n))
        self.beta_mean.resize((2, n)); self.beta_sigma.resize((2, n))
        self.nu.resize((2, n))
        self.dnu_sum.resize((2, n));   self.dnu_sum2.resize((2, n))
        self.dnu_mean.resize((2, n));  self.dnu_sigma.resize((2, n))

        self.beta_sum[:, :] = 0e0;     self.beta_sum2[:, :] = 0e0
        self.dnu_sum[:, :] = 0e0;      self.dnu_sum2[:, :] = 0e0

    def get_stats(self, bpm_data, lin_opt):
        dbeta_max = 5.0; dnu_max = 0.05
        prt = False

        if prt:
            printf('\n bpm                A'
                   '                                   dnu\n')
        self.beta_mean  = np.zeros((2, bpm_data.n_bpm))
        self.beta_sigma = np.zeros((2, bpm_data.n_bpm))
        self.dnu_mean   = np.zeros((2, bpm_data.n_bpm))
        self.dnu_sigma  = np.zeros((2, bpm_data.n_bpm))
        for j in range(bpm_data.n_bpm):
            for k in range(2):
                [self.beta_mean[k, j], self.beta_sigma[k, j]] = \
                    get_m_s(self.n_stats,
                            self.beta_sum[k, j], self.beta_sum2[k, j])
                [self.dnu_mean[k, j], self.dnu_sigma[k, j]] = \
                    get_m_s(self.n_stats,
                            self.dnu_sum[k, j], self.dnu_sum2[k, j])

            if prt:
                printf('%3d  [%6.3f+/-%6.3f, %6.3f+/-%6.3f]'
                       '  [%6.3f+/-%5.3f, %6.3f+/-%5.3f]\n',
                       j+1, self.beta_mean[X_, j], self.beta_sigma[X_, j],
                       self.beta_mean[Y_, j], self.beta_sigma[Y_, j],
                       self.dnu_mean[X_, j], self.dnu_sigma[X_, j],
                       self.dnu_mean[Y_, j], self.dnu_sigma[Y_, j])

        outf = open('tbt.out', 'w')
        dbeta = np.zeros(2); dnu = np.zeros(2)
        fprintf(outf, '\n# bpm  s [m]                 beta [m]'
                '                           nu\n')
        for j in range(bpm_data.n_bpm):
            loc = bpm_data.loc[j]
            for k in range(2):
                # dbeta[k] = beta_mean[k, j] - lin_opt.beta[k, loc]
                dbeta[k] = self.beta_mean[k, j]
                if self.beta_sigma[k, j] > dbeta_max:
                    dbeta[k] = 0e0; self.beta_sigma[k, j] = 0e0

                dnu[k] = self.dnu_mean[k, j] \
                         - (lin_opt.nu[k, loc]-int(lin_opt.nu[k, loc]))
                if self.dnu_sigma[k, j] > dnu_max:
                    printf('\nBPM # %d excluded, plane = %d\n', j, k)
                    dnu[k] = 0e0; self.dnu_sigma[k, j] = 0e0

            fprintf(outf, '%4d %8.3f %7.3f +/- %5.3f %7.3f +/- %5.3f'
                    '%7.3f +/- %5.3f %7.3f +/- %5.3f %8.3f %8.3f\n',
                    j+1, lin_opt.s[loc],
                    dbeta[X_], self.beta_sigma[X_, j],
                    dbeta[Y_], self.beta_sigma[Y_, j],
                    self.dnu_mean[X_, j], self.dnu_sigma[X_, j],
                    self.dnu_mean[Y_, j], self.dnu_sigma[Y_, j],
                    lin_opt.beta[X_, loc], lin_opt.beta[Y_, loc])
        outf.close()


def DFT(x, sgn):
    n = len(x); I = complex(0e0, 1e0); X = np.zeros(n/2+1, dtype=complex)
    for j in range(n/2+1):
	for k in range(0, n):
	    X[j] += x[k]*cmath.exp(sgn*I*2e0*np.pi*float(k*j)/float(n))
    return X


def FFT1(x, window):
    n = len(x); x1 = np.zeros(n)
    for i in range(n):
	if window == 1:
	    # Rectangular.
	    x1[i] = x[i]
	elif window == 2:
	    # Sine.
	    x1[i] = math.sin(float(i)/float(n-1)*np.pi)*x[i]
	elif window == 3:
	    # Sine^2.
	    x1[i] = sqr(math.sin(float(i)/float(n-1)*np.pi))*x[i]
	else:
	    printf('FFT1: not implemented\n')
	    exit(1)
    x1 = np.fft.rfft(x1)
    # Scale FFT by 2/n for compability with 'four1' from Numerical Recipes.
    return [np.abs(x1)*2e0/n, np.angle(x1)]


def FFT2(x, window):
    n = len(x); x1 = np.zeros(n)
    for i in range(n):
	if window == 1:
	    # Rectangular.
	    x1[i] = x[i]
	elif window == 2:
	    # Sine.
	    x1[i] = math.sin(float(i)/float(n-1)*np.pi)*x[i]
	elif window == 3:
	    # Sine^2.
	    x1[i] = sqr(math.sin(float(i)/float(n-1)*np.pi))*x[i]
	else:
	    cout << 'FFT2: not implemented' << '\n'
	    exit(1)
    return np.fft.rfft(x1)


def get_ind(n, k):
    # Spectrum for real signal is mirror symmetric at k = (0, n/2).
    if k == 0:
	ind1 = 1; ind3 = 1
    elif k == n/2:
	ind1 = n/2 - 1; ind3 = n/2 - 1
    else:
	ind1 = k - 1; ind3 = k + 1
    return [ind1, ind3]


def get_nu1(n, A, k, window):
    # E. Asseo, J. Bengtsson, M. Chanel "LEAR Beam Stability Improvements
    # Using FFT Analysis" EPAC 1988.
    nu = 0e0
    [ind1, ind3] = get_ind(n, k)
    if A[ind3] > A[ind1]:
	A1 = A[k]; A2 = A[ind3]; ind = k
    else:
	A1 = A[ind1]; A2 = A[k]
	# Special case for 0 frequency.
	if k != 0:
            ind = ind1
        else:
            ind = -1
    # Avoid division by zero.
    if A1+A2 != 0e0:
	if (window == 1):
	    nu = (ind+A2/(A1+A2))/n
	elif window == 2:
	    nu = (ind-0.5e0+2e0*A2/(A1+A2))/n
	elif window == 3:
	    nu = (ind-1e0+3e0*A2/(A1+A2))/n
	else:
	    cout << 'get_nu1: not defined\n'
    else:
	nu = 0e0
    return nu


def sinc(omega):
    if omega != 0e0:
        return math.sin(omega)/omega
    else:
        return 1e0


def get_A(n, A, nu, k, window):
    corr = 0e0
    if window == 1:
	corr = sinc(np.pi*(k-nu*n))
    elif window == 2:
	corr = (sinc(np.pi*(k+0.5e0-nu*n))+sinc(np.pi*(k-0.5e0-nu*n)))/2e0
    elif window == 3:
	cout << 'get_A: not implemented\n'
	exit(1)
    else:
	cout << 'get_A: not defined\n'
    return A[k]/corr


def get_alpha(n, X, nu, k):
    # M. Bertocco, C. Offelli, D. Petri "Analysis of Damped Sinusoidal Signals
    # via a Frequency-Domain Interpolation Algorithm" IEEE 43, 245-250 (1994).
    I = complex(0e0, 1e0)
    [ind1, ind3] = get_ind(n, k)
    if abs(X[ind3]) > abs(X[ind1]):
	d = 1; rho = X[ind3]/X[k]
    else:
	d = -1; rho = X[ind1]/X[k]
    z = (1e0-rho)/(1e0-rho*cmath.exp(-I*2e0*np.pi*float(d)/float(n)))
    delta = n*cmath.phase(z)/(2e0*np.pi)
    alpha = n*math.log(abs(z))/(2e0*np.pi)
    return [delta, alpha]


def get_peak(n, A):
    k = 0
    peak = 0e0
    for ind2 in range(n/2+1):
	[ind1, ind3] = get_ind(n, ind2)
	if (A[ind2] > peak) and (A[ind1] < A[ind2]) and (A[ind2] > A[ind3]):
	    peak = A[ind2]
	    k = ind2
    return k


def get_phi(n,  k,  nu, phi):
    phi_nu = phi[k] - (n*nu-k)*np.pi
    if phi_nu > np.pi:
	phi_nu -= 2e0*np.pi
    elif phi_nu < -np.pi:
	phi_nu += 2e0*np.pi
    return phi_nu


def get_nu2(n, x, window):
    [A, phi] = FFT1(x, window)
    k = get_peak(n, A)
    nu = get_nu1(n, A, k, window)
    A_nu = get_A(n, A, nu, k, window)
    phi_nu = get_phi(n, k, nu, phi)
    # Rectangular window.
    x = FFT2(x, 1)
    [delta, alpha] = get_alpha(n, x, nu, k)
    return [nu, A_nu, phi_nu, delta, alpha]


def rm_mean(n, x): mean = sum(x); mean /= n; x -= mean


def get_nus(outf, cut, n,  window, bpm_data,  lin_opt,  est_lin_opt):
    prt = False
    sgn = [1, -1]

    tune_sum  = np.zeros(2); tune_sum2  = np.zeros(2)
    alpha_sum = np.zeros(2); alpha_sum2 = np.zeros(2)
    twoJ_sum  = np.zeros(2); twoJ_sum2  = np.zeros(2)
    phi0_sum  = np.zeros(2); phi0_sum2  = np.zeros(2)

    x = np.zeros(n)
    tunes = np.zeros((bpm_data.n_bpm, 2)); As = np.zeros((bpm_data.n_bpm, 2))
    phis = np.zeros((bpm_data.n_bpm, 2))
    delta = np.zeros(2); alpha = np.zeros(2)
    nus = np.zeros((bpm_data.n_bpm, 2)); phi0 = np.zeros(2)
    printf('\n');
    for i in range(bpm_data.n_bpm):
	loc = bpm_data.loc[i]
	for j in range(2):
	    x = bpm_data.data[j, i, cut:n+cut]; rm_mean(n, x)
 
	    [tunes[i, j], As[i, j], phis[i, j], delta[j], alpha[j]] = \
                get_nu2(n, x, window)

	    if sgn[j] < 0: phis[i, j] = -phis[i, j]
	    if phis[i, j] < 0e0: phis[i, j] += 2e0*np.pi
	    nus[i, j] = phis[i, j]/(2e0*np.pi)

	    tune_sum[j] += tunes[i, j]; tune_sum2[j] += sqr(tunes[i, j])
	    alpha_sum[j] += alpha[j]; alpha_sum2[j] += sqr(alpha[j])

	    twoJ = sqr(As[i, j])/lin_opt.beta[j, loc]
	    twoJ_sum[j] += twoJ; twoJ_sum2[j] += sqr(twoJ)

	    phi0[j] = (nus[i, j]-(lin_opt.nu[j, loc]
                                  -int(lin_opt.nu[j, loc])))*2e0*np.pi
	    if phi0[j] < 0e0: phi0[j] += 2e0*np.pi
	    phi0_sum[j] += phi0[j]; phi0_sum2[j] += sqr(phi0[j])

	# if (prt) printf('[%8.6f, %8.6f]\n', tunes[i, X_],
	# tunes[i, Y_])

    twoJ_mean = twoJ_sum/bpm_data.n_bpm
    twoJ_sigma = np.sqrt((bpm_data.n_bpm*twoJ_sum2-sqr(twoJ_sum))
		         /(bpm_data.n_bpm*(bpm_data.n_bpm-1e0)))

    phi0_mean = phi0_sum/bpm_data.n_bpm
    phi0_sigma = np.sqrt((bpm_data.n_bpm*phi0_sum2-np.square(phi0_sum))
		         /(bpm_data.n_bpm*(bpm_data.n_bpm-1e0)))

    printf('\ntwoJ  = [%9.3e+/-%9.3e, %9.3e+/-%9.3e]'
           ', phi0 = [%5.3f+/-%5.3f, %5.3f+/-%5.3f]\n',
           twoJ_mean[X_], twoJ_sigma[X_], twoJ_mean[Y_], twoJ_sigma[Y_],
           phi0_mean[X_], phi0_sigma[X_], phi0_mean[Y_], phi0_sigma[Y_])
    printf('A0    = [%5.3f, %5.3f] mm\n',
           1e3*math.sqrt(twoJ_mean[X_]*est_lin_opt.beta_pinger[X_]),
           1e3*math.sqrt(twoJ_mean[Y_]*est_lin_opt.beta_pinger[Y_]))

    # Normalize.
    if prt:
        printf('\n bpm        A               nu            nu (model)\n')
    dnu = np.zeros(2)
    for i in range(bpm_data.n_bpm):
	loc = bpm_data.loc[i]
	for j in range(2):
	    beta = sqr(As[i, j])/twoJ_mean[j]

	    nus[i, j] -= phi0_mean[j]/(2e0*np.pi)
	    if nus[i, j] < 0e0:	nus[i, j] += 1e0

	    dnu[j] = nus[i, j]-(lin_opt.nu[j, loc] -int(lin_opt.nu[j, loc]))
	    if dnu[j] < -0.5e0: dnu[j] += 1e0
	    if dnu[j] > 0.5e0:	dnu[j] -= 1e0

	    est_lin_opt.beta_sum[j, i]  += beta
	    est_lin_opt.beta_sum2[j, i] += sqr(beta)
	    est_lin_opt.dnu_sum[j, i]   += dnu[j]
	    est_lin_opt.dnu_sum2[j, i]  += sqr(dnu[j])

	fprintf(outf, '%4d %7.3f %8.3f %8.3f\n',
                i+1, lin_opt.s[loc], dnu[X_], dnu[Y_])

	if prt:
	    printf('%3d  [%6.3e, %5.3e]  [%6.3e, %5.3e]  [%6.3f, %5.3f]\n',
                   i+1, 1e3*As[i, X_], 1e3*As[i, Y_],
                   nus[i, X_], nus[i, Y_],
                   lin_opt.nu[X_, loc]-int(lin_opt.nu[X_, loc]),
                   lin_opt.nu[Y_, loc]-int(lin_opt.nu[Y_, loc]))

    for j in range(2):
	est_lin_opt.tune_mean[j] = tune_sum[j]/bpm_data.n_bpm
	if sgn[j] < 0:
	    est_lin_opt.tune_mean[j] = 1e0 - est_lin_opt.tune_mean[j]
	est_lin_opt.tune_sigma[j] = \
            math.sqrt((bpm_data.n_bpm*tune_sum2[j]-sqr(tune_sum[j]))
		   /(bpm_data.n_bpm*(bpm_data.n_bpm-1e0)))

	est_lin_opt.alpha_mean[j] = alpha_sum[j]/bpm_data.n_bpm
	est_lin_opt.alpha_sigma[j] = \
	    math.sqrt((bpm_data.n_bpm*alpha_sum2[j]-sqr(alpha_sum[j]))
		   /(bpm_data.n_bpm*(bpm_data.n_bpm-1e0)))

    printf('\nnu    = [%9.6f+/-%8.6f, %9.6f+/-%8.6f]\n',
           est_lin_opt.tune_mean[X_], est_lin_opt.tune_sigma[X_],
           est_lin_opt.tune_mean[Y_], est_lin_opt.tune_sigma[Y_])
    printf('alpha = [%9.6f+/-%8.6f, %9.6f+/-%8.6f]\n',
           est_lin_opt.alpha_mean[X_], est_lin_opt.alpha_sigma[X_],
           est_lin_opt.alpha_mean[Y_], est_lin_opt.alpha_sigma[Y_])

    printf('%8.5f %8.5f\n', nus[6, X_]-nus[5, X_], nus[6, Y_]-nus[5, Y_])


def get_m_s(n, sum, sum2):
    mean = sum/n; sigma = math.sqrt((n*sum2-sqr(sum))/(n*(n-1e0)))
    return [mean, sigma]


def prt_FFT(cut, xy, window):
    n = len(xy[X_])
    outf = open('sls.out', 'w')
    x1 = xy[:, cut:n+cut]
    for j in range(cut, n+cut):
	fprintf(outf, '%5d %11.3e %11.3e\n', j+1, xy[X_, j], xy[Y_, j])
    outf.close()

    A = np.zeros((2, n/2+1)); phi = np.zeros((2, n/2+1))
    for k in range(0, 2):
	[A[k], phi[k]] = FFT1(x1[k], window)

    outf = open('sls_fft.out', 'w')
    for k in range(n/2+1):
        fprintf(outf, '%5d %9.3e %9.3e %9.3e\n',
                k+1, float(k)/float(n), A[X_, k], A[Y_, k])
    outf.close()


def get_b1ob2_dnu(n, ps1, ps2):
    # Estimate beta_1/beta_2 and dnu by tracking data from two adjacent BPMs.
    printf('\n')
    b1ob2 = np.zeros(2); dnu = np.zeros(2)
    for k in range(2):
	x1_sqr = np.sum(sqr(ps1[k])); x2_sqr = np.sum(sqr(ps2[k]))
	x1x2 = np.sum(ps1[k]*ps2[k])
	x1_sqr /= n; x2_sqr /= n; x1x2 /= n
	b1ob2[k] = x1_sqr/x2_sqr
	dnu[k] = math.acos(x1x2/math.sqrt(x1_sqr*x2_sqr))/(2e0*np.pi)

	printf('b1ob2 = %9.3e, dnu = %5.3f\n', b1ob2[k], dnu[k])
    return [b1ob2, dnu]


def ss_est(cut, n, bpm1, bpm2, bpm_data, lin_opt):
    ps1 = np.zeros((2, n)); ps2 = np.zeros((2, n))
    ps1[:, :] = bpm_data.data[:, bpm1-1, cut:n+cut]
    ps2[:, :] = bpm_data.data[:, bpm2-1, cut:n+cut]

    # Get estimated linear optics parameters.
    [b1ob2, dnu] = get_b1ob2_dnu(n, ps1, ps2)

    # Get reference linear optics parameters.
    loc1 = bpm_data.loc[bpm1-1]; loc2 = bpm_data.loc[bpm2-1]
    beta1 = lin_opt.beta[:, loc1]
    b1ob2 = lin_opt.beta[:, loc1]/lin_opt.beta[:, loc2]
    dnu   = lin_opt.nu[:, loc2]-lin_opt.nu[:, loc1]
    printf('\n')
    for k in range(2):
        printf('b1ob2 = %9.3e, dnu = %5.3f\n', b1ob2[k], dnu[k])

    outf = open('tbt_phase_space.out', 'w')
    ps = np.zeros(4); twoJ = np.zeros(2)
    for j in range(n):
        for k in range(2):
            ps[2*k] = ps1[k, j]/math.sqrt(beta1[k])
            ps[2*k+1] = (math.sqrt(b1ob2[k])*ps2[k, j] \
                         -ps1[k, j]*math.cos(2e0*math.pi*dnu[k])) \
                        /(math.sqrt(beta1[k])*math.sin(2e0*math.pi*dnu[k]))
            twoJ[k] = sqr(ps[2*k]) + sqr(ps[2*k+1])
        fprintf(outf, '%4d', j+1)
        fprintf(outf, '%13.5e %13.5e %13.5e %13.5e %13.5e %13.5e\n',
                ps[x_], ps[px_], ps[y_], ps[py_], twoJ[X_], twoJ[Y_])
    outf.close()
    

def prt_name(outf, name):
    strlen = len(name); j = 0
    while (j < strlen) and (name[j] != ' '):
	fprintf(outf, '%c' % (name[j]))
	j += 1
    fprintf(outf, ',')
    for k in range(j, strlen):
	fprintf(outf, '%c' % (name[k]))


def main():
    bpm_data    = bpm_data_type()
    lin_opt     = lin_opt_type()
    est_lin_opt = est_lin_opt_type()

    home_dir = sys.argv[1] + '/'

    # sls_ri_f6cwo_20.435_8.737_gset7
    file_name = '/linlat_maxlab.out'

    lin_opt.rd_data(home_dir+file_name)

    # Turn-by-turn data.
    window = 2; cut = 0*5; n_turn = 2*1024; bpm1 = 6

    bpm_data.rd_tbt(home_dir+'tbt_090513_215959.log', lin_opt)
    prt_FFT(cut, bpm_data.data[:, bpm1-1], window)

    # Linear optics.
    window = 2; cut = 0; n_turn = 2048

    est_lin_opt.zero(len(bpm_data.loc))
    est_lin_opt.beta_pinger = (6.92e0, 6.76e0)

    outf = open('tbt_optics.out', 'w')

    est_lin_opt.n_stats = 1
    bpm_data.rd_tbt(home_dir+'tbt_090513_215619.log', lin_opt)
    get_nus(outf, cut, n_turn, window, bpm_data, lin_opt, est_lin_opt)

    est_lin_opt.n_stats += 1
    bpm_data.rd_tbt(home_dir+'tbt_090513_215631.log', lin_opt)
    get_nus(outf, cut, n_turn, window, bpm_data, lin_opt, est_lin_opt)

    est_lin_opt.n_stats += 1
    bpm_data.rd_tbt(home_dir+'tbt_090513_215652.log', lin_opt)
    get_nus(outf, cut, n_turn, window, bpm_data, lin_opt, est_lin_opt)

    outf.close()

    est_lin_opt.get_stats(bpm_data, lin_opt)

    # Phase space.
    cut = 10; n_turn = 1024; bpm1 = 5; bpm2 = 6

    bpm_data.rd_tbt(home_dir+'tbt_090513_220010.log', lin_opt)
    ss_est(cut, n_turn, bpm1, bpm2, bpm_data, lin_opt)


main()
