#define NO 1

#include "tracy_lib.h"

int no_tps = NO;

extern bool  freq_map;
extern int   N_Fam, Q_Fam[];

// dynamic aperture
const double  delta  = 3e-2; // delta for off-momentum aperture

double  R0 = 25e-3;


const bool  ID = false, get_DA = true, tune_shift = false, DW100 = true;

const int wp = 0;
// const int wp = 7;
// const int wp = 8
// const int wp = 9;
// const int wp = 10;
// Baseline Lattice
// const int wp = 12;
// Commissioning
// const int wp = 13;

const int get_thor = 0;

const int  max_Fnum = 10, max_quads = 200, max_loc = 10;

bool      IDs;
int       Fnum[max_Fnum], loc[max_loc], n_iter, n_prm, n_loc;
int       quads1[max_quads], quads2[max_quads], sext_scheme;
double    alpha_mp[2][2], beta_mp[2][2], dnu_mp[2], b2s[max_quads];
double    Jx, Jy;
Vector2   nu;
FILE      *fp_lat;
std::ofstream  tune_out;


const char file_lat[] = "get_bn.lat";
const char home_dir[] = "/home/bengtsson";

const int  n_peaks = 10;


void set_roll_design(const int Fnum, const int Knum, const double roll)
{
  // Set design roll [deg].

  Cell[Elem_GetPos(Fnum, Knum)].Elem.M->PdTpar = roll;

  Mpole_SetdT(Fnum, Knum);
}


void set_roll_design_fam(const int Fnum, const double roll)
{
  // Set design roll [deg].
  int k;

  for (k = 1; k <= GetnKid(Fnum); k++)
    set_roll_design(Fnum, k, roll);
}


double get_PS_vol(const ss_vect<tps> A, const ss_vect<double> ps,
		  double twoJ[])
{
  /* Note:

            T                  -1
       M S M  = S,    M = A R A  ,    A, R symplectic

        -1       T           |  0  I |        T   | beta   -alpha |
       A   = -S A  S,    S = |       |,    A A  = |               |
                             | -I  0 |            | -alpha  gamma |

          T             T        T  T
       A A  -> M A (M A)  = M A A  M

     Transform to Floquet Space:

        -1         T
       A   x = -S A  S x,

                               -1    T  -1            T    T
       2J_1 + 2J_2 + 2J_3 = ( A   x )  A   x = ( S x )  A A  ( S x )

  */

  int              i, j;
  double           V;
  ss_vect<double>  dps;
  ss_vect<tps>     Id, J, u, v;

  Id.identity(); J.zero();
  J[x_] = Id[px_]; J[px_] = -Id[x_]; J[y_] = Id[py_]; J[py_] = -Id[y_];
  J[ct_] = Id[delta_]; J[delta_] = -Id[ct_];

  // relative to fixed point
  dps = ps - globval.CODvect;

  // relative to delta dependent fixed point
  for (i = 0; i <= 1; i++) {
    dps[2*i] -= Cell[0].Eta[i]*ps[delta_];
    dps[2*i+1] -= Cell[0].Etap[i]*ps[delta_];
  }

  u = J*dps;

  for (i = 0; i <= 5; i++) {
    v[i] = 0.0;
    for (j = 0; j <= 5; j++)
      v[i] += A[j][i]*u[j];
  }

  V = 0.0;
  for (i = 0; i <= 2; i++) {
    twoJ[i] = sqr(v[2*i].cst()) + sqr(v[2*i+1].cst());
    if (i < 2 || globval.Cavity_on) V += twoJ[i];
  }

  return V;
}


void get_lin_inv(const double Ax, const double Ay, const double delta,
		const int n)
{
  int              i;
  long int         lastpos;
  double           twoJ[3], V;
  ss_vect<double>  ps;
  ss_vect<tps>     A;
  std::fstream          outf;

  globval.Cavity_on = true;

  Ring_GetTwiss(true, 0.0); printglob();

  putlinmat(6, globval.Ascr, A);

  outf.open("lin_inv.out", std::ios::out);
  ps.zero(); ps[x_] = Ax; ps[y_] = Ay; ps[delta_] = delta;
  for (i = 1; i <= n; i++) {
    Cell_Pass(0, globval.Cell_nLoc, ps, lastpos); V = get_PS_vol(A, ps, twoJ);

    outf << std::setprecision(3) << std::scientific
	 << std::setw(3) << i << std::setw(10) << 1e6*twoJ[X_]
	 << std::setw(10) << 1e6*twoJ[Y_] << std::setw(10) << 1e3*1e2*twoJ[Z_]
	 << std::setw(10) << 1e6*V << std::endl;
  }
  outf.close();
}


void get_sympl_form(const double Ax, const double Ay, const double delta,
		    const int n)
{
  int              i;
  long int         lastpos;
  double           omega[3], V;
  ss_vect<double>  ps1, ps2, ps3, dq, dp;
  std::fstream          outf;

  const double  eps = 1e-10;

  globval.Cavity_on = true;

  outf.open("sympl_form.out", std::ios::out);
  ps1.zero(); ps1[x_] = Ax; ps1[y_] = Ay; ps1[delta_] = delta;
  ps2 = ps1; ps3 = ps1;
  for (i = 0; i <= 5; i++)
    if (i % 2 == 0)
      ps2[i] += eps;
    else
      ps3[i] += eps;
  for (i = 1; i <= n; i++) {
    Cell_Pass(0, globval.Cell_nLoc, ps1, lastpos);
    Cell_Pass(0, globval.Cell_nLoc, ps2, lastpos);
    Cell_Pass(0, globval.Cell_nLoc, ps3, lastpos);
    dq = ps2 - ps1; dp = ps3 - ps1;
    omega[X_] = dq[x_]*dp[px_] - dp[x_]*dq[px_];
    omega[Y_] = dq[y_]*dp[py_] - dp[y_]*dq[py_];
    omega[Z_] = -dq[ct_]*dp[delta_] + dp[ct_]*dq[delta_];
    V = omega[X_] + omega[Y_] + omega[Z_];

    outf << std::setprecision(3) << std::scientific
	 << std::setw(3) << i << std::setw(11) << omega[X_]
	 << std::setw(11) << omega[Y_] << std::setw(11) << omega[Z_] << std::setw(11) << V
	 << std::endl;
  }
  outf.close();
}


void get_fixed_points(const double Ax, const double Ay, const int n)
{
  bool             lost;
  int              i, j, k, l, n1;
  long int         lastpos;
  double           A[2], twoJ2[2], twoJ_sum[2], phi0[2], phi1[2];
  double           dphi, dphi2[2], dphi_sum[2], twoJ[2];
  double           r_max, s_twoJ[2], s_phi[2];
  ss_vect<double>  ps, ps_Fl;
  std::fstream          outf;

  const int     n_turn  = 5;
  const double  A_min = 0.1e-3, twoJ_max = 1e-5, phi_max = 2e-1*0.5*2.0*M_PI;

  outf.open("fixed_points.out", std::ios::out);

  r_max = sqrt(sqr(Ax)+sqr(Ay));
  n1 = 0;
  for (i = -n; i <= n; i++) {
    A[X_] = i*Ax/n;
    if (A[X_] == 0.0) A[X_] = A_min;
    for (j = -n; j <= n; j++) {
      n1++;
      A[Y_] = j*Ay/n;
      if (A[Y_] == 0.0) A[Y_] = A_min;

      ps.zero(); ps[x_] = A[X_]; ps[y_] = A[Y_];
      ps_Fl = ps; LinTrans(4, globval.Ascrinv, ps_Fl);

      for (l = 0; l <= 1; l++) {
	phi0[l] = -atan2(ps_Fl[2*l+1], ps_Fl[2*l]);

	twoJ2[l] = 0.0; twoJ_sum[l] = 0.0;
	dphi2[l] = 0.0; dphi_sum[l] = 0.0;
      }

      for (k = 1; k <= n_turn; k++) {
	Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
	lost = lastpos != globval.Cell_nLoc;
	if (!lost) {
	  ps_Fl = ps; LinTrans(4, globval.Ascrinv, ps_Fl);
	  for (l = 0; l <= 1; l++) {
	    twoJ[l] = sqr(ps_Fl[2*l]) + sqr(ps_Fl[2*l+1]);
	    phi1[l] = -atan2(ps_Fl[2*l+1], ps_Fl[2*l]);

	    dphi = phi1[l] - phi0[l];
	    if (dphi < 0.0) dphi += 2.0*M_PI;

	    twoJ2[l] += sqr(twoJ[l]); twoJ_sum[l] += twoJ[l];
	    dphi2[l] += sqr(dphi); dphi_sum[l] += dphi;

	    phi0[l] = phi1[l];
	  }
	} else
	  break;
      }

      for (l = 0; l <= 1; l++)
	if (!lost) {
	  s_twoJ[l] = sqrt((n_turn*twoJ2[l]-sqr(twoJ_sum[l]))
                      /(n_turn*(n_turn-1.0)));
	  s_twoJ[l] = min(s_twoJ[l], twoJ_max);

	  s_phi[l] = sqrt((n_turn*dphi2[l]-sqr(dphi_sum[l]))
                     /(n_turn*(n_turn-1.0)));
	  s_phi[l] = min(s_phi[l], phi_max);
	} else {
	  s_twoJ[l] = twoJ_max; s_phi[l] = phi_max;
	}

      outf << std::fixed << std::setprecision(3)
	   << std::setw(3) << n1
	   << std::setw(8) << 1e3*A[X_] << std::setw(8) << 1e3*A[Y_]
	   << std::scientific
	   << std::setw(10) << 1e6*s_twoJ[X_] << std::setw(10) << 1e6*s_twoJ[Y_]
	   << std::fixed << std::setprecision(5)
	   << std::setw(10) << s_phi[X_] << std::setw(10) << s_phi[Y_]
	   << std::scientific << std::setprecision(3)
	   << std::setw(10) << log(1.0+sqr(s_twoJ[X_])+sqr(s_twoJ[Y_]))
	   << std::setw(10) << min(s_phi[X_], s_phi[Y_])
	   << std::endl;
    }
    outf << std::endl;
  }
  outf.close();
}


void get_phi_stat(const double Ax, const double Ay, const int n)
{
  bool             lost;
  int              i, j, k, l, n1;
  long int         lastpos;
  double           A[2], twoJ2[2], twoJ_sum[2], phi0[2], phi1[2];
  double           dphi, dphi2[2], dphi_sum[2], twoJ[2];
  double           s_twoJ[2], s_phi[2];
  ss_vect<double>  ps, ps_Fl;
  std::fstream          outf;

  const int     n_turn  = 15;
  const double  A_min = 0.1e-3, twoJ_max = 1e-5, phi_max = M_PI;

  outf.open("phi_stat.out", std::ios::out);

  n1 = 0;
  for (i = -n; i <= n; i++) {
    A[X_] = i*Ax/n;
    if (A[X_] == 0.0) A[X_] = A_min;
    for (j = -n; j <= n; j++) {
      n1++;
      A[Y_] = j*Ay/n;
      if (A[Y_] == 0.0) A[Y_] = A_min;

      ps.zero(); ps[x_] = A[X_]; ps[y_] = A[Y_];
      ps_Fl = ps; LinTrans(4, globval.Ascrinv, ps_Fl);

      for (l = 0; l <= 1; l++) {
	phi0[l] = -atan2(ps_Fl[2*l+1], ps_Fl[2*l]);

	twoJ2[l] = 0.0; twoJ_sum[l] = 0.0;
	dphi2[l] = 0.0; dphi_sum[l] = 0.0;
      }

      for (k = 1; k <= n_turn; k++) {
	Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
	lost = lastpos != globval.Cell_nLoc;
	if (!lost) {
	  ps_Fl = ps; LinTrans(4, globval.Ascrinv, ps_Fl);
	  for (l = 0; l <= 1; l++) {
	    twoJ[l] = sqr(ps_Fl[2*l]) + sqr(ps_Fl[2*l+1]);
	    phi1[l] = -atan2(ps_Fl[2*l+1], ps_Fl[2*l]);

	    dphi = phi1[l] - phi0[l];
	    if (dphi < 0.0) dphi += 2.0*M_PI;

	    twoJ2[l] += sqr(twoJ[l]); twoJ_sum[l] += twoJ[l];
	    dphi2[l] += sqr(dphi); dphi_sum[l] += dphi;

	    phi0[l] = phi1[l];
	  }
	} else
	  break;
      }

      for (l = 0; l <= 1; l++)
	if (!lost) {
	  s_twoJ[l] = sqrt((n_turn*twoJ2[l]-sqr(twoJ_sum[l]))
                      /(n_turn*(n_turn-1.0)));
	  s_twoJ[l] = min(s_twoJ[l], twoJ_max);

	  s_phi[l] = sqrt((n_turn*dphi2[l]-sqr(dphi_sum[l]))
                     /(n_turn*(n_turn-1.0)));
	  s_phi[l] = min(s_phi[l], phi_max);
	} else {
	  s_twoJ[l] = twoJ_max; s_phi[l] = phi_max;
	}

      outf << std::fixed << std::setprecision(3)
	   << std::setw(3) << n1
	   << std::setw(8) << 1e3*A[X_] << std::setw(8) << 1e3*A[Y_]
	   << std::scientific
	   << std::setw(10) << 1e6*s_twoJ[X_] << std::setw(10) << 1e6*s_twoJ[Y_]
	   << std::fixed << std::setprecision(5)
	   << std::setw(10) << s_phi[X_] << std::setw(10) << s_phi[Y_]
	   << std::endl;
    }
    outf << std::endl;
  }
  outf.close();
}


bool eng_tol(const bool config, const int n_orbit)
{
  bool  cod;
  int   i;

  if (config) {
    globval.bpm = ElemIndex("bpm");
    globval.hcorr = ElemIndex("chv"); globval.vcorr = ElemIndex("chv");

    // Clear trim setpoints
    for (i = 1; i <= GetnKid(globval.vcorr); i++)
      SetKLpar(globval.vcorr, i, -Dip, 0.0);
    for (i = 1; i <= GetnKid(globval.hcorr); i++)
      SetKLpar(globval.hcorr, i, Dip, 0.0);
  }

  misalign_rms_type(Quad, 100e-6, 100e-6, 0.0, true);

  cod = orb_corr(n_orbit);

  set_bnr_rms_type(Quad, Quad, 0.5e-3, 0.0, true);

  return cod;
}


void ID_correct(void)
{

  N_Fam = 0;

  if (true) {
    Q_Fam[N_Fam++] = ElemIndex("qh1");
    Q_Fam[N_Fam++] = ElemIndex("qh2");
    Q_Fam[N_Fam++] = ElemIndex("qh3");

    Q_Fam[N_Fam++] = ElemIndex("ql1");
    Q_Fam[N_Fam++] = ElemIndex("ql2");
    Q_Fam[N_Fam++] = ElemIndex("ql3");
  } else {
    Q_Fam[N_Fam++] = ElemIndex("qh1sp01");
    Q_Fam[N_Fam++] = ElemIndex("qh2sp01");
    Q_Fam[N_Fam++] = ElemIndex("qh3sp01");

    Q_Fam[N_Fam++] = ElemIndex("ql1sp01");
    Q_Fam[N_Fam++] = ElemIndex("ql2sp01");
    Q_Fam[N_Fam++] = ElemIndex("ql3sp01");

    Q_Fam[N_Fam++] = ElemIndex("qh1sp02");
    Q_Fam[N_Fam++] = ElemIndex("qh2sp02");
    Q_Fam[N_Fam++] = ElemIndex("qh3sp02");

    Q_Fam[N_Fam++] = ElemIndex("ql1sp02");
    Q_Fam[N_Fam++] = ElemIndex("ql2sp02");
    Q_Fam[N_Fam++] = ElemIndex("ql3sp02");
  }

  printf("\n");
  printf("No of quadrupole families: %d\n", N_Fam);

  ini_ID_corr(true);

  if (false) orb_corr(5);

  // 3 DWs.
  ID_corr(10, 2, true);
  // 5 DWs.
  // ID_corr(10, 3, true);
}


void prt_dlat(const char *fname, const bool all)
{
  long int      i = 0;
  double        code = 0.0;
  FILE          *outf;

  outf = file_write(fname);

  for (i = 0; i <= globval.Cell_nLoc; i++) {
    if (all ||
	((Cell[i].Elem.Pkind == Mpole) &&
	 (Cell[i].Elem.M->n_design == Sext))) {
      switch (Cell[i].Elem.Pkind) {
      case drift:
	code = 0.0;
	break;
      case Mpole:
	if (Cell[i].Elem.M->Pirho != 0.0)
	  code = 0.5;
	else if (Cell[i].Elem.M->PBpar[Quad+HOMmax] != 0)
	  code = sgn(Cell[i].Elem.M->PBpar[Quad+HOMmax]);
	else if (Cell[i].Elem.M->PBpar[Sext+HOMmax] != 0)
	  code = 1.5*sgn(Cell[i].Elem.M->PBpar[Sext+HOMmax]);
	else if (Cell[i].Fnum == globval.bpm)
	  code = 2.0;
	else
	  code = 0.0;
	break;
      default:
	code = 0.0;
	break;
      }

      fprintf(outf, "%4ld %15s %6.2f %4.1f %6.2f %7.4f %6.2f %7.4f\n",
	      i, Cell[i].Elem.PName, Cell[i].S, code,
	      1e2*(Cell[i].Beta[X_]-betas0_[i][X_])/betas0_[i][X_],
	      Cell[i].Nu[X_]-nus0_[i][X_],
	      1e2*(Cell[i].Beta[Y_]-betas0_[i][Y_])/betas0_[i][Y_],
	      Cell[i].Nu[Y_]-nus0_[i][Y_]);
    }
  }

  fclose(outf);
}


void prt_ZAP(const int n)
{
  long int  k;
  FILE      *outf;

  outf = file_write("ZAPLAT.DAT");

  fprintf(outf, "%ld %7.5f\n",
	  globval.Cell_nLoc+1, n*Cell[globval.Cell_nLoc].S);
  fprintf(outf, "One super period\n");

  for (k = 0; k <= globval.Cell_nLoc; k++)
    fprintf(outf, "%10.5f %8.5f %9.6f %8.5f %7.3f %8.5f %8.5f %7.5f\n",
	    Cell[k].S,
	    Cell[k].Beta[X_], Cell[k].Alpha[X_],
	    Cell[k].Beta[Y_], Cell[k].Alpha[Y_],
	    Cell[k].Eta[X_], Cell[k].Etap[X_], Cell[k].maxampl[X_][1]);

  fprintf(outf, "0\n");

  fclose(outf);
}


void chk_IBS(void)
{
  int       k;
  double    chi;
  std::ofstream  outf;

  const double  chi_min = 1e-9, chi_max = 1e-5, fact = 2.0;


  file_wr(outf, "chk_IBS.out");

  k = 0; chi_m = chi_min;
  while (chi_m <= chi_max) {
    k++;

    chi = f_IBS(chi_m);

    outf << std::scientific << std::setprecision(5)
	 << std::setw(4) << k << std::setw(12) << chi_m << std::setw(12) << chi;

    chi = get_int_IBS();

    outf << std::scientific << std::setprecision(5)
	 << std::setw(12) << chi << std::endl;

    chi_m *= fact;
  }

  outf.close();
}


void get_Touschek(void)
{
  long int      k;
  int           j;
  double        gamma_z, eps[3];
  double        sum_delta[globval.Cell_nLoc+1][2];
  double        sum2_delta[globval.Cell_nLoc+1][2];
  FILE          *fp;

  const double  Qb = 1.3e-9;

  globval.eps[X_] = 2.040e-9;
  globval.eps[Y_] = 8e-12;
  globval.eps[Z_] = 1.516e-6;

  globval.alpha_z = 2.502e-02; globval.beta_z = 5.733;

  if (true) {
    printf("\n");
    printf("alpha: %11.3e %11.3e %11.3e\n",
	   globval.alpha_rad[X_], globval.alpha_rad[Y_],
	   globval.alpha_rad[Z_]);

    for (j = 0; j < 3; j++)
      eps[j] = globval.eps[j];

    for (j = 1; j <= 0; j++)
      IBS(Qb, globval.eps, eps, true, true);

    for (j = 0; j < 3; j++)
      eps[j] = globval.eps[j];

    for (j = 1; j <= 5; j++)
      IBS_BM(Qb, globval.eps, eps, true, true);
  }

  if (true) {
    globval.delta_RF = 3.0e-2;
//     Touschek(Qb, globval.delta_RF, 0.85e-9, 8e-12, 0.84e-3, 4.4e-3);
    gamma_z = (1.0+sqr(globval.alpha_z))/globval.beta_z;
//     Touschek(Qb, globval.delta_RF, eps[X_], eps[Y_],
// 	     sqrt(gamma_z*eps[Z_]), sqrt(globval.beta_z*eps[Z_]));
    Touschek(Qb, 3.03e-2, 2.446e-9, 9.595e-12, 0.580e-3, 3.33e-3);

    if (false) {
      // initialize momentum aperture arrays
      for(k = 0; k <= globval.Cell_nLoc; k++){
	sum_delta[k][0] = 0.0; sum_delta[k][1] = 0.0;
	sum2_delta[k][0] = 0.0; sum2_delta[k][1] = 0.0;
      }

      Touschek(1.3e-9, globval.delta_RF, false,
	       0.85e-9, 0.008e-9, 0.84e-3, 4.4e-3,
	       512, true, sum_delta, sum2_delta);

      fp = file_write("mom_aper.out");
      for(k = 0; k <= globval.Cell_nLoc; k++)
	fprintf(fp, "%4ld %7.2f %5.3f %6.3f\n",
		k, Cell[k].S, 1e2*sum_delta[k][0], 1e2*sum_delta[k][1]);
    }
  }
}

void get_MM_TP_M(double **M, ss_vect<tps> &map)
{
  int  i, j, k, l, m, n;

  m = 0;
  for (i = 0; i < 6; i++)
    for (j = i; j < 6; j++) {
      m++;

      n = 0;
      for (k = 0; k < 6; k++)
	for (l = k; l < 6; l++) {
	  n++;

	  if (k == l)
	    M[m][n]  = map[j][k]*map[i][l];
	  else
	    M[m][n]  = map[j][k]*map[i][l] + map[i][k]*map[j][l];

	  }
    }
}


void get_x(const ss_vect<tps> &Sigma, double *s)
{
  // Symmetric NxN matrix -> (N+1)xN/2 vector.
  int  i, j, n;

  n = 0;
  for (i = 0; i < 6; i++)
    for (j = i; j < 6; j++) {
      n++;
      s[n] = Sigma[i][j];
    }
}


void get_M(const double *s, ss_vect<tps> &M)
{
  // (N+1)xN/2 vector -> symmetric NxN matrix.
  int  i, j, n, jj[ss_dim];

  for (i = 0; i < ss_dim; i++)
    jj[i] = 0;

  n = 0;
  for (i = 0; i < 6; i++)
    for (j = i; j < 6; j++) {
      n++;
      jj[j] = 1; M[i].pook(jj, s[n]); jj[j] = 0;
      jj[i] = 1; M[j].pook(jj, s[n]); jj[i] = 0;
    }
}


void get_Sigma()
{
  long int      i, lastpos;
  int           j, k, jj[ss_dim];
  double        **M, **U, **V, *w, *s, *ds1, *ds0, norm[3];
  Vector        V1, V2;
  Matrix        S_mat;
  ss_vect<tps>  I, A, Sigma0, Sigma1, dSigma, map, map_tp, map2, S;

  const int     n_DOF = 3, n_ps = 2*n_DOF, n_lin = (n_ps+1)*n_ps/2;
  const double  s_cut = 1e-10;

  s = dvector(1, n_lin); ds0 = dvector(1, n_lin); ds1 = dvector(1, n_lin);
  M = dmatrix(1, n_lin, 1, n_lin); U = dmatrix(1, n_lin, 1, n_lin);
  V = dmatrix(1, n_lin, 1, n_lin); w = dvector(1, n_lin);

  I.identity();

  globval.Cavity_on = true;
  globval.radiation = true;

  Ring_GetTwiss(true, 0.0); printglob();

  getcod(0.0, lastpos);
  printf("\n");
  printf("COD: %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e\n",
	 globval.CODvect[x_], globval.CODvect[px_],
	 globval.CODvect[y_], globval.CODvect[py_],
	 globval.CODvect[delta_], globval.CODvect[ct_]);

  A.zero(); putlinmat(n_ps, globval.Ascr, A); A += globval.CODvect;

  Sigma0 = A*tp_S(n_DOF, A);

  prt_lin_map(n_DOF, Sigma0);

  // M Sigma M^T = (M (M Sigma)^T)^T = M (M Sigma)^T, since Sigma is symmetric.
  Sigma1 = Sigma0;
  for (i = 0; i <= globval.Cell_nLoc; i++) {
    Cell_Pass(i, i, Sigma1, lastpos);
    Sigma1 = tp_S(n_DOF, Sigma1);
    Cell_Pass(i, i, Sigma1, lastpos);
  }

  printf("\n");
  printf("By tracking:\n");
  prt_lin_map(3, Sigma1-Sigma0);

  map.identity(); map += globval.CODvect;
  Cell_Pass(0, globval.Cell_nLoc, map, lastpos);

  get_MM_TP_M(M, map);

  for (j = 1; j <= n_lin; j++)
    M[j][j] -= 1e0;

  // Transport matrix for Sigma is singular, i.e., gamma = (1+sqr(alpha))/beta.

  for (j = 1; j <= n_lin; j++)
    for (k = 1; k <= n_lin; k++)
      U[j][k] = M[j][k];

  dsvdcmp(U, n_lin, n_lin, w, V);

  printf("\n");
  printf("singular values:\n");
  printf("\n");
  for (j = 1; j <= n_lin; j++) {
    printf("%11.3e", w[j]);
    if (w[j] < s_cut) {
      w[j] = 0.0;
      printf(" (zeroed)");
    }
    if (j % 7 == 0) printf("\n");
  }
  if (n_lin % 7 != 0) printf("\n");

  Sigma0.identity();
  for (j = 1; j <= 1; j++) {
    // M Sigma M^T = (M (M Sigma)^T)^T = M (M Sigma)^T
    Sigma1 = Sigma0;
    for (i = 0; i <= globval.Cell_nLoc; i++) {
      Cell_Pass(i, i, Sigma1, lastpos);
      Sigma1 = tp_S(n_DOF, Sigma1);
      Cell_Pass(i, i, Sigma1, lastpos);
    }

    get_x(Sigma1-Sigma0, ds1);

    // Sigma0 -= (M-I, jj)^-1 * dSigma.
    dsvbksb(U, w, V, n_lin, n_lin, ds1, ds0);

    get_x(Sigma0, s); dvsub(s, n_lin, ds0, s); get_M(s, Sigma0);
  }
  get_M(s, Sigma0);

  // Normalize.
  S = get_S(n_DOF); getlinmat(n_ps, S, S_mat);

  for (k = 0; k < ss_dim; k++)
    jj[k] = 0;

  for (j = 0; j < n_DOF; j++) {
    for (k = 0; k < ss_dim; k++) {
      V1[k] = Sigma0[2*j][k]; V2[k] = Sigma0[2*j+1][k];
    }

    LinTrans(n_ps, S_mat, V2);

    norm[j] = 0e0;
    for (k = 0; k < ss_dim; k++)
      norm[j] += V1[k]*V2[k];

    norm[j] = sqrt(norm[j]);
  }

  for (j = 0; j < n_ps; j++)
    for (k = 0; k < n_DOF; k++) {
      jj[2*k] = 1; Sigma0[j].pook(jj, Sigma0[j][2*k]/norm[k]);
      jj[2*k] = 0; jj [2*k+1] = 1;
      Sigma0[j].pook(jj, Sigma0[j][2*k+1]/norm[k]);
      jj[2*k+1] = 0;
    }

  printf("\n");
  printf("Sigma0:\n");
  prt_lin_map(n_DOF, Sigma0);

  Sigma1 = Sigma0;
  for (i = 0; i <= globval.Cell_nLoc; i++) {
    Cell_Pass(i, i, Sigma1, lastpos);
    Sigma1 = tp_S(n_DOF, Sigma1);
    Cell_Pass(i, i, Sigma1, lastpos);
  }

  printf("\n");
  printf("Sigma1:\n");
//   prt_lin_map(n_DOF, Sigma1-Sigma0);
  prt_lin_map(n_DOF, Sigma1);

  free_dvector(s, 1, n_lin); free_dvector(ds0, 1, n_lin);
  free_dvector(ds1, 1, n_lin);
  free_dmatrix(M, 1, n_lin, 1, n_lin); free_dmatrix(U, 1, n_lin, 1, n_lin);
  free_dmatrix(V, 1, n_lin, 1, n_lin); free_dvector(w, 1, n_lin);
}


void FitTune2(double nux, double nuy)
{
  long      i, Fnum;
  iVector2  nq = {0, 0};
  Vector2   nu = {0.0, 0.0};
  fitvect   qfbuf, qdbuf;

  Fnum = ElemIndex("qh2sp01");
  for (i = 1; i <= GetnKid(Fnum); i++)
    qfbuf[nq[0]++] = Elem_GetPos(Fnum, i);

  Fnum = ElemIndex("qh2sp02");
  for (i = 1; i <= GetnKid(Fnum); i++)
    qfbuf[nq[0]++] = Elem_GetPos(Fnum, i);

  Fnum = ElemIndex("qh3sp01");
  for (i = 1; i <= GetnKid(Fnum); i++)
    qdbuf[nq[1]++] = Elem_GetPos(Fnum, i);

  Fnum = ElemIndex("qh3sp02");
  for (i = 1; i <= GetnKid(Fnum); i++)
    qdbuf[nq[1]++] = Elem_GetPos(Fnum, i);

  nu[X_] = nux; nu[Y_] = nuy;

  Ring_Fittune(nu, nueps, nq, qfbuf, qdbuf, nudkL, nuimax);
}


void set_lat(void)
{
  char    str1[max_str], str2[max_str];
  int     n, N, Fnum, qf, qd;
  double  nu[2];

  sprintf(str1, "%s%s", home_dir, "/projects/src/quad-K1_1.txt");
  sprintf(str2, "%s%s", home_dir, "/projects/src/quad-K1_2.txt");

  switch (wp) {
  case 3:
  case 4:
    n = 81;
    break;
  case 5:
    n = 60;
    break;
  case 6:
    n = 75;
    break;
  case 7:
    n = 97;
    break;
  case 8:
    n = 0;
    break;
  case 9:
    n = 0;
    break;
  case 10:
    n = 0;
    break;
  case 11:
    n = 0;
    break;
  case 12:
    n = 0;
    break;
  case 13:
    n = 0;
    break;
  default:
    std::cout << "set_lat: undef. working point" << std::endl;
    exit(1);
    break;
  }

  if (n != 0) {
    set_tune(str1, str2, n);
    Ring_GetTwiss(true, 0.0); printglob();
  }

  if (true) {
    N = 0;
    if (DW100) {
      Fnum = ElemIndex("b1"); N += GetnKid(Fnum);
    } else {
      Fnum = ElemIndex("b1sp01"); N += GetnKid(Fnum);
      Fnum = ElemIndex("b1sp02"); N += GetnKid(Fnum);
    }
    N /= 4;

    printf("\n");
    printf("N = %d\n", N);

    switch (wp) {
    case 3:
      // need to avoid nu_x + 4 nu_y = 96
      nu[X_] = N*33.47/15.0; nu[Y_] = N*15.68/15.0;
      break;
    case 5:
//    nu[X_] = N*33.15/15.0; nu[Y_] = N*15.69/15.0;
      // need to avoid nu_x + 4 nu_y = 96
      nu[X_] = N*33.15/15.0; nu[Y_] = N*15.68/15.0;
      break;
    case 6:
//    nu[X_] = N*33.15/15.0; nu[Y_] = N*16.3/15.0;
      // need to avoid 3 nu_x + 2 nu_y = 132
      nu[X_] = N*33.14/15.0; nu[Y_] = N*16.27/15.0;
      break;
    case 7:
//    nu[X_] = N*33.45/15.0; nu[Y_] = N*16.35/15.0;
      // test of working points
//    nu[X_] = N*33.47/15.0; nu[Y_] = N*16.33/15.0;
//    nu[X_] = N*33.43/15.0; nu[Y_] = N*16.335/15.0;
//    nu[X_] = N*33.46/15.0; nu[Y_] = N*16.39/15.0;
      // new working point
      nu[X_] = N*33.45/15.0; nu[Y_] = N*16.37/15.0;
//    nu[X_] = N*33.46/15.0; nu[Y_] = N*16.36/15.0;
//    nu[X_] = N*33.46/15.0; nu[Y_] = N*16.35/15.0;
      break;
    case 8:
//      nu[X_] = N*33.42/15.0; nu[Y_] = N*16.32/15.0;
//      nu[X_] = N*33.43/15.0; nu[Y_] = N*16.32/15.0;
//      nu[X_] = N*33.47/15.0; nu[Y_] = N*16.34/15.0;
      nu[X_] = N*33.47/15.0; nu[Y_] = N*16.33/15.0;
      break;
    case 9:
      nu[X_] = N*33.14/15.0; nu[Y_] = N*16.31/15.0;
      break;
    case 10:
//      nu[X_] = N*33.12/15.0; nu[Y_] = N*16.18/15.0;
//      nu[X_] = N*33.15/15.0; nu[Y_] = N*16.19/15.0;
//      nu[X_] = N*33.15/15.0; nu[Y_] = N*16.20/15.0;
      nu[X_] = N*33.16/15.0; nu[Y_] = N*16.19/15.0;
      break;
    case 11:
      // May, 2010
//      nu[X_] = N*33.12/15.0; nu[Y_] = N*16.19/15.0;
      // July, 2010
//      nu[X_] = N*33.08/15.0; nu[Y_] = N*16.19/15.0;
//      nu[X_] = N*33.08/15.0; nu[Y_] = N*(16.19-0.02)/15.0;
//      nu[X_] = N*33.08/15.0; nu[Y_] = N*16.19/15.0;
      // Sep, 2010
      // #1
//      nu[X_] = N*33.08/15.0; nu[Y_] = N*16.19/15.0;
      // #2
//      nu[X_] = N*33.15/15.0; nu[Y_] = N*16.22/15.0;
      // #3
//      nu[X_] = N*33.07/15.0; nu[Y_] = N*16.22/15.0;
      // Jan, 2011
//      nu[X_] = N*33.15/15.0; nu[Y_] = N*16.22/15.0;
//      nu[X_] = N*33.13/15.0; nu[Y_] = N*16.22/15.0;
//      nu[X_] = N*33.13/15.0; nu[Y_] = N*16.21/15.0;
      // Dec, 2011
      nu[X_] = N*33.11/15.0; nu[Y_] = N*16.2/15.0;
      // Aug, 2012
//      nu[X_] = N*33.10/15.0; nu[Y_] = N*16.18/15.0;
     break;
    case 12:
      // Baseline Lattice.
      nu[X_] = N*33.166/15.0; nu[Y_] = N*16.271/15.0;
      // Tweaked.
      // nu[X_] = N*33.18/15.0; nu[Y_] = N*16.271/15.0;
      break;
    case 13:
      // nu[X_] = N*33.21/15.0; nu[Y_] = N*16.28/15.0;
      // nu[X_] = N*33.21/15.0; nu[Y_] = N*16.35/15.0;
      // nu[X_] = N*33.21/15.0; nu[Y_] = N*16.34/15.0;
      // nu[X_] = N*33.21/15.0; nu[Y_] = N*16.335/15.0;
      // nu[X_] = N*33.21/15.0; nu[Y_] = N*16.32/15.0;
      // nu[X_] = N*33.19/15.0; nu[Y_] = N*16.31/15.0;
      // nu[X_] = N*33.23/15.0; nu[Y_] = N*16.33/15.0;

      // nu[X_] = N*33.23/15.0; nu[Y_] = N*16.31/15.0;

      // Crossing of skew sextupole resonance.
      // nu[X_] = N*33.17/15.0; nu[Y_] = N*16.35/15.0;
      // nu[X_] = N*33.2/15.0; nu[Y_] = N*16.32/15.0;
      // nu[X_] = N*33.26/15.0; nu[Y_] = N*16.42/15.0;
      nu[X_] = N*33.176/15.0; nu[Y_] = N*16.287/15.0;
      break;
    default:
      std::cout << "set_lat: undef. working point" << std::endl;
      exit(1);
      break;
    }

    if (DW100) {
      qf = ElemIndex("qh2"); qd = ElemIndex("qh3");

      FitTune(qf, qd, nu[X_], nu[Y_]);

      printf("\n");
      printf("nu_x = %7.5f nu_y = %7.5f\n",
	     globval.TotalTune[X_], globval.TotalTune[Y_]);
      printf("  qf = %8.5f   qd = %8.5f\n",
	     GetKpar(qf, 1, Quad), GetKpar(qd, 1, Quad));
    } else {
      FitTune2(nu[X_], nu[Y_]);

      printf("\n");
      printf("nu_x = %7.5f nu_y = %7.5f\n",
	     globval.TotalTune[X_], globval.TotalTune[Y_]);
      printf("  qh2sp01 = %8.5f   qh2sp02 = %8.5f\n",
	     GetKpar(ElemIndex("qh2sp01"), 1, Quad),
	     GetKpar(ElemIndex("qh2sp02"), 1, Quad));
      printf("  qh3sp01 = %8.5f   qh3sp02 = %8.5f\n",
	     GetKpar(ElemIndex("qh3sp01"), 1, Quad),
	     GetKpar(ElemIndex("qh3sp02"), 1, Quad));
    }

    Ring_GetTwiss(true, 0.0); printglob();
  }
}


void set_sext(void)
{
  const int  n_b3 = 9;

  char    str[max_str];
  int     k;
  double  b3s[n_b3];

  sprintf(str, "%s%s", home_dir, "/Thor-2.0/thor/wrk");

  no_sxt();

  switch (get_thor) {
  case 1:
    switch (wp) {
    case 3:
      if (!IDs)
	get_bn(strcat(str, "/sext_3.bare_3"), 0, true);
//      get_bn(strcat(str, "/sext_3.chrom_2"), 0, true);
//      get_bn(strcat(str, "/sext_3.chrom_2_geom_6"), 0, true);
//	get_bn(strcat(str, "/sext_3.chrom_2_geom_7"), 0, true);
//	get_bn(strcat(str, "/sext_3.ds_chrom_2_geom_7"), 0, true);
      else
	get_bn(strcat(str, "/sext_3.DW_3"), 0, true);
      break;
    case 5:
      if (!IDs)
	get_bn(strcat(str, "/sext_3.bare_5"), 0, true);
      else
	get_bn(strcat(str, "/sext_3.DW_5"), 0, true);
      break;
    case 6:
      if (!IDs)
	get_bn(strcat(str, "/sext_3.bare_6"), 0, true);
//	get_bn(strcat(str, "/sext_3.dat"), 0, true);
      else
	get_bn(strcat(str, "/sext_3.DW_6"), 0, true);
//	get_bn(strcat(str, "/sext_3.DW_6_ksi_1.0_4.0"), 0, true);
//	get_bn(strcat(str, "/sext_3.DW_6_ksi_1.0_4.0_dec"), 0, true);
      break;
    case 7:
      if (!IDs)
//	get_bn(strcat(str, "/sext_3.bare_7"), 0, true);
	get_bn(strcat(str, "/sext_3.bare_7_2"), 0, true);
      else
//	get_bn(strcat(str, "/sext_3.DW_7"), 0, true);
	get_bn(strcat(str, "/sext_3.DW_7.2"), 0, true);
      break;
    default:
      std::cout << "set_sext: undef. working point" << std::endl;
      exit(1);
      break;
    }
    break;
  case 2:
    switch (sext_scheme) {
    case 0:
      sprintf(str, "%s%s", home_dir, "/Thor-2.0/thor/wrk/Jun_2012/CDR/");
      break;
    case 1:
      sprintf(str, "%s%s", home_dir,
	      "/Thor-2.0/thor/wrk/Jun_2012/1_DW/ksi_2_2/sext/");
      break;
    case 2:
      sprintf(str, "%s%s", home_dir,
	      "/Thor-2.0/thor/wrk/Jun_2012/1_DW/ksi_2_2/oct/");
      break;
    case 3:
      sprintf(str, "%s%s", home_dir,
	      "/Thor-2.0/thor/wrk/Jun_2012/1_DW/ksi_2_2/b3_b4_b5/25_10_3.5"
	      "/res_1.0/");
      break;
    case 4:
      // sprintf(str, "%s%s", home_dir,
      // 	      "/Results/Commissioning/robust/bare/2/oct/");
      // sprintf(str, "%s%s", home_dir,
      // 	      "/Results/Commissioning/robust/3_DWs/oct/");
      sprintf(str, "%s", "/home/bengtsson/Results/Jan_2014/2nux-nuy/2/bare/");
      break;
    case 5:
      strcpy(str, "");
      break;
    }

    get_bn(strcat(str, "sext.dat"), 0, true);
    break;
  case 3:
    std::cout << "b3s?>";
    for (k = 0; k < n_b3; k++)
      scanf("%lf", &b3s[k]);
    set_bnL_design_fam(ElemIndex("sl1"), Sext, b3s[0], 0.0);
    set_bnL_design_fam(ElemIndex("sl2"), Sext, b3s[1], 0.0);
    set_bnL_design_fam(ElemIndex("sl3"), Sext, b3s[2], 0.0);
    set_bnL_design_fam(ElemIndex("sm1a"), Sext, b3s[3], 0.0);
    set_bnL_design_fam(ElemIndex("sm1b"), Sext, b3s[4], 0.0);
    set_bnL_design_fam(ElemIndex("sm2"), Sext, b3s[5], 0.0);
    set_bnL_design_fam(ElemIndex("sh1"), Sext, b3s[6], 0.0);
    set_bnL_design_fam(ElemIndex("sh2"), Sext, b3s[7], 0.0);
    set_bnL_design_fam(ElemIndex("sh3"), Sext, b3s[8], 0.0);
    break;
  default:
    std::cout << "set_sext: undef. mode" << std::endl;
    exit(1);
    break;
  }
}


void get_Dbn_rms_tol(const char *name, const int n, const double Dbn_max)
{
  char    str[max_str];
  bool    cod;
  int     k, Fnum;
  double  x_aper[n_aper], y_aper[n_aper], DA, Dbn, Dan;
  FILE    *DAf, *outf;

  const int     n_orbit = 3, n_scale = 1, n_step = 25;

  const char    home_dir[] = "/home/bengtsson/projects/in/";

  sprintf(ae_file, "%s%s", home_dir, "Align_Err/AlignErr.dat");
  sprintf(fe_file, "%s%s", home_dir, "Field_Err/tracy_baseline_2_2_10.txt");

  cod = CorrectCOD_N(ae_file, n_orbit, n_scale, 1);

  if (!cod) {
    printf("get_Dbn_rms_tol: orbit correction failed\n");
    exit(1);
  }

  LoadFieldErr(fe_file, false, 1.0, true);

  Fnum = ElemIndex(name);
  if (Fnum == 0) {
    printf("get_Dbn_rms_tol: undef. element %s\n", name);
    exit(1);
  }

  set_bnL_rms_fam(Fnum, abs(n), 0.0, 0.0, true);

  sprintf(str, "Dbn_rms_tol_%+d.%s", n, name);
  outf = file_write(str);

  sprintf(str, "dynap_rms_%+d.%s", n, name);
  DAf = file_write(str);

  printf("\n");
  globval.Cavity_on = true; n_track = 512;
  for (k = 0; k <= n_step; k++) {
    Dbn = (n > 0)? k*Dbn_max/n_step : 0.0;
    Dan = (n < 0)? k*Dbn_max/n_step : 0.0;
    set_bnL_rms_fam(Fnum, abs(n), Dbn, Dan, false);

    dynap(DAf, 5e-3, 0.0, 0.1e-3, n_aper, n_track,
	  x_aper, y_aper, false, true, false);
    DA = get_aper(n_aper, x_aper, y_aper);

    printf("DA = %10.3e [mm^2]\n", 1e6*DA);
    fprintf(outf, "%3d %10.3e %10.3e %10.3e\n", k, Dbn, Dan, 1e6*DA);
    fflush(outf);
  }

  fclose(DAf); fclose(outf);
}


void get_Dbn_sys_tol(const char *name, const int n, const double Dbn_max)
{
  char    str[max_str];
  bool    cod;
  int     k, Fnum;
  double  x_aper[n_aper], y_aper[n_aper], DA, Dbn, Dan;
  FILE    *DAf, *outf;

  const int     n_orbit = 3, n_scale = 1, n_step = 25;

  const char    home_dir[] = "/home/bengtsson/projects/in/";

  sprintf(ae_file, "%s%s", home_dir, "Align_Err/AlignErr.dat");
  sprintf(fe_file, "%s%s", home_dir, "Field_Err/tracy_baseline_2_2_10.txt");

  cod = CorrectCOD_N(ae_file, n_orbit, n_scale, 1);

  if (!cod) {
    printf("get_Dbn_sys_tol: orbit correction failed\n");
    exit(1);
  }

  LoadFieldErr(fe_file, false, 1.0, true);

  Fnum = ElemIndex(name);
  if (Fnum == 0) {
    printf("get_Dbn_sys_tol: undef. element %s\n", name);
    exit(1);
  }

  set_bnL_sys_fam(Fnum, abs(n), 0.0, 0.0);

  sprintf(str, "Dbn_sys_tol_%+d.%s", n, name);
  outf = file_write(str);

  sprintf(str, "dynap_sys_%+d.%s", n, name);
  DAf = file_write(str);

  printf("\n");
  globval.Cavity_on = true; n_track = 512;
  trace = true;
  for (k = 0; k <= n_step; k++) {
    Dbn = (n > 0)? k*Dbn_max/n_step : 0.0;
    Dan = (n < 0)? k*Dbn_max/n_step : 0.0;
    set_bnL_sys_fam(Fnum, abs(n), Dbn, Dan);

    dynap(DAf, 5e-3, 0.0, 0.1e-3, n_aper, n_track,
	  x_aper, y_aper, false, true, false);
    DA = get_aper(n_aper, x_aper, y_aper);

    printf("DA = %10.3e [mm^2]\n", 1e6*DA);
    fprintf(outf, "%3d %10.3e %10.3e %10.3e\n", k, Dbn, Dan, 1e6*DA);
    fflush(outf);
  }

  fclose(DAf); fclose(outf);
}


void get_ID_tol(const int i)
{

  n_aper = 15;
  switch (i) {
  case 1:
    get_Dbn_rms_tol("dw_90", 3, 1.0);
    break;
  case 2:
    get_Dbn_rms_tol("dw_90", -3, 1.0);
    break;
  case 3:
    get_Dbn_rms_tol("dw_90", 4, 25.0);
    break;
  case 4:
    get_Dbn_rms_tol("dw_90", -4, 25.0);
    break;

  case 5:
    get_Dbn_rms_tol("ivu_20", 3, 4.0);
    break;
  case 6:
    get_Dbn_rms_tol("ivu_20", -3, 4.0);
    break;
  case 7:
    get_Dbn_rms_tol("ivu_20", 4, 100.0);
    break;
  case 8:
    get_Dbn_rms_tol("ivu_20", -4, 100.0);
    break;

  case 9:
    get_Dbn_rms_tol("dpwg5b", 3, 1.0);
    break;
  case 10:
    get_Dbn_rms_tol("dpwg5b", -3, 1.0);
    break;
  case 11:
    get_Dbn_rms_tol("dpwg5b", 4, 100.0);
    break;
  case 12:
    get_Dbn_rms_tol("dpwg5b", -4, 100.0);
    break;

  case 13:
    get_Dbn_sys_tol("dpwg5b", 3, 1.0);
    break;
  case 14:
    get_Dbn_sys_tol("dpwg5b", -3, 1.0);
    break;
  case 15:
    get_Dbn_sys_tol("dpwg5b", 4, 100.0);
    break;
  case 16:
    get_Dbn_sys_tol("dpwg5b", -4, 100.0);
    break;
  };
}


void get_sxt_rnd(const char file_name[], const int n)
{
  const int  line_max = 200, n_prm_max = 20;

  char      line[line_max];
  int       j;
  double    b3L[n_prm_max], ksi2 = 0.0;
  tps       K_re, K_im;
  std::ifstream  inf;

  file_rd(inf, file_name);

  // skip blank line
  inf.getline(line, line_max);
  while (!inf.getline(line, line_max).eof()) {
    sscanf(line,
	   "%d %*s = %lf %*s = %lf %*s = %lf %*s = %lf %*s = %lf"
	   " %*s = %lf %*s = %lf %*s = %lf %*s = %lf %*s = %lf",
	   &j, &ksi2, &b3L[0], &b3L[1], &b3L[2], &b3L[3], &b3L[4],
	   &b3L[5], &b3L[6], &b3L[7], &b3L[8]);

    if (j == n) {
      printf("\n");
      printf("%3d %7.1e %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f"
	     " %8.5f\n",
	     j, ksi2, b3L[0], b3L[1], b3L[2], b3L[3], b3L[4], b3L[5],
	     b3L[6], b3L[7], b3L[8]);

      set_bnL_design_fam(ElemIndex("sl1"), Sext, b3L[0], 0.0);
      set_bnL_design_fam(ElemIndex("sl2"), Sext, b3L[1], 0.0);
      set_bnL_design_fam(ElemIndex("sl3"), Sext, b3L[2], 0.0);
      set_bnL_design_fam(ElemIndex("sh1"), Sext, b3L[3], 0.0);
      set_bnL_design_fam(ElemIndex("sh3"), Sext, b3L[4], 0.0);
      set_bnL_design_fam(ElemIndex("sh4"), Sext, b3L[5], 0.0);
      set_bnL_design_fam(ElemIndex("sm1a"), Sext, b3L[6], 0.0);
      set_bnL_design_fam(ElemIndex("sm1b"), Sext, b3L[7], 0.0);
      set_bnL_design_fam(ElemIndex("sm2"), Sext, b3L[8], 0.0);

      break;
    }

    inf.getline(line, line_max);
  }
}


void get_num_map(void)
{
  long int         lastpos;
  int              i, j;
  double           x, y;
  ss_vect<double>  ps;
  ss_vect<tps>     map;
  std::ofstream         outf;

  const int     nx = 15, ny = 6;
  const double  dx = 1e-3, dy = 0.5e-3, x_max = nx*dx, y_max = ny*dy;


  map.identity(); Cell_Pass(0, globval.Cell_nLoc, map, lastpos);

  std::cout << std::endl;
  std::cout << std::fixed << std::setprecision(5)
       << "L [m] = " << std::setw(7) << map[ct_].cst() << std::endl;

  prt_lin_map(3, map);

  file_wr(outf, "num_map_x.dat");

  x = -x_max;
  for (i = -nx; i <= nx; i++) {
    ps.zero(); ps[x_] = x;
    outf << std::scientific << std::setprecision(5) << std::setw(13) << ps[x_];
    Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
    outf << std::scientific << std::setprecision(5) << std::setw(13) << ps;
    x += dx;
  }

  outf.close();

  file_wr(outf, "num_map_y.dat");

  y = -y_max;
  for (i = -ny; i <= ny; i++) {
    ps.zero(); ps[y_] = y;
    outf << std::scientific << std::setprecision(5) << std::setw(13) << ps[y_];
    Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
    outf << std::scientific << std::setprecision(5) << std::setw(13) << ps;
    y += dy;
  }

  outf.close();

  file_wr(outf, "num_map_xy.dat");

  x = -x_max;
  for (i = -nx; i <= nx; i++) {
    y = -y_max;
    for (j = -ny; j <= ny; j++) {
      ps.zero(); ps[x_] = x; ps[y_] = y;
      outf << std::scientific << std::setprecision(5)
	   << std::setw(13) << ps[x_] << std::setw(13) << ps[y_];
      Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
      outf << std::scientific << std::setprecision(5)
	   << std::setw(13) << ps[x_]-x << std::setw(13) << ps[px_]
	   << std::setw(13) << ps[y_]-y << std::setw(13) << ps[py_] << std::endl;
      y += dy;
    }
    outf << std::endl;
    x += dx;
  }

  outf.close();
}


void get_kick_map(const char *file_name, const int nx, const int ny,
		  const long int loc1, const long int loc2,
		  const double Ax, const double Ay)
{
  long int         lastpos;
  int              i1, i2;
  ss_vect<double>  ps0, ps1, dps_map;
  std::ofstream         outf;

  const double  Brho = globval.Energy*1e9/c0;

  file_wr(outf, file_name);

  dps_map[x_] = 2.0*Ax/(nx-1.0); dps_map[y_] = 2.0*Ay/(ny-1.0);

  outf << "# Author:" << std::endl;
  outf << "# Title" << std::endl;
  outf << "# Cell Length [m]" << std::endl;
  outf << std::fixed << std::setprecision(5) << Cell[loc2].S-Cell[loc1].S << std::endl;
  outf << "# Number of Horizontal Points" << std::endl;
  outf << nx << std::endl;
  outf << "# Number of Vertical Points" << std::endl;
  outf << ny << std::endl;

  outf << "# Horizontal 2nd Order Kick [T2m2]" << std::endl;
  outf << "START" << std::endl;

  for (i1 = 0; i1 < nx; i1++)
    outf << std::scientific << std::setprecision(5) << std::setw(13) << -Ax+i1*dps_map[x_];
  outf << std::endl;

  std::cout << std::endl;
  ps0.zero();
  for (i1 = 0; i1 < ny; i1++) {
    ps0[y_] = Ay - i1*dps_map[y_];

    std::cout << ".";

    outf << std::scientific << std::setprecision(5)
	 << std::setw(13) << ps0[y_];

    for (i2 = 0; i2 < nx; i2++) {
    ps0[x_] = -Ax + i2*dps_map[x_];

      ps1 = ps0;
      Cell_Pass(loc1, loc2, ps1, lastpos);

      if (lastpos == loc2) {
	ps1 -= ps0;

	outf << std::scientific << std::setprecision(5)
	     << std::setw(13) << sqr(Brho)*ps1[px_];
      } else
	outf << std::scientific << std::setprecision(5) << std::setw(13) << NAN;
    }

    outf << std::endl;
  }
  std::cout << std::endl;

  outf << "# Vertical 2nd Order Kick [T2m2]" << std::endl;
  outf << "START" << std::endl;

  for (i1 = 0; i1 < nx; i1++)
    outf << std::scientific << std::setprecision(5) << std::setw(13) << -Ax+i1*dps_map[x_];
  outf << std::endl;

  std::cout << std::endl;
  ps0.zero();
  for (i1 = 0; i1 < ny; i1++) {
    ps0[y_] = Ay - i1*dps_map[y_];

    std::cout << ".";

    outf << std::scientific << std::setprecision(5)
	 << std::setw(13) << ps0[y_];

    for (i2 = 0; i2 < nx; i2++) {
    ps0[x_] = -Ax + i2*dps_map[x_];

      ps1 = ps0;
      Cell_Pass(loc1, loc2, ps1, lastpos);

      if (lastpos == loc2) {
	ps1 -= ps0;

	outf << std::scientific << std::setprecision(5)
	     << std::setw(13) << sqr(Brho)*ps1[py_];
      } else
	outf << std::scientific << std::setprecision(5) << std::setw(13) << NAN;
    }

    outf << std::endl;
  }
  std::cout << std::endl;

  outf.close();
}


void get_map_2D(const char *file_name, const int nx, const int ny,
		const long int loc1, const long int loc2,
		const double Ax, const double Ay)
{
  long int         lastpos;
  int              i1, i2;
  ss_vect<double>  ps0, ps1, dps_map;
  std::ofstream         outf;

  file_wr(outf, file_name);

  dps_map[x_] = 2.0*Ax/(nx-1.0); dps_map[y_] = 2.0*Ay/(ny-1.0);

  outf << std::scientific << std::setprecision(5)
       << "# nx = " << nx << "     Ax = " << Ax
       << "     dx = " << dps_map[x_] << std::endl;
  outf << std::scientific << std::setprecision(5)
       << "# ny = " << ny << "     Ay = " << Ay
       << "     dy = " << dps_map[y_] << std::endl;
  outf << "#" << std::endl;

  std::cout << std::endl;
  ps0.zero();
  for (i1 = 0; i1 < nx; i1++) {
    ps0[x_] = -Ax + i1*dps_map[x_];
    for (i2 = 0; i2 < ny; i2++) {
      ps0[y_] = -Ay + i2*dps_map[y_];

//       cout << std::setw(3) << i1+1 << std::setw(3) << i2+1
// 	   << scientific << setprecision(5)
// 	   << std::setw(13) << ps0[x_] << std::setw(13) << ps0[y_] << endl;

      ps1 = ps0;
      Cell_Pass(loc1, loc2, ps1, lastpos);

      ps1 -= ps0;

      outf << std::scientific << std::setprecision(5)
	   << std::setw(13) << ps0[x_] << std::setw(13) << ps0[y_]
	   << std::setw(13) << ps1[px_] << std::setw(13) << ps1[py_] << std::endl;
    }
    outf << std::endl;
  }

  outf.close();
}


ss_vect<tps> get_fix_point(const int Fnum, const int Knum)
{
  long int     loc;
  int          j;
  FieldMapType *FM;
  std::ofstream     outf;

  const int   n_iter = 3;

  trace = true;

  loc = Elem_GetPos(Fnum, Knum); FM = Cell[loc].Elem.FM;

  Cell[loc].Elem.PL = 2.62; FM->phi = 6.0*M_PI/48.0; FM->n_step = 10;

  FM->Ld = 0.0; FM->L1 = 0.0;
  for (j = 1; j <= n_iter; j++) {
    if (j == n_iter) file_wr(outf, "track.dat");

    FieldMap_Pass(Cell[loc], map);

    if (j == n_iter) outf.close();

    FM->phi -= map[px_].cst(); FM->Lm = map[ct_].cst();
    FM->Ld += 2.0*map[x_].cst()/FM->phi; FM->L1 += Cell[loc].Elem.PL - FM->Lm;
  }

  FieldMap_Pass(Cell[loc], map);

  std::cout << std::endl;
  std::cout << std::scientific << std::setprecision(3)
       << "get_fix_point:" << std::setw(11) << map.cst();
  std::cout << std::fixed << std::setprecision(5)
       << "phi [deg] = " << std::setw(7) << FM->phi*180.0/M_PI
       << ", L [m] = " << std::setw(7) << Cell[loc].Elem.PL << std::endl;
  std::cout << std::fixed << std::setprecision(5)
       << "Lr [m] = " << std::setw(7) << FM->Lr
       << ", Lm [m] = " << std::setw(7) << FM->Lm
       << ", Ld [m] = " << std::setw(7) << FM->Ld
       << ", L1 [m] = " << std::setw(7) << FM->L1 << std::endl;

  return map;
}


void err_and_corr(const char *param_file)
{
  bool      cod, cav, rad, aper;

  // save state for dynamics
  cav = globval.Cavity_on; rad = globval.radiation; aper = globval.Aperture_on;

  // set state for correction
  globval.Cavity_on   = false; globval.radiation = false;
  globval.Aperture_on = false;

  get_param(param_file);

  Ring_GetTwiss(true, 0.0); printglob();

  prt_lat("linlat.out", globval.bpm, true);

  get_bare();

//   if (N_calls > 0) ini_ID_corr(false);

  if (strcmp(ae_file, "") != 0) {
    gcmat(1); gcmat(2);
  }

  // Coupling correction and dispersion wave initialization.
  // if (n_lin > 0) ini_skew_cor(disp_wave_y);

  if (strcmp(ae_file, "") != 0) {
    cod = CorrectCOD_N(ae_file, n_orbit, 1, 1);
    printf("\n");
    if (cod)
      printf("orbit correction completed\n");
    else
      chk_cod(cod, "error_and_correction");
  }

  GetEmittance(ElemIndex("cav"), true);

  // if (n_lin > 0) corr_eps_y();

  if (strcmp(ap_file, "") != 0) LoadApers(ap_file, 1.0, 1.0);

  if (strcmp(fe_file, "") != 0) LoadFieldErr(fe_file, false, 1.0, true);

  Ring_GetTwiss(true, 0.0); printglob();

  GetEmittance(ElemIndex("cav"), true);

  prt_beamsizes();

  prt_lat("linlat_err.out", globval.bpm, true);
  prtmfile("flat_file_err.dat");

  globval.Cavity_on = cav; globval.radiation = rad; globval.Aperture_on = aper;
}


void get_oc(std::string file_name, std::string names[])
{
  char     line[max_str], name[max_str];
  int      k;
  std::ifstream inf;

  file_rd(inf, file_name.c_str());
  k = 0;
  while (inf.getline(line, max_str)) {
    k++;
    sscanf(line, "%s", name);
    names[k-1] = name;
  }
  inf.close();
  std::cout << "get_oc: read " << k << " names" << std::endl;
}


void lwr_case(std::string &str)
{
  int k;

  for (k = 0; k < (int)str.length(); k++)
    str[k] = tolower(str[k]);
}


void get_mpoles(std::ofstream &fp, const std::string &fam_name)
{
  long int k;
  int      k_num, sm1;
  std::string   name, name_lwr;
  double   bn, an;

  k_num = 0;
  for (k = 0; k < globval.Cell_nLoc; k++) {
    if (Cell[k].Elem.Pkind == Mpole) {
      name = std::string(Cell[k].Elem.PName);
      name_lwr = name; lwr_case(name_lwr);
      if (k_num == 0) {
	if (fam_name.compare("sm1") == 0)
	  sm1 = 1;
	else if (fam_name.compare("sm3") == 0)
	  sm1 = 3;
	else
	  sm1 = 0;
      }

      if (((sm1 == 0) && (name_lwr.compare(0, 3, fam_name) == 0)) ||
	  ((sm1 == 1) && (name_lwr.compare(0, 3, "sm1") == 0) &&
	   (name_lwr[8] == 'a')) ||
	  ((sm1 == 3) && (name_lwr.compare(0, 3, "sm1") == 0) &&
	   (name_lwr[8] == 'b'))) {

	k_num++;
	get_bn_design_elem(Cell[k].Fnum, Cell[k].Knum,
			   Cell[k].Elem.M->n_design, bn, an);
	fp << std::scientific << std::setprecision(16) << std::setw(11) << name;
	if (sm1 == 0)
	  fp << std::setw(4) << fam_name;
	else if (sm1 == 1)
	  fp << std::setw(4) << "sm1a";
	if (sm1 == 3)
	  fp << std::setw(4) << "sm1b";

	fp << std::setw(4) << k_num << std::setw(3) << Cell[k].Elem.M->n_design
	   << std::setw(24) << bn << std::endl;
      }
    }
  }
  fp << "#" << std::endl;
}


void get_mpoles(std::string file_name)
{
  std::ofstream fp;

  fp.open(file_name.c_str(), std::ios::out);

  get_mpoles(fp, "qh1"); get_mpoles(fp, "qh2"); get_mpoles(fp, "qh3");
  get_mpoles(fp, "qm1"); get_mpoles(fp, "qm2");
  get_mpoles(fp, "ql1"); get_mpoles(fp, "ql2"); get_mpoles(fp, "ql3");
  get_mpoles(fp, "sh1"); get_mpoles(fp, "sh3"); get_mpoles(fp, "sh4");
  get_mpoles(fp, "sm1"); get_mpoles(fp, "sm3");
  get_mpoles(fp, "sm2");
  get_mpoles(fp, "sl1"); get_mpoles(fp, "sl2"); get_mpoles(fp, "sl3");

  fp.close();
}


void set_mpoles(std::string file_name)
{
  std::string       line, name, fam_name;
  int          n, f_num, k_num;
  double       bn, an;
  std::stringstream str;
  std::ifstream     inf;

  const bool prt = false;

  printf("\n");
  printf("set_mpoles\n");

  inf.open(file_name.c_str(), std::ios::in);

  an = 0e0;
  while (getline(inf, line)) {
    if (line.compare(0, 1, "#") != 0) {
      str.clear(); str.str(""); str << line;
      str >> name >> fam_name >> k_num >> n >> bn;
      f_num = ElemIndex(fam_name);
      if (prt)
	printf("%10s %5s %3d %3d %11.3e %3d\n",
	       name.c_str(), fam_name.c_str(), k_num, n, bn, GetnKid(f_num));
      set_bn_design_elem(f_num, k_num, n, bn, an);
    } else
      if (prt) printf("\n");
  }

  inf.close();
}


void ps_bend(const char *bend, const int n,
	     const double x, const double px, const double y, const double py,
	     const double delta)
{
  long int        loc, lastpos;
  int             k;
  ss_vect<double> ps;
  std::ofstream        outf1, outf2;

  no_sxt();

  outf1.open("ps_bend_1.out", std::ios::out);
  outf2.open("ps_bend_2.out", std::ios::out);

  loc = Elem_GetPos(ElemIndex(bend), 1);

  ps[x_] = x; ps[px_] = px; ps[y_] = y; ps[py_] = py;
  ps[delta_] = delta; ps[ct_] = 0e0;

  for (k = 1; k <= n; k++) {
    Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
    outf1 << std::scientific << std::setprecision(5)
	  << std::setw(5) << k << std::setw(13) << ps << std::endl;
    outf2 << std::scientific << std::setprecision(5)
	  << std::setw(5) << k << std::setw(13) << Cell[loc].BeamPos << std::endl;
  }

  outf1.close();
  outf2.close();
}


int main(int argc, char *argv[])
{
  char          ch, line[max_str], fam_name[max_str];
  long int      loc, lastn, lastpos, n;
  int           i, k, Fnum, Knum, prm_type, sf, sd;
  double        bn, an, bn0, an0, b3[2], a3, ps[4][1024], nu[2][1], A[2][1];
  Matrix        M;
  ss_vect<tps>  map;

  const long  seed = 1121;

#define OC 3
#if OC == 1
  const int   n_bpm_Fam = 1, n_hcorr_Fam = 1, n_vcorr_Fam = 1;

  const std::string  bpm_names[n_bpm_Fam]     = { "BPM" };
  const std::string  hcorr_names[n_hcorr_Fam] = { "CHV" };
  const std::string  vcorr_names[n_vcorr_Fam] = { "CHV" };
#elif OC == 2
  const int   n_bpm_Fam = 5, n_hcorr_Fam = 9, n_vcorr_Fam = 9;

  const std::string  bpm_names[n_bpm_Fam] =
    { "PH1", "PH2", "PM1", "PL2", "PL1" };

  const std::string  hcorr_names[n_hcorr_Fam] =
    { "CXSQ1H", "CXH2", "CXM1", "CXM1G4C30B", "CXL2",
      "CXL1", "CXSQ2H", "CXM1G4C01B", "CXH1" };

  const std::string  vcorr_names[n_vcorr_Fam] =
    { "CYSQ1H", "CYH2",   "CYM1",   "CYM1G4C30B", "CYL2",
      "CYL1",   "CYSQ2H", "CYM1G4C01B", "CYH1" };
#elif OC == 3
  const int   n_bpm_Fam = 12, n_hcorr_Fam = 12, n_vcorr_Fam = 12;

  const std::string  bpm_names[n_bpm_Fam] =
    { "PH1G2C30A", "PH2G2C30A", "PM1G4C30A", "PM1G4C30B", "PL2G6C30B",
      "PL1G6C30B", "PL1G2C01A", "PL2G2C01A", "PM1G4C01A", "PM1G4C01B",
      "PH2G6C01B", "PH1G6C01B" };

  const std::string  hcorr_names[n_hcorr_Fam] =
    { "PH1G2C30A", "PH2G2C30A", "PM1G4C30A", "PM1G4C30B", "PL2G6C30B",
      "PL1G6C30B", "PL1G2C01A", "PL2G2C01A", "PM1G4C01A", "PM1G4C01B",
      "PH2G6C01B", "PH1G6C01B" };

  const std::string  vcorr_names[n_vcorr_Fam] =
    { "PH1G2C30A", "PH2G2C30A", "PM1G4C30A", "PM1G4C30B", "PL2G6C30B",
      "PL1G6C30B", "PL1G2C01A", "PL2G2C01A", "PM1G4C01A", "PM1G4C01B",
      "PH2G6C01B", "PH1G6C01B" };
#elif OC == 4
  const int n_bpm_Fam = 180, n_hcorr_Fam = 180, n_vcorr_Fam = 180;

  std::string bpm_names[n_bpm_Fam];
  std::string hcorr_names[n_hcorr_Fam], vcorr_names[n_vcorr_Fam];
#endif

  iniranf(seed); setrancut(1.0);

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  // 3: Oleg I, 4: Oleg II
  FieldMap_filetype = 4; sympl = true;

  freq_map = false;

  n_x = 50; n_y = 30; n_dp = 25; n_tr = 2064;
  x_max_FMA = 20e-3; y_max_FMA = 6e-3; delta_FMA = 3e-2;

  // disable from TPSALib- and LieLib log messages
//  idprset(-1);

  if (true)
    Read_Lattice(argv[1]);
  else {
    globval.Energy = 3e0;
    rdmfile(argv[1]);
  }

  Ring_GetTwiss(true, 0.0); printglob();

  if (false) {
    get_mpoles("sr1115_mpoles.txt");

    exit(0);
  }

  if (false) set_mpoles("sr1115_mpoles.txt");

  if (false) {
    // Design roll: 90 degrees; for thick multipoles only.
    set_roll_design_fam(ElemIndex("sm3a"), 90e0);
    set_roll_design_fam(ElemIndex("sm3b"), 90e0);
    set_roll_design_fam(ElemIndex("sm3c"), 90e0);
    set_roll_design_fam(ElemIndex("sm3d"), 90e0);
  }

  if (get_thor != 0) sscanf(argv[2], "%d", &sext_scheme);

//  get_sxt_rnd("/home/bengtsson/Thor-2.0/thor/wrk/sext_rnd.dat", 1360);
//  get_sxt_rnd("/home/bengtsson/Thor-2.0/thor/wrk/sext_rnd.dat", 5805);
//  get_sxt_rnd("/home/bengtsson/Thor-2.0/thor/wrk/sext_rnd.dat", 4653);
//  get_sxt_rnd("/home/bengtsson/Thor-2.0/thor/wrk/sext_rnd.dat", 3873);
//  get_sxt_rnd("/home/bengtsson/Thor-2.0/thor/wrk/sext_rnd.dat", 3593);
//  get_sxt_rnd("/home/bengtsson/Thor-2.0/thor/wrk/sext_rnd.dat", 2490);

  globval.EPU = false;

  if (false) {
    loc = Elem_GetPos(ElemIndex("qm1"), 1);

    Cell[loc].Elem.M->PdTpar = 0.5e-3;
    Mpole_SetdT(Cell[loc].Fnum, Cell[loc].Knum);
  }

  if (false) {
    map = get_fix_point(ElemIndex("B1"), 1);

    prt_lin_map(3, map);

    getlinmat(6, map, M);
    std::cout << std::endl;
    std::cout << std::scientific << std::setprecision(3)
	 << "1-Det: " << std::setw(9) << 1-DetMat(6, M) << std::endl;

    exit(0);
  }

  if (globval.H_exact) bend_cal();

  if (ID) {
    if (DW100)
      IDs = GetnKid(ElemIndex("dw100")) > 0;
    else
      IDs = GetnKid(ElemIndex("dw90")) > 0;
//    if (!IDs) IDs = GetnKid(ElemIndex("ivu_20_g59")) > 0;
//    if (!IDs) IDs = GetnKid(ElemIndex("u55")) > 0;
  } else
    IDs = false;

  printf("\n");
  printf("IDs = %d\n", IDs);

  if (IDs) {
    get_IDs(); set_IDs(0.0);

    Ring_GetTwiss(true, 0.0); printglob();
  }

  // no_sxt();

  if (wp != 0) {
    set_lat();

    Ring_GetTwiss(true, 0.0); printglob();
  }

  if (sext_scheme != 0) {
    set_sext();

    Ring_GetTwiss(true, 0.0); printglob();
  }

  if (false) {
    sf = ElemIndex("sf"); sd = ElemIndex("sd");
    FitChrom(sf, sd, 0e0, 0e0);
    get_bn_design_elem(sf, 1, Sext, b3[0], a3);
    get_bn_design_elem(sd, 1, Sext, b3[1], a3);
    std::cout << std::endl;
    std::cout << std::fixed << std::setprecision(5)
	 << "sf = " << b3[0] << ", sd = " << b3[1] << std::endl;

    Ring_GetTwiss(true, 0.0); printglob();
  }

  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);

  prt_chrom_lat();

//   prt_ZAP(15);

  if (false) {
#if OC == 4
  cout << endl;
  get_oc("/home/bengtsson/Yoshi/lattice/bpm_names.txt", bpm_names);
  get_oc("/home/bengtsson/Yoshi/lattice/hcor_names.txt", hcorr_names);
  get_oc("/home/bengtsson/Yoshi/lattice/vcor_names.txt", vcorr_names);
#endif

    ini_COD_corr(n_bpm_Fam, bpm_names,
		 n_hcorr_Fam, hcorr_names, n_vcorr_Fam, vcorr_names, true);

    prt_gcmat(1); prt_gcmat(2);
  }

  get_bare();

  if (IDs) {
    set_IDs(1.0);

    Ring_GetTwiss(true, 0.0); printglob();

    prt_lat("linlat_IDs.out", globval.bpm, true);

    ID_correct();

    Ring_GetTwiss(true, 0.0); printglob();

    prt_lat("linlat_IDs_corr.out", globval.bpm, true);
    prt_dlat("dlinlat.out", true);
  }

  if (false) {
    err_and_corr("param.dat");

    exit(0);
  }

  Ring_GetTwiss(true, 0.0); printglob();

  prtmfile("flat_file.dat");

  GetEmittance(ElemIndex("cav"), true);

  // ps_bend("b1", 512, 15e-3, 0e0, 5.75e-3, 0e0, 0e0);

//  sscanf(argv[2], "%d", &i); get_ID_tol(i);

  if (false) {
    std::cout << std::endl;
    std::cout << "computing kick map" << std::endl;

//     get_num_map();

    get_kick_map("super_cell.dat",
		 51, 51, 0, globval.Cell_nLoc, 20e-3, 6.5e-3);
    get_map_2D("super_cell_2D.dat",
	       51, 51, 0, globval.Cell_nLoc, 20e-3, 6.5e-3);
  }

  if (tune_shift) {
    std::cout << std::endl;
    std::cout << "computing tune shifts" << std::endl;
    dnu_dA(20e-3, 7.5e-3, 0.0, 25); get_ksi2(3.0e-2);
  }

  if (freq_map) {
    fmap(n_x, n_y, n_tr, x_max_FMA, y_max_FMA, 0.0, true, true);
    fmapdp(n_x, n_dp, n_tr, x_max_FMA, -delta_FMA, 1e-3, true, true);
  }

//  chk_IBS();

//   get_Touschek();

//   get_Sigma();

  if (false) {
    Fnum = ElemIndex("ivu_20");

    for (i = 1; i <= GetnKid(Fnum); i++)
      set_ID_scl(Fnum, i, ranf());

    if (false) {
      track("track.out", 8e-3, 0e0, 1.5e-3, 0e0, 0e0, 512, lastn, lastpos,
	    0, 0.0);

      GetTrack("track.out", &n, ps[x_], ps[px_], ps[y_], ps[py_]);
      printf("\n");
      printf("Read %ld turns\n", n);
      sin_FFT(n, ps[x_]); sin_FFT(n, ps[y_]);
      GetPeaks(n, ps[x_], 1, nu[X_], A[X_]);
      GetPeaks(n, ps[y_], 1, nu[Y_], A[Y_]);
      printf("\n");
      printf("nu_x = %7.5f, 1-nu_x = %7.5f, nu_y = %7.5f, 1-nu_y = %7.5f\n",
	     nu[X_][0], 1.0-nu[X_][0], nu[Y_][0], 1.0-nu[Y_][0]);
    }

    // globval.Cavity_on = true;
    // globval.radiation = true;

    // 10 IVUs, same strengths
    //   track("track.out", 6.92e-3, 0e0, 1.5e-3, 0e0, 0e0, 512, lastn, lastpos,
    // 	    0, 0.0);
    // 10 IVUs, random strengths
    // track("track.out", 10e-3, 0e0, 1.5e-3, 0e0, 0e0, 20000, lastn, lastpos,
    // 	  0, 0.0);
  }

//  Ring_GetTwiss(true, 0.0);
//  get_fixed_points(40e-3, 25e-3, 25);
  // grid size for DW kick map
//  get_fixed_points(30e-3, 5e-3, 25);
//  get_lin_inv(5e-3, 5e-3, 0e-2, 512);
//  get_sympl_form(10e-3, 6e-3, 0e-2, 1024);

//  get_phi_stat(40e-3, 25e-3, 25);

  if (get_DA) {
    if (false) {
      err_and_corr("param.dat");

      GetEmittance(ElemIndex("cav"), true);
    }

    globval.Cavity_on = true;
    n_aper = 25;
    get_dynap(delta, true);
  }
}
