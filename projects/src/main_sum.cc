#define NO 1

#include "tracy_lib.h"

int no_tps = NO;

extern bool  freq_map;
extern int   N_Fam, Q_Fam[];

// dynamic aperture
const double  delta  = 2.5e-2; // delta for off-momentum aperture

double  R0 = 25e-3;

const int     N_cell = 15;
// for quadrupoles (even in delta)
const double  Ax_tol = 15e-3, delta_tol = 2.5e-2, dnu_x_tol = 0.005;
// for sextupoles (odd in delta)
//const double  Ax_tol = 15e-3, delta_tol = 2.5e-2, dnu_x_tol = 0.0025;

//const int     n_orbit = 3;

#if ORDER > 1
const double  max_Ax = 20e-3, max_Ay = 10e-3, max_delta = 2.5e-2;
#endif

const bool tune_scn = false;               ;

//const int wp = 7;
const int wp = 0;
const int get_thor = 0;

const int  max_Fnum = 10, max_quads = 200, max_loc = 10;

bool      IDs;
int       Fnum[max_Fnum], loc[max_loc], n_iter, n_prm, n_loc;
int       quads1[max_quads], quads2[max_quads];
double    alpha_mp[2][2], beta_mp[2][2], dnu_mp[2], b2s[max_quads];
double    beta0[Cell_nLocMax][2], nu0[Cell_nLocMax][2];
double    Jx, Jy;
Vector2   nu;
FILE      *fp_lat;
ofstream  tune_out;


const char file_lat[] = "get_bn.lat";
const char home_dir[] = "/home/bengtsson";

const int  n_peaks = 10;


#if ORDER > 1


void get_map_normal_form()
{
  int       i;
  long int  lastpos;

  daeps_(1e-15);
  danot_(no_tps-1);
  getcod(0.0, lastpos);

  cout << endl;
  cout << "COD" << endl;
  for (i = 0; i < 2*nd_tps; i++)
    if (true) 
      cout << scientific << setprecision(5)
	   << setw(13) << globval.CODvect[i];
    else
      cout << scientific << setprecision(16)
	   << setw(24) << globval.CODvect[i];
  cout << endl;

  map.identity(); map += globval.CODvect;
  Cell_Pass(0, globval.Cell_nLoc, map, lastpos);
  map -= globval.CODvect;

  danot_(no_tps);
  MNF = MapNorm(map, no_tps);
}

#endif


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


void get_sigma()
{
  int              i, j, k, n, n_iter;
  long             lastpos;
  double           dx_abs = 0.0;
  iVector          jj;
  ss_vect<double>  x0, x1, dx, sigma[6];
  ss_vect<tps>     I, dx0, map;
  
  const double  i_max =    10;
  const double  eps   = 1e-10;

  n = (globval.Cavity_on)? 6 : 4;

  for (j = 0; j < 2*nd_tps; j++)
    jj[j] = (j < n)? 1 : 0;

  x0.zero();

  danot_(2);
  n_iter = 0; I.identity();
  for (i = 0; i < i_max; i++) {
    map.identity(); map += x0; Cell_Pass(0, globval.Cell_nLoc, map, lastpos); 

    x1 = map.cst(); dx = x0 - x1; dx0 = PInv(map-I-x1, jj)*dx;

    for (j = 0; j <= 5; j++)
      for (k = 0; k <= 5; k++)
	sigma[j][k] = map[j][k];

    dx_abs = xabs(n, dx); x0 += dx0.cst();

    printf("%3d, err=%10.3e (%10.3e)\n", i, dx_abs, eps);
    printf(" %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	   x0[x_], x0[px_], x0[y_], x0[py_], x0[delta_], x0[ct_]);
    if ((dx_abs < eps) || (lastpos != globval.Cell_nLoc)) break;
  }
}


double get_psv(const ss_vect<tps> A, const ss_vect<double> ps, double twoJ[])
{
  /* Note:

        T                      -1
       M  J M = J,    M = A R A  ,    A, R symplectic

        -1       T           |  0  I |        T   | beta   -alpha |
       A   = -J A  J,    J = |       |,    A A  = |               |
                             | -I  0 |            | -alpha  gamma |

          T             T        T  T
       A A  -> M A (M A)  = M A A  M

     Transform to Floquet Space:

        -1         T
       A   x = -J A  J x,

                               -1    T    -1              T    T
       2I_1 + 2I_2 + 2I_3 = ( A   x )  ( A   x ) = ( J x )  A A  ( J x )

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
  fstream          outf;

  globval.Cavity_on = true;

  Ring_GetTwiss(true, 0.0); printglob();

  putlinmat(6, globval.Ascr, A);

  outf.open("lin_inv.out", ios::out);
  ps.zero(); ps[x_] = Ax; ps[y_] = Ay; ps[delta_] = delta;
  for (i = 1; i <= n; i++) {
    Cell_Pass(0, globval.Cell_nLoc, ps, lastpos); V = get_psv(A, ps, twoJ);

    outf << setprecision(3) << scientific
	 << setw(3) << i << setw(10) << 1e6*twoJ[X_]
	 << setw(10) << 1e6*twoJ[Y_] << setw(10) << 1e3*1e2*twoJ[Z_]
	 << setw(10) << 1e6*V << endl;
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
  fstream          outf;

  const double  eps = 1e-10;

  globval.Cavity_on = true;

  outf.open("sympl_form.out", ios::out);
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

    outf << setprecision(3) << scientific
	 << setw(3) << i << setw(11) << omega[X_]
	 << setw(11) << omega[Y_] << setw(11) << omega[Z_] << setw(11) << V
	 << endl;
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
  Matrix           R;
  ss_vect<double>  ps, ps_Fl;
  fstream          outf;

  const int     n_turn  = 15*1;
  const double  A_min = 0.1e-3, twoJ_max = 1e-5, phi_max = 0.5*2.0*M_PI;
  const double  phi_Fl[] = { 0.0*M_PI/180.0, 0.0*M_PI/180.0 };

  outf.open("fixed_points.out", ios::out);

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

//      ZeroMat(4, R);
//      R[x_][x_]  = cos(phi_Fl[X_]); R[x_][px_]  = sin(phi_Fl[X_]);
//      R[px_][x_] = sin(phi_Fl[X_]); R[px_][px_] = cos(phi_Fl[X_]);
//      R[y_][y_]  = cos(phi_Fl[Y_]); R[y_][py_]  = sin(phi_Fl[Y_]);
//      R[py_][y_] = sin(phi_Fl[Y_]); R[py_][py_] = cos(phi_Fl[Y_]);

//      LinTrans(4, R, ps_Fl);

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

      outf << fixed << setprecision(3)
	   << setw(3) << n1
	   << setw(8) << 1e3*A[X_] << setw(8) << 1e3*A[Y_]
	   << scientific
	   << setw(10) << 1e6*s_twoJ[X_] << setw(10) << 1e6*s_twoJ[Y_]
	   << fixed << setprecision(5)
	   << setw(10) << s_phi[X_] << setw(10) << s_phi[Y_]
	   << scientific << setprecision(3)
	   << setw(10) << log(1.0+sqr(s_twoJ[X_])+sqr(s_twoJ[Y_]))
	   << setw(10) << log(1.0+sqr(s_phi[X_])+sqr(s_phi[Y_]))
	   << endl;
    }
    outf << endl;
  }
  outf.close();
}


bool track(const double x, const double px, const double y, const double py,
	   const double delta, const long int n, const double f_rf,
	   const bool prt)
{
  long int         i, lastpos;
  ss_vect<double>  ps;
  ofstream         os;

  ps[x_] = x; ps[px_] = px; ps[y_] = y; ps[py_] = py;
  ps[delta_] = delta; ps[ct_] = 0.0;

  if (prt) {
    os.open("track.out", ios::out);
    os << "# Tracking with Thor" << endl;
    os << "#" << endl;
    os << "#  n       x           p_x          y            p_y  "
       << "       delta         cdt" << endl;
    if (f_rf == 0.0) {
      os << "#         [mm]        [mrad]       [mm]         [mrad]"
	 << "                    [mm]" << endl;
      os << scientific << setprecision(16)
	 << setw(4) << 0
	 << setw(24) << 1e3*ps[x_] << setw(24) << 1e3*ps[px_]
	 << setw(24) << 1e3*ps[y_] << setw(24) << 1e3*ps[py_]
	 << setw(24) << 1e2*ps[delta_] 
	 << setw(24) << 1e3*ps[ct_] << endl;
    } else {
      os << "#         [mm]        [mrad]       [mm]         [mrad]"
	 << "                    [deg]" << endl;
      os << scientific << setprecision(16)
	 << setw(4) << 0
	 << setw(24) << 1e3*ps[x_] << setw(24) << 1e3*ps[px_]
	 << setw(24) << 1e3*ps[y_] << setw(24) << 1e3*ps[py_]
	 << setw(24) << 1e2*ps[delta_] 
	 << setw(24) << 2.0*f_rf*180.0*ps[ct_]/c0 << endl;
    }
    os << "#" << endl;
  }

  for (i = 1; i <= n; i++) {
    Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
    if (lastpos == globval.Cell_nLoc) {
      if (prt)
	if (f_rf == 0.0)
	  os << scientific << setprecision(16)
	     << setw(4) << i
	     << setw(24) << 1e3*ps[x_] << setw(24) << 1e3*ps[px_]
	     << setw(24) << 1e3*ps[y_] << setw(24) << 1e3*ps[py_]
	     << setw(24) << 1e2*ps[delta_] 
	     << setw(24) << 1e3*ps[ct_] << endl;
	else
	  os << scientific << setprecision(16)
	     << setw(4) << i
	     << setw(24) << 1e3*ps[x_] << setw(24) << 1e3*ps[px_]
	     << setw(24) << 1e3*ps[y_] << setw(24) << 1e3*ps[py_]
	     << setw(24) << 1e2*ps[delta_] 
	     << setw(24) << 2.0*f_rf*180.0*ps[ct_]/c0 << endl;
    } else
      return false;
  }
  if (prt) os.close();

  return true;
}


void get_r_stable(double &r, const double phi, const double delta,
		  const long int n, const double eps)
{
  /* Binary search for dynamical aperture. */
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


double get_dynap(const double r, const double delta, const int n,
		 const double eps, const int n_pts,
		 double x_min[], double x_max[])
{
  /* Determine the dynamical aperture by tracking.
     Assumes mid-plane symmetry.                    */

  int           i, j;
  double        r1, phi, x0[2], x1[2], x2[2], DA;
  ofstream      os;

  os.open("dynap.dat", ios::out);
  os << "# Dynamical Aperture:" << endl;
  os << "#    x      y" << endl;
  os << "#   [mm]   [mm]" << endl;
  os << "#" << endl;

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
    for (j = 0; j < 2; j++) {
      x_min[j] = min(x2[j], x_min[j]); x_max[j] = max(x2[j], x_max[j]);
    }
    if (i == 0) {
      x0[X_] = x2[X_]; x0[Y_] = x2[Y_];
    } else
      DA += x1[X_]*x2[Y_] - x2[X_]*x1[Y_];
    os << fixed << setprecision(2)
	 << setw(8) << 1e3*x2[X_] << setw(8) << 1e3*x2[Y_] << endl;
    x1[X_] = x2[X_]; x1[Y_] = x2[Y_];
  }
  DA += x2[X_]*x0[Y_] - x0[X_]*x2[Y_];
  // factor of 2 from mid-plane symmetry
  DA = fabs(DA)/sqrt(Cell[globval.Cell_nLoc].Beta[X_]
       *Cell[globval.Cell_nLoc].Beta[Y_]);

  cout << fixed << setprecision(1)
       << "DA^ = " << setw(6) << 1e6*DA
       << ", x^ = " << setw(5) << 1e3*x_min[X_] << " - "
       << setw(4) << 1e3*x_max[X_] << " mm"
       << ", y^ = " << setw(5) << 1e3*x_min[Y_] << " - "
       << setw(4) << 1e3*x_max[Y_] << " mm" << endl;

  return DA;
} 


void prt_scan(const int n, const double nu_x, const double nu_y,
	      const double DA, const double x_min[], const double x_max[])
{

  tune_out << fixed << setprecision(5)
	   << "n = " << setw(4) << n
	   << ", nu_x= " << nu_x << ", nu_y= " << nu_y
	   << setprecision(1)
	   << ", DA^= " << setw(6) << 1e6*DA
	   << ", x^= " << setw(5) << 1e3*x_min[X_] << " - "
	   << setw(4) << 1e3*x_max[X_] << " mm"
	   << ", y^= " << setw(5) << 1e3*x_min[Y_] << " - "
	   << setw(4) << 1e3*x_max[Y_] << " mm" << endl;
}


bool eng_tol(const bool config, const int n_orbit)
{
  bool  cod;
  int   i;

  if (config) {
    globval.bpm = ElemIndex("bpm");
    globval.hcorr = ElemIndex("ch"); globval.vcorr = ElemIndex("cv");

    // Clear trim setpoints
    for (i = 1; i <= GetnKid(globval.vcorr); i++) 
      SetKLpar(globval.vcorr, i, -Dip, 0.0);
    for (i = 1; i <= GetnKid(globval.hcorr); i++) 
      SetKLpar(globval.hcorr, i, Dip, 0.0);
    
    gcmat(globval.bpm, globval.hcorr, 1); gcmat(globval.bpm, globval.vcorr, 2);
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
//    Q_Fam[N_Fam++] = ElemIndex("sh1");
    Q_Fam[N_Fam++] = ElemIndex("qh1");
    Q_Fam[N_Fam++] = ElemIndex("qh2");
    Q_Fam[N_Fam++] = ElemIndex("qh3");

    Q_Fam[N_Fam++] = ElemIndex("ql1");
    Q_Fam[N_Fam++] = ElemIndex("ql2");
    Q_Fam[N_Fam++] = ElemIndex("ql3");
  }

  if (false) {
    Q_Fam[N_Fam++] = ElemIndex("sh1");
//    Q_Fam[N_Fam++] = ElemIndex("twk");
//    Q_Fam[N_Fam++] = -ElemIndex("qh2_dw");
//    Q_Fam[N_Fam++] = ElemIndex("qh1");
    Q_Fam[N_Fam++] = ElemIndex("qh2");
    Q_Fam[N_Fam++] = ElemIndex("qh3");
    Q_Fam[N_Fam++] = ElemIndex("ql1");
    Q_Fam[N_Fam++] = ElemIndex("ql2");
    Q_Fam[N_Fam++] = ElemIndex("ql3");
  }

  printf("\n");
  printf("No of quadrupole families: %d\n", N_Fam);

  ini_ID_corr();

  set_IDs(1.0);

  if (false) orb_corr(5);

  ID_corr(7, 4);
}


void get_DA(const char *file_dir, const int n,
	    const double nu_x, const double nu_y, const bool prt)
{
  bool    cod;
  char    str[max_str];
  double  DA, DA_sum, x_min[2], x_max[2];

  const double  delta   = 2.5e-2;

  sprintf(str, "%s%s", file_dir, "/quad.dat"); get_bn(str, n, false);

//  no_sxt();

  Ring_GetTwiss(false, 0.0);

  if (IDs) {
    get_IDs(); set_IDs(0.0);
  }

  cod = eng_tol(true, n_orbit);
  sprintf(str, "%s%s", file_dir, "/sext.dat"); get_bn(str, n, false);

  globval.Cavity_on = true;

  cout << endl;
  cout << fixed << setprecision(5)
       << "#n = " << setw(4) << n
       << ", nu_x = " << nu_x << ", nu_y = " << nu_y << ", "; 
  if (cod)
    DA = get_dynap(10e-3, -delta, n_track, 0.1e-3, n_aper, x_min, x_max);
  else
    DA = 0.0;
  DA_sum = DA;
  if (prt) {
    tune_out << "#";
    prt_scan(n, nu_x, nu_y, DA, x_min, x_max);
  }

  cout << fixed << setprecision(5)
       << "#n = " << setw(4) << n
       << ", nu_x = " << nu_x << ", nu_y = " << nu_y << ", "; 
  if (cod)
    DA = get_dynap(10e-3, delta, n_track, 0.1e-3, n_aper, x_min, x_max); 
  else
    DA = 0.0;
  DA_sum += DA;
  if (prt) {
    tune_out << "#";
    prt_scan(n, nu_x, nu_y, DA, x_min, x_max);
  }
	
  cout << fixed << setprecision(5)
       << " n = " << setw(4) << n
       << ", nu_x = " << nu_x << ", nu_y = " << nu_y << ", "; 
  if (cod)
    DA = get_dynap(10e-3, 0.0, n_track, 0.1e-3, n_aper, x_min, x_max); 
  else
    DA = 0.0;
  DA_sum += DA; DA = DA_sum/3.0;
  if (prt) {
    tune_out << " ";
    prt_scan(n, nu_x, nu_y, DA, x_min, x_max);
  }
}


void tune_scan(const char *file_dir)
{
  char      line[max_str], str[max_str];
  int       n;
  double    nu_x, nu_y;
  ifstream  inf;

  sprintf(str, "%s%s", file_dir, "/tune_scan.dat");
  inf.open(str, ios::in);

  tune_out.open("tune_scan.dat", ios::out);

  cout << endl;
  while (inf.getline(line, max_str) != NULL) {
    if (strstr(line, "#") == NULL) {
      if (strcmp(line, "") != 0) {
	globval.Cavity_on = false;

	sscanf(line, "%*s %*s %d, nu_x= %lf, nu_y= %lf", &n, &nu_x, &nu_y);

	get_DA(file_dir, n, nu_x, nu_y, true);

      } else
	tune_out << endl;
    }
  }

  inf.close(); tune_out.close();
}


void get_curly_H(long int k1, long int k2)
{
  long int      k, lastpos;
  double        dnu[2];
  ss_vect<tps>  Id, ps, A, ps_Floq, M;
  ofstream      outf;

  Ring_GetTwiss(true, 0.0);

  file_wr(outf, "../out/curly_H.out");

  Id.identity(); ps.zero();
  ps[x_] = Cell[0].Eta[x_]*Id[delta_]; ps[px_] = Cell[0].Etap[px_]*Id[delta_];
  ps[delta_] = Id[delta_];
  A.identity(); putlinmat(4, globval.Ascr, A);

  if (false) {
    cout << "A.A^tp" << endl;
    prt_lin_map(2, A*tp_S(2, A));
    exit(0);
  }

  Cell_Pass(0, k1-1, ps, lastpos); Cell_Pass(0, k1-1, A, lastpos);

  A = get_A_CS(A, dnu);
  
  for (k = k1-1; k <= k2; k++) {
    Cell_Pass(k-1, k, ps, lastpos); Cell_Pass(k-1, k, A, lastpos);

    A = get_A_CS(A, dnu);

    ps_Floq.zero();
    ps_Floq[x_] = ps[x_][delta_]; ps_Floq[px_] = ps[px_][delta_];
    ps_Floq = Inv(A)*ps_Floq;

    outf << scientific << setprecision(3)
	 << setw(4) << k
	 << setw(11) << ps_Floq[x_].cst() << setw(11) << ps_Floq[px_].cst()
	 << setw(11) << (sqr(ps_Floq[x_])+sqr(ps_Floq[px_])).cst()
	 << endl;
  }

  outf.close();

  prt_lat("linlat.out", globval.bpm, true);
}


void prt_mom_map(const double m2[][n_m2])
{
  int  i, j;

  cout << endl;
  for (i = 0; i < n_m2; i++) {
    for (j = 0; j < n_m2; j++)
//      cout << scientific << setprecision(3)
//	   << setw(11) << m2[i][j];
      cout << scientific << setprecision(16)
	   << setw(24) << m2[i][j];
    cout << endl;
  }
}


void get_moments()
{
  int              i, j, k, l, i1, k1;
  long int         lastpos;
  double           twoJ[2], mm2[n_m2][n_m2];
  tps              m2[n_m2];
  tps              sigma;
  ss_vect<double>  ps;
  ss_vect<tps>     A, Id, ps1, dM1, sgma;
  ofstream         outf;

  const int  n_turn = 10;

  Id.identity();

  globval.radiation = false; globval.Cavity_on = false;

  getcod(0.0, lastpos);

  get_map();

  prt_lin_map(3, map);

  prt_lin_map(3, tp_S(3, map)*map);

  if (true) {
    cout << endl;
    i1 = 0;
    for (i = 0; i < 2*nd_tps; i++)
      for (j = i; j < 2*nd_tps; j++) {
        k1 = 0;
        for (k = 0; k < 2*nd_tps; k++) {
          for (l = k; l < 2*nd_tps; l++) {
            mm2[i1][k1] = map[i][k]*map[j][l];
            k1++;
          }
        }
        i1++;
      }

//    prt_mom_map(mm2);
  }

  for (i = 0; i < 6; i++)
    cout << scientific << setprecision(16)
       << setw(24)
       << mm2[i][0]+mm2[i][6]+mm2[i][11]+mm2[i][15]+mm2[i][18]+mm2[i][20]
       << setw(24);

  if (false) {
    ps1.identity(); ps1 += globval.CODvect;
    Cell_Pass(0, globval.Cell_nLoc, ps1, lastpos);
    get_m2(ps1, m2);

    for (i = 0; i < 2; i++)
      cout << scientific << setprecision(3)
	   << setw(11) << m2[i];
  }

  danot_(1);
  sgma.identity(); sgma += globval.CODvect;
  for (i = 1; i <= globval.Cell_nLoc; i++) {
    dM1.identity(); Cell_Pass(i-1, i, dM1, lastpos);
    sgma = tp_S(3, dM1)*sgma*dM1;
  }

  for (i = 0; i < 2; i++)
    cout << sgma[i];

  danot_(2);
  sigma = 0.0;
  for (i = 0; i < 6; i++)
    sigma += sqr(Id[i]);
  Cell_Pass(0, globval.Cell_nLoc, sigma, lastpos);
  cout << sigma;

  exit(0);

  file_wr(outf, "../out/moments.out");

//  globval.radiation  = true;
//  globval.Cavity_on = true;

  Ring_GetTwiss(true, 0.0);

  A.zero(); putlinmat(6, globval.Ascr, A);

  cout << endl;
  ps.zero(); ps[x_] = 1e-3; ps[y_] = 1e-3;
  for (i = 1; i <= n_turn; i++) {
    Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);

    get_twoJ(2, ps, A, twoJ);

    cout << scientific << setprecision(3)
	 << "i = " << setw(4) << i
	 << setw(11) << twoJ[X_] << setw(11) << twoJ[Y_] << endl;
  }

  outf.close();

  prt_lat("linlat.out", globval.bpm, true);

  exit(0);
}


tps SLS_mag(const double bres, const double tlin, const double sfac,
	    const double Isat, const double aexp, const double nexp)
{
  int           k;
  iVector       jj;
  tps           Iabs, I;
  ss_vect<tps>  b;

  Iabs = tps(0.0, 1);

//  b[0] = bres + Iabs*tlin/pow(1.0-sfac*(pow(Iabs/Isat, aexp)), nexp);

  cout << b[0];

  for (k = 0; k < nv_tps; k++)
    jj[k] = 0;

  jj[0] = 1;

  return (PInv(b, jj))[0];
}


void get_slb(void)
{
  long int  loc;
  ofstream  outf;

  outf.open("slb.dat", ios::out);

  loc = Elem_GetPos(ElemIndex("slb"), 1);
  outf << scientific << setprecision(5)
       << setw(13) << Cell[loc].Alpha[X_] << setw(13) << Cell[loc].Alpha[Y_]
       << "\t // alpha" << endl;
  outf << scientific << setprecision(5)
       << setw(13) << Cell[loc].Beta[X_] << setw(13) << Cell[loc].Beta[Y_]
       << "\t // beta" << endl;
  outf << scientific << setprecision(5)
       << setw(13) << 2.0*Cell[loc].Nu[X_] << setw(13) << 2.0*Cell[loc].Nu[Y_]
       << "\t // dnu" << endl;

  outf.close();
}


void get_A1(const double alpha_x, const double beta_x,
	    const double alpha_y, const double beta_y)
{
  ss_vect<tps>  Id;

  map.zero(); MNF.A0.zero(); MNF.A1.zero(); MNF.g = 0.0; Id.identity();

  MNF.A1[x_]  = sqrt(beta_x)*Id[x_];
  MNF.A1[px_] = -alpha_x/sqrt(beta_x)*Id[x_] + 1.0/sqrt(beta_x)*Id[px_];

  MNF.A1[y_]  = sqrt(beta_y)*Id[y_];
  MNF.A1[py_] = -alpha_y/sqrt(beta_y)*Id[y_] + 1.0/sqrt(beta_y)*Id[py_];
}


ss_vect<tps> get_A(ss_vect<tps> &map)
{
  // Compute the Twiss parameters from a periodic map with mid-plane symmetry.

  int           k;
  double        c, s;
  Vector2       alpha, beta, eta, etap, nu;
  iVector       jj;
  ss_vect<tps>  Id, ps;

  Id.identity();

  for (k = 0; k < 2*nd_tps; k++)
    jj[k] = (k < 4)? 1 : 0;

  for (k = 0; k < 2*nd_tps; k++)
    ps[k] = map[k][delta_];

  ps = PInv(Id-map, jj)*ps;

  for (k = 0; k <= 1; k++) {
    eta[k] = ps[2*k].cst(); etap[k] = ps[2*k+1].cst();

    beta[k] = 2.0*map[2*k][2*k+1]
      /sqrt(4.0-sqr(map[2*k][2*k]+map[2*k+1][2*k+1]));

    c = (map[2*k][2*k]+map[2*k+1][2*k+1])/2.0; s = map[2*k][2*k+1]/beta[k];
    nu[k] = atan2(s, c)/(2.0*M_PI);
    if (nu[k] < 0.0) nu[k] += 1.0;

    alpha[k] = (map[2*k][2*k]-map[2*k+1][2*k+1])/(2.0*s);
  }

  return get_A(alpha, beta, eta, etap);
}


void get_chrom(const long int i1, const long int i2)
{
  long  int     lastpos;
  int           k;
  Vector2       nu1, nu2;
  Matrix        M;
  ss_vect<tps>  map;

  globval.emittance = false; globval.radiation = false;

  map.identity(); map[delta_] = -globval.dPcommon/2.0;
  Cell_Pass(i1, i2, map, lastpos); getlinmat(6, map, M);
  // get tune from symplectic matrix
  GetNu(nu1, M);
  map.identity(); map[delta_] = globval.dPcommon/2.0;
  Cell_Pass(i1, i2, map, lastpos); getlinmat(6, map, M);
  // get tune from symplectic matrix
  GetNu(nu2, M);

  for (k = 0; k <= 1; k++)
    globval.Chrom[k] = (nu2[k]-nu1[k])/globval.dPcommon;
}


void get_twiss(const bool chrom, const bool emit)
{
  long  int     lastpos;
  int           k;
  Matrix        R;
  ss_vect<tps>  map, A;

  if (chrom) get_chrom(0, globval.Cell_nLoc);

  globval.emittance = emit; globval.radiation = emit;

  map.identity(); Cell_Pass(0, globval.Cell_nLoc, map, lastpos);
  getlinmat(6, map, globval.OneTurnMat);
  GetNu(globval.TotalTune, globval.OneTurnMat);

  GDiag(4, Cell[globval.Cell_nLoc].S, globval.Ascr, globval.Ascrinv, R,
        globval.OneTurnMat, globval.Omega, globval.Alphac);

  putlinmat(6, globval.Ascr, A); Cell_Pass(0, globval.Cell_nLoc, A, lastpos);

  if (emit)
    for (k = 0; k <= 2; k++)
      globval.eps[k] = -globval.D_rad[k]/(4.0*globval.alpha_rad[k]);
}


void get_ab(const Matrix &A, double alpha[], double beta[],
	    double eta[], double etap[])
{
  ss_vect<tps>  A_map;

  A_map.zero();
  putlinmat(6, A, A_map); get_ab(A_map, alpha, beta, eta, etap);
}


void get_dnu(const Matrix &A, double dnu[])
{
  ss_vect<tps>  A_map;

  A_map.zero();
  putlinmat(6, A, A_map); get_dnu(A_map, dnu);
}


void calc_twiss(const long int i0, const long int i1, const bool prt)
{
  long int  j;
  int       k;
  Vector2   dnu0, dnu;

  if (prt) printf("\n");
  Cell[i0].Nu[X_] = 0.0; Cell[i0].Nu[Y_] = 0.0;
  for (j = i0; j <= i1; j++) {
    get_ab(Cell[j].A, Cell[j].Alpha, Cell[j].Beta, Cell[j].Eta, Cell[j].Etap);

    if (j > i0) {
      get_dnu(Cell[j-1].A, dnu0); get_dnu(Cell[j].A, dnu);
      for (k = 0; k <= 1; k++) {
	dnu[k] -= dnu0[k];
	if ((dnu[k] < 0.0) && (Cell[j].Elem.PL > 0.0)) dnu[k] += 1.0;
	Cell[j].Nu[k] = Cell[j-1].Nu[k] + dnu[k];
      }

      Cell[j].S = Cell[j-1].S + Cell[j].Elem.PL;
    }

    if (prt)
      printf("%4ld %15s %6.2f"
	     " %7.3f %6.3f %6.3f %6.3f %6.3f"
	     " %7.3f %6.3f %6.3f %6.3f %6.3f\n",
	     j, Cell[j].Elem.PName, Cell[j].S,
	     Cell[j].Alpha[X_], Cell[j].Beta[X_], Cell[j].Nu[X_],
	     Cell[j].Eta[X_], Cell[j].Etap[X_],
	     Cell[j].Alpha[Y_], Cell[j].Beta[Y_], Cell[j].Nu[Y_],
	     Cell[j].Eta[Y_], Cell[j].Etap[Y_]);
  }

  for (k = 0; k <= 1; k++)
    globval.TotalTune[k] = Cell[globval.Cell_nLoc].Nu[k];
}


void prt_twiss(void)
{
  long int  j;
  int       k;
  Vector2   dnu0, dnu;

  Cell[0].Nu[X_] = 0.0; Cell[0].Nu[Y_] = 0.0;
  for (j = 0; j <= globval.Cell_nLoc; j++) {
    get_ab(Cell[j].A, Cell[j].Alpha, Cell[j].Beta, Cell[j].Eta, Cell[j].Etap);

    if (j > 0) {
      get_dnu(Cell[j-1].A, dnu0); get_dnu(Cell[j].A, dnu);
      for (k = 0; k <= 1; k++) {
	dnu[k] -= dnu0[k];
	if (dnu[k] < 0.0) dnu[k] += 1.0;
	Cell[j].Nu[k] = Cell[j-1].Nu[k] + dnu[k];
      }
    }

    printf("%4ld %15s %6.2f"
	   " %7.3f %6.3f %6.3f %6.3f %6.3f"
	   " %7.3f %6.3f %6.3f %6.3f %6.3f\n",
	   j, Cell[j].Elem.PName, Cell[j].S,
	   Cell[j].Alpha[X_], Cell[j].Beta[X_], Cell[j].Nu[X_],
	   Cell[j].Eta[X_], Cell[j].Etap[X_],
	   Cell[j].Alpha[Y_], Cell[j].Beta[Y_], Cell[j].Nu[Y_],
	   Cell[j].Eta[Y_], Cell[j].Etap[Y_]);
  }

  for (k = 0; k <= 1; k++)
    globval.TotalTune[k] = Cell[globval.Cell_nLoc].Nu[k];
}


float f_cell(float prms[])
{
  int     i;
  float   f;

  const int     n_prt       = 10;
  const double  eps_x       = 2.0e-9;
  const double  max_chrom[] = { -7.6, -3.1 };
  const double  beta_ss[]   = {  1.0,  1.0 };
  const double  beta_ls[]   = { 10.0,  3.0 };
  const double  max_eta_x   = 0.45;

  n_iter++;

  for (i = 1; i <= n_prm; i++)
    set_bn_design_fam(Fnum[i-1], Quad, prms[i], 0.0);

  get_twiss(true, true);

  for (i = 0; i < n_loc; i++)
    get_ab(Cell[loc[i]].A, Cell[loc[i]].Alpha, Cell[loc[i]].Beta,
	 Cell[loc[i]].Eta, Cell[loc[i]].Etap);

  if (globval.stable)
    f = 1e25*sqr(globval.eps[X_]-eps_x)
      + 1e10*sqr(globval.TotalTune[X_]-fract(nu[X_]))
      + 1e10*sqr(globval.TotalTune[Y_]-fract(nu[Y_]))
      + 1e3*sqr(globval.Chrom[X_]-max_chrom[X_])
      + 1e3*sqr(globval.Chrom[X_]-max_chrom[Y_])
      + 1e8*sqr(Cell[loc[0]].Alpha[X_]) + 1e8*sqr(Cell[loc[0]].Alpha[Y_])
      + 1e2*sqr(Cell[loc[0]].Eta[X_]-max_eta_x)
      + 1e2*sqr(Cell[loc[1]].Eta[X_]) + 1e2*sqr(Cell[loc[1]].Etap[X_])
      + 1e2*sqr(Cell[loc[2]].Beta[X_]-beta_ss[X_])
      + 1e2*sqr(Cell[loc[2]].Beta[Y_]-beta_ss[Y_])
      + 1e2*sqr(Cell[loc[3]].Beta[X_]-beta_ls[X_])
      + 1e2*sqr(Cell[loc[3]].Beta[Y_]-beta_ls[Y_]);
  else
    f = 1e30;

  if (n_iter % n_prt == 0) {
    printf("\n");
    printf("  eps_x = %9.3e nu = (%7.5f, %7.5f) ksi = (%6.3f, %6.3f)\n",
	   globval.eps[X_], globval.TotalTune[X_], globval.TotalTune[Y_],
	   globval.Chrom[X_], globval.Chrom[Y_]);
    printf("  alpha_ss = (%8.1e, %8.1e) beta_ss = (%6.3f, %5.3f)\n",
	   Cell[loc[2]].Alpha[X_], Cell[loc[2]].Alpha[Y_],
	   Cell[loc[2]].Beta[X_], Cell[loc[2]].Beta[Y_]);
    printf("  alpha_ls = (%8.1e, %8.1e) beta_ls = (%6.3f, %5.3f)\n",
	   Cell[loc[3]].Alpha[X_], Cell[loc[3]].Alpha[Y_],
	   Cell[loc[3]].Beta[X_], Cell[loc[3]].Beta[Y_]);
    printf("  alpha_mp = (%8.1e, %8.1e) beta_mp = (%6.3f, %5.3f)"
	   " eta_x^ = %5.3f\n",
	   Cell[loc[0]].Alpha[X_], Cell[loc[0]].Alpha[Y_],
	   Cell[loc[0]].Beta[X_], Cell[loc[0]].Beta[Y_],
	   Cell[loc[0]].Eta[X_]);
    printf("  eta_x = %8.1e eta\'_x = %8.1e\n",
	   Cell[loc[1]].Eta[X_], Cell[loc[1]].Etap[X_]);

    printf("\n");
    printf("%4d %11.3e", n_iter, f);
    for (i = 1; i <= n_prm; i++)
      printf(" %8.5f", prms[i]);
    printf("\n");
  }

  return f;
}


void fit_cell(const double nu_x, const double nu_y)
{
  int       iter, i, j;
  float     *prms, **xi, fret;
  double    bn, an;

  const float  ftol  = 1e-10;

  cout << endl;
  cout << "fit_cell" << endl;

  n_prm = 0;
  Fnum[n_prm++] = ElemIndex("qh1"); Fnum[n_prm++] = ElemIndex("qh2");
  Fnum[n_prm++] = ElemIndex("qh3");
  Fnum[n_prm++] = ElemIndex("ql1"); Fnum[n_prm++] = ElemIndex("ql2");
  Fnum[n_prm++] = ElemIndex("qh3");
  Fnum[n_prm++] = ElemIndex("qm1"); Fnum[n_prm++] = ElemIndex("qm2");

  prms = vector(1, n_prm); xi = matrix(1, n_prm, 1, n_prm);

  n_loc = 0;
  loc[n_loc++] = Elem_GetPos(ElemIndex("sm2"), 1);
  loc[n_loc++] = Elem_GetPos(ElemIndex("b1"), 2);
  loc[n_loc++] = Elem_GetPos(ElemIndex("ss"), 1);
  loc[n_loc++] = Elem_GetPos(ElemIndex("ls"), 1);

  nu[X_] = nu_x; nu[Y_] = nu_y;

  for (i = 1; i <= n_prm; i++) {
    get_bn_design_elem(Fnum[i-1], 1, Quad, bn, an); prms[i] = bn;
    for (j = 1; j <= n_prm; j++)
	xi[i][j] = (i == j)? 1e-3 : 0.0;
  }

  cout << endl;
  n_iter = 0;
  powell(prms, xi, n_prm, ftol, &iter, &fret, f_cell);
  n_iter = -1; f_cell(prms);

  free_vector(prms, 1, n_prm); free_matrix(xi, 1, n_prm, 1, n_prm);
}


double get_ii00m(const int i, const int m, const int N)
{
  long int  k;
  int       n;
  double    b_nL, h_ii00m;

  n = 2*i + m; h_ii00m = 0.0;
  for (k = 0; k <= globval.Cell_nLoc; k++)
    if (Cell[k].Elem.Pkind == Mpole) {
      if (true)
	// relative
	b_nL = Cell[k].Elem.M->PBpar[HOMmax+N];
      else
	// actual
	b_nL = Cell[k].Elem.M->PBsys[HOMmax+n];
      if (Cell[k].Elem.PL != 0.0) b_nL *= Cell[k].Elem.PL;

      h_ii00m += b_nL*pow(Cell[k].Beta[X_], i)*pow(Cell[k].Eta[X_], m);
    }

  h_ii00m *= -N_cell*nok(n, m)*nok(2*i, i)/(n*pow(2.0, 2*i));

  return h_ii00m;
}


void get_tol(const int N, const int n, const double twoJ_x)
{
  int     i, k, k0;
  double  h_ii00m, sum;

  k0 = (n % 2 != 0)? 1 : 2;
  sum = 0.0;
  for (k = k0; k <= n-2; k += 2) {
    i = (n-k)/2;
    h_ii00m = get_ii00m(i, k, N);
    sum += i*h_ii00m*pow(twoJ_x, i-1)*pow(delta_tol, k);
    printf("\n");
    printf("h_%d%d00%d    = %12.5e dnu_x = %12.5e (2J_x)^%d delta^%d\n",
	   i, i, k, h_ii00m, -i*h_ii00m/M_PI, i-1, k);
    printf("dnu_x      = %12.5e\n",
	   -i*h_ii00m/M_PI*pow(twoJ_x, i-1)*pow(delta_tol, k));
  }

  printf("\n");
  printf("dB_%d/B_%d   = %8.5f\n", n, N,
	 -1e4*pow(R0, n-N)*M_PI*dnu_x_tol/sum);
}


void get_h_sys(void)
{
  double  twoJ_x;

  twoJ_x = sqr(Ax_tol)/Cell[globval.Cell_nLoc].Beta[X_];

  printf("\n");
  printf("For full lattice (N = %d)\n", N_cell);
  printf("\n");
  printf("Quadrupoles:\n");
  get_tol(Quad, 6, twoJ_x);
  printf("\n");
  get_tol(Quad, 10, twoJ_x);
  printf("\n");
  get_tol(Quad, 14, twoJ_x);

  printf("\n");
  printf("Sextupoles:\n");
  get_tol(Sext, 9, twoJ_x);
  printf("\n");
  get_tol(Sext, 15, twoJ_x);
  printf("\n");
  get_tol(Sext, 21, twoJ_x);
}


void get_cod_rms(const double dx, const double dy,
		 const int n_seed, const bool all)
{
  const int  n_elem = 4000;

  bool     cod;
  int      i, j, k, n;
  double   x1[n_elem][6], x2[n_elem][6];
  double   x_mean[n_elem][6], x_sigma[n_elem][6];
  FILE     *fp;

  if (globval.Cell_nLoc+1 > n_elem) {
    printf("n_elem exceeded: %d (%ld)\n", n_elem, globval.Cell_nLoc+1);
    exit_(0);
  }

  for (j = 0; j <= globval.Cell_nLoc; j++)
    for (k = 0; k < 6; k++) {
      x1[j][k] = 0.0; x2[j][k] = 0.0;
    }
  
  fp = file_write("cod_rms.out");
  
  for (i = 0; i < n_seed; i++) {
    misalign_rms_fam(ElemIndex("dip"), dx, dy, 0.0, true);
    misalign_rms_fam(ElemIndex("qf"), dx, dy, 0.0, true);
    misalign_rms_fam(ElemIndex("qfend"), dx, dy, 0.0, true);
    misalign_rms_fam(ElemIndex("qdend"), dx, dy, 0.0, true);
    
    cod = orb_corr(3);

    if (cod) {
      n = 0;
      for (j = 0; j <= globval.Cell_nLoc; j++)
	if (all || ((Cell[j].Elem.Pkind == Mpole) &&
		    (Cell[j].Elem.M->n_design == Sext))) {
	  n++;
	  for (k = 0; k < 6; k++) {
	    x1[n-1][k] += Cell[j].BeamPos[k];
	    x2[n-1][k] += sqr(Cell[j].BeamPos[k]);
	  }
	}
    } else
      printf("orb_corr: failed\n");

    // reset orbit trims
    set_bn_design_fam(globval.hcorr, Dip, 0.0, 0.0);
    set_bn_design_fam(globval.vcorr, Dip, 0.0, 0.0);
  }
  
  n = 0;
  for (j = 0; j <= globval.Cell_nLoc; j++)
    if (all || ((Cell[j].Elem.Pkind == Mpole) &&
		(Cell[j].Elem.M->n_design == Sext))) {
      n++;
      for (k = 0; k < 6; k++) {
	x_mean[n-1][k] = x1[n-1][k]/n_seed;
	x_sigma[n-1][k] = sqrt((n_seed*x2[n-1][k]-sqr(x1[n-1][k]))
			       /(n_seed*(n_seed-1.0)));
      }
      fprintf(fp, "%8.3f %6.2f %10.3e +/- %10.3e %10.3e +/- %10.3e\n",
	      Cell[j].S, get_code(Cell[j]),
	      1e3*x_mean[n-1][x_], 1e3*x_sigma[n-1][x_],
	      1e3*x_mean[n-1][y_], 1e3*x_sigma[n-1][y_]);
    } else
      fprintf(fp, "%8.3f %6.2f\n", Cell[j].S, get_code(Cell[j]));
  
  fclose(fp);
}


void corr_stat(const bool hor)
{
  long    i, n, Fnum;
  double  b1L, a1L, sum, sum2, m, s;

  sum = 0.0; sum2 = 0.0;
  Fnum = (hor)? globval.hcorr : globval.vcorr; n = GetnKid(Fnum);

  for (i = 1; i <= n; i++) {
    get_bnL_design_elem(Fnum, i, Dip, b1L, a1L);
    sum += (hor)? b1L : a1L; sum2 += (hor)? sqr(b1L) : sqr(a1L);
  }

  m = (n != 0)? sum/n : 0.0;

  s = (n != 0 && n != 1)? s = (n*sum2-sqr(sum))/(n*(n-1.0)) : 0.0;
  s = (s >= 0.0)? sqrt(s) : 0.0;

  printf("\n");
  if (hor) {
    printf("D(b1L) = %8.1e+/-%7.1e\n", m, s);
  } else {
    printf("D(a1L) = %8.1e+/-%7.1e\n", m, s);
  }
}


void set_bn_s(const int Fnum, const int Knum, const double b2)
{
  long int  k = 0;
  double    a2 = 0.0;

  if (Fnum > 0)
    set_bn_design_elem(Fnum, Knum, Quad, b2, a2);
  else {
    // point to multipole
    k = Elem_GetPos(abs(Fnum), Knum);

    switch (Cell[k-1].Elem.PName[1]) {
    case 'u':
      if (Cell[k+1].Elem.PName[1] == 'd') {
	set_L(Cell[k-1].Fnum, Cell[k-1].Knum, b2);
	set_L(Cell[k+1].Fnum, Cell[k+1].Knum, -b2);
      } else {
	cout << "set_bn_s: configuration error 1 " << Cell[k+1].Elem.PName
	     << " (" << k+1 << ")" << endl;
	exit(0);
      }
      break;
    case 'd':
      if (Cell[k+1].Elem.PName[1] == 'u') {
	set_L(Cell[k-1].Fnum, Cell[k-1].Knum, -b2);
	set_L(Cell[k+1].Fnum, Cell[k+1].Knum, b2);
      } else {
	cout << "set_bn_s: configuration error 2 " << Cell[k+1].Elem.PName
	     << " (" << k+1 << ")" << endl;
	exit(0);
      }
      break;
    default:
      cout << "set_bn_s: configuration error 3 " << Cell[k-1].Elem.PName
	   << " (" << k << ")" << endl;
      exit(0);
      break;
    }
  }
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
	      1e2*(Cell[i].Beta[X_]-beta0[i][X_])/beta0[i][X_],
	      Cell[i].Nu[X_]-nu0[i][X_],
	      1e2*(Cell[i].Beta[Y_]-beta0[i][Y_])/beta0[i][Y_],
	      Cell[i].Nu[Y_]-nu0[i][Y_]);
    }
  }

  fclose(outf);
}


float f_ID_match(float prms[])
{
  long int  lastpos;
  int       i;
  float     f;
  double    dnu[2];

  const int  n_prt = 10;

  n_iter++;

  for (i = 1; i <= n_prm; i++) {
    set_bn_s(sgn(quads1[i-1])*Cell[abs(quads1[i-1])].Fnum,
	     Cell[abs(quads1[i-1])].Knum, prms[i]);
    set_bn_s(sgn(quads2[i-1])*Cell[abs(quads2[i-1])].Fnum,
	     Cell[abs(quads2[i-1])].Knum, prms[i]);
  }

//  Ring_GetTwiss(false, 0.0);
  get_twiss(false, false);
//  get_A1(alpha_mp[0][X_], beta_mp[0][X_], alpha_mp[0][Y_], beta_mp[0][Y_]);
//  Cell_Pass(loc[0], loc[1], MNF.A1, lastpos);
  calc_twiss(loc[0], loc[1], false);

  for (i = 0; i <= 1; i++)
    dnu[i] = fract(Cell[loc[1]].Nu[i]-Cell[loc[0]].Nu[i]);

  if (globval.stable)
    f = 1e6*sqr(dnu[X_]-dnu_mp[X_])
      + 1e6*sqr(dnu[Y_]-dnu_mp[Y_])
      + 1e3*sqr(Cell[loc[0]].Alpha[X_]-alpha_mp[0][X_])
      + 1e3*sqr(Cell[loc[0]].Alpha[Y_]-alpha_mp[0][Y_])
      + 1e0*sqr(Cell[loc[0]].Beta[X_]-beta_mp[0][X_])
      + 1e0*sqr(Cell[loc[0]].Beta[Y_]-beta_mp[0][Y_])
      + 1e3*sqr(Cell[loc[1]].Alpha[X_]-alpha_mp[1][X_])
      + 1e3*sqr(Cell[loc[1]].Alpha[Y_]-alpha_mp[1][Y_])
      + 1e0*sqr(Cell[loc[1]].Beta[X_]-beta_mp[1][X_])
      + 1e0*sqr(Cell[loc[1]].Beta[Y_]-beta_mp[1][Y_]);
  else
    f = 1e30;

  if (n_iter % n_prt == 0) {
    printf("\n");
    printf("  nu    = (%6.3f, %6.3f)\n",
	   dnu[X_]-dnu_mp[X_], dnu[Y_]-dnu_mp[Y_]);
    printf("  alpha = (%6.3f, %6.3f)\n",
	   Cell[loc[0]].Alpha[X_]-alpha_mp[0][X_],
	   Cell[loc[0]].Alpha[Y_]-alpha_mp[0][Y_]);
    printf("  beta  = (%6.3f, %6.3f)\n",
	   Cell[loc[0]].Beta[X_]-beta_mp[0][X_],
	   Cell[loc[0]].Beta[Y_]-beta_mp[0][Y_]);
    printf("  alpha = (%6.3f, %6.3f)\n",
	   Cell[loc[1]].Alpha[X_]-alpha_mp[1][X_],
	   Cell[loc[1]].Alpha[Y_]-alpha_mp[1][Y_]);
    printf("  beta  = (%6.3f, %6.3f)\n",
	   Cell[loc[1]].Beta[X_]-beta_mp[1][X_],
	   Cell[loc[1]].Beta[Y_]-beta_mp[1][Y_]);

    printf("\n");
    printf("%4d %11.3e", n_iter, f);
    for (i = 1; i <= n_prm; i++)
      printf(" %8.5f", prms[i]);
    printf("\n");
  }

  return f;
}


void ID_match(void)
{
  int       iter, i, j;
  float     *prms, **xi, fret;
  double    b2, a2, b2s_init[max_quads];

  const float  ftol  = 1e-3;

  cout << endl;
  cout << "ID_match" << endl;

  n_prm = 0;

//  quads1[n_prm]   = Elem_GetPos(ElemIndex("sh1"), 2);
//  quads2[n_prm++] = Elem_GetPos(ElemIndex("sh1"), 3);
//  quads1[n_prm]   = Elem_GetPos(ElemIndex("qh1"), 2);
//  quads2[n_prm++] = Elem_GetPos(ElemIndex("qh1"), 3);
  quads1[n_prm]   = Elem_GetPos(ElemIndex("qh2"), 2);
  quads2[n_prm++] = Elem_GetPos(ElemIndex("qh2"), 3);
  quads1[n_prm]   = Elem_GetPos(ElemIndex("qh3"), 2);
  quads2[n_prm++] = Elem_GetPos(ElemIndex("qh3"), 3);

  if (true) {
//    quads1[n_prm]   = -Elem_GetPos(ElemIndex("qh1"), 2);
//    quads2[n_prm++] = -Elem_GetPos(ElemIndex("qh1"), 3);
    quads1[n_prm]   = -Elem_GetPos(ElemIndex("qh2"), 2);
    quads2[n_prm++] = -Elem_GetPos(ElemIndex("qh2"), 3);
    quads1[n_prm]   = -Elem_GetPos(ElemIndex("qh3"), 2);
    quads2[n_prm++] = -Elem_GetPos(ElemIndex("qh3"), 3);
  }

  prms = vector(1, n_prm); xi = matrix(1, n_prm, 1, n_prm);

  n_loc = 0;
  loc[n_loc++] = Elem_GetPos(ElemIndex("mp"), 2);
  loc[n_loc++] = Elem_GetPos(ElemIndex("mp"), 3);

  Ring_GetTwiss(false, 0.0);

  if (false) {
    set_IDs(1.0);

    n_prm = 0;

    prms[++n_prm] =  0.23028;
    prms[++n_prm] = -0.59955;
    prms[++n_prm] =  1.50151;
    prms[++n_prm] = -1.80667;
    prms[++n_prm] = -0.32596;
    prms[++n_prm] = -0.42250;

/*    prms[++n_prm] = -0.87503;
    prms[++n_prm] =  1.81526;
    prms[++n_prm] = -1.84911;
    prms[++n_prm] = -0.91296;
    prms[++n_prm] = -0.31610;
    prms[++n_prm] = -0.37166;*/

    for (i = 1; i <= n_prm; i++) {
      set_bn_s(sgn(quads1[i-1])*Cell[abs(quads1[i-1])].Fnum,
	       Cell[abs(quads1[i-1])].Knum, prms[i]);
      set_bn_s(sgn(quads2[i-1])*Cell[abs(quads2[i-1])].Fnum,
	       Cell[abs(quads2[i-1])].Knum, prms[i]);
    }
  }

  for (i = 0; i <= 1; i++) {
    dnu_mp[i] = fract(Cell[loc[1]].Nu[i]-Cell[loc[0]].Nu[i]);
    alpha_mp[0][i] = Cell[loc[0]].Alpha[i];
    beta_mp[0][i] = Cell[loc[0]].Beta[i];
    alpha_mp[1][i] = Cell[loc[1]].Alpha[i];
    beta_mp[1][i] = Cell[loc[1]].Beta[i];
  }

  for (i = 1; i <= n_prm; i++) {
    if (quads1[i-1] > 0) {
      get_bn_design_elem(Cell[quads1[i-1]].Fnum, Cell[quads1[i-1]].Knum,
			 Quad, b2, a2);
      prms[i] = b2;
    } else
      // point to drift adjacent to quad
      prms[i] = get_L(Cell[abs(quads1[i-1])-1].Fnum,
		      Cell[abs(quads1[i-1])-1].Knum);

    b2s_init[i-1] = prms[i];
    for (j = 1; j <= n_prm; j++)
      xi[i][j] = (i == j)? 1e-3 : 0.0;
  }

//  n_prm -= 2;
  set_IDs(1.0);

  cout << endl;
  n_iter = 0;
  powell(prms, xi, n_prm, ftol, &iter, &fret, f_ID_match);
  n_iter = -1; f_ID_match(prms);

  printf("\n");
  printf("grad changes\n");
  printf("\n");
  for (i = 1; i <= n_prm; i++)
    if (quads1[i-1] > 0)
      printf("%10s %6.3f (db2)\n", Cell[abs(quads1[i-1])].Elem.PName,
	     prms[i]-b2s_init[i-1]);
    else
      printf("%10s %6.3f (ds)\n", Cell[abs(quads1[i-1])].Elem.PName,
	     prms[i]-b2s_init[i-1]);

  free_vector(prms, 1, n_prm); free_matrix(xi, 1, n_prm, 1, n_prm);
}


int main(int argc, char *argv[])
{
  char          ch, str1[max_str], str2[max_str];
  int           i, j, n, Fnum, qf, qd, sf, sd, qse, qsf, qsg, N;
  long int      loc, lastpos, lastn;
  double        nu_x, nu_y, sigma_delta, sigma_s, dnu_x, alpha[2], beta[2];
  double        db2L, dx_rms, dbr_n, d, eta[2], etap[2];
  Matrix        M;
  tps           H, H_re, H_im, h, h_re, h_im, K_re, K_im, g_re, g_im;
  ss_vect<tps>  ps, S, Id, nus;
  ifstream      inf;
  ofstream      outf;
  FILE          *fp;


  const long  seed = 1121;

  iniranf(seed); setrancut(1.0);

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  freq_map = false; tune_shift = true;

//  n_x = 100; n_y = 100; n_dp = 100; n_tr = 2064;
//  x_max_FMA = 40e-3; y_max_FMA = 40e-3; delta_FMA = 3e-2;
//  x_max_FMA = 20e-3; y_max_FMA = 15e-3; delta_FMA = 3e-2;
  n_x = 50; n_y = 30; n_dp = 25; n_tr = 2064;
//  x_max_FMA = 40e-3; y_max_FMA = 20e-3; delta_FMA = 3e-2;
//  x_max_FMA = 40e-3; y_max_FMA = 7e-3; delta_FMA = 3e-2;

//  n_x = 20; n_y = 15; n_dp = 15; n_tr = 2064;
//  x_max_FMA = 25e-3; y_max_FMA = 10e-3; delta_FMA = 2.5e-2;
  x_max_FMA = 20e-3; y_max_FMA = 6e-3; delta_FMA = 3e-2;
//  x_max_FMA = 20e-3; y_max_FMA = 3e-3; delta_FMA = 3e-2;



  // disable from TPSALib- and LieLib log messages
//  idprset(-1);

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile("flat_file.dat");

  globval.EPU = false;

  IDs = GetnKid(ElemIndex("dw_90")) > 0;
  if (!IDs) IDs = GetnKid(ElemIndex("ivu_20_g59")) > 0;
  if (!IDs) IDs = GetnKid(ElemIndex("u55")) > 0;

  IDs = false;

  printf("\n");
  printf("IDs = %d\n", IDs);

  if (false) {
    globval.bpm = ElemIndex("bpm");
    globval.hcorr = ElemIndex("ch"); globval.vcorr = ElemIndex("cv");
  }

  if (IDs) {
    get_IDs(); set_IDs(0.0);
  }

//  get_twiss(true, true); prt_twiss(); exit(0);

  if (false) {
    no_sxt();
    fit_cell(2.224, 1.085);
  }

//  no_sxt();

  Ring_GetTwiss(false, 0.0);

  if (false) {
    qf = ElemIndex("qh2"); qd = ElemIndex("qh3");
    FitTune(qf, qd, globval.TotalTune[X_]+0.01, globval.TotalTune[Y_]-0.015);
    printf("\n");
    printf("nu_x = %7.5f nu_y = %7.5f\n",
	   globval.TotalTune[X_], globval.TotalTune[Y_]);
    printf("  qf = %8.5f   qd = %8.5f\n",
	   GetKpar(qf, 1, Quad), GetKpar(qd, 1, Quad));
  }

  if (false) {
//    sf = ElemIndex("sf2h"); sd = ElemIndex("sd1");
//    sf = ElemIndex("sf2"); sd = ElemIndex("sf1");
    sf = ElemIndex("sm2h"); sd = ElemIndex("sm1");
    FitChrom(sf, sd, 0.0, 0.0);
    printf("\n");
    printf("sf = %8.5f sd = %8.5f\n",
	   GetKpar(sf, 1, Sext), GetKpar(sd, 1, Sext));
  }

  // qa
//  cout << SLS_mag(0.074, 0.1968, 1.0, 136.3, 10.6, 0.077);

//  exit(0);

  ch = 'n';
//  cout << "exact Hamiltonian (y/n)? "; scanf("%c", &ch);
  globval.H_exact = (ch == 'y')? true : false;
  if (globval.H_exact) bend_cal();

//  get_curly_H(30, 71);

//  get_moments();

  if (false) {
    Fnum = ElemIndex("epu1");
    for (n = 1; n <= GetnKid(Fnum); n++) {
      loc = Elem_GetPos(Fnum, n);
      Cell[loc].Elem.M->PdTpar = 90.0; Mpole_SetdT(Fnum, n);
    }
  }

  if (false) {
    cout << endl;
    ps.identity(); Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);

    prt_lin_map(3, ps);

    S = get_S(3);

    prt_lin_map(3, ps*S*tp_S(3, ps));

    cout << endl;
    cout << scientific << setprecision(5)
	 << "trace(M_y) = " << setw(13) << ps[y_][y_]+ps[py_][py_] << endl;

    prt_lat("linlat.out", globval.bpm, true);

//    prtmfile("flat_file.dat");

//    exit(0);
  }

  get_map(); prt_lin_map(3, map);

  Ring_GetTwiss(true, 0.0); printglob();

  if (wp != 0) {
//    cout << "n?> "; scanf("%d", &n);

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
    default:
      cout << "undef. working point" << endl;
      exit(1);
      break;
    }

    set_tune(str1, str2, n);

    Ring_GetTwiss(true, 0.0); printglob();

    if (true) {
      N = 0;
      Fnum = ElemIndex("b1"); N += GetnKid(Fnum);
//      Fnum = ElemIndex("b1d"); N += GetnKid(Fnum);
      N /= 4;

      printf("\n");
      printf("N = %d\n", N);

      switch (wp) {
      case 3:
	// need to avoid nu_x + 4 nu_y = 96
	nu_x = N*33.47/15.0; nu_y = N*15.68/15.0;
	break;
      case 5:
//	nu_x = N*33.15/15.0; nu_y = N*15.69/15.0;
	// need to avoid nu_x + 4 nu_y = 96
	nu_x = N*33.15/15.0; nu_y = N*15.68/15.0;
	break;
      case 6:
//	nu_x = N*33.15/15.0; nu_y = N*16.3/15.0;
	// need to avoid 3 nu_x + 2 nu_y = 132
	nu_x = N*33.14/15.0; nu_y = N*16.27/15.0;
	break;
      case 7:
//	nu_x = N*33.45/15.0; nu_y = N*16.35/15.0;
	// test of working points
//	nu_x = N*33.47/15.0; nu_y = N*16.33/15.0;
//	nu_x = N*33.43/15.0; nu_y = N*16.335/15.0;
//	nu_x = N*33.46/15.0; nu_y = N*16.39/15.0;
	// new working point
	nu_x = N*33.45/15.0; nu_y = N*16.37/15.0;
//	nu_x = N*33.46/15.0; nu_y = N*16.36/15.0;
//	nu_x = N*33.46/15.0; nu_y = N*16.35/15.0;
	break;
      default:
	cout << "undef. working point" << endl;
	exit(1);
	break;
      }
	
      qf = ElemIndex("qh2"); qd = ElemIndex("qh3");

      FitTune(qf, qd, nu_x, nu_y);

      printf("\n");
      printf("nu_x = %7.5f nu_y = %7.5f\n",
	     globval.TotalTune[X_], globval.TotalTune[Y_]);
      printf("  qf = %8.5f   qd = %8.5f\n",
	     GetKpar(qf, 1, Quad), GetKpar(qd, 1, Quad));

      Ring_GetTwiss(true, 0.0); printglob();
    }
  }

  if (get_thor) {
    no_sxt();

    sprintf(str1, "%s%s", home_dir, "/Thor-2.0/thor/wrk");
    sprintf(str2, "%s%s", home_dir, "/Thor-2.0/thor/wrk");

    switch (get_thor) {
    case 1:
      switch (wp) {
      case 3:
	if (!IDs)
	  get_bn(strcat(str2, "/sext_3.bare_3"), 0, true);
//          get_bn(strcat(str2, "/sext_3.chrom_2"), 0, true);
//          get_bn(strcat(str2, "/sext_3.chrom_2_geom_6"), 0, true);
//	    get_bn(strcat(str2, "/sext_3.chrom_2_geom_7"), 0, true);
//	    get_bn(strcat(str2, "/sext_3.ds_chrom_2_geom_7"), 0, true);
	else
	  get_bn(strcat(str2, "/sext_3.DW_3"), 0, true);
	break;
      case 5:
	if (!IDs)
	  get_bn(strcat(str2, "/sext_3.bare_5"), 0, true);
	else
	  get_bn(strcat(str2, "/sext_3.DW_5"), 0, true);
	break;
      case 6:
	if (!IDs)
	  get_bn(strcat(str2, "/sext_3.bare_6"), 0, true);
//	    get_bn(strcat(str2, "/sext_3.dat"), 0, true);
	else
	  get_bn(strcat(str2, "/sext_3.DW_6"), 0, true);
//	    get_bn(strcat(str2, "/sext_3.DW_6_ksi_1.0_4.0"), 0, true);
//	    get_bn(strcat(str2, "/sext_3.DW_6_ksi_1.0_4.0_dec"), 0, true);
	break;
      case 7:
	if (!IDs)
//	    get_bn(strcat(str2, "/sext_3.bare_7"), 0, true);
	  get_bn(strcat(str2, "/sext_3.bare_7_2"), 0, true);
	else 
//	    get_bn(strcat(str2, "/sext_3.DW_7"), 0, true);
	  get_bn(strcat(str2, "/sext_3.DW_7.2"), 0, true);
	break;
      default:
	cout << "undef. working point" << endl;
	exit(1);
	break;
      }
      break;
    case 2:
//      get_bn(strcat(str1, "/fit_tune.dat"), 0, true);
//      get_bn(strcat(str2, "/fit_chrom.dat"), 0, true);
//      get_bn(strcat(str2, "/sext_2.dat"), 0, true);
//      get_bn(strcat(str2, "/sext_3.dat"), 0, true);
	get_bn("sext_3.dat", 0, true);
      break;
    otherwise:
      break;
    }

    Ring_GetTwiss(true, 0.0); printglob();
  }

  for (i = 0; i <= globval.Cell_nLoc; i++)
    for (j = 0; j <= 1; j++) {
      beta0[i][j] = Cell[i].Beta[j]; nu0[i][j] = Cell[i].Nu[j];
    }

  prt_lat("linlat.out", globval.bpm, true);

  if (false) {
    d = 2.5e-2;

    sprintf(str1, "cod_dp%+4.1f.out", 1e2*d);
    getcod(d, lastpos); prt_cod(str1, globval.bpm, true);
    sprintf(str1, "cod_dp%+4.1f.out", -1e2*d);
    getcod(-d, lastpos); prt_cod(str1, globval.bpm, true);

    sprintf(str1, "linlat_dp%+4.1f.out", 1e2*d);
    Ring_GetTwiss(true, d); prt_lat(str1, globval.bpm, true);
    sprintf(str1, "linlat_dp%+4.1f.out", -1e2*d);
    Ring_GetTwiss(true, -d); prt_lat(str1, globval.bpm, true);
  }

  prt_chrom_lat();

  if (false) {
    if (false)
      track("track.out", 23.2e-3, 0.0, 0e-3, 0.0,
	    0e-2, 512, lastn, lastpos, 0, 0.0);
    else
      if (true)
	track("dJ.out",
	      sqr(0.6e-3)/Cell[globval.Cell_nLoc].Beta[X_], 0.0,
	      sqr(0.2e-3)/Cell[globval.Cell_nLoc].Beta[Y_], 0.0,
	      0e-2, 20*512, lastn, lastpos, 2, 0.0);
      else if (false)
	// 3nu_x, nu_x + 2nu_y
	track("dJ.out",
	      sqr(16.7e-3)/Cell[globval.Cell_nLoc].Beta[X_], 0.0,
	      sqr(10e-3)/Cell[globval.Cell_nLoc].Beta[Y_], 0.0,
	      0e-2, 20*512, lastn, lastpos, 2, 0.0);
      else if (false)
	// 4nu_x
	track("dJ.out",
	      sqr(23e-3)/Cell[globval.Cell_nLoc].Beta[X_], 0.0,
	      sqr(0e-3)/Cell[globval.Cell_nLoc].Beta[Y_], 0.0,
	      0e-2, 20*512, lastn, lastpos, 2, 0.0);
      else if (false)
	// 3nu_x + 2nu_y
	track("dJ.out",
	      sqr(-15e-3)/Cell[globval.Cell_nLoc].Beta[X_], 0.0,
	      sqr(10e-3)/Cell[globval.Cell_nLoc].Beta[Y_], 0.0,
	      0e-2, 20*512, lastn, lastpos, 2, 0.0);

//    track(0e-3, 0e0, 0.0, 0.0, -3e-2, 512, 0.0, true);

//    getcod(0.0, lastpos);
//    prt_cod("cod.out", globval.bpm, true);

//    exit(0);
  }

  if (IDs) {
    if (false) {
      set_IDs(1.0);

      map.identity(); Fnum = ElemIndex("dw_90");
      Cell_Pass(Elem_GetPos(Fnum, 1), Elem_GetPos(Fnum, GetnKid(Fnum)),
		map, lastpos);
      prt_lin_map(3, map);

      exit(0);
    }

    if (true)
      ID_correct();
    else
      ID_match();

    if (false) {
//      FitTune(qf, qd, nu_x, nu_y);
//      printf("\n");
//      printf("nu_x = %7.5f nu_y = %7.5f\n",
//	         globval.TotalTune[X_], globval.TotalTune[Y_]);
//      printf("  qf = %8.5f   qd = %8.5f\n",
//	         GetKpar(qf, 1, Quad), GetKpar(qd, 1, Quad));

      sf = ElemIndex("sm1"); sd = ElemIndex("sm2h");
      FitChrom(sf, sd, 0.0, 0.0);
      printf("\n");
      printf("sf = %8.5f sd = %8.5f\n",
	     GetKpar(sf, 1, Sext), GetKpar(sd, 1, Sext));
    }

    Ring_GetTwiss(true, 0.0); printglob();

    prt_lat("linlat.out", globval.bpm, true);
    prt_dlat("dlinlat.out", true);
    prt_cod("cod.out", globval.bpm, true);

//    LoadApers("/home/bengtsson/projects/in/Aper/Apertures_DW.dat", 1.0, 1.0);

//    globval.Aperture_on = true;
  }

//   Ring_GetTwiss(true, 0.0); printglob();

  GetEmittance(ElemIndex("cav"), true);

  sigma_s = sqrt(Cell[0].sigma[ct_][ct_]); 
  sigma_delta = sqrt(Cell[0].sigma[delta_][delta_]); 

  if (false)
    LoadFieldErr("/home/bengtsson/projects/in/Field_Err/jb.dat",
		 false, 1.0, true);

//  no_sxt();

//  get_h_sys();

  prtmfile("flat_file.dat");

  //    exit(0);

  if (tune_shift) {
    cout << endl;
    cout << "computing tune shifts" << endl;
    dnu_dA(20e-3, 10e-3, 0.0);
    get_ksi2(3.0e-2);

    //      get_alphac(); get_alphac2();

    //      exit(0);
  }

  if (tune_scn) {
    sprintf(str1, "%s%s", home_dir, "/R3L3SBS4-C6x6Qfgdr/2");
    tune_scan(str1);

    exit(0);
  }

  //    no_sxt();

  Ring_GetTwiss(true, 0.0); printglob();

  //    get_slb();

  //    rdmfile("flat_file.dat");

  if (false) {
    gcmat(globval.bpm, globval.hcorr, 1);
    gcmat(globval.bpm, globval.vcorr, 2);

    get_cod_rms(100e-6, 100e-6, 1000, true);
  }

  if (false) {
    no_sxt(); eng_tol(true, n_orbit); get_bn(str2, n, false);

    Ring_GetTwiss(true, 0.0); printglob();
  }

  GetEmittance(ElemIndex("cav"), true);

  sigma_s = sqrt(Cell[0].sigma[ct_][ct_]); 
  sigma_delta = sqrt(Cell[0].sigma[delta_][delta_]); 

  if (freq_map) {
    fmap(n_x, n_y, n_tr, x_max_FMA, y_max_FMA, 0.0, true, true);
    fmapdp(n_x, n_dp, n_tr, x_max_FMA, -delta_FMA, 1e-3, true, true);
  }

  if (false) {
    //aperture:
    //      ini_aper(-50e-3, 50e-3, -25e-3, 25e-3); 
      
    //      if (false)
    //	set_aper(ElemIndex("d0id"), -25e-3, 25e-3, -2.5e-3, 2.5e-3);
    //	set_aper(ElemIndex("d_u55"), -25e-3, 25e-3, -19e-3/2.0, 19e-3/2.0);
    //      else
    //	set_aper(ElemIndex("d0id"), -25e-3, 25e-3, -5e-3, 5e-3);
      
    //      PrintCh();
      
    //      Touschek(1.3e-9, globval.delta_RF, globval.eps[X_], 0.008e-9,
    //	       sigma_delta, sigma_s);

    //      Touschek(1.3e-9, globval.delta_RF, globval.eps[X_], 1e-11,
    //	       sigma_delta, sigma_s, 512, true);

    globval.delta_RF = 3.0e-2;

//    IBS(1.3e-9, 0.85e-9, 0.008e-9, 0.84e-3, 2*4.4e-3);

    Touschek(1.3e-9, globval.delta_RF, 0.85e-9, 0.008e-9, 0.84e-3, 4.4e-3);
      
    double  sum_delta[globval.Cell_nLoc+1][2];
    double  sum2_delta[globval.Cell_nLoc+1][2];
    double  mean_delta_s[globval.Cell_nLoc+1][2];
    double  sigma_delta_s[globval.Cell_nLoc+1][2];

    // initialize momentum aperture arrays
    for(j = 0; j <= globval.Cell_nLoc; j++){
      sum_delta[j][0] = 0.0; sum_delta[j][1] = 0.0;
      sum2_delta[j][0] = 0.0; sum2_delta[j][1] = 0.0;
    }
 
    Touschek(1.3e-9, globval.delta_RF, false,
	     0.85e-9, 0.008e-9, 0.84e-3, 4.4e-3,
	     512, true, sum_delta, sum2_delta);

    fp = file_write("mom_aper.out"); 
    for(j = 0; j <= globval.Cell_nLoc; j++)
      fprintf(fp, "%4d %7.2f %5.3f %6.3f\n",
	      j, Cell[j].S, 1e2*sum_delta[j][0], 1e2*sum_delta[j][1]);

    exit(0);
  }

  if (false) {
    printf("db2L? "); scanf("%lf", &db2L);
    set_bnr_rms_type(Quad, Quad, db2L, 0.0, true);
  }

  if (false) {
    printf("dx_rms? "); scanf("%lf", &dx_rms);
    misalign_rms_type(Sext, dx_rms, dx_rms, 0.0, true);
  }

//  Ring_GetTwiss(true, 0.0);
//  get_fixed_points(60e-3, 50e-3, 25);
//  get_fixed_points(15e-3, 10e-3, 25);
//  get_lin_inv(5e-3, 5e-3, 0e-2, 512);
//  get_sympl_form(10e-3, 6e-3, 0e-2, 1024);

  n_aper = 25;

  globval.Cavity_on = true;
//  get_dynap(delta);

//  Ring_GetTwiss(true, 0.0); printglob();
//  get_nu(9.9e-3, 25.2e-3, nu_x, nu_y);
//  cout << endl;
//  cout << fixed << setprecision(5)
//       << "nu_x = " << nu_x << ", nu_y = " << nu_y << endl;

#if ORDER == 1
//  get_map(); prt_lin_map(3, map);

#else
//  get_map_normal_form();

//  get_sigma();
#endif
}
