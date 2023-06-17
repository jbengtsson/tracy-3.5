#define NO 5

#include "tracy_lib.h"

int
  no_tps   = NO,
  ndpt_tps = 5;


const double
#if 1
  beta_inj[] = {3.1, 3.0},
#else
  beta_inj[] = {3.1, 2.8},
#endif
  A_max[]    = {3e-3, 1.5e-3},
  delta_max  = 2e-2,
  twoJ[]     = {sqr(A_max[X_])/beta_inj[X_], sqr(A_max[Y_])/beta_inj[Y_]},
  dnu[]      = {0.03, 0.02};


ss_vect<tps> get_map_Fl(const ss_vect<tps> &map)
{
  Matrix       R_mat;
  ss_vect<tps> Id, A0_inv, map1, R;

  const int n_dim = 4, no_delta = 1;

  Id.identity();

  // Find fixed point.
  GoFix(map, MNF.A0, A0_inv, no_delta);

  // Translate to fix point.
  map1 = A0_inv*map*MNF.A0;

  getlinmat(n_dim, map1, globval.OneTurnMat);
  GDiag(n_dim, 1.0, globval.Ascr, globval.Ascrinv, R_mat,
        globval.OneTurnMat, globval.Omega, globval.Alphac);

  MNF.A1 = putlinmat(n_dim, globval.Ascr);
  MNF.A1[ct_] = Id[ct_];
  MNF.A1[delta_] = Id[delta_];

  R = putlinmat(n_dim, R_mat);
  R[ct_] = Id[ct_];
  R[delta_] = Id[delta_];

  // Transform to Floquet space.
  return Inv(MNF.A1)*map1*MNF.A1;
}


ss_vect<tps> get_map_Fl_2(const ss_vect<tps> &map)
{
  const int no_delta = 1;

  MNF = MapNorm(map, no_delta);
  return Inv(MNF.A0*MNF.A1)*map*MNF.A0*MNF.A1;
}


tps get_mns(const tps &a, const int no1, const int no2)
{
  tps b;

  danot_(no1-1);
  b = -a;
  danot_(no2);
  b += a;
  danot_(no_tps);
  return b;
}


ss_vect<tps> get_mns(const ss_vect<tps> &x, const int no1, const int no2)
{
  int          k;
  ss_vect<tps> y;

  for (k = 0; k < nv_tps; k++)
    y[k] = get_mns(x[k], no1, no2);

  return y;
}


void Dragt_Finn_Fact(const ss_vect<tps> &map, ss_vect<tps> &M, tps &h)
{
  /* Dragt-Finn factorization:

       M = exp(:h_no:)*...*exp(:h_4:)*exp(:h_3:)*M                            */

  int          k;
  tps          hn;
  ss_vect<tps> Id, A0_inv, Fn, map1;

  Id.identity();

  danot_(1);
  M = map;

  danot_(no_tps);
  map1 = map*Inv(M);
  h = 0e0;
  for (k = 3; k <= no_tps; k++) {
    Fn = get_mns(map1, k-1, k-1);
    hn = Intd(Fn, -1.0);
    h += hn;
    map1 = map1*LieExp(-hn, Id);
  }
}


void Dragt_Finn_Fact_2(const ss_vect<tps> &map, ss_vect<tps> &M, tps &h)
{
  h = LieFact_DF(map, M);
}


ss_vect<tps> Dragt_Finn_Map(const ss_vect<tps> &M, const tps &h)
{
  int           k;
  ss_vect<tps>  Id, map;

  Id.identity();
  map = M;
  for (k = 3; k <= no_tps; k++)
    map = LieExp(Take(h, k), Id)*map;
  return map;
}


ss_vect<tps> Dragt_Finn_Map_2(const ss_vect<tps> &M, const tps &h)
{
  ss_vect<tps> Id;

  Id.identity();
  return FExpo(h, Id, 3, no_tps, 1)*M;
}


void chk_get_map_Fl(const ss_vect<tps> &map_Fl)
{
  ss_vect<tps> map_Fl_2;

  map_Fl_2 = get_map_Fl_2(map);
  daeps_(1e-6);
  cout << scientific << setprecision(3) << "\nmap_Fl-map_Fl_2\n"
       << setw(11) << map_Fl-map_Fl_2 << "\n";
  daeps_(eps_tps);
}


void chk_Dragt_Finn_Fact(const ss_vect<tps> &map_Fl, const ss_vect<tps> &R,
			 const tps &h)
{
  tps          h_2;
  ss_vect<tps> R_2;

  Dragt_Finn_Fact_2(map_Fl, R_2, h_2);
  daeps_(1e-6);
  cout << scientific << setprecision(3) << "\nR-R_2\n"
       << setw(11) << R-R_2 << "\n";
  cout << scientific << setprecision(3) << "\nh-h_2\n"
       << setw(11) << h-h_2 << "\n";
  daeps_(eps_tps);
}


void chk_Dragt_Finn_Map(const ss_vect<tps> &R, const tps &h,
			const ss_vect<tps> &map_Fl)
{
  ss_vect<tps> map_Fl_2;

  map_Fl_2 = Dragt_Finn_Map_2(R, h);
  danot_(no_tps-1);
  daeps_(1e-6);
  cout << scientific << setprecision(3) << "\nmap_Fl-map_Fl_2\n"
       << setw(11) << map_Fl-map_Fl_2 << "\n";
  daeps_(eps_tps);
}

tps get_k_2(const ss_vect<tps> &R)
{
  int          k;
  double       mu[2];
  tps          k_2, k_2_re, k_2_im;
  ss_vect<tps> Id;

  Id.identity();

  k_2_re = 0e0;
  for (k = 0; k < 2; k++) {
    mu[k] = atan2(R[2*k][2*k+1], R[2*k][2*k]);
    if (mu[k] < 0e0) mu[k] += 2e0*M_PI;
    k_2_re -= mu[k]*Id[2*k]*Id[2*k+1]/2e0;
  }
  printf("\nget_k_2: nu = [%7.5f, %7.5f]\n",
	 mu[X_]/(2e0*M_PI), mu[Y_]/(2e0*M_PI));
  k_2_im = 0e0;
  k_2 = RtoC(k_2_re, k_2_im);
  return k_2;
}


tps get_h(const tps &k_2, const tps &h)
{
  int           k;
  tps           K;
  ss_vect<tps>  Id;

  Id.identity();
  K = k_2;
  for (k = 3; k <= no_tps; k++)
    K = K*LieExp(Take(h, k), Id);
  return K;
}


tps get_h(const ss_vect<tps> &map) { return LieFact(map); }


void analyse(const ss_vect<tps> &map)
{
  tps          h_DF, h, k_2;
  ss_vect<tps> Id, map_Fl, map_Fl_2, R_DF;

  Id.identity();

  danot_(no_tps);

  printf("\nM:\n");
  prt_lin_map(3, map);

  map_Fl = get_map_Fl(map);
  printf("\nmap_Fl:\n");
  prt_lin_map(3, map_Fl);
  if (false) chk_get_map_Fl(map_Fl);

  Dragt_Finn_Fact(map_Fl, R_DF, h_DF);
  if (false) chk_Dragt_Finn_Fact(map_Fl, R_DF, h_DF);

  map_Fl_2 = Dragt_Finn_Map(R_DF, h_DF);
  printf("\nmap_Fl:\n");
  prt_lin_map(3, map_Fl_2);
  if (false) chk_Dragt_Finn_Map(R_DF, h_DF, map_Fl_2);

  k_2 = get_k_2(R_DF);
  cout << scientific << setprecision(3) << "\nk_2:\n"
       << setw(11) << k_2 << "\n";

  h = get_h(map_Fl*Inv(R_DF));
  daeps_(1e-10);
  cout << scientific << setprecision(3) << "\nh-h_DF:\n"
       << setw(11) << h-h_DF << "\n";

  h = h + k_2;
  cout << scientific << setprecision(3) << "\nh:\n"
       << setw(11) << h << "\n";

  if (false) {
    // Doesn't converge for no >= 4;
    // h = get_h(map_Fl*Inv(R));
    h = get_h(map_Fl);
    daeps_(1e-10);
    cout << scientific << setprecision(3) << "\nh:\n"
	 << setw(11) << h << "\n";
    cout << scientific << setprecision(3) << "\nM*h-h:\n"
	 << setw(11) << h*map_Fl-h << "\n";
  }
}


ss_vect<tps> get_A_nl(const tps g)
{
  ss_vect<tps>  Id;

  Id.identity();

  return FExpo(g, Id, 3, no_tps, -1);
}


ss_vect<tps> get_A_nl_inv(const tps g)
{
  ss_vect<tps>  Id;

  Id.identity();

  return FExpo(-g, Id, 3, no_tps, 1);
}


ss_vect<tps> get_R(ss_vect<tps> &map)
{
  Matrix       R_mat;
  ss_vect<tps> Id, A0_inv, R;

  const int n_dim = 4;

  Id.identity();

  // find fixed point
  GoFix(map, MNF.A0, A0_inv, no_tps);

  // translate to fix point
  map = A0_inv*map*MNF.A0;

  getlinmat(n_dim, map, globval.OneTurnMat);
  GDiag(n_dim, 1.0, globval.Ascr, globval.Ascrinv, R_mat,
        globval.OneTurnMat, globval.Omega, globval.Alphac);
  MNF.A1 = putlinmat(n_dim, globval.Ascr);
  MNF.A1[ct_] = Id[ct_];
  MNF.A1[delta_] = Id[delta_];
  R = putlinmat(n_dim, R_mat);
  R[ct_] = Id[ct_];
  R[delta_] = Id[delta_];

  // transform to Floquet space
  map = Inv(MNF.A1)*map*MNF.A1;

  return R;
}


tps get_Ker(const tps &h)
{
  long int      jj[ss_dim];
  int           i, j, k;
  tps           h_Ke;
  ss_vect<tps>  Id;

  for (k = 0; k < nv_tps; k++)
    jj[k] = 0;

  Id.identity(); h_Ke = 0.0;
  for (i = 0; i <= no_tps; i++) {
    jj[x_] = i; jj[px_] = i;
    for (j = 0; j <= no_tps; j++) {
      jj[y_] = j; jj[py_] = j;
      for (k = 0; k <= no_tps; k++) {
	jj[delta_] = k;
	if ((2*i+2*j+k <= no_tps) && ((i != 0) || (j != 0) || (k != 0))) {
	  h_Ke +=
	    h[jj]*pow(Id[x_], i)*pow(Id[px_], i)
	    *pow(Id[y_], j)*pow(Id[py_], j)*pow(Id[delta_], k);
	}
      }
    }
  }

  return h_Ke;
}


tps get_g(const tps nu_x, const tps nu_y, const tps &h)
{
  // Compute g = (1-R)^-1 * h.

  long int     jj1[ss_dim], jj2[ss_dim];
  int          i, j, k, l, m;
  double       re, im;
  tps          h_re, h_im, g_re, g_im, mn1, mn2, cotan;
  ss_vect<tps> Id;

  CtoR(h, h_re, h_im);

  for (k = 0; k < nv_tps; k++) {
    jj1[k] = 0; jj2[k] = 0;
  }

  Id.identity();
  g_re = 0e0;
  g_im = 0e0;
  for (i = 0; i <= no_tps; i++) {
    jj1[x_] = i; jj2[px_] = i;
    for (j = 0; j <= no_tps; j++) {
      jj1[px_] = j; jj2[x_] = j;
      for (k = 0; k <= no_tps; k++) {
	jj1[y_] = k; jj2[py_] = k;
	for (l = 0; l <= no_tps; l++) {
	  jj1[py_] = l; jj2[y_] = l;
	  if ((i+j+k+l <= no_tps) && ((i-j != 0) || (k-l != 0))) {
	    cotan = 1.0/tan(((i-j)*nu_x+(k-l)*nu_y)*M_PI);
	    mn1 =
	      pow(Id[x_], i)*pow(Id[px_], j)*pow(Id[y_], k)*pow(Id[py_], l);
	    mn2 =
	      pow(Id[x_], j)*pow(Id[px_], i)*pow(Id[y_], l)*pow(Id[py_], k);

	    for (m = 0; m <= no_tps-i-j-k-l; m++) {
	      jj1[delta_] = m;
	      jj2[delta_] = m;
	      re = h_re[jj1];
	      im = h_im[jj1];
	      // Compute g.
	      g_re += (re-cotan*im)*(mn1+mn2)*pow(Id[delta_], m)/2e0;
	      g_im += (im+cotan*re)*(mn1-mn2)*pow(Id[delta_], m)/2e0;
	      h_re.pook(jj2, 0e0);
	      h_im.pook(jj2, 0e0);
	    }
	  }
	}
      }
    }
  }

  return RtoC(g_re, g_im);
}


ss_vect<tps> map_norm(const ss_vect<tps> &map)
{
  int          k;
  double       nu0[2];
  tps          hn, hn_re, hn_im, h_ke, gn, Kn, g, K;
  ss_vect<tps> Id, R, A, nus, map1, map2;

  const bool prt = false;

  Id.identity();

  danot_(no_tps-1);

  map1 = map;
  R = get_R(map1);

  danot_(no_tps);

  K = 0e0;
  for (k = 0; k < 2; k++) {
    nu0[k] = atan2(R[2*k][2*k+1], R[2*k][2*k]);
    if (nu0[k] < 0.0) nu0[k] += 2e0*M_PI;
    nu0[k] /= 2e0*M_PI;
    K -= M_PI*nu0[k]*(sqr(Id[2*k])+sqr(Id[2*k+1]));
  }
  if (prt)
    cout << fixed << setprecision(5)
	 << "\nnu0 = (" << nu0[X_] << ", " << nu0[Y_] << ")" << endl;

  // coasting beam
  K += h_ijklm(map1[ct_], 0, 0, 0, 0, 1)*sqr(Id[delta_])/2.0;
  CtoR(K, hn_re, hn_im);
  if (prt)
    cout << scientific << setprecision(3)
	 << "\nK:\n" << hn_re;

  g = 0e0;
  for (k = 3; k <= no_tps; k++) {
    map2 = map1*Inv(R*FExpo(K, Id, 3, k-1, -1));
    hn = Intd(get_mns(map2, k-1, k-1), -1.0);
    gn = get_g(nu0[X_], nu0[Y_], hn);
    g += gn;
    CtoR(hn, hn_re, hn_im);
    Kn = RtoC(get_Ker(hn_re), get_Ker(hn_im)); K += Kn;
    if (prt)
      cout << scientific << setprecision(3)
	   << "\nk = " << k << ":\n" << hn_re;

    A = FExpo(gn, Id, k, k, -1);
    map1 = Inv(A)*map1*A;
  }

  if (prt) cout << scientific << setprecision(3) << "\nhn_re:\n" << hn_re;

  MNF = MapNorm(map, no_tps);

  CtoR(MNF.K-K, hn_re, hn_im);

  if (prt)
    cout << scientific << setprecision(3)
	 << "\nMNF.g-g:\n" << MNF.g-g
	 << "\nMNF.K-K:\n" << MNF.K-K;

  cout << scientific << setprecision(3)
       << "\ng:\n" << g << "\nK:\n" << K;

  return map1;
}


tps get_H(const tps &g, const tps &K)
{
  int          i;
  tps          H, H_d[4];
  ss_vect<tps> Id;

  // Construct generator.
  // K is in Dragt-Finn form but the generators commute.
  Id.identity();
  H = K;
  for (i = no_tps; i >= 3; i--)
    H = H*LieExp(-Take(g, i), Id);
  return H;
}


void H_inv(const ss_vect<tps> &map)
{
  int          i, j, k;
  double       dx[2], x[2];
  tps          H;
  ss_vect<tps> Id, ps, A, A_inv;
  FILE         *outf;

  const int
    n         = 50;
  const double
    x_max[]   = {5e-3, 5e-3};
  const string
    file_name = "H_inv.out";

  danot_(no_tps);

  Id.identity();

  if (true)
    MNF = MapNorm(map, no_tps);
  else
    map_norm(map);

  cout << scientific << setprecision(3)
       << "\nK:\n" << MNF.K << "\ng:\n" << MNF.g;

  A = MNF.A0*MNF.A1*get_A_nl(MNF.g);
  A_inv = Inv(A);

  if (false) {
    danot_(no_tps-1);
    cout << scientific << setprecision(3)
	 << "\nA*e^K*A^-1 - M:\n" << A*LieExp(MNF.K, Id)*A_inv-map;
    danot_(no_tps);
  }

  H = get_H(MNF.g, MNF.K)*Inv(MNF.A1)*Inv(MNF.A0);

  if (false)
    cout << scientific << setprecision(3) << "\nH - M*H:\n" << H*map-H;

  for (k = 0; k < 2; k++)
    dx[k] = x_max[k]/n;

  outf = file_write(file_name.c_str());

  x[X_] = -x_max[X_];
  for (i = -n; i <= n; i++) {
    x[Y_] = -x_max[Y_];
    for (j = -n; j <= n; j++) {
      ps.zero();
      for (k = 0; k < 2; k++)
	ps[2*k] = x[k];
      fprintf(outf, "%8.3f %8.3f %10.3e\n",
	      1e3*x[X_], 1e3*x[Y_], log(-(H*ps).cst()));
      x[Y_] += dx[Y_];
    }
    fprintf(outf, "\n");
    x[X_] += dx[X_];
  }

  fclose(outf);
}


tps g_renorm(const double nu0_x, const double nu0_y,
	     const double nu1_x, const double nu1_y,
	     const tps &g)
{
  // Renormalize g: (1-R^-1)^-1 * h 

  long int     jj1[ss_dim], jj2[ss_dim];
  int          i, j, k, l, m;
  double       re, im, cotan0, cotan1, cotan0_sqr;
  tps          h_re, h_im, g_re, g_im, G_re, G_im, mn1, mn2;
  ss_vect<tps> Id;

  CtoR(g, g_re, g_im);

  for (k = 0; k < ss_dim; k++) {
    jj1[k] = 0; jj2[k] = 0;
  }

  Id.identity(); G_re = 0.0; G_im = 0.0;
  for (i = 0; i <= no_tps; i++) {
    jj1[x_] = i; jj2[px_] = i;
    for (j = 0; j <= i; j++) {
      jj1[px_] = j; jj2[x_] = j;
      for (k = 0; k <= no_tps; k++) {
	jj1[y_] = k; jj2[py_] = k;
	for (l = 0; l <= no_tps; l++) {
	  jj1[py_] = l; jj2[y_] = l;

	  if (i+j+k+l <= no_tps) {
	    cotan0 = 1.0/tan(((i-j)*nu0_x+(k-l)*nu0_y)*M_PI);
	    cotan0_sqr = sqr(cotan0);
	    cotan1 = 1.0/tan(((i-j)*nu1_x+(k-l)*nu1_y)*M_PI);
	    mn1 =
	      pow(Id[x_], i)*pow(Id[px_], j)*pow(Id[y_], k)*pow(Id[py_], l);
	    mn2 =
	      pow(Id[x_], j)*pow(Id[px_], i)*pow(Id[y_], l)*pow(Id[py_], k);

	    for (m = 0; m <= no_tps; m++) {
	      if (i+j+k+l+m <= no_tps) {
		jj1[delta_] = m; jj2[delta_] = m;
		if ((i != j) || (k != l)) {
		  re = g_re[jj1]; im = g_im[jj1];

		  // compute h
		  h_re = (re+cotan0*im)*2.0/(1.0+cotan0_sqr);
		  h_im = (im-cotan0*re)*2.0/(1.0+cotan0_sqr);

		  // renormalize g
		  G_re += (h_re-cotan1*h_im)*(mn1+mn2)*pow(Id[delta_], m)/2.0;
		  G_im += (h_im+cotan1*h_re)*(mn1-mn2)*pow(Id[delta_], m)/2.0;
		  g_re.pook(jj2, 0.0); g_im.pook(jj2, 0.0);
		}
	      }
	    }
	  }
	}
      }
    }
  }

  return RtoC(G_re, G_im);
}


void track_H(const string &file_name, const double Ax, const double Ay, tps &H)
{
  long int        lastpos;
  int             k;
  double          h, J[2], phi[2], nu[2], scl, J2[2], phi2[2];
  tps             scl_tps;
  ss_vect<double> ps, ps_Fl, ps_nl;
  ss_vect<tps>    Id;
  FILE            *outf;

  const double x_ampl[] = {Ax, Ay};

  danot_(no_tps);

  Id.identity();

  ps.zero();
  ps[x_] = x_ampl[X_];
  ps[y_] = x_ampl[Y_];

  ps_Fl = (MNF.A1_inv*MNF.A0_inv*ps).cst();
  ps_nl = (MNF.A_nl_inv*ps_Fl).cst();

  // scale action-angle variables
  for (k = 0; k < 2; k++) {
    nu[k] = (MNF.nus[k]*ps_nl).cst();
    scl = sqrt(MNF.nus[k].cst()/nu[k]);
    ps_Fl[2*k] *= scl; ps_Fl[2*k+1] *= scl;
  }

  outf = file_write(file_name.c_str());

  for (k = 1; k <= 500; k++) {
    Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
    ps_Fl = (MNF.A1_inv*ps).cst();
    J[X_] = (sqr(ps_Fl[x_])+sqr(ps_Fl[px_]))/2.0;
    phi[X_] = atan2(ps_Fl[px_], ps_Fl[x_]);
    J[Y_] = (sqr(ps_Fl[y_])+sqr(ps_Fl[py_]))/2.0;
    phi[Y_] = atan2(ps_Fl[py_], ps_Fl[y_]);

    h = (H*ps_Fl).cst();

    ps_nl = (MNF.A_nl_inv*ps_Fl).cst();
    J2[X_] = (sqr(ps_nl[x_])+sqr(ps_nl[px_]))/2.0;
    phi2[X_] = atan2(ps_nl[px_], ps_nl[x_]);
    J2[Y_] = (sqr(ps_nl[y_])+sqr(ps_nl[py_]))/2.0;
    phi2[Y_] = atan2(ps_nl[py_], ps_nl[y_]);

    fprintf(outf, "%3d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e"
	    " %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	    k, 1e3*ps_Fl[x_], 1e3*ps_Fl[px_], 1e3*ps_Fl[y_], 1e3*ps_Fl[py_],
	    1e6*J[X_], phi[X_], 1e6*J[Y_], phi[Y_], 1e6*J2[X_], phi2[X_],
	    1e6*J2[Y_], phi2[Y_], 1e6*fabs(h));
  }

  fclose(outf);
}


void stability(const ss_vect<tps> &map)
{
  tps          H;
  ss_vect<tps> Id;

  MNF = MapNorm(map, no_tps);

  MNF.nus = dHdJ(MNF.K);

  MNF.A0_inv = Inv(MNF.A0);
  MNF.A1_inv = Inv(MNF.A1);

  // if (false)
  //   g_r = g_renorm(nus[0].cst(), nus[1].cst(), nu[X_], nu[Y_], MNF.g);
  // else
  //   g_r = MNF.g;

  H = get_H(MNF.g, MNF.K);
  MNF.A_nl = get_A_nl(MNF.g);
  MNF.A_nl_inv = get_A_nl_inv(MNF.g);

  // BESSY-III/b3_sf_40Grad_JB.lat.
  track_H("track_H_1.dat", 0.1e-3, 0.1e-3, H);
  track_H("track_H_2.dat", 0.5e-3, 0.5e-3, H);
  track_H("track_H_3.dat", 1.0e-3, 1.0e-3, H);
  track_H("track_H_4.dat", 1.0e-3, 1.0e-3, H);
  track_H("track_H_5.dat", 1.25e-3, 1.0e-3, H);
  track_H("track_H_6.dat", 1.4700475e-3, 1.0e-3, H);
  track_H("track_H_7.dat", 1.5e-3, 1.5e-3, H);
}


void get_sympl_form(FILE *outf, const double Ax, const double Ay,
		    const double delta)
{
  int             j, k;
  long int        lastpos;
  double          omega[3], domega[3], omega2[3], omega_rms[3], V_sum, V_mean;
  ss_vect<double> ps1, ps2, ps3, dq, dp;

  const int
    n = 10;
  const double
    eps = 1e-10,
    cut = 2e-20;

  globval.Cavity_on = false;

  ps1.zero();
  ps1[x_]     = Ax;
  ps1[y_]     = Ay;
  ps1[delta_] = delta;
  ps2 = ps1;
  ps3 = ps1;
  for (k = 0; k < 6; k++)
    if (k % 2 == 0)
      ps2[k] += eps;
    else
      ps3[k] += eps;
  for (k = 0; k < 3; k++)
    omega2[k] = 0e0;
  for (j = 1; j <= n; j++) {
    Cell_Pass(0, globval.Cell_nLoc, ps1, lastpos);
    Cell_Pass(0, globval.Cell_nLoc, ps2, lastpos);
    Cell_Pass(0, globval.Cell_nLoc, ps3, lastpos);
    dq = ps2 - ps1;
    dp = ps3 - ps1;
    V_sum = 0e0;
    for (k = 0; k < 3; k++) {
      omega[k] = dq[2*k]*dp[2*k+1] - dp[2*k]*dq[2*k+1];
      V_sum += omega[k];
      domega[k] = fabs(omega[k]) - sqr(eps);
      if (fabs(domega[k]) > cut) domega[k] = NAN;
      omega2[k] += sqr(domega[k]);
    }
  }

  for (k = 0; k < 3; k++)
    omega_rms[k] = sqrt(omega2[k])/sqr(n);
  V_mean = V_sum/n;

  fprintf(outf, "  %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	  Ax, Ay, delta, omega_rms[X_], omega_rms[Y_], omega_rms[Z_], V_mean);
}


void track_sympl_form(const int n, const double Ax, const double Ay)
{
  int             j, k;
  double          x[2], dx[2];
  FILE            *outf;

  const string
    file_name = "sympl_form.out";
  const double
    x_max[] = {Ax, Ay};

  outf = file_write(file_name.c_str());

  for (k = 0; k < 2; k++)
    dx[k] = x_max[k]/n;

  x[X_] = -x_max[X_];
  for (j = -n; j <= n; j++) {
    x[Y_] = -x_max[Y_];
    for (k = -n; k <= n; k++) {
      get_sympl_form(outf, x[X_], x[Y_], 0e0);
      x[Y_] += dx[Y_];
    }
    fprintf(outf, "\n");
    x[X_] += dx[X_];
  }

  fclose(outf);
}


void track_sympl_form_delta(const int n, const double Ax, const double Ay,
			    const double delta)
{
  int             j, k;
  double          x[2], dx[2];
  FILE            *outf;

  const string
    file_name = "sympl_form_delta.out";
  const double
    x_max[] = {Ax, delta};

  outf = file_write(file_name.c_str());

  for (k = 0; k < 2; k++)
    dx[k] = x_max[k]/n;

  x[X_] = -x_max[X_];
  for (j = -n; j <= n; j++) {
    x[1] = -x_max[1];
    for (k = -n; k <= n; k++) {
      get_sympl_form(outf, x[X_], Ay, x[1]);
      x[1] += dx[1];
    }
    fprintf(outf, "\n");
    x[X_] += dx[X_];
  }

  fclose(outf);
}


int main(int argc, char *argv[])
{
  int          k;
  ss_vect<tps> Id_scl;

  globval.H_exact    = false; globval.quad_fringe    = false;
  globval.Cavity_on  = false; globval.radiation      = false;
  globval.emittance  = false; globval.IBS            = false;
  globval.pathlength = false; globval.bpm            = 0;
  globval.Cart_Bend  = false; globval.dip_edge_fudge = true;
  globval.mat_meth   = false;

  // disable from TPSALib and LieLib log messages
  idprset(-1);

  std::string home_dir = "";

  daeps_(1e-25);

  Id_scl.identity();
  for (k = 0; k < 4; k++)
    Id_scl[k] *= sqrt(twoJ[k/2]);
  Id_scl[delta_] *= delta_max;

  if (!true)
    Read_Lattice((home_dir+argv[1]).c_str());
  else
    rdmfile(argv[1]);

  danot_(no_tps-1);
  get_map(false);
  danot_(no_tps);

  if (!false) analyse(map);

  if (false) H_inv(map);

  if (false) {
    track_sympl_form(30, 5e-3, 5e-3);
    track_sympl_form_delta(30, 5e-3, 0.1e-3, 4e-2);
  }

  if (false) stability(map);
}
