#define NO 4

#include "tracy_lib.h"


// F. Klein's Erlangen Program.

int no_tps   = NO,
    ndpt_tps = 5;


ss_vect<tps> Id(void)
{
  ss_vect<tps> Id;
  Id.identity();
  return Id;
}


tps get_h_k(const tps &h, const int k1, const int k2)
{
  // Get monomials of order [k1..k2].
  long int no;
  tps      h_k;

  no = getno_();
  danot_(k1-1);
  h_k = -h;
  danot_(k2);
  h_k += h;
  danot_(no);
  return h_k;
}


ss_vect<tps> get_map_k(const ss_vect<tps> &x, const int k1, const int k2)
{
  int          k;
  ss_vect<tps> map_k;

  for (k = 0; k < nv_tps; k++)
    map_k[k] = get_h_k(x[k], k1, k2);
  return map_k;
}


ss_vect<tps> get_map_Fl(ss_vect<tps> &map)
{
  Matrix       R_mat;
  ss_vect<tps> Id, R;

  const int n_dim = 4;

  Id.identity();
  // Find fix point.
  GoFix(map, MNF.A0, MNF.A0_inv, no_tps);
  printf("\nA0:\n");
  prt_lin_map(3, MNF.A0);
  // Translate to fix point.
  map = MNF.A0_inv*map*MNF.A0;
  getlinmat(n_dim, map, globval.OneTurnMat);
  GDiag(n_dim, 1e0, globval.Ascr, globval.Ascrinv, R_mat,
        globval.OneTurnMat, globval.Omega, globval.Alphac);
  MNF.A1 = putlinmat(n_dim, globval.Ascr);
  MNF.A1[ct_]    = Id[ct_];
  MNF.A1[delta_] = Id[delta_];
  R = putlinmat(n_dim, R_mat);
  R[ct_]    = Id[ct_];
  R[delta_] = Id[delta_];
  // Transform to Floquet space.
  printf("\nA1:\n");
  prt_lin_map(3, MNF.A1);
  map = Inv(MNF.A1)*map*MNF.A1;
  return R;
}


ss_vect<tps> compute_Dragt_Finn_Map
(const tps &h, const ss_vect<tps> &map, const int k1, const int k2,
 const int reverse)
{
  // Dragt-Finn map:
  // not reverse:
  //   exp(:h_3:) exp(:h_4:) ... exp(:h_no:)
  // reverse:
  //   exp(:h_no:) exp(:h_no-1:) ... exp(:h_3:)
  int          k;
  tps          h_k;
  ss_vect<tps> map1;

  map1.identity();
  for (k = k2; k >= k1; k--) {
    h_k = get_h_k(h, k, k);
    if (!reverse)
      map1 = map1*LieExp(h_k, map);
    else
      map1 = LieExp(h_k, map)*map1;
  }
  return map1;
}


tps LieFact_JB(const ss_vect<tps> &map)
{
  // Dragt-Finn factorization:
  //   M = M_1 exp(:h_3:) ... exp(:h_4:) exp(:h_no:).
  int          k;
  tps          h, h_k;
  ss_vect<tps> Id, A0_inv, M_1, Fn, map1;

  Id.identity();
  map1 = map;
  M_1 = get_map_Fl(map1);
  map1 = map1*Inv(M_1);
  h = 0e0;
  for (k = 3; k <= no_tps; k++) {
    Fn = get_map_k(map1, k-1, k-1);
    h_k = Intd(Fn, -1e0);
    h += h_k;
    map1 = map1*compute_Dragt_Finn_Map(-h_k, Id, k, k, true);
  }
//  cout << h;
  return h;
}


tps get_Ker(const tps &h)
{
  int          i, j, k;
  long int     jj[ss_dim];
  tps          h_Ke;
  ss_vect<tps> Id;

  for (k = 0; k < nv_tps; k++)
    jj[k] = 0;

  Id.identity();
  h_Ke = 0e0;
  for (i = 0; i <= no_tps; i++) {
    jj[x_] = i; jj[px_] = i;
    for (j = 0; j <= no_tps; j++) {
      jj[y_] = j;
      jj[py_] = j;
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
  // Compute g = (1-R)^-1 * h 

  int          i, j, k, l, m;
  long int     jj1[ss_dim], jj2[ss_dim];
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
	    cotan = 1e0/tan(((i-j)*nu_x+(k-l)*nu_y)*M_PI);
	    mn1 =
	      pow(Id[x_], i)*pow(Id[px_], j)*pow(Id[y_], k)*pow(Id[py_], l);
	    mn2 =
	      pow(Id[x_], j)*pow(Id[px_], i)*pow(Id[y_], l)*pow(Id[py_], k);

	    for (m = 0; m <= no_tps-i-j-k-l; m++) {
	      jj1[delta_] = m; jj2[delta_] = m;
	      re = h_re[jj1]; im = h_im[jj1];
	      // compute g
	      g_re += (re-cotan*im)*(mn1+mn2)*pow(Id[delta_], m)/2e0;
	      g_im += (im+cotan*re)*(mn1-mn2)*pow(Id[delta_], m)/2e0;
	      h_re.pook(jj2, 0e0); h_im.pook(jj2, 0e0);
	    }
	  }
	}
      }
    }
  }

  return RtoC(g_re, g_im);
}


ss_vect<tps> compute_map_normal_form(ss_vect<tps> &map)
{
  // Transform map into normal form:
  //   M = (exp(:g_no:) exp(:g_no-1:) ... exp(g_3) A_1 A_0)^-1
  //       exp(:K_2:) exp(:K_3:) ... exp(:K_no:)
  //       exp(:g_no:) exp(:g_no-1:) ... exp(g_3) A_1 A_0
  int          k, n;
  double       nu0[2];
  tps          h_k, h_k_re, h_k_im, h_ke, K, K_k, g, g_k;
  ss_vect<tps> Id, map1, R, A, nus, map2;

  Id.identity();

  danot_(no_tps-1);

  map1 = map;
  R = get_map_Fl(map1);

  danot_(no_tps);

  K = 0e0;
  for (k = 0; k < 2; k++) {
    nu0[k] = atan2(R[2*k][2*k+1], R[2*k][2*k]);
    if (nu0[k] < 0e0) nu0[k] += 2e0*M_PI;
    nu0[k] /= 2e0*M_PI;
    K -= M_PI*nu0[k]*(sqr(Id[2*k])+sqr(Id[2*k+1]));
  }
  cout << endl;
  cout << fixed << setprecision(5)
       << "nu0 = (" << nu0[X_] << ", " << nu0[Y_] << ")" << endl;

  // coasting beam
  K += h_ijklm(map1[ct_], 0, 0, 0, 0, 1)*sqr(Id[delta_])/2e0;
//  CtoR(K, h_k_re, h_k_im);
//  cout << endl << "K:" << h_k_re;

  g = 0e0;
  for (k = 3; k <= no_tps; k++) {
    n = pow(2, k-3);

    map2 = map1*Inv(R*compute_Dragt_Finn_Map(K, Id, 3, k-1, true));
    h_k =
      Intd(get_map_k(map2, k-1, k-1), -1e0);
    g_k = get_g(nu0[X_], nu0[Y_], h_k);
    g += g_k;
    CtoR(h_k, h_k_re, h_k_im);
    K_k = RtoC(get_Ker(h_k_re), get_Ker(h_k_im));
    K += K_k;
//    cout << endl << "k = " << k << h_k_re;

    A = compute_Dragt_Finn_Map(g_k, Id, k, k, true);
    map1 = Inv(A)*map1*A;
  }

#if 1
  MNF = MapNorm(map, no_tps);

  daeps_(1e-8);
  cout << "\nMNF.g-g:\n" << MNF.g-g;
  cout << "\nMNF.K-K:\n" << MNF.K-K;
  daeps_(eps_tps);
#endif

  // MNF.g = g;
  // MNF.K = K;

  return map1;
}


ss_vect<tps> get_A_nl(const tps g)
{
  return compute_Dragt_Finn_Map(g, Id(), 3, no_tps, true);
}


ss_vect<tps> get_A_nl_inv(const tps g)
{
  return compute_Dragt_Finn_Map(-g, Id(), 3, no_tps, false);
}


tps get_H(const tps &g, const tps &K)
{
  return K*get_A_nl_inv(g);
}


tps g_renorm(const double nu0_x, const double nu0_y,
	     const double nu1_x, const double nu1_y,
	     const tps &g)
{
  // Renormalize g: (1-R^-1)^-1 * h.
  long int     jj1[ss_dim], jj2[ss_dim];
  int          i, j, k, l, m;
  double       re, im, cotan0, cotan1, cotan0_sqr;
  tps          h_re, h_im, g_re, g_im, G_re, G_im, mn1, mn2;
  ss_vect<tps> Id;

  CtoR(g, g_re, g_im);

  for (k = 0; k < ss_dim; k++) {
    jj1[k] = 0; jj2[k] = 0;
  }

  Id.identity(); G_re = 0e0; G_im = 0e0;
  for (i = 0; i <= no_tps; i++) {
    jj1[x_] = i; jj2[px_] = i;
    for (j = 0; j <= i; j++) {
      jj1[px_] = j; jj2[x_] = j;
      for (k = 0; k <= no_tps; k++) {
	jj1[y_] = k; jj2[py_] = k;
	for (l = 0; l <= no_tps; l++) {
	  jj1[py_] = l; jj2[y_] = l;

	  if (i+j+k+l <= no_tps) {
	    cotan0 = 1e0/tan(((i-j)*nu0_x+(k-l)*nu0_y)*M_PI);
	    cotan0_sqr = sqr(cotan0);
	    cotan1 = 1e0/tan(((i-j)*nu1_x+(k-l)*nu1_y)*M_PI);
	    mn1 =
	      pow(Id[x_], i)*pow(Id[px_], j)*pow(Id[y_], k)*pow(Id[py_], l);
	    mn2 =
	      pow(Id[x_], j)*pow(Id[px_], i)*pow(Id[y_], l)*pow(Id[py_], k);

	    for (m = 0; m <= no_tps; m++) {
	      if (i+j+k+l+m <= no_tps) {
		jj1[delta_] = m; jj2[delta_] = m;
		if ((i != j) || (k != l)) {
		  re = g_re[jj1]; im = g_im[jj1];

		  // Compute h.
		  h_re = (re+cotan0*im)*2e0/(1e0+cotan0_sqr);
		  h_im = (im-cotan0*re)*2e0/(1e0+cotan0_sqr);

		  // Renormalize g.
		  G_re += (h_re-cotan1*h_im)*(mn1+mn2)*pow(Id[delta_], m)/2e0;
		  G_im += (h_im+cotan1*h_re)*(mn1-mn2)*pow(Id[delta_], m)/2e0;
		  g_re.pook(jj2, 0e0); g_im.pook(jj2, 0e0);
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


void track_H(const string &file_name, const double A_x, const double A_y,
	     tps &H)
{
  long int        lastpos;
  int             k;
  double          h, J[2], phi[2], nu[2], scl, J2[2], phi2[2];
  tps             scl_tps;
  ss_vect<double> ps, ps_Fl, ps_nl;
  FILE            *outf;

  const double x_ampl[] = {A_x, A_y};

  danot_(no_tps);

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
    J[X_] = (sqr(ps_Fl[x_])+sqr(ps_Fl[px_]))/2e0;
    phi[X_] = atan2(ps_Fl[px_], ps_Fl[x_]);
    J[Y_] = (sqr(ps_Fl[y_])+sqr(ps_Fl[py_]))/2e0;
    phi[Y_] = atan2(ps_Fl[py_], ps_Fl[y_]);

    h = (H*ps_Fl).cst();

    ps_nl = (MNF.A_nl_inv*ps_Fl).cst();
    J2[X_] = (sqr(ps_nl[x_])+sqr(ps_nl[px_]))/2e0;
    phi2[X_] = atan2(ps_nl[px_], ps_nl[x_]);
    J2[Y_] = (sqr(ps_nl[y_])+sqr(ps_nl[py_]))/2e0;
    phi2[Y_] = atan2(ps_nl[py_], ps_nl[y_]);

    fprintf(outf, "%3d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e"
	    " %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	    k, 1e3*ps_Fl[x_], 1e3*ps_Fl[px_], 1e3*ps_Fl[y_], 1e3*ps_Fl[py_],
	    1e6*J[X_], phi[X_], 1e6*J[Y_], phi[Y_], 1e6*J2[X_], phi2[X_],
	    1e6*J2[Y_], phi2[Y_], 1e6*fabs(h));
  }

  fclose(outf);
}


void track_H(const ss_vect<tps> &map)
{
  tps H;

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

  track_H("track_H_1.dat", 0.1e-3, 0.1e-3, H);
  track_H("track_H_2.dat", 0.5e-3, 0.5e-3, H);
  track_H("track_H_3.dat", 1.0e-3, 1.0e-3, H);
  track_H("track_H_4.dat", 1.1e-3, 1.0e-3, H);
  track_H("track_H_5.dat", 1.2e-3, 1.0e-3, H);
  track_H("track_H_6.dat", 1.3e-3, 1.0e-3, H);
  track_H("track_H_7.dat", 1.4e-3, 1.0e-3, H);
}


void compute_invariant(ss_vect<tps> &M)
{
  int          k;
  double       nu_int;
  tps          H_2, DH_2, H, DH;
  ss_vect<tps> Id, B;

  const double
    alpha[] = {Cell[0].Alpha[X_], Cell[0].Alpha[Y_]},
    beta[]  = {Cell[0].Beta[X_],  Cell[0].Beta[Y_]},
    gamma[] = {(1e0+sqr(alpha[X_]))/beta[X_], (1e0+sqr(alpha[Y_]))/beta[Y_]},
    eta_x   = Cell[0].Eta[X_],
    etap_x  = Cell[0].Etap[X_];

  Id.identity();

  printf("\n  alpha_x = %6.3f beta_x = %5.3f gamma_x = %5.3f"
	 " eta_x = %10.3e eta'_x = %10.3e\n"
	 "  alpha_y = %6.3f beta_y = %5.3f gamma_y = %5.3f\n",
	 alpha[X_], beta[X_], gamma[X_], eta_x, etap_x,
	 alpha[Y_], beta[Y_], gamma[Y_]);

  danot_(2);

  H_2 = 0e0;
  for (k = 0; k < 2; k++) {
    H_2 +=
      -M_PI*modf(globval.TotalTune[k], &nu_int)
      *(gamma[k]*sqr(Id[2*k])+2e0*alpha[k]*Id[2*k]*Id[2*k+1]
	+beta[k]*sqr(Id[2*k+1]));
  }
  H_2 += M[ct_][delta_]*sqr(Id[delta_])/2e0;

  B.identity();
  B[x_]  += eta_x*Id[delta_];
  B[px_] += etap_x*Id[delta_];
  B[ct_] += -etap_x*Id[x_] - eta_x*Id[px_];

  printf("\nLinear dispersion computed by numerical differentiation.\n");
  printf("\nB:\n");
  prt_lin_map(3, B);

  danot_(no_tps);

  printf("\nLinear dispersion computed by TPSA.\n");
  MNF = MapNorm(M, 1);

  cout << scientific << setprecision(3) << "\ng:\n" << MNF.g;
  printf("\nA0:\n");
  prt_lin_map(3, MNF.A0);

  H_2 = H_2*Inv(MNF.A0);

  daeps_(1e-10);
  H_2 = 1e0*H_2;
  cout << scientific << setprecision(3) << "\nH_2:\n" << H_2;
  daeps_(eps_tps);

  printf("\ne^-H_2:\n");
  prt_lin_map(3, LieExp(H_2, Id));

  // Restore max order after call to LieExp.
  danot_(2);

  DH_2 = H_2*M - H_2;
  daeps_(1e-13);
  DH_2 = 1e0*DH_2;
  cout << scientific << setprecision(3) << "\nH_2*M - H_2:\n" << DH_2;
  daeps_(eps_tps);

  danot_(no_tps);

  MNF = MapNorm(M, no_tps);

  H = get_H(MNF.g, MNF.K)*Inv(MNF.A0*MNF.A1);

  daeps_(1e-7);
  H = 1e0*H;
  cout << scientific << setprecision(3) << "\nH:\n" << H;
  daeps_(eps_tps);

  DH = H*M - H;
  daeps_(1e-7);
  DH = 1e0*DH;
  cout << scientific << setprecision(3) << "\nDH:\n" << DH;
  daeps_(eps_tps);
}


int main(int argc, char *argv[])
{
  long int    lastpos;
  ss_vect<tps> M;

  globval.H_exact    = false; globval.quad_fringe    = false;
  globval.Cavity_on  = false; globval.radiation      = false;
  globval.emittance  = false; globval.IBS            = false;
  globval.pathlength = false; globval.bpm            = 0;
  globval.Cart_Bend  = false; globval.dip_edge_fudge = true;
  globval.mat_meth   = false;

  // disable from TPSALib- and LieLib log messages
  idprset_(-1);

  if (false)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  Ring_GetTwiss(true, 0e0); printglob();

  danot_(no_tps-1);

  M.identity();
  Cell_Pass(0, globval.Cell_nLoc, M, lastpos);
  printf("\nM:\n");
  prt_lin_map(3, M);
  // cout << scientific << setprecision(3) << "\nM:\n" << M << "\n";

  if (!false)
    compute_invariant(M);

  if (false)
    track_H(M);
}
