#define NO 7

#include "tracy_lib.h"

int  no_tps   = NO,
     ndpt_tps = 5;


void prt_lin_map(const int n_DOF, const ss_vect<tps> &map)
{
  int  i, j;

  cout << endl;
  for (i = 1; i <= 2*n_DOF; i++) {
    for (j = 1; j <= 2*n_DOF; j++)
      if (true) 
	cout << scientific << setprecision(6)
	     << setw(15) << getmat(map, i, j);
      else
	cout << scientific << setprecision(16)
	     << setw(24) << getmat(map, i, j);
    cout << endl;
  }
}


tps get_mns(const tps &a, const int no1, const int no2)
{
  tps  b;

  danot_(no1-1); b = -a;
  danot_(no2);   b += a;
  danot_(no_tps);

  return b;
}


ss_vect<tps> get_mns(const ss_vect<tps> &x, const int no1, const int no2)
{
  int           k;
  ss_vect<tps>  y;

  for (k = 0; k < nv_tps; k++)
    y[k] = get_mns(x[k], no1, no2);

  return y;
}


void get_A0(const ss_vect<tps> &map)
{
  ss_vect<tps>  Id, map1, dx, eta, A0, M_inv;

  Id.zero(); Id[delta_] = tps(0.0, delta_+1); dx = map*Id;

  map1 = map - dx; map1[delta_] = 0.0; map1[ct_] = 0.0;

  Id.identity(); eta = Inv(Id-map1)*dx;

  MNF.A0.identity(); MNF.A0 += eta;
}


ss_vect<tps> get_map_Fl(ss_vect<tps> &map)
{
  Matrix        R_mat;
  ss_vect<tps>  A0_inv, R;

  const int  n_dim = 4;

  // find fixed point
  GoFix(map, MNF.A0, A0_inv, no_tps);

  // translate to fix point
  map = A0_inv*map*MNF.A0;

  getlinmat(n_dim, map, globval.OneTurnMat);
  GDiag(n_dim, 1.0, globval.Ascr, globval.Ascrinv, R_mat,
        globval.OneTurnMat, globval.Omega, globval.Alphac);
  MNF.A1.identity(); putlinmat(4, globval.Ascr, MNF.A1);
  R.identity(); putlinmat(4, R_mat, R);

  // transform to Floquet space
  map = Inv(MNF.A1)*map*MNF.A1;

  return R;
}


tps LieFact_JB(const ss_vect<tps> &map)
{
  /* Dragt-Finn factorization:

       M = exp(:h_no:)...exp(:h_4:)exp(:h_3:)R                              */

  int           k;
  tps           h, hn;
  ss_vect<tps>  Id, A0_inv, R, Fn, map1;

  Id.identity();

  map1 = map; R = get_map_Fl(map1); map1 = map1*Inv(R);

  h = 0.0;
  for (k = 3; k <= no_tps; k++) {
    Fn = get_mns(map1, k-1, k-1); hn = Intd(Fn, -1.0); h += hn;
    map1 = map1*FExpo(-hn, Id, k, k, -1);
  }

//  cout << h;

  return h;
}


tps LieFact_SC(const ss_vect<tps> &map)
{
  /* Superconvergent Dragt-Finn factorization:

       M = exp(:h_no,:)...exp(:h_6,9:)exp(:h_4,5:)exp(:h_3:)R                */

  int           j, n;
  tps           h, hn;
  ss_vect<tps>  Id, A0_inv, R, Fn, map1;

  Id.identity();

  map1 = map; R = get_map_Fl(map1); map1 = map1*Inv(R);

  h = 0; n = 1;
  while (n+2 < no_tps) {
    Fn.zero();
    for (j = n+2; j <= 2*n+1; j++)
      Fn += get_mns(map1, j-1, j-1);
    hn = Intd(Fn, -1.0); h += hn;
    map1 = map1*LieExp(-hn, Id);
    n *= 2;
  }

  return h;
}


tps get_Ker(const tps &h)
{
  int           i, j, k;
  iVector       jj;
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


tps get_Img(const tps &h) { return h-get_Ker(h); }


ss_vect<tps> get_R_nl(void)
{
  ss_vect<tps>  Id, R_nl;

  Id.identity(); R_nl = LieExp(MNF.K, Id);

  return R_nl;
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


tps get_g(const tps nu_x, const tps nu_y, const tps &h)
{
  // Compute g = (1-R)^-1 * h 

  int           i, j, k, l, m;
  iVector       jj1, jj2;
  double        re, im;
  tps           h_re, h_im, g_re, g_im, mn1, mn2, cotan;
  ss_vect<tps>  Id;

  CtoR(h, h_re, h_im);

  for (k = 0; k < nv_tps; k++) {
    jj1[k] = 0; jj2[k] = 0;
  }

  Id.identity(); g_re = 0.0; g_im = 0.0;
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
	      jj1[delta_] = m; jj2[delta_] = m;
	      re = h_re[jj1]; im = h_im[jj1];
	      // compute g
	      g_re += (re-cotan*im)*(mn1+mn2)*pow(Id[delta_], m)/2.0;
	      g_im += (im+cotan*re)*(mn1-mn2)*pow(Id[delta_], m)/2.0;
	      h_re.pook(jj2, 0.0); h_im.pook(jj2, 0.0);
	    }
	  }
	}
      }
    }
  }

  return RtoC(g_re, g_im);
}


int pow(const int i, const int n)
{
  int  j, k;

  if (n < 0) {
    cout << "pow ***: neg exponent " << n;
    exit(1);
  }

  j = 1;
  for (k = 1; k <= n; k++)
    j *= i;

  return j;
}


ss_vect<tps> map_norm(void)
{
  int           k, n;
  double        nu0[2];
  tps           hn, hn_re, hn_im, h_ke, gn, Kn, g, K;
  ss_vect<tps>  Id, R, A, nus, map1, map2;

  Id.identity();

  danot_(no_tps-1);

  map1 = map; R = get_map_Fl(map1);

  danot_(no_tps);

  K = 0.0;
  for (k = 0; k < 2; k++) {
    nu0[k] = atan2(R[2*k][2*k+1], R[2*k][2*k]);
    if (nu0[k] < 0.0) nu0[k] += 2.0*M_PI;
    nu0[k] /= 2.0*M_PI;
    K -= M_PI*nu0[k]*(sqr(Id[2*k])+sqr(Id[2*k+1]));
  }
  cout << endl;
  cout << fixed << setprecision(5)
       << "nu0 = (" << nu0[X_] << ", " << nu0[Y_] << ")" << endl;

  // coasting beam
  K += h_ijklm(map1[ct_], 0, 0, 0, 0, 1)*sqr(Id[delta_])/2.0;
//  CtoR(K, hn_re, hn_im);
//  cout << endl << "K:" << hn_re;

  g = 0.0;
  for (k = 3; k <= no_tps; k++) {
    n = pow(2, k-3);

    map2 = map1*Inv(R*FExpo(K, Id, 3, k-1, -1));
    hn = Intd(get_mns(map2, k-1, k-1), -1.0); gn = get_g(nu0[X_], nu0[Y_], hn);
    g += gn;
    CtoR(hn, hn_re, hn_im); Kn = RtoC(get_Ker(hn_re), get_Ker(hn_im)); K += Kn;
//    cout << endl << "k = " << k << hn_re;

    A = FExpo(gn, Id, k, k, -1); map1 = Inv(A)*map1*A;
  }

  MNF = MapNorm(map, no_tps);

  cout << MNF.g-g;
//  CtoR(MNF.K-K, hn_re, hn_im);
//  cout << hn_re;

  return map1;
}


ss_vect<tps> map_norm_SC(void)
{
  int           k, n;
  double        nu0[2];
  tps           hn, hn_re, hn_im, h_ke, gn, Kn, g, K;
  ss_vect<tps>  Id, R, A, nus, map1, map2;

  const bool  sup_conv = true;

  Id.identity();

  danot_(no_tps-1);

  map1 = map; R = get_map_Fl(map1);

  danot_(no_tps);

  K = 0.0;
  for (k = 0; k < 2; k++) {
    nu0[k] = atan2(R[2*k][2*k+1], R[2*k][2*k]);
    if (nu0[k] < 0.0) nu0[k] += 2.0*M_PI;
    nu0[k] /= 2.0*M_PI;
    K -= M_PI*nu0[k]*(sqr(Id[2*k])+sqr(Id[2*k+1]));
  }
  cout << endl;
  cout << fixed << setprecision(5)
       << "nu0 = (" << nu0[X_] << ", " << nu0[Y_] << ")" << endl;

  // coasting beam
  K += h_ijklm(map1[ct_], 0, 0, 0, 0, 1)*sqr(Id[delta_])/2.0;
//  CtoR(K, hn_re, hn_im);
//  cout << endl << "K:" << hn_re;

  g = 0.0;
  for (k = 3; k <= no_tps; k++) {
    n = pow(2, k-3);
    if (n+2 > no_tps) break;

    map2 = map1*Inv(R*FExpo(K, Id, 3, k-1, -1));
    map2 = get_mns(map2, n+1, min(2*n, no_tps-1)); hn = Intd(map2, -1.0);
    // get_g: uses resonance basis
    nus = dHdJ(K);
    gn = get_mns(get_g(nus[3], nus[4], hn), n+2, min(2*n+1, no_tps));
    g += gn;
    CtoR(hn, hn_re, hn_im); Kn = RtoC(get_Ker(hn_re), get_Ker(hn_im)); K += Kn;
//    cout << endl << "k = " << k << Kn << hn_re;

    A = LieExp(gn, Id); map1 = Inv(A)*map1*A;

//    map2 = get_mns(map1, n+1, min(2*n, no_tps-1)); hn = Intd(map2, -1.0);
//    CtoR(hn, hn_re, hn_im);
//    cout << endl << hn_re;
  }

  if (!sup_conv) {
    // transform to regular Map Normal Form
    gn = 0.0;
    for (k = 3; k <= no_tps; k++) {
      n = pow(2, k-3);
      if (n+2 > no_tps) break;
      hn = get_mns(g, n+2, min(2*n+1, no_tps));
      // reversed order
      gn += -LieFact_JB(LieExp(-hn, Id)*R);
    }
    g = gn;

//    cout << g;
  }

  MNF = MapNorm(map, no_tps);

  if (sup_conv) {
    // transform to super convergent Map Normal Form
    map1 = FExpo(-MNF.g, Id, 3, no_tps, 1); MNF.g = -LieFact_SC(map1*R);
  }

  cout << MNF.g-g;
//  CtoR(MNF.K-K, hn_re, hn_im);
//  cout << hn_re;

  return map1;
}


int main(int argc, char *argv[])
{

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  // disable from TPSALib- and LieLib log messages
  idprset_(-1);

  rdmfile("flat_file.dat");

  danot_(1);

//  Ring_GetTwiss(true, 0.0); printglob();

  danot_(no_tps-1);

  get_map();

  prt_lin_map(3, map);

  danot_(no_tps);

//  map_norm();
  map_norm_SC();
}
