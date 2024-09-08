#define NO 4

#include "tracy_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;


tps get_mns(const tps &a, const int no1, const int no2)
{
  tps  b;

  danot_(no1-1);
  b = -a;
  danot_(no2);
  b += a;
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
  MNF.A1 = putlinmat(4, globval.Ascr);
  for (auto k = 4; k < nv_tps; k++)
    MNF.A1[k] = tps(0e0, k+1); 
  R = putlinmat(4, R_mat);
  for (auto k = 4; k < nv_tps; k++)
    R[k] = tps(0e0, k+1); 

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

  map1 = map;
  R = get_map_Fl(map1);
  map1 = map1*Inv(R);

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

       M = exp(:h_no-:)...exp(:h_6-9:)exp(:h_4-5:)exp(:h_3:)R                 */

  int           j, n;
  tps           h, hn;
  ss_vect<tps>  Id, A0_inv, R, Fn, map1;

  Id.identity();

  map1 = map;
  R = get_map_Fl(map1);
  map1 = map1*Inv(R);

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
  long int      jj[ss_dim];
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

#if 0
ss_vect<tps> get_A_nl(const tps g)
{
  ss_vect<tps>  Id;

  Id.identity();

  return FExpo(g, Id, 3, no_tps, -1);
}

#else

ss_vect<tps> get_A_nl(const tps g)
{
  int          j;
  tps          gn;
  ss_vect<tps> Id, A_nl;

  Id.identity(); A_nl = Id;
  for (j = 3; j <= no_tps; j++) {
    gn = Take(g, j);
    A_nl = A_nl*LieExp(gn, Id);
  }
  return A_nl;
}
#endif

tps get_g(const tps nu_x, const tps nu_y, const tps &h)
{
  // Compute g = (1-R)^-1 * h 

  int           i, j, k, l, m;
  long int      jj1[ss_dim], jj2[ss_dim];
  double        re, im;
  tps           h_re, h_im, g_re, g_im, mn1, mn2, cotan;
  ss_vect<tps>  Id;

  CtoR(h, h_re, h_im);

  for (k = 0; k < nv_tps; k++) {
    jj1[k] = 0;
    jj2[k] = 0;
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
	      jj1[delta_] = m;
	      jj2[delta_] = m;
	      re = h_re[jj1];
	      im = h_im[jj1];
	      // Compute g.
	      g_re += (re-cotan*im)*(mn1+mn2)*pow(Id[delta_], m)/2.0;
	      g_im += (im+cotan*re)*(mn1-mn2)*pow(Id[delta_], m)/2.0;
	      h_re.pook(jj2, 0.0);
	      h_im.pook(jj2, 0.0);
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


void map_norm_JB(MNF_struct &MNF)
{
  double        nu_0[2];
  tps           hn, hn_re, hn_im, gn, Kn, g, g_re, g_im, K, K_re, K_im;
  ss_vect<tps>  Id, R, A, nus, map1, map2, A_nl;

  Id.identity();

  danot_(no_tps-1);

  map1 = map;
  R = get_map_Fl(map1);

  printf("\nA0:");
  prt_lin_map(3, MNF.A0);
  printf("\nA1:");
  prt_lin_map(3, MNF.A1);
  printf("\nR:");
  prt_lin_map(3, R);

  danot_(no_tps);

  K = 0e0;
  for (auto k = 0; k < 2; k++) {
    nu_0[k] = atan2(R[2*k][2*k+1], R[2*k][2*k]);
    if (nu_0[k] < 0.0) nu_0[k] += 2.0*M_PI;
    nu_0[k] /= 2.0*M_PI;
    K -= M_PI*nu_0[k]*(sqr(Id[2*k])+sqr(Id[2*k+1]));
  }
  cout << endl;
  cout << fixed << setprecision(5)
       << "nu_0 = [" << nu_0[X_] << ", " << nu_0[Y_] << "]" << "\n";

  // coasting beam
  K += h_ijklm(map1[ct_], 0, 0, 0, 0, 1)*sqr(Id[delta_])/2.0;

  g = 0e0;
  for (auto k = 3; k <= no_tps; k++) {
    map2 = map1*Inv(R*FExpo(K, Id, 3, k-1, -1));
    hn = Intd(get_mns(map2, k-1, k-1), -1.0);
    gn = get_g(nu_0[X_], nu_0[Y_], hn);
    g += gn;
    CtoR(hn, hn_re, hn_im);
    Kn = RtoC(get_Ker(hn_re), get_Ker(hn_im));
    K += Kn;
    A = FExpo(gn, Id, k, k, -1);
    map1 = Inv(A)*map1*A;
  }

  if (!false) {
    MNF = MapNorm(map, no_tps);

    daeps_(1e-8);
    cout << "\ng:" << 1e0*g;
    cout << "\nK:" << 1e0*K;
    cout << "\ng-MNF.g:" << g-MNF.g;
    cout << "\nK-MNF.K:" << K-MNF.K;
    daeps_(1e-30);
  }
}


ss_vect<tps> map_norm_SC(void)
{
  int           k, n;
  double        nu_0[2];
  tps           hn, hn_re, hn_im, gn, Kn, g, K;
  ss_vect<tps>  Id, R, A, nus, map1, map2;

  const bool  sup_conv = true;

  Id.identity();

  danot_(no_tps-1);

  map1 = map;
  R = get_map_Fl(map1);

  danot_(no_tps);

  K = 0.0;
  for (k = 0; k < 2; k++) {
    nu_0[k] = atan2(R[2*k][2*k+1], R[2*k][2*k]);
    if (nu_0[k] < 0.0) nu_0[k] += 2.0*M_PI;
    nu_0[k] /= 2e0*M_PI;
    K -= M_PI*nu_0[k]*(sqr(Id[2*k])+sqr(Id[2*k+1]));
  }
  cout << endl;
  cout << fixed << setprecision(5)
       << "nu_0 = (" << nu_0[X_] << ", " << nu_0[Y_] << ")" << endl;

  // coasting beam
  K += h_ijklm(map1[ct_], 0, 0, 0, 0, 1)*sqr(Id[delta_])/2e0;
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
    CtoR(hn, hn_re, hn_im);
    Kn = RtoC(get_Ker(hn_re), get_Ker(hn_im)); K += Kn;
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
    map1 = FExpo(-MNF.g, Id, 3, no_tps, 1);
    MNF.g = -LieFact_SC(map1*R);
  }

  cout << MNF.g-g;
//  CtoR(MNF.K-K, hn_re, hn_im);
//  cout << hn_re;

  return map1;
}


inline long int* vec2arr(const std::vector<long int> &vec)
{
  static long int jj[ss_dim];
  for (auto k = 0; k < vec.size(); k++)
    jj[k] = vec[k];
  return jj;
 }


void prt_h_ijklm
(const char h, const std::vector<long int> &vec, const tps &h_re,
 const tps &h_im)
{
  std::string index;

  auto jj = vec2arr(vec);
  for (auto k = 0; k < 6; k++)
    index += '0' + jj[k];
  std::cout << std::scientific << std::setprecision(16)
	    << h << "_" << index << " = [" << std::setw(23) << h_re[jj] << ", "
	    << std::setw(23) <<  h_im[jj] << "]\n";
}


void prt_drv_terms(const tps &h, const tps &K)
{
  tps h_re, h_im, K_re, K_im;

  CtoR(h, h_re, h_im);
  CtoR(K, K_re, K_im);

  // Linear chromaticity.
  std::cout << "\n";
  prt_h_ijklm('h', {1, 1, 0, 0, 1, 0, 0}, h_re, h_im);
  prt_h_ijklm('h', {0, 0, 1, 1, 1, 0, 0}, h_re, h_im);

  // First order chromatic terms.
  std::cout << "\n";
  prt_h_ijklm('h', {2, 0, 0, 0, 1, 0, 0}, h_re, h_im);
  prt_h_ijklm('h', {0, 0, 2, 0, 1, 0, 0}, h_re, h_im);
  prt_h_ijklm('h', {1, 0, 0, 0, 2, 0, 0}, h_re, h_im);

  // Normal sextupoles.
  std::cout << "\n";
  prt_h_ijklm('h', {2, 1, 0, 0, 0, 0, 0}, h_re, h_im);
  prt_h_ijklm('h', {3, 0, 0, 0, 0, 0, 0}, h_re, h_im);
  prt_h_ijklm('h', {1, 0, 1, 1, 0, 0, 0}, h_re, h_im);
  prt_h_ijklm('h', {1, 0, 2, 0, 0, 0, 0}, h_re, h_im);
  prt_h_ijklm('h', {1, 0, 0, 2, 0, 0, 0}, h_re, h_im);

  // Amplitude dependent tune shift.
  std::cout << "\n";
  prt_h_ijklm('K', {2, 2, 0, 0, 0, 0, 0}, K_re, K_im);
  prt_h_ijklm('K', {1, 1, 1, 1, 0, 0, 0}, K_re, K_im);
  prt_h_ijklm('K', {0, 0, 2, 2, 0, 0, 0}, K_re, K_im);
}


void prt_ampl_dep_orb_terms(MNF_struct &MNF)
{
  tps g_re, g_im;

  MNF.g = get_mns(MNF.g, 1, 3);
  CtoR(MNF.g, g_re, g_im);
  cout << "\ng_im:\n" << g_im;

  printf("\n  s_21000 = %10.3e\n", h_ijklm(g_im, 2, 1, 0, 0, 0));
  printf("  s_30000 = %10.3e\n", h_ijklm(g_im, 3, 0, 0, 0, 0));
  printf("  s_10110 = %10.3e\n", h_ijklm(g_im, 1, 0, 1, 1, 0));
  printf("  s_10200 = %10.3e\n", h_ijklm(g_im, 1, 0, 2, 0, 0));
  printf("  s_10020 = %10.3e\n", h_ijklm(g_im, 1, 0, 0, 2, 0));
}


void compute_ampl_dep_orb(MNF_struct &MNF)
{
  tps
    x_avg, x_avg_re, x_avg_im, x3_avg, x3_avg_re, x3_avg_im,
    x_y2_avg, x_y2_avg_re, x_y2_avg_im;
  ss_vect<tps>
    Id;

  Id.identity();

  MNF.g = get_mns(MNF.g, 1, 3);

  x_avg = PB(MNF.g, Id[x_]*MNF.A1);
  CtoR(x_avg, x_avg_re, x_avg_im);
  // cout << "\n<x_re>:" << x_avg_re;

  x3_avg = PB(MNF.g, cube(Id[x_])*MNF.A1);
  CtoR(x3_avg, x3_avg_re, x3_avg_im);
  // cout << "\n<x3_re>:" << x3_avg_re;

  x_y2_avg = PB(MNF.g, Id[x_]*sqr(Id[y_])*MNF.A1);
  CtoR(x_y2_avg, x_y2_avg_re, x_y2_avg_im);
  // cout << "\n<x_y^2_re>:" << x_y2_avg_re;

  printf("\n  s_11000 = %10.3e\n", h_ijklm(x_avg_re, 1, 1, 0, 0, 0));
  printf("  s_00110 = %10.3e\n", h_ijklm(x_avg_re, 0, 0, 1, 1, 0));
  printf("  s_22000 = %10.3e\n", h_ijklm(x3_avg_re, 2, 2, 0, 0, 0));
  printf("  s_11110 = %10.3e\n", h_ijklm(x3_avg_re, 1, 1, 1, 1, 0));
  printf("  s_11110 = %10.3e\n", h_ijklm(x_y2_avg_re, 1, 1, 1, 1, 0));
  printf("  s_00220 = %10.3e\n", h_ijklm(x_y2_avg_re, 0, 0, 2, 2, 0));
}


void set_state(void)
{
  globval.H_exact        = false;
  globval.quad_fringe    = false;
  globval.Cavity_on      = false;
  globval.radiation      = false;
  globval.emittance      = false;
  globval.IBS            = false;
  globval.pathlength     = false;
  globval.Aperture_on    = false;
  globval.Cart_Bend      = false;
  globval.dip_edge_fudge = true;
}


int main(int argc, char *argv[])
{

  // disable from TPSALib- and LieLib log messages
  idprset_(-1);

  daeps_(1e-30);

  if (false)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  danot_(1);

//  Ring_GetTwiss(true, 0.0); printglob();

  danot_(no_tps-1);

  get_map(false);

  printf("\nM:");
  prt_lin_map(3, map);

  danot_(no_tps);

  if (!false) {
    tps          h_re, h_im;
    ss_vect<tps> R;

    MNF = MapNorm(map, no_tps);

    auto M_Fl = Inv(MNF.A0*MNF.A1)*map*MNF.A0*MNF.A1;
    auto h = LieFact_DF(M_Fl, R);
    prt_drv_terms(h, MNF.K);

    // h = get_mns(h, 1, 3);
    if (false) {
      CtoR(h, h_re, h_im);
      cout << "\nh_re:" << h_re;
      cout << "\nh_im:" << h_im;
    }
  }
  assert(false);

  map_norm_JB(MNF);
  // map_norm_SC();

  cout << "\ng:" << MNF.g;
  cout << "\nK:" << MNF.K;

  prt_ampl_dep_orb_terms(MNF);
  // compute_ampl_dep_orb(MNF);
}
