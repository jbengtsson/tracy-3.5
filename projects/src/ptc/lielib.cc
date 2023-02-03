#define NO 6

#include "tracy_lib.h"


// F. Klein's Erlangen Program.


int no_tps   = NO,
    ndpt_tps = 5;


tps get_h_k(const tps &h, const int k)
{
  // Take in Forest's F77 LieLib.
  // Get monomials of order k.
  long int no;
  tps      h_k;

  no = getno_();
  danot_(k-1);
  h_k = -h;
  danot_(k);
  h_k += h;
  danot_(no);
  return h_k;
}


ss_vect<tps> get_M_k(const ss_vect<tps> &x, const int k)
{
  // Taked in Forest's F77 LieLib.
  int          i;
  ss_vect<tps> map_k;

  for (i = 0; i < nv_tps; i++)
    map_k[i] = get_h_k(x[i], k);
  return map_k;
}


tps tps_fun
(const tps &a, std::function<double (const long int [])> fun)
{
  // Dacfu in Forest's F77 LieLib.
  // Multiplies mononials I_vec with function f(I_vec).
  char     name[name_len_for+1];
  int      k, n;
  long int ibuf1[bufsize], ibuf2[bufsize], jj[nv_tps];
  double   rbuf[bufsize];
  tps      b;

  a.exprt(rbuf, ibuf1, ibuf2, name);
  n = rbuf[0];
  for (k = 1; k <= n; k++) {
    dehash_(no_tps, nv_tps, ibuf1[k-1], ibuf2[k-1], jj);
    rbuf[k] *= fun(jj);
  }
  b.imprt(n, rbuf, ibuf1, ibuf2);
  return b;
}


double f_int_mon(const long int jj[])
{
  // Integrate monomials:
  //   scl = 1/(|I_vec|+1)
  int    k;
  double scl;

  scl = 0e0;
  for (k = 0; k < 2*nd_tps; k++)
    scl += jj[k];
  scl += 1e0;
  scl = 1e0/scl;
  return scl;
}


tps M_to_h(const ss_vect<tps> &map)
{
  // Intd in Forest's F77 LieLib.
  // E. Forest, M. Berz, J. Irwin "Normal Form Methods for Complicated
  // Periodic Systems: A Complete Solution Using Differential Algebra and Lie
  // Operators" Part. Accel. 24, 91-107 (1989):
  //   Eqs. (34)-(37).
  // Integrate monomials:
  //   M -> exp(:h:)
  int          k;
  tps          f_x, f_px, h;
  ss_vect<tps> Id;

  Id.identity();
  h = 0e0;
  for (k = 0; k < nd_tps; k++) {
    // Integrate monomials.
    f_x = tps_fun(map[2*k+1], f_int_mon)*Id[2*k];
    f_px = tps_fun(map[2*k], f_int_mon)*Id[2*k+1];
    h += f_x - f_px;
  }
  return h;
}


ss_vect<tps> h_to_v(const tps &h)
{
  // Difd in Forest's F77 LieLib:
  // Compute vector flow operator from Lie operator :h:
  //   v = Omega * [del_x H, del_px H]^T
  int          k;
  ss_vect<tps> v;

  for (k = 0; k < nd_tps; k++) {
    v[2*k+1] = Der(h, 2*k+1);
    v[2*k] = -Der(h, 2*k+2);
  }
  return v;
}


tps v_to_tps(const ss_vect<tps> &v, const tps &x)
{
  // Daflo in Forest's F77 LieLib.
  //   y = v * nabla * x
  int k;
  tps y;

  y = 0e0;
  for (k = 0; k < 2*nd_tps; k++)
    y += v[k]*Der(x, k+1);
  return y;
}


tps exp_v_to_tps(const ss_vect<tps> &v, const tps &x, const double eps,
	      const int n_max)
{
  // Expflo in Forest's F77 LieLib:
  //   y = exp(v*nabla) * x
  int    k;
  double eps1;
  tps    y_k, y;

  y_k = y = x;
  for (k = 1; k <= n_max; k++) {
    y_k = v_to_tps(v, y_k/k);
    y += y_k;
    eps1 = abs(y_k);
    if (eps1 < eps)
      break;
  }
  if (eps1 < eps)
    return y;
  else {
    printf("\n*** exp_v_to_tps: did not converge eps = %9.3e (eps = %9.3e)"
	   " n_max = %1d\n", eps1, eps, n_max);
    return NAN;
  }
}


tps exp_v_fac_to_tps(const ss_vect<tps> &v, const tps &x, const int k1,
		      const int k2, const double scl)
{
  // Facflo in Forest's F77 LieLib.
  //   y = exp(D_k1) * exp(D_k1+1) ...  * exp(D_k2) * x
  int          k;
  tps          y;
  ss_vect<tps> v_k;

  const int n_max = 100; 

  y = x;
  for (k = k1; k <= k2; k++) {
    v_k = scl*get_M_k(v, k);
    y = exp_v_to_tps(v_k, y, eps_tps, n_max);
  }
  return y;
}


ss_vect<tps> exp_v_fac_to_M(const ss_vect<tps> &v, const ss_vect<tps> &x,
			    const int k1, const int k2, const double scl)
{
  // Facflod in Forest's F77 LieLib.
  int          k;
  ss_vect<tps> M;

  for (k = 0; k < 2*nd_tps; k++)
    M[k] = exp_v_fac_to_tps(v, x[k], k1, k2, scl);
  return M;
}


ss_vect<tps>M_to_M_fact(const ss_vect<tps> &map)
{
  // Flofac in Forest's F77 LieLib.
  // Factor map:
  //   M = M_2 ... * M_n
  int          j, k;
  tps          y;
  ss_vect<tps> map_lin, v_k, map_res, map_fact;

  const int n_max = 100; 

  map_lin = get_M_k(map, 1);
  map_res = map*Inv(map_lin);
  map_fact.zero();
  for (k = 2; k <= no_tps; k++) {
    map_fact += get_M_k(map_res, k);
    // Loop over phase-space dimensions.
    for (j = 0; j < 2*nd_tps; j++) {
      y = map_res[j];
      v_k = get_M_k(-map_fact, k);
      y = exp_v_to_tps(v_k, y, eps_tps, n_max);
      map_res[j] = y;
    }
  }
  return map_fact;
}


ss_vect<tps>M_to_M_fact2(const ss_vect<tps> &map)
{
  // Obsolete.
  // Flofac in Forest's F77 LieLib.
  // Factor map:
  //   M = M_2 ... * M_n
  int          j, k;
  ss_vect<tps> map_lin, map_res, map_fact;

  map_lin = get_M_k(map, 1);
  map_res = map*Inv(map_lin);
  map_fact.zero();
  for (k = 2; k <= no_tps; k++) {
    map_fact += get_M_k(map_res, k);
    for (j = 0; j < 2*nd_tps; j++)
      map_res[j] = exp_v_fac_to_tps(map_fact, map_res[j], k, k, -1e0);
  }
  return map_fact;
}


tps exp_h_to_tps(const tps &h, const tps &x, const double eps,
		   const int n_max)
{
  // Exp1d in Forest's F77 LieLib.
  //   y = exp(:h:) x
  return exp_v_to_tps(h_to_v(h), x, eps, n_max);
}


ss_vect<tps> exp_h_to_M(const tps &h, const ss_vect<tps> &x, const double eps,
		       const int n_max)
{
  // Expnd2 in Forest's F77 LieLib.
  //   Y = exp(:h:) X
  int          k;
  ss_vect<tps> y;

  y = x;
  for (k = 0; k < 2*nd_tps; k++)
    y[k] = exp_h_to_tps(h, y[k], eps, n_max);
  return y;
}


tps M_to_h_DF(const ss_vect<tps> &map)
{
  // Liefact in Forest's F77 LieLib.
  // A. Dragt, J. Finn "Lie Series and Invariant Functions for Analytic
  // Symplectic maps" J. Math. Phys. 17, 2215-2227 (1976).
  // Dragt-Finn factorization:
  //   M ->  M_lin * exp(:h_3:) * exp(:h_4:) ...  * exp(:h_n:)
  return M_to_h(M_to_M_fact(map));
}


ss_vect<tps> h_DF_to_M
(const tps &h_DF, const ss_vect<tps> &x, const int k1, const int k2)
{
  // Fexpo in Forest's F77 LieLib.
  // Compute map from Dragt-Finn factorisation:
  //   M = exp(:h_3:) * exp(:h_4:) ...  * exp(:h_n:) * X
  ss_vect<tps> v_DF;

  v_DF = h_to_v(h_DF);
  cout << v_DF;
  exit(0);
  return exp_v_fac_to_M(v_DF, x, k1, k2, 1e0);
}


ss_vect<tps> h_DF_to_M2
(const tps &h, const ss_vect<tps> &map, const int k1, const int k2)
{
  // Fexpo in Forest's LieLib.
  // Compute map from Dragt-Finn factorisation:
  //   exp(:h_3:) exp(:h_4:) ... exp(:h_no:)
  int          k;
  tps          h_k;
  ss_vect<tps> map1;

  map1.identity();
  for (k = k2; k >= k1; k--) {
    h_k = get_h_k(h, k);
    map1 = map1*LieExp(h_k, map);
  }
  return map1;
}


int main(int argc, char *argv[])
{
  long int    lastpos;
  tps          h, h_DF;
  ss_vect<tps> Id, M, M2, M_lin;

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

  Id.identity();

  if (!false) {
    danot_(1);
    Ring_GetTwiss(true, 0e0);
    printglob();
  }

  if (false) {
    danot_(no_tps);

    M.identity();
    Cell_Pass(0, globval.Cell_nLoc, M, lastpos);
    printf("\nM:");
    prt_lin_map(3, M);

    h_DF = LieFact_DF(M, M_lin);
    cout << exp_h_to_M(h_DF, Id, eps_tps, 100)-LieExp(h_DF, Id) << "\n";
  }

  if(!false) {
    danot_(no_tps);

    M.identity();
    Cell_Pass(0, globval.Cell_nLoc, M, lastpos);
    printf("\nM:");
    prt_lin_map(3, M);

    cout << M_to_h_DF(M)-LieFact_DF(M, M_lin) << "\n";
  }

  if(false) {
    danot_(no_tps-1);

    M.identity();
    Cell_Pass(0, globval.Cell_nLoc, M, lastpos);
    printf("\nM:");
    prt_lin_map(3, M);

    danot_(no_tps);

    M_lin = get_M_k(M, 1);
    prt_lin_map(3, M_lin);

    h_DF = M_to_h_DF(M);
    cout << "\nh_DF:\n" << h_DF << "\n";

    // M2 = h_DF_to_M(h_DF, Id, 3, no_tps)*M_lin;

    // cout << h_DF_to_M(h_DF, Id, 3, no_tps);
    // cout << FExpo(h_DF, Id, 3, no_tps, -1);
    
    cout << h_to_v(h_DF)-Difd(h_DF, -1e0);

    exit(0);

    danot_(no_tps-1);

    cout << "\nh_M2-M:\n" << 1e0*M2 << "\n";
  }
}
