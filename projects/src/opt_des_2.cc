#define NO 1

#include "tracy_lib.h"

#include "Powell/src/newuoa.h"

#include "drv_terms.cc"

int no_tps = NO;


const bool
  ps_rot          = !false, // Note, needs to be zeroed; after use.
  qf6_rb          = false,
  sp_std          = true,
  pert_dip_cell   = false,
  dphi            = true,
  long_grad_dip[] = {true, true};

/* opt_case:
     opt_mI_std  1,
     match_ls    2,
     opt_mi_sp   3,
     match_ss    4,
     opt_mult    5.                                                           */
const int opt_case = 3;

const int
  n_ic   = 5,
  n_cell = 6;

const double
  ic[n_ic][2] =
  {{-0.0000000016, -0.0000000198},
   {3.9407179898, 0.9733218045},
   {0.0183249406, 0.0000000000},
   {0.0, 0.0}},                    // From Center of Mid Straight:
                                   // alpha, beta, eta, eta'.
 
  eps0_x                = 0.097,

  beta_inj[]            = {10.7, 6.5},
  A_max[]               = {3.5e-3, 1.5e-3},
  delta_max             = 2e-2,
  A_delta_max[]         = {2e-3, 0.1e-3},

  dnu[]                 = {-0.2, -0.05},
  nu[]                  = {65.90/n_cell+dnu[X_], 20.73/n_cell+dnu[Y_]},
  nu_ref[]              = {nu[X_]-(int)nu[X_], nu[Y_]-(int)nu[Y_]},
  beta_ms[]             = { 4.0, 2.0},
  beta_ss[]             = { 2.0, 2.0},
  beta_ls[]             = {12.0, 5.0},
  eta_ms_x              = 19e-3,

  nu_ref_scl            = 0*5e7,

  scl_ksi_1             = 1e0*1e1,
  scl_h_3               = 1e0*1e10,
  scl_h_3_delta         = 1e0*1e10,
  scl_h_4               = 1e0,
  scl_ksi_2             = 1e0*1e5,
  scl_ksi_3             = 0*1e-5,
  scl_chi_2             = 1e0*1e5,
  scl_chi_delta_2       = 1e0*1e5,

  scl_extra             = 1e1,

  twoJ[]                =
    {sqr(A_max[X_])/beta_inj[X_], sqr(A_max[Y_])/beta_inj[Y_]},
  twoJ_delta[]          =
    {sqr(A_delta_max[X_])/beta_inj[X_], sqr(A_delta_max[Y_])/beta_inj[Y_]},
  mI_nu_ref[]           = {1.5, 0.5},
  high_ord_achr_nu[][2] = {{5.0/8.0, 2.0/8.0}, {4.0/8.0, 1.0/8.0}};


// Needs scl_extra.
#include "ctrl_H_2.cc"


void f_der(double *bn, double *df)
{
  lat_prms.set_prm(bn);
  if (lat_constr.phi_scl != 0e0) phi_corr(lat_constr);
  lat_constr.get_dchi2(twoJ, delta_max, twoJ_delta, bn, df);
}


double f_match(double *bn)
{
  double      chi2;
  static bool first = true;

  const int n_prt = 5;

  lat_prms.set_prm(bn);

  eps_x = get_lin_opt(lat_constr);

  if ((lat_constr.high_ord_achr_scl[X_] != 0e0)
      || (lat_constr.high_ord_achr_scl[Y_] != 0e0))
    get_high_ord_achr(lat_constr);

  chi2 =
    lat_constr.get_chi2(twoJ, delta_max, twoJ_delta, lat_prms, bn, false,
			false);

  if (first || (chi2 < lat_constr.chi2)) {
    first = false;
    if (lat_constr.n_iter % n_prt == 0) {
      // Print dchi2.
      lat_constr.get_chi2(twoJ, delta_max, twoJ_delta, lat_prms, bn, false,
			  true);
      prt_f(bn, chi2, lat_constr, lat_prms, false);
    }
    lat_constr.n_iter++;
    lat_constr.chi2 = chi2;
  }

  return chi2;
}


double f_achrom(double *bn)
{
  bool        stable;
  double      chi2;
  static bool first = true;

  const int n_prt = 5;

  lat_prms.set_prm(bn);

  if (lat_constr.phi_scl != 0e0) phi_corr(lat_constr);

  if (lat_constr.ring) stable = get_nu(lat_constr.nu);

  if ((lat_constr.ring && stable) || !lat_constr.ring) {
    eps_x = get_lin_opt(lat_constr);

    fit_ksi1(lat_constr.Fnum_b3, 0e0, 0e0, 1e1);

    chi2 =
      lat_constr.get_chi2(twoJ, delta_max, twoJ_delta, lat_prms, bn, true,
			  false);

    if (first || (chi2 < lat_constr.chi2)) {
      first = false;
      if (lat_constr.n_iter % n_prt == 0) {
	// Print dchi2.
	lat_constr.get_chi2(twoJ, delta_max, twoJ_delta, lat_prms, bn, true,
			    true);
	prt_f(bn, chi2, lat_constr, lat_prms, true);
      }

      lat_constr.n_iter++;
      lat_constr.chi2 = chi2;
    }
  } else
    chi2 = 1e30;

  return chi2;
}


double f_mult(double *bn)
{
  bool        stable;
  double      chi2;
  static bool first = true;

  const int n_prt = 5;

  lat_prms.set_prm(bn);

  fit_ksi1(lat_constr.Fnum_b3, 0e0, 0e0, 1e1);

  chi2 =
    lat_constr.get_chi2(twoJ, delta_max, twoJ_delta, lat_prms, bn, true, false);

  if (first || (chi2 < lat_constr.chi2)) {
    first = false;
    if (lat_constr.n_iter % n_prt == 0) {
      // Print dchi2.
      lat_constr.get_chi2(twoJ, delta_max, twoJ_delta, lat_prms, bn, true,
			  true);
      prt_f(bn, chi2, lat_constr, lat_prms, true);
    }

    lat_constr.n_iter++;
    lat_constr.chi2 = chi2;
  }

  return chi2;
}


double f_nu(double *bn)
{
  bool        stable;
  double      chi2, ksi1[2], s;

  static bool first = true;

  const int    n_prt   = 5;
  const double twoJ1[] = {1e0, 1e0};

  lat_prms.set_prm(bn);

  if (lat_constr.phi_scl != 0e0) phi_corr(lat_constr);

  if (lat_constr.ring) stable = get_nu(lat_constr.nu);

  if ((lat_constr.ring && stable) || !lat_constr.ring) {
    eps_x = get_lin_opt(lat_constr);

    fit_ksi1(lat_constr.Fnum_b3, 0e0, 0e0, 1e1);

    // Get linear chromaticity.
    lat_constr.drv_terms.h_c.clear(); lat_constr.drv_terms.h_s.clear();
    sxt_1(-1e0/(4e0*M_PI), twoJ1, 1e0, 1, 1, 0, 0, 1, ksi1[X_], s, false);
    sxt_1(1e0/(4e0*M_PI), twoJ1, 1e0, 0, 0, 1, 1, 1, ksi1[Y_], s, false);
    lat_constr.drv_terms.h_label.push_back("ksi1   ");
    lat_constr.drv_terms.h_c.push_back(ksi1[X_]);
    lat_constr.drv_terms.h_s.push_back(ksi1[Y_]);

    chi2 =
      lat_constr.get_chi2(twoJ, delta_max, twoJ_delta, lat_prms, bn, false,
			  false);

    if (first || (chi2 < lat_constr.chi2)) {
      first = false;
      if (lat_constr.n_iter % n_prt == 0) {
	// Print dchi2.
	lat_constr.get_chi2(twoJ, delta_max, twoJ_delta, lat_prms, bn, false,
			    true);
	prt_f(bn, chi2, lat_constr, lat_prms, true);
      }

      lat_constr.n_iter++;
      lat_constr.chi2 = chi2;
    }
  } else
    chi2 = 1e30;

  return chi2;
}


void set_dip_cell_std(param_type &prms, constr_type &constr)
{
  std::vector<int>    grad_dip_Fnum, mI_loc;
  std::vector<double> grad_dip_scl;

  // Standard Cell.
  grad_dip_scl.push_back(0.129665);
  grad_dip_scl.push_back(0.149256);
  grad_dip_scl.push_back(0.174737);
  grad_dip_scl.push_back(0.216027);
  grad_dip_scl.push_back(0.330314);

  grad_dip_Fnum.push_back(ElemIndex("dl1a_1"));
  grad_dip_Fnum.push_back(ElemIndex("dl1a_2"));
  grad_dip_Fnum.push_back(ElemIndex("dl1a_3"));
  grad_dip_Fnum.push_back(ElemIndex("dl1a_4"));
  grad_dip_Fnum.push_back(ElemIndex("dl1a_5"));
  if (dphi)
    prms.add_prm(grad_dip_Fnum, grad_dip_scl, -3, -15.0, 15.0, 1.0);
  if (long_grad_dip[0])
    prms.add_prm(grad_dip_Fnum, grad_dip_scl,  2, -15.0, 15.0, 1.0);

  lat_constr.Fnum_b1.push_back(-ElemIndex("dl1a_1"));
  lat_constr.grad_dip_Fnum_b1.push_back(grad_dip_Fnum);

  grad_dip_Fnum.clear();
  grad_dip_Fnum.push_back(ElemIndex("dl2a_1"));
  grad_dip_Fnum.push_back(ElemIndex("dl2a_2"));
  grad_dip_Fnum.push_back(ElemIndex("dl2a_3"));
  grad_dip_Fnum.push_back(ElemIndex("dl2a_4"));
  grad_dip_Fnum.push_back(ElemIndex("dl2a_5"));
  if (dphi)
    prms.add_prm(grad_dip_Fnum, grad_dip_scl, -3, -15.0, 15.0, 1.0);
  if (long_grad_dip[1])
    prms.add_prm(grad_dip_Fnum, grad_dip_scl,  2, -15.0, 15.0, 1.0);

  lat_constr.Fnum_b1.push_back(-ElemIndex("dl2a_1"));
  lat_constr.grad_dip_Fnum_b1.push_back(grad_dip_Fnum);

  if (dphi) {
    // Dipole Cell.
    prms.add_prm("qf4", -3, -15.0, 15.0, 1.0);
    lat_constr.Fnum_b1.push_back(ElemIndex("qf4"));
    if (qf6_rb) {
      prms.add_prm("qf6", -3, -15.0, 15.0, 1.0);
      lat_constr.Fnum_b1.push_back(ElemIndex("qf6"));
    }
    prms.add_prm("qf8", -3, -15.0, 15.0, 1.0);
    lat_constr.Fnum_b1.push_back(ElemIndex("qf8"));

    // Commented out must be defined last.
    // prms.add_prm("dq1", -3, -15.0,   15.0,  1.0);
    lat_constr.Fnum_b1.push_back(ElemIndex("dq1"));
  }

  prms.add_prm("dq1", 2, -15.0, 15.0, 1.0);

  prms.add_prm("qd3", 2, -15.0,  0.0, 1.0);
  prms.add_prm("qf4", 2,   0.0, 15.0, 1.0);
  prms.add_prm("qd5", 2, -15.0,  0.0, 1.0);

  // prms.add_prm("qd3", -1,  -0.4, -0.4, 1.0);
  // prms.add_prm("qd5", -1,  -0.4, -0.4, 1.0);
}


void set_b2_ms_std(param_type &prms)
{
  // Mid-Straight.
  prms.add_prm("qf1",  2,   0.0, 15.0, 1.0);
  prms.add_prm("qd2",  2, -15.0,  0.0, 1.0);

  // Standard-Straight.
  // prms.add_prm("qf6", -1, -0.3,  0.01, 1.0);
  prms.add_prm("qf6",  2,  0.0, 15.0, 1.0);
  prms.add_prm("qf8",  2,  0.0, 15.0, 1.0);
}


void set_constr_std(constr_type &constr)
{
  // constr.add_constr(Elem_GetPos(ElemIndex("qf1"), 1)-1,
  // 		    0e0, 0e0, 0e0, 0e0, 1e8, 1e7,
  // 		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  // constr.add_constr(Elem_GetPos(ElemIndex("qf1"), 2),
  // 		    0e0, 0e0, 0e0, 0e0, 1e8, 1e7,
  // 		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  constr.add_constr(Elem_GetPos(ElemIndex("dl1a_5"), 1)-1,
		    0e0, 0e0, 0e0, 0e0, 1e8, 1e7,
		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  constr.add_constr(Elem_GetPos(ElemIndex("dl1a_5"), 2),
		    0e0, 0e0, 0e0, 0e0, 1e8, 1e7,
		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  constr.add_constr(Elem_GetPos(ElemIndex("ms"), 1),
		    1e7, 1e7, 1e2,         1e2,         1e7,      0e0,
		    0.0, 0.0, beta_ms[X_], beta_ms[Y_], eta_ms_x, 0.0);
  constr.add_constr(Elem_GetPos(ElemIndex("ss"), 1),
		    1e5, 1e5, 1e2,         1e2,         1e7, 1e7,
		    0.0, 0.0, beta_ss[X_], beta_ss[Y_], 0.0, 0.0);
  // Symmetry.
  constr.add_constr(Elem_GetPos(ElemIndex("sf1"), 1),
		    5e2, 1e3, 0e1,  0e2, 0e0, 0e0,
		    0.0, 0.0, 10.0, 1.0, 0.0, 0.0);

  if (false)
    // Increase beta_y.
    constr.add_constr(Elem_GetPos(ElemIndex("sd1"), 1),
		      0e5, 0e5, 5e2, 5e3,  0e7, 0e7,
		      0.0, 0.0, 5.0, 15.0, 0.0, 0.0);
  if (false)
    // Increase beta_y.
    constr.add_constr(Elem_GetPos(ElemIndex("sd2"), 1),
		      0e5, 0e5, 5e2, 5e3,  0e7, 0e7,
		      0.0, 0.0, 5.0, 10.0, 0.0, 0.0);
}


void set_b3_Fam_std(param_type &prms)
{
  std::vector<int> Fnum;

  switch (2) {
  case 1:
    prms.add_prm("sd2", 3, -2e2, 2e2, 1.0);
    break;
  case 2:
    prms.add_prm("of1", 4, -1e4, 1e4, 1.0);
    prms.add_prm("sd2", 3, -2e2, 2e2, 1.0);
    break;
  case 3:
    prms.add_prm("of1", 4, -1e4, 1e4, 1.0);
    prms.add_prm("sd2", 3, -2e2, 2e2, 1.0);
    prms.add_prm("s",   3, -2e2, 2e2, 1.0);
    break;
  case 4:
    prms.add_prm("of1", 4, -1e4, 1e4, 1.0);
    prms.add_prm("sd2", 3, -2e2, 2e2, 1.0);
    prms.add_prm("s",   3, -2e2, 2e2, 1.0);
    prms.add_prm("sh1", 3, -2e2, 2e2, 1.0);
    prms.add_prm("sh2", 3, -2e2, 2e2, 1.0);
    break;
  }

  lat_constr.Fnum_b3.push_back(ElemIndex("sf1"));
  lat_constr.Fnum_b3.push_back(ElemIndex("sd1"));

  if (true) {
    no_sxt();
    set_bn_design_fam(ElemIndex("of1"), Oct, 0.0, 0.0);

    Fnum.push_back(ElemIndex("sf1"));
    Fnum.push_back(ElemIndex("sd1"));
    Fnum.push_back(ElemIndex("sd2"));
    fit_ksi1(Fnum, 0e0, 0e0, 1e1);
  }
}


void set_b3_constr_std(constr_type &constr)
{
  int              k, n;
  std::vector<int> mI_loc;

  mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 1));
  mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 3));
  lat_constr.mI_loc.push_back(mI_loc);

  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ss"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("sf1"), 1));

  n = lat_constr.high_ord_achr_Fnum.size()/2;
  lat_constr.high_ord_achr_dnu.resize(n);
  lat_constr.high_ord_achr_nu.resize(n);
  for (k = 0; k < n; k++) {
    lat_constr.high_ord_achr_dnu[k].resize(2, 0e0);
    lat_constr.high_ord_achr_nu[k].resize(2, 0e0);
  }

  for (k = 0; k < 2; k++)
    lat_constr.high_ord_achr_nu[0][k] = high_ord_achr_nu[0][k];

  for (k = 0; k < 2; k++)
    lat_constr.mI0[k] = mI_nu_ref[k];
}


void set_weights_std(constr_type &constr)
{
  lat_constr.eps_x_scl             = 1e7;

  lat_constr.ksi1_ctrl_scl[0]      = 0e-1;
  lat_constr.ksi1_ctrl_scl[1]      = 0e0;
  lat_constr.ksi1_ctrl_scl[2]      = 0e0;
  lat_constr.mI_scl[X_]            = 1e0*1e6;
  lat_constr.mI_scl[Y_]            = 1e0*1e6;
  lat_constr.high_ord_achr_scl[X_] = 1e-2*1e6;
  lat_constr.high_ord_achr_scl[Y_] = 1e-4*1e6;

  lat_constr.alpha_c_scl           = 1e-6;

  // Super Period.
  lat_constr.phi_scl               = ((dphi)? 1e0 : 0e0);
  lat_constr.phi0                  = 15.0;
  lat_constr.L_scl                 = 0e-10;
  lat_constr.L0                    = 10.0;

  lat_constr.drv_terms.get_h_scl(scl_ksi_1, scl_h_3, scl_h_3_delta, scl_h_4,
				 scl_ksi_2, scl_ksi_3, scl_chi_2,
				 scl_chi_delta_2);
}


void opt_mI_std(param_type &prms, constr_type &constr)
{
  // Parameter Type:
  //   Bend Angle  -3,
  //   Length      -2,
  //   Position    -1,
  //   Quadrupole   2.

  // Set parameters; initialized by optimizer.
  lat_prms.bn_tol = 1e-5; lat_prms.step = 1e0;

  set_dip_cell_std(prms, constr);
  set_b2_ms_std(prms);
  set_b3_Fam_std(prms);

  lat_constr.eps0_x = eps0_x;
  set_constr_std(constr);
  set_b3_constr_std(constr);
  set_weights_std(constr);

  lat_constr.ini_constr(true);

  prt_prms(lat_constr);
}


void set_dip_cell_sp(param_type &prms, constr_type &constr)
{
  std::vector<int>    grad_dip_Fnum;
  std::vector<double> grad_dip_scl;

  // Standard Cell.
  grad_dip_scl.push_back(0.129665);
  grad_dip_scl.push_back(0.149256);
  grad_dip_scl.push_back(0.174737);
  grad_dip_scl.push_back(0.216027);
  grad_dip_scl.push_back(0.330314);

  grad_dip_Fnum.push_back(ElemIndex("dl1a_1"));
  grad_dip_Fnum.push_back(ElemIndex("dl1a_2"));
  grad_dip_Fnum.push_back(ElemIndex("dl1a_3"));
  grad_dip_Fnum.push_back(ElemIndex("dl1a_4"));
  grad_dip_Fnum.push_back(ElemIndex("dl1a_5"));
  if (dphi)
    prms.add_prm(grad_dip_Fnum, grad_dip_scl, -3, -15.0, 15.0, 1.0);
  if (long_grad_dip[0])
    prms.add_prm(grad_dip_Fnum, grad_dip_scl,  2, -15.0, 15.0, 1.0);

  lat_constr.Fnum_b1.push_back(-ElemIndex("dl1a_1"));
  lat_constr.grad_dip_Fnum_b1.push_back(grad_dip_Fnum);

  if (true) {
    // Different longitudinal gradient dipoles.
    grad_dip_Fnum.clear();
    grad_dip_Fnum.push_back(ElemIndex("dl2a_1"));
    grad_dip_Fnum.push_back(ElemIndex("dl2a_2"));
    grad_dip_Fnum.push_back(ElemIndex("dl2a_3"));
    grad_dip_Fnum.push_back(ElemIndex("dl2a_4"));
    grad_dip_Fnum.push_back(ElemIndex("dl2a_5"));
    if (dphi)
      prms.add_prm(grad_dip_Fnum, grad_dip_scl, -3, -15.0, 15.0, 1.0);
    if (long_grad_dip[1])
      prms.add_prm(grad_dip_Fnum, grad_dip_scl,  2, -15.0, 15.0, 1.0);

    lat_constr.Fnum_b1.push_back(-ElemIndex("dl2a_1"));
    lat_constr.grad_dip_Fnum_b1.push_back(grad_dip_Fnum);
  }

  if (dphi) {
    // Dipole Cell.
    prms.add_prm("qf4", -3, -15.0, 15.0, 1.0);
    lat_constr.Fnum_b1.push_back(ElemIndex("qf4"));
    if (qf6_rb) {
      prms.add_prm("qf6", -3, -15.0, 15.0, 1.0);
      lat_constr.Fnum_b1.push_back(ElemIndex("qf6"));
    }
    prms.add_prm("qf8", -3, -15.0, 15.0, 1.0);
    lat_constr.Fnum_b1.push_back(ElemIndex("qf8"));

    // Commented out must be defined last.
    // prms.add_prm("dq1", -3, -15.0,   15.0,  1.0);
    lat_constr.Fnum_b1.push_back(ElemIndex("dq1"));
  }

  prms.add_prm("dq1", 2, -15.0, 15.0, 1.0);

  prms.add_prm("qd3", 2, -15.0,  0.0, 1.0);
  prms.add_prm("qf4", 2,   0.0, 15.0, 1.0);
  prms.add_prm("qd5", 2, -15.0,  0.0, 1.0);
}


void set_b2_ms_std_ls_sp(param_type &prms)
{
  // Mid Straight.
  prms.add_prm("qd2",  2, -15.0, 15.0, 1.0);
  prms.add_prm("qf1",  2, -15.0, 15.0, 1.0);

  // Standard Straight.
  prms.add_prm("qf6", 2, -15.0, 15.0, 1.0);
  prms.add_prm("qf8", 2, -15.0, 15.0, 1.0);

  // Long Straight.
  if (pert_dip_cell)
    prms.add_prm("qd3_c1",   2, -15.0, 15.0, 1.0);
  if (!true)
    prms.add_prm("qd2_c1",   2, -15.0,  0.0, 1.0);
  else
    prms.add_prm("qd2_c1",   2, -15.0, 15.0, 1.0);
  prms.add_prm("qf1_c1",   2, -15.0, 15.0, 1.0);
  prms.add_prm("quad_add", 2, -15.0, 15.0, 1.0);
}


void set_constr_sp(constr_type &constr)
{
  int k;

  // Lattice constraints are: alpha_x,y, beta_x,y, eta_x, eta'_x.
  // Include Matching Sections for Std & Long Straights.
  constr.add_constr(Elem_GetPos(ElemIndex("quad_add"), 1)-1,
		    0e0, 0e0, 0e0, 0e0, 1e8, 1e8,
		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  constr.add_constr(Elem_GetPos(ElemIndex("qf1"), 1),
		    0e0, 0e0, 0e0, 0e0, 1e8, 1e8,
		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  // Include constraint on alpha; in case of using ps_rot.
  for (k = 1; k <= 2; k++)
    constr.add_constr(Elem_GetPos(ElemIndex("ms"), k),
		      1e6, 1e6, 1e2,         1e-1*1e3,         1e8,      1e8,
		      0.0, 0.0, beta_ms[X_], beta_ms[Y_], eta_ms_x, 0.0);
  for (k = 1; k <= 2; k++)
    constr.add_constr(Elem_GetPos(ElemIndex("ss"), k),
		      1e6, 1e6, 1e2,         1e3,         1e8, 1e8,
		      0.0, 0.0, beta_ss[X_], beta_ss[Y_], 0.0, 0.0);
  constr.add_constr(Elem_GetPos(ElemIndex("ls"), 1),
		    1e6, 1e6, 1e2,         1e2,         1e8, 1e8,
		    0.0, 0.0, beta_ls[X_], beta_ls[Y_], 0.0, 0.0);
  // Symmetry & peak dispersion.
  constr.add_constr(Elem_GetPos(ElemIndex("sf1"), 1),
		    5e2, 5e2, 0e3,  0e2, 0e0,       0e0,
		    0.0, 0.0, 10.0, 1.0, 0.0, 0.0);

  if (false)
    // Increase beta_y.
    constr.add_constr(Elem_GetPos(ElemIndex("sd2"), 1),
		      0e5, 0e5, 1e1, 1e3, 0e7, 0e7,
		      0.0, 0.0, 1.0, 8.0, 0.0, 0.0);
  if (false)
    // Increase beta_y.
    constr.add_constr(Elem_GetPos(ElemIndex("qd5"), 1),
		      0e5, 0e5, 0e1, 1e3,  0e7, 0e7,
		      0.0, 0.0, 0.0, 10.0, 0.0, 0.0);
}


void set_b3_Fam_sp(param_type &prms)
{
  std::vector<int> Fnum;

  switch (3) {
  case 1:
    // First control ksi2_x,y.
    prms.add_prm("of1", 4, -1e4, 1e4, 1.0);

    lat_constr.Fnum_b3.push_back(ElemIndex("sf1"));
    lat_constr.Fnum_b3.push_back(ElemIndex("sd1"));
    lat_constr.Fnum_b3.push_back(ElemIndex("sd2"));
    break;
  case 2:
    // Control  [k_22000, k_11110, k_00220] & ksi2_x,y.
    prms.add_prm("sh1", 3, -2e2, 2e2, 1.0);
    prms.add_prm("sh2", 3, -2e2, 2e2, 1.0);
    prms.add_prm("s",   3, -2e2, 2e2, 1.0);
    prms.add_prm("sd2", 3, -2e2, 2e2, 1.0);
    prms.add_prm("of1", 4, -1e4, 1e4, 1.0);

    lat_constr.Fnum_b3.push_back(ElemIndex("sf1"));
    lat_constr.Fnum_b3.push_back(ElemIndex("sd1"));
    break;
  case 3:
    // Control  [k_22000, k_11110, k_00220] & ksi2_x,y.
    prms.add_prm("sf1", 4, -1e4, 1e4, 1.0);
    prms.add_prm("sd1", 4, -1e4, 1e4, 1.0);
    prms.add_prm("sd2", 4, -1e4, 1e4, 1.0);

    prms.add_prm("sf1", 4, -1e4, 1e4, 1.0);
    prms.add_prm("sd1", 4, -1e4, 1e4, 1.0);
    prms.add_prm("sd2", 4, -1e4, 1e4, 1.0);

    if (false) {
      prms.add_prm("sf1", 5, -1e7, 1e7, 1.0);
      prms.add_prm("sd1", 5, -1e7, 1e7, 1.0);
      prms.add_prm("sd2", 5, -1e7, 1e7, 1.0);
    }
    set_bn_design_fam(ElemIndex("sf1"), 5, 0.0, 0.0);
    set_bn_design_fam(ElemIndex("sd1"), 5, 0.0, 0.0);
    set_bn_design_fam(ElemIndex("sd2"), 5, 0.0, 0.0);

    if (true) {
      prms.add_prm("sh1", 3, -2e2, 2e2, 1.0);
      prms.add_prm("sh2", 3, -2e2, 2e2, 1.0);
      prms.add_prm("s",   3, -2e2, 2e2, 1.0);
      prms.add_prm("sd2", 3, -2e2, 2e2, 1.0);
    }

    lat_constr.Fnum_b3.push_back(ElemIndex("sf1"));
    lat_constr.Fnum_b3.push_back(ElemIndex("sd1"));
    break;
  }

  if (!true) {
    no_sxt();
    if (true) set_bn_design_fam(ElemIndex("of1"), Oct, 0.0, 0.0);

    Fnum.push_back(ElemIndex("sf1"));
    Fnum.push_back(ElemIndex("sd1"));
    Fnum.push_back(ElemIndex("sd2"));
    fit_ksi1(Fnum, 0e0, 0e0, 1e1);
  }
}


void set_b3_constr_sp(constr_type &constr)
{
  int              j, k, n;
  std::vector<int> mI_loc;

  lat_constr.eps0_x = eps0_x;

  mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 1));
  mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 3));
  lat_constr.mI_loc.push_back(mI_loc);
  mI_loc.clear();
  mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 5));
  mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 7));
  lat_constr.mI_loc.push_back(mI_loc);
  mI_loc.clear();
  mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 9));
  mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 11));
  lat_constr.mI_loc.push_back(mI_loc);
  mI_loc.clear();
  mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 13));
  mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 15));
  lat_constr.mI_loc.push_back(mI_loc);

  // Std straight.
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("sf1"), 3));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ss"), 1));
  // Long straight.
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ls"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("sf1"), 1));

  n = lat_constr.high_ord_achr_Fnum.size()/2;
  lat_constr.high_ord_achr_dnu.resize(n);
  lat_constr.high_ord_achr_nu.resize(n);
  for (k = 0; k < n; k++) {
    lat_constr.high_ord_achr_dnu[k].resize(2, 0e0);
    lat_constr.high_ord_achr_nu[k].resize(2, 0e0);
  }

  for (j = 0; j < n; j++)
    for (k = 0; k < 2; k++)
      lat_constr.high_ord_achr_nu[j][k] = high_ord_achr_nu[j][k];

  for (k = 0; k < 2; k++)
    lat_constr.mI0[k] = mI_nu_ref[k];
}


void set_weights_sp(constr_type &constr)
{
  lat_constr.eps_x_scl             = 1e0*1e7;

  lat_constr.ksi1_ctrl_scl[0]      = 0e-1;
  lat_constr.ksi1_ctrl_scl[1]      = 0e0;
  lat_constr.ksi1_ctrl_scl[2]      = 0e-1;
  lat_constr.mI_scl[X_]            = 1e-5*1e6;
  lat_constr.mI_scl[Y_]            = 1e-5*1e6;
  lat_constr.high_ord_achr_scl[X_] = 1e-5*1e6;
  lat_constr.high_ord_achr_scl[Y_] = 1e-5*1e6;

  lat_constr.alpha_c_scl           = 1e0*5e-7;

  lat_constr.phi_scl               = ((dphi)? 1e0 : 0e0);
  lat_constr.phi0                  = 60.0;
  lat_constr.L_scl                 = 0e-10;
  lat_constr.L0                    = 10.0;

  lat_constr.nu_ref_scl            = nu_ref_scl;
  lat_constr.nu_ref[X_]            = nu_ref[X_];
  lat_constr.nu_ref[Y_]            = nu_ref[Y_];

  lat_constr.drv_terms.get_h_scl(scl_ksi_1, scl_h_3, scl_h_3_delta, scl_h_4,
				 scl_ksi_2, scl_ksi_3, scl_chi_2,
				 scl_chi_delta_2);
}


void opt_mI_sp(param_type &prms, constr_type &constr)
{
  // Parameter Type:
  //   Bend Angle  -3,
  //   Length      -2,
  //   Position    -1,
  //   Quadrupole   2.

  // Set parameters; initialized by optimizer.
  lat_prms.bn_tol = 1e-5; lat_prms.step = 1.0;

  set_dip_cell_sp(prms, constr);
  set_b2_ms_std_ls_sp(prms);
  set_b3_Fam_sp(prms);

  lat_constr.eps0_x = eps0_x;
  set_constr_sp(constr);
  set_b3_constr_sp(constr);

  set_weights_sp(constr);

  lat_constr.ini_constr(true);

  prt_prms(lat_constr);
}


void set_b2_ls(param_type &prms)
{
  if (pert_dip_cell)
    prms.add_prm("qd3_c1", 2, -15.0, 15.0, 1.0);
  // Long Straight.
  prms.add_prm("qd2_c1",   2, -15.0,  0.0, 1.0);
  prms.add_prm("qf1_c1",   2, -15.0, 15.0, 1.0);
  prms.add_prm("quad_add", 2, -15.0, 15.0, 1.0);

  // Parameters are initialized in optimizer.
}


void set_constr_ls(constr_type &constr)
{
  if (!pert_dip_cell)
    // Eta_x & eta_x' are zero after Dipole Cell.
    constr.add_constr(Elem_GetPos(ElemIndex("ls"), 1),
		      1e4, 1e4, 1e0,         1e-1,        1e-10, 1e-10,
		      0.0, 0.0, beta_ls[X_], beta_ls[Y_], 0.0,   0.0);
  else {
    constr.add_constr(Elem_GetPos(ElemIndex("ls"), 1),
		      1e2, 1e2, 1e0,         1e-1,        1e-10, 1e-10,
		      0.0, 0.0, beta_ls[X_], beta_ls[Y_], 0.0,   0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("quad_add"), 1),
		      0e0, 0e0, 0e0, 0e0, 1e5, 1e5,
		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  }
}


void set_b3_constr_ls(constr_type &constr)
{
  int k, n;
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("sf1"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ls"), 1));

  n = lat_constr.high_ord_achr_Fnum.size() - 1;
  lat_constr.high_ord_achr_dnu.resize(n);
  lat_constr.high_ord_achr_nu.resize(n);
  for (k = 0; k < n; k++) {
    lat_constr.high_ord_achr_dnu[k].resize(2, 0e0);
    lat_constr.high_ord_achr_nu[k].resize(2, 0e0);
  }

  for (k = 0; k < 2; k++)
    lat_constr.high_ord_achr_nu[0][k] = high_ord_achr_nu[0][k];

}


void set_weights_ls(constr_type &constr)
{
  lat_constr.high_ord_achr_scl[X_] = 1e-6*1e6;
  lat_constr.high_ord_achr_scl[Y_] = 1e-6*1e6;
}


void match_ls(param_type &prms, constr_type &constr)
{
  // Parameter Type:
  //   Bend Angle  -3,
  //   Length      -2,
  //   Position    -1,
  //   Quadrupole   2.
  int j, k;

  // Perturbed symmetry at end of Dipole Cell: 
  //   1. Initialize with: [Qf1, Qd2, Qd3].
  //   2. Exclude for 1st pass: pert_dip_cell = false.
  //   3. Include for fine tuning.
  lat_prms.bn_tol = 1e-5; lat_prms.step = 1.0;
  set_b2_ls(prms);

  set_constr_ls(constr);
  set_b3_constr_ls(constr);
  set_weights_ls(constr);

  lat_constr.ini_constr(false);

  for (j = 0; j < n_ic; j++)
    for (k = 0; k < 2; k++)
      lat_constr.ic[j][k] = ic[j][k];
}


void set_b2_ss(param_type &prms)
{
  // Std Straight.
  prms.add_prm("qd2",  2, -15.0, 15.0, 1.0);
  prms.add_prm("qf1",  2, -15.0, 15.0, 1.0);
  // prms.add_prm("qf1", -1,   0.0, -0.5, 1.0);

  // Parameters are initialized in optimizer.
}


void set_constr_ss(constr_type &constr)
{
  constr.add_constr(Elem_GetPos(ElemIndex("ss"), 1),
		    1e2, 1e2, 1e-1,        1e-1,         1e-10, 1e-10,
		    0.0, 0.0, beta_ss[X_], beta_ss[Y_],  0.0,   0.0);
  constr.add_constr(Elem_GetPos(ElemIndex("qf1"), 1),
		    0e0, 0e0, 0e0, 0e0, 1e6, 1e6,
		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
}


void set_b3_constr_ss(constr_type &constr)
{
  int k, n;

  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ms"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ss"), 1));

  n = lat_constr.high_ord_achr_Fnum.size()/2;
  lat_constr.high_ord_achr_dnu.resize(n);
  lat_constr.high_ord_achr_nu.resize(n);
  for (k = 0; k < n; k++) {
    lat_constr.high_ord_achr_dnu[k].resize(2, 0e0);
    lat_constr.high_ord_achr_nu[k].resize(2, 0e0);
  }

  for (k = 0; k < 2; k++)
    lat_constr.high_ord_achr_nu[0][k] = high_ord_achr_nu[0][k]/2e0;
}


void set_weights_ss(constr_type &constr)
{
  lat_constr.high_ord_achr_scl[X_] = 0e0;
  lat_constr.high_ord_achr_scl[Y_] = 0e0;
}


void match_ss(param_type &prms, constr_type &constr)
{
  // Parameter Type:
  //   Bend Angle  -3,
  //   Length      -2,
  //   Position    -1,
  //   Quadrupole   2.
  int j, k;

  // Perturbed symmetry at end of Dipole Cell: 
  //   1. Initialize with: [Qf1, Qd2, Qd3].
  //   2. Exclude for 1st pass: pert_dip_cell = false.
  //   3. Include for fine tuning.
  lat_prms.bn_tol = 1e-5; lat_prms.step = 1.0;
  set_b2_ss(prms);

  set_constr_ss(constr);
  set_b3_constr_ss(constr);
  set_weights_ss(constr);

  lat_constr.ini_constr(false);

  for (j = 0; j < n_ic; j++)
    for (k = 0; k < 2; k++)
      lat_constr.ic[j][k] = ic[j][k];
}


void set_b3_Fam_mult(param_type &prms)
{
  std::vector<int> Fnum;

  switch (3) {
  case 1:
    // Then control [k_22000, k_11110, k_00220].
    prms.add_prm("sh1", 3, -2e2, 2e2, 1.0);
    prms.add_prm("sh2", 3, -2e2, 2e2, 1.0);
    prms.add_prm("s",   3, -2e2, 2e2, 1.0);
    prms.add_prm("sd2", 3, -2e2, 2e2, 1.0);
    break;
  case 2:
    prms.add_prm("of1", 4, -1e4, 1e4, 1.0);
    prms.add_prm("sh1", 3, -2e2, 2e2, 1.0);
    prms.add_prm("sh2", 3, -2e2, 2e2, 1.0);
    prms.add_prm("s",   3, -2e2, 2e2, 1.0);
    prms.add_prm("sd2", 3, -2e2, 2e2, 1.0);
    break;
  case 3:
    // Control  [k_22000, k_11110, k_00220] & ksi2_x,y.
    prms.add_prm("sf1", 4, -1e4, 1e4, 1.0);
    prms.add_prm("sd1", 4, -1e4, 1e4, 1.0);
    prms.add_prm("sd2", 4, -1e4, 1e4, 1.0);

    prms.add_prm("sf1", 4, -1e4, 1e4, 1.0);
    prms.add_prm("sd1", 4, -1e4, 1e4, 1.0);
    prms.add_prm("sd2", 4, -1e4, 1e4, 1.0);

    if (false) {
      prms.add_prm("sf1", 5, -1e7, 1e7, 1.0);
      prms.add_prm("sd1", 5, -1e7, 1e7, 1.0);
      prms.add_prm("sd2", 5, -1e7, 1e7, 1.0);
    }

    if (true) {
      prms.add_prm("sh1", 3, -2e2, 2e2, 1.0);
      prms.add_prm("sh2", 3, -2e2, 2e2, 1.0);
      prms.add_prm("s",   3, -2e2, 2e2, 1.0);
      prms.add_prm("sd2", 3, -2e2, 2e2, 1.0);
    }

    lat_constr.Fnum_b3.push_back(ElemIndex("sf1"));
    lat_constr.Fnum_b3.push_back(ElemIndex("sd1"));
    break;
  }

  lat_constr.Fnum_b3.push_back(ElemIndex("sf1"));
  lat_constr.Fnum_b3.push_back(ElemIndex("sd1"));

  if (!true) {
    no_sxt();

    set_bn_design_fam(ElemIndex("of1"), Oct, 0.0, 0.0);

    Fnum.push_back(ElemIndex("sf1"));
    Fnum.push_back(ElemIndex("sd1"));
    Fnum.push_back(ElemIndex("sd2"));
    fit_ksi1(Fnum, 0e0, 0e0, 1e1);
  }
}


void set_weights_mult(constr_type &constr)
{
  lat_constr.eps_x_scl             = 0e0;

  lat_constr.ksi1_ctrl_scl[0]      = 0e0;
  lat_constr.ksi1_ctrl_scl[1]      = 0e0;
  lat_constr.ksi1_ctrl_scl[2]      = 0e0;
  lat_constr.mI_scl[X_]            = 0e0;
  lat_constr.mI_scl[Y_]            = 0e0;
  lat_constr.high_ord_achr_scl[X_] = 0e0;
  lat_constr.high_ord_achr_scl[Y_] = 0e0;

  lat_constr.drv_terms.get_h_scl(scl_ksi_1, scl_h_3, scl_h_3_delta, scl_h_4,
				 scl_ksi_2, scl_ksi_3, scl_chi_2,
				 scl_chi_delta_2);
}


void opt_mult(param_type &prms, constr_type &constr)
{
  // Set parameters; initialized by optimizer.
  lat_prms.bn_tol = 1e-5; lat_prms.step = 1.0;

  set_b3_Fam_mult(prms);
  set_weights_mult(constr);

  lat_constr.ini_constr(true);
  prt_prms(lat_constr);

  get_nu(lat_constr.nu);
  eps_x = get_lin_opt(lat_constr);
}


void match_als_u(param_type &prms, constr_type &constr)
{
  // Parameter Type:
  //   Bend Angle  -3,
  //   Length      -2,
  //   Position    -1,
  //   Quadrupole   2.
  int                 j, k, n;
  std::vector<int>    grad_dip_Fnum;
  std::vector<double> grad_dip_scl;

  // From Center of Mid Straight: alpha, beta, eta, eta'.
  const int    n_ic        = 4;
  const double ic[n_ic][2] =
    {{0.0, 0.0}, {0.1153496843, 3.2620137331}, {0.0002943186, 0.0}, {0.0, 0.0}};
 
  prms.add_prm("q1",   2, -20.0, 20.0, 1.0);
  prms.add_prm("q2",   2, -20.0, 20.0, 1.0);
  prms.add_prm("q3",   2, -20.0, 20.0, 1.0);
  prms.add_prm("q4",   2, -20.0, 20.0, 1.0);

  prms.add_prm("mqfh", 2, -20.0, 20.0, 1.0);

  // Parameters are initialized in optimizer.

  constr.add_constr(Elem_GetPos(ElemIndex("q1"), 1),
  		    0e0, 0e0, 0e0, 0e0, 1e3, 1e3,
  		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  constr.add_constr(Elem_GetPos(ElemIndex("ssh"), 1),
  		    1e0, 1e0, 1e-5, 1e-5, 1e5, 1e5,
  		    0.0, 0.0, 2.1,  2.3,  0.0, 0.0);

  lat_prms.bn_tol = 1e-5; lat_prms.step = 1.0;

  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ssh"), 1));

  n = lat_constr.high_ord_achr_Fnum.size() - 1;
  lat_constr.high_ord_achr_dnu.resize(n);
  for (k = 0; k < n; k++)
    lat_constr.high_ord_achr_dnu[k].resize(2, 0e0);

  for (k = 0; k < 2; k++)
    lat_constr.high_ord_achr_nu[0][k] = high_ord_achr_nu[0][k];

  lat_constr.high_ord_achr_scl[X_] = 0e0;
  lat_constr.high_ord_achr_scl[Y_] = 0e0;

  lat_constr.ini_constr(false);

  for (j = 0; j < n_ic; j++)
    for (k = 0; k < 2; k++)
      lat_constr.ic[j][k] = ic[j][k];
}


void opt_tba(param_type &prms, constr_type &constr)
{
  // Parameter Type:
  //   Bend Angle  -3,
  //   Length      -2,
  //   Position    -1,
  //   Quadrupole   2.
  int k, n;

  // TBA Cell.
  prms.add_prm("b1",       -3, -20.0,   20.0,  1.0);
  prms.add_prm("b1",       -2,   0.1,    1.2,  1.0);
  // prms.add_prm("b1",       -1,   0.075,  0.07, 1.0);
  prms.add_prm("b1",        2, -20.0,   20.0,  1.0);
  prms.add_prm("b2",       -2,   0.1,    0.7,  1.0);
  prms.add_prm("b2",        2, -20.0,   20.0,  1.0);
  prms.add_prm("qf1a_tba",  2, -20.0,   20.0,  1.0);
  prms.add_prm("qf1a_tba", -2,   0.1,    0.25, 1.0);
  prms.add_prm("qf1b_tba",  2, -20.0,   20.0,  1.0);
  prms.add_prm("qf1b_tba", -2,   0.1,    0.25, 1.0);

  // Mid-Straight.
  prms.add_prm("d8",       -2,   0.075,  0.2,  1.0);
  prms.add_prm("qf1_ms",    2, -20.0,   20.0,  1.0);
  prms.add_prm("qf1_ms",   -2,   0.1,    0.25, 1.0);
  prms.add_prm("d9",       -2,   0.075,  0.2,  1.0);
  prms.add_prm("qd2_ms",    2, -20.0,   20.0,  1.0);
  prms.add_prm("qd2_ms",   -2,   0.1,    0.2,  1.0);

 // Parameters are initialized in optimizer.

  // Lattice constraints are: alpha_x,y, beta_x,y, eta_x, eta'_x.
  constr.add_constr(Elem_GetPos(ElemIndex("b1"), 2),
  		    0e0, 0e0, 0e0, 0e0, 1e5, 1e5,
  		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

  // Include constraint on alpha; in case of using ps_rot.
  constr.add_constr(Elem_GetPos(ElemIndex("ms"), 1),
  		    1e4, 1e4, 1e-10, 1e-10, 0e0, 0e0,
  		    0.0, 0.0, 3.0,   1.5,   0.0, 0.0);

  lat_prms.bn_tol = 1e-6; lat_prms.step = 1.0;

  lat_constr.Fnum_b3.push_back(ElemIndex("sd1a"));
  lat_constr.Fnum_b3.push_back(ElemIndex("sfh"));
  lat_constr.Fnum_b3.push_back(ElemIndex("sd1b"));

  lat_constr.Fnum_b1.push_back(ElemIndex("b1"));
  lat_constr.Fnum_b1.push_back(ElemIndex("b2"));

  lat_constr.high_ord_achr_scl[X_] = 1e3;
  lat_constr.high_ord_achr_scl[Y_] = 1e3;
  for (k = 0; k < 2; k++)
    lat_constr.high_ord_achr_nu[0][k] = high_ord_achr_nu[0][k];

  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ms"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ms"), 2));

  n = lat_constr.high_ord_achr_Fnum.size() - 1;
  lat_constr.high_ord_achr_dnu.resize(n);
  for (k = 0; k < n; k++)
    lat_constr.high_ord_achr_dnu[k].resize(2, 0e0);

  lat_constr.eps_x_scl = 1e2; lat_constr.eps0_x = 0.190;

  lat_constr.mI_scl[X_] = lat_constr.mI_scl[Y_] = 1e-10;
  for (k = 0; k < 2; k++)
    lat_constr.mI0[k] = mI_nu_ref[k];

  // 2 TBA & 1 MS: phi = 7.5.
  lat_constr.phi_scl = 1e0;
  lat_constr.phi0    = 7.5;
  lat_constr.L_scl   = 1e-10;
  lat_constr.L0      = 10.0;

  lat_constr.ini_constr(true);
}


void match_ms(param_type &prms, constr_type &constr)
{
  // Parameter Type:
  //   Bend Angle  -3,
  //   Length      -2,
  //   Position    -1,
  //   Quadrupole   2.
  int                 j, k;
  std::vector<int>    grad_dip_Fnum;
  std::vector<double> grad_dip_scl;

  // From Center of Mid Straight: alpha, beta, eta, eta'.
  const int    n_ic        = 4;
  const double ic[n_ic][2] =
    {{0.0, 0.0}, {7.08276, 3.03915}, {0.0, 0.0}, {0.0, 0.0}};
 
  // Standard Straight.
  grad_dip_scl.push_back(0.129665);
  grad_dip_scl.push_back(0.149256);
  grad_dip_scl.push_back(0.174737);
  grad_dip_scl.push_back(0.216027);
  grad_dip_scl.push_back(0.330314);

  grad_dip_Fnum.push_back(ElemIndex("bl1_1"));
  grad_dip_Fnum.push_back(ElemIndex("bl1_2"));
  grad_dip_Fnum.push_back(ElemIndex("bl1_3"));
  grad_dip_Fnum.push_back(ElemIndex("bl1_4"));
  grad_dip_Fnum.push_back(ElemIndex("bl1_5"));
  prms.add_prm(grad_dip_Fnum, grad_dip_scl, 2, -20.0, 20.0, 1.0);

  grad_dip_Fnum.clear();
  grad_dip_Fnum.push_back(ElemIndex("bl2_1"));
  grad_dip_Fnum.push_back(ElemIndex("bl2_2"));
  grad_dip_Fnum.push_back(ElemIndex("bl2_3"));
  grad_dip_Fnum.push_back(ElemIndex("bl2_4"));
  grad_dip_Fnum.push_back(ElemIndex("bl2_5"));
  prms.add_prm(grad_dip_Fnum, grad_dip_scl, 2, -20.0, 20.0, 1.0);

  // Parameters are initialized in optimizer.

  constr.add_constr(Elem_GetPos(ElemIndex("ms"), 1),
  		    1e0, 1e0, 1e-5, 1e-5, 1e0, 1e0,
  		    0.0, 0.0, 3.0,  1.5,  0.0, 0.0);

  lat_prms.bn_tol = 1e-6; lat_prms.step = 1.0;

  lat_constr.high_ord_achr_scl[X_] = 0e0;
  lat_constr.high_ord_achr_scl[Y_] = 0e0;
  for (k = 0; k < 2; k++)
    lat_constr.high_ord_achr_nu[0][k] = high_ord_achr_nu[0][k]/2.0;

  lat_constr.ini_constr(false);

  for (j = 0; j < n_ic; j++)
    for (k = 0; k < 2; k++)
      lat_constr.ic[j][k] = ic[j][k];
}


void opt_std_cell(param_type &prms, constr_type &constr)
{
  // Parameter Type:
  //   Bend Angle  -3,
  //   Length      -2,
  //   Position    -1,
  //   Quadrupole   2.
  int k, n;

  // TBA Cell.
  prms.add_prm("b1",       -3, -20.0,   20.0,  1.0);
  prms.add_prm("b1",       -2,   0.1,    1.2,  1.0);
  // prms.add_prm("b1",       -1,   0.075,  0.07, 1.0);
  prms.add_prm("b1",        2, -20.0,   20.0,  1.0);
  prms.add_prm("b2",       -2,   0.1,    0.7,  1.0);
  prms.add_prm("b2",        2, -20.0,   20.0,  1.0);
  prms.add_prm("qf1a_tba",  2, -20.0,   20.0,  1.0);
  prms.add_prm("qf1a_tba", -2,   0.1,    0.25, 1.0);
  prms.add_prm("qf1b_tba",  2, -20.0,   20.0,  1.0);
  prms.add_prm("qf1b_tba", -2,   0.1,    0.25, 1.0);

  // Mid-Straight.
  prms.add_prm("d8",       -2,   0.075,  0.2,  1.0);
  prms.add_prm("qf1_ms",    2, -20.0,   20.0,  1.0);
  prms.add_prm("qf1_ms",   -2,   0.1,    0.25, 1.0);
  prms.add_prm("d9",       -2,   0.075,  0.2,  1.0);
  prms.add_prm("qd2_ms",    2, -20.0,   20.0,  1.0);
  prms.add_prm("qd2_ms",   -2,   0.1,    0.2,  1.0);

  // Standard Straight.
  prms.add_prm("d10",    -2,   0.075,  0.4, 1.0);
  prms.add_prm("qd3_ss",  2, -20.0,   20.0, 1.0);
  prms.add_prm("qd3_ss", -2,   0.1,    0.4, 1.0);
  prms.add_prm("d11",    -2,   0.075,  0.4, 1.0);
  prms.add_prm("qf2_ss",  2, -20.0,   20.0, 1.0);
  prms.add_prm("qf2_ss", -2,   0.15,   0.4, 1.0);
  prms.add_prm("d12",    -2,   0.075,  0.4, 1.0);
  prms.add_prm("qd1_ss",  2, -20.0,   20.0, 1.0);
  prms.add_prm("qd1_ss", -2,   0.1,    0.4, 1.0);

  // Parameters are initialized in optimizer.

  // Lattice constraints are: alpha_x,y, beta_x,y, eta_x, eta'_x.
  constr.add_constr(Elem_GetPos(ElemIndex("b2"), 1),
  		    1e3, 1e3, 0e0, 0e-1, 0e0, 1e6,
  		    0.0, 0.0, 0.0, 10.0, 0.0, 0.0);
  constr.add_constr(Elem_GetPos(ElemIndex("b1"), 1)-1,
  		    0e0, 0e0, 0e0, 0e0, 1e6, 1e6,
  		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  constr.add_constr(Elem_GetPos(ElemIndex("b1"), 2),
  		    0e0, 0e0, 0e0, 0e0, 1e6, 1e6,
  		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

  constr.add_constr(Elem_GetPos(ElemIndex("ms"), 1),
  		    1e3, 1e3, 1e-2, 1e-2, 0e0, 0e0,
  		    0.0, 0.0, 3.0,  1.5,  0.0, 0.0);

  constr.add_constr(Elem_GetPos(ElemIndex("ss"), 1),
  		    1e3, 1e3, 1e-2, 1e-2, 0e0, 0e0,
  		    0.0, 0.0, 4.0,  2.5,  0.0, 0.0);

 lat_prms.bn_tol = 1e-6; lat_prms.step = 1.0;

  lat_constr.Fnum_b3.push_back(ElemIndex("sd1a"));
  lat_constr.Fnum_b3.push_back(ElemIndex("sfh"));
  lat_constr.Fnum_b3.push_back(ElemIndex("sd1b"));

  lat_constr.Fnum_b1.push_back(ElemIndex("b1"));
  lat_constr.Fnum_b1.push_back(ElemIndex("b2"));

  lat_constr.high_ord_achr_scl[X_] = 1e3;
  lat_constr.high_ord_achr_scl[Y_] = 1e3;
  for (k = 0; k < 2; k++)
    lat_constr.high_ord_achr_nu[0][k] = high_ord_achr_nu[0][k]/2.0;

  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ms"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("b2"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ss"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("b2"), 3));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("ms"), 2));

  n = lat_constr.high_ord_achr_Fnum.size() - 1;
  lat_constr.high_ord_achr_dnu.resize(n);
  for (k = 0; k < n; k++)
    lat_constr.high_ord_achr_dnu[k].resize(2, 0e0);

  lat_constr.mI_scl[X_] = lat_constr.mI_scl[Y_] = 1e-10;
  for (k = 0; k < 2; k++)
    lat_constr.mI0[k] = mI_nu_ref[k];

  lat_constr.eps_x_scl = 1e3; lat_constr.eps0_x = 0.190;

  // 2 TBA: phi = 15.
  lat_constr.phi_scl = 1e0;
  lat_constr.phi0    = 15.0;
  lat_constr.L_scl   = 1e-10;
  lat_constr.L0      = 10.0;

  lat_constr.ini_constr(true);
}


int main(int argc, char *argv[])
{
  char             buffer[BUFSIZ];
  double           eps_x;

  const double dnu0[] = {0.0, 0.0};

  reverse_elem = !false;

  trace = false;

  // Unbuffered output.
  setvbuf(stdout, buffer, _IONBF, BUFSIZ);

  globval.mat_meth = !false;

  if (!true)
    Read_Lattice(argv[1]);
  else {
    rdmfile(argv[1]);
    globval.dPcommon = 1e-6;
  }

  globval.H_exact    = false; globval.quad_fringe    = false;
  globval.Cavity_on  = false; globval.radiation      = false;
  globval.emittance  = false; globval.IBS            = false;
  globval.pathlength = false; globval.bpm            = 0;
  globval.Cart_Bend  = false; globval.dip_edge_fudge = true;

  if (ps_rot) {
    Ring_GetTwiss(true, 0e0); printglob();
    set_map(ElemIndex("ps_rot"), dnu0);
    if (true) {
      set_map(ElemIndex("ps_rot"), dnu);
      Ring_GetTwiss(true, 0e0); printglob();
    }
  }

  if (false) {
    if (true) {
      Ring_GetTwiss(true, 0e0); printglob();
    } else
      ttwiss(lat_constr.ic[0], lat_constr.ic[1],
	     lat_constr.ic[2], lat_constr.ic[3], 0e0);

    eps_x = get_eps_x1(true, eps_x, sigma_E);
    printf("\neps_x = %5.3f nm.rad\n\n", eps_x);

    prt_lat("linlat1.out", globval.bpm, true);
    prt_lat("linlat.out", globval.bpm, true, 10);
    prt_chrom_lat();
  }

  switch (opt_case) {
  case 1:
    // Optimize Standard Straight: mI & 2*nu.
    opt_mI_std(lat_prms, lat_constr);
    fit_powell(lat_prms, 1e-3, f_achrom);
    // fit_sim_anneal(lat_prms, 1e-3, f_achrom);
    break;
  case 2:
    // Match Long Straight: mI.
    match_ls(lat_prms, lat_constr);
    fit_powell(lat_prms, 1e-3, f_match);
    break;
  case 3:
    opt_mI_sp(lat_prms, lat_constr);
    fit_powell(lat_prms, 1e-3, f_achrom);
    // fit_sim_anneal(lat_prms, 1e-3, f_achrom);
    break;
  case 4:
    // Match Short Straight: mI.
    match_ss(lat_prms, lat_constr);
    fit_powell(lat_prms, 1e-3, f_match);
    break;
  case 5:
    // Optimize multipoles.
    opt_mult(lat_prms, lat_constr);
    fit_powell(lat_prms, 1e-3, f_mult);
    break;
  default:
    printf("\nmain: unknown opt_case %d\n", opt_case);
    exit(1);
    break;
  }

  if (false) {
    // Match ALS-U.
    match_als_u(lat_prms, lat_constr);
    fit_powell(lat_prms, 1e-3, f_match);
  }

  if (false) {
    // Optimize TBA & Mid Straight: Higher-Order-Achromat.
    opt_tba(lat_prms, lat_constr);
    no_sxt();
    fit_powell(lat_prms, 1e-3, f_achrom);
  }

  if (false) {
    // Match Mid Straight.
    match_ms(lat_prms, lat_constr);
    no_sxt();
    fit_powell(lat_prms, 1e-3, f_match);
  }
}
