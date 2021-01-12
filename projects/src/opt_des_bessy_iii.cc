#define NO 1

#include "tracy_lib.h"

#include "Powell/src/newuoa.h"

#include "drv_terms.cc"

int no_tps = NO;


const bool
  ps_rot          = false, // Note, needs to be zeroed; after use.
  pert_dip_cell   = false,
  dphi            = true,
  long_grad_dip[] = {false, false},
  ksi_terms[]     = {!false, false, false},
  drv_terms[]     = {false, false, false},
  tune_fp_terms[] = {false, false};

/* opt_case:
     opt_unit_cell 1,
     match_disp    2,
     match_ss      3,
     opt_mI        4.                                                         */
const int opt_case = 2;

const int
  n_ic   = 5,
  n_cell = 6;

const double
  ic[n_ic][2] =
  {{-0.0000000003, 0.0000000002},
   {6.6214019255, 0.1814730089},
   {0.0380192360, 0.0000000000},
   {-0.0, 0.0}},                   // From Center of Mid Straight:
                                   // alpha, beta, eta, eta'.
 
  eps0_x                = 0.100,
  dnu[]                 = {0.1, 0.0},

  beta_inj[]            = {1.6, 1.6},
  A_max[]               = {6e-3, 3e-3},
  delta_max             = 2e-2,
  A_delta_max[]         = {2e-3, 0.1e-3},

  nu[]                  = {65.90/n_cell+dnu[X_], 20.73/n_cell+dnu[Y_]},
  nu_ref[]              = {nu[X_]-(int)nu[X_], nu[Y_]-(int)nu[Y_]},
  beta_ss[]             = {1.5, 1.5},

  scl_eps_x             = 5e7,
  nu_ref_scl            = 0*5e7,
  alpha_c_scl           = 1e0*5e-7,

  mI_scl                = 1e-6,
  high_ord_achr_scl[]   = {0e-6, 0e-6},

  scl_ksi_1             = (ksi_terms[0])?     1e0*5e-1 : 0e0,
  scl_h_3               = (drv_terms[0])?     1e18 : 0e0,
  scl_h_3_delta         = (drv_terms[1])?     1e-3*1e18 : 0e0,
  scl_h_4               = (drv_terms[2])?     1e18 : 0e0,
  scl_ksi_2             = (ksi_terms[1])?     1e0*1e7  : 0e0,
  scl_ksi_3             = (ksi_terms[2])?     1e0*1e7  : 0e0,
  scl_chi_2             = (tune_fp_terms[0])? 1e0*1e7  : 0e0,
  scl_chi_delta_2       = (tune_fp_terms[1])? 1e5  : 0e0,

  scl_extra             = 0e2,

  twoJ[]                =
    {sqr(A_max[X_])/beta_inj[X_], sqr(A_max[Y_])/beta_inj[Y_]},
  twoJ_delta[]          =
    {sqr(A_delta_max[X_])/beta_inj[X_], sqr(A_delta_max[Y_])/beta_inj[Y_]},
  mI_nu_ref[]           = {1.5, 0.5},
  high_ord_achr_nu[][2] = {{5.0/8.0, 2.0/8.0}, {4.0/8.0, 1.0/8.0}},

  phi_max               = 5.0,
  phi_rb_max            = 0.275,
  b_2_max               = 15.0,
  b_2_dip_max           = 0.7,
  b_3_chrom_max         = 2e3,
  b_3_max               = 2e3,
  b_4_max               = 1e4;


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


void set_unit_cell(param_type &prms, constr_type &constr)
{
  std::vector<int>    grad_dip_Fnum, mI_loc;
  std::vector<double> grad_dip_scl;

  if (false) {
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
      prms.add_prm(grad_dip_Fnum, grad_dip_scl, -3, -phi_max, phi_max, 1.0);
    if (long_grad_dip[0])
      prms.add_prm(grad_dip_Fnum, grad_dip_scl,  2, -b_2_dip_max, b_2_dip_max,
		   1.0);

    lat_constr.Fnum_b1.push_back(-ElemIndex("dl1a_1"));
    lat_constr.grad_dip_Fnum_b1.push_back(grad_dip_Fnum);

    grad_dip_Fnum.clear();
    grad_dip_Fnum.push_back(ElemIndex("dl2a_1"));
    grad_dip_Fnum.push_back(ElemIndex("dl2a_2"));
    grad_dip_Fnum.push_back(ElemIndex("dl2a_3"));
    grad_dip_Fnum.push_back(ElemIndex("dl2a_4"));
    grad_dip_Fnum.push_back(ElemIndex("dl2a_5"));
    if (dphi)
      prms.add_prm(grad_dip_Fnum, grad_dip_scl, -3, -phi_max, phi_max,
		   1.0);
    if (long_grad_dip[1])
      prms.add_prm(grad_dip_Fnum, grad_dip_scl,  2, -b_2_dip_max, b_2_dip_max,
		   1.0);

    lat_constr.Fnum_b1.push_back(-ElemIndex("dl2a_1"));
    lat_constr.grad_dip_Fnum_b1.push_back(grad_dip_Fnum);
  }

  if (dphi) {
    // Dipole Cell.
    prms.add_prm("b1", -3, -phi_max, phi_max, 1.0);
    lat_constr.Fnum_b1.push_back(ElemIndex("b1"));

    // Commented out must be defined last.
    // prms.add_prm("b2", -3, -5.0, 5.0,  1.0);
    lat_constr.Fnum_b1.push_back(ElemIndex("b2"));
  }

  prms.add_prm("b1", 2, -b_2_max, b_2_max, 1.0);
  prms.add_prm("b2", 2, -b_2_max, b_2_max, 1.0);
}


void set_unit_cell_b2(param_type &prms)
{
  // prms.add_prm("qf1",  2,      0.0, b_2_max, 1.0);
  // prms.add_prm("qd2",  2, -b_2_max,     0.0, 1.0);
}


void set_unit_cell_constr(constr_type &constr)
{
  // constr.add_constr(Elem_GetPos(ElemIndex("qf1"), 1)-1,
  // 		    0e0, 0e0, 0e0, 0e0, 1e8, 1e7,
  // 		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  // constr.add_constr(Elem_GetPos(ElemIndex("qf1"), 2),
  // 		    0e0, 0e0, 0e0, 0e0, 1e8, 1e7,
  // 		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  // constr.add_constr(Elem_GetPos(ElemIndex("dl1a_5"), 1)-1,
  // 		    0e0, 0e0, 0e0, 0e0, 1e8, 1e7,
  // 		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  // constr.add_constr(Elem_GetPos(ElemIndex("dl1a_5"), 2),
  // 		    0e0, 0e0, 0e0, 0e0, 1e8, 1e7,
  // 		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  // constr.add_constr(Elem_GetPos(ElemIndex("b0"), 2),
  // 		    1e5, 1e5, 1e2,         1e2,         1e7, 1e7,
  // 		    0.0, 0.0, beta_ss[X_], beta_ss[Y_], 0.0, 0.0);
}


void set_unit_cell_b3_Fam(param_type &prms)
{
  std::vector<int> Fnum;

  prms.add_prm("s1",3, -b_3_chrom_max, b_3_chrom_max, 1.0);
  prms.add_prm("s2",3, -b_3_chrom_max, b_3_chrom_max, 1.0);

  // Control of linear chromaticity.
  lat_constr.Fnum_b3.push_back(ElemIndex("s1"));
  lat_constr.Fnum_b3.push_back(ElemIndex("s2"));

  if (true) {
    no_sxt();
    Fnum.push_back(ElemIndex("s1"));
    Fnum.push_back(ElemIndex("s2"));
    fit_ksi1(Fnum, 0e0, 0e0, 1e1);
  }
}


void set_unit_cell_b3_constr(constr_type &constr)
{
  int              k, n;
  std::vector<int> mI_loc;

  // mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 1));
  // mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 3));
  // lat_constr.mI_loc.push_back(mI_loc);

  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("s0"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("s1"), 1));

  n = lat_constr.high_ord_achr_Fnum.size()/2;
  lat_constr.high_ord_achr_dnu.resize(n);
  lat_constr.high_ord_achr_nu.resize(n);
  for (k = 0; k < n; k++) {
    lat_constr.high_ord_achr_dnu[k].resize(2, 0e0);
    lat_constr.high_ord_achr_nu[k].resize(2, 0e0);
  }

  for (k = 0; k < 2; k++)
    lat_constr.high_ord_achr_nu[0][k] = high_ord_achr_nu[0][k];

  // for (k = 0; k < 2; k++)
  //   lat_constr.mI0[k] = mI_nu_ref[k];
}


void set_unit_cell_weights(constr_type &constr)
{
  lat_constr.eps_x_scl             = 1e7;

  lat_constr.ksi1_ctrl_scl[0]      = 0e-1;
  lat_constr.ksi1_ctrl_scl[1]      = 0e0;
  lat_constr.ksi1_ctrl_scl[2]      = 0e0;
  lat_constr.mI_scl[X_]            = 0e6;
  lat_constr.mI_scl[Y_]            = 0e6;
  lat_constr.high_ord_achr_scl[X_] = high_ord_achr_scl[X_];
  lat_constr.high_ord_achr_scl[Y_] = high_ord_achr_scl[Y_];

  lat_constr.alpha_c_scl           = 1e-6;

  // Super Period.
  lat_constr.phi_scl               = ((dphi)? 1e0 : 0e0);
  lat_constr.phi0                  = 4.5;
  lat_constr.L_scl                 = 0e-10;
  lat_constr.L0                    = 10.0;

  lat_constr.drv_terms.get_h_scl(scl_ksi_1, scl_h_3, scl_h_3_delta, scl_h_4,
				 scl_ksi_2, scl_ksi_3, scl_chi_2,
				 scl_chi_delta_2);
}


void opt_unit_cell(param_type &prms, constr_type &constr)
{
  // Parameter Type:
  //   Bend Angle  -3,
  //   Length      -2,
  //   Position    -1,
  //   Quadrupole   2.

  // Set parameters; initialized by optimizer.
  lat_prms.bn_tol = 1e-5; lat_prms.step = 1e0;

  set_unit_cell(prms, constr);
  set_unit_cell_b2(prms);
  set_unit_cell_b3_Fam(prms);

  lat_constr.eps0_x = eps0_x;
  set_unit_cell_constr(constr);
  // set_unit_cell_b3_constr(constr);
  set_unit_cell_weights(constr);

  lat_constr.ini_constr(true);

  prt_prms(lat_constr);
}


void set_dip_cell_disp(param_type &prms, constr_type &constr)
{
  std::vector<int>    grad_dip_Fnum, mI_loc;
  std::vector<double> grad_dip_scl;

  if (false) {
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
      prms.add_prm(grad_dip_Fnum, grad_dip_scl, -3, -phi_max, phi_max,
		   1.0);
    if (long_grad_dip[0])
      prms.add_prm(grad_dip_Fnum, grad_dip_scl,  2, -b_2_dip_max, b_2_dip_max,
		   1.0);

    lat_constr.Fnum_b1.push_back(-ElemIndex("dl1a_1"));
    lat_constr.grad_dip_Fnum_b1.push_back(grad_dip_Fnum);

    grad_dip_Fnum.clear();
    grad_dip_Fnum.push_back(ElemIndex("dl2a_1"));
    grad_dip_Fnum.push_back(ElemIndex("dl2a_2"));
    grad_dip_Fnum.push_back(ElemIndex("dl2a_3"));
    grad_dip_Fnum.push_back(ElemIndex("dl2a_4"));
    grad_dip_Fnum.push_back(ElemIndex("dl2a_5"));
    if (dphi)
      prms.add_prm(grad_dip_Fnum, grad_dip_scl, -3, -phi_max, phi_max,
		   1.0);
    if (long_grad_dip[1])
      prms.add_prm(grad_dip_Fnum, grad_dip_scl,  2, -b_2_dip_max, b_2_dip_max,
		   1.0);

    lat_constr.Fnum_b1.push_back(-ElemIndex("dl2a_1"));
    lat_constr.grad_dip_Fnum_b1.push_back(grad_dip_Fnum);
  }

  if (dphi) {
    // Dipole Cell.
    prms.add_prm("b1c", -3, -phi_max, phi_max, 1.0);
    lat_constr.Fnum_b1.push_back(ElemIndex("b1c"));
    prms.add_prm("b1e", -3, -phi_max, phi_max, 1.0);
    lat_constr.Fnum_b1.push_back(ElemIndex("b1e"));

    // Commented out must be defined last.
    // prms.add_prm("b2", -3, -5.0, 5.0,  1.0);
    lat_constr.Fnum_b1.push_back(ElemIndex("b2"));
  }

  prms.add_prm("b1c", 2, -b_2_max, b_2_max, 1.0);
  prms.add_prm("b1e", 2, -b_2_max, b_2_max, 1.0);
  prms.add_prm("b2",  2, -b_2_max, b_2_max, 1.0);
}


void set_b2_disp(param_type &prms) { }


void set_constr_disp(constr_type &constr)
{
  // constr.add_constr(Elem_GetPos(ElemIndex("qf1"), 1)-1,
  // 		    0e0, 0e0, 0e0, 0e0, 1e8, 1e7,
  // 		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  // constr.add_constr(Elem_GetPos(ElemIndex("qf1"), 2),
  // 		    0e0, 0e0, 0e0, 0e0, 1e8, 1e7,
  // 		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  // constr.add_constr(Elem_GetPos(ElemIndex("dl1a_5"), 1)-1,
  // 		    0e0, 0e0, 0e0, 0e0, 1e8, 1e7,
  // 		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  // constr.add_constr(Elem_GetPos(ElemIndex("dl1a_5"), 2),
  // 		    0e0, 0e0, 0e0, 0e0, 1e8, 1e7,
  // 		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  constr.add_constr(Elem_GetPos(ElemIndex("b1e"), 2),
  		    1e5, 1e5, 0e2, 0e2, 1e8, 1e7,
  		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  // Symmetry.
  // constr.add_constr(Elem_GetPos(ElemIndex("sf1"), 1),
  // 		    5e2, 1e3, 0e1,  0e2, 0e0, 0e0,
  // 		    0.0, 0.0, 10.0, 1.0, 0.0, 0.0);

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


void set_b3_Fam_disp(param_type &prms)
{
  std::vector<int> Fnum;

  prms.add_prm("s1",  3, -b_3_chrom_max, b_3_chrom_max, 1.0);
  prms.add_prm("s2",  3, -b_3_chrom_max, b_3_chrom_max, 1.0);

  // For control of linear chromaticity.
  lat_constr.Fnum_b3.push_back(ElemIndex("s1"));
  lat_constr.Fnum_b3.push_back(ElemIndex("s2"));

  if (true) {
    no_sxt();

    Fnum.push_back(ElemIndex("s1"));
    Fnum.push_back(ElemIndex("s2"));
    fit_ksi1(Fnum, 0e0, 0e0, 1e1);
  }
}


void set_b3_constr_disp(constr_type &constr)
{
  int              k, n;
  std::vector<int> mI_loc;

  // mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 1));
  // mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 3));
  // lat_constr.mI_loc.push_back(mI_loc);

  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("s1"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("s2"), 1));

  n = lat_constr.high_ord_achr_Fnum.size()/2;
  lat_constr.high_ord_achr_dnu.resize(n);
  lat_constr.high_ord_achr_nu.resize(n);
  for (k = 0; k < n; k++) {
    lat_constr.high_ord_achr_dnu[k].resize(2, 0e0);
    lat_constr.high_ord_achr_nu[k].resize(2, 0e0);
  }

  for (k = 0; k < 2; k++)
    lat_constr.high_ord_achr_nu[0][k] = high_ord_achr_nu[0][k];

  // for (k = 0; k < 2; k++)
  //   lat_constr.mI0[k] = mI_nu_ref[k];
}


void set_weights_disp(constr_type &constr)
{
  lat_constr.eps_x_scl             = 1e5;

  lat_constr.ksi1_ctrl_scl[0]      = 0e-1;
  lat_constr.ksi1_ctrl_scl[1]      = 0e0;
  lat_constr.ksi1_ctrl_scl[2]      = 0e0;
  lat_constr.mI_scl[X_]            = 0e6;
  lat_constr.mI_scl[Y_]            = 0e6;
  lat_constr.high_ord_achr_scl[X_] = high_ord_achr_scl[X_];
  lat_constr.high_ord_achr_scl[Y_] = high_ord_achr_scl[Y_];

  lat_constr.alpha_c_scl           = 1e-8;

  // Super Period.
  lat_constr.phi_scl               = ((dphi)? 1e0 : 0e0);
  lat_constr.phi0                  = 18.0;
  lat_constr.L_scl                 = 0e-10;
  lat_constr.L0                    = 10.0;

  lat_constr.drv_terms.get_h_scl(scl_ksi_1, scl_h_3, scl_h_3_delta, scl_h_4,
				 scl_ksi_2, scl_ksi_3, scl_chi_2,
				 scl_chi_delta_2);
}


void opt_disp(param_type &prms, constr_type &constr)
{
  // Parameter Type:
  //   Bend Angle  -3,
  //   Length      -2,
  //   Position    -1,
  //   Quadrupole   2.

  // Set parameters; initialized by optimizer.
  lat_prms.bn_tol = 1e-5; lat_prms.step = 1e0;

  set_dip_cell_disp(prms, constr);
  set_b2_disp(prms);
  set_b3_Fam_disp(prms);

  lat_constr.eps0_x = eps0_x;
  set_constr_disp(constr);
  // set_b3_constr(constr);
  set_weights_disp(constr);

  lat_constr.ini_constr(true);

  prt_prms(lat_constr);
}


void set_dip_cell(param_type &prms, constr_type &constr)
{
  std::vector<int>    grad_dip_Fnum, mI_loc;
  std::vector<double> grad_dip_scl;

  if (false) {
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
      prms.add_prm(grad_dip_Fnum, grad_dip_scl, -3, -phi_max, phi_max,
		   1.0);
    if (long_grad_dip[0])
      prms.add_prm(grad_dip_Fnum, grad_dip_scl,  2, -b_2_dip_max, b_2_dip_max,
		   1.0);

    lat_constr.Fnum_b1.push_back(-ElemIndex("dl1a_1"));
    lat_constr.grad_dip_Fnum_b1.push_back(grad_dip_Fnum);

    grad_dip_Fnum.clear();
    grad_dip_Fnum.push_back(ElemIndex("dl2a_1"));
    grad_dip_Fnum.push_back(ElemIndex("dl2a_2"));
    grad_dip_Fnum.push_back(ElemIndex("dl2a_3"));
    grad_dip_Fnum.push_back(ElemIndex("dl2a_4"));
    grad_dip_Fnum.push_back(ElemIndex("dl2a_5"));
    if (dphi)
      prms.add_prm(grad_dip_Fnum, grad_dip_scl, -3, -phi_max, phi_max,
		   1.0);
    if (long_grad_dip[1])
      prms.add_prm(grad_dip_Fnum, grad_dip_scl,  2, -b_2_dip_max, b_2_dip_max,
		   1.0);

    lat_constr.Fnum_b1.push_back(-ElemIndex("dl2a_1"));
    lat_constr.grad_dip_Fnum_b1.push_back(grad_dip_Fnum);
  }

  if (dphi) {
    // Dipole Cell.
    prms.add_prm("b1", -3, -phi_max, phi_max, 1.0);
    lat_constr.Fnum_b1.push_back(ElemIndex("b1"));

    // Commented out must be defined last.
    // prms.add_prm("d2", -3, -5.0, 5.0,  1.0);
    lat_constr.Fnum_b1.push_back(ElemIndex("b2"));
  }

  prms.add_prm("b1", 2, -b_2_max, b_2_max, 1.0);
  prms.add_prm("b2", 2, -b_2_max, b_2_max, 1.0);
}


void set_b2(param_type &prms)
{
  prms.add_prm("q1",  2, -b_2_max, b_2_max, 1.0);
  prms.add_prm("q2",  2, -b_2_max, b_2_max, 1.0);
  prms.add_prm("q3",  2, -b_2_max, b_2_max, 1.0);
}


void set_constr(constr_type &constr)
{
  // constr.add_constr(Elem_GetPos(ElemIndex("qf1"), 1)-1,
  // 		    0e0, 0e0, 0e0, 0e0, 1e8, 1e7,
  // 		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  // constr.add_constr(Elem_GetPos(ElemIndex("qf1"), 2),
  // 		    0e0, 0e0, 0e0, 0e0, 1e8, 1e7,
  // 		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  // constr.add_constr(Elem_GetPos(ElemIndex("dl1a_5"), 1)-1,
  // 		    0e0, 0e0, 0e0, 0e0, 1e8, 1e7,
  // 		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  // constr.add_constr(Elem_GetPos(ElemIndex("dl1a_5"), 2),
  // 		    0e0, 0e0, 0e0, 0e0, 1e8, 1e7,
  // 		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  constr.add_constr(Elem_GetPos(ElemIndex("cav"), 1),
  		    1e5, 1e5, 1e2,         1e2,         1e7, 1e7,
  		    0.0, 0.0, beta_ss[X_], beta_ss[Y_], 0.0, 0.0);
  // Symmetry.
  // constr.add_constr(Elem_GetPos(ElemIndex("sf1"), 1),
  // 		    5e2, 1e3, 0e1,  0e2, 0e0, 0e0,
  // 		    0.0, 0.0, 10.0, 1.0, 0.0, 0.0);

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


void set_b3_Fam(param_type &prms)
{
  std::vector<int> Fnum;

  prms.add_prm("s1",  3, -b_3_chrom_max, b_3_chrom_max, 1.0);
  prms.add_prm("s2",  3, -b_3_chrom_max, b_3_chrom_max, 1.0);

  // For control of linear chromaticity.
  lat_constr.Fnum_b3.push_back(ElemIndex("s1"));
  lat_constr.Fnum_b3.push_back(ElemIndex("s2"));

  if (true) {
    no_sxt();

    Fnum.push_back(ElemIndex("s1"));
    Fnum.push_back(ElemIndex("s2"));
    fit_ksi1(Fnum, 0e0, 0e0, 1e1);
  }
}


void set_b3_constr(constr_type &constr)
{
  int              k, n;
  std::vector<int> mI_loc;

  // mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 1));
  // mI_loc.push_back(Elem_GetPos(ElemIndex("sf1"), 3));
  // lat_constr.mI_loc.push_back(mI_loc);

  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("s1"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("s2"), 1));

  n = lat_constr.high_ord_achr_Fnum.size()/2;
  lat_constr.high_ord_achr_dnu.resize(n);
  lat_constr.high_ord_achr_nu.resize(n);
  for (k = 0; k < n; k++) {
    lat_constr.high_ord_achr_dnu[k].resize(2, 0e0);
    lat_constr.high_ord_achr_nu[k].resize(2, 0e0);
  }

  for (k = 0; k < 2; k++)
    lat_constr.high_ord_achr_nu[0][k] = high_ord_achr_nu[0][k];

  // for (k = 0; k < 2; k++)
  //   lat_constr.mI0[k] = mI_nu_ref[k];
}


void set_weights(constr_type &constr)
{
  lat_constr.eps_x_scl             = 1e7;

  lat_constr.ksi1_ctrl_scl[0]      = 0e-1;
  lat_constr.ksi1_ctrl_scl[1]      = 0e0;
  lat_constr.ksi1_ctrl_scl[2]      = 0e0;
  lat_constr.mI_scl[X_]            = 0e6;
  lat_constr.mI_scl[Y_]            = 0e6;
  lat_constr.high_ord_achr_scl[X_] = high_ord_achr_scl[X_];
  lat_constr.high_ord_achr_scl[Y_] = high_ord_achr_scl[Y_];

  lat_constr.alpha_c_scl           = 1e-6;

  // Super Period.
  lat_constr.phi_scl               = ((dphi)? 1e0 : 0e0);
  lat_constr.phi0                  = 18.0;
  lat_constr.L_scl                 = 0e-10;
  lat_constr.L0                    = 10.0;

  lat_constr.drv_terms.get_h_scl(scl_ksi_1, scl_h_3, scl_h_3_delta, scl_h_4,
				 scl_ksi_2, scl_ksi_3, scl_chi_2,
				 scl_chi_delta_2);
}


void opt_mI(param_type &prms, constr_type &constr)
{
  // Parameter Type:
  //   Bend Angle  -3,
  //   Length      -2,
  //   Position    -1,
  //   Quadrupole   2.

  // Set parameters; initialized by optimizer.
  lat_prms.bn_tol = 1e-5; lat_prms.step = 1e0;

  set_dip_cell(prms, constr);
  set_b2(prms);
  set_b3_Fam(prms);

  lat_constr.eps0_x = eps0_x;
  set_constr(constr);
  // set_b3_constr(constr);
  set_weights(constr);

  lat_constr.ini_constr(true);

  prt_prms(lat_constr);
}


void set_b2_ss(param_type &prms)
{
  // Std Straight.
  prms.add_prm("q1",  2, -b_2_max, b_2_max, 1.0);
  prms.add_prm("q2",  2, -b_2_max, b_2_max, 1.0);
  prms.add_prm("q3",  2, -b_2_max, b_2_max, 1.0);

  // Parameters are initialized in optimizer.
}


void set_constr_ss(constr_type &constr)
{
  constr.add_constr(Elem_GetPos(ElemIndex("l2"), 1),
		    1e2, 1e2, 1e-1,        1e-1,         1e-10, 1e-10,
		    0.0, 0.0, beta_ss[X_], beta_ss[Y_],  0.0,   0.0);
  constr.add_constr(Elem_GetPos(ElemIndex("q3"), 1),
		    0e0, 0e0, 0e0, 0e0, 1e6, 1e6,
		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
}


void set_b3_constr_ss(constr_type &constr)
{
  int k, n;

  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("s1"), 1));
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("s2"), 1));

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
    Ring_GetTwiss(true, 0e0); printglob();
    set_map(ElemIndex("ps_rot"), dnu);
    Ring_GetTwiss(true, 0e0); printglob();
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
    // Optimize Unit Cell.
    opt_unit_cell(lat_prms, lat_constr);
    fit_powell(lat_prms, 1e-3, f_achrom);
    // fit_sim_anneal(lat_prms, 1e-3, f_achrom);
    break;
  case 2:
    opt_disp(lat_prms, lat_constr);
    fit_powell(lat_prms, 1e-3, f_achrom);
    // fit_sim_anneal(lat_prms, 1e-3, f_achrom);
    break;
  case 3:
    match_ss(lat_prms, lat_constr);
    fit_powell(lat_prms, 1e-3, f_match);
    break;
  case 4:
    opt_mI(lat_prms, lat_constr);
    fit_powell(lat_prms, 1e-3, f_achrom);
    // fit_sim_anneal(lat_prms, 1e-3, f_achrom);
    break;
  default:
    printf("\nmain: unknown opt_case %d\n", opt_case);
    exit(1);
    break;
  }
}
