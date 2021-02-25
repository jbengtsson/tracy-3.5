#define NO 1

#include "tracy_lib.h"

#include "Powell/src/newuoa.h"

#include "drv_terms.cc"

int no_tps = NO;


const bool
  ps_rot          = !false, // Note, needs to be zeroed; after use.
  pert_dip_cell   = false,
  dphi            = true,
  long_grad_dip   = false,
  ksi_terms[]     = {!false, !false, false},
  drv_terms[]     = {!false, !false, false},
  tune_fp_terms[] = {false, false};

/* opt_case:
     make_uc        1,
     cross_res      2,
     opt_unit_cell  3,
     match_disp     4, no ksi_terms & drv_terms, b1e single block dipole,
     match_straight 5,
     opt_sp         6.                                                        */
const int opt_case = 6
;

const int
  n_ic   = 5,
  n_cell = 4;

const double
  beta_ms[]             = {1.5, 1.5},

// From the center of the straight: alpha, beta, eta, eta'.
  ic[n_ic][2]           =
  {{-0.0000000000, -0.0000000000}, {3.2262890063, 2.4083826129},
   {0.0306340488, 0.0000000000}, {0.0, 0.0}},

#if 1
  eps0_x                = 0.125,
#else
  eps0_x                = 0.050,
#endif

  beta_inj[]            = {2.0, 2.0},
  A_max[]               = {6e-3, 3e-3},
  delta_max             = 2e-2,
  A_delta_max[]         = {2e-3, 0.1e-3},

  dnu[]                 = {0.0, -0.1},

// 3 for opt_unit_cell.
#define NU 1
#if NU == 1
  // 5 unit cells.
  nu_ref[]              = {1.0/2.0, 1.0/2.0},
  high_ord_achr_nu[][2] = {{2.0/5.0, 1.0/10.0}, {0.0, 0.0}},

#elif NU == 2
  // 4 unit cells.
  nu_ref[]              = {3.0/8.0, 1.0/8.0},
  high_ord_achr_nu[][2] = {{3.0/8.0, 1.0/8.0}, {0.0, 0.0}},
#elif NU == 3
  // Unit cell.
  nu_ref[]              = {2.0/5.0, 1.0/10.0},
  high_ord_achr_nu[][2] = {{0.0, 0.0}, {0.0, 0.0}},
#endif

  scl_eps_x             = 5e7,
  nu_ref_scl            = 1e-2*1e7,
  alpha_c_scl           = 1e0*5e-7,

  mI_scl                = 1e-6,

  scl_ksi_1             = (ksi_terms[0])?     1e0*1e0 : 0e0,
  scl_h_3               = (drv_terms[0])?     1e-10*1e14 : 0e0,
  scl_h_3_delta         = (drv_terms[1])?     1e14 : 0e0,
  scl_h_4               = (drv_terms[2])?     1e18 : 0e0,
  scl_ksi_2             = (ksi_terms[1])?     1e-1*1e7 : 0e0,
  scl_ksi_3             = (ksi_terms[2])?     1e0*1e7 : 0e0,
  scl_chi_2             = (tune_fp_terms[0])? 1e1*1e4 : 0e0,
  scl_chi_delta_2       = (tune_fp_terms[1])? 1e5 : 0e0,

  scl_extra             = 0e4,

  twoJ[]                =
    {sqr(A_max[X_])/beta_inj[X_], sqr(A_max[Y_])/beta_inj[Y_]},
  twoJ_delta[]          =
    {sqr(A_delta_max[X_])/beta_inj[X_], sqr(A_delta_max[Y_])/beta_inj[Y_]},

  mI_nu_ref[]           = {1.5, 0.5},

  phi_max               = 3.0,
  b_2_max               = 15.0,
  b_2_dip_max           = 4.0,
  b_3_chrom_max         = 1.5e3,
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
      prt_f(bn, chi2, lat_constr, lat_prms, true, false);
    }
    lat_constr.n_iter++;
    lat_constr.chi2 = chi2;
  }

  return chi2;
}


double f_res(double *bn)
{
  long int    lastpos;
  int         k;
  double      chi2;
  static bool first = true;

  const int n_prt = 5;

  lat_prms.set_prm(bn);

  if (lat_constr.phi_scl != 0e0) phi_corr(lat_constr);

  getcod(0e0, lastpos);

  for (k = 0; k < 2; k++) {
    lat_constr.nu_cos[k] = globval.OneTurnMat[2*k][2*k];
    lat_constr.nu_sin[k] =
      sqrt(fabs(globval.OneTurnMat[2*k][2*k+1]*globval.OneTurnMat[2*k+1][2*k]));
    if (globval.OneTurnMat[2*k][2*k+1] < 0e0) lat_constr.nu_sin[k] *= -1e0;
  }

  globval.Alphac = globval.OneTurnMat[ct_][delta_]/Cell[globval.Cell_nLoc].S;

  chi2 =
    lat_constr.get_chi2(twoJ, delta_max, twoJ_delta, lat_prms, bn, false,
			false);

  if (first || (chi2 < lat_constr.chi2)) {
    first = false;
    if (lat_constr.n_iter % n_prt == 0) {
      printf("\n");
      prtmat(6, globval.OneTurnMat);
      lat_constr.get_chi2(twoJ, delta_max, twoJ_delta, lat_prms, bn, false,
			  true);
      prt_f(bn, chi2, lat_constr, lat_prms, false, false);
    }

    lat_constr.n_iter++;
    lat_constr.chi2 = chi2;
  }

  return chi2;
}


double f_achrom(double *bn)
{
  bool        stable;
  double      chi2, L_b1e, L_b1c, dchi2;
  static bool first = true;

  const bool prt   = false;
  const int  n_prt = 5;

  lat_prms.set_prm(bn);

  if (lat_constr.phi_scl != 0e0) phi_corr(lat_constr);

  if (lat_constr.ring) stable = get_nu(lat_constr.nu);
  if (prt) printf("\nf_achrom: stable = %d\n", stable);

  if ((lat_constr.ring && stable) || !lat_constr.ring) {
    eps_x = get_lin_opt(lat_constr);

    fit_ksi1(lat_constr.Fnum_b3, 0e0, 0e0, 1e1);

    chi2 =
      lat_constr.get_chi2(twoJ, delta_max, twoJ_delta, lat_prms, bn, true,
			  false);

    if (false) {
      L_b1e = get_L(ElemIndex("b1e"), 1);
      L_b1c = get_L(ElemIndex("b1c"), 1);
      dchi2 = scl_extra*sqr(L_b1e-L_b1c);
      chi2 += dchi2;
      printf("\n  L_b1e = %5.3f L_b1c = %5.3f (%9.3e)\n", L_b1e, L_b1c, dchi2);
    }

    if (first || (chi2 < lat_constr.chi2)) {
      first = false;
      if (lat_constr.n_iter % n_prt == 0) {
	// Print dchi2.
	lat_constr.get_chi2(twoJ, delta_max, twoJ_delta, lat_prms, bn, true,
			    true);
	prt_f(bn, chi2, lat_constr, lat_prms, true, true);
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
      prt_f(bn, chi2, lat_constr, lat_prms, true, true);
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
	prt_f(bn, chi2, lat_constr, lat_prms, true, true);
      }

      lat_constr.n_iter++;
      lat_constr.chi2 = chi2;
    }
  } else
    chi2 = 1e30;

  return chi2;
}


void set_b2_uc(param_type &prms)
{
  const double ds_max = 0.3;

  prms.add_prm("q1",  2, -b_2_max, b_2_max, 1.0);
  prms.add_prm("q2",  2, -b_2_max, b_2_max, 1.0);
  prms.add_prm("q3",  2, -b_2_max, b_2_max, 1.0);

  if (false) {
    prms.add_prm("q1", -1, 0.1-ds_max, 0.1-ds_max, 1.0);
    prms.add_prm("q2", -1, 0.1-ds_max, 0.1-ds_max, 1.0);
    prms.add_prm("q3", -1, 0.1-ds_max, 2.5-ds_max, 1.0);
  }

  prms.add_prm("b1e",  2, -b_2_max, b_2_max, 1.0);
  prms.add_prm("b1e", -2, 0.2,      1.0,     1.0);

  prms.add_prm("b1c",  2, -b_2_max, b_2_max, 1.0);
  prms.add_prm("b2",   2, -b_2_max, b_2_max, 1.0);

  // Parameters are initialized in optimizer.
}


void set_constr_uc(constr_type &constr)
{
  // Symmetry.
  constr.add_constr(Elem_GetPos(ElemIndex("s1"), 1),
		    1e5, 1e5, 0e0, 0e0, 0e0, 1e7,
		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  // Center of dipole.
  constr.add_constr(Elem_GetPos(ElemIndex("b1c"), 1),
		    1e5, 1e5, 0e0, 0e0, 0e0, 1e7,
		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
}


void set_b3_constr_uc(constr_type &constr)
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


void set_weights_uc(constr_type &constr)
{
  lat_constr.high_ord_achr_scl[X_] = 0e0;
  lat_constr.high_ord_achr_scl[Y_] = 0e0;
}


void make_uc(param_type &prms, constr_type &constr)
{
  // Parameter Type:
  //   Bend Angle  -3,
  //   Length      -2,
  //   Position    -1,
  //   Quadrupole   2.
  int j, k;

  lat_prms.bn_tol = 1e-5; lat_prms.step = 1.0;
  set_b2_uc(prms);

  set_constr_uc(constr);
  set_b3_constr_uc(constr);
  set_weights_uc(constr);

  lat_constr.ini_constr(false);

  for (j = 0; j < n_ic; j++)
    for (k = 0; k < 2; k++)
      lat_constr.ic[j][k] = ic[j][k];
}


void set_res(param_type &prms, constr_type &constr)
{
  std::vector<int>    grad_dip_Fnum, mI_loc;
  std::vector<double> grad_dip_scl;

  const double ds_max = 0.1;

  if (dphi) {
    // Dipole Cell.
    prms.add_prm("b1c", -3, 0.0, phi_max, 1.0);
    lat_constr.Fnum_b1.push_back(ElemIndex("b1c"));

    // Commented out must be defined last.
    // prms.add_prm("b2", -3, -5.0, 5.0,  1.0);
    lat_constr.Fnum_b1.push_back(ElemIndex("b2"));
  }

  prms.add_prm("b1c", 2, -b_2_max, b_2_max, 1.0);
  prms.add_prm("b2", 2, -b_2_max, b_2_max, 1.0);

  if (!false) {
    prms.add_prm("b2", -1, 0.1-ds_max, 0.1-ds_max, 1.0);
    prms.add_prm("s2", -1, 0.1-ds_max, 0.1-ds_max, 1.0);
  }
}


void set_b2_res(param_type &prms)
{
  // prms.add_prm("s1", 2, -b_2_max, b_2_max, 1.0);
  // prms.add_prm("s2", 2, -b_2_max, b_2_max, 1.0);

  // Parameters are initialized in optimizer.
}


void set_constr_res(constr_type &constr)
{
  // constr.add_constr(Elem_GetPos(ElemIndex("s1"), 1),
  // 		    1e-3*1e5, 1e-3*1e5, 0e0, 0e0, 0e0, 1e-3*1e7,
  // 		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
}


void set_b3_Fam_res(param_type &prms)
{
  std::vector<int> Fnum;

  // Seg. fault if removed.
  prms.add_prm("s1",3, -b_3_chrom_max, b_3_chrom_max, 1.0);
  prms.add_prm("s2",3, -b_3_chrom_max, b_3_chrom_max, 1.0);

  // Control of linear chromaticity.
  lat_constr.Fnum_b3.push_back(ElemIndex("s1"));
  lat_constr.Fnum_b3.push_back(ElemIndex("s2"));
}


void set_b3_constr_res(constr_type &constr)
{
  int k, n;

  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("b1c"), 1)-1);
  lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("b1c"), 2));

  n = lat_constr.high_ord_achr_Fnum.size()/2;
  lat_constr.high_ord_achr_dnu.resize(n);
  lat_constr.high_ord_achr_nu.resize(n);
  for (k = 0; k < n; k++) {
    lat_constr.high_ord_achr_dnu[k].resize(2, 0e0);
    lat_constr.high_ord_achr_nu[k].resize(2, 0e0);
  }

  for (k = 0; k < 2; k++)
    lat_constr.high_ord_achr_nu[0][k] = high_ord_achr_nu[0][k];
}


void set_weights_res(constr_type &constr)
{
  int k;

  lat_constr.nu_cos_ref_scl[X_] = 1e3;
  lat_constr.nu_cos_ref_scl[Y_] = 1e1;
  lat_constr.alpha_c_scl        = 1e-9;

  lat_constr.nu_ref_scl         = nu_ref_scl;

  for (k = 0; k < 2; k++) {
    lat_constr.nu_ref[k] = nu_ref[k];
    lat_constr.nu_cos_ref[k] = cos(2e0*M_PI*lat_constr.nu_ref[k]);
    lat_constr.nu_sin_ref[k] = sin(2e0*M_PI*lat_constr.nu_ref[k]);
  }

  // Super Period.
  lat_constr.phi_scl            = ((dphi)? 1e0 : 0e0);
  lat_constr.phi0               = 4.5;

  // Seg. fault if removed.
  lat_constr.drv_terms.get_h_scl(scl_ksi_1, scl_h_3, scl_h_3_delta, scl_h_4,
				 scl_ksi_2, scl_ksi_3, scl_chi_2,
				 scl_chi_delta_2);
}


void cross_res(param_type &prms, constr_type &constr)
{
  // Parameter Type:
  //   Bend Angle  -3,
  //   Length      -2,
  //   Position    -1,
  //   Quadrupole   2.

  // Set parameters; initialized by optimizer.
  lat_prms.bn_tol = 1e-5; lat_prms.step = 1e0;

  set_res(prms, constr);
  set_b2_res(prms);
  set_b3_Fam_res(prms);

  lat_constr.eps0_x = eps0_x;
  set_constr_res(constr);
  set_b3_constr_res(constr);
  set_weights_res(constr);

  lat_constr.ini_constr(true);

  prt_prms(lat_constr);
}


void set_unit_cell(param_type &prms, constr_type &constr)
{
  std::vector<int>    grad_dip_Fnum, mI_loc;
  std::vector<double> grad_dip_scl;

  // Dipole Cell.
  if (dphi) {
    if (!long_grad_dip) {
      prms.add_prm("b1c", -3, -phi_max, phi_max, 1.0);
      lat_constr.Fnum_b1.push_back(ElemIndex("b1c"));

      // Commented out must be defined last.
      // prms.add_prm("b2", -3, -5.0, 5.0,  1.0);
      lat_constr.Fnum_b1.push_back(ElemIndex("b2"));

      prms.add_prm("b1c",  2, -b_2_max, b_2_max, 1.0);

      prms.add_prm("b2",   2, -b_2_max, b_2_max, 1.0);
    } else {
      prms.add_prm("b1c1", -3, -phi_max, phi_max, 1.0);
      lat_constr.Fnum_b1.push_back(ElemIndex("b1c1"));
      prms.add_prm("b1c2", -3, -phi_max, phi_max, 1.0);
      lat_constr.Fnum_b1.push_back(ElemIndex("b1c2"));
      prms.add_prm("b1c3", -3, -phi_max, phi_max, 1.0);
      lat_constr.Fnum_b1.push_back(ElemIndex("b1c3"));
      prms.add_prm("b1c4", -3, -phi_max, phi_max, 1.0);
      lat_constr.Fnum_b1.push_back(ElemIndex("b1c4"));
      prms.add_prm("b1c5", -3, -phi_max, phi_max, 1.0);
      lat_constr.Fnum_b1.push_back(ElemIndex("b1c5"));

      // Commented out must be defined last.
      // prms.add_prm("b2", -3, -5.0, 5.0,  1.0);
      lat_constr.Fnum_b1.push_back(ElemIndex("b2"));

      prms.add_prm("b1c1",  2, -b_2_max, b_2_max, 1.0);
      prms.add_prm("b1c2",  2, -b_2_max, b_2_max, 1.0);
      prms.add_prm("b1c3",  2, -b_2_max, b_2_max, 1.0);
      prms.add_prm("b1c4",  2, -b_2_max, b_2_max, 1.0);
      prms.add_prm("b1c5",  2, -b_2_max, b_2_max, 1.0);

      prms.add_prm("b2",   2, -b_2_max, b_2_max, 1.0);
    }
  }
}


void set_b2_unit_cell(param_type &prms)
{
  // Parameters are initialized in optimizer.
}


void set_constr_unit_cell(constr_type &constr)
{
  if (!false)
    // Zero eta_x & eta'_x to simplify matching.
    constr.add_constr(Elem_GetPos(ElemIndex("b1c1"), 2),
		      0e0, 0e0, 0e0, 0e0, 1e8, 1e8,
		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
}


void set_b3_Fam_unit_cell(param_type &prms)
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


void set_b3_constr_unit_cell(constr_type &constr)
{
  int k, n;

  if (!long_grad_dip) {
    lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("b1c"), 1));
    lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("b1c"), 3));
  } else {
    lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("b1c1"), 1));
    lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("b1c1"), 3));
  }

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


void set_weights_unit_cell(constr_type &constr)
{
  lat_constr.eps_x_scl             = 1e-3*1e7;

  lat_constr.alpha_c_scl           = 1e-3*1e-6;

  lat_constr.nu_ref_scl            = nu_ref_scl;
  lat_constr.nu_ref[X_]            = nu_ref[X_];
  lat_constr.nu_ref[Y_]            = nu_ref[Y_];

  // Super Period.
  lat_constr.phi_scl               = ((dphi)? 1e0 : 0e0);

#define N_SP 18
#if N_SP == 20
  // 20 super periods 5-BA.
  lat_constr.phi0                  = 18.0/4.0;
#elif N_SP == 18
  // 18 super periods 6-BA.
  lat_constr.phi0                  = 20.0/5.0;
#elif N_SP == 16
  // 16 super periods 6-BA.
  lat_constr.phi0                  = 22.5/5.0;
#endif

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
  set_b2_unit_cell(prms);
  set_b3_Fam_unit_cell(prms);

  lat_constr.eps0_x = eps0_x;
  set_constr_unit_cell(constr);
  set_b3_constr_unit_cell(constr);
  set_weights_unit_cell(constr);

  lat_constr.ini_constr(true);

  prt_prms(lat_constr);
}


void set_dip_disp(param_type &prms, constr_type &constr)
{
  if (dphi) {
    if (!long_grad_dip) {
      prms.add_prm("b1e",  2, -b_2_dip_max, b_2_dip_max, 1.0);
      if (!false) prms.add_prm("b1e", -2, 0.3, 0.7, 1.0);
    } else {
      prms.add_prm("b1e1",  2, -b_2_dip_max, b_2_dip_max, 1.0);
      prms.add_prm("b1e2",  2, -b_2_dip_max, b_2_dip_max, 1.0);
      prms.add_prm("b1e3",  2, -b_2_dip_max, b_2_dip_max, 1.0);
      prms.add_prm("b1e4",  2, -b_2_dip_max, b_2_dip_max, 1.0);
      prms.add_prm("b1e5",  2, -b_2_dip_max, b_2_dip_max, 1.0);
    }
  }
}


void set_b2_disp(param_type &prms)
{
  // Parameters are initialized in optimizer.
}


void set_constr_disp(constr_type &constr)
{
  constr.add_constr(Elem_GetPos(ElemIndex("b1e1"), 1),
		    0e0, 0e0, 0e0, 0e0, 1e6, 1e6,
		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
}


void set_b3_constr_disp(constr_type &constr)
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


void set_weights_disp(constr_type &constr)
{
}


void match_disp(param_type &prms, constr_type &constr)
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
  set_dip_disp(prms, constr);
  set_b2_disp(prms);

  set_constr_disp(constr);
  set_b3_constr_disp(constr);
  set_weights_disp(constr);

  lat_constr.ini_constr(false);

  for (j = 0; j < n_ic; j++)
    for (k = 0; k < 2; k++)
      lat_constr.ic[j][k] = ic[j][k];
}


void set_dip_straight(param_type &prms, constr_type &constr)
{
  // if (dphi) {
  //   if (!long_grad_dip) {
  //     prms.add_prm("b1e",  2, -b_2_dip_max, b_2_dip_max, 1.0);
  //     if (!false) prms.add_prm("b1e", -2, 0.3, 0.7, 1.0);
  //   } else {
  //     prms.add_prm("b1e1",  2, -b_2_dip_max, b_2_dip_max, 1.0);
  //     prms.add_prm("b1e2",  2, -b_2_dip_max, b_2_dip_max, 1.0);
  //     prms.add_prm("b1e3",  2, -b_2_dip_max, b_2_dip_max, 1.0);
  //     prms.add_prm("b1e4",  2, -b_2_dip_max, b_2_dip_max, 1.0);
  //     prms.add_prm("b1e5",  2, -b_2_dip_max, b_2_dip_max, 1.0);
  //   }
  // }
}


void set_b2_straight(param_type &prms)
{
  prms.add_prm("q1",  2, -b_2_max, 0.0,     1.0);
  prms.add_prm("q2",  2,  0.0,     b_2_max, 1.0);
  prms.add_prm("q3",  2, -b_2_max, 0.0,     1.0);

  // Min & max are min distance to next element.
  if (!false) {
    // Introduce a tiny interval.
    prms.add_prm("q1", -1, 0.1, 0.1, 1.0);
    prms.add_prm("q2", -1, 0.1, 0.1, 1.0);
    prms.add_prm("q3", -1, 0.1, 2.1, 1.0);
  }

  // Parameters are initialized in optimizer.
}


void set_constr_straight(constr_type &constr)
{
  constr.add_constr(Elem_GetPos(ElemIndex("l7"), 1),
  		    1e5, 1e5, 0e0, 0e0, 0e0, 0e0,
  		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
}


void set_b3_constr_straight(constr_type &constr)
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


void set_weights_straight(constr_type &constr)
{
}


void match_straight(param_type &prms, constr_type &constr)
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
  set_dip_straight(prms, constr);
  set_b2_straight(prms);

  set_constr_straight(constr);
  set_b3_constr_straight(constr);
  set_weights_straight(constr);

  lat_constr.ini_constr(false);

  for (j = 0; j < n_ic; j++)
    for (k = 0; k < 2; k++)
      lat_constr.ic[j][k] = ic[j][k];
}


void set_dip_sp(param_type &prms, constr_type &constr)
{
  if (dphi) {
    if (!long_grad_dip) {
      prms.add_prm("b0", -3, -phi_max, phi_max, 1.0);
      lat_constr.Fnum_b1.push_back(ElemIndex("b0"));

      prms.add_prm("b1", -3, -phi_max, phi_max, 1.0);
      lat_constr.Fnum_b1.push_back(ElemIndex("b1"));

      prms.add_prm("mb0", -3, -5.0, 5.0,  1.0);
      lat_constr.Fnum_b1.push_back(ElemIndex("mb0"));
      // Commented out must be defined last.
      prms.add_prm("mb1", -3, -5.0, 5.0,  1.0);
      lat_constr.Fnum_b1.push_back(ElemIndex("mb1"));

      prms.add_prm("b0",  2, -b_2_dip_max, b_2_dip_max, 1.0);
      if (false) prms.add_prm("b0", -2, 0.3, 0.7, 1.0);
      prms.add_prm("b1",  2, -b_2_max,     b_2_max,     1.0);

      prms.add_prm("mb0", 2, -b_2_max, b_2_max, 1.0);
      prms.add_prm("mb1", 2, -b_2_max, b_2_max, 1.0);
    } else {
      prms.add_prm("b1e1", -3, -phi_max, phi_max, 1.0);
      lat_constr.Fnum_b1.push_back(ElemIndex("b1e1"));
      prms.add_prm("b1e2", -3, -phi_max, phi_max, 1.0);
      lat_constr.Fnum_b1.push_back(ElemIndex("b1e2"));
      prms.add_prm("b1e3", -3, -phi_max, phi_max, 1.0);
      lat_constr.Fnum_b1.push_back(ElemIndex("b1e3"));
      prms.add_prm("b1e4", -3, -phi_max, phi_max, 1.0);
      lat_constr.Fnum_b1.push_back(ElemIndex("b1e4"));
      prms.add_prm("b1e5", -3, -phi_max, phi_max, 1.0);
      lat_constr.Fnum_b1.push_back(ElemIndex("b1e5"));

      prms.add_prm("b1c1", -3, -phi_max, phi_max, 1.0);
      lat_constr.Fnum_b1.push_back(ElemIndex("b1c1"));
      prms.add_prm("b1c2", -3, -phi_max, phi_max, 1.0);
      lat_constr.Fnum_b1.push_back(ElemIndex("b1c2"));
      prms.add_prm("b1c3", -3, -phi_max, phi_max, 1.0);
      lat_constr.Fnum_b1.push_back(ElemIndex("b1c3"));
      prms.add_prm("b1c4", -3, -phi_max, phi_max, 1.0);
      lat_constr.Fnum_b1.push_back(ElemIndex("b1c4"));
      prms.add_prm("b1c5", -3, -phi_max, phi_max, 1.0);
      lat_constr.Fnum_b1.push_back(ElemIndex("b1c5"));

      // Commented out must be defined last.
      // prms.add_prm("b2", -3, -5.0, 5.0,  1.0);
      lat_constr.Fnum_b1.push_back(ElemIndex("b2"));

      prms.add_prm("b1e1",  2, -b_2_dip_max, b_2_dip_max, 1.0);
      prms.add_prm("b1e2",  2, -b_2_dip_max, b_2_dip_max, 1.0);
      prms.add_prm("b1e3",  2, -b_2_dip_max, b_2_dip_max, 1.0);
      prms.add_prm("b1e4",  2, -b_2_dip_max, b_2_dip_max, 1.0);
      prms.add_prm("b1e5",  2, -b_2_dip_max, b_2_dip_max, 1.0);

      prms.add_prm("b1c1",  2, -b_2_dip_max, b_2_dip_max, 1.0);
      prms.add_prm("b1c2",  2, -b_2_dip_max, b_2_dip_max, 1.0);
      prms.add_prm("b1c3",  2, -b_2_dip_max, b_2_dip_max, 1.0);
      prms.add_prm("b1c4",  2, -b_2_dip_max, b_2_dip_max, 1.0);
      prms.add_prm("b1c5",  2, -b_2_dip_max, b_2_dip_max, 1.0);

      prms.add_prm("b2",   2, -b_2_max, b_2_max, 1.0);
    }
  }
}


void set_b2_sp(param_type &prms)
{
  prms.add_prm("uq1",  2, -b_2_max, b_2_max, 1.0);
  prms.add_prm("uq2",  2, -b_2_max, b_2_max, 1.0);
  prms.add_prm("uq3",  2, -b_2_max, b_2_max, 1.0);

  // Parameters are initialized in optimizer.
}


void set_constr_sp(constr_type &constr)
{
  int k;

  // Beta functions at center of straight
  constr.add_constr(Elem_GetPos(ElemIndex("ul4"), 2),
  		    1e5, 1e5, 1e2,         1e2,         1e7, 1e7,
  		    0.0, 0.0, beta_ms[X_], beta_ms[Y_], 0.0, 0.0);
  // Symmetry.
  constr.add_constr(Elem_GetPos(ElemIndex("s0"), 1),
		    1e5, 1e5, 0e0, 0e0, 0e0, 1e8,
		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  constr.add_constr(Elem_GetPos(ElemIndex("s0"), 3),
		    1e5, 1e5, 0e0, 0e0, 0e0, 1e8,
		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

  if (!long_grad_dip) {
    // Symmetry & eta' at dipole centers.
    constr.add_constr(Elem_GetPos(ElemIndex("b0"), 1),
		      1e5, 1e5, 0e0, 0e0, 0e0, 1e8,
		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("b0"), 3),
		      1e5, 1e5, 0e0, 0e0, 0e0, 1e8,
		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    // Eta & eta' at arc entrance.
    constr.add_constr(Elem_GetPos(ElemIndex("mb0"), 1)-1,
		      0e0, 0e0, 0e0, 0e0, 1e8, 1e8,
		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  } else {
    // Symmetry & eta' at dipole centers.
    constr.add_constr(Elem_GetPos(ElemIndex("b1c1"), 1),
		      1e5, 1e5, 0e0, 0e0, 0e0, 1e8,
		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    constr.add_constr(Elem_GetPos(ElemIndex("b1c1"), 3),
		      1e5, 1e5, 0e0, 0e0, 0e0, 1e8,
		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    // Eta & eta' at arc entrance.
    constr.add_constr(Elem_GetPos(ElemIndex("b1e5"), 1)-1,
		      0e0, 0e0, 0e0, 0e0, 1e8, 1e8,
		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    // Eta at dipole centres.
    for (k = 1; k <= 5; k += 2)
      constr.add_constr(Elem_GetPos(ElemIndex("b1c1"), k),
			0e0, 0e0, 0e0, 0e0, 0e7, 0e0,
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  }
}


void set_b3_Fam_sp(param_type &prms)
{
  std::vector<int> Fnum;

  prms.add_prm("s0",  3, -b_3_chrom_max, b_3_chrom_max, 1.0);
  prms.add_prm("s1",  3, -b_3_chrom_max, b_3_chrom_max, 1.0);

  // For control of linear chromaticity.
  lat_constr.Fnum_b3.push_back(ElemIndex("s0"));
  lat_constr.Fnum_b3.push_back(ElemIndex("s1"));

  if (true) {
    no_sxt();

    Fnum.push_back(ElemIndex("s0"));
    Fnum.push_back(ElemIndex("s1"));
    fit_ksi1(Fnum, 0e0, 0e0, 1e1);
  }
}


void set_b3_constr_sp(constr_type &constr)
{
  int k, n;

  if (!long_grad_dip) {
    lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("b0"), 1));
    lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("b0"), 3));
  } else {
    lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("b1c1"), 1));
    lat_constr.high_ord_achr_Fnum.push_back(Elem_GetPos(ElemIndex("b1c1"), 3));
  }

  n = lat_constr.high_ord_achr_Fnum.size()/2;
  lat_constr.high_ord_achr_dnu.resize(n);
  lat_constr.high_ord_achr_nu.resize(n);
  for (k = 0; k < n; k++) {
    lat_constr.high_ord_achr_dnu[k].resize(2, 0e0);
    lat_constr.high_ord_achr_nu[k].resize(2, 0e0);
  }

  for (k = 0; k < 2; k++)
    lat_constr.high_ord_achr_nu[0][k] = high_ord_achr_nu[0][k];
}


void set_weights_sp(constr_type &constr)
{
  lat_constr.eps_x_scl             = 1e-3*1e7;

  lat_constr.alpha_c_scl           = 1e1*1e-6;

  lat_constr.nu_ref_scl            = nu_ref_scl;
  lat_constr.nu_ref[X_]            = 2*nu_ref[X_];
  lat_constr.nu_ref[Y_]            = 2*nu_ref[Y_];

  lat_constr.high_ord_achr_scl[X_] = 1e1*1e6;
  lat_constr.high_ord_achr_scl[Y_] = 1e1*1e6;

  // Super Period.
  lat_constr.phi_scl               = ((dphi)? 1e0 : 0e0);
#define N_SP 16
#if N_SP == 20
  // 20 super periods 5-BA.
  lat_constr.phi0                  = 2*18.0;
#elif N_SP == 18
  // 18 super periods 6-BA.
  lat_constr.phi0                  = 2*20.0;
#else
  // 16 super periods 6-BA.
  lat_constr.phi0                  = 2*22.5;
#endif
  lat_constr.L_scl                 = 0e-10;
  lat_constr.L0                    = 10.0;

  lat_constr.drv_terms.get_h_scl(scl_ksi_1, scl_h_3, scl_h_3_delta, scl_h_4,
				 scl_ksi_2, scl_ksi_3, scl_chi_2,
				 scl_chi_delta_2);
}


void opt_sp(param_type &prms, constr_type &constr)
{
  // Parameter Type:
  //   Bend Angle  -3,
  //   Length      -2,
  //   Position    -1,
  //   Quadrupole   2.

  // Set parameters; initialized by optimizer.
  lat_prms.bn_tol = 1e-5; lat_prms.step = 1e0;

  set_dip_sp(prms, constr);
  set_b2_sp(prms);
  set_b3_Fam_sp(prms);

  lat_constr.eps0_x = eps0_x;
  set_constr_sp(constr);
  set_b3_constr_sp(constr);
  set_weights_sp(constr);

  lat_constr.ini_constr(true);

  prt_prms(lat_constr);
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
    // Ring_GetTwiss(true, 0e0); printglob();
    // set_map(ElemIndex("ps_rot"), dnu0);
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
    make_uc(lat_prms, lat_constr);
    fit_powell(lat_prms, 1e-3, f_match);
    break;
  case 2:
    cross_res(lat_prms, lat_constr);
    fit_powell(lat_prms, 1e-3, f_res);
    break;
  case 3:
    // Optimize Unit Cell.
    opt_unit_cell(lat_prms, lat_constr);
    fit_powell(lat_prms, 1e-3, f_achrom);
    // fit_sim_anneal(lat_prms, 1e-3, f_achrom);
    break;
  case 4:
    match_disp(lat_prms, lat_constr);
    fit_powell(lat_prms, 1e-3, f_match);
    // fit_sim_anneal(lat_prms, 1e-3, f_achrom);
    break;
  case 5:
    match_straight(lat_prms, lat_constr);
    fit_powell(lat_prms, 1e-3, f_match);
    // fit_sim_anneal(lat_prms, 1e-3, f_achrom);
    break;
  case 6:
    opt_sp(lat_prms, lat_constr);
    fit_powell(lat_prms, 1e-3, f_achrom);
    // fit_sim_anneal(lat_prms, 1e-3, f_achrom);
    break;
  default:
    printf("\nmain: unknown opt_case %d\n", opt_case);
    exit(1);
    break;
  }
}
