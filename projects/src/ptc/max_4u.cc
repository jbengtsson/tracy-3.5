#define NO 4

#include "tracy_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;


const double
  A_max[]   = {6e-3, 3e-3},
  delta_max = 3e-2;


void set_ps_rot(const double dnu_x, const double dnu_y)
{
  const double
    dnu_0[] = {0e0, 0e0},
    dnu[]   = {dnu_x, dnu_y};

  set_map(ElemIndex("ps_rot"), dnu_0);
  set_map(ElemIndex("ps_rot"), dnu);
}


tps get_h_local(const ss_vect<tps> &map, const bool dragt_finn)
{
  ss_vect<tps>  map1, R;

  if (dragt_finn)
    // Dragt-Finn factorization.
    return LieFact_DF(map, R);
  else {
    // Single Lie exponent.
    danot_(1);
    map1 = map;
    danot_(no_tps);
    return LieFact(map*Inv(map1));
  }
}


double h_abs_ijklm
(const tps &h_re, const tps &h_im, const int i, const int j, const int k,
 const int l, const int m)
{
  return sqrt(sqr(h_ijklm(h_re, i, j, k, l, m))
  	      +sqr(h_ijklm(h_im, i, j, k, l, m)));
}


std::vector<double> get_h_ijklm(ss_vect<tps> &Id_scl)
{
  long int            lastpos;
  std::vector<double> h_buf;
  tps                 h, h_re, h_im;
  ss_vect<tps>        map;

  danot_(no_tps-1);
  map.identity();
  Cell_Pass(0, globval.Cell_nLoc, map, lastpos);
  danot_(no_tps);
  h = get_h_local(map, true);
  CtoR(h*Id_scl, h_re, h_im);

  h_buf.push_back(h_abs_ijklm(h_re, h_im, 1, 1, 0, 0, 0));
  h_buf.push_back(h_abs_ijklm(h_re, h_im, 1, 1, 1, 1, 0));
  h_buf.push_back(h_abs_ijklm(h_re, h_im, 0, 0, 1, 1, 0));

  h_buf.push_back(h_abs_ijklm(h_re, h_im, 1, 1, 0, 0, 2));
  h_buf.push_back(h_abs_ijklm(h_re, h_im, 0, 0, 1, 1, 2));

  return h_buf;
}


void prt_h_ijklm
(ostream &outf, const std::vector<double> dnu,
 const std::vector<double> &h_ijklm)
{
  outf << fixed << setprecision(3) << setw(7) << dnu[X_]  << setw(7) << dnu[Y_];
  for (int k = 0; k < h_ijklm.size(); k++)
    outf << scientific << setprecision(3) << setw(11) << h_ijklm[k];
  outf << "\n";
}


double get_chi_2(std::vector<double> &h_ijklm)
{
  double chi_2 = 0e0;

  for (int k = 0; k < h_ijklm.size(); k++)
    chi_2 += sqr(h_ijklm[k]);
  return chi_2;
}

void scan_h_ijklm
(ss_vect<tps> &Id_scl, const int n, const double dnu_max_x,
 const double dnu_max_y)
{
  const double
    dnu_max[] = {dnu_max_x, dnu_max_y};
  const string
    file_name = "h_ijklm.dat"; 

  double              nu_step[2], chi_2, chi_2_min = 1e30;
  std::vector<double> h_ijklm, dnu, dnu_min, h_ijklm_min;
  ofstream            outf;

  file_wr(outf, file_name.c_str());

  for (int k = 0; k < 2; k++) {
    nu_step[k] = dnu_max[k]/n;
    dnu.push_back(-dnu_max[k]);
  }

  for (int k = 0; k < 5; k++)
    h_ijklm_min.push_back(1e30);

  printf("\n");
  for (int j = 0; j < 2*n+1; j++) {
    printf("  j = %d\n", j);
    for (int k = 0; k < 2*n+1; k++) {
      if (k == 0)
	outf << "\n";
      set_ps_rot(dnu[X_], dnu[Y_]);
      h_ijklm = get_h_ijklm(Id_scl);
      chi_2 = get_chi_2(h_ijklm);
      if (chi_2 < chi_2_min) {
	dnu_min = dnu;
	chi_2_min = chi_2;
      }
      for (int l = 0; l < 5; l++)
	h_ijklm_min[l] = min(fabs(h_ijklm[l]), h_ijklm_min[l]);

      prt_h_ijklm(cout, dnu, h_ijklm);
      prt_h_ijklm(outf, dnu, h_ijklm);

      dnu[Y_] += nu_step[Y_];
    }
    dnu[X_] += nu_step[X_];
    dnu[Y_] = -dnu_max[Y_];
  }

  outf.close();

  cout << scientific << setprecision(3) << "\ndh_ijklm rms = "
       << setw(9) << sqrt(chi_2_min) << "\n";
  set_ps_rot(dnu_min[X_], dnu_min[Y_]);
  h_ijklm = get_h_ijklm(Id_scl);
  prt_h_ijklm(cout, dnu_min, h_ijklm);
  cout << "\nh_ijklm_min  = ";
  for (int k = 0; k < h_ijklm_min.size(); k++)
    cout << scientific << setprecision(3) << setw(10) << h_ijklm_min[k];
  cout << "\n";
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
  double       twoJ;
  ss_vect<tps> Id_scl;

  set_state();

  // disable from TPSALib and LieLib log messages
  idprset(-1);

  daeps_(eps_tps);

  if (!true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  if (false) {
    Ring_GetTwiss(true, 0e0);
    printglob();
  }

  if (!false) {
    // Twiss functions are needed for set_ps_rot.
    Ring_GetTwiss(true, 0e0);
    printglob();

    Id_scl.identity();
    for (int k = 0; k < 2; k++) {
      twoJ = sqr(A_max[k])/Cell[0].Beta[k];
      Id_scl[2*k] *= sqrt(twoJ);
      Id_scl[2*k+1] *= sqrt(twoJ);
    }
    Id_scl[delta_] *= delta_max;

    scan_h_ijklm(Id_scl, 10, 1.0, 1.0);
  }

  if (false) {
    // Twiss functions are needed for set_ps_rot.
    Ring_GetTwiss(true, 0e0);
    printglob();
    set_ps_rot(0.1, -0.1);
    Ring_GetTwiss(true, 0e0);
    printglob();
    prtmfile("flat_file.dat");
  }
}
