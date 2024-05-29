#define NO 4

#include "tracy_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;


const double
  beta_inj[] = {5.0, 3.0},
  A_max[]    = {6e-3, 3e-3},
  delta_max  = 3e-2,
  twoJ[]     = {sqr(A_max[X_])/beta_inj[X_], sqr(A_max[Y_])/beta_inj[Y_]};


void set_ps_rot(const double dnu_x, const double dnu_y)
{
  const double
    dnu_0[] = {0e0, 0e0},
    dnu[]   = {dnu_x, dnu_y};

  set_map(ElemIndex("ps_rot"), dnu_0);
  set_map(ElemIndex("ps_rot"), dnu);
}


std::vector<double> get_k_ijklm(void)
{
  long int            lastpos;
  std::vector<double> k_ijklm;
  tps                 K_re, K_im;
  ss_vect<tps>        Id_scl, map;

  Id_scl.identity();
  for (int k = 0; k < 4; k++)
    Id_scl[k] *= sqrt(twoJ[k/2]);
  Id_scl[delta_] *= delta_max;

  danot_(no_tps-1);
  map.identity();
  Cell_Pass(0, globval.Cell_nLoc, map, lastpos);
  danot_(no_tps);
  MNF = MapNorm(map, 1);

  CtoR(MNF.K*Id_scl, K_re, K_im);

  k_ijklm.push_back(h_ijklm(K_re, 1, 1, 0, 0, 0));
  k_ijklm.push_back(h_ijklm(K_re, 1, 1, 1, 1, 0));
  k_ijklm.push_back(h_ijklm(K_re, 0, 0, 1, 1, 0));

  k_ijklm.push_back(h_ijklm(K_re, 1, 1, 0, 0, 2));
  k_ijklm.push_back(h_ijklm(K_re, 0, 0, 1, 1, 2));

  return k_ijklm;
}


void prt_k_ijklm
(ostream &outf, const std::vector<double> dnu,
 const std::vector<double> &k_ijklm)
{
  outf << fixed << setprecision(3) << setw(7) << dnu[X_]  << setw(7) << dnu[Y_];
  for (int k = 0; k < k_ijklm.size(); k++)
    outf << scientific << setprecision(3) << setw(11) << k_ijklm[k];
  outf << "\n";
}


double get_chi_2(std::vector<double> &k_ijklm)
{
  double chi_2 = 0e0;

  for (int k = 0; k < k_ijklm.size(); k++)
    chi_2 += sqr(k_ijklm[k]);
  return chi_2;
}

void scan_k_ijklm(const int n, const double dnu_max_x, const double dnu_max_y)
{
  const double
    dnu_max[] = {dnu_max_x, dnu_max_y};
  const string
    file_name = "k_ijklm.dat"; 

  double              nu_step[2], chi_2, chi_2_min = 1e30;
  std::vector<double> k_ijklm, dnu, dnu_min, k_ijklm_min;
  ofstream            outf;

  file_wr(outf, file_name.c_str());

  for (int k = 0; k < 2; k++) {
    nu_step[k] = dnu_max[k]/n;
    dnu.push_back(-dnu_max[k]);
  }

  for (int k = 0; k < 5; k++)
    k_ijklm_min.push_back(1e30);

  printf("\n");
  for (int j = 0; j < 2*n+1; j++) {
    printf("  j = %d\n", j);
    for (int k = 0; k < 2*n+1; k++) {
      if (k == 0)
	outf << "\n";
      set_ps_rot(dnu[X_], dnu[Y_]);
      k_ijklm = get_k_ijklm();
      chi_2 = get_chi_2(k_ijklm);
      if (chi_2 < chi_2_min) {
	dnu_min = dnu;
	chi_2_min = chi_2;
      }
      for (int l = 0; l < 5; l++)
	k_ijklm_min[l] = min(fabs(k_ijklm[l]), k_ijklm_min[l]);

      prt_k_ijklm(outf, dnu, k_ijklm);

      dnu[Y_] += nu_step[Y_];
    }
    dnu[X_] += nu_step[X_];
    dnu[Y_] = -dnu_max[Y_];
  }

  outf.close();

  cout << scientific << setprecision(3) << "\ndk_ijklm rms = "
       << setw(9) << sqrt(chi_2_min) << "\n";
  set_ps_rot(dnu_min[X_], dnu_min[Y_]);
  k_ijklm = get_k_ijklm();
  prt_k_ijklm(cout, dnu_min, k_ijklm);
  cout << "\nk_ijklm_min  = ";
  for (int k = 0; k < k_ijklm_min.size(); k++)
    cout << scientific << setprecision(3) << setw(10) << k_ijklm_min[k];
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

  if (false) {
    // Twiss functions are needed for set_ps_rot.
    Ring_GetTwiss(true, 0e0);
    printglob();
    scan_k_ijklm(10, 1.0, 1.0);
  }

  if (!false) {
    // Twiss functions are needed for set_ps_rot.
    Ring_GetTwiss(true, 0e0);
    printglob();
    set_ps_rot(0.1, -0.1);
    Ring_GetTwiss(true, 0e0);
    printglob();
    prtmfile("flat_file.dat");
  }
}
