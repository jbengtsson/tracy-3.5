#include <assert.h>

#define NO 4

#include "tracy_lib.h"


int
  no_tps   = NO,
  ndpt_tps = 5;


const double
  A_max[]   = {6e-3, 3e-3},
  delta_max = 4e-2;


void set_ps_rot(const string &fam_name, const double dnu_x, const double dnu_y)
{
  const int
    Fnum = ElemIndex(fam_name.c_str());
  const double
    dnu_0[] = {0e0, 0e0},
    dnu[]   = {dnu_x, dnu_y};

  set_map(Fnum, dnu_0);
  set_map(Fnum, dnu);
}


ss_vect<tps> get_map(void)
{
  long int     lastpos;
  ss_vect<tps> map;

  map.identity();
  Cell_Pass(0, globval.Cell_nLoc, map, lastpos);
  return map;
}


double k_abs_ijklm
(const tps &k_re, const tps &k_im, const int i, const int j, const int k,
 const int l, const int m)
{
  return sqrt(sqr(h_ijklm(k_re, i, j, k, l, m))
  	      +sqr(h_ijklm(k_im, i, j, k, l, m)));
}


std::vector<double> get_h_ijklm(const ss_vect<tps> &Id_scl)
{
  std::vector<double> k_buf;
  tps                 k_re, k_im;
  ss_vect<tps>        map;
  MNF_struct          MNF;

  danot_(no_tps-1);
  map = get_map();
  danot_(no_tps);
  MNF = MapNorm(map, no_tps);

  CtoR(MNF.K*Id_scl, k_re, k_im);

  k_buf.push_back(k_abs_ijklm(k_re, k_im, 1, 1, 0, 0, 0));
  k_buf.push_back(k_abs_ijklm(k_re, k_im, 1, 1, 1, 1, 0));
  k_buf.push_back(k_abs_ijklm(k_re, k_im, 0, 0, 1, 1, 0));

  k_buf.push_back(k_abs_ijklm(k_re, k_im, 1, 1, 0, 0, 2));
  k_buf.push_back(k_abs_ijklm(k_re, k_im, 0, 0, 1, 1, 2));

  return k_buf;
}


double get_chi_2(const std::vector<double> &h_ijklm)
{
  // Compute RMS for the on-momentum tune tune footprint terms.
  double chi_2 = 0e0;

  for (int k = 0; k < 3; k++)
    chi_2 += sqr(h_ijklm[k]);
  return chi_2;
}

void prt_k_ijklm
(ostream &outf, const double nu[], const double dnu[],
 const std::vector<double> &h_ijklm)
{
  outf << fixed << setprecision(5)
       << setw(9) << nu[X_] << setw(9) << nu[Y_]
       << setw(9) << dnu[X_] << setw(9) << dnu[Y_];
  for (int k = 0; k < 3; k++)
    outf << scientific << setprecision(5) << setw(13) << h_ijklm[k];
  outf << scientific << setprecision(5)
       << setw(13) << sqrt(get_chi_2(h_ijklm));
  for (int k = 3; k < 5; k++)
    outf << scientific << setprecision(5) << setw(13) << h_ijklm[k];
  outf << "\n";
}


void scan_nu
(const string &fam_name, const ss_vect<tps> &Id_scl,
 const int n_step, const double nu0_x, const double nu0_y,
 const double dnu_offset_x, const double dnu_offset_y,
 const double dnu_max_x, const double dnu_max_y)
{
  const double
    nu0[]         = {nu0_x, nu0_y},
    dnu_offset[] = {dnu_offset_x, dnu_offset_y},
    dnu_max[]    = {dnu_max_x, dnu_max_y};
  const string
    file_name = "k_ijklm.dat"; 

  double              nu[2], dnu[2];
  std::vector<double> k_ijklm;
  ofstream            outf;

  file_wr(outf, file_name.c_str());

  outf << "\n#   nu_x     nu_y    dnu_x    dnu_y     k_22000      k_11110"
       << "      k_00220       rms      k_11002      k_00112\n";
  cout << "\n";
  for (int j = 0; j <= n_step; j++) {
    cout << "  j = " << j << "\n";
    dnu[X_] = j*dnu_max[X_]/n_step + dnu_offset[X_];
    nu[X_] = nu0[X_] + dnu[X_];
    for (int k = 0; k <= n_step; k++) {
      dnu[Y_] = k*dnu_max[Y_]/n_step + dnu_offset[Y_];
      nu[Y_] = nu0[Y_] + dnu[Y_];
      if ((j !=0) && (k == 0)) {
	outf << "\n";
      }
      // A 1/2 ps_rot at the entrance & exit of the super period for a symmetric
      // approach.
      set_ps_rot(fam_name, dnu[X_]/2e0, dnu[Y_]/2e0);
      k_ijklm = get_h_ijklm(Id_scl);

      prt_k_ijklm(outf, nu, dnu, k_ijklm);
    }
  }

  outf.close();
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

  daeps_(1e-30);

  if (!true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  // Twiss functions required by ps_rot.
  Ring_GetTwiss(true, 0e0);
  printglob();

  Id_scl.identity();
  for (int k = 0; k < 2; k++) {
    twoJ = sqr(A_max[k])/Cell[0].Beta[k];
    Id_scl[2*k] *= sqrt(twoJ);
    Id_scl[2*k+1] *= sqrt(twoJ);
  }
  Id_scl[delta_] *= delta_max;

  if (false)
    scan_nu("ps_rot1", Id_scl, 50, 0.290, 0.180, -0.15, -0.15, 0.3, 0.3);

  if (false)
    scan_nu("ps_rot2", Id_scl, 50, 2.274, 1.403, -0.5, -0.5, 1.0, 1.0);

  if (!false) {
    // A 1/2 ps_rot at the entrance & exit of the super period for a symmetric
    // approach.
    set_ps_rot("ps_rot2", -0.15/2e0, -0.1/2e0);

    prtmfile("flat_file.fit");
  }
}
