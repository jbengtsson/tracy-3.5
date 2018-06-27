#define NO 5

#include "tracy_lib.h"

int  no_tps   = NO,
     ndpt_tps = 5;


const int n_alphac = 3;


void get_alphac(double alphac[])
{
  int i;

  for (i = 0; i < n_alphac; i++)
    alphac[i] = h_ijklm(map[ct_], 0, 0, 0, 0, i+1)/Cell[globval.Cell_nLoc].S;
}


double H_long(const double phi, const double delta,
	      const int h_rf, const double V_rf, const double phi0,
	      const double alphac[])
{
  int    i;
  double H;

  H = V_rf/(globval.Energy*1e9)*(cos(phi+phi0)+phi*sin(phi0));
  for (i = 2; i <= n_alphac+1; i++)
    H += 2.0*pi*h_rf*alphac[i-2]*pow(delta, (double)i)/i;
  return H;
}


void prt_H_long(const int n, const double phi_max, const double delta_max,
		const string &cav_name, const double U0, const bool neg_alphac)
{
  int            loc, i, j, h_rf;
  double         V_rf, alphac[n_alphac], phi, delta, H, delta_rf, phi0;
  std::ofstream  os;

  os.open("H_long.dat", std::ios::out);

  get_alphac(alphac);

  loc = Elem_GetPos(ElemIndex(cav_name.c_str()), 1);
  h_rf = Cell[loc].Elem.C->Ph;
  V_rf = Cell[loc].Elem.C->Pvolt;

  phi0 = - fabs(asin(U0/V_rf));
  if (neg_alphac) phi0 += pi;

  delta_rf = sqrt(-V_rf*cos(pi+phi0)*(2e0-(pi-2e0*(pi+phi0))*tan(pi+phi0))
	     /(alphac[0]*pi*h_rf*globval.Energy*1e9));
  std::cout << std::endl << std::fixed << std::setprecision(1)
	    << "U0     = " << 1e-3*U0 << " keV" << std::endl;
  if (!neg_alphac) 
    std::cout << std::fixed << std::setprecision(2)
	      << "phi0 = " << fabs(phi0)*180e0/pi-180e0
	      << " deg" << std::endl;
  else
    std::cout << std::fixed << std::setprecision(2)
	      << "phi0             = 180 - " << fabs(phi0)*180e0/pi-180e0
	      << " deg" << std::endl;

  for (i = -n; i <= n ; i++) {
    for (j = -n; j <= n ; j++) {
      phi = i*phi_max/n; delta = j*delta_max/n;
      H = H_long(phi, delta, h_rf, V_rf, pi+phi0, alphac);
      os << std::fixed
	 << std::setprecision(2)<< std::setw(8) <<  (phi0+phi)*180e0/pi
	 << std::setprecision(5) << std::setw(10) << 1e2*delta
	 << std::scientific << std::setw(13) << H << std::endl;
    }
    os << std::endl;
  }
  os.close();
}


void prt_alphac()
{
  double alphac[n_alphac];

  get_alphac(alphac);
  std::cout << std::endl << std::scientific << std::setprecision(3)
	    << "alphac = " << std::setw(10) << alphac[0] << " + "
	    << std::setw(10) << alphac[1] << "*delta + "
	    << std::setw(10) << alphac[2] << "*delta^2" << std::endl;
}


void prt_h_K(const double twoJx, const double twoJy, const double delta)
{
  tps           K_re, K_im, h_re, h_im;
  ss_vect<tps>  twoJ, Id_scl, nus;
  std::ofstream outf;

  twoJ[X_] = twoJx; twoJ[Y_] = twoJy;

  Id_scl.identity();
  Id_scl[x_] *= sqrt(twoJ[X_]); Id_scl[px_] *= sqrt(twoJ[X_]);
  Id_scl[y_] *= sqrt(twoJ[Y_]); Id_scl[py_] *= sqrt(twoJ[Y_]);
  Id_scl[delta_] *= delta;

  danot_(no_tps);
  MNF = MapNorm(map, no_tps); nus = dHdJ(MNF.K);
  CtoR(MNF.K, K_re, K_im); CtoR(get_h(), h_re, h_im);

  file_wr(outf, "h.out");
  outf << h_re*Id_scl << h_im*Id_scl;
  outf.close();
  file_wr(outf, "K.out");
  outf << K_re*Id_scl << K_im*Id_scl;
  outf.close();
}


int main(int argc, char *argv[])
{
  const double beta[]  = {3.4, 1.9},
               A_max[] = {6e-3, 4e-3}, delta = 5e-2;

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  // disable from TPSALib- and LieLib log messages
  idprset_(-1);

  if (true)
    Read_Lattice(argv[1]);
  else {
    rdmfile(argv[1]);
  }

  danot_(1);
  Ring_GetTwiss(true, 0e0); printglob();

  danot_(no_tps-1);
  get_map(false);
  prt_h_K(sqr(A_max[X_])/beta[X_], sqr(A_max[Y_])/beta[Y_], delta);

  danot_(no_tps-1);
  globval.Cavity_on = true; globval.radiation = true;
  get_map(false);
  prt_H_long(10, M_PI, 10e-2, "cav", -544.7e3, true);
  prt_alphac();
}
