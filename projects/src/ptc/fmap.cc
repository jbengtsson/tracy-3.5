#define NO 8
//#define NO 15

#include "tracy_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;


ss_vect<tps> A_inv;


void get_map_normal_form()
{
  danot_(no_tps);

  MNF = MapNorm(map, no_tps);
}


void get_A(void)
{
  int          j;
  iVector      jj;
  tps          gn;
  ss_vect<tps> Id, A;

  Id.identity(); A = MNF.A1;
  for (j = no_tps; j >= 3; j--) {
    gn = Take(MNF.g, j); A = A*LieExp(gn, Id);
  }
  
  for (j = 0; j < nv_tps; j++)
    jj[j] = (j < 4)? 1 : 0;

  A_inv = PInv(A, jj);
}


void get_twoJ(const ss_vect<double> &ps, double twoJ[])
{   
  int             j;
  ss_vect<double> z;
  ss_vect<tps>    Id;

  z = (A_inv*ps).cst();

  for (j = 0; j < 2; j++)
    twoJ[j] = sqr(z[2*j]) + sqr(z[2*j+1]);
}
  

void fmap(const double Ax, const double Ay)
{
  const int  n_ampl = 10;

  int             i, j;
  double          twoJ[2], nu[2], A[2];
  ss_vect<double> ps;
  ss_vect<tps>    Id, nus;
  std::ofstream   fmap_out;

  file_wr(fmap_out, "fmap_est.dat");

  danot_(no_tps-1);
  get_map(false);
  danot_(no_tps);
  get_map_normal_form(); nus = dHdJ(MNF.K); get_A();

  Id.zero(); ps.zero();
  for (i = -n_ampl; i <= n_ampl; i++) {
    A[X_] = i*Ax/n_ampl; ps[x_] = A[X_];
    for (j = -n_ampl; j <= n_ampl; j++) {
      A[Y_] = j*Ay/n_ampl; ps[y_] = A[Y_]; get_twoJ(ps, twoJ);

      Id[x_] = sqrt(twoJ[X_]); Id[px_] = sqrt(twoJ[X_]);
      Id[y_] = sqrt(twoJ[Y_]); Id[py_] = sqrt(twoJ[Y_]);

      nu[X_] = (nus[3]*Id).cst(); nu[Y_] = (nus[4]*Id).cst();

      fmap_out << std::fixed << std::setprecision(3)
	       << std::setw(8) << 1e3*A[X_] << std::setw(8) << 1e3*A[Y_]
	       << std::scientific << std::setprecision(5)
	       << std::setw(13) << nu[X_] << std::setw(13) << nu[Y_]
	       << std::fixed << std::setprecision(1)
	       << std::setw(4) << 1.0 << std::endl;

      if (j == n_ampl) fmap_out << std::endl;
    }
  }

  fmap_out.close();
}


void fmapdp(const double Ax, const double delta)
{
  const int  n_ampl = 10;

  int             i, j;
  double          twoJ[2], nu[2], A[2], d;
  ss_vect<double> ps;
  ss_vect<tps>    Id, nus;
  std::ofstream   fmapdp_out;

  file_wr(fmapdp_out, "fmapdp_est.dat");

  danot_(no_tps-1);
  get_map(false);
  danot_(no_tps);
  get_map_normal_form(); nus = dHdJ(MNF.K); get_A();

  Id.zero(); ps.zero();
  for (i = -n_ampl; i <= n_ampl; i++) {
    A[X_] = i*Ax/n_ampl; ps[x_] = A[X_];
    for (j = -n_ampl; j <= n_ampl; j++) {
      d = j*delta/n_ampl; ps[delta_] = d; get_twoJ(ps, twoJ);

      Id[x_] = sqrt(twoJ[X_]); Id[px_] = sqrt(twoJ[X_]);
      Id[y_] = sqrt(twoJ[Y_]); Id[py_] = sqrt(twoJ[Y_]);
      Id[delta_] = d;

      nu[X_] = (nus[3]*Id).cst(); nu[Y_] = (nus[4]*Id).cst();

      fmapdp_out << std::fixed << std::setprecision(3)
	       << std::setw(8) << 1e3*A[X_] << std::setw(8) << 1e2*d
	       << std::scientific << std::setprecision(5)
	       << std::setw(13) << nu[X_] << std::setw(13) << nu[Y_]
	       << std::fixed << std::setprecision(1)
	       << std::setw(4) << 1.0 << std::endl;

      if (j == n_ampl) fmapdp_out << std::endl;
    }
  }

  fmapdp_out.close();
}


int main(int argc, char *argv[])
{
  int  k;

  const double A_max[] = {6.0e-3, 6.0e-3}, delta_max = 5e-2;

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  // disable from TPSALib- and LieLib log messages
  idprset_(-1);

  if (false)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  globval.EPU = true;

  danot_(1);

  Ring_GetTwiss(true, 0.0); printglob();

  prt_lat("linlat.out", globval.bpm, true);

  k = atoi(argv[2]);
  if (k == 1) {
    std::cout << "fmap" << std::endl;
    fmap(A_max[X_], A_max[Y_]);
  } else if (k == 2) {
    std::cout << "fmapdp" << std::endl;
    fmapdp(0.1e-3, delta_max);
  } else {
    std::cout << "bad param.: k = " << k << std::endl;
    exit(1);
  }
}
