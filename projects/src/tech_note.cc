#define NO 1

#include "tracy_lib.h"

int no_tps = NO;

double get_eps_x1(const bool track)
{
  // eps_x [nm.rad].
  long int     lastpos;
  double       eps_x;
  ss_vect<tps> A;

  const bool prt = !false;

  if (track) {
    globval.emittance = true;
    // A = get_A(ic[0], ic[1], ic[2], ic[3]);
    putlinmat(6, globval.Ascr, A);
    Cell_Pass(0, globval.Cell_nLoc, A, lastpos);
    // Cell_Twiss(0, globval.Cell_nLoc, A, false, false, 0e0);
    globval.emittance = false;
  }

  eps_x = 1470e0*sqr(globval.Energy)*I5/(I2-I4);

  if (prt) {
    printf("\neps_x = %5.3f nm.rad\n", eps_x);
    printf("J_x   = %5.3f, J_z = %5.3f\n", 1.0-I4/I2, 2.0+I4/I2);
  }

  return eps_x;
}


int main(int argc, char *argv[])
{
  double       eps_x;
  ss_vect<tps> A;

  const double
    eta0[]     = { 0.00152705,  0.0},
    etap0[]    = { 0.0,         0.0},
    alpha0[]   = { 0.0,         0.0},
    beta0[]    = { 0.0962146,   0.2},
    eta1[]     = { 0.00609655,  0.0},
    etap1[]    = { 0.0252206,   0.0},
    alpha1[]   = {-5.22971,    -1.22356},
    beta1[]    = { 1.82155,     0.73649},
    m_etap1[]  = {-0.0252206,   0.0},
    m_alpha1[] = { 5.22971,     1.22356};
  
  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  globval.H_exact        = false;  globval.quad_fringe = false;
  globval.Cavity_on      = false;  globval.radiation   = false;
  globval.emittance      = false;  globval.IBS         = false;
  globval.pathlength     = false;  globval.bpm         = 0;
  globval.dip_edge_fudge = true;

  if (!false) set_map_per(ElemIndex("Mir"), alpha1, beta1, eta1, etap1);

  if (!true) {
    Ring_GetTwiss(true, 0e0); printglob();
    // trace = true;
    // GetEmittance(ElemIndex("cav"), true);
  } else {
    globval.emittance = true;
    A = get_A(m_alpha1, beta1, eta1, m_etap1);
    Cell_Twiss(0, globval.Cell_nLoc, A, false, false, 0e0);
    eps_x = get_eps_x1(false);
    globval.emittance = false;
  }

  ttwiss(m_alpha1, beta1, eta1, m_etap1, 0e0);

  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);
  prtmfile("flat_file.dat");
}
