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
    eta1[]     = { 0.0177041,   0.0},
    etap1[]    = { 0.0566527,   0.0},
    alpha1[]   = {-6.45968,     2.46939},
    beta1[]    = { 2.73766,     1.32365},
    m_etap1[]  = {-etap1[X_],  -etap1[Y_]},
    m_alpha1[] = {-alpha1[X_], -alpha1[Y_]};
  
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

  if (true) {
    Ring_GetTwiss(true, 0e0); printglob();
    GetEmittance(ElemIndex("cav"), true);
    Ring_GetTwiss(true, 0e0);
  } else {
    // globval.emittance = true;
    // // A = get_A(m_alpha1, beta1, eta1, m_etap1);
    // A = get_A(m_alpha1, beta1, eta1, etap1);
    // Cell_Twiss(0, globval.Cell_nLoc, A, false, false, 0e0);
    // eps_x = get_eps_x1(false);
    // globval.emittance = false;

    ttwiss(m_alpha1, beta1, eta1, m_etap1, 0e0);
  }

  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);
  prtmfile("flat_file.dat");
}
