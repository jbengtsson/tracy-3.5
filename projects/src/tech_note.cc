#define NO 1

#include "tracy_lib.h"

int no_tps = NO;

double get_eps_x1(const bool track)
{
  // eps_x [nm.rad].
  long int     lastpos;
  double       I[6], eps_x;
  ss_vect<tps> A;

  const bool   prt     = !false;
  const double C_q_scl = 1e18*C_q/sqr(m_e);

  if (track) {
    globval.emittance = true;
    // A = get_A(ic[0], ic[1], ic[2], ic[3]);
    putlinmat(6, globval.Ascr, A);
    Cell_Pass(0, globval.Cell_nLoc, A, lastpos);
    // Cell_Twiss(0, globval.Cell_nLoc, A, false, false, 0e0);
    globval.emittance = false;
  }

  get_I(I, false);

  eps_x = 1e9*C_q_scl*sqr(globval.Energy)*I[5]/(I[2]-I[4]);

  if (prt) {
    printf("\neps_x = %5.3f nm.rad\n", eps_x);
    printf("J_x   = %5.3f, J_z = %5.3f\n", 1.0-I[4]/I[2], 2.0+I[4]/I[2]);
  }

  return eps_x;
}


void get_map(const double delta)
{
  long int     lastpos;
  ss_vect<tps> map;

  map.identity(); map[delta_] = delta;
  Cell_Pass(0, globval.Cell_nLoc, map, lastpos);

  prt_lin_map(3, map);
}


int main(int argc, char *argv[])
{
  int          j, k;
  double       eps_x, del;
  ss_vect<tps> A;

  const double
#define CASE 5
#if CASE == 1
    eta0[]     = { 0.0078751,  0.0},
    etap0[]    = { 0.0,        0.0},
    alpha0[]   = { 0.0,        0.0},
    beta0[]    = { 0.21821,    5.0},
    eta1[]     = { 0.0328105,  0.0},
    etap1[]    = { 0.063868,   0.0},
    alpha1[]   = {-7.90529,    2.59913},
    beta1[]    = { 5.24392,    1.37611},
    m_etap1[]  = {-etap0[X_],  -etap0[Y_]},
    m_alpha1[] = {-alpha0[X_], -alpha0[Y_]};
#elif CASE == 2
    eta0[]     = {  0.0017490,  0.0},
    etap0[]    = {  0.0,        0.0},
    alpha0[]   = {  0.0,        0.0},
    beta0[]    = {  0.11271,    5.0},
    eta1[]     = {  0.02512,    0.0},
    etap1[]    = {  0.05110,    0.0},
    alpha1[]   = {-14.56328,    2.59911},
    beta1[]    = {  9.36865,    1.37607},
    m_etap1[]  = {-etap0[X_],  -etap0[Y_]},
    m_alpha1[] = {-alpha0[X_], -alpha0[Y_]};
#elif CASE == 3
    eta0[]     = {  0.0017490,  0.0},
    etap0[]    = {  0.0,        0.0},
    alpha0[]   = {  0.0,        0.0},
    beta0[]    = {  3*0.11271,  5.0},
    eta1[]     = {  0.02512,    0.0},
    etap1[]    = {  0.05110,    0.0},
    alpha1[]   = { -5.57908,    2.59911},
    beta1[]    = {  3.88386,    1.37607},
    m_etap1[]  = {-etap0[X_],  -etap0[Y_]},
    m_alpha1[] = {-alpha0[X_], -alpha0[Y_]};
#elif CASE == 4
    eta0[]     = {  0.0017490/2,  0.0},
    etap0[]    = {  0.0,        0.0},
    alpha0[]   = {  0.0,        0.0},
    beta0[]    = {  3*0.11271,  5.0},
    eta1[]     = {  0.02373,    0.0},
    etap1[]    = {  0.04978,    0.0},
    alpha1[]   = { -5.57908,    2.59911},
    beta1[]    = {  3.88386,    1.37607},
    m_etap1[]  = {-etap0[X_],  -etap0[Y_]},
    m_alpha1[] = {-alpha0[X_], -alpha0[Y_]};
#elif CASE == 5
    eta0[]     = { 0.0002943186, 0.0},
    etap0[]    = { 0.0,          0.0},
    alpha0[]   = { 0.0,          0.0},
    beta0[]    = { 0.1153496843, 3.2620137331},
    eta1[]     = { 0.0,          0.0},
    etap1[]    = { 0.0,          0.0},
    alpha1[]   = { 0.0,          0.0},
    beta1[]    = { 1.76416,      2.38834},
    m_etap1[]  = {-etap0[X_],  -etap0[Y_]},
    m_alpha1[] = {-alpha0[X_], -alpha0[Y_]};
#endif

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  globval.H_exact        = false;  globval.quad_fringe = false;
  globval.Cavity_on      = false;  globval.radiation   = false;
  globval.emittance      = false;  globval.IBS         = false;
  globval.pathlength     = false;  globval.bpm         = 0;
  globval.dip_edge_fudge = true;

  if (false) set_map_per(ElemIndex("Mir"), alpha1, beta1, eta1, etap1);

  if (!true) {
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

    ttwiss(alpha0, beta0, eta0, etap0, 0e0);
  }

  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);
  prtmfile("flat_file.dat");
}
