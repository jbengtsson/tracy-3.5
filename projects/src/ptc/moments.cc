#define NO 5

#include "tracy_lib.h"


// F. Klein's Erlangen Program.


int no_tps   = NO,
    ndpt_tps = 5;


void compute_C_S_long(double &alpha_z, double &beta_z)
{
  alpha_z =
    -globval.Ascr[ct_][ct_]*globval.Ascr[delta_][ct_]
    - globval.Ascr[ct_][delta_]*globval.Ascr[delta_][delta_];
  beta_z = sqr(globval.Ascr[ct_][ct_]) + sqr(globval.Ascr[ct_][delta_]);
}


ss_vect<tps> compute_A_A_t(void)
{
  int          k;
  double       alpha_z, beta_z;
  ss_vect<tps> Id, A_A_t;

  Id.identity();

  if (nd_tps > 2) compute_C_S_long(alpha_z, beta_z);

  A_A_t.zero();
  for (k = 0; k < nd_tps; k++) {
    if (k < 2) {
      A_A_t[2*k] = Cell[0].Beta[k]*Id[2*k] - Cell[0].Alpha[k]*Id[2*k+1];
      A_A_t[2*k+1] =
	-Cell[0].Alpha[k]*Id[2*k]
	+ (1e0+sqr(Cell[0].Alpha[k]))/Cell[0].Beta[k]*Id[2*k+1];
    } else {
      A_A_t[ct_] = beta_z*Id[ct_] - alpha_z*Id[delta_];
      A_A_t[delta_] = -alpha_z*Id[ct_]	+ (1e0+sqr(alpha_z))/beta_z*Id[delta_];
    }
  }

  return A_A_t;
}


void tst_moment(void)
{
  long int     lastpos;
  ss_vect<tps> M, A_A_t;

  const int
    nd_tps = 3,
#if 0
    loc   = globval.Cell_nLoc;
#else
    loc   = 10;
#endif

  globval.Cavity_on = !false;
  globval.radiation = !false;

  danot_(1);

  Ring_GetTwiss(true, 0e0);
  printglob();

  danot_(no_tps-1);

  M.identity();
  Cell_Pass(0, loc, M, lastpos);
  printf("\nM:");
  prt_lin_map(nd_tps, M);

  A_A_t = compute_A_A_t();
  printf("\nInitial A_A_t beta = [%5.3f, %5.3f]:",
	 A_A_t[x_][x_], A_A_t[y_][y_]);
  prt_lin_map(nd_tps, A_A_t);

  A_A_t = M*tp_S(2, M*A_A_t);

  printf("\nFinal A_A_t beta = [%5.3f, %5.3f]:",
	 A_A_t[x_][x_], A_A_t[y_][y_]);
  prt_lin_map(nd_tps, A_A_t);

  A_A_t = compute_A_A_t();
  printf("\nInitial A_A_t beta = [%5.3f, %5.3f]:",
	 A_A_t[x_][x_], A_A_t[y_][y_]);
  prt_lin_map(nd_tps, A_A_t);

  Cell_Pass(0, loc, A_A_t, lastpos);
  A_A_t = tp_S(2, A_A_t);

  Cell_Pass(0, loc, A_A_t, lastpos);
  printf("\nFinal A_A_t beta = [%5.3f, %5.3f]:",
	 A_A_t[x_][x_], A_A_t[y_][y_]);
  prt_lin_map(nd_tps, A_A_t);
}


tps compute_twoJ(const double eps[], const ss_vect<tps> &A_A_t)
{
  int          j, k;
  tps          twoJ;
  ss_vect<tps> Id, quad_form;

  const ss_vect<tps> omega = get_S(nd_tps);

  Id.identity();

  quad_form = tp_S(nd_tps, omega)*A_A_t*omega;
  twoJ = 0e0;
  for (j = 0; j < 2*nd_tps; j++)
    for (k = 0; k < 2*nd_tps; k++)
      twoJ += Id[j]*sqrt(eps[j/2])*quad_form[j][k]*sqrt(eps[k/2])*Id[k];
  return twoJ;
}


void compute_Twiss(const tps &twoJ, double alpha[], double beta[])
{
  long int jj[ss_dim];
  int      k;

  for (k = 0; k < ss_dim; k++)
    jj[k] = 0;
  for (k = 0; k < 2; k++) {
    jj[2*k]   = 1;
    jj[2*k+1] = 1;
    alpha[k] = twoJ[jj]/2e0;
    jj[2*k]   = 0;
    jj[2*k+1] = 0;
    jj[2*k+1] = 2;
    beta[k] = twoJ[jj];
    jj[2*k+1] = 0;
  }
}


int main(int argc, char *argv[])
{

  globval.H_exact    = false; globval.quad_fringe    = false;
  globval.Cavity_on  = false; globval.radiation      = false;
  globval.emittance  = false; globval.IBS            = false;
  globval.pathlength = false; globval.bpm            = 0;
  globval.Cart_Bend  = false; globval.dip_edge_fudge = true;
  globval.mat_meth   = false;

  // disable from TPSALib- and LieLib log messages
  idprset_(-1);

  if (false)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  if (!false) {
    danot_(1);
    Ring_GetTwiss(true, 0e0);
    printglob();
  }

  tst_moment();
}
