#define NO 3

#include "tracy_lib.h"

int  no_tps   = NO,
     ndpt_tps = 5;


void get_map_normal_form()
{

  danot_(no_tps);

  MNF = MapNorm(map, no_tps);
}


void get_Ke_Im(const tps &h, tps &h_Ke, tps &h_Im)
{
  int           i, j, k;
  iVector       jj;
  ss_vect<tps>  Id;

  for (k = 0; k < nv_tps; k++)
    jj[k] = 0;

  Id.identity(); h_Ke = 0.0;
  for (i = 0; i <= no_tps; i++) {
    jj[x_] = i; jj[px_] = i;
    for (j = 0; j <= no_tps; j++) {
      jj[y_] = j; jj[py_] = j;
      for (k = 0; k <= no_tps; k++) {
	jj[delta_] = k;
	if (2*i+2*j+k <= no_tps) {
	  h_Ke +=
	    h[jj]*pow(Id[x_], i)*pow(Id[px_], i)
	    *pow(Id[y_], j)*pow(Id[py_], j)*pow(Id[delta_], k);
	}
      }
    }
  }

  h_Im = h - h_Ke;
}


tps renorm_g(const double nu0_x, const double nu0_y,
	     const double nu1_x, const double nu1_y,
	     const tps &g)
{
  int           i, j, k, l, m;
  iVector       jj;
  double        scl;
  tps           g_re, g_im, G_re, G_im;
  ss_vect<tps>  Id;

  CtoR(g, g_re, g_im);

  for (k = 0; k < nv_tps; k++)
    jj[k] = 0;

  Id.identity(); G_re = 0.0; G_im = 0.0;
  for (i = 0; i <= no_tps; i++) {
    jj[x_] = i;
    for (j = i; j <= no_tps; j++) {
      jj[px_] = j;
      for (k = 0; k <= no_tps; k++) {
	jj[y_] = k;
	for (l = k; l <= no_tps; l++) {
	  jj[py_] = l;

	  scl =
	    (1.0-1.0/tan(((i-j)*nu0_x+(k-l)*nu0_y)*M_PI))
	    /(1.0-1.0/tan(((i-j)*nu1_x+(k-l)*nu1_y)*M_PI));

	  for (m = 0; m <= no_tps; m++) {
	    jj[delta_] = m;
	    if ((i != j) || (k != l)) {
	      G_re +=
		g_re[jj]*scl*(pow(Id[x_], i)*pow(Id[px_], j)
			   *pow(Id[y_], k)*pow(Id[py_], l)
			   +pow(Id[x_], j)*pow(Id[px_], i)
			   *pow(Id[y_], l)*pow(Id[py_], k))
		*pow(Id[delta_], m);

	      G_im +=
		g_im[jj]*scl*(pow(Id[x_], i)*pow(Id[px_], j)
			   *pow(Id[y_], k)*pow(Id[py_], l)
			   +pow(Id[x_], j)*pow(Id[px_], i)
			   *pow(Id[y_], l)*pow(Id[py_], k))
		*pow(Id[delta_], m);
	    }
	  }
	}
      }
    }
  }

  cout << g_re << G_re;

  return RtoC(G_re, G_im);
}


void tst_g(void)
{
  tps           G;
  ss_vect<tps>  nus;

  danot_(no_tps-1);
  get_map();
  danot_(no_tps);
  get_map_normal_form(); nus = dHdJ(MNF.K);

  cout << endl;
  cout << fixed << setprecision(5)
       << "nu_x = " << nus[0].cst() << ", nu_y = " << nus[1].cst() << endl;

  G = renorm_g(nus[0].cst(), nus[1].cst(), nus[0].cst(), nus[1].cst(), MNF.g);

//  cout << MNF.g << G;
}


int main(int argc, char *argv[])
{
  // disable from TPSALib- and LieLib log messages
  idprset_(-1);

  if (false)
    Read_Lattice(argv[1]);
  else
    rdmfile("flat_file.dat");

  danot_(1);

  Ring_GetTwiss(true, 0.0); printglob();

  tst_g();
}
