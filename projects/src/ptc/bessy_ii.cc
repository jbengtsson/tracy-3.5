#define NO 5

#include "tracy_lib.h"

int  no_tps   = NO,
     ndpt_tps = 5;


void get_alpha_c(std::vector<double> &alpha_c)
{
  for (int k = 0; k < NO-1; k++)
    alpha_c.push_back
      (h_ijklm(map[ct_], 0, 0, 0, 0, k+1)/Cell[globval.Cell_nLoc].S);
}


void prt_alpha_c(std::vector<double> &alpha_c)
{
  double po2, q, pm;

  cout << "\nalpha_c = ";
  for (int k = 0; k < alpha_c.size(); k++) {
    if (k == 0)
      cout << scientific << setprecision(3)
	   << setw(9) << alpha_c[k];
    else if (k == 1)
      cout << scientific << setprecision(3)
	   << setw(9) << alpha_c[k] << "*delta";
    else if (k > 1)
      cout << scientific << setprecision(3)
	   << setw(9) << alpha_c[k] << "*delta^" << k;
    if (k < alpha_c.size()-1)
      cout << " + ";
  }
  

  if (NO >= 3) {
    cout << fixed << setprecision(2)
	 << "\nFixed points to O(3): " << 0
	 << ", " << -1e2*alpha_c[0]/alpha_c[1] <<"%\n";
  }

  if (NO >= 4) {
    po2 = alpha_c[1]/(2.0*alpha_c[2]); q = alpha_c[0]/alpha_c[2];
    pm = sqrt(sqr(po2)-q);
    cout << fixed << setprecision(2)
	 << "\nFixed points to O(4): " << 0
	 << ", " << -1e2*(po2+pm) << "%, "
	 << -1e2*(po2-pm) << "%\n";
  }
}


double get_f_RF(const int Fnum)
{
  const int loc = Elem_GetPos(Fnum, 1);
  return Cell[loc].Elem.C->f_RF;
}


void set_f_RF(const int Fnum, const double f_RF)
{
  const int loc = Elem_GetPos(Fnum, 1);
  Cell[loc].Elem.C->f_RF = f_RF;
}


void compute_f_0(const int Fnum, const std::vector<double> &alpha_c)
{
  // Validate estimated:
  //   df_RF = -(alpha_c[0]*delta+alpha_c[1]*sqr(delta)+...)*f_RF
  const int
    n[]       = {-15, 10},
    cod_n_max = 10;
  const double
    cod_eps   = 1e-10,
    f_RF_step = 1e3,
    C_0       = Cell[globval.Cell_nLoc].S;

  long int lastpos;
  double   f_RF, df_RF, dct, delta, df_RF_mod;

  trace = false;

  globval.Cavity_on = globval.radiation = globval.pathlength = true;

  danot_(1);

  f_RF = get_f_RF(Fnum);

  printf("\n");
  for (int k = n[0]; k <= n[1]; k++) {
    df_RF = k*f_RF_step;
    set_f_RF(Fnum, f_RF+df_RF);

    GetCOD(cod_n_max, cod_eps, 0e0, lastpos);
    dct = globval.CODvect[ct_];
    delta = globval.CODvect[delta_];

    df_RF_mod = 0e0;
    for (int k = 0; k < alpha_c.size(); k++)
      df_RF_mod -= alpha_c[0]*pow(delta, k+1);
    df_RF_mod *= f_RF;
    printf("  %5.1f %5.1f %8.5f\n", 1e-3*df_RF, 1e-3*df_RF_mod, dct);
    fflush(stdout);
  }
}


void set_lat_state(void)
{
  globval.H_exact     = false;
  globval.quad_fringe = false;
  globval.Cavity_on   = false;
  globval.radiation   = false;
  globval.emittance   = false;
  globval.IBS         = false;
  globval.pathlength  = false;
  globval.bpm         = 0;
}


int main(int argc, char *argv[])
{
  std::vector<double> alpha_c;

  set_lat_state();

  // disable from TPSALib- and LieLib log messages
  idprset_(-1);

  if (!true)
    Read_Lattice(argv[1]);
  else {
    rdmfile(argv[1]);
  }

  if (false) {
    danot_(1);
    Ring_GetTwiss(true, 0e0);
    printglob();
  }

  danot_(no_tps-1);
  get_map(false);

  get_alpha_c(alpha_c);
  prt_alpha_c(alpha_c);

  if (!false)
    compute_f_0(ElemIndex("cavh1t8r"), alpha_c);
}
