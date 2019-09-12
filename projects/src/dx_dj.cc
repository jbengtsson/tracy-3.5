#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


void dx_dJ(const int n_turn, const double Ax, const double Ay)
{
  long int                      lastpos;
  int                           j, k, n;
  ss_vect<double>               ps;
  std::vector<int>              loc;
  std::vector<ss_vect<double> > ps_loc;
  ofstream                      outf;

  file_wr(outf, "dx_dJ.out");

  ps.zero();
  printf("\n");
  for (k = 0; k <= globval.Cell_nLoc; k++)
    if (!false || ((Cell[k].Elem.Pkind == Mpole)
	 && (Cell[k].Elem.M->PBpar[Sext+HOMmax] != 0e0))) {
      loc.push_back(k); ps_loc.push_back(ps);
    }

  n = loc.size();

  for (k = 0; k < n; k++)
    printf("  %3d %7.3f %10s\n",
	   loc[k], Cell[loc[k]].S, Cell[loc[k]].Elem.PName);

  globval.Cavity_on = false; globval.radiation = false;

  ps.zero();
  ps[x_] = Ax; ps[y_] = Ay;
  for (j = 1; j <= n_turn; j++) {
    Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
    for (k = 0; k < n; k++)
      ps_loc[k] += Cell[loc[k]].BeamPos;
  }

  for (k = 0; k < n; k++)
    for (j = 0; j < ss_dim; j++)
#if 1
      ps_loc[k][j] /= n_turn;
#else
    ps_loc[k][j] = fabs(ps_loc[k][j]/n_turn);
#endif

  for (k = 0; k < n; k++)
    outf << setw(4) << k
	 << fixed << setprecision(5) << setw(9) << Cell[loc[k]].S
	 << " " << setw(10) << Cell[loc[k]].Elem.PName
	 << scientific << setprecision(5) << setw(13) << ps_loc[k] << "\n";

  outf.close();
}


int main(int argc, char *argv[])
{

  reverse_elem = !false;

  trace = !true;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  globval.H_exact    = false; globval.quad_fringe    = false;
  globval.Cavity_on  = false; globval.radiation      = false;
  globval.emittance  = false; globval.IBS            = false;
  globval.pathlength = false; globval.bpm            = 0;
  globval.Cart_Bend  = false; globval.dip_edge_fudge = true;

  globval.Cavity_on = false; globval.radiation = false;
  Ring_GetTwiss(true, 0e0); printglob();

  dx_dJ(500, 1e-3, 0e-3);

}
