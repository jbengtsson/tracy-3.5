
#define NO 4

#include "tracy_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;


double       twoJ[2];
ss_vect<tps> Id_scl;

const double beta_inj[] = {10.0, 3.0}, A_max[] = {5e-3,   2e-3},
                          delta_max = 3e-2;


void prt_eta(FILE *fp, const int order, const int k, const ss_vect<tps> &cod)
{
  int i, j;

  fprintf(fp, "%4d %10s %8.3f %4.1f",
	  k, Cell[k].Elem.PName, Cell[k].S, get_code(Cell[k]));
  for (i = 1; i <= order; i++)
    for (j = 0; j < 2; j++)
      fprintf(fp, " %13.5e", h_ijklm(cod[j], 0, 0, 0, 0, i));
    fprintf(fp, "\n");
}


ss_vect<tps> get_A_nl(const tps g)
{
  int          j;
  tps          gn;
  ss_vect<tps> Id, A_nl;

  Id.identity(); A_nl = Id;
  for (j = 3; j <= no_tps; j++) {
    gn = Take(g, j); A_nl = A_nl*LieExp(gn, Id);
  }
  return A_nl;
}


void prt_dbeta_deta(FILE *fp, const int k, const ss_vect<tps> &A)
{
  ss_vect<tps> A_Atp;


  A_Atp = A*tp_S(3, A);

  fprintf(fp, "%4d %10s %8.3f %4.1f %10.5f %10.5f %8.5f %8.5f"
	  " %10.5f %10.5f %8.5f %8.5f",
	  k, Cell[k].Elem.PName, Cell[k].S, get_code(Cell[k]),
	  h_ijklm(A_Atp[x_]*Id_scl, 1, 0, 0, 0, 0),
	  h_ijklm(A_Atp[y_]*Id_scl, 0, 0, 1, 0, 0),
	  h_ijklm(A[x_], 0, 0, 0, 0, 1),
	  h_ijklm(A[px_], 0, 0, 0, 0, 1),
	  h_ijklm(A_Atp[x_]*Id_scl, 1, 0, 0, 0, 0)
	  +h_ijklm(A_Atp[x_]*Id_scl, 1, 0, 0, 0, 1),
	  h_ijklm(A_Atp[y_]*Id_scl, 0, 0, 1, 0, 1)
	  + h_ijklm(A_Atp[y_]*Id_scl, 0, 0, 1, 0, 0),
	  h_ijklm(A[x_], 0, 0, 0, 0, 1)
	  +h_ijklm(A[x_]*delta_max, 0, 0, 0, 0, 2),
	  h_ijklm(A[px_], 0, 0, 0, 0, 1)
	  +h_ijklm(A[px_]*delta_max, 0, 0, 0, 0, 2));
  fprintf(fp, "\n");
}


void get_dbeta_eta()
{
  int          k;
  ss_vect<tps> A;
  FILE         *outf;

  const string file_name = "dbeta_deta.out";

  get_map(false); MNF = MapNorm(map, 1);
  // A = MNF.A0*MNF.A1*get_A_nl(MNF.g);
  A = LieExp(MNF.g, MNF.A0*MNF.A1);

  outf = file_write(file_name.c_str());
  for (k = 0; k <= globval.Cell_nLoc; k++) {
    Elem_Pass(k, A);
    prt_dbeta_deta(outf, k, A);
  }
  fclose(outf);

  // std::cout << std::scientific << std::setprecision(3)
  // 	    << MNF.A0 << "\n";
}


int main(int argc, char *argv[])
{
  int j;

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  // disable from TPSALib- and LieLib log messages
  idprset_(-1);

  if (false)
    Read_Lattice(argv[1]);
  else {
    rdmfile(argv[1]);
  }

  for (j = 0; j < 2; j++)
    twoJ[j] = sqr(A_max[j])/beta_inj[j];
  Id_scl.identity();
  // for (j = 0; j < 4; j++)
  //   Id_scl[j] *= sqrt(twoJ[j/2]);
  Id_scl[delta_] *= delta_max;

  danot_(2);

  Ring_GetTwiss(true, 0.0); printglob();

  if (!false) {
    set_map("ps_rot", 0.12, 0.08);
    Ring_GetTwiss(true, 0e0); printglob();
  }

  danot_(no_tps);

  get_dbeta_eta();
}
