#define NO 1

#include "tracy_lib.h"

int no_tps = NO;

const int n_cell = 6;

#if 1
const double nu[] = {10.87, 3.29};
#else
const double nu[] = {64.91/n_cell, 19.66/n_cell};
#endif


void prt_name(FILE *outf, const char *name, const string &str, const int len)
{
  int j, k;

  j = 0;
  do {
    fprintf(outf, "%c", name[j]);
    j++;
  } while (name[j] != ' ');
  fprintf(outf, str.c_str());
  for (k = j; k < len; k++)
    fprintf(outf, " ");
}


void prt_b2(FILE *outf, const std::vector<int> &Fnum_b2)
{
  int k, loc;
  elemtype *elemp;

  for (k = 0; k < (int)Fnum_b2.size(); k++) {
    loc = Elem_GetPos(Fnum_b2[k], 1);
    elemp = &Cell[loc].Elem;
    prt_name(outf, elemp->PName, ":", 8);
    fprintf(outf, " quadrupole, l = %10.8f, k = %13.8f, n = nquad"
	    ", method = 4;\n",
	    elemp->PL, elemp->M->PBpar[Quad+HOMmax]);
  }
}


void fit_nu(const std::vector<int> &Fnum_b2,
	    const double nu_x, const double nu_y, const double db2L)
{
  int    n_b2, j, k, n_svd;
  double **A, **U, **V, *w, *b, *x, b2L, a3L;
  FILE   *outf;

  const bool   prt = !false;
  const int    m   = 2;
  const double
    nu0[]  = {nu_x, nu_y},
    svd_cut = 1e-10;
  const string file_name = "fit_nu_b2.out";

  n_b2 = Fnum_b2.size();

  A = dmatrix(1, m, 1, n_b2); U = dmatrix(1, m, 1, n_b2);
  V = dmatrix(1, n_b2, 1, n_b2);
  w = dvector(1, n_b2); b = dvector(1, m); x = dvector(1, n_b2);

  for (k = 1; k <= n_b2; k++) {
    set_dbnL_design_fam(Fnum_b2[k-1], Quad, db2L, 0e0);
    Ring_GetTwiss(true, 0e0);
    if (prt)
      printf("\nfit_nu: nu1+ = [%9.5f, %9.5f]\n",
	     globval.TotalTune[X_], globval.TotalTune[Y_]);

    for (j = 1; j <= m; j++)
      A[j][k] = globval.TotalTune[j-1];
    set_dbnL_design_fam(Fnum_b2[k-1], Quad, -2e0*db2L, 0e0);
    Ring_GetTwiss(true, 0e0);
    if (prt)
      printf("fit_nu: nu1- = [%9.5f, %9.5f]\n",
	 globval.TotalTune[X_], globval.TotalTune[Y_]);
    for (j = 1; j <= 2; j++) {
      A[j][k] -= globval.TotalTune[j-1]; A[j][k] /= 2e0*db2L;
    }

    set_dbnL_design_fam(Fnum_b2[k-1], Quad, db2L, 0e0);
  }

  Ring_GetTwiss(true, 0e0);
  if (prt)
    printf("\nfit_nu: nu1  = [%9.5f, %9.5f]\n",
	   globval.TotalTune[X_], globval.TotalTune[Y_]);
  for (j = 1; j <= 2; j++) {
    b[j] = -(globval.TotalTune[j-1]-nu0[j-1]);
    if (prt) printf(" %9.5f", b[j]);
  }
  if (prt) printf("\n");
    

  dmcopy(A, m, n_b2, U); dsvdcmp(U, m, n_b2, w, V);

  printf("\nfit_nu:\n  singular values:\n");
  n_svd = 0;
  for (j = 1; j <= n_b2; j++) {
    printf("  %9.3e", w[j]);
    if (w[j] < svd_cut) {
      w[j] = 0e0;
      printf(" (zeroed)");
    } else {
      if (n_svd > 2) {
	printf("fit_nu: more than 2 non-zero singular values");
	exit(1);
      }
      n_svd++;
    }
    printf("\n");
  }

  dsvbksb(U, w, V, m, n_b2, b, x);

  if (prt) printf("\n");
  for (k = 1; k <= n_b2; k++) {
    set_dbnL_design_fam(Fnum_b2[k-1], Quad, x[k], 0e0);
    if (prt) printf(" %9.5f", x[k]);
  }
  if (prt) printf("\n");


  if (prt) {
    printf("  b2:\n  ");
    for (k = 0; k < n_b2; k++) {
      get_bn_design_elem(Fnum_b2[k], 1, Quad, b2L, a3L);
      printf(" %9.5f", b2L);
    }
    printf("\n");
  }

  outf = file_write(file_name.c_str());
  prt_b2(outf, Fnum_b2);
  fclose(outf);

  free_dmatrix(A, 1, m, 1, n_b2); free_dmatrix(U, 1, m, 1, n_b2);
  free_dmatrix(V, 1, n_b2, 1, n_b2);
  free_dvector(w, 1, n_b2); free_dvector(b, 1, m); free_dvector(x, 1, n_b2);
}


void fit_nu(const double nu[])
{
  std::vector<int> Fam_b2;

  Fam_b2.push_back(ElemIndex("qf1"));
  Fam_b2.push_back(ElemIndex("qd2"));

  Fam_b2.push_back(ElemIndex("qf1_c1"));
  Fam_b2.push_back(ElemIndex("qd2_c1"));
  Fam_b2.push_back(ElemIndex("quad_add"));

  fit_nu(Fam_b2, nu[X_], nu[Y_], 1e-3);

  Ring_GetTwiss(true, 0e0); printglob();
}


int main(int argc, char *argv[])
{

  reverse_elem = !false;

  trace = false;

  globval.mat_meth = false;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  globval.H_exact    = false; globval.quad_fringe    = false;
  globval.Cavity_on  = false; globval.radiation      = false;
  globval.emittance  = false; globval.IBS            = false;
  globval.pathlength = false; globval.bpm            = 0;
  globval.Cart_Bend  = false; globval.dip_edge_fudge = true;

  Ring_GetTwiss(true, 0e0); printglob();

}
