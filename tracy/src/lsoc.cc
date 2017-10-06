/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

*/


static bool              first_h[] = {true, true}, first_v[] = {true, true};
       int               n_bpm_[2], n_corr_[2];
       long unsigned int *bpms_[2], *corrs_[2];
static double            *w_lsoc[2], **A_lsoc[2], **U_lsoc[2], **V_lsoc[2];
static double            *w_lstc[2], **A_lstc[2], **U_lstc[2], **V_lstc[2];


void zero_trims(void)
{
  int      j, k;
  long int loc;

  for (k = 0; k < 2; k++)
    for (j = 1; j <= n_corr_[k]; j++) {
      loc = corrs_[k][j];
      set_bn_design_elem(Lattice.Cell[loc].Fnum, Lattice.Cell[loc].Knum, Dip, 0e0, 0e0);
    }
}


void prt_gcmat(const int plane)
{
  int  i, j, k;
  FILE *outf = NULL;

  k = plane - 1;

  printf("\n");
  printf("no of bpms = %d, no of corrs = %d, plane = %d\n",
	 n_bpm_[k], n_corr_[k], plane);

  if (plane == 1)
    outf = file_write("svdh.out");
  else if (plane == 2)
    outf = file_write("svdv.out");
  else {
    printf("prt_gcmat: undefined plane %d\n", plane);
    exit_(1);
  }

  fprintf(outf,"# total no of monitors:                %d\n", n_bpm_[k]);

  if (plane == 1)
    fprintf(outf,"# total no of horizontal correctors: %d\n", n_corr_[k]);
  else
    fprintf(outf,"# total no of vertical correctors:   %d\n", n_corr_[k]);

  fprintf(outf, "# A[%d][%d] = \n", n_bpm_[k], n_corr_[k]);
  for (i = 1; i <= n_bpm_[k]; i++) {
    for (j = 1; j <= n_corr_[k]; j++)
      fprintf(outf, "% .13e ", A_lsoc[k][i][j]);
    fprintf(outf, "\n");
  }

  fprintf(outf, "# U[%d][%d] = \n", n_bpm_[k], n_corr_[k]);
  for (i = 1; i <= n_bpm_[k]; i++) {
    for (j = 1; j <= n_corr_[k]; j++)
      fprintf(outf, "% .13e ", U_lsoc[k][i][j]);
    fprintf(outf, "\n");
  }

  fprintf(outf, "# w[%d]    = \n", n_corr_[k]);
  for (j = 1; j <= n_corr_[k]; j++)
    fprintf(outf, "% .13e ", w_lsoc[k][j]);
  fprintf(outf, "\n");
  fprintf(outf, "# V[%d][%d] = \n", n_bpm_[k], n_corr_[k]);

  for (i = 1; i <= n_corr_[k]; i++) {
    for (j = 1; j <= n_corr_[k]; j++)
      fprintf(outf, "% .13e ", V_lsoc[k][i][j]);
    fprintf(outf, "\n");
  }

  fclose(outf);
}


void gcmat(const int plane)
{
  /* Get orbit response matrix

                -----------
              \/beta  beta
                    i     j
        A   = ------------- cos(nu pi - 2 pi|nu  - nu |)
         ij   2 sin(pi nu)                     i     j

  */

  int      i, j, k;
  long int loc;
  double   nu, betai, betaj, nui, nuj, spiq;

  const double eps = 1e-4;

  k = plane - 1;

  nu = globval.TotalTune[k]; spiq = sin(M_PI*nu);

  for (i = 1; i <= n_bpm_[k]; i++) {
    loc = bpms_[k][i]; betai = Lattice.Cell[loc].Beta[k]; nui = Lattice.Cell[loc].Nu[k];
    for (j = 1; j <= n_corr_[k]; j++) {
      loc = corrs_[k][j]; betaj = Lattice.Cell[loc].Beta[k]; nuj = Lattice.Cell[loc].Nu[k];
      A_lsoc[k][i][j] =
	sqrt(betai*betaj)/(2.0*spiq)*cos(nu*M_PI-fabs(2.0*M_PI*(nui-nuj)));
    }
  }

  for (i = 1; i <= n_bpm_[k]; i++)
    for (j = 1; j <= n_corr_[k]; j++)
      U_lsoc[k][i][j] = A_lsoc[k][i][j];

  dsvdcmp(U_lsoc[k], n_bpm_[k], n_corr_[k], w_lsoc[k], V_lsoc[k]);

  printf("\n");
  printf("gcmat singular values:\n");
  for (j = 1; j <= n_corr_[k]; j++) {
    printf("%11.3e", w_lsoc[k][j]);
    if (w_lsoc[k][j] < eps) {
      w_lsoc[k][j] = 0e0;
      printf(" (zeroed)");
    }
    if (j % 5 == 0) printf("\n");
  }
  if (n_corr_[k] % 5 != 0) printf("\n");

  if (trace) prt_gcmat(plane);
}


void gcmat(const int n_bpm, const long int bpms[],
	   const int n_corr, const long int corrs[], const int plane,
	   const bool svd)
{
  bool first;
  int  i, k;

  k = plane - 1;

  first = (plane == 1)? first_h[0] : first_v[0];
  if (first) {
    if (plane == 1)
      first_h[0] = false;
    else
      first_v[0] = false;

    bpms_[k] = lvector(1, n_bpm); corrs_[k] = lvector(1, n_corr);

    A_lsoc[k] = dmatrix(1, n_bpm, 1, n_corr);
    U_lsoc[k] = dmatrix(1, n_bpm, 1, n_corr);
    w_lsoc[k] = dvector(1, n_corr);
    V_lsoc[k] = dmatrix(1, n_corr, 1, n_corr);
  }

  for (i = 1; i <= n_bpm; i++)
    bpms_[k][i] = bpms[i-1];

  for (i = 1; i <= n_corr; i++)
    corrs_[k][i] = corrs[i-1];

  n_bpm_[k] = n_bpm; n_corr_[k] = n_corr;

  if (svd) gcmat(plane);
}


void gcmat(const int bpm, const int corr, const int plane)
{
  int i, k;

  k = plane - 1; n_bpm_[k] = Lattice.GetnKid(bpm);
  n_corr_[k] = Lattice.GetnKid(corr);

  long int bpms[n_bpm_[k]], corrs[n_corr_[k]];

  for (i = 1; i <= n_bpm_[k]; i++)
    bpms[i-1] = Lattice.Elem_GetPos(bpm, i);

  for (i = 1; i <= n_corr_[k]; i++)
    corrs[i-1] = Lattice.Elem_GetPos(corr, i);

  gcmat(n_bpm_[k], bpms, n_corr_[k], corrs, plane, true);
}


void lsoc(const int plane, const double scl)
{
  int      j, k;
  long int loc;
  double   *b, *x;

  k = plane - 1;

  b = dvector(1, n_bpm_[k]); x = dvector(1, n_corr_[k]);

  for (j = 1; j <= n_bpm_[k]; j++) {
    loc = bpms_[k][j];
    b[j] = -Lattice.Cell[loc].BeamPos[2*k] + Lattice.Cell[loc].dS[k];
  }
      
  dsvbksb(U_lsoc[k], w_lsoc[k], V_lsoc[k], n_bpm_[k], n_corr_[k], b, x);

  for (j = 1; j <= n_corr_[k]; j++) {
    loc = corrs_[k][j];
    if (plane == 1)
      set_dbnL_design_elem(Lattice.Cell[loc].Fnum, Lattice.Cell[loc].Knum, Dip, -scl*x[j], 0e0);
    else
      set_dbnL_design_elem(Lattice.Cell[loc].Fnum, Lattice.Cell[loc].Knum, Dip, 0e0, scl*x[j]);
  }

  free_dvector(b, 1, n_bpm_[k]); free_dvector(x, 1, n_corr_[k]);
}


void gtcmat(const int plane)
{
  /* Get trajectory response matrix

                -----------
        A   = \/beta  beta  sin(2 pi(nu  - nu ))
         ij         i     j            i     j

  */

  int      i, j, k;
  long int loc_bpm, loc_corr;
  double   betai, betaj, nui, nuj;

  const double eps = 1e-4;

  k = plane - 1;

  for (i = 1; i <= n_bpm_[k]; i++) {
    loc_bpm = bpms_[k][i];
    betai = Lattice.Cell[loc_bpm].Beta[k]; nui = Lattice.Cell[loc_bpm].Nu[k];
    for (j = 1; j <= n_corr_[k]; j++) {
      loc_corr = corrs_[k][j];
      betaj = Lattice.Cell[loc_corr].Beta[k]; nuj = Lattice.Cell[loc_corr].Nu[k];
      if (loc_bpm > loc_corr)
	A_lstc[k][i][j] = sqrt(betai*betaj)*sin(2.0*M_PI*(nui-nuj));
      else
	A_lstc[k][i][j] = 0e0;
    }
  }

  for (i = 1; i <= n_bpm_[k]; i++)
    for (j = 1; j <= n_corr_[k]; j++)
      U_lstc[k][i][j] = A_lstc[k][i][j];

  dsvdcmp(U_lstc[k], n_bpm_[k], n_corr_[k], w_lstc[k], V_lstc[k]);

  printf("\n");
  printf("gtcmat singular values:\n");
  for (j = 1; j <= n_corr_[k]; j++) {
    printf("%11.3e", w_lstc[k][j]);
    if (w_lstc[k][j] < eps) {
      w_lstc[k][j] = 0e0;
      printf(" (zeroed)");
    }
    if (j % 5 == 0) printf("\n");
  }
  if (n_corr_[k] % 5 != 0) printf("\n");

  if (trace) prt_gcmat(plane);
}


void gtcmat(const int n_bpm, const long int bpms[],
	    const int n_corr, const long int corrs[], const int plane,
	    const bool svd)
{
  bool first;
  int  i, k;

  k = plane - 1;

  first = (plane == 1)? first_h[1] : first_v[1];
  if (first) {
    if (plane == 1)
      first_h[1] = false;
    else
      first_v[1] = false;

    bpms_[k] = lvector(1, n_bpm); corrs_[k] = lvector(1, n_corr);

    A_lstc[k] = dmatrix(1, n_bpm, 1, n_corr);
    U_lstc[k] = dmatrix(1, n_bpm, 1, n_corr);
    w_lstc[k] = dvector(1, n_corr);
    V_lstc[k] = dmatrix(1, n_corr, 1, n_corr);
  }

  for (i = 1; i <= n_bpm; i++)
    bpms_[k][i] = bpms[i-1];

  for (i = 1; i <= n_corr; i++)
    corrs_[k][i] = corrs[i-1];

  n_bpm_[k] = n_bpm; n_corr_[k] = n_corr;

  if (svd) gtcmat(plane);
}


void lstc(const int plane, const double scl)
{
  int      j, k;
  long int loc;
  double   *b, *x;

  k = plane - 1;

  b = dvector(1, n_bpm_[k]); x = dvector(1, n_corr_[k]);

  for (j = 1; j <= n_bpm_[k]; j++) {
    loc = bpms_[k][j];
    b[j] = -Lattice.Cell[loc].BeamPos[2*k] + Lattice.Cell[loc].dS[k];

    if (trace) std::cout << std::scientific << std::setprecision(5)
		    << "b[" << std::setw(3) << j << "] = "
		    << std::setw(12) << b[j] << std::endl;
  }
      
  dsvbksb(U_lstc[k], w_lstc[k], V_lstc[k], n_bpm_[k], n_corr_[k], b, x);

  for (j = 1; j <= n_corr_[k]; j++) {
    loc = corrs_[k][j];
    if (plane == 1) {
      if (trace) std::cout << std::scientific << std::setprecision(5)
		      << "(b_1L)[" << std::setw(3) << j << "] = "
		      << std::setw(12)<< -x[j] << std::endl;

      set_dbnL_design_elem(Lattice.Cell[loc].Fnum, Lattice.Cell[loc].Knum, Dip,
			   -scl*x[j], 0e0);
    } else {
      if (trace) std::cout << std::scientific << std::setprecision(5)
		      << "(a_1L)[" << std::setw(3) << j << "] = "
		      << std::setw(12)<< x[j] << std::endl;

      set_dbnL_design_elem(Lattice.Cell[loc].Fnum, Lattice.Cell[loc].Knum, Dip,
			   0e0, scl*x[j]);
    }
  }

  free_dvector(b, 1, n_bpm_[k]); free_dvector(x, 1, n_corr_[k]);
}
