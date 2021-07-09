/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

*/

  
static bool       first_h[] = {true, true}, first_v[] = {true, true};
int               n_bpm_[2], n_corr_[2];
std::vector<int>  bpms_[2], corrs_[2];
static gsl_vector *w_lsoc[2], *w_lstc[2];
static gsl_matrix *A_lsoc[2], *U_lsoc[2], *V_lsoc[2],
                  *A_lstc[2], *U_lstc[2], *V_lstc[2];


void zero_trims(LatticeType &lat)
{
  int      j, k;
  long int loc;

  for (k = 0; k < 2; k++)
    for (j = 0; j < n_corr_[k]; j++) {
      loc = corrs_[k][j];
      set_bn_design_elem(lat, lat.elems[loc]->Fnum, lat.elems[loc]->Knum, Dip,
			 0e0, 0e0);
    }
}


void prt_gcmat(const int plane)
{
  int  i, j;
  FILE *outf = NULL;

  const int k = plane - 1;

  printf("\nno of bpms = %d, no of corrs = %d, plane = %d\n",
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
  for (i = 0; i < n_bpm_[k]; i++) {
    for (j = 0; j < n_corr_[k]; j++)
      fprintf(outf, "% .13e ", gsl_matrix_get(A_lsoc[k], i, j));
    fprintf(outf, "\n");
  }

  fprintf(outf, "# U[%d][%d] = \n", n_bpm_[k], n_corr_[k]);
  for (i = 0; i < n_bpm_[k]; i++) {
    for (j = 0; j < n_corr_[k]; j++)
      fprintf(outf, "% .13e ", gsl_matrix_get(U_lsoc[k], i, j));
    fprintf(outf, "\n");
  }

  fprintf(outf, "# w[%d]    = \n", n_corr_[k]);
  for (j = 0; j < n_corr_[k]; j++)
    fprintf(outf, "% .13e ", gsl_vector_get(w_lsoc[k], j));
  fprintf(outf, "\n# V[%d][%d] = \n", n_bpm_[k], n_corr_[k]);

  for (i = 0; i < n_corr_[k]; i++) {
    for (j = 0; j < n_corr_[k]; j++)
      fprintf(outf, "% .13e ", gsl_matrix_get(V_lsoc[k], i, j));
    fprintf(outf, "\n");
  }

  fclose(outf);
}


void gcmat(LatticeType &lat, const int plane)
{
  /* Get orbit response matrix

                -----------
              \/beta  beta
                    i     j
        A   = ------------- cos(nu pi - 2 pi|nu  - nu |)
         ij   2 sin(pi nu)                     i     j

  */

  int      i, j;
  long int loc;
  double   nu, betai, betaj, nui, nuj, spiq;

  const int    k   = plane - 1;
  const double eps = 1e-4;

  gsl_vector *work = gsl_vector_alloc(n_corr_[k]);

  nu = lat.conf.TotalTune[k]; spiq = sin(M_PI*nu);

  for (i = 0; i < n_bpm_[k]; i++) {
    loc = bpms_[k][i]; betai = lat.elems[loc]->Beta[k];
    nui = lat.elems[loc]->Nu[k];
    for (j = 0; j < n_corr_[k]; j++) {
      loc = corrs_[k][j]; betaj = lat.elems[loc]->Beta[k];
      nuj = lat.elems[loc]->Nu[k];
      gsl_matrix_set
	(A_lsoc[k], i, j,
	 sqrt(betai*betaj)/(2.0*spiq)*cos(nu*M_PI-fabs(2.0*M_PI*(nui-nuj))));
    }
  }

  gsl_matrix_memcpy(U_lsoc[k], A_lsoc[k]);
  gsl_linalg_SV_decomp(U_lsoc[k], V_lsoc[k], w_lsoc[k], work);

  printf("\ngcmat singular values:\n");
  for (j = 0; j < n_corr_[k]; j++) {
    printf("%11.3e", gsl_vector_get(w_lsoc[k], j));
    if (gsl_vector_get(w_lsoc[k], j) < eps) {
      gsl_vector_set(w_lsoc[k], j, 0e0);
      printf(" (zeroed)");
    }
    if ((j+1) % 5 == 0) printf("\n");
  }
  if (n_corr_[k] % 5 != 0) printf("\n");

  if (lat.conf.trace) prt_gcmat(plane);

  gsl_vector_free(work);
}


void gcmat(LatticeType &lat, const int n_bpm, const std::vector<int> bpms,
	   const int n_corr, const std::vector<int> corrs, const int plane,
	   const bool svd)
{
  int  i;

  const bool first = (plane == 1)? first_h[0] : first_v[0];
  const int  k     = plane - 1;

  if (first) {
    if (plane == 1)
      first_h[0] = false;
    else
      first_v[0] = false;

    gsl_matrix *A = gsl_matrix_alloc(n_bpm, n_corr);
    A_lsoc[k] = A;
    gsl_matrix *U = gsl_matrix_alloc(n_bpm, n_corr);
    U_lsoc[k] = U;
    gsl_vector *w = gsl_vector_alloc(n_corr);
    w_lsoc[k] = w;
    gsl_matrix *V = gsl_matrix_alloc(n_corr, n_corr);
    V_lsoc[k] = V;
  }

  bpms_[k].clear();
  corrs_[k].clear();
  for (i = 0; i < n_bpm; i++)
    bpms_[k].push_back(bpms[i]);
  for (i = 0; i < n_corr; i++)
    corrs_[k].push_back(corrs[i]);
  n_bpm_[k] = n_bpm; n_corr_[k] = n_corr;

  if (svd) gcmat(lat, plane);
}


void gcmat(LatticeType &lat, const int bpm, const int corr,
	   const int plane)
{
  int              i;
  std::vector<int> bpms, corrs;

  const int k = plane - 1;

  n_bpm_[k]  = lat.GetnKid(bpm);
  n_corr_[k] = lat.GetnKid(corr);

  bpms_[k].clear();
  corrs_[k].clear();
  for (i = 0; i < n_bpm_[k]; i++)
    bpms.push_back(lat.Elem_GetPos(bpm, i+1));
  for (i = 0; i < n_corr_[k]; i++)
    corrs.push_back(lat.Elem_GetPos(corr, i+1));

  gcmat(lat, n_bpm_[k], bpms, n_corr_[k], corrs, plane, true);
}


void lsoc(LatticeType &lat, const int plane, const double scl)
{
  int      j;
  long int loc;

  const int k = plane - 1;

  gsl_vector *b = gsl_vector_alloc(n_bpm_[k]);
  gsl_vector *x = gsl_vector_alloc(n_corr_[k]);

  for (j = 0; j < n_bpm_[k]; j++) {
    loc = bpms_[k][j];
    gsl_vector_set(b, j, -lat.elems[loc]->BeamPos[2*k]+lat.elems[loc]->dS[k]);
  }
      
  gsl_linalg_SV_solve(U_lsoc[k], V_lsoc[k], w_lsoc[k], b, x);

  for (j = 0; j < n_corr_[k]; j++) {
    loc = corrs_[k][j];
    if (plane == 1)
      set_dbnL_design_elem(lat, lat.elems[loc]->Fnum, lat.elems[loc]->Knum, Dip,
			   -scl*gsl_vector_get(x, j), 0e0);
    else
      set_dbnL_design_elem(lat, lat.elems[loc]->Fnum, lat.elems[loc]->Knum, Dip,
			   0e0, scl*gsl_vector_get(x, j));
  }

  gsl_vector_free(b);
  gsl_vector_free(x);
}


void gtcmat(LatticeType &lat, const int plane)
{
  /* Get trajectory response matrix

                -----------
        A   = \/beta  beta  sin(2 pi(nu  - nu ))
         ij         i     j            i     j

  */

  int      i, j;
  long int loc_bpm, loc_corr;
  double   betai, betaj, nui, nuj;

  const int    k   = plane - 1;
  const double eps = 1e-4;

  gsl_vector *work = gsl_vector_alloc(n_corr_[k]);

  for (i = 0; i < n_bpm_[k]; i++) {
    loc_bpm = bpms_[k][i];
    betai = lat.elems[loc_bpm]->Beta[k]; nui = lat.elems[loc_bpm]->Nu[k];
    for (j = 0; j < n_corr_[k]; j++) {
      loc_corr = corrs_[k][j];
      betaj = lat.elems[loc_corr]->Beta[k]; nuj = lat.elems[loc_corr]->Nu[k];
      if (loc_bpm > loc_corr)
	gsl_matrix_set(A_lstc[k], i, j,
		       sqrt(betai*betaj)*sin(2.0*M_PI*(nui-nuj)));
      else
	gsl_matrix_set(A_lstc[k], i, j, 0e0);
    }
  }

  gsl_matrix_memcpy(U_lsoc[k], A_lsoc[k]);
  gsl_linalg_SV_decomp(U_lstc[k], V_lstc[k], w_lstc[k], work);

  printf("\ngtcmat singular values:\n");
  for (j = 0; j < n_corr_[k]; j++) {
    printf("%11.3e", gsl_vector_get(w_lstc[k], j));
    if (gsl_vector_get(w_lstc[k], j) < eps) {
      gsl_vector_set(w_lstc[k], j, 0e0);
      printf(" (zeroed)");
    }
    if ((j+1) % 5 == 0) printf("\n");
  }
  if (n_corr_[k] % 5 != 0) printf("\n");

  if (lat.conf.trace) prt_gcmat(plane);
}


void gtcmat(LatticeType &lat, const int n_bpm, const std::vector<int> bpms,
	    const int n_corr, const std::vector<int> corrs, const int plane,
	    const bool svd)
{
  int  i;

  const bool first = (plane == 1)? first_h[1] : first_v[1];
  const int  k     = plane - 1;

  if (first) {
    if (plane == 1)
      first_h[1] = false;
    else
      first_v[1] = false;

    gsl_matrix *A = gsl_matrix_alloc(n_bpm, n_corr);
    A_lstc[k] = A;
    gsl_matrix *U = gsl_matrix_alloc(n_bpm, n_corr);
    U_lstc[k] = U;
    gsl_vector *w = gsl_vector_alloc(n_corr);
    w_lstc[k] = w;
    gsl_matrix *V = gsl_matrix_alloc(n_corr, n_corr);
    V_lstc[k] = V;
  }

  bpms_[k].clear();
  corrs_[k].clear();
  for (i = 0; i < n_bpm; i++)
    bpms_[k].push_back(bpms[i]);
  for (i = 0; i < n_corr; i++)
    corrs_[k].push_back(corrs[i]);
  n_bpm_[k] = n_bpm; n_corr_[k] = n_corr;

  if (svd) gtcmat(lat, plane);
}


void lstc(LatticeType &lat, const int plane, const double scl)
{
  int      j;
  long int loc;

  const int k = plane - 1;

  gsl_vector *b = gsl_vector_alloc(n_bpm_[k]);
  gsl_vector *x = gsl_vector_alloc(n_corr_[k]);

  for (j = 0; j < n_bpm_[k]; j++) {
    loc = bpms_[k][j];
    gsl_vector_set(b, j, -lat.elems[loc]->BeamPos[2*k]+lat.elems[loc]->dS[k]);

    if (lat.conf.trace)
      std::cout << std::scientific << std::setprecision(5)
		<< "b[" << std::setw(3) << j << "] = "
		<< std::setw(12) << gsl_vector_get(b, j) << std::endl;
  }

  gsl_linalg_SV_solve(U_lstc[k], V_lstc[k], w_lstc[k], b, x);

  for (j = 0; j < n_corr_[k]; j++) {
    loc = corrs_[k][j];
    if (plane == 1) {
      if (lat.conf.trace)
	std::cout << std::scientific << std::setprecision(5)
		  << "(b_1L)[" << std::setw(3) << j << "] = "
		  << std::setw(12)<< -gsl_vector_get(x, j) << std::endl;

      set_dbnL_design_elem(lat, lat.elems[loc]->Fnum, lat.elems[loc]->Knum, Dip,
			   -scl*gsl_vector_get(x, j), 0e0);
    } else {
      if (lat.conf.trace)
	std::cout << std::scientific << std::setprecision(5)
		  << "(a_1L)[" << std::setw(3) << j << "] = "
		  << std::setw(12) << gsl_vector_get(x, j) << std::endl;

      set_dbnL_design_elem(lat, lat.elems[loc]->Fnum, lat.elems[loc]->Knum, Dip,
			   0e0, scl*gsl_vector_get(x, j));
    }
  }

  gsl_vector_free(b);
  gsl_vector_free(x);
}
