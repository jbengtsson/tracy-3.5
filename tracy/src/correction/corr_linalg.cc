// Shared SVD-with-singular-value-cut substrate. See correction/corr_linalg.h.


void corr::svd_decomp_cut(double **M, const int m, const int n, double **U,
			  double *w, double **V, const double s_cut,
			  const bool prt, const int n_prt)
{
  int i, j;

  for (i = 1; i <= m; i++)
    for (j = 1; j <= n; j++)
      U[i][j] = M[i][j];

  dsvdcmp(U, m, n, w, V);

  if (prt) {
    printf("\n");
    printf("singular values: s_cut = %10.3e\n", s_cut);
    printf("\n");
  }

  for (i = 1; i <= n; i++) {
    if (prt) printf("%11.3e", w[i]);
    if (w[i] < s_cut) {
      w[i] = 0e0;
      if (prt) printf(" (zeroed)");
    }
    if (prt && (i % n_prt == 0)) printf("\n");
  }
  if (prt && (n % n_prt != 0)) printf("\n");
}


void corr::svd_backsub(double **U, double *w, double **V, const int m,
		       const int n, double b[], double x[])
{
  dsvbksb(U, w, V, m, n, b, x);
}
