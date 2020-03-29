
// Poincare Map Object.
// Johan Bengtsson 29/03/20.


struct PoincareMap {

private:
public:
  int n_DOF;      // Degrees-of-Freedom.
  double
    C,            // Circumference.
    nu[3],        // Tunes; in Floquet Space.
    ksi_1[2],     // Linear Chromaticity.
    ksi_2[2],     // Second Order Chromaticity.
    dnu_dJ[3],    // Anharmonic terms.
    D[3];         // Diffusion Coefficients.
  ss_vect<double>
  ps_cod;         // Fixed Point.
  ss_vect<tps>
    M_lat,        // Linear map.
    A,            // Transformation to Floquet Space.
    A_tp,         // Transpose ditto.
    R,            // Floquet Space Rotation.
    M_diff,       // Diffusion Matrix.
    M_Chol_tp;    // Cholesky Decomposition of Diffusion Matrix.

  void GetM(const bool cav, const bool rad);
  void GetA(void);
  void GetM_diff(void);
  void GetM_Chol(void);
};


void PrtMap(const double n_DOF, ss_vect<tps> map)
{
  int    i, j;

  const int    n_dec = 6;
  const double eps   = 1e-20;

  for (i = 0; i < 2*n_DOF; i++) {
    for (j = 0; j < 2*n_DOF; j++) {
      if (!true)
	printf("%*.*e", n_dec+8, n_dec, map[i][j]);
      else
	printf("%*.*e",
	       n_dec+8, n_dec, (fabs(map[i][j]) > eps)? map[i][j] : 0e0);
    }
    printf("\n");
  }
}


void PoincareMap::GetM(const bool cav, const bool rad)
{
  int long lastpos;
    
  const bool prt = !false;

  globval.Cavity_on = cav; globval.radiation = rad;

  C = Cell[globval.Cell_nLoc].S;

  n_DOF = (!cav)? 2 : 3;

  getcod(0e0, lastpos);
  ps_cod = globval.CODvect;

  M_lat.identity(); M_lat += globval.CODvect;
  Cell_Pass(0, globval.Cell_nLoc, M_lat, lastpos);
  M_lat -= globval.CODvect;

  if (prt) {
    printf("\nCavity %d, Radiation %d\n", cav, rad);
    cout << scientific << setprecision(6)
	 << "\nCOD:\n" << setw(14) << ps_cod << "\n";
    printf("\nM:\n");
    PrtMap(3, M_lat);
  }

  GetA();

  if (rad) GetM_diff();
}


void PoincareMap::GetA(void)
{
  int    k;
  double alpha_c, dnu[3];
  Matrix M_mat, A_mat, A_inv_mat, A_tp_mat, R_mat;

  const bool prt = !false;

  getlinmat(2*n_DOF, M_lat, M_mat);
  GDiag(2*n_DOF, C, A_mat, A_inv_mat, R_mat, M_mat, nu[Z_], alpha_c);

  A = get_A_CS(n_DOF, putlinmat(2*n_DOF, A_mat), dnu);
  if (n_DOF < 3) {
    A[ct_] = tps(0e0, ct_+1); A[delta_] = tps(0e0, delta_+1);
  }
  // Note, for radiation, A is not symplectic.
  getlinmat(2*n_DOF, A, A_tp_mat); TpMat(2*n_DOF, A_tp_mat);
  A_tp = putlinmat(2*n_DOF, A_tp_mat);
  R = putlinmat(2*n_DOF, R_mat);
  for (k = 0; k < 2; k++)
    nu[k] = atan2(R[2*k][2*k+1], R[2*k][2*k])/(2e0*M_PI);

  if (prt) {
    printf("\nnu = [");
    for (k = 0; k < n_DOF; k++) {
      printf("%7.5f", nu[k]);
      if (k < n_DOF-1) printf(", ");
    }
    printf("]\n");
    printf("\nA:\n");
    PrtMap(n_DOF, A);
    printf("\nR:\n");
    PrtMap(n_DOF, R);
  }
}


void PoincareMap::GetM_diff(void)
{
  int long     lastpos;
  int          k;
  ss_vect<tps> A1;

  const bool prt = !false;

  globval.emittance = true;

  A1 = A + globval.CODvect;
  Cell_Pass(0, globval.Cell_nLoc, A1, lastpos);
  for (k = 0; k < DOF; k++)
    D[k] = globval.D_rad[k];

  globval.emittance = false;

  M_diff.identity();
  for (k = 0; k < 2*n_DOF; k++)
    M_diff[k] *= D[k/2];
  M_diff = M_diff*A*A_tp;

  if (prt) {
    printf("\nD = [");
    for (k = 0; k < n_DOF; k++) {
      printf("%9.3e", D[k]);
      if (k < n_DOF-1) printf(", ");
    }
    printf("]\n");
    printf("\nM_diff:\n");
    PrtMap(n_DOF, M_diff);
  }

  GetM_Chol();
}


ss_vect<tps> Mat2Map(const int n, double **M)
{
  int          j, k;
  ss_vect<tps> map;

  map.zero();
  for (j = 1; j <= n; j++)
    for (k = 1; k <= n; k++)
      map[j-1] += M[j][k]*tps(0e0, k);
  return map;
}


void PoincareMap::GetM_Chol(void)
{
  int    j, k;
  double *diag, **d1, **d2, **d2_tp;
  ss_vect<tps> L, L_tp;

  const bool prt = !false;
  const int  n   = 2*n_DOF;

  diag = dvector(1, n); d1 = dmatrix(1, n, 1, n);
  d2 = dmatrix(1, n, 1, n); d2_tp = dmatrix(1, n, 1, n);

  for (j = 1; j <= n; j++)
    for (k = 1; k <= n; k++) {
      d1[j][k] = M_diff[j-1][k-1];
    }

  dcholdc(d1, n, diag);
  for (j = 1; j <= n; j++)
    for (k = 1; k <= j; k++)
      d2[j][k] = (j == k)? diag[j] : d1[j][k];
  dmtranspose(d2, n, n, d2_tp);

  M_Chol_tp = Mat2Map(n, d2_tp);

  if (prt) {
  printf("\nM_Chol_tp:\n");
  PrtMap(n_DOF, M_Chol_tp);
  }

  free_dvector(diag, 1, n); free_dmatrix(d1, 1, n, 1, n);
  free_dmatrix(d2, 1, n, 1, n); free_dmatrix(d2_tp, 1, n, 1, n);
}


void get_Poincare_Map(void)
{
  PoincareMap map;

  if (!false) no_sxt();

  if (false) {
    map.GetM(false, false);
    map.GetM(true, false);
    map.GetM(true, true);
  }

  map.GetM(true, true);
}
