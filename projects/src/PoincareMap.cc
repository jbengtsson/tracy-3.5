
// Poincare Map Object.
// Johan Bengtsson 29/03/20.


struct PoincareMap {

private:

public:
  int n_DOF;           // Degrees-of-Freedom.
  bool
    cav_on,            // RF Cavity on/off.
    rad_on;            // Classical Radiation on/off.
  double
    C,                 // Circumference.
    nu[3],             // Tunes; in Floquet Space.
    ksi_1[2],          // Linear Chromaticity.
    ksi_2[2],          // Second Order Chromaticity.
    dnu_dJ[3],         // Anharmonic terms.
    D[3];              // Diffusion Coefficients.
  ss_vect<double>
  ps_cod;              // Fixed Point.
  ss_vect<tps>
    M_lat,             // Linear map.
    A, A_tp,           // Transformation to Floquet Space and transpose.
    R,                 // Floquet Space Rotation.
    M_diff, M_Chol_tp; // Diffusion Matrix & Cholesky Decomposition.

  void GetM_lat(const bool cav, const bool rad);
  void GetA(void);
  void GetM_diff(void);
  void GetM_Chol(void);
  void print(void);
};


void PrtMap(const double n_DOF, ss_vect<tps> map)
{
  int i, j;

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


void PoincareMap::GetM_lat(const bool cav, const bool rad)
{
  int long lastpos;
  
  cav_on = cav; this->rad_on = rad;
  
  globval.Cavity_on = cav_on; globval.radiation = rad_on;

  C = Cell[globval.Cell_nLoc].S;

  n_DOF = (!cav)? 2 : 3;

  getcod(0e0, lastpos);
  ps_cod = globval.CODvect;

  M_lat.identity(); M_lat += globval.CODvect;
  Cell_Pass(0, globval.Cell_nLoc, M_lat, lastpos);
  M_lat -= globval.CODvect;

  GetA();

  if (rad) GetM_diff();
}


void PoincareMap::GetA(void)
{
  int    k;
  double alpha_c, dnu[3];
  Matrix M_mat, A_mat, A_inv_mat, A_tp_mat, R_mat;

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
}


void PoincareMap::GetM_diff(void)
{
  int long     lastpos;
  int          k;
  ss_vect<tps> A1;

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

  const int n = 2*n_DOF;

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

  free_dvector(diag, 1, n); free_dmatrix(d1, 1, n, 1, n);
  free_dmatrix(d2, 1, n, 1, n); free_dmatrix(d2_tp, 1, n, 1, n);
}


void PoincareMap::print(void)
{
  int k;

  printf("\nCavity %d, Radiation %d\n", cav_on, rad_on);
  cout << scientific << setprecision(6)
       << "\nCOD:\n" << setw(14) << ps_cod << "\n";
  printf("\nM:\n");
  PrtMap(3, M_lat);

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

  printf("\nD = [");
  for (k = 0; k < n_DOF; k++) {
    printf("%9.3e", D[k]);
    if (k < n_DOF-1) printf(", ");
  }
  printf("]\n");
  printf("\nM_diff:\n");
  PrtMap(n_DOF, M_diff);

  printf("\nM_Chol_tp:\n");
  PrtMap(n_DOF, M_Chol_tp);
}


void get_Poincare_Map(void)
{
  PoincareMap map;

  if (!false) no_sxt();

  if (false) {
    map.GetM_lat(false, false);
    map.GetM_lat(true, false);
    map.GetM_lat(true, true);
  }

  map.GetM_lat(true, true);
  map.print();
}
