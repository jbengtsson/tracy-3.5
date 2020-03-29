
// Class for Poincare map.
// Johan Bengtsson 29/03/20.


struct PoincareMap {

private:
public:
  int n_DOF;             // Degrees-of-Freedom.
  double
    nu[3],               // Tunes; in Floquet Space.
    ksi_1[2],            // Linear Chromaticity.
    ksi_2[2],            // Second Order Chromaticity.
    dnu_dJ[3];           // Anharmonic terms.
  ss_vect<double> x_cod; // Fixed Point.
  ss_vect<tps>
    M,                   // Linear map.
    A,                   // Transformation to Floquet Space.
    R;                   // Floquet Space Rotation.

  void GetA(const int n_DOF, const double S, ss_vect<tps> &M);
};


void PrtMap(const double n_DOF, ss_vect<tps> map)
{
  int    i, j;
  double val;

  const int    n_dec = 6;
  const double eps   = 1e-15;

  for (i = 0; i < 2*n_DOF; i++) {
    for (j = 0; j < 2*n_DOF; j++) {
      val = map[i][j];
      printf("%*.*e", n_dec+8, n_dec, (fabs(val) > eps)? val : 0e0);
    }
    printf("\n");
  }
}


void PoincareMap::GetA(const int n_DOF, const double S, ss_vect<tps> &map)
{
  int    k;
  double alpha_c, dnu[3];
  Matrix M_mat, A_mat, A_inv_mat, R_mat;

  const bool prt = !false;

  M = map;

  getlinmat(2*n_DOF, M, M_mat);
  GDiag(2*n_DOF, S, A_mat, A_inv_mat, R_mat, M_mat, nu[Z_], alpha_c);

  A = get_A_CS(n_DOF, putlinmat(2*n_DOF, A_mat), dnu);
  if (n_DOF < 3) {
    A[ct_] = tps(0e0, ct_+1);
    A[delta_] = tps(0e0, delta_+1);
  }
  R = putlinmat(2*n_DOF, R_mat);
  for (k = 0; k < 2; k++)
    nu[k] = globval.TotalTune[k];

  if (prt) {
    printf("\nnu = [");
    for (k = 0; k < n_DOF; k++) {
      printf("%7.5f", nu[k]);
      if (k < n_DOF-1) printf(", ");
    }
    printf("]\n");
    printf("M:\n");
    PrtMap(n_DOF, M);
    printf("A:\n");
    PrtMap(n_DOF, A);
    printf("R:\n");
    PrtMap(n_DOF, R);
  }
}
