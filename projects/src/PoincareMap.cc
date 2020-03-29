
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

  void GetA(const int n_DOF, const double S, Matrix &M);
};


void PrtMap(const double n_DOF, ss_vect<tps> map)
{
  int    i, j;
  double val;

  const int    n_dec = 6;
  const double eps   = 1e-15;

  cout << endl;
  for (i = 1; i <= 2*n_DOF; i++) {
    for (j = 1; j <= 2*n_DOF; j++) {
      val = map[i][j];
      cout << scientific << setprecision(n_dec)
	   << setw(n_dec+8) << ((fabs(val) > eps)? val : 0e0);
    }
    cout << endl;
  }
}


void PoincareMap::GetA(const int n_DOF, const double S, Matrix &M)
{
  int    k;
  double alpha_c, nu[3];
  Matrix A, A_inv, R;

  const bool prt = !false;

  this->M = putlinmat(2*n_DOF, M);

  GDiag(2*n_DOF, S, A, A_inv, R, M, this->nu[Z_], alpha_c);
  this->A.identity();
  this->A = get_A_CS(n_DOF, putlinmat(2*n_DOF, A), nu);
  this->R = putlinmat(2*n_DOF, R);
  for (k = 0; k < 2; k++)
    this->nu[k] = globval.TotalTune[k];

  if (prt) {
    printf("\nGetA:\n  nu:");
    for (k = 0; k < n_DOF; k++)
      printf(" %7.5f", this->nu[k]);
    printf("\n");
  }
}
