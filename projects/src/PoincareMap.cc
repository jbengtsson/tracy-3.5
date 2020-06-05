
// Poincare Map Object.
//
// Computes the Poincare map for fast tracking, for a general lattice of
// arbitray complexity, e.g. with: linear coupling, imperfections, and related
// corrections. It is partitioned so that nonlinear effects can be included to
// arbitrary order by a parametric approach. In particular: momentum compaction,
// chromaticity, and anharmonic terms.
//
// Use Cases include lattice modelling for collective phenomena.
//
// Johan Bengtsson 29/03/20.

#include <random>


const int n_DOF = 3;

struct PoincareMapType;


struct BeamType {
private:
public:
  int
    n_part;            // No of particles.
  std::vector< ss_vect<double> >
    ps;                // Phase space coordinates. 
  ss_vect<double>
    mean;              // Average.
  ss_vect<tps>
    Sigma;             // Variances.

  void BeamInit_dist(const int n_part, const double eps_x, const double eps_y,
		     const double eps_z, const ss_vect<tps> &A);
  void BeamInit_Sigma(const double eps_x, const double eps_y,
		      const double eps_z, const ss_vect<tps> &A);
  void BeamStats(void);
  void print(const PoincareMapType &map);
};


struct PoincareMapType {
private:
public:
  bool
    cav_on,            // RF Cavity on/off.
    rad_on;            // Classical Radiation on/off.
  double
    C,                 // Circumference.
    E0,                // Beam Energy.
    alpha[3], beta[3], // Twiss Parameters.
    nu[3],             // Tunes; in Floquet Space.
    ksi1[2],           // Linear Chromaticity.
    ksi2[2],           // Second Order Chromaticity.
    alpha_c,           // Linear Momentum Compaction.
    dnu_dJ[3],         // Anharmonic terms.
    U0,                // Energy loss per turn.
    delta_cav,         // Momentum change due Energy Loss.
    delta_tau,         // Momentum change due Radiation Damping.
    V_RF,              // RF Voltage.
    f_RF,              // RF frequency.
    phi0,              // Synchronous phase.
    tau[3],            // Damping times.
    D[3],              // Diffusion Coefficients.
    eps[3];            // Emittance.
  ss_vect<double>
    ps_cod,            // Fixed Point.
    eta;               // Linear Dispersion.
  ss_vect<tps>
    M_num,             // numerical Poincare Map.
    M_lat,             // Linear Map for Lattice.
    M_cav,             // Linear Map for RF Cavity.
    A,                 // Transf. to Floquet Space.
    A0,                // Transf. to remove delta dependent Fixed Point.
    M_Fl,              // Map in Floquet Space.
    M_delta,           // Radiation Loss Map.
    M_tau,             // Radiation Damping Map.
    M,                 // Poincare Map.
    M_diff, M_Chol_tp;    // Diffusion Matrix & Cholesky Decomposition.

  void GetMap(void);
  void GetDisp(void);
  void GetM_cav(void);
  void GetM_delta(void);
  void GetM_tau(void);
  void GetM_diff(void);
  void GetM_Chol_tp(void);
  void GetM_rad(void);
  void GetM(const bool cav, const bool rad);
  void propagate(const int n, BeamType &beam) const;
  void propagate(const int n, ss_vect<tps> &Sigma) const;
  void print(void);
};


double DetMap(ss_vect<tps> A)
{
  Matrix A_mat;

  getlinmat(2*n_DOF, A, A_mat);
  return DetMat(2*n_DOF, A_mat);
}


ss_vect<tps> TpMap(ss_vect<tps> A)
{
  Matrix A_mat;

  getlinmat(2*n_DOF, A, A_mat);
  TpMat(2*n_DOF, A_mat);
  return putlinmat(2*n_DOF, A_mat);
}


void PrtVec(const string &str, ss_vect<double> vec)
{
  int k;

  const int n_dec = 6;

  printf("%s\n", str.c_str());
  for (k = 0; k < 2*n_DOF; k++)
    printf("%*.*e", n_dec+8, n_dec, vec[k]);
  printf("\n");
}


void PrtMap(const string &str, ss_vect<tps> map)
{
  int i, j;

  const int n_dec = 6+6;

  printf("%s\n", str.c_str());
  for (i = 0; i < 2*n_DOF; i++) {
    for (j = 0; j < 2*n_DOF; j++)
      printf("%*.*e", n_dec+8, n_dec, map[i][j]);
    printf("\n");
  }
}


void prt_vec(const string &str, const int n, double *a)
{
  int i;

  const int n_prt = 5;

  printf("%s\n", str.c_str());
  for (i = 1; i <= n; i++) {
    printf("%11.3e", a[i]);
    if (i % n_prt == 0) printf("\n");
  }
  if (n % n_prt != 0) printf("\n");
}


double get_curly_H_x(const double alpha_x, const double beta_x,
		     const ss_vect<double> &eta)
{
  double gamma_x;

  gamma_x = (1e0+sqr(alpha_x))/beta_x;
  return
    gamma_x*sqr(eta[x_])+2e0*alpha_x*eta[x_]*eta[px_]*beta_x*sqr(eta[px_]);
}


ss_vect<double> get_eps(const ss_vect<tps> Sigma, const ss_vect<tps> A)
{
  int             k;
  ss_vect<double> eps;
  ss_vect<tps>    diag;

  diag = Inv(A)*Sigma*Inv(TpMap(A));
  for (k = 0; k < 2*n_DOF; k++) {
    eps[k] = diag[k][k];
  }
  return eps;
}


void prt_Sigma(const ss_vect<tps> &Sigma, const PoincareMapType &map)
{
  int             k;
  ss_vect<double> eps;

  eps = get_eps(Sigma, map.A);

  PrtVec("\neps:", eps);
  PrtMap("\nSigma:", Sigma);
  printf("\nsigma_kk:\n ");
  for (k = 0; k < 2*n_DOF; k++)
    printf(" %11.5e", sqrt(Sigma[k][k]));
  printf("\n");
}


ss_vect<tps> get_map(double **A)
{
  int          j, k;
  ss_vect<tps> Id, B;

  Id.identity(); B.zero();
  for (j = 0; j < 2*n_DOF; j++)
    for (k = 0; k < 2*n_DOF; k++)
      B[j] += A[j+1][k+1]*Id[k];
  return B;
}


void get_mat(const ss_vect<tps> &A, double **B)
{
  int i, j;

  for (i = 0; i < 2*n_DOF; i++)
    for (j = 0; j < 2*n_DOF; j++)
      B[i+1][j+1] = A[i][j];
}


void map2vec(const ss_vect<tps> &M, double *M_vec)
{
  // Matrix vectorization.
  int i, j, k;

  k = 0;
  for (i = 0; i < 2*n_DOF; i++)
    for (j = 0; j < 2*n_DOF; j++) {
      k++;
      M_vec[k] = M[j][i];
    }
}


ss_vect<tps> vec2map(double *M_vec)
{
  // Inverse matrix vectorization.
  int          i, j, k;
  ss_vect<tps> Id, M;

  Id.identity();
  M.zero(); k = 0;
  for (i = 0; i < 2*n_DOF; i++)
    for (j = 0; j < 2*n_DOF; j++) {
      k++;
      M[j] += M_vec[k]*Id[i];
    }
  return M;
}


void Kronecker_prod(const int n, double **A, double **B, double **C)
{
  int i, j, k, l;

  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++)
      for (k = 1; k <= n; k++)
	for (l = 1; l <= n; l++)
	  C[(i-1)*n+k][(j-1)*n+l] = A[i][j]*B[k][l];
}


ss_vect<tps> get_emit(ss_vect<tps> &M, ss_vect<tps> &D)
{
  int          i, j;
  double       *D_vec, *Sigma_vec, **M_mat, **M_M_mat, **Id, **MmI, **MmI_inv;
  ss_vect<tps> Sigma;

  const int
    mat_dim     = 2*n_DOF,
    mat_vec_dim = sqr(mat_dim);

  D_vec = dvector(1, mat_vec_dim);
  Sigma_vec = dvector(1, mat_vec_dim);
  M_mat = dmatrix(1, mat_dim, 1, mat_dim);
  M_M_mat = dmatrix(1, mat_vec_dim, 1, mat_vec_dim);
  Id = dmatrix(1, mat_vec_dim, 1, mat_vec_dim);
  MmI = dmatrix(1, mat_vec_dim, 1, mat_vec_dim);
  MmI_inv = dmatrix(1, mat_vec_dim, 1, mat_vec_dim);

  map2vec(D, D_vec);
  get_mat(M, M_mat);

  Kronecker_prod(mat_dim, M_mat, M_mat, M_M_mat);

  for (i = 1; i <= mat_vec_dim; i++)
    for (j = 1; j <= mat_vec_dim; j++)
      Id[i][j] = (i == j)? 1e0 : 0e0;
 
  dmsub(Id, mat_vec_dim, mat_vec_dim, M_M_mat, MmI);
  dinverse(MmI, mat_vec_dim, MmI_inv);  
  dmvmult(MmI_inv, mat_vec_dim, mat_vec_dim, D_vec, mat_vec_dim, Sigma_vec);
  Sigma = vec2map(Sigma_vec);

  free_dvector(D_vec, 1, mat_vec_dim);
  free_dvector(Sigma_vec, 1, mat_vec_dim);
  free_dmatrix(M_mat, 1, mat_dim, 1, mat_dim);
  free_dmatrix(M_M_mat, 1, mat_vec_dim, 1, mat_vec_dim);
  free_dmatrix(Id, 1, mat_vec_dim, 1, mat_vec_dim);
  free_dmatrix(MmI, 1, mat_vec_dim, 1, mat_vec_dim);
  free_dmatrix(MmI_inv, 1, mat_vec_dim, 1, mat_vec_dim);

  return Sigma;
}


ss_vect<tps> GetSigma(const double C, const double tau[], const double D[],
	      ss_vect<tps> &A)
{
  int          k;
  double       eps[3];
  ss_vect<tps> Sigma;

  for (k = 0; k < n_DOF; k++)
    eps[k] = D[k]*tau[k]*c0/(2e0*C);
  for (k = 0; k < 2*n_DOF; k++)
    Sigma[k] = eps[k/2]*tps(0e0, k+1);

  return A*Sigma*TpMap(A);
}


void GetA(const int n_DOF, const double C, const ss_vect<tps> &M,
	  ss_vect<tps> &A, ss_vect<tps> &R, double tau[])
{
  int    k;
  double nu_z, alpha_c;
  Matrix M_mat, A_mat, A_inv_mat, R_mat;

  getlinmat(2*n_DOF, M, M_mat);
  GDiag(2*n_DOF, C, A_mat, A_inv_mat, R_mat, M_mat, nu_z, alpha_c);

  for (k = 0; k < n_DOF; k++)
    tau[k] = -C/(c0*globval.alpha_rad[k]);

  A = putlinmat(2*n_DOF, A_mat);
  if (n_DOF < 3) {
    // Coasting longitudinal plane.
    A[ct_] = tps(0e0, ct_+1);
    A[delta_] = tps(0e0, delta_+1);
  }
  R = putlinmat(2*n_DOF, R_mat);
}


ss_vect<tps> GetA0(const ss_vect<double> &eta)
{
  // Canonical Transfomation to delta dependent Fixed Point.
  int          k;
  ss_vect<tps> Id, A0;

  Id.identity(); A0.identity();
  for (k = 0; k < 2; k++)
    A0[k] += eta[k]*Id[delta_];
  // Symplectic flow.
  A0[ct_] += eta[px_]*Id[x_] - eta[x_]*Id[px_];
  return A0;
}


ss_vect<tps> GetA1(const double alpha[], const double beta[])
{
  // Canonical Transformation to Floquet Space.
  int          k;
  ss_vect<tps> Id, A;

  Id.identity(); A.zero();
  for (k = 0; k < n_DOF; k++) {
    if (k < 2) {
      A[2*k] += sqrt(beta[k])*Id[2*k];
      A[2*k+1] += -alpha[k]/sqrt(beta[k])*Id[2*k] + 1e0/sqrt(beta[k])*Id[2*k+1];
    } else {
      A[ct_] += sqrt(beta[Z_])*Id[ct_];
      A[delta_] +=
	-alpha[Z_]/sqrt(beta[Z_])*Id[ct_] + 1e0/sqrt(beta[Z_])*Id[delta_];
    }
  }
  return A;
}


ss_vect<tps> GetR(const double nu[])
{
  int          k;
  ss_vect<tps> Id, R;

  Id.identity(); R.zero();
  for (k = 0; k < n_DOF; k++) {
    R[2*k] = cos(2e0*M_PI*nu[k])*Id[2*k] + sin(2e0*M_PI*nu[k])*Id[2*k+1];
    R[2*k+1] = -sin(2e0*M_PI*nu[k])*Id[2*k] + cos(2e0*M_PI*nu[k])*Id[2*k+1];
  }
  R[ct_] = cos(2e0*M_PI*nu[Z_])*Id[ct_] + sin(2e0*M_PI*nu[Z_])*Id[delta_];
  R[delta_] = -sin(2e0*M_PI*nu[Z_])*Id[ct_] + cos(2e0*M_PI*nu[Z_])*Id[delta_];
  return R;
}


ss_vect<double> GetEta(const ss_vect<tps> &M)
{
  long int        jj[ss_dim];
  int             k;
  ss_vect<double> eta;
  ss_vect<tps>    Id, M_x;

  Id.identity();
  for (k = 0; k < ss_dim; k++)
    jj[k] = 0;
  // Only use 2x2 part of hor. matrix.
  M_x.zero();
  for (k = 0; k < 2; k++) {
    M_x[k] = M[k][x_]*Id[x_] + M[k][px_]*Id[px_];
    jj[k] = 1;
  }
  eta.zero();
  for (k = 0; k < 2; k++)
    eta[k] = M[k][delta_];
  eta = (PInv(Id-M_x, jj)*eta).cst();
  return eta;
}


void GetTwiss(const ss_vect<tps> &A, const ss_vect<tps> &R,
	      double alpha[], double beta[], ss_vect<double> &eta, double nu[])
{
  long int     jj[ss_dim];
  int          k;
  ss_vect<tps> Id, A_t, scr, A_A_tp;

  Id.identity();

  A_t.zero();
  for (k = 4; k < 6; k++)
    A_t[k] = A[k][ct_]*Id[ct_] + A[k][delta_]*Id[delta_];
  for (k = 0; k < ss_dim; k++)
    jj[k] = 0;
  jj[ct_] = 1; jj[delta_] = 1;
  scr = A*PInv(A_t, jj);
  eta.zero();
  for (k = 0; k < 2; k++)
    eta[k] = scr[k][delta_];
 
  // A_t.identity();
  // for (k = 4; k < 6; k++)
  //   A_t[k] = A[k][ct_]*Id[ct_] + A[k][delta_]*Id[delta_];
  // scr = A*Inv(A_t);
  // eta.zero();
  // for (k = 0; k < 4; k++)
  //   eta[k] = scr[k][delta_];

  for (k = 0; k < n_DOF; k++) {
    if (k < 2) {
      // Phase space rotation by betatron motion.
      alpha[k] = -(A[2*k][2*k]*A[2*k+1][2*k] + A[2*k][2*k+1]*A[2*k+1][2*k+1]);
      beta[k] = sqr(A[2*k][2*k]) + sqr(A[2*k][2*k+1]);
      nu[k] = atan2(R[2*k][2*k+1], R[2*k][2*k])/(2e0*M_PI);
    } else {
      alpha[Z_] =
	-(A[ct_][ct_]*A[delta_][ct_] + A[ct_][delta_]*A[delta_][delta_]);
      beta[Z_] = sqr(A[ct_][ct_]) + sqr(A[ct_][delta_]);
      nu[Z_] = atan2(R[ct_][delta_], R[ct_][ct_])/(2e0*M_PI);
    }
  }
}


void PoincareMapType::GetM_cav(void)
{
  const int loc = Elem_GetPos(ElemIndex("cav"), 1);

  V_RF = Cell[loc].Elem.C->Pvolt;
  f_RF = Cell[loc].Elem.C->Pfreq;
  phi0 = asin(delta_cav*E0/V_RF);

  M_cav.identity();
  M_cav[delta_] -=
    V_RF*2e0*M_PI*f_RF*cos(phi0)/(1e9*globval.Energy*c0)*tps(0e0, ct_+1);
}


void PoincareMapType::GetM_delta(void)
{
  int k;

  M_delta.identity();
  for (k = 0; k < n_DOF; k++) {
    M_delta[2*k] *= 1e0 + delta_tau;
    M_delta[2*k+1] /= 1e0 + delta_tau;
  }
}


void PoincareMapType::GetM_tau(void)
{
  int k;
 
  M_tau.zero();
  for (k = 0; k < n_DOF; k++) {
    M_tau[2*k] = exp(-C/(c0*tau[k]))*tps(0e0, 2*k+1);
    M_tau[2*k+1] = exp(-C/(c0*tau[k]))*tps(0e0, 2*k+2);
  }
}


void PoincareMapType::GetM_diff(void)
{
  int long     lastpos;
  int          k;
  ss_vect<tps> As, D_diag;

  globval.Cavity_on = true; globval.radiation = true;
  globval.emittance = true;

  As = A + ps_cod;
  Cell_Pass(0, globval.Cell_nLoc, As, lastpos);

  globval.emittance = false;

  for (k = 0; k < n_DOF; k++) {
    D[k] = globval.D_rad[k];
    eps[k] = D[k]*tau[k]*c0/(2e0*C);
  }

  D_diag.zero();
  for (k = 0; k < 2*n_DOF; k++)
    D_diag[k] = D[k/2]*tps(0e0, k+1);
  M_diff = A*D_diag*TpMap(A);

  GetM_Chol_tp();
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


void PoincareMapType::GetM_Chol_tp(void)
{
  int          j, k, j1, k1;
  double       *diag, **d1, **d2;

  const int n = 2*n_DOF;

  diag = dvector(1, n); d1 = dmatrix(1, n, 1, n);
  d2 = dmatrix(1, n, 1, n);

  // Remove vertical plane.
  for (j = 1; j <= 4; j++) {
    j1 = (j < 3)? j : j+2;
    for (k = 1; k <= 4; k++) {
      k1 = (k < 3)? k : k+2;
      d1[j][k] = M_diff[j1-1][k1-1];
    }
  }

  dcholdc(d1, 4, diag);

  for (j = 1; j <= 4; j++)
    for (k = 1; k <= j; k++)
      d2[j][k] = (j == k)? diag[j] : d1[j][k];

  M_Chol_tp.zero();
  for (j = 1; j <= 4; j++) {
     j1 = (j < 3)? j : j+2;
     for (k = 1; k <= 4; k++) {
       k1 = (k < 3)? k : k+2;
       M_Chol_tp[j1-1] += d2[j][k]*tps(0e0, k1);
     }
  }

  free_dvector(diag, 1, n); free_dmatrix(d1, 1, n, 1, n);
  free_dmatrix(d2, 1, n, 1, n);
}


void PoincareMapType::GetMap(void)
{
  int long lastpos;

  globval.Cavity_on = cav_on; globval.radiation = rad_on;

  getcod(0e0, lastpos);

  M_num = putlinmat(2*n_DOF, globval.OneTurnMat);
  ps_cod = globval.CODvect;

  U0 = globval.dE*E0;
  delta_cav = U0/E0;
  alpha_c = M_num[ct_][delta_]/C;
}


void PoincareMapType::GetM(const bool cav, const bool rad)
{
  cav_on = cav; rad_on = rad;

  C = Cell[globval.Cell_nLoc].S; E0 = 1e9*globval.Energy;

  GetMap();
  GetA(n_DOF, C, M_num, A, M_Fl, tau);
  delta_tau = exp(-C/(c0*tau[Z_])) - 1e0;
  GetM_tau();
  M_tau = A*M_tau*Inv(A);
  GetM_delta();
  GetM_cav();

  M = M_num;

  // Nota Bene: the Twiss parameters are not used; i.e., included for
  // convenience, e.g. cross checks, etc.
  GetA(n_DOF, C, M_num, A, M_Fl, tau);
  GetTwiss(A, M_Fl, alpha, beta, eta, nu);
  eta = GetEta(M);

  if (rad) GetM_diff();
}


void get_stats(const int n, const ss_vect<double> &sum,
	       const ss_vect<double> &sum2, ss_vect<double> &m,
	       ss_vect<double> &s)
{
  int k;

  for (k = 0; k < 2*n_DOF; k++) {
    m[k] = sum[k]/n;
    s[k] = (n*sum2[k]-sqr(sum[k]))/(n*(n-1e0));
  }
}


void prt_eps(FILE *outf, const int n, const ss_vect<tps> Sigma,
	     const PoincareMapType &map)
{
  int             k;
  ss_vect<double> eps;

  eps = get_eps(Sigma, map.A);
  fprintf(outf, "%7d", n);
  for (k = 0; k < 2*n_DOF; k++)
    fprintf(outf, "%13.5e", eps[k]);
  fprintf(outf, "%13.5e %13.5e\n",
	  sqrt(Sigma[delta_][delta_]), sqrt(Sigma[ct_][ct_]));
}


void PoincareMapType::propagate(const int n, BeamType &beam) const
{
  int             i, j = 0, k, n_rnd;
  double          rnd;
  ss_vect<double> X, sum, sum2, m, s;
  FILE            *outf;

  typedef minstd_rand0 default_random_engine;

  std::default_random_engine       rand;
  std::normal_distribution<double> norm_ranf(0e0, 1e0);

  const int    n_prt     = 10;
  const string file_name = "propagate.out";

  outf = file_write(file_name.c_str());

  sum.zero(); sum2.zero(); n_rnd = 0;
  for (i = 1; i <= n; i++) {
    for (j = 0; j < (int)beam.ps.size(); j++) {
      n_rnd++;
      for (k = 0; k < 2*n_DOF; k++) {
	rnd = norm_ranf(rand); sum[k] += rnd; sum2[k] += sqr(rnd);
	X[k] = rnd;
      }
      beam.ps[j] = (M*beam.ps[j]+M_Chol_tp*X).cst();
    }

    beam.BeamStats();

    if (i % n_prt == 0) prt_eps(outf, i, beam.Sigma, *this);
  }
  if (n % n_prt != 0) prt_eps(outf, n, beam.Sigma, *this);

  fclose(outf);

  get_stats(n_rnd, sum, sum2, m, s);
  PrtVec("\nmean_X:", m);
  PrtVec("sigma_X:", s);
}


void PoincareMapType::propagate(const int n, ss_vect<tps> &Sigma) const
{
  int             j, k, n_rnd;
  double          rnd;
  ss_vect<double> sum, sum2, m, s;
  ss_vect<tps>    X;
  FILE            *outf;

  const int          n_prt     = 1;
  const string       file_name = "propagate.out";
  const ss_vect<tps> M_Chol    = TpMap(M_Chol_tp);

  std::default_random_engine       rand;
  std::normal_distribution<double> norm_ranf(0e0, 1e0);

  outf = file_write(file_name.c_str());

  n_rnd = 0;
  for (j = 1; j <= n; j++) {
    n_rnd++;
    for (k = 0; k < 2*n_DOF; k++) {
      rnd = norm_ranf(rand); sum[k] += rnd; sum2[k] += sqr(rnd);
      X[k] = rnd*tps(0e0, k+1);
    }
    Sigma = M*TpMap(M*Sigma) + M_Chol_tp*sqr(X)*M_Chol;

    if (j % n_prt == 0) prt_eps(outf, j, Sigma, *this);
  }
  if (n % n_prt != 0) prt_eps(outf, n, Sigma, *this);

  fclose(outf);

  get_stats(n_rnd, sum, sum2, m, s);
  PrtVec("\nmean_X:", m);
  PrtVec("sigma_X:", s);
}


ss_vect<tps> GetM_track(const bool cav_on, const bool rad_on,
			const ss_vect<tps> &ps_cod)
{
  int long     lastpos;
  ss_vect<tps> M;

  globval.Cavity_on = cav_on; globval.radiation = rad_on;
  M.identity(); M += ps_cod;
  Cell_Pass(0, globval.Cell_nLoc, M, lastpos);
  return M;
}


ss_vect<tps> GetMTwiss(const double alpha[], const double beta[],
		       const ss_vect<double> eta, const double nu[])
{
  ss_vect<tps> Id, A, M;

  Id.identity();
  A = GetA0(eta)*GetA1(alpha, beta);
  M = A*GetR(nu)*Inv(A);
  return M;
}


void PoincareMapType::print(void)
{
  printf("\nCavity %d, Radiation %d\n", cav_on, rad_on);
  printf("\nC [m]      = %7.5f\n", C);

  printf("alpha      = [%9.5f, %9.5f, %9.5f]\n",
	 alpha[X_], alpha[Y_], alpha[Z_]);
  printf("beta       = [%9.5f, %9.5f, %9.5f]\n",
	 beta[X_], beta[Y_], beta[Z_]);
  printf("eta        = [%9.5f, %9.5f]\n", eta[x_], eta[px_]);
  printf("nu         = [%9.5f, %9.5f, %9.5f]\n", nu[X_], nu[Y_], nu[Z_]);

  printf("V_RF       = %7.5e\n", V_RF);
  printf("f_RF [MhZ] = %7.5f\n", 1e-6*f_RF);
  printf("phi0       = %7.5f\n", phi0*180.0/M_PI);

  printf("E0 [GeV]   = %7.5f\n", 1e-9*E0);
  printf("U0 [keV]   = %7.5f\n", 1e-3*U0);
  printf("delta_tau  = %11.5e\n", delta_tau);
  printf("delta_cav  = %11.5e\n", delta_cav);
  printf("alpha_c    = %11.5e\n", alpha_c);
  printf("tau        = [%11.5e, %11.5e, %11.5e]\n", tau[X_], tau[Y_], tau[Z_]);
  printf("D          = [%11.5e, %11.5e, %11.5e]\n", D[X_], D[Y_], D[Z_]);
  printf("eps        = [%11.5e, %11.5e, %11.5e]\n", eps[X_], eps[Y_], eps[Z_]);

  cout << scientific << setprecision(6) << "\nCOD:\n" << setw(14) << ps_cod
       << "\n";

  PrtMap("\nM_num (input map):", M_num);
  printf("\nDet{M_num}-1 = %10.3e\n", DetMap(M_num)-1e0);

  PrtMap("\nM = A*M_Fl*A^-1:", A*M_Fl*Inv(A));
  printf("\nDet{A*M_Fl*A^-1}-1 = %10.3e\n",
	 DetMap(A*M_Fl*Inv(A))-1e0);

  PrtMap("\nA:", A);
  printf("\nDet{A}-1 = %10.3e\n", DetMap(A)-1e0);

  PrtMap("\nM_lat (w/o Cavity):", M_lat);
  printf("\nDet{M_lat}-1 = %10.3e\n", DetMap(M_lat)-1e0);

  if (!false) {
    ss_vect<tps> M1;
    M1 = GetM_track(false, true, ps_cod);
    PrtMap("\nM_lat; cross check from tracking:", M1);
    printf("\nDet{M_lat}-1 = %10.3e\n", DetMap(M1)-1e0);
  }

  PrtMap("\nM_cav:", M_cav);

  if (rad_on) {
    PrtMap("\nM_delta (radiation loss):", M_delta);
    printf("\nDet{M_delta}-1 = %10.3e\n", DetMap(M_delta)-1e0);
    PrtMap("\nM_tau (radiation damping):", M_tau);
    printf("\nDet{M_tau}-1 = %10.3e\n", DetMap(M_tau)-1e0);

    PrtMap("\nM_diff:", M_diff);
    PrtMap("\nM_Chol_tp:", M_Chol_tp);
  }
}


void BeamType::BeamInit_dist(const int n_part, const double eps_x,
			     const double eps_y, const double eps_z,
			     const ss_vect<tps> &A)
{
  int             j, k;
  double          phi;
  ss_vect<double> ps_Fl;

  this->n_part = n_part;

  const double eps[] = {eps_x, eps_y, eps_z};

  for (j = 0; j < n_part; j++) {
    ps_Fl.zero();
    for (k = 0; k < n_DOF; k++) {
      phi = 2e0*M_PI*ranf();
      // Single particle amplitude.
      ps_Fl[2*k] = sqrt(2e0*eps[k])*cos(phi);
      ps_Fl[2*k+1] = -sqrt(2e0*eps[k])*sin(phi);
    }
    this->ps.push_back((A*ps_Fl).cst());
  }
}


void BeamType::BeamInit_Sigma(const double eps_x, const double eps_y,
			      const double eps_z, const ss_vect<tps> &A)
{
  int          k;

  const double eps[] = {eps_x, eps_y, eps_z};

  Sigma.zero();
  for (k = 0; k < 2*n_DOF; k++)
    Sigma[k] += eps[k/2]*tps(1e0, k+1);
  Sigma = A*Sigma*TpMap(A);
}


void get_stats(const int n, const ss_vect<double> &sum,
	       const ss_vect<tps> &sum2, ss_vect<double> &mean,
	       ss_vect<tps> &Sigma)
{
  int j, k;

  Sigma.zero();
  for (j = 0; j < 2*n_DOF; j++) {
    mean[j] = sum[j]/n;
    for (k = 0; k < 2*n_DOF; k++) {
      Sigma[j] += (n*sum2[j][k]-sum[j]*sum[k])/(n*(n-1e0))*tps(0e0, k+1);
    }
  }
}


void BeamType::BeamStats(void)
{
  int             i, j, k;
  ss_vect<double> sum;
  ss_vect<tps>    sum2;

  sum.zero(); sum2.zero();
  for (i = 0; i < n_part; i++) {
    sum += ps[i];
    for (j = 0; j < 2*n_DOF; j++)
      for (k = 0; k < 2*n_DOF; k++)
	sum2[j] += ps[i][j]*ps[i][k]*tps(0e0, k+1);
  }
  get_stats(n_part, sum, sum2, mean, Sigma);
}


void BeamType::print(const PoincareMapType &map) { prt_Sigma(Sigma, map); }


void BenchMark(const int n_part, const int n_turn, const PoincareMapType &map)
{
  BeamType beam;

  printf("\nBenchMark:\n");
  beam.BeamInit_dist(n_part, 0e-9, 1e-9, 0e-3, map.A);
  beam.BeamStats(); beam.print(map);
  map.propagate(n_turn, beam);
  beam.BeamStats(); beam.print(map);
}


void BenchMark(const int n_turn, const PoincareMapType &map)
{
  BeamType beam;
 
  printf("\nBenchMark:\n");
  beam.BeamInit_Sigma(0e-9, 1e-9, 0e-9, map.A);
  beam.print(map);
  map.propagate(n_turn, beam.Sigma);
  beam.print(map);
}


void tst_Chol_decomp(const PoincareMapType map)
{
  int             i, j, k;
  double          rnd;
  ss_vect<double> ps, X, sum, sum2_vec, m, s, mean;
  ss_vect<tps>    sum2, Sigma;

  std::default_random_engine       rand;
  std::normal_distribution<double> norm_ranf(0e0, 1e0);

  const int n_rnd = 10000;

  if (false) {
    PrtMap("\nL^T*L:", map.M_Chol_tp*TpMap(map.M_Chol_tp));
    PrtMap("D:", map.M_diff);
  }

  X.zero(); sum.zero(); sum2_vec.zero();
  for (j = 0; j < n_rnd; j++)
    for (k = 0; k < 2*n_DOF; k++) {
      rnd = norm_ranf(rand);
      sum[k] += rnd; sum2_vec[k] += sqr(rnd);
      X[k] += rnd;
    }
  for (k = 0; k < 2*n_DOF; k++)
    X[k] /= n_rnd;

  get_stats(n_rnd, sum, sum2_vec, m, s);
  PrtVec("\nmean:", m);
  PrtVec("sigma:", s);

  sum.zero(); sum2.zero();
  for (i = 0; i < n_rnd; i++) {
    for (k = 0; k < 2*n_DOF; k++) {
      rnd = norm_ranf(rand);
      X[k] = rnd;
    }
    ps = (map.M_Chol_tp*X).cst();
    sum += ps;
    for (j = 0; j < 2*n_DOF; j++)
      for (k = 0; k < 2*n_DOF; k++)
	sum2[j] += ps[j]*ps[k]*tps(0e0, k+1);
  }

  get_stats(n_rnd, sum, sum2, mean, Sigma);
  PrtVec("\nmean:", mean);
  PrtMap("Sigma:", Sigma);
  PrtMap("D:", map.M_diff);
}


void chk_map(const PoincareMapType &map)
{
  double          tau[3], alpha[3], beta[3], nu[3], alpha_c;
  ss_vect<double> eta, eta1;
  ss_vect<tps>    M, M_Fl, A, A0, A1;

  M = map.M_num;
  PrtMap("\nM:", M);
  M = Inv(map.M_tau)*M;
  PrtMap("\nInv(M_tau)*M:", M);
  M = Inv(map.M_cav)*M;
  PrtMap("\nInv(M_cav)*Inv(M_tau)*M:", M);
  M = Inv(map.M_delta)*M;
  PrtMap("\nInv(M_delta)*Inv(M_cav)*Inv(M_tau)*M:", M);
  eta = GetEta(M);
  eta1 = eta; eta1[delta_] = 1e0;
  cout << scientific << setprecision(12)
       << "\neta:  \n" << setw(20) << eta1 << "\n";
  cout << scientific << setprecision(12)
       << "M*eta:\n" << setw(20) << (M*eta1).cst() << "\n";
  A0 = GetA0(eta);
  M = Inv(A0)*M*A0;
  PrtMap("\nInv(A0)*M*A0:", M);

  GetA(2, map.C, M, A, M_Fl, tau);
  GetTwiss(A, M_Fl, alpha, beta, eta1, nu);
  alpha_c = M[ct_][delta_]/map.C;

  printf("\nReconstruct Lattice:\n");
  A0 = GetA0(eta);
  A1 = GetA1(alpha, beta);
  M = GetR(nu);
  M[ct_] += alpha_c*map.C*tps(0e0, delta_+1);
  PrtMap("\nR:", M);
  M = A1*M*Inv(A1);
  PrtMap("\nInv(A1)*M*A1:", M);
  M = A0*M*Inv(A0);
  PrtMap("\nA0*A1*M*Inv(A1)*Inv(A0):", M);
  M = map.M_tau*map.M_cav*map.M_delta*M;
  PrtMap("\nM_rec:", M);
  PrtMap("\nM_num:", map.M_num);
}


void get_Poincare_Map(void)
{
  ss_vect<tps>    Sigma;
  PoincareMapType map;

  if (!false) no_sxt();

  map.GetM(true, true); map.print();

  if (false) tst_Chol_decomp(map);

  if (false) {
    Sigma = get_emit(map.M, map.M_diff);
    prt_Sigma(Sigma, map);
    Sigma = GetSigma(map.C, map.tau, map.D, map.A);
    prt_Sigma(Sigma, map);
  }

  if (!false) chk_map(map);

  if (false) {
    if (true)
      BenchMark(30000, map);
    else
      BenchMark(10000, 30000, map);
    Sigma = GetSigma(map.C, map.tau, map.D, map.A);
    prt_Sigma(Sigma, map);
  }
}
