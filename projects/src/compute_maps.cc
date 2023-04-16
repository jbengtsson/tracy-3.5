#define NO 1

#include "tracy_lib.h"

int no_tps = NO;

void prt_vec(const int n_dof, const string &str, const ss_vect<double> vec)
{
  const int n_dec = 6;

  printf("%s\n", str.c_str());
  for (int k = 0; k < 2*n_dof; k++)
    printf("%*.*e", n_dec+8, n_dec, vec[k]);
  printf("\n");
}


void prt_map(const int n_dof, const string &str, const ss_vect<tps> map)
{
  const int n_dec = 6;

  printf("%s\n", str.c_str());
  for (int i = 0; i < 2*n_dof; i++) {
    for (int j = 0; j < 2*n_dof; j++)
      printf("%*.*e", n_dec+8, n_dec, map[i][j]);
    printf("\n");
  }
}


ss_vect<tps> tp_map(const int n_dof, const ss_vect<tps> &A)
{
  // Matrix transpose.
  Matrix A_mat;

  getlinmat(2*n_dof, A, A_mat);
  TpMat(2*n_dof, A_mat);
  return putlinmat(2*n_dof, A_mat);
}


ss_vect<tps> compute_M(void)
{
  // Compute the Poincaré map.
  long int     lastpos;
  ss_vect<tps> M;

  M.identity();
  Cell_Pass(0, globval.Cell_nLoc, M, lastpos);
  return M;
}


ss_vect<tps> compute_M(const ss_vect<double> fixed_point)
{
  // Compute the Poincaré map for the fix point.
  long int     lastpos;
  ss_vect<tps> M;

  M.identity();
  M += fixed_point;
  Cell_Pass(0, globval.Cell_nLoc, M, lastpos);
  M -= fixed_point;
  return M;
}


ss_vect<double> compute_eta(const ss_vect<tps> &M)
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


ss_vect<tps> compute_A0(const ss_vect<double> &eta)
{
  // Canonical Transfomation to delta dependent Fixed Point.
  int          k;
  ss_vect<tps> Id, A0;

  Id = A0.identity();
  for (k = 0; k < 2; k++)
    A0[k] += eta[k]*Id[delta_];
  // Symplectic flow.
  A0[ct_] += eta[px_]*Id[x_] - eta[x_]*Id[px_];
  return A0;
}


ss_vect<tps> compute_A
(const double C, const ss_vect<tps> &M, double nu[], double alpha_rad[],
 const bool cav_on)
{
  // Poincaré map Diagonalization.
  // The damping coeffients alpha[] are obtained from the eigen values.
  int          k;
  double       nu_s, alpha_c;
  Matrix       M_mat, A_mat, A_inv_mat, R_mat;
  ss_vect<tps> A;

  const bool prt   = !false;
  const int  n_dof = (cav_on)? 3 : 2;

  getlinmat(2*n_dof, M, M_mat);
  GDiag(2*n_dof, C, A_mat, A_inv_mat, R_mat, M_mat, nu_s, alpha_c);
  for (int k = 0; k < n_dof; k++)
    alpha_rad[k] = globval.alpha_rad[k];

  if (prt) {
    printf("\ncompute_A:\n  nu = [");
    for (k = 0; k < nd_tps; k++) {
      nu[k] = GetAngle(R_mat[2*k][2*k], R_mat[2*k][2*k+1])/(2e0*M_PI);
      printf("%21.15e%s", nu[k], (k < 2)? ", " : "]\n");
    }
  }

  A = putlinmat(2*n_dof, A_mat);
  if (n_dof != 3) {
    // Coasting longitudinal plane.
    A[ct_] = tps(0e0, ct_+1);
    A[delta_] = tps(0e0, delta_+1);
  }
  return A;
}


void compute_tau
(const double C, const double alpha_rad[], double tau[])
{
  // Compute the damping times.
  
  for (int k = 0; k < nd_tps; k++)
    tau[k] = -C/(c0*alpha_rad[k]);
}


ss_vect<tps> compute_M_tau(const double alpha[])
{
  int          k;
  ss_vect<tps> M_tau;
 
  M_tau.zero();
  for (k = 0; k < nd_tps; k++) {
    M_tau[2*k] = exp(alpha[k])*tps(0e0, 2*k+1);
    M_tau[2*k+1] = exp(alpha[k])*tps(0e0, 2*k+2);
  }
  return M_tau;
}


void compute_D
(const ss_vect<double> &fixed_point, const ss_vect<tps> &A, double D[])
{
  // Compute the diffusion coefficients.

  int long     lastpos;
  ss_vect<tps> As, D_diag;

  globval.Cavity_on = globval.radiation = globval.emittance = true;

  As = A + fixed_point;
  Cell_Pass(0, globval.Cell_nLoc, As, lastpos);

  globval.emittance = false;

  for (int k = 0; k < nd_tps; k++)
    D[k] = globval.D_rad[k];
}

void compute_eps
(const double C, const double tau[], const double D[], double eps[])
{
  // Compute the eigen emittances.

  for (int k = 0; k < nd_tps; k++)
    eps[k] = D[k]*tau[k]*c0/(2e0*C);
}


void compute_bunch_size
(const ss_vect<tps> &A, const double eps[], double &sigma_s,
 double &sigma_delta)
{
  // Compute the bunch size: sigma_s & sigma_delta.
  double alpha_s, beta_s, gamma_s;

  alpha_s = -A[ct_][ct_]*A[delta_][ct_] - A[ct_][delta_]*A[delta_][delta_];
  beta_s = sqr(A[ct_][ct_]) + sqr(A[ct_][delta_]);
  gamma_s = (1e0+sqr(alpha_s))/beta_s;

  sigma_s = sqrt(beta_s*eps[Z_]);
  sigma_delta = sqrt(gamma_s*eps[Z_]);
}


ss_vect<tps> compute_D_mat(const ss_vect<tps> &A, const double D[])
{
  // Compute the diffusion matrix.
  ss_vect<tps> D_diag;

  D_diag.zero();
  for (int k = 0; k < 2*nd_tps; k++)
    D_diag[k] = D[k/2]*tps(0e0, k+1);
  return A*D_diag*tp_map(nd_tps, A);
}


void prt_map(const string file_name, const ss_vect<tps> map)
{
  ofstream outf;

  file_wr(outf, file_name.c_str());
  outf << scientific << setprecision(3) << setw(11) << map;
  outf.close();
}


void prt_mat(const int n, const string &str, double **M)
{
  int i, j;

  const int n_dec = 8;

  printf("%s\n", str.c_str());
  for (i = 1; i <= n; i++) {
    for (j = 1; j <= n; j++)
      printf("%*.*e", n_dec+8, n_dec, M[i][j]);
    printf("\n");
  }
}


ss_vect<tps> compute_M_Chol(const ss_vect<tps> &M_diff)
{
  // Compute the Cholesky decomposition: D = L^T L.
  int          j, k, j1, k1;
  double       *diag, **d1, **d2;
  ss_vect<tps> M_Chol_t;

  const int n = 2*nd_tps;

  diag = dvector(1, n);
  d1 = dmatrix(1, n, 1, n);
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
    for (k = 1; k <= 4; k++)
      if (k <= j)
	d2[j][k] = (j == k)? diag[j] : d1[j][k];
      else
	d2[j][k] = 0e0;

  M_Chol_t.zero();
  for (j = 1; j <= 4; j++) {
    j1 = (j < 3)? j : j+2;
    for (k = 1; k <= 4; k++) {
      k1 = (k < 3)? k : k+2;
      M_Chol_t[j1-1] += d2[j][k]*tps(0e0, k1);
    }
  }

  free_dvector(diag, 1, n);
  free_dmatrix(d1, 1, n, 1, n);
  free_dmatrix(d2, 1, n, 1, n);

  return tp_map(3, M_Chol_t);
}


ss_vect<tps> compute_M_cav
(const string &cav_name, double phi0)
{
  tps          ct0;
  ss_vect<tps> M;
  CavityType   *C;

  const int loc = Elem_GetPos(ElemIndex(cav_name.c_str()), 1);

  C = Cell[loc].Elem.C;

  ct0 = tps(0e0, ct_+1);
  M.identity();
#if 1
  C->phi_RF = phi0;
  M[ct_] += ct0.cst();
  Cav_Pass(Cell[loc], M);
#else
  M[delta_] += -C->V_RF/(1e9*globval.Energy)*sin(2e0*M_PI*C->f_RF*ct0/c0+phi0);
#endif
  M[delta_] -= M[delta_].cst();
  return M;
}


void compute_M_delta(const double alpha[])
{
  int          k;
  ss_vect<tps> M_delta;

  M_delta.identity();
  for (k = 0; k < nd_tps; k++) {
    M_delta[2*k] *= 1e0 + alpha[k];
    M_delta[2*k+1] /= 1e0 + alpha[k];
  }
}


double DetMap(ss_vect<tps> A)
{
  Matrix A_mat;

  getlinmat(2*nd_tps, A, A_mat);
  return DetMat(2*nd_tps, A_mat);
}


void compute_map
(const ss_vect<tps> M_2D, const ss_vect<tps> M_cav,
 const ss_vect<tps> M_3D_no_rad, const ss_vect<tps> A_3D_no_rad,
 const ss_vect<tps> M_tau, const ss_vect<tps> M_3D, const ss_vect<tps> A_3D,
 const double alpha[], const double U0)
{
  ss_vect<tps> M, M1;

  printf("\nCompute map:\n");
  M = A_3D*M_tau*Inv(A_3D_no_rad)*M_2D*M_cav*A_3D_no_rad*Inv(A_3D);
  prt_map(nd_tps,
	  "\nM_3D = A_3D*M_tau*A_3D_no_rad^-1*M_2D*M_cav*A_3D_no_rad*A_3D^1:",
	  M);
  prt_map(nd_tps, "\nM_3D:", M_3D);
  prt_map(nd_tps, "\nM_tau:", M_tau);
  M = A_3D*Inv(A_3D_no_rad);
  prt_map(nd_tps, "\nA_3D.A_3D_no_rad^-1:", M);
  M1.identity();
  for (int k = 0; k < 2*nd_tps; k++)
    M1[k] *= (k % 2 == 0)? pow(M_tau[k][k], 2) : 1e0/pow(M_tau[k][k], 2);
  prt_map(nd_tps, "\nM1:", M1);
  prt_map(nd_tps, "\nM1.A_3D_no_rad:", M1*A_3D_no_rad);
  prt_map(nd_tps, "\nA_3D:", A_3D);
}


void compute_A_2D(const ss_vect<double> &eta, const double alpha[],
		  const double beta[], ss_vect<tps> &A0, ss_vect<tps> &A1)
{
  int          k;
  ss_vect<tps> Id;

  const double n_dof = 2;

  Id.identity();
  A0 = compute_A0(eta);

  A1.identity();
  for (k = 0; k < n_dof; k++) {
    A1[2*k]   = sqrt(beta[k])*Id[2*k]; 
    A1[2*k+1] = (-alpha[k]*Id[2*k] + Id[2*k+1])/sqrt(beta[k]);
  }
}


ss_vect<tps> compute_R_2D(const double Circ, const double alpha_c,
			  const double nu[], ss_vect<tps> &A)
{
  int          k;
  double       dct;
  ss_vect<tps> Id, R;

  const double n_dof = 2;

  Id = R.identity();
  for (k = 0; k < n_dof; k++) {
    R[2*k] = cos(2e0*M_PI*nu[k])*Id[2*k] + sin(2e0*M_PI*nu[k])*Id[2*k+1];
    R[2*k+1] = -sin(2e0*M_PI*nu[k])*Id[2*k] + cos(2e0*M_PI*nu[k])*Id[2*k+1];
  }
  // Compute time-of-flight contribution from dispersion.
  dct = (A*R*Inv(A))[ct_][delta_];
  R[ct_] += (-dct + Circ*alpha_c)*Id[delta_];
  return R;
}


void compute_M_2D(void)
{
  ss_vect<double> eta;
  ss_vect<tps>    M, A, A0, A1, R;

  globval.radiation = globval.Cavity_on = false;
  Ring_GetTwiss(true, 0e0);
  M = putlinmat(2*nd_tps, globval.OneTurnMat);
  // prt_map(nd_tps, "\nM:", M);

  const double
    Circ    = Cell[globval.Cell_nLoc].S,
    alpha_c = globval.Alphac,
    alpha[] = {Cell[0].Alpha[X_], Cell[0].Alpha[Y_]},
    beta[]  = {Cell[0].Beta[X_], Cell[0].Beta[Y_]},
    nu[]    = {globval.TotalTune[X_], globval.TotalTune[Y_]};

  eta.zero();
  eta[x_] = Cell[0].Eta[X_];
  eta[px_] = Cell[0].Etap[X_];

  compute_A_2D(eta, alpha, beta, A0, A1);
  A = A0*A1;
  R = compute_R_2D(Circ, alpha_c, nu, A);
  M = A*R*Inv(A);

  prt_map(nd_tps, "\nM:", M);
}


ss_vect<tps> compute_A_3D(const double alpha[], const double beta[])
{
  int          k;
  ss_vect<tps> Id, A;

  Id = A.identity();
  for (k = 0; k < nd_tps; k++) {
    if (k < 2) {
      A[2*k] = sqrt(beta[k])*Id[2*k]; 
      A[2*k+1] = (-alpha[k]*Id[2*k] + Id[2*k+1])/sqrt(beta[k]);
    } else {
      A[ct_] = sqrt(beta[k])*Id[ct_]; 
      A[delta_] = (-alpha[k]*Id[ct_] + Id[delta_])/sqrt(beta[k]);
    }
  }
  return A;
}


ss_vect<tps> compute_R_3D(const double nu[])
{
  int          k;
  ss_vect<tps> Id, R;

  Id = R.identity();
  for (k = 0; k < nd_tps; k++) {
    if (k < 2) {
      R[2*k]= cos(2e0*M_PI*nu[k])*Id[2*k] + sin(2e0*M_PI*nu[k])*Id[2*k+1];
      R[2*k+1] = -sin(2e0*M_PI*nu[k])*Id[2*k] + cos(2e0*M_PI*nu[k])*Id[2*k+1];
    } else {
      R[delta_]= cos(2e0*M_PI*nu[k])*Id[delta_] + sin(2e0*M_PI*nu[k])*Id[ct_];
      R[ct_] = -sin(2e0*M_PI*nu[k])*Id[delta_] + cos(2e0*M_PI*nu[k])*Id[ct_];
    }
  }
  return R;
}


void compute_M_3D(const double phi0, const ss_vect<tps> &M_tau)
{
  double          dnu[3];
  ss_vect<double> eta;
  ss_vect<tps>    M, A, R, A0, A1, M_cav;

  globval.radiation = false;
  globval.Cavity_on = true;
  Ring_GetTwiss(true, 0e0);
  M = putlinmat(2*nd_tps, globval.OneTurnMat);
  A = get_A_CS(nd_tps, putlinmat(2*nd_tps, globval.Ascr), dnu);

  eta = compute_eta(M);

  prt_map(nd_tps, "\nM:", M);
  prt_map(nd_tps, "\nA:", A);

  globval.alpha_z =
    -globval.Ascr[ct_][ct_]*globval.Ascr[delta_][ct_]
    - globval.Ascr[ct_][delta_]*globval.Ascr[delta_][delta_];
  globval.beta_z = sqr(globval.Ascr[ct_][ct_]) + sqr(globval.Ascr[ct_][delta_]);

  const double
    Circ    = Cell[globval.Cell_nLoc].S,
    alpha_c = globval.Alphac,
    alpha[] = {Cell[0].Alpha[X_], Cell[0].Alpha[Y_], globval.alpha_z},
    beta[]  = {Cell[0].Beta[X_], Cell[0].Beta[Y_], globval.beta_z},
    nu[]    =
    {globval.TotalTune[X_], globval.TotalTune[Y_], globval.Omega};

#if 1
  A = compute_A_3D(alpha, beta);
  R = compute_R_3D(nu);
  // R = M_tau*R;
  M = A*R*Inv(A);
#else
  compute_A_2D(eta, alpha, beta, A0, A1);
  A = A0*A1;
  R = compute_R_2D(Circ, alpha_c, nu, A);
  R = M_tau*R;
  M = A*R*Inv(A);

  globval.Cavity_on = true;
  M_cav = compute_M_cav("cav", phi0);
  printf("\nphi0 [deg] = 180 - %4.2f\n", fabs(phi0*180e0/M_PI));
  M = M*M_cav;
#endif

  prt_map(nd_tps, "\nM:", M);
}


ss_vect<tps> compute_transp(const int n_dof, const ss_vect<tps> &A)
{
  // Compute the transpose of a matrix.
  Matrix A_mat;

  getlinmat(2*n_dof, A, A_mat);
  TpMat(2*n_dof, A_mat);
  return putlinmat(2*n_dof, A_mat);
}


double compute_trace(const int n_dof, const ss_vect<tps> &A)
{
  int    k;
  double trace;

  trace = 0e0;
  for (k = 0; k < 2*n_dof; k++)
    trace += A[k][k];
  return trace;
}


double compute_det(const int n_dof, const ss_vect<tps> &A)
{
  // Compute the determinant of a matrix.
  Matrix A_mat;

  getlinmat(2*n_dof, A, A_mat);
  return DetMat(2*n_dof, A_mat);
}


ss_vect<tps> compute_sympl_conj(const ss_vect<tps> &A)
{
  // Symplectic conjugate for a 2x2 matrix.
  ss_vect<tps> Id, B;

  Id.identity();
  B.zero();
  B[x_]  =  A[px_][px_]*Id[x_] - A[x_][px_]*Id[px_];
  B[px_] = -A[px_][x_]*Id[x_]  + A[x_][x_]*Id[px_];
  return B;
}


void compute_normal_mode_form(const ss_vect<tps> &T)
{
  int
    j, k;
  double
    cos_mu1_m_cos_mu2, cs_2phi, sn_2phi, cs_phi, sn_phi, sgn_D_det;
  ss_vect<tps>
    Id, M, N, m, n, D, D_inv, R;

  const int plane = 4;

  Id.identity();
  M = N = m = n.zero();
  for (j = 0; j < 2; j++)
    for (k = 0; k < 2; k++) {
      M[j] += T[j][k]*Id[k];
      N[j] += T[j+plane][k+plane]*Id[k];
      m[j] += T[j+plane][k]*Id[k];
      n[j] += T[j][k+plane]*Id[k];
    }

  if (!false) {
    prt_map(nd_tps, "\nT:", T);
    prt_map(nd_tps, "\nM:", M);
    prt_map(nd_tps, "\nN:", N);
    prt_map(nd_tps, "\nm:", m);
    prt_map(nd_tps, "\nn:", n);
    printf("\n|n| = %11.5e\n", compute_det(1, n));
  }

  cos_mu1_m_cos_mu2 =
    1e0/2e0*compute_trace(1, M-N)
    *sqrt(1e0+(2e0*compute_det(1, m)+compute_trace(1, n*m))
	  /sqr(1e0/2e0*compute_trace(1, M-N)));
  printf("\n  cos(mu1) - cos(mu2) = %11.5e\n", cos_mu1_m_cos_mu2);

  cs_2phi = 1e0/2e0*compute_trace(1, M-N)/cos_mu1_m_cos_mu2;

  if (fabs(cs_2phi) < 1e0) {
    printf("compute_normal_mode_form: not tested!");
    exit(1);

    sn_2phi =
      fabs(sqrt(2e0*compute_det(1, m)+compute_trace(1, n*m))/cos_mu1_m_cos_mu2);
    cs_phi = sqrt((1e0+cs_2phi)/2e0);
    sn_phi = sn_2phi/sqrt(2e0*(1e0-cs_2phi));
    printf("\n  cos(2 phi) = %11.5e\n", cs_2phi);
    printf("  sin(2 phi) = %11.5e\n", sn_2phi);
    printf("  sin(2 phi) = %11.5e\n", sin(acos(cs_2phi)));
  } else {
    sn_2phi =
      fabs(sqrt(-(2e0*compute_det(1, m)+compute_trace(1, n*m)))
	   /cos_mu1_m_cos_mu2);
    cs_phi = sqrt((1e0+cs_2phi)/2e0);
    sn_phi = sn_2phi/sqrt(2e0*(1e0+cs_2phi));
    printf("\n  cosh(2 phi) = %11.5e\n", cs_2phi);
    printf("  sinh(2 phi) = %11.5e\n", sn_2phi);
  }

  D =
    -(m+compute_transp(1, get_S(1))*compute_transp(1, n)*get_S(1));
  // Supported operators issue.
  D *= 1e0/(cos_mu1_m_cos_mu2*sn_2phi);
  for (k = y_; k < 2*nd_tps; k++)
    D[k]= Id[k];
  D_inv = Inv(D);
  sgn_D_det = compute_det(1, D);
  printf("\n  |D| = %11.5e\n", sgn_D_det);

  prt_map(nd_tps, "\nD:", D);

  R.identity();
  R[x_]  =
    Id[x_]*cs_phi
    + sgn_D_det*(D_inv[x_][x_]*Id[plane]+D_inv[x_][px_]*Id[plane+1])*sn_phi;
  R[px_] =
    Id[px_]*cs_phi
    + sgn_D_det*(D_inv[px_][x_]*Id[plane]+D_inv[px_][px_]*Id[plane+1])*sn_phi;
  R[plane]  =
    -(D[x_][x_]*Id[x_]+D[x_][px_]*Id[px_])*sn_phi + Id[plane]*cs_phi;
  R[plane+1] =
    -(D[px_][x_]*Id[x_]+D[px_][px_]*Id[px_])*sn_phi + Id[plane+1]*cs_phi;

  prt_map(nd_tps, "\nR:", R);
  prt_map(nd_tps, "\nR^-1.T.R:", Inv(R)*T*R);
  prt_map(nd_tps, "\nT:", T);

  const ss_vect<tps>
    A = M + sgn_D_det*Inv(D)*m*(sn_phi/cs_phi),
    B = N + D*n*(sn_phi/cs_phi);
  M = A*sqr(cs_phi) + sgn_D_det*Inv(D)*B*D*sqr(sn_phi);
  N = B*sqr(cs_phi) + sgn_D_det*D*A*Inv(D)*B*D*sqr(sn_phi);
  m = -(D*A-B*D)*sn_phi*cs_phi;
  n = -(sgn_D_det*A*Inv(D)-sgn_D_det*Inv(D)*B)*sn_phi*cs_phi;
  prt_map(nd_tps, "\nA:", A);
  prt_map(nd_tps, "\nB:", B);
  prt_map(nd_tps, "\nM:", M);
  prt_map(nd_tps, "\nN:", N);
  prt_map(nd_tps, "\nm:", m);
  prt_map(nd_tps, "\nn:", n);
}


void compute_maps(void)
{
  long int
    lastpos;
  double
    U0, nu[3], alpha_rad[3], dnu[3], phi0, dummy[3], D[3], tau[3], eps[3],
    sigma_s, sigma_delta;
  ss_vect<double>
    fixed_point, fixed_point_no_rad, eta;
  ss_vect<tps>
    M_3D, A_3D, M_3D_no_rad, A_3D_no_rad, M_tau, M_cav, M_2D, M, D_mat, M_Chol,
    A0;
   CavityType
     *C;

  const string
    file_name = "compute_maps",
    cav_name  = "cav";
  const int
    loc       = Elem_GetPos(ElemIndex(cav_name.c_str()), 1);
  const double
    Circ      = Cell[globval.Cell_nLoc].S,
    E0        = 1e9*globval.Energy;

  C = Cell[loc].Elem.C;

  globval.radiation = false;
  globval.Cavity_on = true;

  C->phi_RF = 0e0;
  getcod(0e0, lastpos);
  fixed_point = globval.CODvect;
  cout << scientific << setprecision(3)
       << "\nFixed point =" << setw(11) << fixed_point << "\n";
  U0 = globval.dE*E0;
  M_3D = compute_M(fixed_point);
  prt_map(nd_tps, "\nM_3D:", M_3D);
  A_3D = compute_A(Circ, M_3D, nu, alpha_rad, true);
  A_3D = get_A_CS(nd_tps, A_3D, dnu);
  prt_map(nd_tps, "\nA_3D:", A_3D);

  prt_map(nd_tps, "\nR = Inv(A_3D)*M_3D*A_3D:", Inv(A_3D)*M_3D*A_3D);
  prt_map(file_name+"_M_3D.dat", M_3D);
  prt_map(file_name+"_A_3D.dat", A_3D);

  compute_normal_mode_form(M_3D);

  exit(0);

  if (!false)
    compute_M_2D();

  M_tau = compute_M_tau(alpha_rad);
  prt_map(nd_tps, "\nM_tau:", M_tau);
  prt_map(file_name+"_M_tau.dat", M_tau);

  phi0 = asin(U0/C->V_RF);
  C->phi_RF = phi0;
  compute_M_3D(phi0, M_tau);

  exit(0);

  globval.radiation = false;
  globval.Cavity_on = true;

  phi0 = asin(U0/C->V_RF);
  C->phi_RF = phi0;
  printf("\nphi0 [deg] = 180 - %4.2f\n", fabs(phi0*180e0/M_PI));
  getcod(0e0, lastpos);
  fixed_point_no_rad = globval.CODvect;
  cout << scientific << setprecision(3)
       << "\nFixed point =" << setw(11) << fixed_point_no_rad << "\n";
  M_3D_no_rad = compute_M(fixed_point_no_rad);
  prt_map(nd_tps, "\nM_3D_no_rad:", M_3D_no_rad);
  A_3D_no_rad = compute_A(Circ, M_3D_no_rad, nu, dummy, true);
  A_3D_no_rad = get_A_CS(nd_tps, A_3D_no_rad, dnu);
  prt_map(nd_tps, "\nA_3D_no_rad:", A_3D_no_rad);

  M_3D_no_rad = A_3D_no_rad*Inv(M_tau)*Inv(A_3D)*M_3D*A_3D*Inv(A_3D_no_rad);
  prt_map(nd_tps,
	  "\nM_3D_no_rad = "
	  "A_3D_no_rad.M_tau^-1.A_3D^-1.M_3D*A_3D.A_3D_no_rad^-1:",
	  M_3D_no_rad);
  prt_map(nd_tps, "\nM_3D_no_rad^t.Omega.M_3D_no_rad - Omega:",
	  tp_map(3, M_3D_no_rad)*get_S(3)*M_3D_no_rad-get_S(3));

  exit(0);

  M_3D_no_rad = A_3D_no_rad*Inv(M_tau)*Inv(A_3D)*M_3D*A_3D*Inv(A_3D_no_rad);
  prt_map(nd_tps,
	  "\nM_3D_no_rad = "
	  "A_3D_no_rad.M_tau^-1.A_3D^-1.M_3D*A_3D.A_3D_no_rad^-1:",
	  M_3D_no_rad);
  prt_map(nd_tps, "\nM_3D_no_rad^t.Omega.M_3D_no_rad - Omega:",
	  tp_map(3, M_3D_no_rad)*get_S(3)*M_3D_no_rad-get_S(3));
  compute_A(Circ, M_3D_no_rad, nu, dummy, true);
  prt_map(file_name+"_A_3D_no_rad.dat", A_3D_no_rad);

  phi0 = asin(U0/C->V_RF);
  M_cav = compute_M_cav("cav", phi0);
  printf("\nphi0 [deg] = 180 - %4.2f\n", fabs(phi0*180e0/M_PI));
  prt_map(nd_tps, "\nM_cav:", M_cav);

  M_2D = M_3D_no_rad*Inv(M_cav);
  prt_map(nd_tps, "\nM_2D = M_3D_no_rad.M_cav^-1:", M_2D);
  prt_map(nd_tps, "\nM_2D^t.Omega.M_2D - Omega:",
	  tp_map(3, M_2D)*get_S(3)*M_2D-get_S(3));
  // Compute nu.
  compute_A(Circ, M_2D, nu, dummy, false);
  prt_map(file_name+"_M_2D.dat", M_2D);

  globval.radiation = globval.Cavity_on = false;

  M = compute_M();
  prt_map(nd_tps, "\nM_2D:", M);
 // Compute nu.
  compute_A(Circ, M, nu, dummy, false);

  if (false) {
    // Remove linear dispersion.
    printf("\nRemoving linear dispersion:\n");
    eta = compute_eta(M);
    A0 = compute_A0(eta);
    M = Inv(A0)*M*A0;
    prt_map(nd_tps, "\nM:", M);
  }

  if (!false)
    compute_map
      (M_2D, M_cav, M_3D_no_rad, A_3D_no_rad, M_tau, M_3D, A_3D, alpha_rad, U0);

  compute_D(fixed_point, A_3D, D);
  compute_tau(Circ, alpha_rad, tau);
  compute_eps(Circ, tau, D, eps);
  compute_bunch_size(A_3D, eps, sigma_s, sigma_delta);
  D_mat = compute_D_mat(A_3D, D);
  prt_map(nd_tps, "\nDiffusion Matrix:", D_mat);

  M_Chol = compute_M_Chol(D_mat);
  prt_map(nd_tps, "\nCholesky Decomposition:", M_Chol);
  prt_map(file_name+"_M_Chol.dat", M_Chol);

  globval.CODvect[ct_] /= c0;
  cout << scientific << setprecision(3)
       << "\nFixed point         = " << setw(11) << globval.CODvect << "\n";
  globval.CODvect[ct_] *= c0;
  cout << scientific << setprecision(3)
       << "eta                 ="  << setw(11) << eta << "\n";
  printf("U0 [keV]            = %3.1f\n", 1e-3*U0);
  printf("U0/E0               = %10.3e\n", U0/E0);
  printf("phi0 [deg]          = 180 - %4.2f\n", fabs(phi0*180e0/M_PI));
  printf("tau [msec]          = [%5.3f, %5.3f, %5.3f]\n",
	 1e3*tau[X_], 1e3*tau[Y_], 1e3*tau[Z_]);
  printf("D                   = [%9.3e, %9.3e, %9.3e]\n",
	 D[X_], D[Y_], D[Z_]);
  printf("eps                 = [%9.3e, %9.3e, %9.3e]\n",
	 eps[X_], eps[Y_], eps[Z_]);
  printf("sigma_x [micro m]   = %5.3f\n", 1e6*sqrt(eps[X_]*Cell[0].Beta[X_]));
  printf("sigma_x'[micro rad] = %5.3f\n",
	 1e6*sqrt(eps[X_]*(1e0+sqr(Cell[0].Alpha[X_]))/Cell[0].Beta[X_]));
  printf("sigma_t [ps]        = %5.3f\n", 1e12*sigma_s/c0);
  printf("sigma_delta         = %9.3e\n", sigma_delta);
}


void set_state(void)
{
  globval.H_exact        = false;
  globval.quad_fringe    = false;
  globval.Cavity_on      = false;
  globval.radiation      = false;
  globval.emittance      = false;
  globval.IBS            = false;
  globval.pathlength     = false;
  globval.Aperture_on    = false;
  globval.Cart_Bend      = false;
  globval.dip_edge_fudge = true;
}


int main(int argc, char *argv[])
{

  reverse_elem     = true;
  globval.mat_meth = false;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  set_state();

  Ring_GetTwiss(true, 0e0);

  printglob();
  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);
  prtmfile("flat_file.dat");

  daeps_(1e-15);
  compute_maps();
}
