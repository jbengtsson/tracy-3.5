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


void prt_map(const string file_name, const ss_vect<tps> map)
{
  ofstream outf;

  file_wr(outf, file_name.c_str());
  outf << scientific << setprecision(3) << setw(11) << map;
  outf.close();
}


ss_vect<tps> compute_transp(const int n_dof, const ss_vect<tps> &A)
{
  // Matrix transpose.
  Matrix A_mat;

  getlinmat(2*n_dof, A, A_mat);
  TpMat(2*n_dof, A_mat);
  return putlinmat(2*n_dof, A_mat);
}


double compute_trace(const int n_dof, const ss_vect<tps> &A)
{
  // Matrix trace.
  int    k;
  double trace;

  trace = 0e0;
  for (k = 0; k < 2*n_dof; k++)
    trace += A[k][k];
  return trace;
}


double compute_det(const int n_dof, const ss_vect<tps> &A)
{
  // Matrix determinant.
  Matrix A_mat;

  getlinmat(2*n_dof, A, A_mat);
  return DetMat(2*n_dof, A_mat);
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
  // Compute the Poincaré map for the fixed point.
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
  // Canonical Transfomation to delta dependent fixed point.
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
(const double Circ, const double alpha_rad[], double tau[])
{
  // Compute the damping times.
  
  for (int k = 0; k < nd_tps; k++)
    tau[k] = -Circ/(c0*alpha_rad[k]);
}


void compute_eps
(const double Circ, const double tau[], const double D[], double eps[])
{
  // Compute the eigen emittances.

  for (int k = 0; k < nd_tps; k++)
    eps[k] = D[k]*tau[k]*c0/(2e0*Circ);
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

void compute_Twiss_long(void)
{
  const ss_vect<tps> A = putlinmat(2*nd_tps, globval.Ascr);

  globval.alpha_z =
    -A[ct_][ct_]*A[delta_][ct_] - A[ct_][delta_]*A[delta_][delta_];
  globval.beta_z =
    sqr(A[ct_][ct_]) + sqr(A[ct_][delta_]);
}


void compute_bunch_size
(const double eps[], double &sigma_s, double &sigma_delta)
{
  // Compute the bunch size: sigma_s & sigma_delta.
  double gamma_s;

  compute_Twiss_long();
  gamma_s = (1e0+sqr(globval.alpha_z))/globval.beta_z;

  sigma_s = sqrt(globval.beta_z*eps[Z_]);
  sigma_delta = sqrt(gamma_s*eps[Z_]);
}


ss_vect<tps> compute_D_mat(const ss_vect<tps> &A, const double D[])
{
  // Compute the diffusion matrix.
  ss_vect<tps> D_diag;

  D_diag.zero();
  for (int k = 0; k < 2*nd_tps; k++)
    D_diag[k] = D[k/2]*tps(0e0, k+1);
  return A*D_diag*compute_transp(nd_tps, A);
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

  return compute_transp(3, M_Chol_t);
}


ss_vect<tps> compute_M_cav(const string &cav_name, double phi_RF)
{
  tps          ct0;
  ss_vect<tps> M;
  CavityType   *C;

  const int loc = Elem_GetPos(ElemIndex(cav_name.c_str()), 1);

  C = Cell[loc].Elem.C;

  ct0 = tps(0e0, ct_+1);
  M.identity();
#if 1
  C->phi_RF = phi_RF;
  M[ct_] += ct0.cst();
  Cav_Pass(Cell[loc], M);
#else
  M[delta_] +=
    -C->V_RF/(1e9*globval.Energy)*sin(2e0*M_PI*C->f_RF*ct0/c0+phi_RF);
#endif
  M[delta_] -= M[delta_].cst();
  return M;
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


void get_Twiss(const int n_dof, double alpha[], double beta[], double nu[])
{
  int k;

  for (k = 0; k < n_dof; k++) {
    if (k < 2) {
      alpha[k] = Cell[k].Alpha[k];
      beta[k] = Cell[k].Beta[k];
      nu[k] = globval.TotalTune[k];
    } else {
      alpha[k] = globval.alpha_z;
      beta[k] = globval.beta_z;
      nu[k] = globval.Omega;
    }
  }
}


ss_vect<tps> compute_M_2D(void)
{
  double          alpha[2], beta[2], nu[2];
  ss_vect<double> eta;
  ss_vect<tps>    M, A, A0, A1, R;

  const double
    Circ    = Cell[globval.Cell_nLoc].S,
    alpha_c = globval.Alphac;

  globval.radiation = globval.Cavity_on = false;

  Ring_GetTwiss(true, 0e0);
  M = putlinmat(2*nd_tps, globval.OneTurnMat);
  prt_map(nd_tps, "\nM:", M);

  get_Twiss(2, alpha, beta, nu);

  eta.zero();
  eta[x_] = Cell[0].Eta[X_];
  eta[px_] = Cell[0].Etap[X_];

  compute_A_2D(eta, alpha, beta, A0, A1);
  A = A0*A1;
  R = compute_R_2D(Circ, alpha_c, nu, A);
  M = A*R*Inv(A);
  prt_map(nd_tps, "\nM:", M);

  return M;
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


ss_vect<tps> compute_normal_mode_form(const ss_vect<tps> &T)
{
  int
    j, k;
  double
    cos_mu1_m_cos_mu2, cs_2phi, sn_2phi, cs_phi, sn_phi, sgn_D_det;
  ss_vect<tps>
    Id, M, N, m, n, D, D_inv, R, A, B;

  const bool prt  = false;
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

  if (prt)
    prt_map(nd_tps, "\nT:", T);

  cos_mu1_m_cos_mu2 =
    1e0/2e0*compute_trace(1, M-N)
    *sqrt(1e0+(2e0*compute_det(1, m)+compute_trace(1, n*m))
	  /sqr(1e0/2e0*compute_trace(1, M-N)));
  printf("\n  cos(mu1) - cos(mu2) = %12.5e\n", cos_mu1_m_cos_mu2);

  cs_2phi = 1e0/2e0*compute_trace(1, M-N)/cos_mu1_m_cos_mu2;

  if (fabs(cs_2phi) < 1e0) {
    printf("\n*** compute_normal_mode_form: not tested!\n");
    exit(1);

    sn_2phi =
      fabs(sqrt(2e0*compute_det(1, m)+compute_trace(1, n*m))/cos_mu1_m_cos_mu2);
    cs_phi = sqrt((1e0+cs_2phi)/2e0);
    sn_phi = sn_2phi/sqrt(2e0*(1e0-cs_2phi));
    printf("  cos(2 phi)          = %12.5e\n", cs_2phi);
    printf("  sin(2 phi)          = %12.5e\n", sn_2phi);
  } else {
    sn_2phi =
      fabs(sqrt(-(2e0*compute_det(1, m)+compute_trace(1, n*m)))
	   /cos_mu1_m_cos_mu2);
    cs_phi = sqrt((1e0+cs_2phi)/2e0);
    sn_phi = sn_2phi/sqrt(2e0*(1e0+cs_2phi));
    printf("  cosh(2 phi)         = %12.5e\n", cs_2phi);
    printf("  sinh(2 phi)         = %12.5e\n", sn_2phi);
  }

  D =
    -(m+compute_transp(1, get_S(1))*compute_transp(1, n)*get_S(1));
  // Supported operators issue.
  D *= 1e0/(cos_mu1_m_cos_mu2*sn_2phi);
  for (k = y_; k < 2*nd_tps; k++)
    D[k]= Id[k];
  D_inv = Inv(D);
  sgn_D_det = compute_det(1, D);
  printf("  |D|                 = %12.5e\n", sgn_D_det);

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

  if (prt) {
    prt_map(nd_tps, "\nT:", T);
    prt_map(nd_tps, "\nR:", R);
    prt_map(nd_tps, "\nR^-1.T.R:", Inv(R)*T*R);

    A = M + sgn_D_det*Inv(D)*m*(sn_phi/cs_phi);
    B = N + D*n*(sn_phi/cs_phi);
    M = A*sqr(cs_phi) + sgn_D_det*Inv(D)*B*D*sqr(sn_phi);
    N = B*sqr(cs_phi) + sgn_D_det*D*A*Inv(D)*B*D*sqr(sn_phi);
    m = -(D*A-B*D)*sn_phi*cs_phi;
    n = -sgn_D_det*(A*Inv(D)-Inv(D)*B)*sn_phi*cs_phi;

    prt_map(nd_tps, "\nD:", D);
    prt_map(nd_tps, "\nA:", A);
    prt_map(nd_tps, "\nB:", B);
    prt_map(nd_tps, "\nM:", M);
    prt_map(nd_tps, "\nN:", N);
    prt_map(nd_tps, "\nm:", m);
    prt_map(nd_tps, "\nn:", n);
  }

  return R;
}


void compute_M_3D
(const string &file_name, const string &cav_name, double &phi_RF,
 ss_vect<double> &fixed_point, double &U0, double alpha_rad[],
 ss_vect<tps> &A_3D)
{
  long int        lastpos;
  int             k;
  double          alpha[3], beta[3], nu[3];
  ss_vect<double> eta;
  ss_vect<tps>    M_3D_no_rad, M_3D, A_CS, A_sb, M_tau, R, M;
  CavityType      *C;

  const int
    loc = Elem_GetPos(ElemIndex(cav_name.c_str()), 1);
  const double
    E0  = 1e9*globval.Energy;

  C = Cell[loc].Elem.C;

  // Compute radiation effects.
  globval.radiation = globval.Cavity_on = true;

  getcod(0e0, lastpos);
  fixed_point = globval.CODvect;
  U0 = globval.dE*E0;
  phi_RF = asin(U0/C->V_RF);
  C->phi_RF = phi_RF;
  printf("\nphi_RF [deg] = 180 - %4.2f\n", fabs(phi_RF*180e0/M_PI));

  // With RF cavity but without radiation.
  globval.radiation = false;
  globval.Cavity_on = true;

  Ring_GetTwiss(true, 0e0);
  compute_Twiss_long();
  M_3D_no_rad = putlinmat(2*nd_tps, globval.OneTurnMat);

  get_Twiss(nd_tps, alpha, beta, nu);

  A_CS = compute_A_3D(alpha, beta);
  A_sb = compute_normal_mode_form(M_3D_no_rad);

  R = compute_R_3D(nu);
  M = A_sb*A_CS*R*Inv(A_sb*A_CS);
  prt_map(nd_tps, "\nM - with RF cavity & no radiation:", M);
  prt_map(nd_tps, "\nValidation:", M-M_3D_no_rad);

  // With RF cavity & radiation.
  globval.radiation = globval.Cavity_on = true;

  C->phi_RF = 0e0;
  printf("\nphi_RF [deg] = 180 - %4.2f\n", fabs(C->phi_RF*180e0/M_PI));
  Ring_GetTwiss(true, 0e0);
  for (k = 0; k < nd_tps; k++)
    alpha_rad[k] = globval.alpha_rad[k];
  compute_Twiss_long();
  M_3D = putlinmat(2*nd_tps, globval.OneTurnMat);

  get_Twiss(nd_tps, alpha, beta, nu);

  M_tau = compute_M_tau(globval.alpha_rad);

  A_CS = compute_A_3D(alpha, beta);
  A_sb = compute_normal_mode_form(M_3D_no_rad);
  A_3D = A_sb*A_CS;

  R = compute_R_3D(nu);
  M = A_3D*M_tau*R*Inv(A_3D);
  prt_map(nd_tps, "\nM - w/ RF cavity & radiation:", M);
  prt_map(nd_tps, "\nValidation:", M-M_3D);

  prt_map(file_name+"_R.dat", R);
  prt_map(file_name+"_M_tau.dat", M_tau);
  prt_map(file_name+"_A_CS.dat", A_CS);
  prt_map(file_name+"_A_sb.dat", A_sb);
}


void prt_summary
(ss_vect<double> &fixed_point, const double U0, const double phi_RF,
 const double tau[], const double D[], const double eps[], const double sigma_s,
 const double sigma_delta)
{
  int k;

  const double E0 = 1e9*globval.Energy;

  fixed_point[ct_] /= c0;
  cout << scientific << setprecision(3)
       << "\n  Fixed point         =" << setw(11) << fixed_point << "\n";
  fixed_point[ct_] *= c0;
  printf("  U0 [keV]            = %3.1f\n", 1e-3*U0);
  printf("  U0/E0               = %10.3e\n", U0/E0);
  printf("  phi_RF [deg]        = 180 - %4.2f\n", fabs(phi_RF*180e0/M_PI));
  printf("  tau [msec]          = [");
  for (k = 0; k < nd_tps; k++)
    printf("%5.3f%s", 1e3*tau[k], (k < 2)? ", " : "]\n");
  printf("  D                   = [");
  for (k = 0; k < nd_tps; k++)
    printf("%9.3e%s", D[k], (k < 2)? ", " : "]\n");
  printf("  eps                 = [");
  for (k = 0; k < nd_tps; k++)
    printf("%9.3e%s", eps[k], (k < 2)? ", " : "]\n");
  printf("  sigma_x [micro m]   = %5.3f\n", 1e6*sqrt(eps[X_]*Cell[0].Beta[X_]));
  printf("  sigma_x'[micro rad] = %5.3f\n",
	 1e6*sqrt(eps[X_]*(1e0+sqr(Cell[0].Alpha[X_]))/Cell[0].Beta[X_]));
  printf("  sigma_t [ps]        = %5.3f\n", 1e12*sigma_s/c0);
  printf("  sigma_delta         = %9.3e\n", sigma_delta);
}


void compute_maps(void)
{
  double
    phi_RF, U0, alpha_rad[3], D[3], tau[3], eps[3], sigma_s, sigma_delta;
  ss_vect<double>
    fixed_point;
  ss_vect<tps>
    M_3D, A_3D, D_mat, M_Chol;

  const string
    file_name = "compute_maps",
    cav_name  = "cav";
  const double
    Circ      = Cell[globval.Cell_nLoc].S;

  compute_M_3D(file_name, cav_name, phi_RF, fixed_point, U0, alpha_rad, A_3D);

  compute_D(fixed_point, A_3D, D);
  compute_tau(Circ, alpha_rad, tau);
  compute_eps(Circ, tau, D, eps);
  compute_bunch_size(eps, sigma_s, sigma_delta);
  D_mat = compute_D_mat(A_3D, D);
  prt_map(nd_tps, "\nDiffusion Matrix:", D_mat);

  M_Chol = compute_M_Chol(D_mat);
  prt_map(nd_tps, "\nCholesky Decomposition:", M_Chol);
  prt_map(file_name+"_M_Chol.dat", M_Chol);

  prt_summary(fixed_point, U0, phi_RF, tau,  D,  eps, sigma_s, sigma_delta);
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
