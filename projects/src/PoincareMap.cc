
// Poincare Map Object.
// Johan Bengtsson 29/03/20.


struct BeamType {
private:
public:
  int                            n_part;
  std::vector< ss_vect<double> > ps;
  ss_vect<double>                mean;
  ss_vect<tps>                   sigma;

  void BeamInit(const int n_part, const double eps_x, const double eps_y,
		const double sigma_delta);
  void BeamStats(void);
  void print(void);
};


struct PoincareMap {
private:
public:
  int n_DOF;           // Degrees-of-Freedom.
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
    delta_rad,         // Momentum change.
    V_RF,              // RF Voltage.
    f_RF,              // RF frequency.
    phi0,              // Synchronous phase.
    tau[3],            // Damping times.
    D[3];              // Diffusion Coefficients.
  ss_vect<double>
    ps_cod,            // Fixed Point.
    eta1;              // Linear Dispersion.
  ss_vect<tps>
    map_num,           // numerical Poincare Map.
    M_lat,             // Linear Map for Lattice.
    M_cav,             // Linear Map for RF Cavity.
    A,                 // Transformation to Floquet Space.
    R,                 // Phase Space Rotation.
    M_delta,           // Radiation Loss Map.
    M_tau,             // Radiation Damping Map.
    M,                 // Poincare Map.
    M_diff, M_Chol;    // Diffusion Matrix & Cholesky Decomposition.

  void GetMap(void);
  void GetM_lat(void);
  void GetDisp(void);
  void GetA(void);
  void GetM_cav(void);
  void GetM_delta(void);
  void GetM_tau(void);
  void GetM_diff(void);
  void GetM_Chol(void);
  void GetM_rad(void);
  void GetM(const bool cav, const bool rad);
  void propagate(const int n, BeamType &beam);
  void print(void);
};

void PrtMap(const string &str, const double n_DOF, ss_vect<tps> map)
{
  int i, j;

  const int    n_dec = 6;
  const double eps   = 1e-20;

  printf("%s\n", str.c_str());
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


double DetMap(const int n, ss_vect<tps> A)
{
  Matrix A_mat;

  getlinmat(2*n, A, A_mat);
  return DetMat(2*n, A_mat);
}


ss_vect<tps> TpMap(const int n, ss_vect<tps> A)
{
  Matrix A_mat;

  getlinmat(2*n, A, A_mat);
  TpMat(2*n, A_mat);
  return putlinmat(2*n, A_mat);
}


void PoincareMap::GetM_lat(void)
{
  int          k;
  double       mu[2], gamma[2];
  ss_vect<tps> Id;

  Id.identity();

  for (k = 0; k < 2; k++) {
    mu[k] = 2e0*M_PI*nu[k];
    gamma[k] = (1e0+sqr(alpha[k]))/beta[k];
  }

  M_lat.identity();
  for (k = 0; k < 2; k++) {
    M_lat[2*k] =
      (cos(mu[k])+alpha[k]*sin(mu[k]))*Id[2*k] + beta[k]*sin(mu[k])*Id[2*k+1];
    M_lat[2*k+1] =
      -gamma[k]*sin(mu[k])*Id[2*k]
      + (cos(mu[k])-alpha[k]*sin(mu[k]))*Id[2*k+1];
  }
  M_lat[ct_] += alpha_c*C*Id[delta_];

  M_lat[x_] += globval.OneTurnMat[x_][delta_]* Id[delta_];
  M_lat[px_] += globval.OneTurnMat[px_][delta_]*Id[delta_];
  M_lat[ct_] +=
    (M_lat[x_][x_]*M_lat[px_][delta_]-M_lat[px_][x_]*M_lat[x_][delta_])*Id[x_]
    + (M_lat[x_][px_]*M_lat[px_][delta_]-M_lat[px_][px_]*M_lat[x_][delta_])
    *Id[px_];
}


void PoincareMap::GetDisp(void)
{
  long int     lastpos, jj[6];
  int          k;
  ss_vect<tps> Id, map;

  Id.identity();

  globval.Cavity_on = false; globval.radiation = false;
  getcod(0e0, lastpos);

  map = putlinmat(6, globval.OneTurnMat);
  eta1.zero(); eta1[x_] = map[x_][delta_]; eta1[px_] = map[px_][delta_];
  for (k = 0; k < 6; k++)
    jj[k] = 0;
  jj[x_] = 1; jj[px_] = 1;
  eta1 = (PInv(Id-map, jj)*eta1).cst();
}


void GetTwiss(ss_vect<tps> A, double alpha[], double beta[])
{
  int k;

  for (k = 0; k < 3; k++) {
    alpha[k] = -(A[2*k][2*k]*A[2*k+1][2*k]+A[2*k][2*k+1]*A[2*k+1][2*k+1]);
    beta[k] = sqr(A[2*k][2*k]) + sqr(A[2*k][2*k+1]);
  }
}


void PoincareMap::GetA(void)
{
  int          k;
  double       a_c;
  ss_vect<tps> A1;
  Matrix       M_mat, A_mat, A_inv_mat, R_mat;

  getlinmat(6, map_num, M_mat);
  GDiag(2*n_DOF, C, A_mat, A_inv_mat, R_mat, M_mat, nu[Z_], a_c);

  for (k = 0; k < n_DOF; k++)
    tau[k] = -C/(c0*globval.alpha_rad[k]);

  A = putlinmat(2*n_DOF, A_mat);
  if (n_DOF < 3) {
    A[ct_] = tps(0e0, ct_+1); A[delta_] = tps(0e0, delta_+1);
  }
  R = putlinmat(2*n_DOF, R_mat);

  GetTwiss(A, alpha, beta);

  // GetDisp();
}


void PoincareMap::GetM_cav(void)
{
  const int loc = Elem_GetPos(ElemIndex("cav"), 1);

  V_RF = Cell[loc].Elem.C->Pvolt;
  f_RF = Cell[loc].Elem.C->Pfreq;
  phi0 = asin(U0/V_RF);

  M_cav.identity();
  M_cav[delta_] -=
    V_RF*2e0*M_PI*f_RF*cos(phi0)/(1e9*globval.Energy*c0)*tps(0e0, ct_+1);
}


void PoincareMap::GetM_delta(void)
{
  ss_vect<tps> Id;

  Id.identity();

  M_delta.identity();
  M_delta[ct_] /= 1e0 + delta_rad;
  M_delta[delta_] *= 1e0 + delta_rad;
}


void PoincareMap::GetM_tau(void)
{
  int          k;
  ss_vect<tps> Id;

  Id.identity();

  M_tau.zero();
  for (k = 0; k < n_DOF; k++) {
    M_tau[2*k] = exp(-C/(c0*tau[k]))*Id[2*k];
    M_tau[2*k+1] = exp(-C/(c0*tau[k]))*Id[2*k+1];
  }
}


ss_vect<tps> Get_A_A_tp(const int n, const double alpha[], const double beta[])
{
  int          k;
  double       gamma[3];
  ss_vect<tps> Id, A_A_tp;

  Id.identity();

  for (k = 0; k < n; k++)
    gamma[k] = (1e0+sqr(alpha[k]))/beta[k];
  A_A_tp.zero();
  for (k = 0; k < n; k++) {
    A_A_tp[2*k] += beta[k]*Id[2*k] - alpha[k]*Id[2*k+1];
    A_A_tp[2*k+1] += -alpha[k]*Id[2*k] + gamma[k]*Id[2*k+1];
  }
  return A_A_tp;
}


void PoincareMap::GetM_diff(void)
{
  int long     lastpos;
  int          k;
  ss_vect<tps> As, A_A_tp;

  globval.Cavity_on = true; globval.radiation = true;
  globval.emittance = true;

  As = A + ps_cod;
  Cell_Pass(0, globval.Cell_nLoc, As, lastpos);

  globval.emittance = false;

  for (k = 0; k < n_DOF; k++)
    D[k] = globval.D_rad[k];

  M_diff.identity();
  for (k = 0; k < 2*n_DOF; k++)
    // Diff. Coeffs. are for the Actions.
    M_diff[k] *= 2e0*D[k/2];
  // A_A_tp = Get_A_A_tp(n_DOF, alpha, beta);
  // PrtMap("A_A_tp:", n_DOF, A_A_tp);
  // M_diff = M_diff*A_A_tp;
  M_diff = M_diff*A*TpMap(n_DOF, A);

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
  M_Chol.zero();
  for (j = 1; j <= 4; j++) {
     j1 = (j < 3)? j : j+2;
     for (k = 1; k <= 4; k++) {
       k1 = (k < 3)? k : k+2;
       M_Chol[j1-1] += d2[j][k]*tps(0e0, k1);
     }
  }

  free_dvector(diag, 1, n); free_dmatrix(d1, 1, n, 1, n);
  free_dmatrix(d2, 1, n, 1, n);
}


void PoincareMap::GetMap(void)
{
  int long lastpos;

  globval.Cavity_on = cav_on; globval.radiation = rad_on;

  getcod(0e0, lastpos);

  map_num = putlinmat(2*n_DOF, globval.OneTurnMat);
  ps_cod = globval.CODvect;
  printf("\nDet{M}-1 = %10.3e\n", DetMap(n_DOF, map)-1e0);

  U0 = globval.dE*E0;
  delta_rad = U0/E0;
  alpha_c = map_num[ct_][delta_]/C;
}


void PoincareMap::GetM(const bool cav, const bool rad)
{
  int k;

  cav_on = cav; rad_on = rad;

  C = Cell[globval.Cell_nLoc].S; E0 = 1e9*globval.Energy;
  n_DOF = (!cav_on)? 2 : 3;

  GetMap();
  GetA();
  GetM_tau();
  R = Inv(M_tau)*R;
  for (k = 0; k < n_DOF; k++)
    nu[k] = atan2(R[2*k][2*k+1], R[2*k][2*k])/(2e0*M_PI);

  GetM_lat();
  GetM_cav();
  M = M_cav*M_lat;
  if (rad) {
    GetM_delta();
    M = M_tau*M_delta*M;
    GetM_diff();
  }
}


void PoincareMap::propagate(const int n, BeamType &beam)
{
  int             i, j, k;
  ss_vect<double> X;
  ofstream        outf;

  const int    n_prt     = 10;
  const string file_name = "propagate.out";

  file_wr(outf, file_name.c_str());

  for (i = 1; i <= n; i++) {
    for (j = 0; j < (int)beam.ps.size(); j++) {
      for (k = 0; k < 6; k++)
	X[k] = normranf();
      beam.ps[j] = (M_lat*beam.ps[j]+M_Chol*X).cst();
    }

    beam.BeamStats();

    if (i % n_prt == 0) {
      outf << setw(7) << i;
      for (j = 0; j < 6; j++)
	outf << scientific << setprecision(5) << setw(13) << beam.sigma[j][j];
      outf << "\n";
    }
  }
  if (n % n_prt != 0) {
    outf << setw(7) << n;
    for (j = 0; j < 6; j++)
      outf << scientific << setprecision(5) << setw(13) << beam.sigma[j][j];
    outf << "\n";
  }

  outf.close();
}


void PoincareMap::print(void)
{
  printf("\nCavity %d, Radiation %d\n", cav_on, rad_on);
  printf("\nC [m]      = %7.5f\n", C);

  printf("alpha      = [%7.5f, %7.5f]\n", alpha[X_], alpha[Y_], alpha[Z_]);
  printf("beta       = [%7.5f, %7.5f]\n", beta[X_], beta[Y_], beta[Z_]);
  printf("eta        = [%7.5f, %7.5f]\n", eta1[x_], eta1[px_]);
  printf("nu         = [%7.5f, %7.5f, %7.5f]\n", nu[X_], nu[Y_], nu[Z_]);

  printf("V_RF       = %7.5e\n", V_RF);
  printf("f_RF [MhZ] = %7.5f\n", 1e-6*f_RF);
  printf("phi0       = %7.5f\n", phi0*180.0/M_PI);

  printf("E0 [GeV]   = %7.5f\n", 1e-9*E0);
  printf("U0 [keV]   = %7.5f\n", 1e-3*U0);
  printf("delta_rad  = %11.5e\n", delta_rad);
  printf("alpha_c    = %11.5e\n", alpha_c);
  printf("tau        = [%11.5e, %11.5e, %11.5e]\n", tau[X_], tau[Y_], tau[Z_]);
  printf("D          = [%11.5e, %11.5e, %11.5e]\n", D[X_], D[Y_], D[Z_]);

  cout << scientific << setprecision(6) << "\nCOD:\n" << setw(14) << ps_cod
       << "\n";

  PrtMap("\nA:", 3, A);
  printf("\nDet{A}-1 = %10.3e\n", DetMap(n_DOF, A)-1e0);

  PrtMap("\nM_lat:", 3, M_lat);
  printf("\nDet{M_lat}-1 = %10.3e\n", DetMap(n_DOF, M_lat)-1e0);
  PrtMap("\nM_cav:", 3, M_cav);

  PrtMap("\nM:", 3, M);
  printf("\nDet{M}-1 = %10.3e\n", DetMap(n_DOF, M)-1e0);

  // if (rad_on) {
  //   PrtMap("\nM_delta:", 3, M_delta);
  //   printf("\nDet{M_delta}-1 = %10.3e\n", DetMap(n_DOF, M_delta)-1e0);
  //   PrtMap("\nM_tau:", 3, M_tau);
  //   printf("\nDet{M_tau}-1 = %10.3e\n", DetMap(n_DOF, M_tau)-1e0);

  //   PrtMap("\nM_diff:", n_DOF, M_diff);
  //   PrtMap("\nM_Chol:", n_DOF, M_Chol);
  // }
}


void BeamType::BeamInit(const int n_part, const double eps_x,
			const double eps_y, const double sigma_s)
{
  int             j, k;
  double          phi;
  ss_vect<double> ps;

  this->n_part = n_part;

  globval.Cavity_on = false; globval.radiation = false;
  Ring_GetTwiss(true, 0e0);

  const double
    eps[]   = {eps_x, eps_y},
    alpha[] = {Cell[0].Alpha[X_], Cell[0].Alpha[Y_]},
    beta[]  = {Cell[0].Beta[X_],  Cell[0].Beta[Y_]};

  for (j = 0; j < n_part; j++) {
    ps.zero();
    for (k = 0; k < 2; k++) {
      phi = 2e0*M_PI*ranf();
      ps[2*k] = sqrt(beta[k]*eps[k])*cos(phi);
      ps[2*k+1] = -sqrt(eps[k]/beta[k])*(sin(phi)-alpha[k]*cos(phi));
      ps[ct_] = sigma_s*(2e0*ranf()-1e0);
      // ps_delta = ;
    }
    this->ps.push_back(ps);
  }
}


void BeamType::BeamStats(void)
{
  int             i, j, k;
  double          val;
  ss_vect<double> sum;
  ss_vect<tps>    sum2;

  sum.zero(); sum2.zero();
  for (i = 0; i < n_part; i++) {
    sum += ps[i];
    for (j = 0; j < 6; j++)
      for (k = 0; k < 6; k++)
	sum2[j] += ps[i][j]*ps[i][k]*tps(0e0, k+1);
  }

  sigma.zero();
  for (j = 0; j < 6; j++) {
    mean[j] = sum[j]/n_part;
    for (k = 0; k < 6; k++) {
      val = (n_part*sum2[j][k]-sum[j]*sum[k])/(n_part*(n_part-1e0));
      if (val >= 0e0) sigma[j] += sqrt(val)*tps(0e0, k+1);
    }
  }
}

void BeamType::print(void)
{
  int k;

  printf("\nmean: ");
  cout << scientific << setprecision(5) << setw(13) << mean << "\n";
  printf("sigma:");
  for (k = 0; k < 6; k++)
    printf(" %12.5e", sigma[k][k]);
  printf("\n");
}


void get_Poincare_Map(void)
{
  PoincareMap map;
  BeamType    beam;

  const int long seed = 1121;

  if (!false) no_sxt();

  if (false) {
    // map.GetM(false, false); map.print();
    // map.GetM(true, false); map.print();
    map.GetM(true, true); map.print();
    exit(0);
  }

  map.GetM(true, true); map.print();

  ss_vect<tps> M1, M2, M3;
  PrtMap("\nmap_num", 3, map.map_num);
  printf("\nDet{map_num}-1 = %10.3e\n", DetMap(3, map.map_num)-1e0);
  PrtMap("\nM_tau", 3, map.M_tau);
  M1 = map.M_tau*map.A*map.R*Inv(map.A);
  PrtMap("\nA*M_tau*R*Inv(A)", 3, M1);
  M2 = map.A*map.R*Inv(map.A);
  PrtMap("\nA*R*Inv(A)", 3, M2);
  M3 = map.M_cav*map.M_lat;
  PrtMap("\nM_cav*M_lat", 3, M3);
  exit(0);

  if (false) {
    PrtMap("\nM", 3, map.map_num);
    printf("\Det{M}-1 = %10.3e\n", DetMap(3, map.map_num)-1e0);
    PrtMap("A*M*A^-1", 3, Inv(map.A)*map.map_num*map.A);
    ss_vect<tps> A_A_tp = map.A*TpMap(3, map.A);
    PrtMap("A*A^tp", 3, A_A_tp);
    printf("\Det{A_A_tp}-1 = %10.3e\n", DetMap(3, A_A_tp)-1e0);
    A_A_tp = map.map_num*A_A_tp*TpMap(3, map.map_num);
    PrtMap("M*A*A_tp*M^tp", 3, A_A_tp);
    printf("\Det{A_A_tp}-1 = %10.3e\n", DetMap(3, A_A_tp)-1e0);
    exit(0);
  }

  if (!false) {
    Matrix       L_tp_mat;
    ss_vect<tps> L, L_tp, Id;
    Id.identity();
    L.zero();
    L[x_] = 3.6592e-06*Id[x_];
    L[px_] = 1.7514e-08*Id[x_] + 3.1063e-07*Id[px_];
    L[delta_] = 2.1727e-05*Id[x_] - 4.4389e-06*Id[px_] + 2.6758e-05*Id[delta_];
    L[ct_] =
      8.4551e-07*Id[x_] - 1.1983e-07*Id[px_]
      + 1.2481e-06*Id[delta_] + 8.7481e-07*Id[ct_];

    getlinmat(6, L, L_tp_mat); TpMat(6, L_tp_mat); L_tp = putlinmat(6, L_tp_mat);

    PrtMap("\nL:", 3, L);
    PrtMap("\nL*L_tp:", 3, L*L_tp);

    exit(0);
  }

  iniranf(seed); setrancut(5e0);

  printf("\ntracking:\n");
  beam.BeamInit(10000, 0e-9, 6e-9, 0e-3);
  beam.BeamStats(); beam.print();
  map.propagate(20000, beam);
  beam.BeamStats(); beam.print();

  double eps_x, sigma_delta, U_0, J[3];
  get_eps_x(eps_x, sigma_delta, U_0, J);
  printf("\nsigma_x, sigma_px: %12.5e %12.5e\n",
	 sqrt(Cell[0].Beta[X_]*eps_x
	      +sqr(map.M_lat[x_][delta_]*sigma_delta)),
	 sqrt((1e0+sqr(Cell[0].Alpha[X_]))/Cell[0].Beta[X_]*eps_x
	      +sqr(map.M_lat[px_][delta_]*sigma_delta)));
}
