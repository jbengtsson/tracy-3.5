
struct PoincareMap {
  // Poincare map for tracking.

private:
public:
  double          nu[3], ksi1[2], ksi2[2], dnudJ[3];
  ss_vect<double> x_cod;
  ss_vect<tps>    M, A;

  ss_vect<tps> GetLat(void);
  ss_vect<tps> GetDelta(const double delta);
  ss_vect<tps> GetCav(void);
  ss_vect<tps> GetRad(void);
  ss_vect<tps> GetMap(const bool cav, const bool rad);
};


ss_vect<tps> PoincareMap::GetLat(void)
{
  int          k;
  double       mu[2], alpha[2], beta[2], gamma[2];
  ss_vect<tps> Id, M;

  const bool prt = !false;
  const double
    C       = Cell[globval.Cell_nLoc].S,
    alpha_c = globval.Alphac;

  Id.identity();

  for (k = 0; k < 2; k++) {
    mu[k] = 2e0*M_PI*globval.TotalTune[k];
    alpha[k] = Cell[0].Alpha[k];
    beta[k] = Cell[0].Beta[k];
    gamma[k] = (1e0+sqr(alpha[k]))/beta[k];
  }

  M.identity();
  for (k = 0; k < 2; k++) {
    M[2*k] =
      (cos(mu[k])+alpha[k]*sin(mu[k]))*Id[2*k] + beta[k]*sin(mu[k])*Id[2*k+1];
    M[2*k+1] =
      -gamma[k]*sin(mu[k])*Id[2*k]
      + (cos(mu[k])-alpha[k]*sin(mu[k]))*Id[2*k+1];
  }
  M[ct_] += alpha_c*C*Id[delta_];

  M[x_] += globval.OneTurnMat[x_][delta_]* Id[delta_];
  M[px_] += globval.OneTurnMat[px_][delta_]*Id[delta_];
  M[ct_] +=
    (M[x_][x_]*M[px_][delta_]-M[px_][x_]*M[x_][delta_])*Id[x_]
    +(M[x_][px_]*M[px_][delta_]-M[px_][px_]*M[x_][delta_])*Id[px_];

  if (prt) {
    printf("\nGetLat:\n  C, alpha_c: %7.3f %12.5e\n", C, alpha_c);
    prt_lin_map(3, M);
  }

  return M;
}


ss_vect<tps> PoincareMap::GetDelta(const double delta)
{
  ss_vect<tps> Id, M;

  const bool prt = !false;

  Id.identity();

  M.identity();
  M[ct_] /= 1e0 + delta;
  M[delta_] *= 1e0 + delta;

  if (prt) {
    printf("\nGetDelta:\n  delta %12.5e\n", delta);
    prt_lin_map(3, M);
  }

  return M;
}


ss_vect<tps> PoincareMap::GetCav(void)
{
  ss_vect<tps> Id, M;

  const bool prt = !false;
  const int  loc = Elem_GetPos(ElemIndex("cav"), 1);
  const double
    V_RF = Cell[loc].Elem.C->Pvolt,
    f_RF = Cell[loc].Elem.C->Pfreq,
    phi0 = asin(globval.U0/V_RF);

  Id.identity();

  M.identity();
  M[delta_] -= V_RF*2e0*M_PI*f_RF*cos(phi0)/(1e9*globval.Energy*c0)*Id[ct_];

  if (prt) {
    printf("\nGetCav:\n  V_RF, U0, phi0: %10.3e %10.3e %9.6f\n",
	   V_RF, globval.U0, phi0*180.0/M_PI);
    prt_lin_map(3, M);
  }

  return M;
}


ss_vect<tps> PoincareMap::GetRad(void)
{
  int          k;
  double       tau[3];
  ss_vect<tps> Id, M;

  const bool   prt = !false;
  const double C   = Cell[globval.Cell_nLoc].S;

  Id.identity();

  M.zero();
  for (k = 0; k < 3; k++) {
    tau[k] = -C/(c0*globval.alpha_rad[k]);
    M[2*k] = exp(-C/(c0*tau[k]))*Id[2*k];
    M[2*k+1] = exp(-C/(c0*tau[k]))*Id[2*k+1];
  }

  if (prt) {
    printf("\nGetRad:\n  tau = [%11.5e, %11.5e, %11.5e]\n",
	   tau[X_], tau[Y_], tau[Z_]);
    prt_lin_map(3, M);
  }

  return M;
}


ss_vect<tps> PoincareMap::GetMap(const bool cav, const bool rad)
{
  double       delta_rad, C, tau_z, delta_rad2;
  ss_vect<tps> M;

  globval.U0 = globval.dE*1e9*globval.Energy;
  delta_rad = globval.U0/(1e9*globval.Energy);

  C     = Cell[globval.Cell_nLoc].S;
  tau_z = -C/(c0*globval.alpha_rad[Z_]);
  delta_rad2 = exp(-C/(c0*tau_z)) - 1e0;

  // printf("\n  delta_rad delta_rad2 delta_cod %12.5e %12.5e %12.5e\n",
  // 	 delta_rad, delta_rad2, globval.CODvect[delta_]);

  M = GetLat();
  
  if (rad) M = GetDelta(delta_rad)*M;
  if (cav) M = GetCav()*M;
  if (rad) M = GetRad()*M;

  return M;
}
