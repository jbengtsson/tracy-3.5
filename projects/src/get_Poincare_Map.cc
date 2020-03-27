ss_vect<tps> get_M_lat(void)
{
  int          k;
  double       alpha_c, mu[2], alpha[2], beta[2], gamma[2];
  ss_vect<tps> Id, M;

  const bool   prt = false;
  const double C = Cell[globval.Cell_nLoc].S;

  Id.identity();

  alpha_c = globval.Alphac;

  for (k = 0; k < 2; k++) {
    mu[k] = 2e0*M_PI*globval.TotalTune[k];
    alpha[k] = Cell[0].Alpha[k]; beta[k] = Cell[0].Beta[k];
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
  // Reverse order for [ct, delta].
  M[ct_] += alpha_c*C*Id[delta_];

  M[x_] += globval.OneTurnMat[x_][delta_]* Id[delta_];
  M[px_] += globval.OneTurnMat[px_][delta_]*Id[delta_];
  M[ct_] +=
    (M[x_][x_]*M[px_][delta_]-M[px_][x_]*M[x_][delta_])*Id[x_]
    +(M[x_][px_]*M[px_][delta_]-M[px_][px_]*M[x_][delta_])*Id[px_];

  if (prt) {
    printf("\nget_M_lat:\n  C, alpha_c: %7.3f %12.5e\n", C, alpha_c);
    prt_lin_map(3, M);
  }

  return M;
}


ss_vect<tps> get_M_cav(const double delta)
{
  int          loc;
  double       V_RF, f_RF, phi0;
  ss_vect<tps> Id, M;

  const bool prt = !false;

  Id.identity();

  loc = Elem_GetPos(ElemIndex("cav"), 1);
  V_RF = Cell[loc].Elem.C->Pvolt;
  f_RF = Cell[loc].Elem.C->Pfreq;
  globval.U0 = globval.dE*1e9*globval.Energy;
  phi0 = globval.CODvect[ct_]*f_RF/c0*2e0*M_PI;

  M.identity();
  M[ct_] *= 1e0 + 2e0*delta;
  M[delta_] /= 1e0 + 2e0*delta;
  M[delta_] += -V_RF*(2e0*M_PI*f_RF/c0)/(1e9*globval.Energy)*Id[ct_];

  if (prt) {
    printf("\nget_M_cav:\n  V_RF, U0, phi0, delta: %10.3e %10.3e %10.3e\n",
	   V_RF, globval.U0, phi0*180.0/M_PI);
    cout << scientific << setprecision(3)
	 << "  cod:" << setw(11) << globval.CODvect << "\n";
    prt_lin_map(3, M);
  }

  return M;
}


ss_vect<tps> get_M_rad(void)
{
  int          k;
  double       tau[3];
  ss_vect<tps> Id, M;

  const bool   prt = false;
  const double C = Cell[globval.Cell_nLoc].S;

  Id.identity();

  M.zero();
  for (k = 0; k < 3; k++) {
    tau[k] = -C/(c0*globval.alpha_rad[k]);
    M[2*k] = exp(-C/(c0*tau[k]))*Id[2*k];
    M[2*k+1] = exp(-C/(c0*tau[k]))*Id[2*k+1];
  }

  if (prt) {
    printf("\nget_M_rad:\n  tau = [%11.5e, %11.5e, %11.5e]\n",
	   tau[X_], tau[Y_], tau[Z_]);
    prt_lin_map(3, M);
  }

  return M;
}


ss_vect<tps> get_Poincare_Map(const bool cav, const bool rad)
{
  ss_vect<tps> M;

  putlinmat(6, globval.OneTurnMat, M);
  prt_lin_map(3, Inv(get_M_rad())*M*Inv(get_M_lat()));

  M = get_M_lat();
  if (cav) M = get_M_cav(globval.CODvect[delta_])*M;
  // if (rad) M = get_M_rad()*M;

  return M;
}
