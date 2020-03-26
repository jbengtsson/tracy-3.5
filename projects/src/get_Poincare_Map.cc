

ss_vect<tps> get_Poincare_Map(const bool cav, const bool rad)
{
  int          k;
  double       C, alpha_c, mu[3], alpha[3], beta[3], gamma[3], tau[3];
  ss_vect<tps> Id, Id_can, M, M_rad, A, A_A_tp;

  const int n_DOF = (!cav)? 2 : 3;

  Id.identity();

  C = Cell[globval.Cell_nLoc].S;
  alpha_c = globval.Alphac;

  for (k = 0; k < 2; k++) {
    mu[k] = 2e0*M_PI*globval.TotalTune[k];
    alpha[k] = Cell[0].Alpha[k]; beta[k] = Cell[0].Beta[k];
    gamma[k] = (1e0+sqr(alpha[k]))/beta[k];
  }

  if (cav) {
    mu[Z_] = -2e0*M_PI*globval.Omega;

    if (!false) {
      putlinmat(2*n_DOF, globval.Ascr, A);
      // Canonical coordinates are [ct, -delta].
      Id_can.identity(); Id_can[delta_] *= -1e0;
      A = A*Id_can;
      A_A_tp = A*tp_S(n_DOF, A);
      prt_lin_map(3, A);
      prt_lin_map(3, A_A_tp);

      putlinmat(2*n_DOF, globval.OneTurnMat, M);
      prt_lin_map(3, Inv(A)*M*A);

      alpha[Z_] = A_A_tp[ct_][delta_];
      beta[Z_] = -A_A_tp[ct_][ct_];
      gamma[Z_] = A_A_tp[delta_][delta_];

      printf("\n  nu_z, alpha_z, beta_z: %12.5e %12.5e %12.5e\n",
	     -globval.Omega, alpha[Z_], beta[Z_]);
    }

    alpha[Z_] =
      -globval.Ascr[ct_][ct_]*globval.Ascr[delta_][ct_]
      - globval.Ascr[ct_][delta_]*globval.Ascr[delta_][delta_];
    beta[Z_] = sqr(globval.Ascr[ct_][ct_]) + sqr(globval.Ascr[ct_][delta_]);
    gamma[Z_] = (1e0+sqr(alpha[Z_]))/beta[Z_];
  }

  M.identity();
  for (k = 0; k < 3; k++) {
    M[2*k] =
      (cos(mu[k])+alpha[k]*sin(mu[k]))*Id[2*k] + beta[k]*sin(mu[k])*Id[2*k+1];
    M[2*k+1] =
      -gamma[k]*sin(mu[k])*Id[2*k]
      + (cos(mu[k])-alpha[k]*sin(mu[k]))*Id[2*k+1];

    if (k == 2) {
      // Reverse order for [ct, delta].
      if (!cav)
	M[ct_] += alpha_c*C*Id[delta_];
      else {
#if 1
	// Exact.
	M[2*k+1] =
	  (cos(mu[k])+alpha[k]*sin(mu[k]))*Id[2*k+1]
	  + beta[k]*sin(mu[k])*Id[2*k];
	M[2*k] =
	  -gamma[k]*sin(mu[k])*Id[2*k+1]
	  + (cos(mu[k])-alpha[k]*sin(mu[k]))*Id[2*k];
#else
	// Leading order.
	M[2*k+1] = Id[2*k+1] + globval.OneTurnMat[ct_][delta_]*Id[2*k];
	M[2*k] =
	  -sqr(mu[Z_])/globval.OneTurnMat[ct_][delta_]*Id[2*k+1]
	  + (1e0-sqr(mu[Z_]))*Id[2*k];
#endif
      }
    }
  }

  M[x_] += globval.OneTurnMat[x_][delta_]* Id[delta_];
  M[px_] += globval.OneTurnMat[px_][delta_]*Id[delta_];
  M[ct_] +=
    (M[x_][x_]*M[px_][delta_]-M[px_][x_]*M[x_][delta_])*Id[x_]
    +(M[x_][px_]*M[px_][delta_]-M[px_][px_]*M[x_][delta_])*Id[px_];
  if (cav)
    M[delta_] +=
      -sqr(mu[Z_])*M[ct_][x_]/M[ct_][delta_]*Id[x_]
      -sqr(mu[Z_])*M[ct_][px_]/M[ct_][delta_]*Id[px_];

  if (rad) {
    M_rad.zero();
    for (k = 0; k < 3; k++) {
      tau[k] = -C/(c0*globval.alpha_rad[k]);
      M_rad[2*k] = exp(-C/(c0*tau[k]))*Id[2*k];
      M_rad[2*k+1] = exp(-C/(c0*tau[k]))*Id[2*k+1];
    }
    M = M_rad*M;
  }

  printf("\n  nu_z, alpha_z, beta_z: %12.5e %12.5e %12.5e\n",
	 -globval.Omega, alpha[Z_], beta[Z_]);
  printf("  tau =                   [%11.5e, %11.5e, %11.5e]\n",
	 tau[X_], tau[Y_], tau[Z_]);

  return M;
}
