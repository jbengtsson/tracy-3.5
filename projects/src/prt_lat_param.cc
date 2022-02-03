void get_bend(const char *name, double &L, double &B, double &B_2)
{
  int    loc;
  double rho_inv, b_2;

  const double Brho = globval.Energy*1e9/c0;

  loc = Elem_GetPos(ElemIndex(name), 1);
  L = Cell[loc].Elem.PL;
  rho_inv = Cell[loc].Elem.M->Pirho;
  b_2 = Cell[loc].Elem.M->PBpar[Quad+HOMmax];
  B = Brho*rho_inv;
  B_2 = Brho*b_2;
}


void get_RB(const char *name, int &n_RB, double &B_2, double &dx)
{
  double L, B;

  get_bend(name, L, B, B_2);
  n_RB = GetnKid(ElemIndex(name));
  dx = B/B_2;
}


void get_max_bn(const int n, int &n_bn, double &bn_max)
{
  int k;

  n_bn = 0; bn_max = 0e0;
  for (k = 0; k <= globval.Cell_nLoc; k++)
    if (Cell[k].Elem.Pkind == Mpole)
      if (Cell[k].Elem.M->PBpar[n+HOMmax] != 0e0) {
	n_bn++;
	bn_max = max(fabs(Cell[k].Elem.M->PBpar[n+HOMmax]), bn_max);
      }
}


void set_sxt(void)
{
  int    k;
  double b_n;

  for (k = 0; k < globval.Elem_nFam; k++)
    if ((ElemFam[k].ElemF.Pkind == Mpole)
	&& (ElemFam[k].ElemF.M->n_design == Sext)) {
      b_n = ElemFam[k].ElemF.M->PBpar[Sext+HOMmax];
      set_bn_design_fam(k+1, Sext, b_n, 0e0);
    }
}


void prt_lat_param(char *file_name, char *cav_name, const double eps_s,
		   const double beta_s)
{
  int    loc, n_bn, n_RB, n_RB_tot;
  double bn_max, B_2, dx, dx_max, L, B, eps_x, sigma_delta, U_0, J[3], tau[3];
  double I[6], curly_H_x, eps_eff_x;
  FILE   *outf;

  const double Brho = globval.Energy*1e9/c0;

  outf = file_write(file_name);

  fprintf(outf, "\n\nBasic parameters\n");

  fprintf(outf, "C [m]\t%5.3f\n", Cell[globval.Cell_nLoc].S);

  loc = Elem_GetPos(ElemIndex(cav_name), 1);
  fprintf(outf, "Harmonic no\t%1d\n", Cell[loc].Elem.C->Ph);
  fprintf(outf, "RF frequency [MHz]\t%8.6f\n", 1e-6*Cell[loc].Elem.C->Pfreq);

  fprintf(outf, "LS [m]\t%4.2f\n",
	 2e0*Cell[Elem_GetPos(ElemIndex("quad_add"), 1)-1].S);
  fprintf(outf, "SS [m]\t%4.2f\n",
	 Cell[Elem_GetPos(ElemIndex("qf1"), 2)-1].S
	 -Cell[Elem_GetPos(ElemIndex("qf1"), 1)].S);
  fprintf(outf, "MS [m]\t%4.2f\n",
	 Cell[Elem_GetPos(ElemIndex("sh2"), 2)-1].S
	 -Cell[Elem_GetPos(ElemIndex("sh2"), 1)].S);
  fprintf(outf, "\n\n\n");

  get_max_bn(Quad, n_bn, bn_max);
  fprintf(outf, "B_2^ [T/m]\t%1.0f\n", Brho*bn_max);
  get_max_bn(Sext, n_bn, bn_max);
  fprintf(outf, "b_3^ [T/mˆ2]\t%1.0f\n", Brho*bn_max);
  get_max_bn(Oct, n_bn, bn_max);
  fprintf(outf, "b_4^ [T/mˆ3]\t%1.0f\n", Brho*bn_max);
  fprintf(outf, "n_b_4\t%1d\n", n_bn);

  get_RB("qf4", n_RB_tot, B_2, dx);
  dx_max = max(fabs(dx), dx_max);
  get_RB("qf8", n_RB, B_2, dx);
  n_RB_tot += n_RB;
  dx_max = max(fabs(dx), dx_max);
  fprintf(outf, "dx^ [mm]\t%3.1f\n", 1e3*dx_max);
  fprintf(outf, "No of RB\t%1d\n", n_RB_tot);
  get_RB("qf4", n_RB_tot, B_2, dx);
  fprintf(outf, "QF4 RB [T/m], [mm]\t%3.1f, %4.2f\n", B_2, 1e3*dx);
  fprintf(outf, "QF4 RB [T/m], [mm]\t%3.1f, %4.2f\n", B_2, 1e3*dx);
  fprintf(outf, "QF4 RB [T/m], [mm]\t%3.1f, %4.2f\n", B_2, 1e3*dx);
  get_RB("qf6", n_RB_tot, B_2, dx);
  fprintf(outf, "QF6 RB [T/m], [mm]\t%3.1f, %4.2f\n", B_2, 1e3*dx);
  get_RB("qf8", n_RB_tot, B_2, dx);
  fprintf(outf, "QF8 RB [T/m], [mm]\t%3.1f, %4.2f\n", B_2, 1e3*dx);

  get_bend("dl1a_1", L, B, B_2);
  fprintf(outf, "\nDL1A_1 [m], [T], [T/m]\t%5.3f, %5.3f, %4.2f\n", L, B, B_2);
  get_bend("dl1a_2", L, B, B_2);
  fprintf(outf, "DL1A_2 [m], [T], [T/m]\t%5.3f, %5.3f, %4.2f\n", L, B, B_2);
  get_bend("dl1a_3", L, B, B_2);
  fprintf(outf, "DL1A_3 [m], [T], [T/m]\t%5.3f, %5.3f, %4.2f\n", L, B, B_2);
  get_bend("dl1a_4", L, B, B_2);
  fprintf(outf, "DL1A_4 [m], [T], [T/m]\t%5.3f, %5.3f, %4.2f\n", L, B, B_2);
  get_bend("dl1a_5", L, B, B_2);
  fprintf(outf, "DL1A_5 [m], [T], [T/m]\t%5.3f, %5.3f, %4.2f\n", L, B, B_2);

  get_bend("dl2a_1", L, B, B_2);
  fprintf(outf, "DL2A_1 [m], [T], [T/m]\t%5.3f, %5.3f, %4.2f\n", L, B, B_2);
  get_bend("dl2a_2", L, B, B_2);
  fprintf(outf, "DL2A_2 [m], [T], [T/m]\t%5.3f, %5.3f, %4.2f\n", L, B, B_2);
  get_bend("dl2a_3", L, B, B_2);
  fprintf(outf, "DL2A_3 [m], [T], [T/m]\t%5.3f, %5.3f, %4.2f\n", L, B, B_2);
  get_bend("dl2a_4", L, B, B_2);
  fprintf(outf, "DL2A_4 [m], [T], [T/m]\t%5.3f, %5.3f, %4.2f\n", L, B, B_2);
  get_bend("dl2a_5", L, B, B_2);
  fprintf(outf, "DL2A_5 [m], [T], [T/m]\t%5.3f, %5.3f, %4.2f\n", L, B, B_2);

  get_bend("dq1", L, B, B_2);
  fprintf(outf, "DQ [m], [T], [T/m]\t%5.3f, %5.3f, %4.2f\n", L, B, B_2);

  no_sxt();
  Ring_GetTwiss(true, 0e0);  printglob();
  fprintf(outf, "\nnu\t%4.2f, %4.2f\n",
	 globval.TotalTune[X_], globval.TotalTune[Y_]);
  fprintf(outf, "ksi\t%1.0f, %1.0f\n",
	  globval.Chrom[X_], globval.Chrom[Y_]);
  set_sxt();
  Ring_GetTwiss(true, 0e0);  printglob();
  fprintf(outf, "ksi corr.\t%4.2f, %4.2f\n",
	  globval.Chrom[X_], globval.Chrom[Y_]);

  fprintf(outf, "alpha_c (w/o IDs) [1e-4]\t%4.2f\n", 1e4*globval.Alphac);
  fprintf(outf, "alpha_c (w/ IDs) [1e-4]\t%4.2f\n", NAN);

  globval.Cavity_on = true; globval.radiation = true;
  Ring_GetTwiss(true, 0e0);
  fprintf(outf, "nu_s (w/o IDs)\t%7.5f\n", fabs(globval.Omega));
  fprintf(outf, "nu_s (w/ IDs)\t%7.5f\n", NAN);
  globval.Cavity_on = false; globval.radiation = false;
  Ring_GetTwiss(true, 0e0);

  fprintf(outf, "\nEmittance & energy spread\n");

  get_eps_x(eps_x, sigma_delta, U_0, J, tau, I, false);
  fprintf(outf, "eps_x (w/o IDs) [pm.rad]\t%1.0f\n", 1e12*eps_x);
  loc = Elem_GetPos(ElemIndex("ms"), 1);
  curly_H_x = sqr(Cell[loc].Eta[X_])/Cell[loc].Beta[X_];
  eps_eff_x = eps_x*sqrt(1e0+curly_H_x*sqr(sigma_delta)/eps_x);
  fprintf(outf, "epsˆ*_x MS [pm.rad]\t%1.0f\n", 1e12*eps_eff_x);
  fprintf(outf, "\n\n");
  fprintf(outf, "eps_x (w/ IDs) [pm.rad]\t%1.0f\n", NAN);
  fprintf(outf, "\n\n");
  fprintf(outf, "sigma_delta (w/o IDs) [1e-3]\t%4.2f\n", 1e3*sigma_delta);
  fprintf(outf, "\n\n");
  fprintf(outf, "sigma_delta (w/ IDs) [1e-3]\t%4.2f\n", NAN);
  fprintf(outf, "\n\n");
  fprintf(outf, "sigma_s (w/o IDs) [mm]\t%4.2f\n", 1e3*sqrt(beta_s*eps_s));
  fprintf(outf, "\n\n\n\n\n");

  fprintf(outf, "\nU_0 (w/o IDs) [MeV]\t%5.3f\n", 1e-6*U_0);
  fprintf(outf, "U_0 (w/ IDs) [MeV]\t%4.2f\n", NAN);
  fprintf(outf, "J_x (w/o IDs) [MeV]\t%4.2f\n", J[X_]);
  fprintf(outf, "J_x (w/ IDs) [MeV]\t%4.2f\n", NAN);
  fprintf(outf, "J_s (w/o IDs) [MeV]\t%4.2f\n", J[Z_]);
  fprintf(outf, "J_s (w/ IDs) [MeV]\t%4.2f\n", NAN);

  fprintf(outf, "I1\t%5.3f\n", I[1]);
  fprintf(outf, "I2\t%5.3f\n", I[2]);
  fprintf(outf, "I3\t%5.3f\n", I[3]);
  fprintf(outf, "I4\t%5.3f\n", I[4]);
  fprintf(outf, "I5\t%10.8f\n", I[5]);

  loc = Elem_GetPos(ElemIndex("ls"), 1);
  fprintf(outf, "beta_x, beta_y, eta_x (LS)\t%3.1f, %3.1f, %3.1f\n",
	 Cell[loc].Beta[X_], Cell[loc].Beta[Y_], Cell[loc].Eta[X_]);
  loc = Elem_GetPos(ElemIndex("ss"), 1);
  fprintf(outf, "beta_x, beta_y, eta_x (SS)\t%3.1f, %3.1f, %3.1f\n",
	 Cell[loc].Beta[X_], Cell[loc].Beta[Y_], Cell[loc].Eta[X_]);
  loc = Elem_GetPos(ElemIndex("ms"), 1);
  fprintf(outf, "beta_x, beta_y, eta_x (MS)\t%3.1f, %3.1f, %5.3f\n",
	 Cell[loc].Beta[X_], Cell[loc].Beta[Y_], Cell[loc].Eta[X_]);

  fprintf(outf, "tau (w/o IDs) [msec]\t%3.1f, %3.1f, %3.1f\n",
	 1e3*tau[X_], 1e3*tau[Y_], 1e3*tau[Z_]);
  fprintf(outf, "tau (w/ IDs)\t%3.1f, %3.1f, %3.1f\n", NAN, NAN, NAN);

  fprintf(outf, "\nLifetime\n");

  fclose(outf);
}
