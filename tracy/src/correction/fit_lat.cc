// Linear-lattice tune / chromaticity fits — see correction/fit_lat.h.

namespace corr {

// Drop knob directions whose singular value falls below this.
static const double fit_s_cut = 1e-10;


// Exact, case-insensitive family lookup; 0 on a miss, so the caller can name
// the offender. ElemF.PName is space-padded and need not be NUL-terminated.
static long fam_index(const std::string &name)
{
  std::string key = name;

  for (auto &c : key)
    c = tolower(c);
  for (long Fnum = 1; Fnum <= globval.Elem_nFam; Fnum++) {
    const char *PName = ElemFam[Fnum-1].ElemF.PName;
    size_t     len    = SymbolLength;

    while ((len > 0) && ((PName[len-1] == ' ') || (PName[len-1] == '\0')))
      len--;
    if ((key.length() == len) && (strncmp(key.c_str(), PName, len) == 0))
      return Fnum;
  }
  return 0;
}


// Resolve family names to indices; the error names the keyword that set them.
static bool resolve_fams(const char *keyword,
			 const std::vector<std::string> &fams,
			 std::vector<long> &Fnum)
{
  Fnum.clear();
  if (fams.empty()) {
    printf("\ncorr: no fit families, set %s\n", keyword);
    return false;
  }
  for (const auto &fam : fams) {
    const long k = fam_index(fam);

    if (k == 0) {
      printf("\ncorr: fit family '%s' not found, set %s\n", fam.c_str(),
	     keyword);
      return false;
    }
    Fnum.push_back(k);
  }
  return true;
}


enum fit_obs { fit_nu, fit_xi };


// Re-evaluate the observable. Returns false if the closed-orbit finder failed
// or the ring went unstable, which leaves val meaningless.
static bool get_obs(const fit_obs obs, double val[])
{
  if (obs == fit_nu)
    Ring_GetTwiss(false, 0e0);
  else
    Ring_Getchrom(0e0);

  if (!status.codflag || !globval.stable)
    return false;

  for (int k = 0; k < 2; k++)
    val[k] = (obs == fit_nu)? globval.TotalTune[k] : globval.Chrom[k];
  return true;
}


// Spread db_nL over the family's n_kid members — so db_nL is a whole-family
// step — and accumulate it in db_net, so a failure can undo what was applied.
static void apply_dbnL(const long Fnum, const int n, const double db_nL,
		       double &db_net)
{
  set_dbnL_design_fam(Fnum, n, db_nL/GetnKid(Fnum), 0e0);
  db_net += db_nL;
}


// Fit a 2-vector observable to target[] with the knob families in Fnum, whose
// order-n integrated strengths are the free parameters.
//
// The Jacobian A[j][k] = d(obs_j)/d(db_nL of family k) is built once by central
// differencing and re-used every iteration (chord iteration). It is 2 x n_knob
// and solved by SVD, so any knob count works.
//
// On instability or a lost closed orbit the knobs are rolled back to their entry
// values; on running out of steps the best iterate is kept. False either way.
static bool fit_lat(const char *what, const fit_obs obs,
		    const std::vector<long> &Fnum, const int n,
		    const double target[], const double db_nL,
		    const double eps, const int imax)
{
  bool     valid, converged = false;
  int      i, j, k;
  double   val[2] = {0e0, 0e0}, val_0[2] = {0e0, 0e0};
  double   val_best[2] = {0e0, 0e0}, res_best = 0e0;
  double   val_p[2], val_m[2], res, b_n, a_n;
  double   **A, **U, **V, *w, *dval, *db;

  const int m = 2, n_knob = Fnum.size();

  std::vector<double> db_net(n_knob, 0e0), db_best(n_knob, 0e0);

  A    = dmatrix(1, m, 1, n_knob);
  U    = dmatrix(1, m, 1, n_knob);
  V    = dmatrix(1, n_knob, 1, n_knob);
  w    = dvector(1, n_knob);
  dval = dvector(1, m);
  db   = dvector(1, n_knob);

  printf("\n%s: target [%9.5f, %9.5f], %d knob(s), db_%dL = %9.3e\n",
	 what, target[0], target[1], n_knob, n, db_nL);

  valid = get_obs(obs, val_0);

  // Jacobian by central differencing, once.
  for (k = 1; valid && (k <= n_knob); k++) {
    apply_dbnL(Fnum[k-1], n, db_nL, db_net[k-1]);
    if (!(valid = get_obs(obs, val_p))) break;
    apply_dbnL(Fnum[k-1], n, -2e0*db_nL, db_net[k-1]);
    if (!(valid = get_obs(obs, val_m))) break;
    apply_dbnL(Fnum[k-1], n, db_nL, db_net[k-1]);

    for (j = 1; j <= m; j++)
      A[j][k] = (val_p[j-1]-val_m[j-1])/(2e0*db_nL);
    if (trace)
      printf("  %-*.*s probe + [%9.5f, %9.5f], - [%9.5f, %9.5f]\n",
	     SymbolLength, SymbolLength, ElemFam[Fnum[k-1]-1].ElemF.PName,
	     val_p[0], val_p[1], val_m[0], val_m[1]);
  }

  if (valid) {
    corr::svd_decomp_cut(A, m, n_knob, U, w, V, fit_s_cut, trace);
    if (trace) dmdump(stdout, "\n  A:", A, m, n_knob, "%11.3e");
  }

  // Chord iteration: re-evaluate, solve A*db = target - current, apply.
  for (i = 0; valid; i++) {
    if (!(valid = get_obs(obs, val))) break;

    res = sqrt(sqr(target[0]-val[0])+sqr(target[1]-val[1]));
    if ((i == 0) || (res < res_best)) {
      res_best = res;
      db_best  = db_net;
      memcpy(val_best, val, sizeof(val_best));
    }
    if (trace)
      printf("  it %2d: [%9.5f, %9.5f], residual %9.3e\n", i, val[0], val[1],
	     res);
    if (res < eps) {
      converged = true;
      break;
    }
    if (i == imax) break;

    for (j = 1; j <= m; j++)
      dval[j] = target[j-1]-val[j-1];
    corr::svd_backsub(U, w, V, m, n_knob, dval, db);
    if (trace) dvdump(stdout, "\n  db_n:", db, n_knob, "%11.3e");

    for (k = 1; k <= n_knob; k++)
      apply_dbnL(Fnum[k-1], n, db[k], db_net[k-1]);
  }

  if (converged)
    printf("%s: [%9.5f, %9.5f] -> [%9.5f, %9.5f]\n", what, val_0[0], val_0[1],
	   val[0], val[1]);
  else if (!valid) {
    // Undo everything; adding -db_net zeroes db_net in passing.
    for (k = 0; k < n_knob; k++)
      apply_dbnL(Fnum[k], n, -db_net[k], db_net[k]);
    printf("%s: FAILED, unstable or no closed orbit; knobs fully restored,"
	   " back at [%9.5f, %9.5f]\n", what, val_0[0], val_0[1]);
  } else {
    for (k = 0; k < n_knob; k++)
      apply_dbnL(Fnum[k], n, db_best[k]-db_net[k], db_net[k]);
    printf("%s: not converged in %d step(s), keeping the best iterate:"
	   " [%9.5f, %9.5f] -> [%9.5f, %9.5f], residual %9.3e\n", what, imax,
	   val_0[0], val_0[1], val_best[0], val_best[1], res_best);
  }

  if (trace && valid) {
    printf("\n  b_%d of member 1 of each family:\n", n);
    for (k = 0; k < n_knob; k++) {
      get_bn_design_elem(Fnum[k], 1, n, b_n, a_n);
      printf("    %-*.*s %10.5f\n", SymbolLength, SymbolLength,
	     ElemFam[Fnum[k]-1].ElemF.PName, b_n);
    }
  }

  free_dmatrix(A, 1, m, 1, n_knob);
  free_dmatrix(U, 1, m, 1, n_knob);
  free_dmatrix(V, 1, n_knob, 1, n_knob);
  free_dvector(w, 1, n_knob);
  free_dvector(dval, 1, m);
  free_dvector(db, 1, n_knob);

  return converged;
}


bool fit_tune(const std::vector<std::string> &fams, const double nu_x,
	      const double nu_y, const double db_2L, const double eps,
	      const int imax)
{
  bool              ok;
  std::vector<long> Fnum;

  // Absolute, against the raw globval.TotalTune — not wrapped mod 1.
  const double target[] = {nu_x, nu_y};

  if (!resolve_fams("tune_fams", fams, Fnum))
    return false;

  ok = fit_lat("corr::fit_tune", fit_nu, Fnum, Quad, target, db_2L, eps, imax);

  Ring_GetTwiss(true, 0e0);
  printglob();

  return ok;
}


bool fit_chrom(const std::vector<std::string> &fams, const double chrom_x,
	       const double chrom_y, const double db_3L, const double eps,
	       const int imax)
{
  bool              ok, rad, cav;
  std::vector<long> Fnum;

  const double target[] = {chrom_x, chrom_y};

  if (!resolve_fams("chrom_fams", fams, Fnum))
    return false;

  // Ring_Getchrom is a linear-optics measurement, and get_DA_real calls this
  // right after GetEmittance turns radiation and the cavity on.
  rad = globval.radiation;
  cav = globval.Cavity_on;
  globval.radiation = false;
  globval.Cavity_on = false;

  ok = fit_lat("corr::fit_chrom", fit_xi, Fnum, Sext, target, db_3L, eps, imax);

  globval.radiation = rad;
  globval.Cavity_on = cav;

  Ring_GetTwiss(true, 0e0);
  printglob();

  return ok;
}

}  // namespace corr
