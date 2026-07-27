// Linear-lattice tune / chromaticity fits — see correction/fit_lat.h.

namespace corr {

// Exact, case-insensitive element-family lookup. ElemIndex does the same, but
// exit_()s on a miss with no context; returning 0 lets the caller name both the
// offending family and the keyword that set it. ElemF.PName is a partsName,
// i.e. a space-padded char[SymbolLength] that need not be NUL-terminated.
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
    printf("corr: no fit families — set %s\n", keyword);
    return false;
  }
  for (const auto &fam : fams) {
    const long k = fam_index(fam);

    if (k == 0) {
      printf("corr: fit family '%s' not found — set %s\n", fam.c_str(),
	     keyword);
      return false;
    }
    Fnum.push_back(k);
  }
  return true;
}


// Cell positions of a family's members, for the legacy fitters below.
// TRANSITIONAL: fitvect is a fixed long[fitvectmax] the fitters write into
// without a bound check, so a large family smashes the stack — hence the guard.
// Both go away with Ring_Fittune / Ring_Fitchrom in the next commit.
static bool fam_positions(const long Fnum, long buf[], long &n_mem)
{
  n_mem = GetnKid(Fnum);
  if (n_mem > fitvectmax) {
    printf("corr: fit family '%s' has %ld members, fitvect holds %d\n",
	   ElemFam[Fnum-1].ElemF.PName, n_mem, fitvectmax);
    return false;
  }
  for (long k = 1; k <= n_mem; k++)
    buf[k-1] = Elem_GetPos(Fnum, k);
  return true;
}


// TRANSITIONAL: Ring_Fittune / Ring_Fitchrom take exactly two families, so only
// the first two are used until the N-knob SVD solver lands.
static bool two_fams(const std::vector<long> &Fnum, iVector2 &n_mem,
		     long buf_h[], long buf_v[])
{
  if (Fnum.size() != 2) {
    printf("corr: %zu fit families given, the fitter takes exactly 2\n",
	   Fnum.size());
    return false;
  }
  return (fam_positions(Fnum[0], buf_h, n_mem[0])
	  && fam_positions(Fnum[1], buf_v, n_mem[1]));
}


bool fit_tune(const std::vector<std::string> &fams, const double nu_x,
	      const double nu_y)
{
  double            TotalTuneX, TotalTuneY;
  iVector2          nq;
  Vector2           nu;
  fitvect           qfbuf, qdbuf;
  std::vector<long> Fnum;

  const double dk = 1e-3;

  printf("\ncorr::fit_tune: fitting nu.\n");
  if (!resolve_fams("tune_fams", fams, Fnum)
      || !two_fams(Fnum, nq, qfbuf, qdbuf))
    return false;

  nu[0] = nu_x;
  nu[1] = nu_y;

  printf("Fittune: nq[0]=%ld nq[1]=%ld\n", nq[0], nq[1]);
  TotalTuneX = globval.TotalTune[0];
  TotalTuneY = globval.TotalTune[1];
  Ring_Fittune(nu, (double)1e-4, nq, qfbuf, qdbuf, dk, 50L);
  printf("Fittune: nux= %f dnux= %f nuy= %f dnuy= %f\n",
	 globval.TotalTune[0], globval.TotalTune[0] - TotalTuneX,
	 globval.TotalTune[1], globval.TotalTune[1] - TotalTuneY);

  Ring_GetTwiss(true, 0.0);
  printglob();

  return true;
}


bool fit_chrom(const std::vector<std::string> &fams, const double chrom_x,
	       const double chrom_y)
{
  double            ChromaX, ChromaY;
  iVector2          ns;
  Vector2           si;
  fitvect           sfbuf, sdbuf;
  std::vector<long> Fnum;

  const double dks = 1e-3;

  printf("\ncorr::fit_chrom: fitting chi^(1)\n");
  if (!resolve_fams("chrom_fams", fams, Fnum)
      || !two_fams(Fnum, ns, sfbuf, sdbuf))
    return false;

  si[0] = chrom_x;
  si[1] = chrom_y;

  printf("Fitchrom: ns[0]=%ld ns[1]=%ld\n", ns[0], ns[1]);
  ChromaX = globval.Chrom[0];
  ChromaY = globval.Chrom[1];
  Ring_Fitchrom(si, 1e-4, ns, sfbuf, sdbuf, dks, 50L);
  printf("Fitchrom: six= %f dsix= %f siy= %f dsiy= %f\n",
	 globval.Chrom[0], globval.Chrom[0] - ChromaX, globval.Chrom[1],
	 globval.Chrom[1] - ChromaY);

  Ring_GetTwiss(true, 0.0);
  printglob();

  return true;
}

}  // namespace corr
