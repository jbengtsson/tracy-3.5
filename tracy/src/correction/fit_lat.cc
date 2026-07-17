// Linear-lattice tune / chromaticity fits — see correction/fit_lat.h.


void corr::fit_tune(const std::string &fam_h, const std::string &fam_v,
		    const double nu_x, const double nu_y)
{
  double   TotalTuneX, TotalTuneY, dk;
  iVector2 nq;
  Vector2  nu;
  fitvect  qfbuf, qdbuf;
  long     i;

  printf("\ncorr::fit_tune: fitting nu.\n");
  dk = 1e-3;
  nq[0] = nq[1] = 0;
  nu[0] = nu_x;
  nu[1] = nu_y;
  for (i = 0; i <= globval.Cell_nLoc; i++) {
    if (Cell[i].Elem.Pkind == Mpole) {
      if (strncmp(Cell[i].Elem.PName, fam_h.c_str(), fam_h.length()) == 0) {
	qfbuf[nq[0]] = i;
	nq[0]++;
      }
      if (strncmp(Cell[i].Elem.PName, fam_v.c_str(), fam_v.length()) == 0) {
	qdbuf[nq[1]] = i;
	nq[1]++;
      }
    }
  }

  printf("Fittune: nq[0]=%ld nq[1]=%ld\n", nq[0], nq[1]);
  TotalTuneX = globval.TotalTune[0];
  TotalTuneY = globval.TotalTune[1];
  Ring_Fittune(nu, (double)1e-4, nq, qfbuf, qdbuf, dk, 50L);
  printf("Fittune: nux= %f dnux= %f nuy= %f dnuy= %f\n",
	 globval.TotalTune[0], globval.TotalTune[0] - TotalTuneX,
	 globval.TotalTune[1], globval.TotalTune[1] - TotalTuneY);

  Ring_GetTwiss(true, 0.0);
  printglob();
}


void corr::fit_chrom(const std::string &fam_h, const std::string &fam_v,
		     const double chrom_x, const double chrom_y)
{
  double   ChromaX, ChromaY, dks;
  iVector2 ns;
  Vector2  si;
  fitvect  sfbuf, sdbuf;
  long     i;

  printf("\ncorr::fit_chrom: fitting chi^(1)\n");
  dks = 1e-3;
  ns[0] = ns[1] = 0;
  si[0] = chrom_x;
  si[1] = chrom_y;
  for (i = 0; i <= globval.Cell_nLoc; i++) {
    if (Cell[i].Elem.Pkind == Mpole) {
      if (strncmp(Cell[i].Elem.PName, fam_h.c_str(), fam_h.length()) == 0) {
	sfbuf[ns[0]] = i;
	ns[0]++;
      }
      if (strncmp(Cell[i].Elem.PName, fam_v.c_str(), fam_v.length()) == 0) {
	sdbuf[ns[1]] = i;
	ns[1]++;
      }
    }
  }

  printf("Fitchrom: ns[0]=%ld ns[1]=%ld\n", ns[0], ns[1]);
  ChromaX = globval.Chrom[0];
  ChromaY = globval.Chrom[1];
  Ring_Fitchrom(si, 1e-4, ns, sfbuf, sdbuf, dks, 50L);
  printf("Fitchrom: six= %f dsix= %f siy= %f dsiy= %f\n",
	 globval.Chrom[0], globval.Chrom[0] - ChromaX, globval.Chrom[1],
	 globval.Chrom[1] - ChromaY);

  Ring_GetTwiss(true, 0.0);
  printglob();
}
