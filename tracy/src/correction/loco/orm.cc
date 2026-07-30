// Measured orbit response matrix. See correction/loco/orm.h.


void corr::measure_orm_column(const int fnum, const int knum,
			      const int kick_type, const double kick,
			      const int bpm_loc[], const int n_bpm,
			      const int read_coord, double resp[])
{
  bool     cod;
  int      i;
  long int lastpos;
  double   *orbitP, *orbitN;

  orbitP = dvector(1, n_bpm); orbitN = dvector(1, n_bpm);

  SetdKLpar(fnum, knum, kick_type, kick);
  cod = getcod(0e0, lastpos); chk_cod(cod, "measure_orm_column");
  for (i = 1; i <= n_bpm; i++)
    orbitP[i] = Cell[bpm_loc[i-1]].BeamPos[read_coord];

  SetdKLpar(fnum, knum, kick_type, -2e0*kick);
  cod = getcod(0e0, lastpos); chk_cod(cod, "measure_orm_column");
  for (i = 1; i <= n_bpm; i++)
    orbitN[i] = Cell[bpm_loc[i-1]].BeamPos[read_coord];

  SetdKLpar(fnum, knum, kick_type, kick);

  for (i = 1; i <= n_bpm; i++)
    resp[i] = (orbitP[i]-orbitN[i])*0.5/kick;

  free_dvector(orbitP, 1, n_bpm); free_dvector(orbitN, 1, n_bpm);
}
