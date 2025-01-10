void prt_lat_maxlab(const char *fname, const int Fnum, const bool all)
{
  // Generate CSV file for the linear optics.
  long int i = 0;
  FILE     *outf;

  outf = fopen(fname, "w");
  fprintf(outf, "#        name               s      type"
	  "    alphax    betax      nux       etax     etapx");
  fprintf(outf, "      alphay    betay      nuy      etay      etapy\n");
  fprintf(outf, "#                          [m]"
	  "                        [m]                 [m]");
  fprintf(outf, "                            [m]                [m]\n");
  fprintf(outf, "#\n");

  for (i = 0; i <= globval.Cell_nLoc; i++) {
    if (all || (Cell[i].Fnum == Fnum)) {
      fprintf(outf, "%4ld, ", i);
      prt_name(outf, Cell[i].Elem.PName);
      fprintf(outf, " %9.5f, %4.1f,"
           " %9.5f, %8.5f, %8.5f, %8.5f, %8.5f,"
           " %9.5f, %8.5f, %8.5f, %8.5f, %8.5f\n",
           Cell[i].S, get_code(Cell[i]),
           Cell[i].Alpha[X_], Cell[i].Beta[X_], Cell[i].Nu[X_],
           Cell[i].Eta[X_], Cell[i].Etap[X_],
           Cell[i].Alpha[Y_], Cell[i].Beta[Y_], Cell[i].Nu[Y_],
           Cell[i].Eta[Y_], Cell[i].Etap[Y_]);
    }
  }

  fclose(outf);
}
