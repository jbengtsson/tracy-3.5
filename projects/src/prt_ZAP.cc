void prt_ZAP(const int n)
{
  long int k;
  FILE     *outf;

  outf = file_write("ZAPLAT.DAT");

  fprintf(outf, "%ld %7.5f\n",
	  globval.Cell_nLoc+1, n*Cell[globval.Cell_nLoc].S);
  fprintf(outf, "One super period\n");

  for (k = 0; k <= globval.Cell_nLoc; k++)
    fprintf(outf, "%10.5f %8.5f %9.6f %8.5f %7.3f %8.5f %8.5f %7.5f\n",
	    Cell[k].S,
	    Cell[k].Beta[X_], Cell[k].Alpha[X_],
	    Cell[k].Beta[Y_], Cell[k].Alpha[Y_],
	    Cell[k].Eta[X_], Cell[k].Etap[X_], Cell[k].maxampl[X_][1]);

  fprintf(outf, "0\n");

  fclose(outf);
}
