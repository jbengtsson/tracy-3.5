/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        


   To generate a lattice flat file.

   Type codes:

     marker     -1
     drift	 0
     multipole   1
     cavity      2
     thin kick   3
     wiggler     4
     kick_map    6

   Integration methods:

     fourth order symplectic integrator   4

   Format:

     name, family no, kid no, element no
     type code, integration method, no of integration steps
     apertures: xmin, xmax, ymin, ymax

   The following lines follows depending on element type.

     type

     drift:	 L

     multipole:  hor., ver. displacement, roll angle (design),
                                          roll angle (error)
                 L, 1/rho, entrance angle, exit angle
                 apertures[4]
		 no of nonzero multipole coeff., n design
		 n, b , a
		     n   n
		    ...

     wiggler:    L [m], lambda [m]
                 no of harmonics
                 harm no, kxV [1/m], BoBrhoV [1/m], kxH, BoBrhoH, phi
                    ...

     cavity:	 cavity voltage/beam energy [eV], omega/c, beam energy [eV]

     thin kick:	 hor., ver. displacement, roll angle (total)
		 no of nonzero multipole coeff.
		 n, b , a
		     n   n
		    ...

     kick_map:   order <file name>

*/


#define snamelen 10

// numerical type codes
#define marker_   -1
#define drift_     0
#define mpole_     1
#define cavity_    2
#define thinkick_  3
#define wiggler_   4
#define kick_map_  6


void prtName(FILE *fp, const int i,
	     const int type, const int method, const int N)
{
  fprintf(fp, "%-15s %4d %4d %4d\n",
	  Cell[i].Elem.PName, Cell[i].Fnum, Cell[i].Knum, i);
  fprintf(fp, " %3d %3d %3d\n", type, method, N);
  fprintf(fp, " %23.16e %23.16e %23.16e %23.16e\n",
	  Cell[i].maxampl[X_][0], Cell[i].maxampl[X_][1],
	  Cell[i].maxampl[Y_][0], Cell[i].maxampl[Y_][1]);
}


void prtHOM(FILE *fp, const int n_design, const mpolArray PB, const int Order)
{
  int i, nmpole;
  
  nmpole = 0;
  for (i = 1; i <= Order; i++)
    if ((PB[HOMmax-i] != 0.0) || (PB[HOMmax+i] != 0.0)) nmpole++;
  fprintf(fp, "  %2d %2d\n", nmpole, n_design);
  for (i = 1; i <= Order; i++) {
    if ((PB[HOMmax-i] != 0.0) || (PB[HOMmax+i] != 0.0))
      fprintf(fp, "%3d %23.16e %23.16e\n", i, PB[HOMmax+i], PB[HOMmax-i]);
  }
}


void prtmfile(const char mfile_dat[])
{
  int     i, j;
  FILE    *mfile;

  mfile = file_write(mfile_dat);
  for (i = 0; i <= globval.Cell_nLoc; i++) {
    switch (Cell[i].Elem.Pkind) {
    case drift:
      prtName(mfile, i, drift_, 0, 0);
      fprintf(mfile, " %23.16e\n", Cell[i].Elem.PL);
      break;
    case Mpole:
      if (Cell[i].Elem.PL != 0.0) {
	prtName(mfile, i, mpole_, Cell[i].Elem.M->Pmethod, Cell[i].Elem.M->PN);
	fprintf(mfile, " %23.16e %23.16e %23.16e %23.16e\n",
		Cell[i].dS[X_], Cell[i].dS[Y_],
		Cell[i].Elem.M->PdTpar,
		Cell[i].Elem.M->PdTsys
		+Cell[i].Elem.M->PdTrms*Cell[i].Elem.M->PdTrnd);
	fprintf(mfile, " %23.16e %23.16e %23.16e %23.16e %23.16e\n",
		Cell[i].Elem.PL, Cell[i].Elem.M->Pirho,
		Cell[i].Elem.M->PTx1, Cell[i].Elem.M->PTx2,
		Cell[i].Elem.M->Pgap);
	prtHOM(mfile, Cell[i].Elem.M->n_design, Cell[i].Elem.M->PB,
	       Cell[i].Elem.M->Porder);
      } else {
	prtName(mfile, i, thinkick_, Cell[i].Elem.M->Pmethod,
		Cell[i].Elem.M->PN);
	fprintf(mfile, " %23.16e %23.16e %23.16e\n",
		Cell[i].dS[X_], Cell[i].dS[Y_],
		Cell[i].Elem.M->PdTsys
		+Cell[i].Elem.M->PdTrms*Cell[i].Elem.M->PdTrnd);
	prtHOM(mfile, Cell[i].Elem.M->n_design, Cell[i].Elem.M->PB,
	       Cell[i].Elem.M->Porder);
      }
      break;
    case Wigl:
      prtName(mfile, i, wiggler_, Cell[i].Elem.W->Pmethod, Cell[i].Elem.W->PN);
      fprintf(mfile, " %23.16e %23.16e\n",
	      Cell[i].Elem.PL, Cell[i].Elem.W->lambda);
      fprintf(mfile, "%2d\n", Cell[i].Elem.W->n_harm);
      for (j = 0; j < Cell[i].Elem.W->n_harm; j++) {
	fprintf(mfile, "%2d %23.16e %23.16e %23.16e %23.16e %23.16e\n",
		Cell[i].Elem.W->harm[j],
		Cell[i].Elem.W->kxV[j], Cell[i].Elem.W->BoBrhoV[j],
		Cell[i].Elem.W->kxH[j], Cell[i].Elem.W->BoBrhoH[j],
		Cell[i].Elem.W->phi[j]);
      }
      break;
    case Cavity:
      prtName(mfile, i, cavity_, 0, 0);
      fprintf(mfile, " %23.16e %23.16e %d %23.16e\n",
	      Cell[i].Elem.C->Pvolt/(1e9*globval.Energy),
	      2.0*M_PI*Cell[i].Elem.C->Pfreq/c0, Cell[i].Elem.C->Ph,
	      1e9*globval.Energy);
      break;
    case marker:
      prtName(mfile, i, marker_, 0, 0);
      break;
    case Insertion:
      prtName(mfile, i, kick_map_, Cell[i].Elem.ID->Pmethod,
	      Cell[i].Elem.ID->PN);
      if (Cell[i].Elem.ID->firstorder)
	fprintf(mfile, " %3.1lf %1d %s\n",
		Cell[i].Elem.ID->scaling, 1, Cell[i].Elem.ID->fname1);
      else
	fprintf(mfile, " %3.1lf %1d %s\n",
		Cell[i].Elem.ID->scaling, 2, Cell[i].Elem.ID->fname2);
      break;
    default:
      printf("prtmfile: unknown type %d\n", Cell[i].Elem.Pkind);
      exit(1);
      break;
    }
  }
  fclose(mfile);
}
