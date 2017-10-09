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

     cavity:	 cavity voltage/beam energy [eV], omega/c, beam energy [eV],
                 phi

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
	  Lattice.Cell[i].Elem.Name,
	  Lattice.Cell[i].Fnum, Lattice.Cell[i].Knum, i);
  fprintf(fp, " %3d %3d %3d\n", type, method, N);
  fprintf(fp, " %23.16e %23.16e %23.16e %23.16e\n",
	  Lattice.Cell[i].maxampl[X_][0], Lattice.Cell[i].maxampl[X_][1],
	  Lattice.Cell[i].maxampl[Y_][0], Lattice.Cell[i].maxampl[Y_][1]);
}


void prtHOM(FILE *fp, const int n_design, const mpolArray B, const int Order)
{
  int i, nmpole;
  
  nmpole = 0;
  for (i = 1; i <= Order; i++)
    if ((B[HOMmax-i] != 0.0) || (B[HOMmax+i] != 0.0)) nmpole++;
  fprintf(fp, "  %2d %2d\n", nmpole, n_design);
  for (i = 1; i <= Order; i++) {
    if ((B[HOMmax-i] != 0.0) || (B[HOMmax+i] != 0.0))
      fprintf(fp, "%3d %23.16e %23.16e\n", i, B[HOMmax+i], B[HOMmax-i]);
  }
}


void Lattice_Type::prtmfile(const char mfile_dat[])
{
  int     i, j;
  FILE    *mfile;

  mfile = file_write(mfile_dat);
  for (i = 0; i <= globval.Cell_nLoc; i++) {
    switch (Lattice.Cell[i].Elem.Kind) {
    case drift:
      prtName(mfile, i, drift_, 0, 0);
      fprintf(mfile, " %23.16e\n", Lattice.Cell[i].Elem.L);
      break;
    case Mpole:
      if (Lattice.Cell[i].Elem.L != 0.0) {
	prtName(mfile, i, mpole_, Lattice.Cell[i].Elem.M->method,
		Lattice.Cell[i].Elem.M->N);
	fprintf(mfile, " %23.16e %23.16e %23.16e %23.16e\n",
		Lattice.Cell[i].dS[X_], Lattice.Cell[i].dS[Y_],
		Lattice.Cell[i].Elem.M->dTpar,
		Lattice.Cell[i].Elem.M->dTsys
		+Lattice.Cell[i].Elem.M->dTrms*Lattice.Cell[i].Elem.M->dTrnd);
	fprintf(mfile, " %23.16e %23.16e %23.16e %23.16e %23.16e\n",
		Lattice.Cell[i].Elem.L, Lattice.Cell[i].Elem.M->irho,
		Lattice.Cell[i].Elem.M->Tx1, Lattice.Cell[i].Elem.M->Tx2,
		Lattice.Cell[i].Elem.M->gap);
	prtHOM(mfile, Lattice.Cell[i].Elem.M->n_design,
	       Lattice.Cell[i].Elem.M->B, Lattice.Cell[i].Elem.M->order);
      } else {
	prtName(mfile, i, thinkick_, Lattice.Cell[i].Elem.M->method,
		Lattice.Cell[i].Elem.M->N);
	fprintf(mfile, " %23.16e %23.16e %23.16e\n",
		Lattice.Cell[i].dS[X_], Lattice.Cell[i].dS[Y_],
		Lattice.Cell[i].Elem.M->dTsys
		+Lattice.Cell[i].Elem.M->dTrms*Lattice.Cell[i].Elem.M->dTrnd);
	prtHOM(mfile, Lattice.Cell[i].Elem.M->n_design,
	       Lattice.Cell[i].Elem.M->B, Lattice.Cell[i].Elem.M->order);
      }
      break;
    case Wigl:
      prtName(mfile, i, wiggler_, Lattice.Cell[i].Elem.W->method,
	      Lattice.Cell[i].Elem.W->N);
      fprintf(mfile, " %23.16e %23.16e\n",
	      Lattice.Cell[i].Elem.L, Lattice.Cell[i].Elem.W->lambda);
      fprintf(mfile, "%2d\n", Lattice.Cell[i].Elem.W->n_harm);
      for (j = 0; j < Lattice.Cell[i].Elem.W->n_harm; j++) {
	fprintf(mfile, "%2d %23.16e %23.16e %23.16e %23.16e %23.16e\n",
		Lattice.Cell[i].Elem.W->harm[j],
		Lattice.Cell[i].Elem.W->kxV[j],
		Lattice.Cell[i].Elem.W->BoBrhoV[j],
		Lattice.Cell[i].Elem.W->kxH[j],
		Lattice.Cell[i].Elem.W->BoBrhoH[j],
		Lattice.Cell[i].Elem.W->phi[j]);
      }
      break;
    case Cavity:
      prtName(mfile, i, cavity_, 0, 0);
      fprintf(mfile, " %23.16e %23.16e %d %23.16e %23.16e\n",
	      Lattice.Cell[i].Elem.C->volt/(1e9*globval.Energy),
	      2.0*M_PI*Lattice.Cell[i].Elem.C->freq/c0,
	      Lattice.Cell[i].Elem.C->h,
	      1e9*globval.Energy, Lattice.Cell[i].Elem.C->phi);
      break;
    case marker:
      prtName(mfile, i, marker_, 0, 0);
      break;
    case Insertion:
      prtName(mfile, i, kick_map_, Lattice.Cell[i].Elem.ID->method,
	      Lattice.Cell[i].Elem.ID->N);
      if (Lattice.Cell[i].Elem.ID->firstorder)
	fprintf(mfile, " %3.1lf %1d %s\n",
		Lattice.Cell[i].Elem.ID->scaling, 1,
		Lattice.Cell[i].Elem.ID->fname1);
      else
	fprintf(mfile, " %3.1lf %1d %s\n",
		Lattice.Cell[i].Elem.ID->scaling, 2,
		Lattice.Cell[i].Elem.ID->fname2);
      break;
    default:
      printf("prtmfile: unknown type %d\n", Lattice.Cell[i].Elem.Kind);
      exit(1);
      break;
    }
  }
  fclose(mfile);
}
