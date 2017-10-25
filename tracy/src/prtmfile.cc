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
	  Lattice.Cell[i]->Name,
	  Lattice.Cell[i]->Fnum, Lattice.Cell[i]->Knum, i);
  fprintf(fp, " %3d %3d %3d\n", type, method, N);
  fprintf(fp, " %23.16e %23.16e %23.16e %23.16e\n",
	  Lattice.Cell[i]->maxampl[X_][0], Lattice.Cell[i]->maxampl[X_][1],
	  Lattice.Cell[i]->maxampl[Y_][0], Lattice.Cell[i]->maxampl[Y_][1]);
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


void LatticeType::prtmfile(const char mfile_dat[])
{
  int           i, j;
  MpoleType     *M;
  WigglerType   *W;
  CavityType    *C;
  InsertionType *ID;
  FILE          *mfile;

  mfile = file_write(mfile_dat);
  for (i = 0; i <= Lattice.param.Cell_nLoc; i++) {
    printf("%8s\n", Lattice.Cell[i]->Name);
    switch (Lattice.Cell[i]->Kind) {
    case drift:
      prtName(mfile, i, drift_, 0, 0);
      fprintf(mfile, " %23.16e\n", Lattice.Cell[i]->L);
      break;
    case Mpole:
      M = static_cast<MpoleType*>(Lattice.Cell[i]);
      if (Lattice.Cell[i]->L != 0.0) {
	prtName(mfile, i, mpole_, M->method, M->N);
	fprintf(mfile, " %23.16e %23.16e %23.16e %23.16e\n",
		Lattice.Cell[i]->dS[X_], Lattice.Cell[i]->dS[Y_],
		M->dTpar, M->dTsys+M->dTrms*M->dTrnd);
	fprintf(mfile, " %23.16e %23.16e %23.16e %23.16e %23.16e\n",
		Lattice.Cell[i]->L, M->irho, M->Tx1, M->Tx2,	M->gap);
	prtHOM(mfile, M->n_design, M->B, M->order);
      } else {
	prtName(mfile, i, thinkick_, M->method,	M->N);
	fprintf(mfile, " %23.16e %23.16e %23.16e\n",
		Lattice.Cell[i]->dS[X_], Lattice.Cell[i]->dS[Y_],
		M->dTsys+M->dTrms*M->dTrnd);
	prtHOM(mfile, M->n_design, M->B, M->order);
      }
      break;
    case Wigl:
      W = static_cast<WigglerType*>(Lattice.Cell[i]);
      prtName(mfile, i, wiggler_, W->method, W->N);
      fprintf(mfile, " %23.16e %23.16e\n",
	      Lattice.Cell[i]->L, W->lambda);
      fprintf(mfile, "%2d\n", W->n_harm);
      for (j = 0; j < W->n_harm; j++) {
	fprintf(mfile, "%2d %23.16e %23.16e %23.16e %23.16e %23.16e\n",
		W->harm[j], W->kxV[j], W->BoBrhoV[j], W->kxH[j], W->BoBrhoH[j],
		W->phi[j]);
      }
      break;
    case Cavity:
      C = static_cast<CavityType*>(Lattice.Cell[i]);
      prtName(mfile, i, cavity_, 0, 0);
      fprintf(mfile, " %23.16e %23.16e %d %23.16e %23.16e\n",
	      C->volt/(1e9*Lattice.param.Energy), 2.0*M_PI*C->freq/c0, C->h,
	      1e9*Lattice.param.Energy, C->phi);
      break;
    case marker:
      prtName(mfile, i, marker_, 0, 0);
      break;
    case Insertion:
      ID = static_cast<InsertionType*>(Lattice.Cell[i]);
      prtName(mfile, i, kick_map_, ID->method, ID->N);
      if (ID->firstorder)
	fprintf(mfile, " %3.1lf %1d %s\n", ID->scaling, 1, ID->fname1);
      else
	fprintf(mfile, " %3.1lf %1d %s\n", ID->scaling, 2, ID->fname2);
      break;
    default:
      printf("prtmfile: unknown type %d\n", Lattice.Cell[i]->Kind);
      exit(1);
      break;
    }
  }
  fclose(mfile);
}
