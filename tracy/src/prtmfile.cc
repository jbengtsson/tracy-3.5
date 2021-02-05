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
     map         7

   Integration methods:

     fourth order symplectic integrator   4

   Format:

     name, family no, kid no, element no
     type code, integration method, no of integration steps, reverse
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

     wiggler:    L [m], Lambda [m]
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

     map:

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
#define map_       7


void prtName(FILE *fp, const int i,
	     const int type, const int method, const int N, const bool reverse)
{
  fprintf(fp, "%-15s %4d %4d %4d\n",
	  lat.elems[i]->PName, lat.elems[i]->Fnum, lat.elems[i]->Knum, i);
  fprintf(fp, " %3d %3d %3d %4d\n", type, method, N, reverse);
  fprintf(fp, " %23.16e %23.16e %23.16e %23.16e\n",
	  lat.elems[i]->maxampl[X_][0], lat.elems[i]->maxampl[X_][1],
	  lat.elems[i]->maxampl[Y_][0], lat.elems[i]->maxampl[Y_][1]);
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
  int           i, j, k;
  MpoleType     *M;
  WigglerType   *W;
  CavityType    *C;
  InsertionType *ID;
  MapType       *Mapp;
  FILE          *mfile;

  const int n_ps = 6;

  mfile = file_write(mfile_dat);
  for (i = 0; i <= globval.Cell_nLoc; i++) {
    switch (lat.elems[i]->Pkind) {
    case drift:
      prtName(mfile, i, drift_, 0, 0, 0);
      fprintf(mfile, " %23.16e\n", lat.elems[i]->PL);
      break;
    case Mpole:
      M = dynamic_cast<MpoleType*>(lat.elems[i]);
      if (lat.elems[i]->PL != 0.0) {
	prtName(mfile, i, mpole_, M->Pmethod, M->PN,
		lat.elems[i]->Reverse);
	fprintf(mfile, " %23.16e %23.16e %23.16e %23.16e\n",
		lat.elems[i]->dS[X_], lat.elems[i]->dS[Y_],
		M->PdTpar, M->PdTsys+M->PdTrms*M->PdTrnd);
	fprintf(mfile, " %23.16e %23.16e %23.16e %23.16e %23.16e\n",
		lat.elems[i]->PL, M->Pirho, M->PTx1, M->PTx2, M->Pgap);
	prtHOM(mfile, M->n_design, M->PB, M->Porder);
      } else {
	prtName(mfile, i, thinkick_, M->Pmethod, M->PN, lat.elems[i]->Reverse);
	fprintf(mfile, " %23.16e %23.16e %23.16e\n",
		lat.elems[i]->dS[X_], lat.elems[i]->dS[Y_],
		M->PdTsys+M->PdTrms*M->PdTrnd);
	prtHOM(mfile, M->n_design, M->PB, M->Porder);
      }
      break;
    case Wigl:
      W = dynamic_cast<WigglerType*>(lat.elems[i]);
      prtName(mfile, i, wiggler_, W->Pmethod, W->PN, lat.elems[i]->Reverse);
      fprintf(mfile, " %23.16e %23.16e\n", lat.elems[i]->PL, W->Lambda);
      fprintf(mfile, "%2d\n", W->n_harm);
      for (j = 0; j < W->n_harm; j++) {
	fprintf(mfile, "%2d %23.16e %23.16e %23.16e %23.16e %23.16e\n",
		W->harm[j], W->kxV[j], W->BoBrhoV[j], W->kxH[j], W->BoBrhoH[j],
		W->phi[j]);
      }
      break;
    case Cavity:
      C = dynamic_cast<CavityType*>(lat.elems[i]);
      prtName(mfile, i, cavity_, 0, 0, 0);
      fprintf(mfile, " %23.16e %23.16e %d %23.16e %23.16e\n",
	      C->Pvolt/(1e9*globval.Energy), 2.0*M_PI*C->Pfreq/c0, C->Ph,
	      1e9*globval.Energy, C->phi);
      break;
    case marker:
      prtName(mfile, i, marker_, 0, 0, 0);
      break;
    case Insertion:
      ID = dynamic_cast<InsertionType*>(lat.elems[i]);
      prtName(mfile, i, kick_map_, ID->Pmethod, ID->PN, lat.elems[i]->Reverse);
      if (ID->firstorder)
	fprintf(mfile, " %3.1lf %1d %s\n", ID->scaling, 1, ID->fname1);
      else
	fprintf(mfile, " %3.1lf %1d %s\n", ID->scaling, 2, ID->fname2);
      break;
    case Map:
      Mapp = dynamic_cast<MapType*>(lat.elems[i]);
      prtName(mfile, i, map_, 0, 0, 0);
      for (j = 0; j < n_ps; j++) {
	for (k = 0; k < n_ps; k++)
	  fprintf(mfile, " %23.16le", Mapp->M[j][k]);
	fprintf(mfile, "\n");
      }
      break;
    default:
      printf("prtmfile: unknown type %d\n", lat.elems[i]->Pkind);
      exit(1);
      break;
    }
  }
  fclose(mfile);
}
