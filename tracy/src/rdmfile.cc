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

   Integration methods:

     linear, matrix style (obsolete)              0
     2nd order symplectic integrator (obsolete)   2
     4th order symplectic integrator              4

   Format:

     name, family no, kid no, element no
     type code, integration method, no of integration steps
     apertures: xmin, xmax, ymin, ymax

   The following lines follows depending on element type.

     type

     drift:	 L

     multipole:  hor., ver. displ., roll angle (design), roll angle (error)
                 L, 1/rho, entrance angle, exit angle
		 no of nonzero multipole coeff., n design
		 n, b , a
		     n   n
		     .
		     .
		     .

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
		     .
		     .
		     .

     kick_map:   scale order <file name>

*/


#define snamelen   10

// numerical type codes
#define marker_   -1
#define drift_     0
#define mpole_     1
#define cavity_    2
#define thinkick_  3
#define wiggler_   4
#define insertion_ 6


std::ifstream  inf;


void get_kind(const int kind, elemtype &Elem)
{

  switch (kind) {
  case marker_:
    Elem.Kind = PartsKind(marker);
    break;
  case drift_:
    Elem.Kind = PartsKind(drift);
    Drift_Alloc(&Elem);
    break;
  case mpole_:
    Elem.Kind = PartsKind(Mpole);
    Mpole_Alloc(&Elem);
    Elem.M->thick = pthicktype(thick);
    break;
  case cavity_:
    Elem.Kind = PartsKind(Cavity);
    Cav_Alloc(&Elem);
    break;
  case thinkick_:
    Elem.Kind = PartsKind(Mpole);
    Mpole_Alloc(&Elem);
    Elem.M->thick = pthicktype(thin);
    break;
  case wiggler_:
    Elem.Kind = PartsKind(Wigl);
    Wiggler_Alloc(&Elem);
    break;
  case insertion_:
    Elem.Kind = PartsKind(Insertion);
    Insertion_Alloc(&Elem);
    break;
  default:
    std::cout << "get_kind: unknown type " << kind << " "
	      << Elem.Name << std::endl;
    exit_(1);
    break;
  }
}


void LatticeType::rdmfile(const char *mfile_dat)
{
  char      line[max_str], file_name[max_str];
  int       j, k, nmpole, kind, method, n;
  long int  i;
  double    dTerror;

  bool  prt = false;

  std::cout << std::endl;
  std::cout << "reading machine file: " << mfile_dat << std::endl;

  file_rd(inf, mfile_dat);

  while (inf.getline(line, max_str)) {
    if (prt) printf("%s\n", line);
    sscanf(line, "%*s %*d %*d %ld", &i);

    Lattice.Cell[i].dS[X_] = 0.0; Lattice.Cell[i].dS[Y_] = 0.0;
    Lattice.Cell[i].dT[X_] = 1.0; Lattice.Cell[i].dT[Y_] = 0.0;

    sscanf(line, "%s %d %d", Lattice.Cell[i].Elem.Name,
	   &Lattice.Cell[i].Fnum, &Lattice.Cell[i].Knum);

    // For compability with lattice parser.
    k = 0;
    while (Lattice.Cell[i].Elem.Name[k] != '\0')
      k++;
    for (j = k; j < SymbolLength; j++)
      Lattice.Cell[i].Elem.Name[j] = ' ';

    if (Lattice.Cell[i].Knum == 1) {
      strcpy(Lattice.ElemFam[Lattice.Cell[i].Fnum-1].ElemF.Name,
	     Lattice.Cell[i].Elem.Name);
      Lattice.param.Elem_nFam = max(Lattice.Cell[i].Fnum, Lattice.param.Elem_nFam);
    }

    if (i > 0) {
      Lattice.ElemFam[Lattice.Cell[i].Fnum-1].KidList[Lattice.Cell[i].Knum-1]
	= i;
      Lattice.ElemFam[Lattice.Cell[i].Fnum-1].nKid =
	max(Lattice.Cell[i].Knum, Lattice.ElemFam[Lattice.Cell[i].Fnum-1].nKid);
    }

    inf.getline(line, max_str);
    if (prt) printf("%s\n", line);
    sscanf(line, "%d %d %d", &kind, &method, &n);
    get_kind(kind, Lattice.Cell[i].Elem);
    if (i > 0)
      Lattice.ElemFam[Lattice.Cell[i].Fnum-1].ElemF.Kind
	= Lattice.Cell[i].Elem.Kind;

    inf.getline(line, max_str);
    if (prt) printf("%s\n", line);
    sscanf(line, "%lf %lf %lf %lf",
	   &Lattice.Cell[i].maxampl[X_][0], &Lattice.Cell[i].maxampl[X_][1],
	   &Lattice.Cell[i].maxampl[Y_][0], &Lattice.Cell[i].maxampl[Y_][1]);

    Lattice.Cell[i].Elem.L = 0.0;

    switch (Lattice.Cell[i].Elem.Kind) {
    case undef:
      std::cout << "rdmfile: unknown type " << i << std::endl;
      exit_(1);
      break;
    case marker:
      break;
    case drift:
      inf.getline(line, max_str);
      if (prt) printf("%s\n", line);
      sscanf(line, "%lf", &Lattice.Cell[i].Elem.L);
      break;
    case Cavity:
      inf.getline(line, max_str);
      if (prt) printf("%s\n", line);
      sscanf(line, "%lf %lf %d %lf %lf",
	     &Lattice.Cell[i].Elem.C->volt, &Lattice.Cell[i].Elem.C->freq,
	     &Lattice.Cell[i].Elem.C->h, &Lattice.param.Energy,
	     &Lattice.Cell[i].Elem.C->phi);
      Lattice.param.Energy *= 1e-9;
      Lattice.Cell[i].Elem.C->volt *= Lattice.param.Energy*1e9;
      Lattice.Cell[i].Elem.C->freq *= c0/(2.0*M_PI);
     break;
    case Mpole:
      Lattice.Cell[i].Elem.M->method = method; Lattice.Cell[i].Elem.M->N = n;

      if (Lattice.Cell[i].Elem.M->thick == thick) {
	inf.getline(line, max_str);
	if (prt) printf("%s\n", line);
	sscanf(line, "%lf %lf %lf %lf",
	       &Lattice.Cell[i].dS[X_], &Lattice.Cell[i].dS[Y_],
	       &Lattice.Cell[i].Elem.M->dTpar, &dTerror);
	Lattice.Cell[i].dT[X_]
	  = cos(dtor(dTerror+Lattice.Cell[i].Elem.M->dTpar));
	Lattice.Cell[i].dT[Y_]
	  = sin(dtor(dTerror+Lattice.Cell[i].Elem.M->dTpar));
	Lattice.Cell[i].Elem.M->dTrms
	  = dTerror - Lattice.Cell[i].Elem.M->dTpar;
	Lattice.Cell[i].Elem.M->dTrnd = 1e0;

	inf.getline(line, max_str);
	if (prt) printf("%s\n", line);
	sscanf(line, "%lf %lf %lf %lf %lf",
	       &Lattice.Cell[i].Elem.L, &Lattice.Cell[i].Elem.M->irho,
	       &Lattice.Cell[i].Elem.M->Tx1, &Lattice.Cell[i].Elem.M->Tx2,
	       &Lattice.Cell[i].Elem.M->gap);
	if (Lattice.Cell[i].Elem.M->irho != 0.0)
	  Lattice.Cell[i].Elem.M->order = 1;
      } else {
	inf.getline(line, max_str);
	if (prt) printf("%s\n", line);
	sscanf(line, "%lf %lf %lf",
	       &Lattice.Cell[i].dS[X_], &Lattice.Cell[i].dS[Y_], &dTerror); 
	Lattice.Cell[i].dT[X_] = cos(dtor(dTerror));
	Lattice.Cell[i].dT[Y_] = sin(dtor(dTerror));
	Lattice.Cell[i].Elem.M->dTrms = dTerror;
	Lattice.Cell[i].Elem.M->dTrnd = 1e0;
      }

      Lattice.Cell[i].Elem.M->c0
	= sin(Lattice.Cell[i].Elem.L*Lattice.Cell[i].Elem.M->irho/2.0);
      Lattice.Cell[i].Elem.M->c1
	= cos(dtor(Lattice.Cell[i].Elem.M->dTpar))
	                    *Lattice.Cell[i].Elem.M->c0;
      Lattice.Cell[i].Elem.M->s1 = sin(dtor(Lattice.Cell[i].Elem.M->dTpar))
	                    *Lattice.Cell[i].Elem.M->c0;

      inf.getline(line, max_str);
      if (prt) printf("%s\n", line);
      sscanf(line, "%d %d", &nmpole, &Lattice.Cell[i].Elem.M->n_design);
      for (j = 1; j <= nmpole; j++) {
	inf.getline(line, max_str);
	if (prt) printf("%s\n", line);
	sscanf(line, "%d", &n);
	sscanf(line, "%*d %lf %lf",
	       &Lattice.Cell[i].Elem.M->B[HOMmax+n],
	       &Lattice.Cell[i].Elem.M->B[HOMmax-n]);
	Lattice.Cell[i].Elem.M->Bpar[HOMmax+n]
	  = Lattice.Cell[i].Elem.M->B[HOMmax+n];
	Lattice.Cell[i].Elem.M->Bpar[HOMmax-n]
	  = Lattice.Cell[i].Elem.M->B[HOMmax-n];
	Lattice.Cell[i].Elem.M->order = max(n, Lattice.Cell[i].Elem.M->order);
      }
      break;
    case Wigl:
      Lattice.Cell[i].Elem.W->method = method; Lattice.Cell[i].Elem.W->N = n;

      inf.getline(line, max_str);
      if (prt) printf("%s\n", line);
      sscanf(line, "%lf %lf",
	     &Lattice.Cell[i].Elem.L, &Lattice.Cell[i].Elem.W->lambda);

      inf.getline(line, max_str);
      if (prt) printf("%s\n", line);
      sscanf(line, "%d", &Lattice.Cell[i].Elem.W->n_harm);

      if (Lattice.Cell[i].Knum == 1)
	Wiggler_Alloc(&Lattice.ElemFam[Lattice.Cell[i].Fnum-1].ElemF);
      for (j = 0; j < Lattice.Cell[i].Elem.W->n_harm; j++) {
	inf.getline(line, max_str);
	if (prt) printf("%s\n", line);
	sscanf(line, "%d %lf %lf %lf %lf %lf",
	       &Lattice.Cell[i].Elem.W->harm[j],
	       &Lattice.Cell[i].Elem.W->kxV[j],
	       &Lattice.Cell[i].Elem.W->BoBrhoV[j],
	       &Lattice.Cell[i].Elem.W->kxH[j],
	       &Lattice.Cell[i].Elem.W->BoBrhoH[j],
	       &Lattice.Cell[i].Elem.W->phi[j]);
	Lattice.ElemFam[Lattice.Cell[i].Fnum-1].ElemF.W->BoBrhoV[j]
	  = Lattice.Cell[i].Elem.W->BoBrhoV[j];
	Lattice.ElemFam[Lattice.Cell[i].Fnum-1].ElemF.W->BoBrhoH[j]
	  = Lattice.Cell[i].Elem.W->BoBrhoH[j];
      }
      break;
    case Insertion:
      Lattice.Cell[i].Elem.ID->method = method;
      Lattice.Cell[i].Elem.ID->N = n;

      inf.getline(line, max_str);
      if (prt) printf("%s\n", line);
      sscanf(line, "%lf %d %s",
	     &Lattice.Cell[i].Elem.ID->scaling, &n, file_name);

      if (n == 1) {
	Lattice.Cell[i].Elem.ID->firstorder = true;
	Lattice.Cell[i].Elem.ID->secondorder = false;

	strcpy(Lattice.Cell[i].Elem.ID->fname1, file_name);
	Read_IDfile(Lattice.Cell[i].Elem.ID->fname1, Lattice.Cell[i].Elem.L,
		    Lattice.Cell[i].Elem.ID->nx, Lattice.Cell[i].Elem.ID->nz,
		    Lattice.Cell[i].Elem.ID->tabx,
		    Lattice.Cell[i].Elem.ID->tabz,
		    Lattice.Cell[i].Elem.ID->thetax1,
		    Lattice.Cell[i].Elem.ID->thetaz1,
		    Lattice.Cell[i].Elem.ID->long_comp,
		    Lattice.Cell[i].Elem.ID->B2,
		    Lattice.Cell[i].Elem.ID->linear);
      } else if (n == 2) {
	Lattice.Cell[i].Elem.ID->firstorder = false;
	Lattice.Cell[i].Elem.ID->secondorder = true;

	strcpy(Lattice.Cell[i].Elem.ID->fname2, file_name);
	Read_IDfile(Lattice.Cell[i].Elem.ID->fname2, Lattice.Cell[i].Elem.L,
		    Lattice.Cell[i].Elem.ID->nx, Lattice.Cell[i].Elem.ID->nz,
		    Lattice.Cell[i].Elem.ID->tabx,
		    Lattice.Cell[i].Elem.ID->tabz,
		    Lattice.Cell[i].Elem.ID->thetax,
		    Lattice.Cell[i].Elem.ID->thetaz,
		    Lattice.Cell[i].Elem.ID->long_comp,
		    Lattice.Cell[i].Elem.ID->B2,
		    Lattice.Cell[i].Elem.ID->linear);
      } else {
	std::cout << "rdmfile: undef order " << n << std::endl;
	exit_(1);
      }

      if (Lattice.Cell[i].Elem.ID->method == 1)
	Lattice.Cell[i].Elem.ID->linear = true;
      else
	Lattice.Cell[i].Elem.ID->linear = false;

      if (!Lattice.Cell[i].Elem.ID->linear) {
	Lattice.Cell[i].Elem.ID->tx = dmatrix(1, Lattice.Cell[i].Elem.ID->nz,
				      1, Lattice.Cell[i].Elem.ID->nx);
	Lattice.Cell[i].Elem.ID->tz = dmatrix(1, Lattice.Cell[i].Elem.ID->nz,
				      1, Lattice.Cell[i].Elem.ID->nx);
	Lattice.Cell[i].Elem.ID->tab1
	  = (double *)malloc((Lattice.Cell[i].Elem.ID->nx)*sizeof(double));
	Lattice.Cell[i].Elem.ID->tab2
	  = (double *)malloc((Lattice.Cell[i].Elem.ID->nz)*sizeof(double));
	Lattice.Cell[i].Elem.ID->f2x = dmatrix(1, Lattice.Cell[i].Elem.ID->nz,
				       1, Lattice.Cell[i].Elem.ID->nx);
	Lattice.Cell[i].Elem.ID->f2z = dmatrix(1, Lattice.Cell[i].Elem.ID->nz,
				       1, Lattice.Cell[i].Elem.ID->nx);
	Matrices4Spline(Lattice.Cell[i].Elem.ID);
      }

/*      free_matrix(tx, 1, nz, 1, nx); free_matrix(tz, 1, nz, 1, nx);
      free(tab1); free(tab2);
      free_matrix(f2x, 1, nz, 1, nx); free_matrix(f2z, 1, nz, 1, nx); */
      break;
    case FieldMap:
      break;
    default:
      std::cout << "rdmfile: unknown type" << std::endl;
      exit_(1);
      break;
    }

    if (i == 0)
      Lattice.Cell[i].S = 0.0;
    else
      Lattice.Cell[i].S = Lattice.Cell[i-1].S + Lattice.Cell[i].Elem.L;
  }
  
  Lattice.param.Cell_nLoc = i;
 
  Lattice.param.dPcommon = 1e-8;
  Lattice.param.CODeps = 1e-14;
  Lattice.param.CODimax = 40;

  SI_init();

  std::cout << std::endl;
  std::cout  << std::fixed << std::setprecision(5)
	<< "rdmfile: read " << Lattice.param.Cell_nLoc << " elements, C = "
	<< Lattice.Cell[Lattice.param.Cell_nLoc].S << std::endl;

  inf.close();
}
