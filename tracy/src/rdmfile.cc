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


void get_kind(const int kind, CellType *Cell)
{
  MpoleType *M;

  switch (kind) {
  case marker_:
    Cell->Kind = PartsKind(marker);
    Cell = new MarkerType();
    break;
  case drift_:
    Cell->Kind = PartsKind(drift);
    Cell = new DriftType();
    break;
  case mpole_:
    Cell->Kind = PartsKind(Mpole);
    M = new MpoleType();
    M->thick = pthicktype(thick);
    Cell = M;
    break;
  case cavity_:
    Cell->Kind = PartsKind(Cavity);
    Cell = new CavityType();
    break;
  case thinkick_:
    Cell->Kind = PartsKind(Mpole);
    M = new MpoleType();
    M->thick = pthicktype(thin);
    Cell = M;
    break;
  case wiggler_:
    Cell->Kind = PartsKind(Wigl);
    Cell = new WigglerType();
    break;
  case insertion_:
    Cell->Kind = PartsKind(Insertion);
    Cell = new InsertionType();
    break;
  default:
    std::cout << "get_kind: unknown type " << kind << " "
	      << Cell->Name << std::endl;
    exit_(1);
    break;
  }
}


void LatticeType::rdmfile(const char *mfile_dat)
{
  char          line[max_str], file_name[max_str];
  int           j, k, nmpole, kind, method, n;
  long int      i;
  double        dTerror;
  MpoleType     *M;
  WigglerType   *W;
  CavityType    *C;
  InsertionType *ID;

  bool prt = false;

  std::cout << std::endl;
  std::cout << "reading machine file: " << mfile_dat << std::endl;

  file_rd(inf, mfile_dat);

  while (inf.getline(line, max_str)) {
    if (prt) printf("%s\n", line);
    sscanf(line, "%*s %*d %*d %ld", &i);

    Lattice.Cell[i]->dS[X_] = 0.0; Lattice.Cell[i]->dS[Y_] = 0.0;
    Lattice.Cell[i]->dT[X_] = 1.0; Lattice.Cell[i]->dT[Y_] = 0.0;

    sscanf(line, "%s %d %d", Lattice.Cell[i]->Name,
	   &Lattice.Cell[i]->Fnum, &Lattice.Cell[i]->Knum);

    // For compability with lattice parser.
    k = 0;
    while (Lattice.Cell[i]->Name[k] != '\0')
      k++;
    for (j = k; j < SymbolLength; j++)
      Lattice.Cell[i]->Name[j] = ' ';

    if (Lattice.Cell[i]->Knum == 1) {
      strcpy(Lattice.ElemFam[Lattice.Cell[i]->Fnum-1].CellF->Name,
	     Lattice.Cell[i]->Name);
      Lattice.param.Elem_nFam =
	max(Lattice.Cell[i]->Fnum, Lattice.param.Elem_nFam);
    }

    if (i > 0) {
      Lattice.ElemFam[Lattice.Cell[i]->Fnum-1].KidList[Lattice.Cell[i]->Knum-1]
	= i;
      Lattice.ElemFam[Lattice.Cell[i]->Fnum-1].nKid =
	max(Lattice.Cell[i]->Knum,
	    Lattice.ElemFam[Lattice.Cell[i]->Fnum-1].nKid);
    }

    inf.getline(line, max_str);
    if (prt) printf("%s\n", line);
    sscanf(line, "%d %d %d", &kind, &method, &n);
    get_kind(kind, Lattice.Cell[i]);
    if (i > 0)
      Lattice.ElemFam[Lattice.Cell[i]->Fnum-1].CellF->Kind
	= Lattice.Cell[i]->Kind;

    inf.getline(line, max_str);
    if (prt) printf("%s\n", line);
    sscanf(line, "%lf %lf %lf %lf",
	   &Lattice.Cell[i]->maxampl[X_][0], &Lattice.Cell[i]->maxampl[X_][1],
	   &Lattice.Cell[i]->maxampl[Y_][0], &Lattice.Cell[i]->maxampl[Y_][1]);

    Lattice.Cell[i]->L = 0.0;

    switch (Lattice.Cell[i]->Kind) {
    case undef:
      std::cout << "rdmfile: unknown type " << i << std::endl;
      exit_(1);
      break;
    case marker:
      break;
    case drift:
      inf.getline(line, max_str);
      if (prt) printf("%s\n", line);
      sscanf(line, "%lf", &Lattice.Cell[i]->L);
      break;
    case Cavity:
      C = static_cast<CavityType*>(Lattice.Cell[i]);
      inf.getline(line, max_str);
      if (prt) printf("%s\n", line);
      sscanf(line, "%lf %lf %d %lf %lf",
	     &C->volt, &C->freq, &C->h, &Lattice.param.Energy, &C->phi);
      Lattice.param.Energy *= 1e-9;
      C->volt *= Lattice.param.Energy*1e9;
      C->freq *= c0/(2.0*M_PI);
     break;
    case Mpole:
      M = static_cast<MpoleType*>(Lattice.Cell[i]);
      M->method = method; M->N = n;

      if (M->thick == thick) {
	inf.getline(line, max_str);
	if (prt) printf("%s\n", line);
	sscanf(line, "%lf %lf %lf %lf",
	       &Lattice.Cell[i]->dS[X_], &Lattice.Cell[i]->dS[Y_],
	       &M->dTpar, &dTerror);
	Lattice.Cell[i]->dT[X_] = cos(dtor(dTerror+M->dTpar));
	Lattice.Cell[i]->dT[Y_] = sin(dtor(dTerror+M->dTpar));
	M->dTrms = dTerror - M->dTpar;M->dTrnd = 1e0;

	inf.getline(line, max_str);
	if (prt) printf("%s\n", line);
	sscanf(line, "%lf %lf %lf %lf %lf",
	       &Lattice.Cell[i]->L, &M->irho, &M->Tx1, &M->Tx2, &M->gap);
	if (M->irho != 0.0)
	  M->order = 1;
      } else {
	inf.getline(line, max_str);
	if (prt) printf("%s\n", line);
	sscanf(line, "%lf %lf %lf",
	       &Lattice.Cell[i]->dS[X_], &Lattice.Cell[i]->dS[Y_], &dTerror); 
	Lattice.Cell[i]->dT[X_] = cos(dtor(dTerror));
	Lattice.Cell[i]->dT[Y_] = sin(dtor(dTerror));
	M->dTrms = dTerror; M->dTrnd = 1e0;
      }

      M->c0 = sin(Lattice.Cell[i]->L*M->irho/2.0);
      M->c1 = cos(dtor(M->dTpar))*M->c0;
      M->s1 = sin(dtor(M->dTpar))*M->c0;

      inf.getline(line, max_str);
      if (prt) printf("%s\n", line);
      sscanf(line, "%d %d", &nmpole, &M->n_design);
      for (j = 1; j <= nmpole; j++) {
	inf.getline(line, max_str);
	if (prt) printf("%s\n", line);
	sscanf(line, "%d", &n);
	sscanf(line, "%*d %lf %lf",
	       &M->B[HOMmax+n], &M->B[HOMmax-n]);
	M->Bpar[HOMmax+n] = M->B[HOMmax+n]; M->Bpar[HOMmax-n] = M->B[HOMmax-n];
	M->order = max(n, M->order);
      }
      break;
    case Wigl:
      W = static_cast<WigglerType*>(Lattice.Cell[i]);
      W->method = method; W->N = n;

      inf.getline(line, max_str);
      if (prt) printf("%s\n", line);
      sscanf(line, "%lf %lf", &Lattice.Cell[i]->L, &W->lambda);

      inf.getline(line, max_str);
      if (prt) printf("%s\n", line);
      sscanf(line, "%d", &W->n_harm);

      // if (Lattice.Cell[i]->Knum == 1)
      // 	W->Wiggler_Alloc
      // 	  (&Lattice.ElemFam[Lattice.Cell[i]->Fnum-1].ElemF);
      for (j = 0; j < W->n_harm; j++) {
	inf.getline(line, max_str);
	if (prt) printf("%s\n", line);
	sscanf(line, "%d %lf %lf %lf %lf %lf",
	       &W->harm[j], &W->kxV[j], &W->BoBrhoV[j], &W->kxH[j],
	       &W->BoBrhoH[j], &W->phi[j]);
	W =
	  static_cast<WigglerType*>
	  (Lattice.ElemFam[Lattice.Cell[i]->Fnum-1].CellF);
	W->BoBrhoV[j] = W->BoBrhoV[j];W->BoBrhoH[j] = W->BoBrhoH[j];
      }
      break;
    case Insertion:
      ID = static_cast<InsertionType*>(Lattice.Cell[i]);
      ID->method = method;
      ID->N = n;

      inf.getline(line, max_str);
      if (prt) printf("%s\n", line);
      sscanf(line, "%lf %d %s", &ID->scaling, &n, file_name);

      if (n == 1) {
	ID->firstorder = true;
	ID->secondorder = false;

	strcpy(ID->fname1, file_name);
	Read_IDfile(ID->fname1, Lattice.Cell[i]->L, ID->nx, ID->nz,
		    ID->tabx, ID->tabz, ID->thetax1,
		    ID->thetaz1, ID->long_comp, ID->B2, ID->linear);
      } else if (n == 2) {
	ID->firstorder = false;
	ID->secondorder = true;

	strcpy(ID->fname2, file_name);
	Read_IDfile(ID->fname2, Lattice.Cell[i]->L, ID->nx, ID->nz,
		    ID->tabx, ID->tabz, ID->thetax, ID->thetaz,
		    ID->long_comp, ID->B2, ID->linear);
      } else {
	std::cout << "rdmfile: undef order " << n << std::endl;
	exit_(1);
      }

      if (ID->method == 1)
	ID->linear = true;
      else
	ID->linear = false;

      if (!ID->linear) {
	ID->tx = dmatrix(1, ID->nz, 1, ID->nx);
	ID->tz = dmatrix(1, ID->nz, 1, ID->nx);
	ID->tab1
	  = (double *)malloc((ID->nx)*sizeof(double));
	ID->tab2
	  = (double *)malloc((ID->nz)*sizeof(double));
	ID->f2x = dmatrix(1, ID->nz, 1, ID->nx);
	ID->f2z = dmatrix(1, ID->nz, 1, ID->nx);
	Matrices4Spline(ID);
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
      Lattice.Cell[i]->S = 0.0;
    else
      Lattice.Cell[i]->S = Lattice.Cell[i-1]->S + Lattice.Cell[i]->L;
  }
  
  Lattice.param.Cell_nLoc = i;
 
  Lattice.param.dPcommon = 1e-8;
  Lattice.param.CODeps = 1e-14;
  Lattice.param.CODimax = 40;

  SI_init();

  std::cout << std::endl;
  std::cout  << std::fixed << std::setprecision(5)
	<< "rdmfile: read " << Lattice.param.Cell_nLoc << " elements, C = "
	<< Lattice.Cell[Lattice.param.Cell_nLoc]->S << std::endl;

  inf.close();
}
