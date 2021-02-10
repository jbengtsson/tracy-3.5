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
     type code, integration method, no of integration steps, reverse
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
#define kick_map   6
#define map_       7


const int line_max = 200;


ElemType* elem_alloc(const int kind)
{
  ElemType   *Elem;
  MpoleType  *M;
  CavityType *C;

  switch (kind) {
  case marker_:
    Elem = Marker_Alloc();
    Elem->Pkind = PartsKind(marker);
    break;
  case drift_:
    Elem = Drift_Alloc();
    Elem->Pkind = PartsKind(drift);
    break;
  case mpole_:
    Elem = Mpole_Alloc();
    Elem->Pkind = PartsKind(Mpole);
    M = dynamic_cast<MpoleType*>(Elem);
    M->Pthick = pthicktype(thick);
    break;
  case cavity_:
    Elem = Cavity_Alloc();
    Elem->Pkind = PartsKind(Cavity);
    C = dynamic_cast<CavityType*>(Elem);
    break;
  case thinkick_:
    Elem = Mpole_Alloc();
    Elem->Pkind = PartsKind(Mpole);
    M = dynamic_cast<MpoleType*>(Elem);
    M->Pthick = pthicktype(thin);
    break;
  case wiggler_:
    Elem = Wiggler_Alloc();
    Elem->Pkind = PartsKind(Wigl);
    break;
  case kick_map:
    Elem = Insertion_Alloc();
    Elem->Pkind = PartsKind(Insertion);
    break;
  case map_:
    Elem = Map_Alloc();
    Elem->Pkind = PartsKind(Map);
    break;
  default:
    printf("elem_alloc: unknown type %d", kind);
    Elem = NULL;
    exit_(1);
    break;
  }
  return Elem;
}


void get_elem(std::ifstream &inf, char *line, long int &i, int &kind)
{
  char file_name[line_max];

  const int n_ps = 6;

  int           Fnum, Knum, j, k, nmpole, method, n, reverse;
  double        dTerror, val[n_ps];
  ss_vect<tps>  Id;
  partsName     name;
  MpoleType     *M;
  CavityType    *C;
  WigglerType   *W;
  InsertionType *ID;
  MapType       *Mapp;

  const bool prt = false;

  if (prt) printf("%s\n", line);

  sscanf(line, "%*s %*d %*d %ld", &i);
  sscanf(line, "%s %d %d", name, &Fnum, &Knum);

  // For compability with lattice parser.
  k = 0;
  while (name[k] != '\0')
    k++;
  for (j = k; j < SymbolLength; j++)
    name[j] = ' ';

  inf.getline(line, line_max);
  if (prt) printf("%s\n", line);
  sscanf(line, "%d %d %d %d", &kind, &method, &n, &reverse);

  lat.elems[i] = elem_alloc(kind);
 
  memcpy(lat.elems[i]->PName, name, sizeof(partsName));
  lat.elems[i]->Fnum = Fnum; lat.elems[i]->Knum = Knum;

  lat.elems[i]->dS[X_] = 0e0; lat.elems[i]->dS[Y_] = 0e0;
  lat.elems[i]->dT[X_] = 1e0; lat.elems[i]->dT[Y_] = 0e0;

  inf.getline(line, line_max);
  if (prt) printf("%s\n", line);
  sscanf(line, "%lf %lf %lf %lf",
	 &lat.elems[i]->maxampl[X_][0], &lat.elems[i]->maxampl[X_][1],
	 &lat.elems[i]->maxampl[Y_][0], &lat.elems[i]->maxampl[Y_][1]);

  lat.elems[i]->PL = 0e0;
  lat.elems[i]->Reverse = (reverse == 1);

  switch (lat.elems[i]->Pkind) {
  case undef:
    std::cout << "rdmfile: unknown type " << i << std::endl;
    exit_(1);
    break;
  case marker:
    break;
  case drift:
    inf.getline(line, line_max);
    if (prt) printf("%s\n", line);
    sscanf(line, "%lf", &lat.elems[i]->PL);
    break;
  case Cavity:
    C = dynamic_cast<CavityType*>(lat.elems[i]);
    inf.getline(line, line_max);
    if (prt) printf("%s\n", line);
    sscanf(line, "%lf %lf %d %lf %lf",
	   &C->Pvolt, &C->Pfreq, &C->Ph, &globval.Energy, &C->phi);
    globval.Energy *= 1e-9;
    C->Pvolt *= globval.Energy*1e9;
    C->Pfreq *= c0/(2.0*M_PI);
    break;
  case Mpole:
    M = dynamic_cast<MpoleType*>(lat.elems[i]);
    M->Pmethod = method; M->PN = n;

    if (M->Pthick == thick) {
      inf.getline(line, line_max);
      if (prt) printf("%s\n", line);
      sscanf(line, "%lf %lf %lf %lf",
	     &lat.elems[i]->dS[X_], &lat.elems[i]->dS[Y_],
	     &M->PdTpar, &dTerror);
      lat.elems[i]->dT[X_] = cos(dtor(dTerror+M->PdTpar));
      lat.elems[i]->dT[Y_] = sin(dtor(dTerror+M->PdTpar));
      M->PdTrms = dTerror - M->PdTpar;
      M->PdTrnd = 1e0;

      inf.getline(line, line_max);
      if (prt) printf("%s\n", line);
      sscanf(line, "%lf %lf %lf %lf %lf",
	     &lat.elems[i]->PL, &M->Pirho, &M->PTx1, &M->PTx2, &M->Pgap);
      if (M->Pirho != 0e0) M->Porder = 1;
    } else {
      inf.getline(line, line_max);
      if (prt) printf("%s\n", line);
      sscanf(line, "%lf %lf %lf",
	     &lat.elems[i]->dS[X_], &lat.elems[i]->dS[Y_], &dTerror); 
      lat.elems[i]->dT[X_] = cos(dtor(dTerror));
      lat.elems[i]->dT[Y_] = sin(dtor(dTerror));
      M->PdTrms = dTerror; M->PdTrnd = 1e0;
    }

    M->Pc0 = sin(lat.elems[i]->PL*M->Pirho/2.0);
    M->Pc1 = cos(dtor(M->PdTpar))*M->Pc0;
    M->Ps1 = sin(dtor(M->PdTpar))*M->Pc0;

    inf.getline(line, line_max);
    if (prt) printf("%s\n", line);
    sscanf(line, "%d %d", &nmpole, &M->n_design);
    for (j = 1; j <= nmpole; j++) {
      inf.getline(line, line_max);
      if (prt) printf("%s\n", line);
      sscanf(line, "%d", &n);
      sscanf(line, "%*d %lf %lf", &M->PB[HOMmax+n], &M->PB[HOMmax-n]);
      M->PBpar[HOMmax+n] = M->PB[HOMmax+n];
      M->PBpar[HOMmax-n] = M->PB[HOMmax-n];
      M->Porder = max(n, M->Porder);
    }

    if (globval.mat_meth && (M->Pthick == thick))
      M->M_lin = get_lin_map(lat.elems[i], 0e0);
    break;
  case Wigl:
    W = dynamic_cast<WigglerType*>(lat.elems[i]);
    W->Pmethod = method; W->PN = n;

    inf.getline(line, line_max);
    if (prt) printf("%s\n", line);
    sscanf(line, "%lf %lf", &lat.elems[i]->PL, &W->Lambda);

    inf.getline(line, line_max);
    if (prt) printf("%s\n", line);
    sscanf(line, "%d", &W->n_harm);

    for (j = 0; j < W->n_harm; j++) {
      inf.getline(line, line_max);
      if (prt) printf("%s\n", line);
      sscanf(line, "%d %lf %lf %lf %lf %lf", &W->harm[j],
	     &W->kxV[j], &W->BoBrhoV[j], &W->kxH[j], &W->BoBrhoH[j],
	     &W->phi[j]);
    }
    break;
  case Insertion:
    ID = dynamic_cast<InsertionType*>(lat.elems[i]);
    ID->Pmethod = method; ID->PN = n;

    inf.getline(line, line_max);
    if (prt) printf("%s\n", line);
    sscanf(line, "%lf %d %s", &ID->scaling, &n, file_name);

    if (n == 1) {
      ID->firstorder = true;
      ID->secondorder = false;

      strcpy(ID->fname1, file_name);
      Read_IDfile(ID->fname1, lat.elems[i]->PL, ID->nx, ID->nz,
		  ID->tabx, ID->tabz, ID->thetax1, ID->thetaz1,
		  ID->long_comp, ID->B2);
    } else if (n == 2) {
      ID->firstorder = false;
      ID->secondorder = true;

      strcpy(ID->fname2, file_name);
      Read_IDfile(ID->fname2, lat.elems[i]->PL, ID->nx, ID->nz,
		  ID->tabx, ID->tabz, ID->thetax, ID->thetaz,
		  ID->long_comp, ID->B2);
    } else {
      std::cout << "rdmfile: undef order " << n << std::endl;
      exit_(1);
    }

    if (ID->Pmethod == 1)
      ID->linear = true;
    else
      ID->linear = false;

    if (!ID->linear) {
      ID->tx = dmatrix(1, ID->nz,
		       1, ID->nx);
      ID->tz = dmatrix(1, ID->nz,
		       1, ID->nx);
      ID->tab1 = (double *)malloc((ID->nx)*sizeof(double));
      ID->tab2 = (double *)malloc((ID->nz)*sizeof(double));
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
  case Map:
    Mapp = dynamic_cast<MapType*>(lat.elems[i]);
    Id.identity(); Mapp->M.zero();
    for (j = 0; j < n_ps; j++) {
      inf.getline(line, line_max);
      if (prt) printf("%s\n", line);
      sscanf(line, "%lf %lf %lf %lf %lf %lf",
	     &val[0], &val[1], &val[2], &val[3], &val[4], &val[5]);
      for (k = 0; k < n_ps; k++)
	Mapp->M[j] += val[k]*Id[k];
    }
    break;
  default:
    std::cout << "rdmfile: unknown type" << std::endl;
    exit_(1);
    break;
  }
}


void get_elemf(const int i, const int kind)
{
  if (i == 0)
    lat.elems[i]->S = 0e0;
  else
    lat.elems[i]->S = lat.elems[i-1]->S + lat.elems[i]->PL;

  if (i > 0) {
    lat.elemf[lat.elems[i]->Fnum-1].ElemF = elem_alloc(kind);

    if (lat.elems[i]->Knum == 1)
      *lat.elemf[lat.elems[i]->Fnum-1].ElemF = *lat.elems[i];

    lat.elemf[lat.elems[i]->Fnum-1].KidList[lat.elems[i]->Knum-1] = i;
    lat.elemf[lat.elems[i]->Fnum-1].nKid =
      max(lat.elems[i]->Knum, lat.elemf[lat.elems[i]->Fnum-1].nKid);

    globval.Elem_nFam = max((long)lat.elems[i]->Fnum, globval.Elem_nFam);
  }
}


void rdmfile(const char *mfile_dat)
{
  char          line[line_max];
  long int      i;
  int           kind;
  std::ifstream inf;

  std::cout << std::endl;
  std::cout << "reading machine file: " << mfile_dat << std::endl;

  file_rd(inf, mfile_dat);

  while (inf.getline(line, line_max)) {
    get_elem(inf, line, i, kind);
    get_elemf(i, kind);
  }
  
  globval.Cell_nLoc = i;
 
  globval.dPcommon = 1e-8; globval.CODeps = 1e-14; globval.CODimax = 40;

  SI_init();

  printf("\nrdmfile: read %ld elements, C = %7.5f\n",
	 globval.Cell_nLoc, lat.elems[globval.Cell_nLoc]->S);

  inf.close();
}
