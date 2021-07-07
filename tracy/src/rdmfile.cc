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


void get_elem(std::ifstream &inf, std::vector<ElemType*> &elems,
	      ConfigType &conf, char *line, long int &i, int &kind)
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

  elems[i] = elem_alloc(kind);
 
  memcpy(elems[i]->PName, name, sizeof(partsName));
  elems[i]->Fnum = Fnum; elems[i]->Knum = Knum;

  elems[i]->dS[X_] = 0e0; elems[i]->dS[Y_] = 0e0;
  elems[i]->dT[X_] = 1e0; elems[i]->dT[Y_] = 0e0;

  inf.getline(line, line_max);
  if (prt) printf("%s\n", line);
  sscanf(line, "%lf %lf %lf %lf",
	 &elems[i]->maxampl[X_][0], &elems[i]->maxampl[X_][1],
	 &elems[i]->maxampl[Y_][0], &elems[i]->maxampl[Y_][1]);

  elems[i]->PL = 0e0;
  elems[i]->Reverse = (reverse == 1);

  switch (elems[i]->Pkind) {
  case undef:
    std::cout << "rdmfile: unknown type " << i << std::endl;
    exit_(1);
    break;
  case marker:
    break;
  case drift:
    inf.getline(line, line_max);
    if (prt) printf("%s\n", line);
    sscanf(line, "%lf", &elems[i]->PL);
    break;
  case Cavity:
    C = dynamic_cast<CavityType*>(elems[i]);
    inf.getline(line, line_max);
    if (prt) printf("%s\n", line);
    sscanf(line, "%lf %lf %d %lf %lf",
	   &C->Pvolt, &C->Pfreq, &C->Ph, &conf.Energy, &C->phi);
    conf.Energy *= 1e-9;
    C->Pvolt *= conf.Energy*1e9;
    C->Pfreq *= c0/(2.0*M_PI);
    break;
  case Mpole:
    M = dynamic_cast<MpoleType*>(elems[i]);
    M->Pmethod = method; M->PN = n;

    if (M->Pthick == thick) {
      inf.getline(line, line_max);
      if (prt) printf("%s\n", line);
      sscanf(line, "%lf %lf %lf %lf",
	     &elems[i]->dS[X_], &elems[i]->dS[Y_],
	     &M->PdTpar, &dTerror);
      elems[i]->dT[X_] = cos(degtorad(dTerror+M->PdTpar));
      elems[i]->dT[Y_] = sin(degtorad(dTerror+M->PdTpar));
      M->PdTrms = dTerror - M->PdTpar;
      M->PdTrnd = 1e0;

      inf.getline(line, line_max);
      if (prt) printf("%s\n", line);
      sscanf(line, "%lf %lf %lf %lf %lf",
	     &elems[i]->PL, &M->Pirho, &M->PTx1, &M->PTx2, &M->Pgap);
      if (M->Pirho != 0e0) M->Porder = 1;
    } else {
      inf.getline(line, line_max);
      if (prt) printf("%s\n", line);
      sscanf(line, "%lf %lf %lf",
	     &elems[i]->dS[X_], &elems[i]->dS[Y_], &dTerror); 
      elems[i]->dT[X_] = cos(degtorad(dTerror));
      elems[i]->dT[Y_] = sin(degtorad(dTerror));
      M->PdTrms = dTerror; M->PdTrnd = 1e0;
    }

    M->Pc0 = sin(elems[i]->PL*M->Pirho/2.0);
    M->Pc1 = cos(degtorad(M->PdTpar))*M->Pc0;
    M->Ps1 = sin(degtorad(M->PdTpar))*M->Pc0;

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

    if (conf.mat_meth && (M->Pthick == thick))
      M->M_transp = get_transp_mat(elems[i], 0e0);
    break;
  case Wigl:
    W = dynamic_cast<WigglerType*>(elems[i]);
    W->Pmethod = method; W->PN = n;

    inf.getline(line, line_max);
    if (prt) printf("%s\n", line);
    sscanf(line, "%lf %lf", &elems[i]->PL, &W->Lambda);

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
    ID = dynamic_cast<InsertionType*>(elems[i]);
    ID->Pmethod = method; ID->PN = n;

    inf.getline(line, line_max);
    if (prt) printf("%s\n", line);
    sscanf(line, "%lf %d %s", &ID->scaling, &n, file_name);

    if (n == 1) {
      ID->firstorder = true;
      ID->secondorder = false;

      strcpy(ID->fname1, file_name);
      Read_IDfile(ID->fname1, conf, ID);
    } else if (n == 2) {
      ID->firstorder = false;
      ID->secondorder = true;

      strcpy(ID->fname2, file_name);
      Read_IDfile(ID->fname2, conf, ID);
    } else {
      std::cout << "rdmfile: undef order " << n << std::endl;
      exit_(1);
    }

    if (ID->Pmethod == 1)
      ID->linear = true;
    else
      ID->linear = false;

    if (!ID->linear) {
      // ID->tx = dmatrix(1, ID->nz,
      // 		       1, ID->nx);
      // ID->tz = dmatrix(1, ID->nz,
      // 		       1, ID->nx);
      // ID->tab1 = (double *)malloc((ID->nx)*sizeof(double));
      // ID->tab2 = (double *)malloc((ID->nz)*sizeof(double));
      // ID->f2x = dmatrix(1, ID->nz, 1, ID->nx);
      // ID->f2z = dmatrix(1, ID->nz, 1, ID->nx);
      // Matrices4Spline(ID);
    }

    /*      free_matrix(tx, 1, nz, 1, nx); free_matrix(tz, 1, nz, 1, nx);
	    free(tab1); free(tab2);
	    free_matrix(f2x, 1, nz, 1, nx); free_matrix(f2z, 1, nz, 1, nx); */
    break;
  case FieldMap:
    break;
  case Map:
    Mapp = dynamic_cast<MapType*>(elems[i]);
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


void get_elemf(std::vector<ElemFamType> &elemf, std::vector<ElemType*> &elems,
	       ConfigType &conf, const int i, const int kind)
{
  if (i == 0)
    elems[i]->S = 0e0;
  else
    elems[i]->S = elems[i-1]->S + elems[i]->PL;

  if (i > 0) {
    elemf[elems[i]->Fnum-1].ElemF = elem_alloc(kind);

    if (elems[i]->Knum == 1)
      *elemf[elems[i]->Fnum-1].ElemF = *elems[i];

    elemf[elems[i]->Fnum-1].KidList[elems[i]->Knum-1] = i;
    elemf[elems[i]->Fnum-1].nKid =
      max(elems[i]->Knum, elemf[elems[i]->Fnum-1].nKid);

    conf.Elem_nFam = max((long)elems[i]->Fnum, conf.Elem_nFam);
  }
}


void LatticeType::rdmfile(const char *mfile_dat)
{
  char          line[line_max];
  long int      i;
  int           kind;
  std::ifstream inf;

  std::cout << std::endl;
  std::cout << "reading machine file: " << mfile_dat << std::endl;

  file_rd(inf, mfile_dat);

  while (inf.getline(line, line_max)) {
    get_elem(inf, elems, conf, line, i, kind);
    get_elemf(elemf, elems, conf, i, kind);
  }
  
  conf.Cell_nLoc = i;
 
  conf.dPcommon = 1e-8; conf.CODeps = 1e-14; conf.CODimax = 40;

  SI_init();

  printf("\nrdmfile: read %ld elements, C = %7.5f\n",
	 conf.Cell_nLoc, elems[conf.Cell_nLoc]->S);

  inf.close();
}
