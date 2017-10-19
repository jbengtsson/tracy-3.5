
MarkerType::MarkerType(void)
{
  this->Kind = PartsKind(marker);
  this->dT[0] = 1e0; this->dT[1] = 0e0;
  this->dS[0] = 0e0; this->dS[1] = 0e0;
}


void MarkerType::Init(int Fnum)
{
  int         i;
  ElemFamType *elemfamp;
  CellType    *cellp;

  elemfamp = &Lattice.ElemFam[Fnum-1];
  for (i = 0; i < elemfamp->nKid; i++) {
    cellp = &Lattice.Cell[elemfamp->KidList[i]];

    strncpy(cellp->Name, elemfamp->Name, NameLength);
    cellp->Kind = elemfamp->Kind; cellp->Reverse = elemfamp->Reverse;
  }
} 

DriftType::DriftType(void)
{
  this->Kind = PartsKind(drift);
  this->dT[0] = 1e0; this->dT[1] = 0e0;
  this->dS[0] = 0e0; this->dS[1] = 0e0;
}


void DriftType::Init(int Fnum)
{
  int         i;
  ElemFamType *elemfamp;
  CellType    *cellp;

  elemfamp = &Lattice.ElemFam[Fnum-1];
  for (i = 1; i <= elemfamp->nKid; i++) {
    cellp = &Lattice.Cell[elemfamp->KidList[i-1]];

    strncpy(cellp->Name, elemfamp->Name, NameLength);
    cellp->L = elemfamp->L; cellp->Kind = elemfamp->Kind;
    cellp->Reverse = elemfamp->Reverse;
  }
}


MpoleType::MpoleType(void)
{
  int j;

  this->method = Meth_Fourth; this->N = 0;
  for (j = 0; j <= 1; j++) {
    this->dSsys[j] = 0e0; this->dSrms[j] = 0e0; this->dSrnd[j] = 0e0;
  }
  this->dTpar = 0e0; this->dTsys = 0e0; this->dTrms = 0e0; this->dTrnd = 0e0;
  for (j = -HOMmax; j <= HOMmax; j++) {
    this->B[j+HOMmax]    = 0e0; this->Bpar[j+HOMmax] = 0e0;
    this->Bsys[j+HOMmax] = 0e0; this->Brms[j+HOMmax] = 0e0;
    this->Brnd[j+HOMmax] = 0e0;
  }
  this->order = 0; this->n_design = 0;
  this->irho = 0e0; this->Tx1 = 0e0; this->Tx2 = 0e0; this->gap = 0e0;

  this->c0 = 0e0; this->c1 = 0e0; this->s1 = 0e0;

  this->Kind = PartsKind(Mpole);

  this->dT[0] = cos(dtor(this->dTpar));
  this->dT[1] = sin(dtor(this->dTpar));
  this->dS[0] = 0e0; this->dS[1] = 0e0;

  if (this->L != 0e0 || this->irho != 0e0) {
    this->thick = pthicktype(thick);
    this->c0 = sin(this->L*this->irho/2e0);
    this->c1 = this->dT[0]*this->c0;
    this->s1 = this->dT[1]*this->c0;
  } else
    this->thick = pthicktype(thin);
}


void MpoleType::Init(int Fnum)
{
  int         i;
  double      phi;
  ElemFamType *elemfamp;
  CellType    *cellp;
  MpoleType   *M;

  elemfamp = &Lattice.ElemFam[Fnum-1];
  M = static_cast<MpoleType*>(&elemfamp->CellF);
  memcpy(M->B, M->Bpar, sizeof(mpolArray));
  M->order = Updateorder(elemfamp->CellF);
  for (i = 1; i <= elemfamp->nKid; i++) {
    cellp = &Lattice.Cell[elemfamp->KidList[i-1]];

    strncpy(cellp->Name, elemfamp->Name, NameLength);
    cellp->L = elemfamp->L; cellp->Kind = elemfamp->Kind;
    cellp->Reverse = elemfamp->Reverse;
    if (reverse_elem && (cellp->Reverse == true)) {
      // Swap entrance and exit angles.
      printf("Swapping entrance and exit angles for %8s %2d\n",
	     cellp->Name, i);
      phi = M->Tx1; M->Tx1 = M->Tx2; M->Tx2 = phi; 
    }
  }
}


CavityType::CavityType(void)
{

  this->volt = 0e0; this->freq = 0e0; this->phi = 0e0; this->h = 0;
  this->entry_focus = false; this->exit_focus = false;

  this->dT[0] = 1e0; this->dT[1] = 0e0;
  this->dS[0] = 0e0; this->dS[1] = 0e0;
}


void CavityType::Init(int Fnum)
{
  int         i;
  ElemFamType *elemfamp;
  CellType    *cellp;

  elemfamp = &Lattice.ElemFam[Fnum-1];
  for (i = 0; i < elemfamp->nKid; i++) {
    cellp = &Lattice.Cell[elemfamp->KidList[i]];

    strncpy(cellp->Name, elemfamp->Name, NameLength);
    cellp->L = elemfamp->L; cellp->Kind = elemfamp->Kind;
    cellp->Reverse = elemfamp->Reverse;
  }
}


WigglerType::WigglerType(void)
{
  int j;

  this->method = Meth_Linear; this->N = 0;
  for (j = 0; j <= 1; j++) {
    this->dSsys[j] = 0e0; this->dSrnd[j] = 0e0;
  }
  this->dTpar = 0e0; this->dTsys = 0e0; this->dTrnd = 0e0;
  this->n_harm = 0; this->lambda = 0e0;
  for (j = 0; j < n_harm_max; j++) {
    this->BoBrhoV[j] = 0e0; this->BoBrhoH[j] = 0e0;
    this->kxV[j] = 0e0; this->kxH[j] = 0e0;
    this->phi[j] = 0e0;
  }
  for (j = 0; j <= HOMmax; j++)
    this->BW[j+HOMmax] = 0e0;
  this->order = 0;

  this->dT[0] = cos(dtor(this->dTpar));
  this->dT[1] = sin(dtor(this->dTpar));
  this->dS[0] = 0e0; this->dS[1] = 0e0;
}


void WigglerType::Init(int Fnum)
{
  int         i;
  ElemFamType *elemfamp;
  CellType    *cellp;
  WigglerType *W;

  elemfamp = &Lattice.ElemFam[Fnum-1];
  W = static_cast<WigglerType*>(&elemfamp->CellF);
  /* ElemF.M^.B := ElemF.M^.Bpar; */
  W->order = Quad;
  for (i = 1; i <= elemfamp->nKid; i++) {
    cellp = &Lattice.Cell[elemfamp->KidList[i-1]];

    strncpy(cellp->Name, elemfamp->Name, NameLength);
    cellp->L = elemfamp->L; cellp->Kind = elemfamp->Kind;
    cellp->Reverse = elemfamp->Reverse;
  }
}


FieldMapType::FieldMapType(void)
{

  this->n_step = 0; this->n[X_] = 0; this->n[Y_] = 0; this->n[Z_] = 0;
  this->scl = 1e0;
  this->phi = 0e0; this->Ld = 0e0; this->L1 = 0e0; this->cut = 0;
  this->x0 = 0e0;

  this->dT[0] = 1e0; this->dT[1] = 0e0;
  this->dS[X_] = 0e0; this->dS[Y_] = 0e0;
}


void FieldMapType::Init(int Fnum)
{
  int         i;
  ElemFamType *elemfamp;
  CellType    *cellp;

  elemfamp = &Lattice.ElemFam[Fnum-1];
  for (i = 1; i <= elemfamp->nKid; i++) {
    cellp = &Lattice.Cell[elemfamp->KidList[i-1]];

    strncpy(cellp->Name, elemfamp->Name, NameLength);
    cellp->L = elemfamp->L; cellp->Kind = elemfamp->Kind;
    cellp->Reverse = elemfamp->Reverse;
  }
}


InsertionType::InsertionType(void)
{
  int i = 0, j = 0;

  this->method = Meth_Linear; this->N = 0;
  this->nx = 0; this->nz = 0;

  /* Initialisation thetax and thetaz to 0*/

  // first order kick map
  if (this->firstorder){
    for (i = 0; i < IDZMAX; i++){
      for (j = 0; j < IDXMAX; j++) {
	this->thetax1[i][j] = 0e0; this->thetaz1[i][j] = 0e0;
	this->B2[i][j] = 0e0;
      }
    }
  }

  // second order kick map
  if (this->secondorder) {
    for (i = 0; i < IDZMAX; i++) {
      for (j = 0; j < IDXMAX; j++) {
          this->thetax[i][j] = 0e0; this->thetaz[i][j] = 0e0;
	  this->B2[i][j] = 0e0;
      }
    }
  }

  // stuffs for interpolation
  for (j = 0; j < IDXMAX; j++)
    this->tabx[j] = 0e0;

  for (j = 0; j < IDZMAX; j++)
    this->tabz[j] = 0e0;

  // filenames
  strcpy(this->fname1,""); strcpy(this->fname2,"");

//  this->kx = 0e0;
  for (j = 0; j <= 1; j++) {
    this->dSsys[j] = 0e0; this->dSrnd[j] = 0e0;
  }
  this->dTpar = 0e0; this->dTsys = 0e0; this->dTrnd = 0e0;
//  for (j = 0; j <= HOMmax; j++)
//    this->BW[j+HOMmax] = 0e0;
  this->order = 0;

  this->dT[0] = cos(dtor(this->dTpar));
  this->dT[1] = sin(dtor(this->dTpar));
  this->dS[0] = 0e0; this->dS[1] = 0e0;
}


void InsertionType::Init(int Fnum)
{
  int         i;
  ElemFamType *elemfamp;
  CellType    *cellp;

  elemfamp = &Lattice.ElemFam[Fnum-1];
//  elemfamp->CellF->order = order;
//  x = elemfamp->CellF->BW[Quad + HOMmax];
  for (i = 1; i <= elemfamp->nKid; i++) {
    cellp = &Lattice.Cell[elemfamp->KidList[i-1]];

    strncpy(cellp->Name, elemfamp->Name, NameLength);
    cellp->L = elemfamp->L; cellp->Kind = elemfamp->Kind;
    cellp->Reverse = elemfamp->Reverse;
  }
}


SpreaderType::SpreaderType(void)
{
  int k;

  for (k = 0; k < Spreader_max; k++)
    this->Cell_ptrs[k] = NULL;

  this->dT[0] = 1e0; this->dT[1] = 0e0;
  this->dS[0] = 0e0; this->dS[1] = 0e0;
}


void SpreaderType::Init(int Fnum)
{
  int         i;
  ElemFamType *elemfamp;
  CellType    *cellp;

  elemfamp = &Lattice.ElemFam[Fnum-1];
  for (i = 1; i <= elemfamp->nKid; i++) {
     cellp = &Lattice.Cell[elemfamp->KidList[i-1]];

    strncpy(cellp->Name, elemfamp->Name, NameLength);
    cellp->Kind = elemfamp->Kind;
    cellp->Reverse = elemfamp->Reverse;
  }
}


RecombinerType::RecombinerType(void)
{
  this->dT[0] = 1e0; this->dT[1] = 0e0;
  this->dS[0] = 0e0; this->dS[1] = 0e0;
}


void RecombinerType::Init(int Fnum)
{
  int         i;
  ElemFamType *elemfamp;
  CellType    *cellp;

  elemfamp = &Lattice.ElemFam[Fnum-1];
  for (i = 1; i <= elemfamp->nKid; i++) {
    cellp = &Lattice.Cell[elemfamp->KidList[i-1]];

    strncpy(cellp->Name, elemfamp->Name, NameLength);
    cellp->Kind = elemfamp->Kind;
    cellp->Reverse = elemfamp->Reverse;
  }
}


SolenoidType::SolenoidType(void)
{
  int j;

  this->N = 0;
  for (j = 0; j <= 1; j++) {
    this->dSsys[j] = 0e0; this->dSrms[j] = 0e0; this->dSrnd[j] = 0e0;
  }
  this->dTpar = 0e0; this->dTsys = 0e0; this->dTrnd = 0e0;

  this->dT[0] = 1e0; this->dT[1] = 0e0;
  this->dS[0] = 0e0; this->dS[1] = 0e0;
}


void SolenoidType::Init(int Fnum)
{
  int         i;
  ElemFamType *elemfamp;
  CellType    *cellp;

  elemfamp = &Lattice.ElemFam[Fnum-1];
  for (i = 1; i <= elemfamp->nKid; i++) {
    cellp = &Lattice.Cell[elemfamp->KidList[i-1]];

    strncpy(cellp->Name, elemfamp->Name, NameLength);
    cellp->L = elemfamp->L;
    cellp->Kind = elemfamp->Kind;
    cellp->Reverse = elemfamp->Reverse;
  }
}


void Fam_Init(int Fnum)
{
  int         i;
  ElemFamType *elemfamp;
  CellType    *cellp;

  elemfamp = &Lattice.ElemFam[Fnum-1];
  for (i = 1; i <= elemfamp->nKid; i++) {
    cellp = &Lattice.Cell[elemfamp->KidList[i-1]];

    strncpy(cellp->Name, elemfamp->Name, NameLength);
    cellp->L = elemfamp->L;
    cellp->Kind = elemfamp->Kind;
    cellp->Reverse = elemfamp->Reverse;
  }
}
