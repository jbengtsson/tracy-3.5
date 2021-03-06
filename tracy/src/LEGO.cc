
MarkerType::MarkerType(void) : CellType()
{
  this->Elem.Kind = PartsKind(marker);

  this->dR[X_] = 1e0; this->dR[Y_] = 0e0;
  this->dS[X_] = 0e0; this->dS[Y_] = 0e0;
}


DriftType::DriftType(void) : CellType()
{
  this->Elem.Kind = PartsKind(drift);

  this->dR[X_] = 1e0; this->dR[Y_] = 0e0;
  this->dS[X_] = 0e0; this->dS[Y_] = 0e0;
}


MpoleType::MpoleType(void) : CellType()
{
  int j;

  this->Elem.Kind = PartsKind(Mpole);

  this->method = Meth_Fourth; this->N = 0;
  for (j = 0; j <= 1; j++) {
    this->dSsys[j] = 0e0; this->dSrms[j] = 0e0; this->dSrnd[j] = 0e0;
  }
  this->dRpar = 0e0; this->dRsys = 0e0; this->dRrms = 0e0; this->dRrnd = 0e0;
  for (j = -HOMmax; j <= HOMmax; j++) {
    this->B[j+HOMmax]    = 0e0; this->Bpar[j+HOMmax] = 0e0;
    this->Bsys[j+HOMmax] = 0e0; this->Brms[j+HOMmax] = 0e0;
    this->Brnd[j+HOMmax] = 0e0;
  }
  this->order = 0; this->n_design = 0;
  this->irho = 0e0; this->Tx1 = 0e0; this->Tx2 = 0e0; this->gap = 0e0;

  this->c0 = 0e0; this->c1 = 0e0; this->s1 = 0e0;

  this->dR[X_] = cos(dtor(this->dRpar));
  this->dR[Y_] = sin(dtor(this->dRpar));
  this->dS[X_] = 0e0; this->dS[Y_] = 0e0;

  if (this->L != 0e0 || this->irho != 0e0) {
    this->thick = thicktype(thick_);
    this->c0 = sin(this->L*this->irho/2e0);
    this->c1 = this->dR[X_]*this->c0;
    this->s1 = this->dR[Y_]*this->c0;
  } else
    this->thick = thicktype(thin_);
}


WigglerType::WigglerType(void) : CellType()
{
  int j;

  this->Elem.Kind = PartsKind(Wigl);

  this->method = Meth_Linear; this->N = 0;
  for (j = 0; j <= 1; j++) {
    this->dSsys[j] = 0e0; this->dSrnd[j] = 0e0;
  }
  this->dRpar = 0e0; this->dRsys = 0e0; this->dRrnd = 0e0;
  this->n_harm = 0; this->lambda = 0e0;
  for (j = 0; j < n_harm_max; j++) {
    this->BoBrhoV[j] = 0e0; this->BoBrhoH[j] = 0e0;
    this->kxV[j] = 0e0; this->kxH[j] = 0e0;
    this->phi[j] = 0e0;
  }
  for (j = 0; j <= HOMmax; j++)
    this->BW[j+HOMmax] = 0e0;
  this->order = 0;

  this->dR[X_] = cos(dtor(this->dRpar));
  this->dR[Y_] = sin(dtor(this->dRpar));
  this->dS[X_] = 0e0; this->dS[Y_] = 0e0;
}


InsertionType::InsertionType(void) : CellType()
{
  int i = 0, j = 0;

  this->Elem.Kind = PartsKind(Insertion);

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
  this->dRpar = 0e0; this->dRsys = 0e0; this->dRrnd = 0e0;
//  for (j = 0; j <= HOMmax; j++)
//    this->BW[j+HOMmax] = 0e0;
  this->order = 0;

  this->dR[X_] = cos(dtor(this->dRpar));
  this->dR[Y_] = sin(dtor(this->dRpar));
  this->dS[X_] = 0e0; this->dS[Y_] = 0e0;
}


FieldMapType::FieldMapType(void) : CellType()
{

  this->Elem.Kind = PartsKind(FieldMap);

  this->n_step = 0; this->n[X_] = 0; this->n[Y_] = 0; this->n[Z_] = 0;
  this->scl = 1e0;
  this->phi = 0e0; this->Ld = 0e0; this->L1 = 0e0; this->cut = 0;
  this->x0 = 0e0;

  this->dR[X_] = 1e0; this->dR[Y_] = 0e0;
  this->dS[X_] = 0e0; this->dS[Y_] = 0e0;
}


CavityType::CavityType(void) : CellType()
{

  this->Elem.Kind = PartsKind(Cavity);

  this->volt = 0e0; this->freq = 0e0; this->phi = 0e0; this->h = 0;
  this->entry_focus = false; this->exit_focus = false;

  this->dR[X_] = 1e0; this->dR[Y_] = 0e0;
  this->dS[X_] = 0e0; this->dS[Y_] = 0e0;
}


SpreaderType::SpreaderType(void) : CellType()
{
  int k;

  this->Elem.Kind = PartsKind(Spreader);

  for (k = 0; k < Spreader_max; k++)
    this->Cell_ptrs[k] = NULL;

  this->dR[X_] = 1e0; this->dR[Y_] = 0e0;
  this->dS[X_] = 0e0; this->dS[Y_] = 0e0;
}


RecombinerType::RecombinerType(void) : CellType()
{
  this->Elem.Kind = PartsKind(Recombiner);

  this->dR[X_] = 1e0; this->dR[Y_] = 0e0;
  this->dS[X_] = 0e0; this->dS[Y_] = 0e0;
}


SolenoidType::SolenoidType(void) : CellType()
{
  int j;

  this->Elem.Kind = PartsKind(Solenoid);

  this->N = 0;
  for (j = 0; j <= 1; j++) {
    this->dSsys[j] = 0e0; this->dSrms[j] = 0e0; this->dSrnd[j] = 0e0;
  }
  this->dRpar = 0e0; this->dRsys = 0e0; this->dRrnd = 0e0;

  this->dR[X_] = 1e0; this->dR[Y_] = 0e0;
  this->dS[X_] = 0e0; this->dS[Y_] = 0e0;
}


void LatticeType::Fam_Init(const int Fnum)
{
  int         i;
  ElemFamType *elemfamp;
  CellType    *cellp;

  elemfamp = &this->ElemFam[Fnum-1];
  for (i = 1; i <= elemfamp->nKid; i++) {
    cellp = this->Cell[elemfamp->KidList[i-1]];

    strncpy(cellp->Name, elemfamp->CellF->Name, NameLength);
    cellp->L = elemfamp->CellF->L;
    cellp->Elem.Kind = elemfamp->CellF->Elem.Kind;
    cellp->Elem.Reverse = elemfamp->CellF->Elem.Reverse;
  }
}


long LatticeType::Elem_Index(const std::string &name)
{
  long        i, j;
  std::string name1;

  name1 = name;
  j = (signed)name.length();
  for (i = 0; i < j; i++)
    name1[i] = tolower(name1[i]);
  for (i = j; i < SymbolLength; i++)
    name1 += ' ';

  if (trace) {
    std::cout << std::endl;
    std::cout << "Elem_Index: " << name << " (";
    for (i = 0; i < (signed)name1.length(); i++)
      std::cout << std::setw(4) << (int)name1[i];
    std::cout << std::setw(4) << (int)name1[name1.length()] << " )"
	      << std::endl;
    std::cout << std::endl;
  }

  if (this->param.Elem_nFam > Elem_nFamMax) {
    printf("Elem_Index: Elem_nFamMax exceeded: %ld(%d)\n",
           this->param.Elem_nFam, Elem_nFamMax);
    exit_(1);
  }

  i = 1;
  while (i <= this->param.Elem_nFam) {
    if (trace) {
      std::cout << std::setw(2) << (name1 == this->ElemFam[i-1].CellF->Name)
	   << " " << name1 << " " << this->ElemFam[i-1].CellF->Name << " (";
      for (j = 0; j < SymbolLength; j++)
	std::cout << std::setw(4) << (int)this->ElemFam[i-1].CellF->Name[j];
      std::cout  << " )" << std::endl;
    }

    if (name1 == this->ElemFam[i-1].CellF->Name) break;

    i++;
  }

  if (name1 != this->ElemFam[i-1].CellF->Name) {
    std::cout << "ElemIndex: undefined element " << name << std::endl;
    exit_(1);
  }

  return i;
}


std::ostream& MarkerType::Show(std::ostream &str) const
{
  str << "      Marker |" << this->Name << "|" << "\n";
  return str;
}


std::ostream& DriftType::Show(std::ostream &str) const
{
  str << std::fixed << std::setprecision(3)
      << "      Drift  |" << this->Name << "|"
      << " L = " << setw(6) << this->L << "\n";
  return str;
}


std::ostream& MpoleType::Show(std::ostream &str) const
{
  str << std::fixed << std::setprecision(3)
      << "      Mpole  |" << this->Name << "|"
      << " L = " << setw(6) << this->L
      << ", Method = " << this->method << ", N = " << this->N << "\n";
  return str;
}


std::ostream& CavityType::Show(std::ostream &str) const
{
  str << "      Cavity |" << this->Name << "|" << "\n";
  return str;
}


std::ostream& CellType::Show(std::ostream &str) const
{
  str << "      CellType\n";
  return str;
};


std::ostream& LatticeType::Show_ElemFam(std::ostream &str) const
{
  int       j, k;
  MpoleType *M;

  const int n_prt = 10;

  str << "\nElemFam:\n";
  for (j = 0; j < this->param.Elem_nFam; j++) {
    str << std::fixed << std::setprecision(3)
	<< "\n  " << std::setw(3) << j+1
	<< " " << std::setw(3) << std::setw(3) << this->ElemFam[j].nKid
	<< " " << std::setw(5) << this->ElemFam[j].CellF->L
	<< " " << std::setw(2) << this->ElemFam[j].CellF->Elem.Kind;
    if (this->ElemFam[j].CellF->Elem.Kind == PartsKind(Mpole)) {
      M = static_cast<MpoleType*>(this->ElemFam[j].CellF);
      str << " " << std::setw(2) << M->n_design;
    } else 
      str << "   ";

    str << " |" << this->ElemFam[j].CellF->Name << "|";
    for (k = 1; k <= this->ElemFam[j].nKid; k++) {
      str << " " << std::setw(4) << this->ElemFam[j].KidList[k-1];
      if (k % n_prt == 0) str << "\n                                       ";
    }
    str << "\n";
    // this->ElemFam[j].CellF->Show(str);
  }

  return str;
}


std::ostream& LatticeType::Show(std::ostream &str) const
{
  int       k;
  MpoleType *M;

  str << "\nLattice:\n";
  for (k = 0; k <= this->param.Cell_nLoc; k++) {
    str << std::fixed << std::setprecision(3)
	<< "  " << std::setw(4) << k
	<< " " << std::setw(3) << this->Cell[k]->Fnum
	<< " " << std::setw(3) << this->Cell[k]->Knum
	<< " " << std::setw(2) << this->Cell[k]->Elem.Kind
	<< " " << std::setw(1) << this->Cell[k]->Elem.Reverse;
    if (this->Cell[k]->Elem.Kind == PartsKind(Mpole)) {
      M = static_cast<MpoleType*>(this->Cell[k]);
      str << std::fixed << std::setprecision(3)
	  << " " << std::setw(2) << M->n_design
	  << " " << std::setw(2) << M->order
	  << " " << std::setw(7) << M->Bpar[Quad+HOMmax]
	  << " " << std::setw(7) << M->B[Quad+HOMmax];
    } else 
      str << "   ";
    str << " |" << this->Cell[k]->Name << "|\n"; 
  }

  return str;
}
