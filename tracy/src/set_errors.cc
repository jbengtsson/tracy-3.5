

const static bool normal = true; // Normal or rectangular distribution.

const static long int
  k_ = 19,
  c_ = 656329L,
  m_ = 100000001;

const static int
  maxiter = 100;

static long int
  rseed0,
  rseed;
static double
  normcut_;

void iniranf(const long i) { rseed0 = i; rseed = i; }

void newseed(void)
{
  rseed0 = (k_*rseed0+c_) % m_; rseed = (rseed0+54321) % m_;
}

// Random number generator with rectangular distribution.
double ranf(void) { rseed = (k_*rseed+c_) % m_; return (rseed/1e8); }

void setrancut(const double cut)
{
  printf("\nsetrancut: cut set to %3.1f\n", cut);
  normcut_ = cut;
}

double normranf(void)
{  
  int    i, j;
  double f, w;

  j = 0;
  do {
    j++;
    w = 0.0;
    for (i = 1; i <= 12; i++)
      w += ranf();
    f = w - 6.0;
  }
  while (fabs(f) > fabs(normcut_) && j <= maxiter);

  if (j > maxiter) {
    printf("normranf: algorithm did not converge\n");
    exit(1);
  }
  return f;
}


// Misalignments.

void CheckAlignTol(LatticeType &lat, const char *OutputFile)
  // check aligment errors of individual magnets on giders
  // the dT and roll angle are all printed out
{
  int          i, j;
  int          n_girders;
  int          gs_Fnum, ge_Fnum;
  int          gs_nKid, ge_nKid;
  int          dip_Fnum,dip_nKid;
  int          loc, loc_gs, loc_ge;
  char         *name;
  double       s;
  double       PdSsys[2], PdSrms[2], PdSrnd[2], dS[2], dT[2];
  MpoleType    *M;
  std::fstream fout;

  gs_Fnum = lat.conf.gs;   gs_nKid = lat.GetnKid(gs_Fnum);
  ge_Fnum = lat.conf.ge;   ge_nKid = lat.GetnKid(ge_Fnum);
  if (gs_nKid == ge_nKid)
    n_girders= gs_nKid;
  else {
    std::cout << " The numbers of GS and GE not same. " << std::endl;
    exit (1);
  }

  fout.open(OutputFile,std::ios::out);
  if(!fout) {
    std::cout << "error in opening the file  " << std::endl;
    exit_(0);
  }

  fout << "Girders, Quads, Sexts:  " << std::endl;
  for (i = 1; i <= n_girders; i++){
    fout << i << ":" << std::endl;
    loc_gs = lat.Elem_GetPos(gs_Fnum, i); loc_ge = lat.Elem_GetPos(ge_Fnum, i);

    loc = loc_gs;
    M = dynamic_cast<MpoleType*>(lat.elems[loc]);
    PdSsys[X_] = M->PdSsys[X_];
    PdSsys[Y_] = M->PdSsys[Y_];
    PdSrms[X_] = M->PdSrms[X_];
    PdSrms[Y_] = M->PdSrms[Y_];
    PdSrnd[X_] = M->PdSrnd[X_];
    PdSrnd[Y_] = M->PdSrnd[Y_];
    dS[X_] = lat.elems[loc]->dS[X_]; dS[Y_] = lat.elems[loc]->dS[Y_];
    dT[0] = lat.elems[loc]->dT[0]; dT[1] = lat.elems[loc]->dT[1];
    s = lat.elems[loc]->S; name = lat.elems[loc]->PName;
    fout << "  " << name << "  " << loc << "   " << s
	 << "  " <<  PdSsys[X_] << "  " <<  PdSsys[Y_]
	 << "   " << PdSrms[X_] << "  " <<  PdSrms[Y_]
	 << "   " << PdSrnd[X_] << "  " <<  PdSrnd[Y_]
         << "   " << M->PdTrms << "  "
	 << M->PdTrnd << "   " << dS[X_]     << "  " <<  dS[Y_]
	 << "   " << atan2( dT[1], dT[0] )  << std::endl;

    for (j = loc_gs+1; j < loc_ge; j++) {
      loc = j;
      M = dynamic_cast<MpoleType*>(lat.elems[loc]);
      if ((lat.elems[j]->Pkind == Mpole) &&
	  (M->n_design >= Quad || M->n_design >= Sext)) {
	PdSsys[X_] = M->PdSsys[X_];
	PdSsys[Y_] = M->PdSsys[Y_];
	PdSrms[X_] = M->PdSrms[X_];
	PdSrms[Y_] = M->PdSrms[Y_];
	PdSrnd[X_] = M->PdSrnd[X_];
	PdSrnd[Y_] = M->PdSrnd[Y_];
	dS[X_] = lat.elems[loc]->dS[X_]; dS[Y_] = lat.elems[loc]->dS[Y_];
	dT[0] = lat.elems[loc]->dT[0];   dT[1] = lat.elems[loc]->dT[1];
	s = lat.elems[loc]->S; name=lat.elems[loc]->PName;
	fout << "  " << name << "  " << loc << "   " << s
	     << "  " <<  PdSsys[X_] << "  " <<  PdSsys[Y_]
	     << "   " << PdSrms[X_] << "  " <<  PdSrms[Y_]
	     << "   " << PdSrnd[X_] << "  " <<  PdSrnd[Y_]
	     << "   " << M->PdTrms << "  "
	     << M->PdTrnd
	     << "   " << dS[X_] << "  " <<  dS[Y_]
	     << "   " << atan2( dT[1], dT[0] )  << std::endl;
      }
    }

    loc = loc_ge;
    M = dynamic_cast<MpoleType*>(lat.elems[loc]);
    PdSsys[X_] = M->PdSsys[X_];
    PdSsys[Y_] = M->PdSsys[Y_];
    PdSrms[X_] = M->PdSrms[X_];
    PdSrms[Y_] = M->PdSrms[Y_];
    PdSrnd[X_] = M->PdSrnd[X_];
    PdSrnd[Y_] = M->PdSrnd[Y_];
    dS[X_] = lat.elems[loc]->dS[X_]; dS[Y_] = lat.elems[loc]->dS[Y_];
    dT[0] = lat.elems[loc]->dT[0]; dT[1] = lat.elems[loc]->dT[1];
    s=lat.elems[loc]->S; name=lat.elems[loc]->PName;
    fout << "  " << name << "  " << loc << "   " << s
	 << "  " <<  PdSsys[X_] << "  " <<  PdSsys[Y_]
	 << "   " << PdSrms[X_] << "  " <<  PdSrms[Y_]
	 << "   " << PdSrnd[X_] << "  " <<  PdSrnd[Y_]
         << "   " << M->PdTrms
	 << "  " << M->PdTrnd
         << "   " << dS[X_]     << "  " <<  dS[Y_]
	 << "   " << atan2( dT[1], dT[0] )  << std::endl;

  }

  fout << "  " << std::endl;
  fout << "Dipoles:  " << std::endl;
  dip_Fnum = ElemIndex("B1"); dip_nKid = lat.GetnKid(dip_Fnum);
  for (i = 1; i <= dip_nKid; i++){
    loc = lat.Elem_GetPos(dip_Fnum, i);
    M = dynamic_cast<MpoleType*>(lat.elems[loc]);
    PdSsys[X_] = M->PdSsys[X_];
    PdSsys[Y_] = M->PdSsys[Y_];
    PdSrms[X_] = M->PdSrms[X_];
    PdSrms[Y_] = M->PdSrms[Y_];
    PdSrnd[X_] = M->PdSrnd[X_];
    PdSrnd[Y_] = M->PdSrnd[Y_];
    dS[X_] = lat.elems[loc]->dS[X_]; dS[Y_] = lat.elems[loc]->dS[Y_];
    dT[0] = lat.elems[loc]->dT[0]; dT[1] = lat.elems[loc]->dT[1];
    s = lat.elems[loc]->S; name = lat.elems[loc]->PName;
    fout << "  " << name << "  " << loc << "   " << s
	 << "  " <<  PdSsys[X_] << "  " <<  PdSsys[Y_]
	 << "   " << PdSrms[X_] << "  " <<  PdSrms[Y_]
	 << "   " << PdSrnd[X_] << "  " <<  PdSrnd[Y_]
	 << "   " << M->PdTrms
	 << "  " << M->PdTrnd
	 << "   " << dS[X_]     << "  " <<  dS[Y_]
	 << "   " << atan2( dT[1], dT[0] )  << std::endl;
  }

  fout.close();
}


void misalign_rms_elem(LatticeType &lat, const int Fnum, const int Knum,
		       const double dx_rms, const double dy_rms,
		       const double dr_rms, const bool new_rnd)
{
  long int  loc;
  MpoleType *mp;

  loc = lat.Elem_GetPos(Fnum, Knum);
  mp = dynamic_cast<MpoleType*>(lat.elems[loc]);

  mp->PdSrms[X_] = dx_rms; mp->PdSrms[Y_] = dy_rms; mp->PdTrms = dr_rms;
  if (new_rnd) {
    if (normal) {
      mp->PdSrnd[X_] = normranf(); mp->PdSrnd[Y_] = normranf();
      mp->PdTrnd = normranf();
    } else {
      mp->PdSrnd[X_] = ranf(); mp->PdSrnd[Y_] = ranf();
      mp->PdTrnd = ranf();
    }
  }

  lat.SetdS(Fnum, Knum); lat.SetdT(Fnum, Knum);
}

void misalign_sys_elem(LatticeType &lat, const int Fnum, const int Knum,
		       const double dx_sys, const double dy_sys,
		       const double dr_sys)
{
  long int  loc;
  MpoleType *mp;

  loc = lat.Elem_GetPos(Fnum, Knum);
  mp  = dynamic_cast<MpoleType*>(lat.elems[loc]);

  mp->PdSsys[X_] = dx_sys; mp->PdSsys[Y_] = dy_sys; mp->PdTsys = dr_sys;

  lat.SetdS(Fnum, Knum); lat.SetdT(Fnum, Knum);
}

void misalign_rms_fam(LatticeType &lat, const int Fnum,
		      const double dx_rms, const double dy_rms,
		      const double dr_rms, const bool new_rnd)
{
  int  i;

  for (i = 1; i <= lat.GetnKid(Fnum); i++)
    misalign_rms_elem(lat, Fnum, i, dx_rms, dy_rms, dr_rms, new_rnd);
}

void misalign_sys_fam(LatticeType &lat, const int Fnum,
		      const double dx_sys, const double dy_sys,
		      const double dr_sys)
{
  int  i;

  for (i = 1; i <= lat.GetnKid(Fnum); i++)
    misalign_sys_elem(lat, Fnum, i, dx_sys, dy_sys, dr_sys);
}

void misalign_rms_type(LatticeType &lat, const int type,
		       const double dx_rms, const double dy_rms,
		       const double dr_rms, const bool new_rnd)
{
  long int  k;
  MpoleType *M;

  if ((type >= All) && (type <= HOMmax)) {
    for (k = 1; k <= lat.conf.Cell_nLoc; k++) {
      M = dynamic_cast<MpoleType*>(lat.elems[k]);
      if ((lat.elems[k]->Pkind == Mpole) &&
	  ((type == M->n_design) || ((type == All) &&
	   ((lat.elems[k]->Fnum != lat.conf.gs)
	    && (lat.elems[k]->Fnum != lat.conf.ge))))) {
	// if all: skip girders
	misalign_rms_elem(lat, lat.elems[k]->Fnum, lat.elems[k]->Knum,
			  dx_rms, dy_rms, dr_rms, new_rnd);
      }
    }
  } else {
    printf("misalign_rms_type: incorrect type %d\n", type); exit_(1);
  }
}

void misalign_sys_type(LatticeType &lat, const int type,
		       const double dx_sys, const double dy_sys,
		       const double dr_sys)
{
  long int  k;
  MpoleType *M;

  if ((type >= All) && (type <= HOMmax)) {
    for (k = 1; k <= lat.conf.Cell_nLoc; k++) {
      M = dynamic_cast<MpoleType*>(lat.elems[k]);
      if ((lat.elems[k]->Pkind == Mpole) &&
	  ((type == M->n_design) || ((type == All) &&
	   ((lat.elems[k]->Fnum != lat.conf.gs)
	    && (lat.elems[k]->Fnum != lat.conf.ge))))) {
	// if all: skip girders
	misalign_sys_elem(lat, lat.elems[k]->Fnum, lat.elems[k]->Knum,
			  dx_sys, dy_sys, dr_sys);
      }
    }
  } else {
    printf("misalign_sys_type: incorrect type %d\n", type); exit_(1);
  }
}

void misalign_rms_girders(LatticeType &lat, const int gs, const int ge,
			  const double dx_rms, const double dy_rms,
			  const double dr_rms, const bool new_rnd)
{
  int       i, k, n_girders, n_ge, n_gs;
  long int  loc_gs, loc_ge, j;
  double    s_gs, s_ge, dx_gs[2], dx_ge[2], s;
  MpoleType *Mgs, *Mge, *Mj;

  n_gs = lat.GetnKid(gs); n_ge = lat.GetnKid(ge);

  if (n_gs == n_ge)
    n_girders = n_gs;
  else {
    std::cout << "set_girders: no of GS != no of GE" << std::endl;
    exit (1);
  }

  misalign_rms_fam(lat, gs, dx_rms, dy_rms, dr_rms, new_rnd);
  misalign_rms_fam(lat, ge, dx_rms, dy_rms, dr_rms, new_rnd);

  for (i = 1; i <= n_girders; i++) {
    loc_gs = lat.Elem_GetPos(gs, i);
    loc_ge = lat.Elem_GetPos(ge, i);
    s_gs = lat.elems[loc_gs]->S;
    s_ge = lat.elems[loc_ge]->S;

    // roll for a rigid boby
    // Note, girders needs to be introduced as gs->ge pairs
    Mgs = dynamic_cast<MpoleType*>(lat.elems[loc_gs]);
    Mge = dynamic_cast<MpoleType*>(lat.elems[loc_ge]);
    Mge->PdTrnd = Mgs->PdTrnd;
    lat.SetdT(ge, i);

    for (k = 0; k <= 1; k++) {
      dx_gs[k] = lat.elems[loc_gs]->dS[k]; dx_ge[k] = lat.elems[loc_ge]->dS[k];
    }

    // move elements onto mis-aligned girder
    for (j = loc_gs+1; j < loc_ge; j++) {
      if ((lat.elems[j]->Pkind == Mpole)
	  || (lat.elems[j]->Fnum == lat.conf.bpm)) {
	Mj = dynamic_cast<MpoleType*>(lat.elems[j]);
        s = lat.elems[j]->S;
	for (k = 0; k <= 1; k++)
	  Mj->PdSsys[k] = dx_gs[k] + (dx_ge[k]-dx_gs[k])*(s-s_gs)/(s_ge-s_gs);
	Mj->PdTsys = Mgs->PdTrms*Mgs->PdTrnd;
      }
    }
  }
}


void misalign_sys_girders(LatticeType &lat, const int gs, const int ge,
			  const double dx_sys, const double dy_sys,
			  const double dr_sys)
{
  int       i, k, n_girders, n_ge, n_gs;
  long int  loc_gs, loc_ge, j;
  double    s_gs, s_ge, dx_gs[2], dx_ge[2], s;
  MpoleType *Mgs, *Mge, *Mj;

  n_gs = lat.GetnKid(gs); n_ge = lat.GetnKid(ge);

  if (n_gs == n_ge)
    n_girders = n_gs;
  else {
    std::cout << "set_girders: no of GS != no of GE" << std::endl;
    exit (1);
  }

  misalign_sys_fam(lat, gs, dx_sys, dy_sys, dr_sys);
  misalign_sys_fam(lat, ge, dx_sys, dy_sys, dr_sys);

  for (i = 1; i <= n_girders; i++) {
    loc_gs = lat.Elem_GetPos(gs, i); loc_ge = lat.Elem_GetPos(ge, i);
    s_gs = lat.elems[loc_gs]->S; s_ge = lat.elems[loc_ge]->S;

    // roll for a rigid boby
    // Note, girders needs to be introduced as gs->ge pairs
    Mgs = dynamic_cast<MpoleType*>(lat.elems[loc_gs]);
    Mge = dynamic_cast<MpoleType*>(lat.elems[loc_ge]);
    Mge->PdTrnd = Mgs->PdTrnd;
    lat.SetdT(ge, i);

    for (k = 0; k <= 1; k++) {
      dx_gs[k] = lat.elems[loc_gs]->dS[k]; dx_ge[k] = lat.elems[loc_ge]->dS[k];
    }

    // move elements onto mis-aligned girder
    for (j = loc_gs+1; j < loc_ge; j++) {
      if ((lat.elems[j]->Pkind == Mpole)
	  || (lat.elems[j]->Fnum == lat.conf.bpm)) {
	Mj = dynamic_cast<MpoleType*>(lat.elems[j]);
        s = lat.elems[j]->S;
	for (k = 0; k <= 1; k++)
	  Mj->PdSsys[k] = dx_gs[k] + (dx_ge[k]-dx_gs[k])*(s-s_gs)/(s_ge-s_gs);
	Mj->PdTsys = Mgs->PdTrms*Mgs->PdTrnd;
      }
    }
  }
}


// Apertures.

void set_aper_elem(LatticeType &lat, const int Fnum, const int Knum,
		   const double Dxmin, const double Dxmax,
		   const double Dymin, const double Dymax)
{
  int k;

    k = lat.Elem_GetPos(Fnum, Knum);
    lat.elems[k]->maxampl[X_][0] = Dxmin; lat.elems[k]->maxampl[X_][1] = Dxmax;
    lat.elems[k]->maxampl[Y_][0] = Dymin; lat.elems[k]->maxampl[Y_][1] = Dymax;
 }

void set_aper_fam(LatticeType &lat, const int Fnum,
		  const double Dxmin, const double Dxmax,
		  const double Dymin, const double Dymax)
{
  int k;

  for (k = 1; k <= lat.GetnKid(Fnum); k++)
    set_aper_elem(lat, Fnum, k, Dxmin, Dxmax, Dymin, Dymax);
}

void set_aper_type(LatticeType &lat, const int type, const double Dxmin,
		   const double Dxmax, const double Dymin, const double Dymax)
{
  long int  k;
  MpoleType *M;

  if (type >= All && type <= HOMmax) {
    for(k = 1; k <= lat.conf.Cell_nLoc; k++) {
      M = dynamic_cast<MpoleType*>(lat.elems[k]);
      if (((lat.elems[k]->Pkind == Mpole) && (M->n_design == type))
	  || (type == All))
	set_aper_elem(lat, lat.elems[k]->Fnum, lat.elems[k]->Knum, Dxmin, Dxmax,
		      Dymin, Dymax);
    }
  } else
    printf("set_aper_type: bad design type %d\n", type);
}


// Miscellaneous.

double get_L(LatticeType &lat, const int Fnum, const int Knum)
{
  return lat.elems[lat.Elem_GetPos(Fnum, Knum)]->PL;
}


void set_L(LatticeType &lat, const int Fnum, const int Knum, const double L)
{
  long int  loc;
  double    phi;
  ElemType  *elemp;
  MpoleType *M;

  loc = lat.Elem_GetPos(Fnum, Knum);
  elemp = lat.elems[loc];
  if (elemp->Pkind == Mpole) {
    M = dynamic_cast<MpoleType*>(elemp);
    if (M->Pirho != 0e0) {
      // Phi is constant.
      phi = elemp->PL*M->Pirho; M->Pirho = phi/L;
      // M->Pc0 = sin(phi/2e0);
    }
  }
  elemp->PL = L;
}


void set_L(LatticeType &lat, const int Fnum, const double L)
{
  int k;

  for (k = 1; k <= lat.GetnKid(Fnum); k++)
    set_L(lat, Fnum, k, L);
}


void set_dL(LatticeType &lat, const int Fnum, const int Knum, const double dL)
{
  lat.elems[lat.Elem_GetPos(Fnum, Knum)]->PL += dL;
}


void set_dL(LatticeType &lat, const int Fnum, const double dL)
{
  int k;

  for (k = 1; k <= lat.GetnKid(Fnum); k++)
    set_dL(lat, Fnum, k, dL);
}


// Multipoles.

void LatticeType::set_b_n(const set_mpole set_, const int Fnum, const int Knum,
			 const int n, const double b_n)
{
  //       / b_n, n > 0 
  // : set |
  //       \ a_n, n < 0

  ElemType  *elemp;
  MpoleType *M;

  if ((-HOMmax <= n) && (n <= HOMmax)) {
    elemp = elems[Elem_GetPos(Fnum, Knum)];
    M = dynamic_cast<MpoleType*>(elemp);

    switch (set_) {
    case set_mpole(b_n_):
      M->PBpar[HOMmax+n] = b_n;
      break;
    case set_mpole(db_n_):
      M->PBpar[HOMmax+n] += b_n;
      break;
    case set_mpole(b_nL_):
      M->PBpar[HOMmax+n] = (elemp->PL != 0e0)? b_n/elemp->PL : b_n;
      break;
    case set_mpole(db_nL_):
      M->PBpar[HOMmax+n] += (elemp->PL != 0e0)? b_n/elemp->PL : b_n;
      break;
    default:
      printf("\nset_b_n: unknown set_ %d\n", set_);
      exit(1);
      break;
    }

   SetPB(Fnum, Knum, n);
  } else {
    printf("set_b_n: n < 1 (%d)\n", n);
    exit(1);
  }
}

//------------------------------------------------------------------------------

void get_bn_design_elem(LatticeType &lat, const int Fnum, const int Knum,
			const int n, double &bn, double &an)
{
  ElemType  *elem;
  MpoleType *M;

  if (n < 1) {
    std::cout << "get_bn_design_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  elem = lat.elems[lat.Elem_GetPos(Fnum, Knum)];
  M = dynamic_cast<MpoleType*>(elem);

  bn = M->PBpar[HOMmax+n]; an = M->PBpar[HOMmax-n];
}


void get_bnL_design_elem(LatticeType &lat, const int Fnum, const int Knum,
			 const int n, double &bnL, double &anL)
{
  ElemType  *elem;
  MpoleType *M;

  if (n < 1) {
    std::cout << "get_bnL_design_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  elem = lat.elems[lat.Elem_GetPos(Fnum, Knum)];
  M = dynamic_cast<MpoleType*>(elem);

  bnL = M->PBpar[HOMmax+n]; anL = M->PBpar[HOMmax-n];

  if (elem->PL != 0e0) {
    bnL *= elem->PL; anL *= elem->PL;
  }
}


void set_bn_design_elem(LatticeType &lat, const int Fnum, const int Knum,
			const int n, const double bn, const double an)
{
  ElemType  *elem;
  MpoleType *M;

  if (n < 1) {
    std::cout << "set_bn_design_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  elem = lat.elems[lat.Elem_GetPos(Fnum, Knum)];
  M = dynamic_cast<MpoleType*>(elem);

  M->PBpar[HOMmax+n] = bn; M->PBpar[HOMmax-n] = an;

  lat.SetPB(Fnum, Knum, n); lat.SetPB(Fnum, Knum, -n);
}


void set_bn_design_fam(LatticeType &lat, const int Fnum, const int n,
		       const double bn, const double an)
{
  int k;

  if (n < 1) {
    std::cout << "set_bn_design_fam: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  for (k = 1; k <= lat.GetnKid(Fnum); k++)
    set_bn_design_elem(lat, Fnum, k, n, bn, an);
}


void set_dbn_design_elem(LatticeType &lat, const int Fnum, const int Knum,
			 const int n, const double dbn, const double dan)
{
  ElemType  *elem;
  MpoleType *M;

  if (n < 1) {
    std::cout << "set_dbn_design_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  elem = lat.elems[lat.Elem_GetPos(Fnum, Knum)];
  M = dynamic_cast<MpoleType*>(elem);

  M->PBpar[HOMmax+n] += dbn; M->PBpar[HOMmax-n] += dan;

  lat.SetPB(Fnum, Knum, n); lat.SetPB(Fnum, Knum, -n);
}


void set_dbn_design_fam(LatticeType &lat, const int Fnum, const int n,
			const double dbn, const double dan)
{
  int k;

  if (n < 1) {
    std::cout << "set_dbn_design_fam: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  for (k = 1; k <= lat.GetnKid(Fnum); k++)
    set_dbn_design_elem(lat, Fnum, k, n, dbn, dan);
}


void set_bnL_design_elem(LatticeType &lat, const int Fnum, const int Knum,
			 const int n, const double bnL, const double anL)
{
  ElemType  *elem;
  MpoleType *M;

  if (n < 1) {
    std::cout << "set_bnL_design_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  elem = lat.elems[lat.Elem_GetPos(Fnum, Knum)];
  M = dynamic_cast<MpoleType*>(elem);

  if (elem->PL != 0e0) {
    M->PBpar[HOMmax+n] = bnL/elem->PL;
    M->PBpar[HOMmax-n] = anL/elem->PL;
  } else {
    // thin kick
    M->PBpar[HOMmax+n] = bnL; M->PBpar[HOMmax-n] = anL;
  }

  lat.SetPB(Fnum, Knum, n); lat.SetPB(Fnum, Knum, -n);
}


void set_dbnL_design_elem(LatticeType &lat, const int Fnum, const int Knum,
			  const int n, const double dbnL, const double danL)
{
  ElemType  *elem;
  MpoleType *M;

  if (n < 1) {
    std::cout << "set_dbnL_design_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  elem = lat.elems[lat.Elem_GetPos(Fnum, Knum)];
  M = dynamic_cast<MpoleType*>(elem);

  if (elem->PL != 0e0) {
    M->PBpar[HOMmax+n] += dbnL/elem->PL;
    M->PBpar[HOMmax-n] += danL/elem->PL;
  } else {
    // thin kick
    M->PBpar[HOMmax+n] += dbnL; M->PBpar[HOMmax-n] += danL;
  }

  lat.SetPB(Fnum, Knum, n); lat.SetPB(Fnum, Knum, -n);
}


void set_dbnL_design_fam(LatticeType &lat, const int Fnum, const int n,
			 const double dbnL, const double danL)
{
  int k;

  if (n < 1) {
    std::cout << "set_dbnL_design_fam: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  for (k = 1; k <= lat.GetnKid(Fnum); k++)
    set_dbnL_design_elem(lat, Fnum, k, n, dbnL, danL);
}


void set_bnL_design_fam(LatticeType &lat, const int Fnum, const int n,
			const double bnL, const double anL)
{
  int k;

  if (n < 1) {
    std::cout << "set_bnL_design_fam: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  for (k = 1; k <= lat.GetnKid(Fnum); k++)
    set_bnL_design_elem(lat, Fnum, k, n, bnL, anL);
}

//------------------------------------------------------------------------------

void set_bnL_design_type(LatticeType &lat, const int type, const int n,
			 const double bnL, const double anL)
{
  long int  k;
  MpoleType *M;

  if (n < 1) {
    std::cout << "set_bnL_design_type: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  if ((type >= Dip) && (type <= HOMmax)) {
    for (k = 1; k <= lat.conf.Cell_nLoc; k++) {
      M = dynamic_cast<MpoleType*>(lat.elems[k]);
      if ((lat.elems[k]->Pkind == Mpole) && (M->n_design == type))
	set_bnL_design_elem(lat, lat.elems[k]->Fnum, lat.elems[k]->Knum, n, bnL,
			    anL);
    }
  } else
    printf("Bad type argument to set_bnL_design_type()\n");
}

//------------------------------------------------------------------------------

void set_bnL_sys_elem(LatticeType &lat, const int Fnum, const int Knum,
		      const int n, const double bnL, const double anL)
{
  ElemType  *elem;
  MpoleType *M;

  if (n < 1) {
    std::cout << "set_bnL_sys_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  elem = lat.elems[lat.Elem_GetPos(Fnum, Knum)];
  M = dynamic_cast<MpoleType*>(elem);

  if (elem->PL != 0e0) {
    M->PBsys[HOMmax+n] = bnL/elem->PL;
    M->PBsys[HOMmax-n] = anL/elem->PL;
  } else {
    // thin kick
    M->PBsys[HOMmax+n] = bnL; M->PBsys[HOMmax-n] = anL;
  }

  lat.SetPB(Fnum, Knum, n); lat.SetPB(Fnum, Knum, -n);
}


void set_bnL_sys_fam(LatticeType &lat, const int Fnum, const int n,
		     const double bnL, const double anL)
{
  int k;

  if (n < 1) {
    std::cout << "set_bnL_sys_fam: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  for (k = 1; k <= lat.GetnKid(Fnum); k++)
    set_bnL_sys_elem(lat, Fnum, k, n, bnL, anL);
}


void set_bnL_sys_type(LatticeType &lat, const int type, const int n,
		      const double bnL, const double anL)
{
  long int  k;
  MpoleType *M;

  if (n < 1) {
    std::cout << "set_bnL_sys_type: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  if (type >= Dip && type <= HOMmax) {
    for(k = 1; k <= lat.conf.Cell_nLoc; k++) {
      M = dynamic_cast<MpoleType*>(lat.elems[k]);
      if ((lat.elems[k]->Pkind == Mpole) && (M->n_design == type))
	set_bnL_sys_elem(lat, lat.elems[k]->Fnum, lat.elems[k]->Knum, n, bnL,
			 anL);
    }
  } else
    printf("Bad type argument to set_bnL_sys_type()\n");
}


void set_bnL_rms_elem(LatticeType &lat, const int Fnum, const int Knum,
		      const int n, const double bnL, const double anL,
		      const bool new_rnd)
{
  ElemType  *elem;
  MpoleType *M;

  bool prt = false;

  if (n < 1) {
    std::cout << "set_bnL_rms_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  elem = lat.elems[lat.Elem_GetPos(Fnum, Knum)];
  M = dynamic_cast<MpoleType*>(elem);

  if (elem->PL != 0e0) {
    M->PBrms[HOMmax+n] = bnL/elem->PL;
    M->PBrms[HOMmax-n] = anL/elem->PL;
  } else {
    // thin kick
    M->PBrms[HOMmax+n] = bnL; M->PBrms[HOMmax-n] = anL;
  }

  if(new_rnd){
    if (normal) {
      M->PBrnd[HOMmax+n] = normranf();
      M->PBrnd[HOMmax-n] = normranf();
    } else {
      M->PBrnd[HOMmax+n] = ranf(); M->PBrnd[HOMmax-n] = ranf();
    }
  }

  if (prt)
    printf("set_bnL_rms_elem:  Fnum = %d, Knum = %d"
	   ", bnL = %e, anL = %e %e %e\n",
	   Fnum, Knum, bnL, anL,
	   M->PBrms[HOMmax+n], M->PBrms[HOMmax-n]);

  lat.SetPB(Fnum, Knum, n); lat.SetPB(Fnum, Knum, -n);
}


void set_bnL_rms_fam(LatticeType &lat, const int Fnum, const int n,
		     const double bnL, const double anL, const bool new_rnd)
{
  int k;

  if (n < 1) {
    std::cout << "set_bnL_rms_fam: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  for (k = 1; k <= lat.GetnKid(Fnum); k++)
    set_bnL_rms_elem(lat, Fnum, k, n, bnL, anL, new_rnd);
}


void set_bnL_rms_type(LatticeType &lat, const int type, const int n,
		      const double bnL, const double anL, const bool new_rnd)
{
  long int  k;
  MpoleType *M;

  if (n < 1) {
    std::cout << "get_bnL_rms_type: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  if (type >= Dip && type <= HOMmax) {
    for(k = 1; k <= lat.conf.Cell_nLoc; k++) {
      M = dynamic_cast<MpoleType*>(lat.elems[k]);
      if ((lat.elems[k]->Pkind == Mpole) && (M->n_design == type))
	set_bnL_rms_elem(lat, lat.elems[k]->Fnum, lat.elems[k]->Knum, n, bnL,
			 anL, new_rnd);
    }
  } else
    printf("Bad type argument to set_bnL_rms_type()\n");
}


void set_bnr_sys_elem(LatticeType &lat, const int Fnum, const int Knum,
		      const int n, const double bnr, const double anr)
{
  int       nd;
  MpoleType *M;
  bool      prt = false;

  if (n < 1) {
    std::cout << "set_bnr_sys_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  M = dynamic_cast<MpoleType*>(lat.elems[lat.Elem_GetPos(Fnum, Knum)]);
  nd = M->n_design;
  // errors are relative to design values for (Dip, Quad, Sext, ...)
  M->PBsys[HOMmax+n] = bnr*M->PBpar[HOMmax+nd];
  M->PBsys[HOMmax-n] = anr*M->PBpar[HOMmax+nd];

  lat.SetPB(Fnum, Knum, n); lat.SetPB(Fnum, Knum, -n);

  if (prt)
    printf("set the n=%d component of %s to %e %e %e\n",
	   n, lat.elems[lat.Elem_GetPos(Fnum, Knum)]->PName,
	   bnr, M->PBpar[HOMmax+nd], M->PBsys[HOMmax+n]);
}


void set_bnr_sys_fam(LatticeType &lat, const int Fnum, const int n,
		     const double bnr, const double anr)
{
  int k;

  if (n < 1) {
    std::cout << "set_bnr_sys_fam: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  for (k = 1; k <= lat.GetnKid(Fnum); k++)
    set_bnr_sys_elem(lat, Fnum, k, n, bnr, anr);
}


void set_bnr_sys_type(LatticeType &lat, const int type, const int n,
		      const double bnr, const double anr)
{
  long int  k;
  MpoleType *M;

  if (n < 1) {
    std::cout << "set_bnr_sys_type: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  if (type >= Dip && type <= HOMmax) {
    for(k = 1; k <= lat.conf.Cell_nLoc; k++) {
      M = dynamic_cast<MpoleType*>(lat.elems[k]);
      if ((lat.elems[k]->Pkind == Mpole) && (M->n_design == type))
	set_bnr_sys_elem(lat, lat.elems[k]->Fnum, lat.elems[k]->Knum, n, bnr,
			 anr);
    }
  } else
    printf("Bad type argument to set_bnr_sys_type()\n");
}


void set_bnr_rms_elem(LatticeType &lat, const int Fnum, const int Knum,
		      const int n, const double bnr, const double anr,
		      const bool new_rnd)
{
  int       nd;
  MpoleType *M;

  bool prt = false;

  if (n < 1) {
    std::cout << "set_bnr_rms_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  M = dynamic_cast<MpoleType*>(lat.elems[lat.Elem_GetPos(Fnum, Knum)]);
  nd = M->n_design;
  // errors are relative to design values for (Dip, Quad, Sext, ...)
  if (nd == Dip) {
    M->PBrms[HOMmax+n] = bnr*M->Pirho; M->PBrms[HOMmax-n] = anr*M->Pirho;
  } else {
    M->PBrms[HOMmax+n] = bnr*M->PBpar[HOMmax+nd];
    M->PBrms[HOMmax-n] = anr*M->PBpar[HOMmax+nd];
  }

  if(new_rnd){
    if (normal) {
      M->PBrnd[HOMmax+n] = normranf(); M->PBrnd[HOMmax-n] = normranf();
    } else {
      M->PBrnd[HOMmax+n] = ranf(); M->PBrnd[HOMmax-n] = ranf();
    }
  }

  lat.SetPB(Fnum, Knum, n); lat.SetPB(Fnum, Knum, -n);

  if (prt) {
    printf("set_bnr_rms_elem:  Fnum = %d, Knum = %d, n = %d, n_design = %d"
	   ", new_rnd = %d, r_# = (%e, %e)\n",
	   Fnum, Knum, n, nd, new_rnd,
	   M->PBrnd[HOMmax+n], M->PBrnd[HOMmax-n]);
    printf("  (bnr, anr) = (%e, %e), PBrms = (%e, %e), PB = (%e, %e)\n",
	   bnr, anr, M->PBrms[HOMmax+n], M->PBrms[HOMmax-n],
	   M->PB[HOMmax+n], M->PB[HOMmax-n]);
  }
}


void set_bnr_rms_fam(LatticeType &lat, const int Fnum, const int n,
		     const double bnr, const double anr, const bool new_rnd)
{
  int k;

  if (n < 1) {
    std::cout << "set_bnr_rms_fam: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  for (k = 1; k <= lat.GetnKid(Fnum); k++)
    set_bnr_rms_elem(lat, Fnum, k, n, bnr, anr, new_rnd);
}


void set_bnr_rms_type(LatticeType &lat, const int type, const int n,
		      const double bnr, const double anr, const bool new_rnd)
{
  long int  k;
  MpoleType *M;

  if (n < 1) {
    std::cout << "set_bnr_rms_type: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  if (type >= Dip && type <= HOMmax) {
    for(k = 1; k <= lat.conf.Cell_nLoc; k++) {
      M = dynamic_cast<MpoleType*>(lat.elems[k]);
      if ((lat.elems[k]->Pkind == Mpole) && (M->n_design == type))
	set_bnr_rms_elem(lat, lat.elems[k]->Fnum, lat.elems[k]->Knum, n, bnr,
			 anr, new_rnd);
    }
  } else
    printf("Bad type argument to set_bnr_rms_type()\n");
}
