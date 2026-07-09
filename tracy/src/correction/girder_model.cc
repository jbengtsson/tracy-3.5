// Cormisal (girder-based) alignment error model -- see
// correction/girder_model.h. Bodies extracted verbatim from param_data_type.

// GirderSetup/SetCorMis reporting + plot dumps (were #defines in param.h).
static const bool reportflag = true, plotflag = true;

namespace corr {

void girder_model::GirderSetup() {
  bool     giropen, ismag;
  double   s0, s1, s2, circ;
  long     ngir, i0, ic, i, countmag;
  CellType cell;
  char     elem[SymbolLength+1];
  FILE     *outf;
  char     fname[30];

  printf("Girder Setup \n");


// allocate the girders (if any)
// also enter mid pos and angle in lattice structure}
  ngir = 0;
  circ = 0; giropen=false;  
  for (i = 0; i <= globval.Cell_nLoc; i++) {
 
    if (i == ilatmax) {
      printf("i %ld exceeds %d\n", i, ilatmax-1);
      exit(1);
    }

    Lattice[i].igir=-1;
 
    getelem(i, &cell);
    
    circ=circ+cell.Elem.PL;

    if ((cell.Elem.PName[0] == 'g') && (cell.Elem.PName[1]=='t')
	&& (cell.Elem.PName[2] == 'y')
    && (cell.Elem.PName[3]=='p')
	&& ((cell.Elem.PName[4] == '0') || (cell.Elem.PName[4]=='1'))) {
      if (giropen) {
// if girder is open, close it:
        giropen = false;
        if (cell.Elem.PName[4]=='0')
	  Girder[ngir-1].gco[1]=0;
	else
	  {Girder[ngir-1].gco[1]=1;}
        Girder[ngir-1].gsp[1]=circ;  
        Girder[ngir-1].ilat[1]=i;
      } else {
// if girder is not open, open a new girder:
        ngir++;
	if (ngir == igrmax)  {
          printf("ngir %ld exceeds %d\n", ngir, igrmax-1); exit(1);
        }
        giropen=true;
        Girder[ngir-1].gdx[0]=0; 
        Girder[ngir-1].gdx[1]=0; 
        Girder[ngir-1].gdy[0]=0; 
        Girder[ngir-1].gdy[1]=0; 
        Girder[ngir-1].gdt=0;
        if (cell.Elem.PName[4]=='0')
	  Girder[ngir-1].gco[0]=0;
	else
	  Girder[ngir-1].gco[0]=1;
        Girder[ngir-1].gsp[0]=circ;  
        Girder[ngir-1].ilat[0]=i;
        Girder[ngir-1].igir[0]=-1; 
        Girder[ngir-1].igir[1]=-1;
        Girder[ngir-1].level=1;
      }
    }
    Lattice[i].smid=circ-cell.Elem.PL/2;
  }//for

  for (i=0;i<ngir;i++)
    for (ic=Girder[i].ilat[0]; ic<=Girder[i].ilat[1]; ic++)
      Lattice[ic].igir=i;

  NGirderLevel[0]=ngir;
  
  // find compounds, i.e. elements which are to be treated as one block w.r.t.
  // misalignment two types: bracketed by girder type 2,3 or series of magnets
  // w/o space between. first select all compound elements, defined by bracket
  // of type 2,3 girders:
   
  s1=0; s2=0; giropen=false;
  for (i = 0; i <= globval.Cell_nLoc; i++) {

    getelem(i, &cell);
    s2=s1+cell.Elem.PL;
    if ((cell.Elem.PName[0]=='g') && (cell.Elem.PName[1]=='t')
	&& (cell.Elem.PName[2]=='y') && (cell.Elem.PName[3]=='p')
	&& ((cell.Elem.PName[4]=='2')||(cell.Elem.PName[4]=='3'))) {

      if (giropen) {
	// if compound is open, close it:
        giropen=false;
        Girder[ngir-1].gsp[1]=s2;
        Girder[ngir-1].ilat[1]=i;
        Girder[ngir-1].igir[1]=Lattice[i].igir; 
        if (cell.Elem.PName[4]=='2')
	  Girder[ngir-1].gco[1]=2;
	else
	  Girder[ngir-1].gco[1]=3; 
      } else {
	// if compound is not open, open a new one:
        ngir++;
        giropen=true;
        Girder[ngir-1].gsp[0]=s2;
        Girder[ngir-1].ilat[0]=i;
        Girder[ngir-1].igir[0]=Lattice[i].igir; 
        Girder[ngir-1].gco[0]=0;
        Girder[ngir-1].level=2;
        if (cell.Elem.PName[4]=='2')
	  Girder[ngir-1].gco[0]=2;
	else
	  Girder[ngir-1].gco[0]=3; 
      }
    }
    s1=s2;
  }//for


  for (i=NGirderLevel[0];i<ngir;i++)
    for (ic=Girder[i].ilat[0];ic<=Girder[i].ilat[1];ic++)
      Lattice[ic].igir=i;
  NGirderLevel[1]=ngir;

  // make a compound element if we have a series of magnets with no gap between,
  // i.e. sext|ch|cv|sext
  s0=0; s1=0; s2=0; i0=0;
  giropen=false; countmag=0;
  for (i = 0; i <= globval.Cell_nLoc; i++) {
    getelem(i, &cell);
    s2=s1+cell.Elem.PL;
 
    ismag= (cell.Elem.Pkind==Mpole);

    if (giropen) {
      // keep s0, was set when merging started
      if  ((ismag) || (fabs(cell.Elem.PL)< seps) ) {  // continue merge
        if (ismag) { countmag++;};
      } else { //stop merge
        if (countmag>1) {
          ngir++;
          Girder[ngir-1].gsp[0]=s0;
          Girder[ngir-1].ilat[0]=i0;
          Girder[ngir-1].igir[0]=Lattice[i0].igir;
          Girder[ngir-1].gco[0]=0;
          Girder[ngir-1].gsp[1]=s1;
          Girder[ngir-1].ilat[1]=i-1;
          Girder[ngir-1].igir[1]=Lattice[i-1].igir;
          Girder[ngir-1].gco[1]=0;
          Girder[ngir-1].level=3;
	  //          countmag=0;
        } 
        countmag=0;
        giropen=false;
      } 
    } else {
      if (ismag) {
        giropen=true; //start a new merge
        s0=s1; i0=i;
        countmag=1;
      }
    }
    s1=s2;
  }//for


  for (i=NGirderLevel[1];i<ngir;i++){
    for (ic=Girder[i].ilat[0];ic<=Girder[i].ilat[1];ic++){Lattice[ic].igir=i;}
  }
  NGirderLevel[2]=ngir;

  for (i = 0; i <= globval.Cell_nLoc; i++) {
    getelem(i, &cell);
    if (!(cell.Elem.Pkind==Mpole)) { Lattice[i].igir=-1;}
  }

  if (reportflag) {
    if (ngir>0) {
      for (i=0;i<ngir;i++) {
        printf( "gir %ld lev %ld igir %ld %ld co %ld %ld sp %f %f \n", 
		i, Girder[i].level,Girder[i].igir[0],Girder[i].igir[1],
		Girder[i].gco[0],Girder[i].gco[1],Girder[i].gsp[0],
		Girder[i].gsp[1]);
      }
    } 
    for (i = 0; i <= globval.Cell_nLoc; i++) {
      getelem(i, &cell);
      TracyStrcpy( elem, cell.Elem.PName);
     if (Lattice[i].igir > -1){
       printf("pos %ld %s %f gir %ld lev %ld",
	      i, elem, Lattice[i].smid, Lattice[i].igir,
	      Girder[Lattice[i].igir].level);
        if (Girder[Lattice[i].igir].level >=2){ 
          printf(" --> %ld %ld \n", Girder[Lattice[i].igir].igir[0],
		 Girder[Lattice[i].igir].igir[1]);
        } else { printf("\n");}
      } else {printf("pos %ld %s --- ---- \n", i,elem);}
    }

    if (plotflag) {
      strcpy(fname, "gsetup.plt");
      outf = fopen(fname,"w" );
 
      if (ngir>0) {
        fprintf(outf, "%ld %ld %ld %ld %f3 \n",
		globval.Cell_nLoc, NGirderLevel[0], NGirderLevel[1],
		NGirderLevel[2], s2);
        for (i=0;i<ngir;i++)
          fprintf(outf, "%ld %ld %ld %ld %ld %ld %f3 %f3 \n", 
		  i, Girder[i].level, Girder[i].igir[0],Girder[i].igir[1],
		  Girder[i].gco[0],Girder[i].gco[1],Girder[i].gsp[0],
		  Girder[i].gsp[1]);
      } 

      s1=0;
      for (i = 0; i <= globval.Cell_nLoc; i++) {
        getelem(i, &cell);
        s2=s1+cell.Elem.PL;
        ismag= (cell.Elem.Pkind==Mpole);
        if (ismag) {
          TracyStrcpy( elem, cell.Elem.PName);
          fprintf(outf, "%ld %ld %f3 %f3 %s\n",
		  i, Lattice[i].igir, s1, s2, elem);
        } 
	s1=s2;
      }
      if (outf != NULL) fclose(outf);
      outf = NULL;
    }
  }

}


void girder_model::SetCorMis(double gxrms, double gyrms, double gtrms,
				double jxrms, double jyrms, double exrms,
				double eyrms, double etrms, double rancutx,
				double rancuty, double rancutt, long iseed)
{
  double   jdx, jdy, ggxrms, ggyrms, r, g3dx, g3dy, g3dt, gelatt, att, dx, dy;
  double   dt, s1, s2;
  long     i, isup;
  CellType cell;
  char     elem[SymbolLength+1];
  FILE     *outf;
  char     fname[30];

  printf("SetCorMis: initializing seed %ld\n", iseed);
  iniranf(iseed);

  // TODO(girder-zero-span): the girder-support interpolations below all divide
  // by a girder's span (Girder[isup].gsp[1]-Girder[isup].gsp[0]). A zero-span
  // girder makes r = inf and 0*(1-inf) = NaN, corrupting the misalignment of
  // every supported element even at zero error amplitude -> lost beam.
  // GirderSetup can produce zero-span level-3 girders from runs of consecutive
  // zero-length magnets (e.g. adjacent zero-length correctors). Guard the span
  // (skip support / treat as free end when span <= seps) as part of the owed
  // physics validation of the girder->element translation. Not fixed here to
  // keep this a behavior-preserving extraction.

  /*
     set misalignments to girder ends:
     simple shortcut for joints: just use prev girder and add +/- joint play 
     if prev-girder-end and this-girder-start both have link flag gco=1
     no further options like in OPA
  */

  for (i=0;i< NGirderLevel[0];i++){
    setrancut(rancutx);
    Girder[i].gdx[0] = gxrms*normranf();
    Girder[i].gdx[1] = gxrms*normranf();
    setrancut(rancuty);
    Girder[i].gdy[0] = gyrms*normranf();
    Girder[i].gdy[1] = gyrms*normranf();
    setrancut(rancutt);
    Girder[i].gdt    = gtrms*normranf();
    if ((Girder[i].gco[0]==1) && (i>0)) {
      if (Girder[i-1].gco[1]==1) {
        setrancut(rancutx);
        jdx=jxrms*normranf();
        setrancut(rancuty);
        jdy=jyrms*normranf();
        Girder[i].gdx[0] = Girder[i-1].gdx[1]+jdx;
        Girder[i].gdy[0] = Girder[i-1].gdy[1]+jdy;
        Girder[i-1].gdx[1] -= jdx;
        Girder[i-1].gdy[1] -= jdy;
      }
    }
  }

  /*
  set misalignment for level 2 girder, which are supported by other girder.
  no, use joint play for connection of level 2 girder to level 1 girder
  [no further error applied for level 2, since the error is given by the
  supporting girders,
  and elements may receive additional individual errors later]
  */
  ggxrms=jxrms; 
  ggyrms=jyrms;

  for (i=NGirderLevel[0]; i<NGirderLevel[1];i++) {
    isup= Girder[i].igir[0]; 
    if (isup > -1) {
      r = (Girder[i].gsp[0]-Girder[isup].gsp[0])
	/(Girder[isup].gsp[1]-Girder[isup].gsp[0]);
      setrancut(rancutx);
      Girder[i].gdx[0] = Girder[isup].gdx[0]*(1-r)+Girder[isup].gdx[1]*r
	+ ggxrms*normranf();
      setrancut(rancuty);
      Girder[i].gdy[0] = Girder[isup].gdy[0]*(1-r)+Girder[isup].gdy[1]*r
	+ ggyrms*normranf();
      if (Girder[i].gco[0]==3) {
	Girder[i].gdt=Girder[isup].gdt;
      } else {
	setrancut(rancutt);Girder[i].gdt=gtrms*normranf();
      }
      // printf("\nGir up %ld %f %f %f %f %f \n",
      // 	     i, Girder[i].gsp[0], Girder[i].gsp[1], Girder[i].gdx[0]*1e6,
      // 	     Girder[i].gdx[1]*1e6, r);
      // printf("  supp %ld %f %f %f %f \n",
      // 	     isup, Girder[isup].gsp[0], Girder[isup].gsp[1],
      // 	     Girder[isup].gdx[0]*1e6, Girder[isup].gdx[1]*1e6);

      /* 
	 contact 3 (2-point) transmits roll error from supporting girder,
	 contact 2 (1-point) is free.
	 if contact 2 -> set gdt, but will be overwritten if other end is
	 contact 3
	 if other end is also contact 2, this value is taken, because gdt then
	 is arbitrary
	 if ends are free, treat like contact 0
	 if end 1 only is free, then check if gdt may have been set at end 0
      */
    } else {
      setrancut(rancutx);
      Girder[i].gdx[0] = gxrms*normranf();
      setrancut(rancuty);
      Girder[i].gdy[0] = gyrms*normranf();
      setrancut(rancutt);
      Girder[i].gdt    = gtrms*normranf();
    }
    isup= Girder[i].igir[1]; 
    if (isup > -1) {
      r=(Girder[i].gsp[1]-Girder[isup].gsp[0])
	/(Girder[isup].gsp[1]-Girder[isup].gsp[0]);
      setrancut(rancutx);
      Girder[i].gdx[1] = Girder[isup].gdx[0]*(1-r)+Girder[isup].gdx[1]*r
	+ ggxrms*normranf();
      setrancut(rancuty);
      Girder[i].gdy[1] = Girder[isup].gdy[0]*(1-r)+Girder[isup].gdy[1]*r
	+ ggyrms*normranf();
      if (Girder[i].gco[1] ==3) {Girder[i].gdt = Girder[isup].gdt;}
      // printf("Gir dn %ld %f %f %f %f %f \n",
      // 	     i, Girder[i].gsp[0], Girder[i].gsp[1], Girder[i].gdx[0]*1e6,
      // 	     Girder[i].gdx[1]*1e6, r);
      // printf("  supp %ld %f %f %f %f \n", isup, Girder[isup].gsp[0],
      // 	     Girder[isup].gsp[1], Girder[isup].gdx[0]*1e6,
      // 	     Girder[isup].gdx[1]*1e6);
    } else {
      setrancut(rancutx);
      Girder[i].gdx[1] = gxrms*normranf();
      setrancut(rancuty);
      Girder[i].gdy[1] = gyrms*normranf();
      if (Girder[i].igir[0] == -1) {
	setrancut(rancutt); Girder[i].gdt = gtrms*normranf();
      }
    }
  }

  // set misalignment for level 3 girder, which are compound elements, which
  // have common element displacement error.
  
  gelatt=1.0;

  for (i=NGirderLevel[1]; i< NGirderLevel[2];i++) {
    setrancut(rancutx);
    g3dx = exrms*normranf()*gelatt;
    setrancut(rancuty);
    g3dy = eyrms*normranf()*gelatt;
    setrancut(rancutt);
    g3dt = etrms*normranf()*gelatt;
    isup= Girder[i].igir[0]; 
    if (isup > -1) {
      r=(Girder[i].gsp[0]-Girder[isup].gsp[0])
	/(Girder[isup].gsp[1]-Girder[isup].gsp[0]);
      Girder[i].gdx[0] =Girder[isup].gdx[0]*(1-r)+Girder[isup].gdx[1]*r + g3dx;
      Girder[i].gdy[0] =Girder[isup].gdy[0]*(1-r)+Girder[isup].gdy[1]*r + g3dy;
      Girder[i].gdt    =Girder[isup].gdt+g3dt; // presume contact 3, rigid
                                               // connection
    } else {
      setrancut(rancutx);
      Girder[i].gdx[0] =exrms*normranf();
      setrancut(rancuty);
      Girder[i].gdy[0] =eyrms*normranf();
      setrancut(rancutt);
      Girder[i].gdt    =etrms*normranf();
    }
    isup= Girder[i].igir[1]; 
    if (isup > -1) {
      r =(Girder[i].gsp[1]-Girder[isup].gsp[0])
	/(Girder[isup].gsp[1]-Girder[isup].gsp[0]);
      Girder[i].gdx[1] =Girder[isup].gdx[0]*(1-r)+Girder[isup].gdx[1]*r + g3dx;
      Girder[i].gdy[1] =Girder[isup].gdy[0]*(1-r)+Girder[isup].gdy[1]*r + g3dy;
      Girder[i].gdt    =Girder[isup].gdt+g3dt; // should be on same girder and
                                               // give same result
    } else {
      setrancut(rancutx);
      Girder[i].gdx[1] =exrms*normranf();
      setrancut(rancuty);
      Girder[i].gdy[1] =eyrms*normranf();
      if (isup ==-1) {setrancut(rancutt);Girder[i].gdt =etrms*normranf();}
    }
  }


// set misalignments of elements on girders:
  for (i = 0; i <= globval.Cell_nLoc; i++) {
    getelem(i, &cell);
      
    if (cell.Elem.Pkind==Mpole) {
      if ((cell.Fnum != globval.hcorr) && (cell.Fnum != globval.vcorr)) {
        setrancut(rancutx);
	dx =exrms*normranf();
        setrancut(rancuty);
	dy =eyrms*normranf();
        setrancut(rancutt);
	dt =etrms*normranf();
	isup=Lattice[i].igir;

	if (isup > -1) {
	  if (Girder[isup].level >= 2) {att=0;} else {att=gelatt;}
	  r =(Lattice[i].smid-Girder[isup].gsp[0])
	    /(Girder[isup].gsp[1]-Girder[isup].gsp[0]);
	  dx =dx*att+Girder[isup].gdx[0]*(1-r)+Girder[isup].gdx[1]*r;
	  dy =dy*att+Girder[isup].gdy[0]*(1-r)+Girder[isup].gdy[1]*r;
	  dt =dt*att+Girder[isup].gdt;
	}

	cell.Elem.M->PdSsys[0] = dx;
	cell.Elem.M->PdSsys[1] = dy;
	cell.Elem.M->PdTsys    = dt;

	putelem(i, &cell);

	Mpole_SetdS(cell.Fnum, cell.Knum);
	Mpole_SetdT(cell.Fnum, cell.Knum);
      }
    }
  }

  if (plotflag) {

    snprintf(fname, sizeof(fname), "cormis_%ld.plt",iseed);
    outf = fopen(fname,"w" );

    if (NGirderLevel[2] > 0) {
      fprintf(outf, "%ld %ld %ld %ld \n",
	      globval.Cell_nLoc, NGirderLevel[0], NGirderLevel[1],
	      NGirderLevel[2]);
      for (i=0;i<NGirderLevel[2];i++) {
        fprintf(outf, "%ld %f %f %f %f %f %f %f  \n", 
         Girder[i].level, Girder[i].gsp[0],Girder[i].gsp[1],
		Girder[i].gdx[0]*1e6,Girder[i].gdx[1]*1e6,
		Girder[i].gdy[0]*1e6,Girder[i].gdy[1]*1e6, Girder[i].gdt*1e6);
      }
    } 


    s1=0;
    for (i = 0; i <= globval.Cell_nLoc; i++) {
      getelem(i, &cell);
      s2=s1+cell.Elem.PL;
      if (cell.Elem.Pkind==Mpole) {
        dx = cell.dS[0];
        dy = cell.dS[1];
        dt = atan(cell.dT[1]/cell.dT[0]);
        TracyStrcpy( elem, cell.Elem.PName);
        fprintf(outf,"%ld %f %f %f %f %f %s \n", i,  s1, s2, dx*1e6, dy*1e6, dt*1e6, elem); 
      }
      s1=s2;
    }
    if (outf != NULL) fclose(outf);
    outf = NULL;
  }
}


void CorMis_in(double *gdxrms, double *gdzrms, double *gdarms, double *jdxrms, double *jdzrms, double *edxrms, double *edzrms, double *edarms, double *bdxrms, double *bdzrms, double *bdarms, double *rancutx, double *rancuty, double *rancutt, long *iseed, long *iseednr)
{
  char a;
  long i;
  bool includeMON;
  FILE *cinf;
  
  cinf = fopen( "cormis.dat" , "r");

  printf("Apply errors also to BPMs (name=MON)? (Y/n) \n");
  fscanf(cinf, "%c%*[^\n]", &a);
  getc(cinf);
  if (a == '\n')
    a = ' ';
  includeMON = (a != 'n' && a != 'N');
  if (includeMON)
    printf("BPMs with errors\n");
  else
    printf("BPMs without errors\n");

  printf("\nInput of error amplitudes for gaussian errors\n");
  printf("----------------------------------------------------------\n");
  printf("Give rms errors for displacements in micron, horizontal and"
	 " vertical:\n");
  printf("__ Absolute displacement of girder joints and ends :\n");
  fscanf(cinf, "%lg%lg%lg%*[^\n]", gdxrms, gdzrms, gdarms);
  getc(cinf);
  printf("% .5E% .5E% .5E\n", *gdxrms, *gdzrms, *gdarms);
  printf("__ Relative displacement WITHIN girder joints (joint play):\n");
  fscanf(cinf, "%lg%lg%*[^\n]", jdxrms, jdzrms);
  getc(cinf);
  printf("% .5E% .5E\n", *jdxrms, *jdzrms);
  printf("__ Relative displacement of elements ON a girder:\n");
  fscanf(cinf, "%lg%lg%lg%*[^\n]", edxrms, edzrms, edarms);
  getc(cinf);
  printf("% .5E% .5E% .5E\n", *edxrms, *edzrms, *edarms);
  if (includeMON) {
    printf("__ Relative displacement of BPMs:\n");
    fscanf(cinf, "%lg%lg%lg%*[^\n]", bdxrms, bdzrms, bdarms);
    getc(cinf);
    printf("% .5E% .5E% .5E\n", *bdxrms, *bdzrms, *bdarms);
  } else {
    (*bdxrms)=(*bdzrms)=(*bdarms)=0.;
  }
  printf("__ Gaussian cut:\n");
  fscanf(cinf, "%lg%lg%lg%*[^\n]", rancutx, rancuty, rancutt);
  getc(cinf);
  printf("% .5E% .5E% .5E\n", *rancutx, *rancuty, *rancutt);
  printf("__ init seed values:\n");
  fscanf(cinf, "%ld", iseednr);
  if (*iseednr > iseednrmax) {
    printf("Iseednr %ld exceeds %d\n", *iseednr, iseednrmax); exit(1);
  }
  for (i=0; i<*iseednr; i++) {
    fscanf(cinf, "%ld", &iseed[i]);
    printf("%ld ", iseed[i]);
  }
  fscanf(cinf,"%*[^\n]"); printf("\n\n");
  getc(cinf);

  printf("rms Girder error  : dx=%5.0f um, dz=%5.0f um, da=%5.0f udeg\n",
	 *gdxrms, *gdzrms, *gdarms);
  printf("rms Joint  error  : dx=%5.0f um, dz=%5.0f um\n", *jdxrms, *jdzrms);
  printf("rms Element error : dx=%5.0f um, dz=%5.0f um, da=%5.0f udeg\n",
	 *edxrms, *edzrms, *edarms);
  if (includeMON)
    printf("rms BBA error     : dx=%5.0f um, dz=%5.0f um, da=%5.0f udeg\n",
	   *bdxrms, *bdzrms, *bdarms);
  printf("\n");
  printf("Gaussian cut      : cutx=%5.0f, cuty=%5.0f, cutt=%5.0f sigma\n",
	 *rancutx, *rancuty, *rancutt);
  printf("init seed values  : %ld seeds= ", *iseednr);
  for (i=0; i<*iseednr; i++)
    printf("%ld ",iseed[i]);
  printf("\n\n");

  *gdxrms=(*gdxrms)*1E-6;  *gdzrms=(*gdzrms)*1E-6; *gdarms=(*gdarms)*1E-6;
  *jdxrms=(*jdxrms)*1E-6;  *jdzrms=(*jdzrms)*1E-6;
  *edxrms=(*edxrms)*1E-6;  *edzrms=(*edzrms)*1E-6; *edarms=(*edarms)*1E-6;
  *bdxrms=(*bdxrms)*1E-6;  *bdzrms=(*bdzrms)*1E-6; *bdarms=(*bdarms)*1E-6;
}

}  // namespace corr
