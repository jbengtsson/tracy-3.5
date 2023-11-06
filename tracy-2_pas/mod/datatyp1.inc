  const	gvecmax=1024; strlmax=80; matdim=6; HOMmax=10; Namelength=15;
	on=true; off=false; lasfnb=13;
	drift=0; wigl=1; mpole=2; cavity=3; marker=4; undef=5;
	Dipole=1; Quadrupole=2; Sextupole=3;
	DBnamelen=39;

  	nvmax=6;

  type	ivector2	= array [1..2] of integer;
	vector2		= array [1..2] of double;
	vector		= array [1..matdim] of double;
	matrix		= array [1..matdim] of vector;

	ivector		= array [1..nvmax] of integer;
	DAvect		= array [0..nvmax] of double;
	DAmap		= array [1..nvmax] of DAvect;

	lasftyp		= array [1..lasfnb] of integer;

	graphvect	= array [1..gvecmax] of double;

	strbuf		= array [1..strlmax] of char;
	tstring		= record
			    len : integer;
			    str : strbuf;
			  end;

	partsname	= array [1..Namelength] of char;
	DBnameType	= array [1..DBNamelen] of char;
	mpolArray	= Array [-HOMmax..HOMmax] of double;

	MpoleType	= Record
			    Pmethod	: integer;	{ Integration Method }
			    PN		: integer;	{ number of integration steps }
			    { Displacement Errors }
			    PdSsys  	: Vector2;	{ systematic [m]  }
			    PdSrms  	: Vector2;	{ rms [m] }
			    PdSrnd  	: Vector2;	{ random number }
			    { Tilt angle }
			    PdTpar	: double;	{ design [deg] }
			    PdTsys	: double;	{ systematic [deg] }
			    PdTrms	: double;	{ rms [deg] }
			    PdTrnd	: double;	{ random number }
			    { Multipole strengths }
			    PBpar	: mpolArray;	{ design }
			    PBsys	: mpolArray;	{ systematic }
			    PBrms	: mpolArray;	{ rms }
			    PBrnd	: mpolArray;	{ random number }
			    PB		: mpolArray;	{ total }
			    Porder	: Integer;	{ The highest order in PB } 
			    Pthick 	: integer;
			      PTx1		: double; { horizontal entrance angle [deg]}
			      PTx2		: double; { horizontal exit angle [deg]}
			      Pirho		: double; { 1/rho [1/m] }
			      Pa1, Pb1		: double; { corrections for tilt error of bend }
			      AU55, AD55	: matrix;
			  end;

	ElementType	= record
			    PName	: partsName;	{ Element name }
			    PL		: double;	{ Length[m] }
			    Pkind	: integer;
			    M		: MpoleType;
			  end;

	CellType	= Record
			    Fnum    : integer; { Element Family # }
			    Knum    : integer; { Element Kid # }
			    S       : double;  { distance from the beginning }
			    dS      : Vector2;
			    dT      : Vector2; { cos(dT), sin(dT) }
			    Elem    : ElementType;
			    Nu      : Vector2;
			    Alpha   : Vector2;
			    Beta    : Vector2;
			    Eta     : Vector2;
			    Etap     : Vector2;
			    BeamPos : Vector;
			  end;

	globvalrec	= record
			    dPcommon			: double;
			    dPparticle			: double;
			    maxampl			: vector2;
			    TotalTune	 		: vector2;
			    Omega	 		: double;
			    Chrom			: vector2;
			    Energy 			: double;
			    Cell_nLoc			: integer;
			    Element_nFam		: integer;
			    QforTune			: ivector2;
			    SforChrom			: ivector2;
			    QforDisp			: integer;
			    PosforDisp			: integer;
			    CODimax			: integer;
			    CODeps			: double;
			    CODvect			: vector;
			    cohv			: ivector2;
			    bpm				: integer;
			    OneTurnMat, Ascr, Ascrinv	: matrix;
			    MatMeth, Cavity_on		: boolean;
			  end;

	var	globval	: globvalrec;
		pi	: double;
