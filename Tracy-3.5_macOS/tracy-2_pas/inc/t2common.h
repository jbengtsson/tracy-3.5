  { interface }

  const Dip = 1; Quad = 2; Sext = 3;
	Horizontal = 1; Vertical = 2;
	graphvectmax = 4096;
	no = 1; nv = 6; nd2 = 6;
	fitvectmax = 200;

  type	iVector2	= Array [1..2] of integer;
	Vector2		= Array [1..2] of double;
	Vector3		= Array [1..3] of double;
	graphvect	= array [1..graphvectmax] of real;
	fitvect		= array [1..fitvectmax] of integer;

	globvalrec	= record
			    dPcommon			: double;
			    dPparticle			: double;
			    maxampl			: vector2;
			    TotalTune	 		: vector2;
			    Omega	 		: double;
			    Alphac	 		: double;
			    Chrom			: vector2;
			    Energy 			: double;
			    Cell_nLoc			: integer;
			    Elem_nFam			: integer;
			    CODimax			: integer;
			    CODeps			: double;
			    CODvect			: vector;
			    bpm				: integer;
			    OneTurnMat, Ascr, Ascrinv	: matrix;
			    MatMeth, Cavity_on		: boolean;
			    radiation, emittance	: boolean;
			    dE				: double;
			    rad, qfluct			: array [1..3] of double;
			    pathlength			: boolean;
			  end;

	statusrec	= record
			    tuneflag	: boolean;
			    chromflag	: boolean;
			    codflag	: boolean;
			    mapflag	: boolean;
			  end;

  var	globval			: [external] globvalrec;
  	status			: [external] statusrec;
	trace, cellconcat	: [external] boolean;
