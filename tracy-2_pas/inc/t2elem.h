  { interface }

  const	HOMmax	=  15;

	Meth_Linear = 0; Meth_First = 1; Meth_Second = 2; Meth_Fourth = 4;
	Meth_genfun = 5;

	{ maximum number of families for Elem_NFam }
	Elem_nFamMax = 200;

	{ maximum number of elements for (Cell_nLoc)}
	Cell_nLocMax = 2500;

	{ maximum number of kids }
	nKidMax = 500;

	DBNameLen = 39;

  type	partsName	= packed array [1..Namelength] of char;
	pthicktype	= (thick, thin);
	partsKind	= (drift, Wigl, Mpole, Cavity, marker, undef);
	DBNameType      = packed array [1..DBNameLen] of char;

	mpolArray	= Array [-HOMmax..HOMmax] of double;

	DriftType	= Record
			    D55: matrix;
			  end;
	DriftPtr	= ^DriftType;

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
			    case Pthick : pthicktype of
			      Thin : ({ nothing });
			      Thick: ({ Bending Angles}
				      PTx1		: double; { horizontal entrance angle [deg]}
				      PTx2		: double; { horizontal exit angle [deg]}
				      Pgap		: double; { total magnet gap [m] }
				      Pirho		: double; { 1/rho [1/m] }
				      Pc0, Pc1, Ps1	: double; { corrections for tilt error of bend }
				      AU55, AD55	: matrix);
			  end;
	MpolePtr	= ^MpoleType;

	WigglerType	= Record
			    Pmethod	: integer;	{ Integration Method }
			    PN		: integer;	{ number of integration steps }
			    { Displacement Error }
			    PdSsys	: Vector2;	{ systematic [m]  }
			    PdSrms	: Vector2;	{ rms [m] }
			    PdSrnd	: Vector2;	{ random number }
			    { Tilt angle }
			    PdTpar	: double;	{ design [deg] }
			    PdTsys	: double;	{ systematic [deg] }
			    PdTrms	: double;	{ rms [deg] }
			    PdTrnd	: double;	{ random number }
			    { Strength }
			    Plperiod	: double;	{ Length Period [m] }
			    Pnperiod	: integer;	{ Number of periods }
			    PBoBrho	: double;	{ B/Brho }
                            PKx         : double;       { kx }
			    PBW		: mpolArray;
			    W55		: matrix;
			    Porder	: Integer;	{ The highest order in PB } 
			  end;
	WigglerPtr	= ^WigglerType;

	CavityType	= Record
			    Pvolt	: double;	{ Vrf [V] }
			    Pfreq	: double;	{ Vrf [Hz] }
			    Ph		: integer;	{ harmonic number }
			  end;
	CavityPtr	= ^CavityType;

	elemtype	= RECORD
			    PName	: partsName;	{ Element name }
			    PL		: double;	{ Length[m] }
			    case Pkind : PartsKind of
			      Drift : (D : DriftPtr);
			      Mpole : (M : MpolePtr);
			      Wigl  : (W : WigglerPtr);
			      Cavity: (C : CavityPtr);
			      Marker: ({ nothing });
			  end;
   
	ElemFamType	= Record
			    ElemF	: elemtype;
			    nKid	: integer;
			    KidList	: Array [1..nKidmax] of integer;
			    NoDBN	: integer;
			    DBNlist	: Array [1..nKidmax] of DBNameType;
			  end;

	CellType	= Record
			    Fnum	: integer; { Element Family # }
			    Knum	: integer; { Element Kid # }
			    S		: double;  { position }
			    dS		: Vector2;
			    dT		: Vector2; { cos(dT), sin(dT) }
			    Elem	: elemtype;
			    Nu		: Vector2;
			    Alpha	: Vector2;
			    Beta	: Vector2;
			    Eta		: Vector2;
			    Etap	: Vector2;
			    BeamPos	: Vector;
			  end;

  var	ElemFam					: [external] Array [1..Elem_nFamMax] of ElemFamType;
	Cell					: [external] Array [0..Cell_nLocmax] of CellType; 
	Fdrift1, Fkick1, Fdrift2,  Fkick2	: [external] double;
	crad, cfluc				: [external] double;
 

  function GetnKid(Fnum1 : integer) : integer; external;

  function Elem_GetPos(Fnum1, Knum1 : integer) : integer; external;

  procedure getelem(i : integer; var cellrec : celltype); external;

  procedure putelem(i : integer; var cellrec : celltype); external;

  procedure GtoL(var X : vector; var S, T : vector2;
		 c0, c1, s1 : double); external;
  { from Global to Local coordinate }

  procedure GtoL_M(var X : matrix; var T : vector2); external;
  { from Global to Local coordinates }

  procedure DAgtoL(var map : DAmap; var S, T : vector2;
		   c0, c1, s1 : double); external;
  { from Global to Local coordinate }

  procedure LtoG(var X : vector; var S, T : vector2;
		 c0, c1, s1 : double); external;
  { from Local to Global coordinate }

  procedure LtoG_M(var X : matrix; var T : vector2); external;
  { from Local to Global coordinate }

  procedure DALtoG(var map : DAmap; var S, T : vector2;
		   c0, c1, s1 : double); external;
  { from Local to Global coordinate }

  function GetNofDBN(elem : integer) : integer; external;
  { Returns the number of DBName defined for the elem-th element family }
 
  procedure GetDBN(elem, kid : integer; var Name : DBNameType); external;
  { Returns the kid-th database name of the elem-th element family }

  procedure InitFkick; external;
  { Init constants for 4-th order integrator }


  procedure Elem_Print(var f : text; Fnum1 : integer); external;
  
  { If Fnum1=0 then print all the elements }

  function Elem_GetKval(Fnum1, Knum1, Order : integer) : double; external;

  procedure Elem_GetBend(Fnum1, Knum1 : integer;
				     var Tx, Tx1, Tx2 : double); external;

  function Elem_GetOrd(Fnum1, Knum1 : integer) : integer; external;


  procedure MulLsMat(var a, b : matrix); external;

  procedure LinsTrans(var a : matrix; var b : vector); external;


  procedure Drift_Alloc(var Elem : elemtype); external;

  procedure Drift_Print(var f : text; Fnum1 : integer); external;

  PROCEDURE Drift_SetMatrix(Fnum1, Knum1 : integer); external;

  Procedure Drift_Init(Fnum1 : integer); external;

  Procedure Drift_Pass(var Cell : CEllType; var x : vector); external;

  Procedure Drift_Pass_M(var Cell : CEllType; var xref : vector;
				  var x : matrix); external;

  Procedure Drift_DApass(var Cell : CEllType; var map : DAmap); external;


  procedure thinkick(Order : integer; var MB : mpolarray;
		     L, irho : double; pthick : pthicktype;
		     var x : vector); external;

  procedure thinkick_M(Order : integer; var MB : mpolarray;
		       L, irho : double; pthick : pthicktype;
		       var xref : vector;
		       var x : matrix); external;

  procedure DAthinkick(Order : integer; var MB : mpolarray;
		       L, irho : double; pthick : pthicktype;
		       var map : DAmap); external;


  procedure Mpole_Alloc(var Elem:elemtype); external;

  PROCEDURE Mpole_SetMatrix(Fnum1, Knum1 : integer; K : double); external;

  procedure Mpole_Init(Fnum1 : integer); external;

  PROCEDURE Mpole_SetIntegrator(Fnum1, Method, Segment : integer ); external;

  procedure Mpole_Print(var f : text; Fnum1 : integer); external;

  procedure Mpole_GetBend(Fnum1, Knum1 : integer;
				   var Tx, Tx1, Tx2 : double); external;

  PROCEDURE Mpole_SetPB(Fnum1, Knum1, Order : integer); external;

  function Mpole_GetPB(Fnum1, Knum1, Order : integer) : double; external;

  function Mpole_GetOrder(Fnum1, Knum1 : integer) : integer; external;

  procedure Mpole_GetPmeth(Fnum1, Knum1 : integer;
				    var Pmethodt, PNt : integer;
				    var Fdrift1t, Fkick1t, Fdrift2t,
				    Fkick2t : double); external;


  PROCEDURE Mpole_SetdS(Fnum1, Knum1 : integer); external;

  PROCEDURE Mpole_SetdT(Fnum1, Knum1 : integer); external;

  procedure Mpole_Pass(var Cell : CEllType; var x : vector); external;

  procedure Mpole_Pass_M(var Cell : CEllType; var xref : vector;
				  var x : matrix); external;

  procedure Mpole_DApass(var Cell : CEllType; var map : DAmap); external;


  procedure Wiggler_Alloc(var Elem : elemtype);external;

  procedure Wiggler_Print(var f : text; Fnum1 : integer);external;

  PROCEDURE Wiggler_SetMatrix(Fnum1, Knum1 : integer;
				       L, lambda, k0, kx : double); external;

  PROCEDURE Wiggler_SetPB(Fnum1, Knum1, Order : integer);external;

  PROCEDURE Wiggler_Init(Fnum1 : integer);external ;

  procedure wiggler_Pass(var Cell : CEllType; var X : vector);external;

  procedure wiggler_Pass_M(var Cell : CEllType; var xref : vector;
				    var x : matrix);external;

  procedure wiggler_DApass(var Cell : CEllType; var map : DAmap);external;


  PROCEDURE Cav_Alloc(var Elem:elemtype); external;

  PROCEDURE Cav_Set(Fnum1, Knum1 : integer; Freq, Volt : double;
		    h : integer); external;

  PROCEDURE Cav_Get(Fnum1, Knum1 : integer; var Freq, Volt : double;
		    var h : integer); external;

  Procedure Cav_Init(Fnum1 : integer); external;

  PROCEDURE Cav_Pass(var Cell : CEllType; var X : vector); external;

  PROCEDURE Cav_DApass(var Cell : CEllType; var map : DAmap); external;

  Procedure Marker_Init(Fnum1 : integer); external;

  PROCEDURE Marker_Pass(var Cell : CEllType; var X : vector); external;
