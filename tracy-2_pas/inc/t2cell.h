  { interface }

  const maxtransfmat=500; maxkicks=5;
        debug=false;

  var   ntransfmat      : [external] integer;
        transfmat       : [external] array [1..maxtransfmat] of matrix;
        kicks           : [external] array [1..maxtransfmat, 1..maxkicks] of integer;

  procedure Cell_SetdP(dP:double); external;

  procedure Cell_Init; external;

  procedure Elem_Pass(i : integer; var x : vector); external;

  procedure Elem_Pass_M(i : integer; var xref : vector;
				    var x : matrix); external;

  procedure Elem_DApass(i : integer; var map : DAmap); external;

  procedure Cell_Pass(i0, i1 : integer; var x : vector;
		      var lastpos : integer); external;
     { Track particle from i0 to i1, set x[5]:=dP }

  procedure Cell_Pass_M(i0, i1 : integer; var xref : vector;
		      var mat : matrix; var lastpos : integer); external;
     { Track matrix from i0 to i1 around ref. orbit }

  procedure Cell_DApass(i0, i1 : integer; var map : DAmap;
				 var lastpos : integer); external;
     { Track matrix from i0 to i1 around ref. orbit }

  procedure Cell_Concat(dP : double); external;

  procedure Cell_fPass(var x : vector; var lastpos : integer); external;

  procedure Cell_fPass_M(var xref : vector; var mat : matrix;
				  var lastpos : integer); external;

  procedure Cell_GetCOD(imax : integer; eps, dP : double;
				 var lastpos : integer); external;
