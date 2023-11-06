  { interface }

{
  var   stable  : [external] boolean;
}

  PROCEDURE Cell_GetABGN(var M : matrix; var alpha, beta, gamma,
				  nu : vector2); external;

  procedure Cell_MatTwiss(i0, i1 : integer; var Ascr : matrix;
			  chroma, ring : boolean; dP : double); external;

  procedure Cell_DATwiss(i0, i1 : integer; var Ascr : DAmap;
			 chroma, ring : boolean; dP : double); external;

  procedure Ring_Getchrom(dP : double); external;

  procedure Ring_GetTwiss(chroma : boolean; dP : double); external;

  PROCEDURE Ring_Fittune(var nu : vector2; eps : double;
				  var nq : ivector2; var qf, qd : fitvect;
				  dkL : double; imax : integer); external;

  PROCEDURE Ring_Fitchrom(var ksi : vector2; eps : double;
				   var ns : ivector2; var sf, sd : fitvect;
				   dkpL : double; imax : integer); external;

  PROCEDURE Ring_FitDisp(pos : integer; eta, eps : double;
				  nq : integer; var q : fitvect;
				  dkL : double; imax : integer); external;

