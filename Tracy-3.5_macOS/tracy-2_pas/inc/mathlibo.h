  const matdim	= 6;

  type{	double		= real;}
	vector		= array [1..matdim] of double;
	matrix		= array [1..matdim] of vector;

{  var	pi		: [external] double;
        rseed0, rseed   : [external] integer;
        normcut_        : [external] double;
}

  { Extensions to the standard functions }

{  FUNCTION dble(x : real) : double; external;}

{  FUNCTION sngl(x : double) : real; external;}

  FUNCTION min_(x1, x2 : double) : double; external;

  FUNCTION max_(x1, x2 : double) : double; external;

  FUNCTION pwr(x, y : double) : double; external;

  FUNCTION tan_(x : double) : double; external;

  FUNCTION cosh_(x : double) : double; external;

  FUNCTION sinh_(x : double) : double; external;

  FUNCTION tanh_(x : double) : double; external;

  PROCEDURE iniranf(i : integer); external;

  PROCEDURE newseed; external;

  FUNCTION ranf : double; external;

  PROCEDURE setrancut(cut : double); external;

  FUNCTION normranf : double; external;

  { Conversion routines }

  FUNCTION dtor (d : double) : double; external;

  FUNCTION sign(x : double) : integer; external;

  Function GetAngle(x, y : double) : double; external;

  { Matrix routines }

  PROCEDURE UnitMat(n : integer; VAR a : matrix); external; 

  PROCEDURE CopyVec(n : integer; VAR a, b : vector); external;

  PROCEDURE CopyMat(n : integer; VAR a, b : matrix); external;

  PROCEDURE AddVec(n : integer; VAR a, b : vector); external;

  PROCEDURE SubVec(n : integer; VAR a, b : vector); external;

  PROCEDURE AddMat(n : integer; VAR a, b : matrix); external;

  PROCEDURE SubMat(n : integer; VAR a, b : matrix); external;

  PROCEDURE LinTrans(n : integer; VAR a : matrix; VAR x : vector); external;

  PROCEDURE MulLMat(n : integer; VAR a, b : matrix); external;

  PROCEDURE MulRMat(n : integer; VAR a, b : matrix); external;

  FUNCTION TrMat(n : integer; VAR a : matrix) : double; external;

  PROCEDURE TpMat(n : integer; VAR a : matrix); external;

  FUNCTION DetMat(n : integer; VAR a : matrix) : double; external;

  function InvMat(n : integer; VAR a : matrix) : boolean; external;

  procedure prtmat(n : integer; var a : matrix); external;
