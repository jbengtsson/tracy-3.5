  { Interface }

  const	nvmax=6; nomax=1; DAnamlen=10;

  type	DAvect		= array [0..nvmax] of double;
	ivector		= array [1..nvmax] of integer;
	DAnambuf	= packed array [1..DAnamlen] of char;
	funnambuf	= packed array [1..4] of char;
	DAmap		= array [1..nvmax] of DAvect;

  var	DAdim, DAord, ndim2	: [external] integer;
	DAeps1			: [external] double;


  procedure DAeps(eps : double); external;

  procedure DAini(o, n, u : integer); external;

  procedure lieini(no, nv, nd2i : integer); external;

  procedure DAall(var x : DAvect; nd2 : integer;
                  DAname : DAnambuf; no, nv : integer); external;

  procedure DAdal(var x : Davect; DAdim : integer); external;

  procedure DAmapall(var x : DAmap; nd2 : integer;
                     DAname : DAnambuf; no, nv : integer); external;

  procedure DAmapdal(var x : DAmap; nd2 : integer); external;

  procedure DAvar(var x : DAvect; r : double; i : integer); external;

  procedure DAcon(var x : DAvect; r : double); external;

  procedure ETini(var map : DAmap); external;

  procedure DApek(var x : DAvect; var jj : ivector; var r : double); external;

  procedure DApok(var x : DAvect; var jj : ivector; r : double); external;

  function getmat(var map : DAmap; i, j : integer) : double; external;

  procedure putmat(var map : DAmap; i, j : integer; r : double); external;

  procedure getlinmat(nv : integer; var map : DAmap; var mat : matrix); external;

  procedure putlinmat(nv : integer; var mat : matrix; var map : DAmap); external;

  procedure DAcop(var x, z : DAvect); external;

  procedure CopyMap(var map1, map2 : DAmap); external;

  procedure DAadd(var x, y, z : DAvect); external;

  procedure DAsub(var x, y, z : DAvect); external;

  procedure DAmul(var x, y, z : DAvect); external;

  procedure DAdiv(var x, y, z : DAvect); external;

  procedure DAcad(var x : DAvect; y : double; var z : DAvect); external;

  procedure DAcsu(var x : DAvect; y : double; var z : DAvect); external;

  procedure DAcmu(var x : DAvect; y : double; var z : DAvect); external;

  procedure DAsuc(var x : DAvect; y : double; var z : DAvect); external;

  procedure DAcdi(var x : DAvect; y : double; var z : DAvect); external;

  procedure DAdic(var x : DAvect; y : double; var z : DAvect); external;

  procedure DApos(var x, z : DAvect); external;

  procedure DAsqr(var x, z : DAvect); external;

  procedure DAcma(var x, y : DAvect; rb : double; var z : DAvect); external;

  procedure DAlin(var x : DAvect; ra : double; var y : DAvect; rb : double;
		  var z : DAvect); external;

  procedure DAfun(fun : funnambuf; var x, z : DAvect); external;

  procedure dacct(var x : DAmap; i : integer; var y : DAmap; j : integer;
		  var z : DAmap; k : integer); external;

  procedure dainv(var x : DAmap; i : integer; var z : DAmap; k : integer); external;

  procedure Rotmap(n : integer; var map : DAmap; var R : matrix); external;

  procedure DAwrite(var x : DAvect; var mapfil : text); external;

  procedure prtmap(var map : DAmap; var mapfil : text); external;

  procedure prtdamat(n : integer; var map : DAmap); external;
