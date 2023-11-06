  module t2elem(input, output);

  %include 'mathlib.def'
  %include 'dab.def'
  %include 'pascommon.def'
  %include 't2common.def'
  %include 'B2perp.def'
  %include 'etwigg.def'
  %include 'gfwigg.def'
  %include 't2elemo.def'

  { interface }

  var	ElemFam					: [global] Array [1..Elem_nFamMax] of ElemFamType;
	Cell					: [global] Array [0..Cell_nLocmax] of CellType; 
	Fdrift1, Fkick1, Fdrift2, Fkick2	: [global] double;
	crad, cfluc				: [global] double;
 

  [global] function GetnKid(Fnum1 : integer) : integer; forward;

  [global] function Elem_GetPos(Fnum1, Knum1 : integer) : integer; forward;

  [global] procedure getelem(i : integer; var cellrec : celltype); forward;

  [global] procedure putelem(i : integer; var cellrec : celltype); forward;

  [global] procedure GtoL(var X : vector; var S, T : vector2;
			  c0, c1, s1 : double); forward;
  { from Global to Local coordinate }

  [global] procedure GtoL_M(var X : matrix; var T : vector2); forward;
  { from Global to Local coordinates }

  [global] procedure DAgtoL(var map : DAmap; var S, T : vector2;
			    c0, c1, s1 : double); forward;
  { from Global to Local coordinate }

  [global] procedure LtoG(var X : vector; var S, T : vector2;
			  c0, c1, s1 : double); forward;
  { from Local to Global coordinate }

  [global] procedure LtoG_M(var X : matrix; var T : vector2); forward;
  { from Local to Global coordinate }

  [global] procedure DALtoG(var map : DAmap; var S, T : vector2;
			    c0, c1, s1 : double); forward;
  { from Local to Global coordinate }

  [global] function GetNofDBN(elem : integer) : integer; forward;
  { Returns the number of DBName defined for the elem-th element family }
 
  [global] procedure GetDBN(elem, kid : integer; var Name : DBNameType); forward;
  { Returns the kid-th database name of the elem-th element family }

  [global] procedure InitFkick; forward;
  { Init constants for 4-th order integrator }


  [global] procedure Elem_Print(var f : text; Fnum1 : integer); forward;
  
  { If Fnum1=0 then print all the elements }

  [global] function Elem_GetKval(Fnum1, Knum1, Order : integer) : double; forward;

  [global] procedure Elem_GetBend(Fnum1, Knum1 : integer;
				  var Tx, Tx1, Tx2 : double); forward;

  [global] function Elem_GetOrd(Fnum1, Knum1 : integer) : integer; forward;


  [global] procedure MulLsMat(var a, b : matrix); forward;

  [global] procedure LinsTrans(var a : matrix; var b : vector); forward;


  [global] procedure Drift_Alloc(var Elem : elemtype); forward;

  [global] procedure Drift_Print(var f : text; Fnum1 : integer); forward;

  [global] PROCEDURE Drift_SetMatrix(Fnum1, Knum1 : integer); forward;

  [global] Procedure Drift_Init(Fnum1 : integer); forward;

  [global] Procedure Drift_Pass(var Cell : CEllType; var x : vector); forward;

  [global] Procedure Drift_Pass_M(var Cell : CEllType; var xref : vector;
				  var x : matrix); forward;

  [global] Procedure Drift_DApass(var Cell : CEllType; var map : DAmap); forward;


  [global] procedure thinkick(Order : integer; var MB : mpolarray;
			      L, irho : double; pthick : pthicktype;
			      var x : vector); forward;

  [global] procedure thinkick_M(Order : integer; var MB : mpolarray;
			        L, irho : double; pthick : pthicktype;
			        var xref : vector;
			        var x : matrix); forward;

  [global] procedure DAthinkick(Order : integer; var MB : mpolarray;
			        L, irho : double; pthick : pthicktype;
			        var map : DAmap); forward;


  [global] procedure Mpole_Alloc(var Elem : elemtype); forward;

  [global] PROCEDURE Mpole_SetMatrix(Fnum1, Knum1 : integer; K : double); forward;

  [global] procedure Mpole_Init(Fnum1 : integer); forward;

  [global] PROCEDURE Mpole_SetIntegrator(Fnum1, Method, Segment : integer ); forward;

  [global] procedure Mpole_Print(var f : text; Fnum1 : integer); forward;

  [global] procedure Mpole_GetBend(Fnum1, Knum1 : integer;
				   var Tx, Tx1, Tx2 : double); forward;

  [global] PROCEDURE Mpole_SetPB(Fnum1, Knum1, Order : integer); forward;

  [global] function Mpole_GetPB(Fnum1, Knum1, Order : integer) : double; forward;

  [global] function Mpole_GetOrder(Fnum1, Knum1 : integer) : integer; forward;

  [global] procedure Mpole_GetPmeth(Fnum1, Knum1 : integer;
				    var Pmethodt, PNt : integer;
				    var Fdrift1t, Fkick1t, Fdrift2t, 
				    Fkick2t : double); forward;

  [global] PROCEDURE Mpole_SetdS(Fnum1, Knum1 : integer); forward;

  [global] PROCEDURE Mpole_SetdT(Fnum1, Knum1 : integer); forward;


  [global] procedure Mpole_Pass(var Cell : CEllType; var x : vector); forward;

  [global] procedure Mpole_Pass_M(var Cell : CEllType; var xref : vector;
				  var x : matrix); forward;

  [global] procedure Mpole_DApass(var Cell : CEllType; var map : DAmap); forward;


  [global] procedure Wiggler_Alloc(var Elem : elemtype); forward;

  [global] procedure Wiggler_Print(var f : text; Fnum1 : integer); forward;

  [global] PROCEDURE Wiggler_SetMatrix(Fnum1, Knum1 : integer;
				       L, lambda, k0, kx : double); forward;

  [global] PROCEDURE Wiggler_SetPB(Fnum1, Knum1, Order : integer); forward;

  [global] PROCEDURE Wiggler_Init(Fnum1 : integer); forward ;

  [global] procedure wiggler_Pass(var Cell : CEllType; var X : vector); forward;

  [global] procedure wiggler_Pass_M(var Cell : CEllType; var xref : vector;
				    var x : matrix); forward;

  [global] procedure wiggler_DApass(var Cell : CEllType; var map : DAmap); forward;


  [global] PROCEDURE Cav_Alloc(var Elem:elemtype); forward;

  [global] PROCEDURE Cav_Set(Fnum1, Knum1 : integer; Freq, Volt : double;
			     h : integer); forward;

  [global] PROCEDURE Cav_Get(Fnum1, Knum1 : integer; var Freq, Volt : double;
			     var h : integer); forward;

  [global] Procedure Cav_Init(Fnum1 : integer); forward;

  [global] PROCEDURE Cav_Pass(var Cell : CEllType; var X : vector); forward;

  [global] PROCEDURE Cav_DApass(var Cell : CEllType; var map : DAmap); forward;

  [global] procedure Marker_Init(Fnum1 : integer); forward;

  [global] PROCEDURE Marker_Pass(var Cell : CEllType; var X : vector); forward;

  { implementation }

  function GetnKid{Fnum1 : integer}{ : integer};
  begin
    GetnKid := ElemFam[Fnum1].nKid;
  end;

  function Elem_GetPos{Fnum1, Knum1 : integer}{ : integer};
  begin
    Elem_GetPos := ElemFam[Fnum1].KidList[Knum1];
  end;

  procedure getelem{i : integer; var cellrec : celltype};
  begin
    cellrec := cell[i];
  end;

  procedure putelem{i : integer; var cellrec : celltype};
  begin
    cell[i] := cellrec;
  end;

  procedure GtoL{var X : vector; var S, T : vector2; c0, c1, s1, : double};
  var	x1	: vector;
  begin
    { Simplified rotated prot }
    x[2] := x[2] + c1; x[4] := x[4] + s1;
    { Translate }
    x[1] := x[1] - S[1]; x[3] := x[3] - S[2];
    { Rotate }
    copyvec(4, x, x1);
    x[1] := T[1]*x1[1] + T[2]*x1[3]; x[2] := T[1]*x1[2] + T[2]*x1[4];
    x[3] := -T[2]*x1[1] + T[1]*x1[3]; x[4] := -T[2]*x1[2] + T[1]*x1[4];
    { Simplified prot }
    x[2] := x[2] - c0;
  end;

  procedure GtoL_M{var X : matrix; var T : vector2};
  var	R	: matrix;
  begin
    { Rotate }
    R[1, 1] := T[1]; R[1, 2] := 0.0; R[1, 3] := T[2]; R[1, 4] := 0.0;
    R[2, 1] := 0.0; R[2, 2] := T[1]; R[2, 3] := 0.0; R[2, 4] := T[2];
    R[3, 1] := -T[2]; R[3, 2] := 0.0; R[3, 3] := T[1]; R[3, 4] := 0.0;
    R[4, 1] := 0.0; R[4, 2] := -T[2]; R[4, 3] := 0.0; R[4, 4] := T[1];
    MulLmat(4, R, X);
  end;

  procedure DAGtoL{var map : DAmap; var S, T : vector2; c0, c1, s1 : double};
  var	map1		: DAmap;
  begin
    damapall(map1, 6, 'map1      ', 6, 6);
    { Simplified rotated prot }
    DAcad(map[2], c1, map[2]);
    DAcad(map[4], s1, map[4]);
    { Translate }
    DAcsu(map[1], S[1], map[1]);
    DAcsu(map[3], S[2], map[3]);
    { Rotate }
    Copymap(map, map1);
    DAlin(map1[1], T[1], map1[3], T[2], map[1]);
    DAlin(map1[2], T[1], map1[4], T[2], map[2]);
    DAlin(map1[1], -T[2], map1[3], T[1], map[3]);
    DAlin(map1[2], -T[2], map1[4], T[1], map[4]);
    { Simplified prot }
    DAcsu(map[2], c0, map[2]);
    damapdal(map1, 6);
  end;

  procedure LtoG{var X : vector; var S, T : vector2; c0, c1, s1 : double};
  var	x1	: vector;
  begin
    { Simplified prot }
    x[2] := x[2] - c0;
    { Rotate }
    copyvec(4, x, x1);
    x[1] := T[1]*x1[1] - T[2]*x1[3]; x[2] := T[1]*x1[2] - T[2]*x1[4];
    x[3] := T[2]*x1[1] + T[1]*x1[3]; x[4] := T[2]*x1[2] + T[1]*x1[4];
    { Translate }
    x[1] := x[1] + S[1]; x[3] := x[3] + S[2];
    { prot rotated }
    x[2] := x[2] + c1; x[4] := x[4] + s1;
  end;

  procedure LtoG_M{var X : matrix; var T : vector2};
  var	R	: matrix;
  begin
    { Rotate }
    R[1, 1] := T[1]; R[1, 2] := 0.0; R[1, 3] := -T[2]; R[1, 4] := 0.0;
    R[2, 1] := 0.0; R[2, 2] := T[1]; R[2, 3] := 0.0; R[2, 4] := -T[2];
    R[3, 1] := T[2]; R[3, 2] := 0.0; R[3, 3] := T[1]; R[3, 4] := 0.0;
    R[4, 1] := 0.0; R[4, 2] := T[2]; R[4, 3] := 0.0; R[4, 4] := T[1];
    MulLmat(4, R, X);
  end;

  procedure DALtoG{var map : DAmap; var S, T : vector2; c0, c1, s1 : double};
  var	map1	: DAmap;
  begin
    damapall(map1, 6, 'map1      ', 6, 6);
    { Simplified prot }
    DAcsu(map[2], c0, map[2]);
    { Rotate }
    Copymap(map, map1);
    DAlin(map1[1], T[1], map1[3], -T[2], map[1]);
    DAlin(map1[2], T[1], map1[4], -T[2], map[2]);
    DAlin(map1[1], T[2], map1[3], T[1], map[3]);
    DAlin(map1[2], T[2], map1[4], T[1], map[4]);
    { Translate }
    DAcad(map[1], S[1], map[1]);
    DAcad(map[3], S[2], map[3]);
    { prot rotated }
    DAcad(map[2], c1, map[2]);
    DAcad(map[4], s1, map[4]);
    damapdal(map1, 6);
  end;

  function GetNofDBN{elem : integer}{ : integer};
  begin
    GetNofDBN := ElemFam[elem].NoDBN;
  end;

  procedure GetDBN{elem, kid : integer; var Name : DANameType};
  var i	: integer;
  begin
    with ElemFam[elem] do 
    begin
      if kid > NoDBN then 
      begin
        for i:=1 to DBNameLen do
	  Name[i] := ' ';
      end
      else
	Name := DBNlist[kid];
    end;
  end;

  procedure InitFkick;

  function thirdroot(a : double) : double;
  { By substitution method }
  var	i	: integer;
	x	: double;
  begin
    x := 1.0; i := 0;
    repeat
      i := succ(i);
      x := (x+a)/(sqr(x)+1d0);
    until i = 250;
    thirdroot := x;
  end;

  begin
    Fdrift1 :=  1d0/(2d0*(2d0-thirdroot(2d0))); Fdrift2 :=  0.5d0 - Fdrift1;
    Fkick1  :=  2d0*Fdrift1; Fkick2  :=  1d0 - 2d0*Fkick1;
  end;

  procedure Elem_Print{f : text; Fnum1 : integer};
  var	i	: integer;
  begin
    if Fnum1 = 0 then
      for i:=1 to globval.Elem_NFam do
	Elem_Print(f, i)
    else
      case ElemFam[Fnum1].ElemF.Pkind of 
        drift: Drift_Print(f, Fnum1);
        mpole: Mpole_Print(f, Fnum1);
        wigl : Wiggler_Print(f, Fnum1);
        otherwise { nothing } ;
      end;
    end;

  function Elem_GetKval{Fnum1, Knum1, Order : integer) : double};
  begin
    with ElemFam[Fnum1].ElemF do
    begin
      case Pkind of 
        drift: Elem_GetKval := 0.0;
        mpole: begin
		 if M^.Pthick = thick then
		   Elem_GetKval := PL*Mpole_GetPB(Fnum1, Knum1, Order)
		 else
		   Elem_GetKval := Mpole_GetPB(Fnum1, Knum1, Order);
	       end;
        wigl : Elem_GetKval := PL
		 *sqrt(2.0*Cell[ElemFam[Fnum1].KidList[Knum1]].Elem.W^.PBW[Order]);
        otherwise { nothing } ;
      end;
    end;
  end;

  procedure Elem_GetBend{Fnum1, Knum1 : integer; var Tx, Tx1, Tx2 : double};
  begin
    case ElemFam[Fnum1].Elemf.Pkind of 
      drift: begin
	       Tx := 0.0; Tx1 := 0.0; Tx2 := 0.0;
	     end;
      mpole: Mpole_GetBend(Fnum1, Knum1, Tx, Tx1, Tx2);
      wigl : begin
	       Tx := 0.0; Tx1 := 0.0; Tx2 := 0.0;
	     end;
      otherwise { nothing } ;
    end;
  end;

  function Elem_GetOrd{Fnum1, Knum1 : integer) : integer};
  begin
    case ElemFam[Fnum1].Elemf.Pkind of 
      drift: Elem_GetOrd := 0;
      mpole: with Cell[ElemFam[Fnum1].KidList[Knum1]].Elem.M^ do
	       Elem_GetOrd := Porder;
      wigl : Elem_GetOrd := 0;
      otherwise { nothing } ;
    end;
  end;

  procedure make3by3(var A : matrix; a11, a12, a13, 
				     a21, a22, a23, 
				     a31, a32, a33 : double);
  begin
    UnitMat(matdim, A);
    a[1, 1] := a11; a[1, 2] := a12; a[1, 3] := a13;
    a[2, 1] := a21; a[2, 2] := a22; a[2, 3] := a23;
    a[3, 1] := a31; a[3, 2] := a32; a[3, 3] := a33;
  end;

  procedure make4by5(var A : matrix; a11, a12, a15, 
				     a21, a22, a25, 
				     a33, a34, a35, 
				     a43, a44, a45 : double);
  begin
    UnitMat(matdim, A);
    a[1, 1] := a11; a[1, 2] := a12; a[1, 5] := a15;
    a[2, 1] := a21; a[2, 2] := a22; a[2, 5] := a25;
    a[3, 3] := a33; a[3, 4] := a34; a[3, 5] := a35;
    a[4, 3] := a43; a[4, 4] := a44; a[4, 5] := a45;
  end;

  procedure mergeto4by5(var A, AH, AV : matrix);
  var	i, j	: integer;
  begin
    UnitMat(matdim, A);
    for i:=1 to 2 do
    begin
      A[i, 5] := AH[i, 3];
      A[i+2, 5] := AV[i, 3];
      for j:=1 to 2 do 
      begin 
        A[i, j] := AH[i, j];
        A[i+2, j+2] := AV[i, j];
      end;
    end;
  end;

  procedure LinsTrans{var a : matrix; var b : vector};
  const	n=4;
  VAR	j	: integer;
	c	: vector;
  begin
    CopyVec(n, b, c);
    LinTrans(n, a, c);
    for j:=1 to n do
      c[j] := c[j] + a[j, n+1]*b[n+1] + a[n+1, j];
    CopyVec(n, c, b);
  end;

  procedure MulLsMat{var a, b : matrix};
  const	n=4;
  VAR	i, k	: integer;
	c	: matrix;
  begin
    CopyMat(n, b, c);
    MulLMat(n, a, c);
    for i:=1 to n do
    begin
      c[i, n+1] := a[i, n+1]; c[n+1, i] := 0.0;
      for k:=1 to n do
      begin
        c[i, n+1] := c[i, n+1] + a[i, k]*b[k, n+1];
        c[n+1, i] := c[n+1, i] + a[i, k]*b[n+1, k];
      end;
    end;
    c[n+1, n+1] := 1.0;
    CopyMat(n+1, c, b);
  end;

  procedure Drift_Alloc;
  begin
    New(Elem.D);
  end;

  procedure Drift_Print{Fnum1 : integer};
  begin
    with ElemFam[Fnum1], ElemF do
    begin
      writeln(f, 'Element[', Fnum1:3, ' ] ');
      writeln(f, '   Name: ', pName, ',  Kind:   drift,  L=', PL:15);
      writeln(f, '   nKid:', nKid:3);
      writeln(f);
    end;
  end;

  PROCEDURE Drift_SetMatrix{Fnum1, Knum1 : integer};
  var	L	: double;
  BEGIN
    if ElemFam[Fnum1].nKid > 0 then
      with Cell[ ElemFam[Fnum1].Kidlist[Knum1] ], Elem, D^ do
      begin
        L := PL/(1+globval.dPparticle);
        make4by5(D55, 1, l, 0, 
                      0, 1, 0, 
                      1, l, 0, 
                      0, 1, 0);
      end;
  END;

  procedure Drift_Init{Fnum1 : integer};
  var	i	: integer;
  begin
    with ElemFam[Fnum1] do
    begin
      for i:=1 to nKid do
        with Cell[KidList[i]] do
        begin
	  Drift_Alloc(Elem);
	  Elem.Pname := Elemf.Pname;
	  Elem.PL := Elemf.PL;
	  Elem.Pkind := Elemf.Pkind;
	  Elem.D^ := Elemf.D^;
	  dT[1] := 1e0; dT[2] := 0.0; dS[1] := 0.0; dS[2] := 0.0;
          Drift_SetMatrix(Fnum1, i);
        end;
    end;
  end;

  procedure Drft(L, d : double; var x : vector);
  {		        L
		d = ---------
		    1 + delta
  }
  begin
    x[1] := x[1] + x[2]*d; x[3] := x[3] + x[4]*d;
    x[6] := x[6] + d*(sqr(x[2])+sqr(x[4]))/(2*(1+x[5]));
    if globval.pathlength then x[6] := x[6] + L;
  end;

  Procedure Drift_Pass{var Cell : CEllType; var X : vector};
  begin
    With Cell, Elem do
    begin
      Drft(PL, PL/(1+x[5]), x); CopyVec(6, x, BeamPos);
    end;
  end;

  Procedure Drift_Pass_M{var Cell : CEllType; var xref : vector;
			 var x : matrix};
  begin
    With Cell, Elem do
    begin
      MulLMat(5, D^.D55, x); Drft(PL, PL/(1+xref[5]), xref);
    end;
  end;

  Procedure DAdrft(L : double; d : davect; var map : DAmap);
  var	x1, x2, x3	: DAvect;
	map1		: DAmap;
  begin
    daall(x1, 1, 'x1        ', 6, 6);
    daall(x2, 1, 'x2        ', 6, 6);
    daall(x3, 1, 'x3        ', 6, 6);
    damapall(map1, 6, 'map1      ', 6, 6);
    CopyMap(map, map1);
    DAmul(map1[2], d, x1); DAadd(map1[1], x1, map[1]);
    DAmul(map1[4], d, x1); DAadd(map1[3], x1, map[3]);
    DAsqr(map1[2], x1); DAsqr(map1[4], x2); DAadd(x1, x2, x3);
    DAmul(x3, d, x2); DAcdi(x2, 2.0, x1);
    DAcad(map1[5], 1.0, x2); DAdiv(x1, x2, x3);
    DAadd(map1[6], x3, map[6]);
    if globval.pathlength then DAcad(map[6], L, map[6]);
    dadal(x1, 1); dadal(x2, 1); dadal(x3, 1); damapdal(map1, 6);
  end;

  Procedure Drift_DApass{var Cell : CEllType; var map : DAmap};
  var	x1, d1	: DAvect;
  begin
    daall(x1, 1, 'x1        ', 6, 6);
    daall(d1, 1, 'd1        ', 6, 6);
    DAcad(map[5], 1.0, x1); DAdic(x1, Cell.Elem.PL, d1);
    Dadrft(Cell.Elem.PL, d1, map);
    dadal(x1, 1); dadal(d1, 1);
  end;

  procedure thinkick{Order : integer; var MB : mpolarray; L, irho : double;
            	     pthick : pthicktype; var x : vector};
  { Calculate multipole kick. The kick is given by

                   e L      L delta      L x
	theta  = - --- B  + -------  -  -----  , 
             x     p    y     rho           2
                    0                    rho

                 e L
	theta  = --- B
             y   p    x
                  0

    where

                           ====
                           \                       n-1
	(B + iB  ) = B rho  >   (ia  + b ) (x + iy)
	  y    x           /       n    n
	                   ====

    where

			e      1
			-- = -----
			p    B rho
			 0
  }
  var	j					: integer;
	BxoBrho, ByoBrho, ByoBrho1, x1, x3, x5	: double;
	B2, xp, yp, psi, psf			: double;
	B					: vector3;
  begin
    x1 := x[1]; x3 := x[3]; x5 := x[5];
    if (1 <= Order) and (Order <= HOMmax) then
    begin
      ByoBrho := MB[order]; BxoBrho := MB[-order];
      for j:=order-1 downto 1 do
      begin
        ByoBrho1 := x1*ByoBrho - x3*BxoBrho + MB[ j];
        BxoBrho  := x3*ByoBrho + x1*BxoBrho + MB[-j];
        ByoBrho := ByoBrho1;
      end;
      if globval.radiation and (pthick = thick) then  
      begin 
	psi := 1d0 + x5; xp := x[2]/psi; yp := x[4]/psi;
	B[1] := BxoBrho; B[2] := ByoBrho + irho; B[3] := 0d0; 
        B2 := B2perp(irho, B, x1, xp, yp);
        x[5] := x[5] - crad*sqr(psi)*B2*(1d0+x1*irho+(sqr(xp)+sqr(yp))/2d0)*L;
	psf := 1d0 + x[5]; x[2] := xp*psf; x[4] := yp*psf;
      end;
      x[2] := x[2] - L*(ByoBrho-(x5-x1*irho)*irho); x[4] := x[4] + L*BxoBrho;
      x[6] := x[6] + L*irho*x1;
    end;
  end;

  procedure thinkick_M{Order : integer; var MB : mpolarray; L, irho : double;
            	       pthick : pthicktype; var xref : vector; var x : matrix};
  var	i	: integer;
	MMB	: MpolArray;
	z	: vector;
	Mk	: matrix;
  begin
    if (2 <= Order) and (Order <= HOMmax) then
    begin
      for i:=2 to Order do
      begin
        MMB[i-1] := (i-1)*MB[i]; MMB[-(i-1)] := (i-1)*MB[-i];
      end;
      z[1] := xref[1]; z[2] := 0.0; z[3] := xref[3]; z[4] := 0.0;
      z[5] := 0.0; z[6] := 0.0;
      thinkick(Order-1, MMB, L, 0.0, pthick, z);
      z[2] := z[2] - L*sqr(irho);
      Unitmat(5, Mk);
      Mk[2, 1] := z[2]; Mk[2, 3] := z[4]; Mk[4, 1] := z[4]; Mk[4, 3] := -z[2];
      MulLMat(5, Mk, x);
    end;
  end;

  procedure DAthinkick{Order : integer; var MB : mpolarray; L, irho : double;
            	       var map : DAmap};
  var	j					: integer;
        delta2, b30cf				: double;
	BxoBrho, ByoBrho, ByoBrho1		: DAvect;
	xn, x1, x2, x3, den, psi, psf, B2	: DAvect;
	map1, xp				: DAmap;
	B, e					: array [1..3] of DAvect;
  begin
    daall(x1, 1, 'x1        ', 6, 6); daall(x2, 1, 'x2        ', 6, 6);
    daall(x3, 1, 'x3        ', 6, 6);
    daall(BxoBrho, 1, 'BxoBrho ', 6, 6); daall(ByoBrho, 1, 'ByoBrho ', 6, 6);
    daall(ByoBrho1, 1, 'ByoBrho1', 6, 6);
    damapall(map1, 6, 'map1      ', 6, 6);
    CopyMap(map, map1);
    if (1 <= Order) and (Order <= HOMmax) then
    begin
      DAcon(ByoBrho, MB[order]); DAcon(BxoBrho, MB[-order]);
      for j:=order-1 downto 1 do
      begin
        DAmul(map1[1], ByoBrho, x1); DAmul(map1[3], BxoBrho, x2);
	DAsub(x1, x2, x3); DAcad(x3, MB[j], ByoBrho1);
        DAmul(map1[3], ByoBrho, x1); DAmul(map1[1], BxoBrho, x2);
	DAadd(x1, x2, x3); DAcad(x3, MB[-j], BxoBrho);
        DAcop(ByoBrho1, ByoBrho);
      end;
      if globval.radiation and (pthick = thick) then
      begin
	daall(xn, 1, 'xn        ', 6, 6);
        for j:=1 to 3 do
	begin
	  daall(e[j], 1, 'e         ', 6, 6);
	  daall(B[j], 1, 'B         ', 6, 6);
	end;
        daall(den, 1, 'den       ', 6, 6);
        damapall(xp, 6, 'xp        ', 6, 6);
        daall(psi, 1, 'psi       ', 6, 6); daall(psf, 1, 'psf       ', 6, 6);
        daall(B2, 1, 'B2        ', 6, 6);
        copymap(map, xp);
        dacad(map1[5], 1d0, psi);
        dadiv(map1[2], psi, xp[2]); dadiv(map1[4], psi, xp[4]);
        dacop(BxoBrho, B[1]); dacad(ByoBrho, irho, B[2]); dacon(B[3], 0d0);
	{ (1d0+x1*irho+0.5d0*(xp*xp+yp*yp))*L---- -> x3 }
	dasqr(xp[2], x1); dasqr(xp[4], x2); daadd(x1, x2, x3);
	dacmu(map1[1], irho, x2); dacad(x2, 1d0, x1); dasqr(x1, x2);
	daadd(x2, x3, x2); dafun('SQRT', x2, xn); dadic(xn, 1d0, xn);
	damul(xp[2], xn, e[1]); damul(xp[4], xn, e[2]); damul(x1, xn, e[3]);
	damul(B[2], e[3], x1); damul(B[3], e[2], x2); dasub(x1, x2, xn);
	dasqr(xn, B2);
	damul(B[2], e[1], x1); damul(B[1], e[2], x2); dasub(x1, x2, xn);
	dasqr(xn, x1);
	daadd(x1, B2, B2);
	damul(B[1], e[3], x1); damul(B[3], e[1], x2); dasub(x1, x2, xn);
	dasqr(xn, x1);
	daadd(x1, B2, B2);
        dacmu(x3, 0.5d0, x1); dacma(x1, map1[1], irho, x2); dacad(x2, 1d0, x1);
        dacmu(x1, L, x3);
	{ end of (1d0+x1*irho+0.5d0*(xp*xp+yp*yp))*L---- -> x3 }
        dasqr(psi, x1); damul(x1, x3, x2); damul(B2, x2, x1);
        dacmu(x1, crad, den); dasub(map1[5], den, map[5]);
        dacad(map[5], 1d0, psf);
	damul(xp[2], psf, map[2]); damul(xp[4], psf, map[4]);
        if globval.emittance then
        begin
          if B2[0] > 0d0 then 
          begin
            { synchrotron integrals }
            b30cf := pwr(B2[0], 1.5d0)*cfluc;
            dasqr(psi, psi); dasqr(psi, psi);
            damul(x3, psi, x1); dacmu(x1, b30cf, x2); delta2 := x2[0];
            dainv(xp, 6, xp, 6);
            for j:=1 to 3 do
              globval.qfluct[j] := globval.qfluct[j]
		 + ( pwr(xp[2*j-1, 5], 2d0)+pwr(xp[2*j, 5], 2d0) )*delta2;
          end;
        end;
        dadal(xn, 1);
        for j:=1 to 3 do
	begin
          dadal(e[j], 1); dadal(B[j], 1);
	end;
        dadal(B2, 1); dadal(den, 1); damapdal(xp, 6);
	dadal(psi, 1); dadal(psf, 1);
      end;      
      DAcma(map1[5], map1[1], -irho, x1); DAcma(ByoBrho, x1, -irho, x2);
      DAcma(map[2], x2, -L, map[2]); DAcma(map[4], BxoBrho, L, map[4]);
      DAcma(map1[6], map1[1], L*irho, map[6]);
    end;
    dadal(x1, 1); dadal(x2, 1); dadal(x3, 1);
    dadal(BxoBrho, 1); dadal(ByoBrho, 1); dadal(ByoBrho1, 1);
    damapdal(map1, 6);
  end;

  procedure Mpole_Alloc;
  var	j	: integer;
  begin
    New(Elem.M);
    with Elem, M^ do
    begin
      Pmethod := Meth_Linear; pn := 0;
      for j:=1 to 2 do
      begin
        PdSsys[j] := 0.0; PdSrnd[j] := 0.0;
      end;
      PdTpar := 0.0; PdTsys := 0.0; PdTrnd := 0.0;
      for j:=-HOMmax to HOMmax do
      begin
        PB[j] := 0.0; PBpar[j] := 0.0; PBsys[j] := 0.0;
	PBrms[j] := 0.0; PBrnd[j] := 0.0;
      end;
      Porder := 0; Pirho := 0.0; Ptx1 := 0.0; Ptx2 := 0.0; Pgap := 0.0;
    END;
  end;

  PROCEDURE driftmat(VAR ah : Matrix; L : double);
  BEGIN
    L := L/(1+globval.dPparticle);
    make4by5(ah, 1, L, 0, 
                 0, 1, 0, 
                 1, L, 0, 
                 0, 1, 0);
  END;

  PROCEDURE quadmat(VAR ahv : Matrix; L, k : double);
  {	output	Ah, Av:	transfer matrix
	input	L:	length [m]
		K:	gradient } 
  VAR	t, sk, sk0, s, c	: double;
	a, ah, av		: matrix;
  BEGIN
    if k <> 0.0 then
    begin{ focusing }
      sk0 := sqrt(abs(k));
      t := L*sk0/sqrt(1+globval.dPparticle); c := cos(t); s := sin(t);
      sk := sk0*sqrt(1+globval.dPparticle);
      make3by3(a,  c,   s/sk, 0, 
                 -sk*s,  c,   0, 
                   0,    0,   1);
      if k > 0.0 then
	CopyMat(3, a, ah)
      else
	CopyMat(3, a, av);
      c := cosh_(t); s := sinh_(t); sk := sk0*sqrt(1+globval.dPparticle);
      make3by3(a, c,   s/sk, 0, 
                 sk*s,  c,   0, 
                  0,    0,   1);
      if k > 0.0 then
	CopyMat(3, a, av)
      else
	CopyMat(3, a, ah);
      mergeto4by5(AHV, AH, AV);      
    end
    else
      driftmat(AHV, L);        
  END;

  function psi(irho, phi, gap : double) : double;
  { Correction for magnet gap }
  const	K1 = 0.5; K2 = 0.0;
  begin
    psi := K1*gap*irho*(1+sqr(sin(dtor(phi))))/cos(dtor(phi))
        *(1-K1*K2*gap*irho*tan_(dtor(phi)));
  end;

  PROCEDURE bendmat(VAR M : Matrix; L, irho, phi1, phi2, gap, k : double);
  {	output	Ah:	transfer matrix
	input	L:	length [m]
		irho:	1/rho [1/m]
		phi1:	entrance edge angle [degres]
		phi2:	exit edge angle [degres]
		K:	gradient = n/Rho
  } 
  VAR	r, s, c, sk, p, fk, afk	: double;
	edge, ah, av		: Matrix;
	coef, scoef		: double;
  BEGIN
    if irho = 0.0 then 
      quadmat(M, L, k)
    else
    BEGIN
      coef := 1+globval.dPparticle; scoef := sqrt(coef); r := L*irho/scoef;
      IF k = 0.0 THEN
      BEGIN
	c := cos(r); s := sin(r);
        make3by3(ah,    c,          s/(irho*scoef), (1-c)/irho, 
                    -s*scoef*irho,      c,          s*scoef, 
                        0,              0,             1);
        make3by3(av, 1, L/coef, 0, 
                     0,   1,    0, 
                     0,   0,    1);
      END
      else
      BEGIN { gradient bend, k= n/rho^2 }
	fk := -k - sqr(irho); afk := abs(fk); sk := sqrt(afk); p := L*sk/scoef;
	if fk < 0.0 then 
        BEGIN { horizontally focusing }
          c := cos(p); s := sin(p);
	  make3by3(ah,    c,        s/sk/scoef, irho*(1.0-c)/(coef*afk), 
		      -scoef*sk*s,      c,           -scoef*sk*s, 
		          0,            0,                1);
           sk := sqrt(abs(k)); p := L*sk/scoef; c := cosh_(p); s := sinh_(p);
	   make3by3(av,   c,        s/sk/scoef, 0, 
		       sk*s*scoef,     c,      0, 
		          0,           0,      1);
	END
	else 
        BEGIN
	  c := cosh_(p); s := sinh_(p);
	  make3by3(ah,   c,        s/sk/scoef,  (c-1)*irho/afk, 
	              scoef*s*sk,      c,      scoef*s*irho/sk, 
		         0,            0,            1);
          sk := sqrt(abs(k)); p := L*sk/scoef; c := cos(p);s := sin(p);
	  make3by3(av,     c,       s/sk/scoef, 0, 
		      -sk*s*scoef,      c,      0, 
		           0,           0,      1)
	END;
      END;
      { Edge focusing, no effect due to gap between AU and AD }
      if (phi1 <> 0.0) or (gap > 0.0) then
      BEGIN
        UnitMat(3, edge);
	edge[2, 1] := irho*tan_(dtor(phi1));
        MulRMat(3, ah, edge);
	edge[2, 1] := -irho*tan_(dtor(phi1)-psi(irho, phi1, gap));
        MulRMat(3, av, edge);
      END
      else if (phi2 <> 0.0) or (gap < 0.0) then
      BEGIN
        UnitMat(3, edge);
	edge[2, 1] := irho*tan_(dtor(phi2));
        MulLMat(3, edge, ah);
	edge[2, 1] := -irho*tan_(dtor(phi2)-psi(irho, phi2, abs(gap)));
        MulLMat(3, edge, av);
      END;
      mergeto4by5(M, AH, AV);      
    END;
  END;

  function UpdatePorder(Elem : elemtype) : integer;
  var	i, order	: integer;
  begin
    with Elem, M^ do
    begin
      if Pirho <> 0.0 then
        order := 1
      else
        order := 0;
      for i := -HOMmax to HOMmax do
        if (PB[i] <> 0.0) and (abs(i) > order) then order := abs(i);
    end;
    UpdatePorder := order;
  end;

  PROCEDURE Mpole_SetMatrix{Fnum1, Knum1 : integer; K : double};
  begin
    if ElemFam[Fnum1].nKid > 0 then
      with Cell[ ElemFam[Fnum1].Kidlist[Knum1] ], Elem, M^ do
      begin
	bendmat(au55, PL/2d0, Pirho, Ptx1, 0.0, Pgap, k);
	bendmat(ad55, PL/2d0, Pirho, 0.0, Ptx2, -Pgap, k);
      end;
  end;

  procedure Mpole_Init{Fnum1 : integer};
  var	x	: double;
	i	: integer;
  begin
    with ElemFam[Fnum1] do
    begin
      ElemF.M^.PB := ElemF.M^.PBpar;
      ElemF.M^.Porder := UpdatePorder(ElemF);
      x := ElemF.M^.PBpar[quad];  
      for i:=1 to nKid do
        with Cell[KidList[i]] do
        begin
	  Mpole_Alloc(Elem);
	  Elem.Pname := Elemf.Pname;
	  Elem.PL := Elemf.PL;
	  Elem.Pkind := Elemf.Pkind;
	  Elem.M^ := Elemf.M^;
	  dT[1] := 1e0; dT[2] := 0.0; dS[1] := 0.0; dS[2] := 0.0;
          with Elem do
          if (PL <> 0.0) or (M^.Pirho <> 0.0) then
	  begin
	    M^.Pthick := Thick; M^.Pc0 := 0.0; M^.Pc1 := 0.0; M^.Ps1 := 0.0;
            Mpole_SetMatrix(Fnum1, i, x);
	  end
          else
	    M^.Pthick := Thin;
        end;
    end;
  end;

  PROCEDURE Mpole_SetIntegrator{Fnum1 : integer; Method : integer;
				Segment : integer};
  var	i	: integer ;
  begin{1}
    with ElemFam[Fnum1] do
    begin 
      ElemF.M^.Pmethod := Method;  ElemF.M^.PN := Segment; 
      for i:=1 to nKid do
	with Cell[ KidList[i] ].Elem.M^ do
	begin
	  Pmethod := Method; PN := Segment; 
	  if Method = Meth_linear then Mpole_SetMatrix(Fnum1, i, PB[quad]);
	end;
    end;
  end;{1}

  procedure Mpole_Print{Fnum1 : integer};
  begin
    with ElemFam[Fnum1].ElemF, M^ do
    begin
      writeln(f, 'Element[', Fnum1:3, ' ] ');
      writeln(f, '   Name: ', pName, ',  Kind:   mpole,  L=', PL:15);
      writeln(f, '   Method: ', Pmethod:0, ', N=', PN:4);
    end;
  end;

  procedure Mpole_GetBend{Fnum1, Knum1 : integer; var Tx, Tx1, Tx2 : double};
  begin
    with Cell[ElemFam[Fnum1].KidList[Knum1]].Elem, M^ do
      if Pthick = thick then
      begin
        Tx := PL*Pirho*180/pi; Tx1 := PTx1; Tx2 := PTx2;
      end
      else
      begin
        Tx := 0.0; Tx1 := 0.0; Tx2 := 0.0;
      end;
  end;

  PROCEDURE Mpole_SetPB{Fnum1, Knum1, Order : integer};
  begin
    with Cell[ ElemFam[Fnum1].KidList[Knum1] ], Elem, M^ do
    begin
      PB[Order] := PBpar[Order] + PBsys[Order] + PBrms[Order]*PBrnd[Order];
      if (abs(Order) > Porder) and (PB[Order] <> 0.0) then Porder := abs(Order);
      if (Pmethod = Meth_linear) and (Order = 2) then
        Mpole_SetMatrix(Fnum1, Knum1, PB[Order]);
    end;
    cellconcat := false;
  end;

  function Mpole_GetPB{Fnum1, Knum1, Order : integer}{ : double};
  begin
    with Cell[ElemFam[Fnum1].KidList[Knum1]].Elem.M^ do
      Mpole_GetPB := PB[Order];
  end;

  function Mpole_GetOrder{Fnum1, Knum1 : integer}{ : integer};
  begin
    with Cell[ElemFam[Fnum1].KidList[Knum1]].Elem.M^ do
      Mpole_GetOrder := Porder;
  end;

  procedure Mpole_GetPmeth{Fnum1, Knum1 : integer; var Pmethodt, PNt : integer;
			   var Fdrift1t, Fkick1t, Fdrift2t, Fkick2t : double};
  begin
    with Cell[ElemFam[Fnum1].KidList[Knum1]].Elem.M^ do
    begin
      Pmethodt := ord(Pmethod); PNt := PN;
    end;
    Fdrift1t := Fdrift1; Fkick1t := Fkick1;
    Fdrift2t := Fdrift2; Fkick2t := Fkick2;
  end;

  procedure Mpole_SetdS{Fnum1, Knum1 : integer};
  var	j	: integer;
  begin
    with Cell[ ElemFam[Fnum1].KidList[Knum1] ], Elem, M^ do
    begin
      for j:=1 to 2 do
	dS[j] := PdSsys[j] + PdSrms[j]*PdSrnd[j];
    end;
    cellconcat := false;
  end;

  procedure Mpole_SetdT{Fnum1, Knum1 : integer};
  begin
    with Cell[ ElemFam[Fnum1].KidList[Knum1] ], Elem, M^ do
    begin
      dT[1] := cos( dtor( PdTpar+PdTsys+PdTrms*PdTrnd ) );
      dT[2] := sin( dtor( PdTpar+PdTsys+PdTrms*PdTrnd ) );
      { Calculate simplified prots }
      Pc0 := sin( PL*Pirho/2d0 );
      Pc1 := cos( dtor(PdTpar) )*Pc0; Ps1 := sin( dtor(PdTpar) )*Pc0;
    end;
    cellconcat := false;
  end;

  Procedure EdgeFocus(irho, phi, gap : double; var x : vector);
  begin
    x[2] := x[2] + irho*tan_(dtor(phi))*x[1];
    x[4] := x[4] - irho*tan_(dtor(phi)-psi(irho, phi, gap))*x[3];
  end;

  procedure Mpole_Pass{var Cell : CEllType; var X  : vector};
  var	seg				: integer;
	k, L, L1, L2, D1, D2, K1, K2	: double;
  begin
    with Cell, Elem, M^ do
    begin
      { Global -> Local }
      GtoL(x, dS, dT, Pc0, Pc1, Ps1);
      case Pmethod of
        Meth_Linear, Meth_First:
	  begin { Tracy integrator }
	    if Pthick = thick then
	    begin
              { thick element }
              { First Linear  }
	      LinTrans(5, AU55, x);
	      k := PB[quad]; PB[quad] := 0;
	      { Kick }
	      thinkick(Porder, PB, PL, 0, pthick, x);
	      PB[quad] := k;
              { Second Linear }
	      LinTrans(5, AD55, x);
            end
	    else
	    begin
              { thin kick }
              { Kick }
	      thinkick(Porder, PB, 1, 0, pthick, x);
            end;
            { Save Beam Pos }
            CopyVec(5, x, BeamPos);
          end;
        Meth_Second:
	  begin { second order integrator }
            if Pthick = thick then
	    begin
              { thick element }
              if (PTx1 <> 0.0) or (Pgap <> 0.0) then
		EdgeFocus(Pirho, PTx1, Pgap, x);
              L := PL/PN; L1 := 0.5*L; D1 := 0.5*L/(1+x[5]);
              for seg := 1 to PN do
              begin 
                Drft(L1, D1, x);
                thinkick(Porder, PB, L, Pirho, pthick, x);
		if globval.radiation then D1 := 0.5*L/(1+x[5]);
                Drft(L1, D1, x);
              end;
              if (PTx2 <> 0.0) or (Pgap <> 0.0) then
		EdgeFocus(Pirho, PTx2, Pgap, x);
            end
	    else
	    begin
              { thin kick }
	      thinkick(Porder, PB, 1, 0, pthick, x);
            end;
            { Save Beam Pos }
            CopyVec(6, x, BeamPos);
          end;
        Meth_Fourth:
	  begin { 4-th order integrator }
            if Pthick = thick then
	    begin
              { thick element }
              if (PTx1 <> 0.0) or (Pgap <> 0.0) then
		EdgeFocus(Pirho, PTx1, Pgap, x);
              L  := PL/PN; L1 := Fdrift1*L; L2 := Fdrift2*L;
              D1 := Fdrift1*L/(1+x[5]); D2 := Fdrift2*L/(1+x[5]);
	      K1 := Fkick1*L; K2 := Fkick2*L;
              for seg := 1 to PN do
              begin 
                Drft(L1, D1, x);
                thinkick(Porder, PB, K1, Pirho, pthick, x);
		if globval.radiation then D2 := fdrift2*L/(1+x[5]);
                Drft(L2, D2, x);
                thinkick(Porder, PB, K2, Pirho, pthick, x);
		if globval.radiation then D2 := fdrift2*L/(1+x[5]);
                Drft(L2, D2, x);
                thinkick(Porder, PB, K1, Pirho, pthick, x);
		if globval.radiation then D1 := fdrift1*L/(1+x[5]);
                Drft(L1, D1, x);
              end;
              if (PTx2 <> 0.0) or (Pgap <> 0.0) then
		EdgeFocus(Pirho, PTx2, Pgap, x);
            end
	    else
	    begin
              { thin kick }
	      thinkick(Porder, PB, 1, 0, pthick, x);
            end;
            { Save Beam Pos }
            CopyVec(6, x, BeamPos);
          end;
        otherwise { Nothing }
      end;
      { Local -> Global }
      LtoG(x, dS, dT, Pc0, Pc1, Ps1);
    end;
  end;{1}

  procedure Mpole_Pass_M{var Cell : CEllType; var xref : vector;
			 var x : matrix};
  var	k	: double;
  begin
    with Cell, Elem, M^ do
    begin
      { Global -> Local }
      GtoL_M(x, dT); GtoL(xref, dS, dT, Pc0, Pc1, Ps1);
      case Pmethod of
        Meth_Linear, Meth_First:
	  begin { Tracy integrator }
	    if Pthick = thick then 
	    begin
	      { thick element }
              { First Linear }
	      MulLmat(5, AU55, x); LinTrans(5, AU55, xref);
	      k := PB[quad]; PB[quad] := 0;
	      { Kick }
	      thinkick_M(Porder, PB, PL, 0, pthick, xref, x);
	      thinkick(Porder, PB, PL, 0, pthick, xref);
	      PB[quad] := k;
              { Second Linear }
              MulLmat(5, AD55, x); LinTrans(5, AD55, xref);
	    end
	    else
	    begin
	      { thin kick }
	      thinkick_M(Porder, PB, 1, 0, pthick, xref, x);
	      thinkick(Porder, PB, 1, 0, pthick, xref);
	    end;
	  end;
        otherwise { Nothing }
      end;
      { Local -> Global }
      LtoG_M(x, dT); LtoG(xref, dS, dT, Pc0, Pc1, Ps1);
    end;
  end;

  Procedure DAedgeFocus(irho, phi, gap : double; var map : DAmap);
  var	map1	: DAmap;
  begin
    damapall(map1, 6, 'map1      ', 6, 6);
    CopyMap(map, map1);
    DAcma(map1[2], map1[1], irho*tan_(dtor(phi)), map[2]);
    DAcma(map1[4], map1[3], -irho*tan_(dtor(phi)-psi(irho, phi, gap)), 
	  map[4]);
    damapdal(map1, 6);
  end;

  procedure Mpole_DApass{var Cell : CEllType; var map : DAmap};
  var	seg			: integer;
	L, L1, L2, K1, K2	: double;
	D1, D2, x1		: DAvect;
  begin{1}
    daall(D1, 1, 'D1        ', 6, 6);
    daall(D2, 1, 'D2        ', 6, 6);
    daall(x1, 1, 'x1        ', 6, 6);
    with Cell, Elem, M^ do
    begin
      { Global -> Local }
      DAgtoL(map, dS, dT, Pc0, Pc1, Ps1);
      case Pmethod of
        Meth_Second:
	  begin { second order integrator }
            if Pthick = thick then
	    begin
              { thick element }
              if (PTx1 <> 0) or (Pgap <> 0.0) then
		DAedgeFocus(Pirho, PTx1, Pgap, map);
              L := PL/PN; L1 := 0.5*L;
	      DAcad(map[5], 1.0, x1); DAdic(x1, 0.5*L, D1);
              for seg := 1 to PN do
              begin 
                DAdrft(L1, D1, map);
                DAthinkick(Porder, PB, L, Pirho, pthick, map);
		if globval.radiation then
		begin
		  DAcad(map[5], 1.0, x1); DAdic(x1, 0.5*L, D1);
		end;
                DAdrft(L1, D1, map);
              end;
              if (PTx2 <> 0) or (Pgap <> 0.0) then
		DAedgeFocus(Pirho, PTx2, Pgap, map);
            end
	    else
	    begin
              { thin kick }
	      DAthinkick(Porder, PB, 1, 0, pthick, map);
            end;
          end;
        Meth_Fourth:
	  begin { 4-th order integrator }
            if Pthick = thick then
	    begin
              { thick element }
              if (PTx1 <> 0) or (Pgap <> 0.0) then
		DAEdgeFocus(Pirho, PTx1, Pgap, map);
              L  := PL/PN; L1 := Fdrift1*L; L2 := Fdrift2*L;
              DAcad(map[5], 1.0, x1);
	      DAdic(x1, Fdrift1*L, D1); K1 := Fkick1*L;
              DAdic(x1, Fdrift2*L, D2); K2 := Fkick2*L;
              for seg := 1 to PN do
              begin 
                DADrft(L1, D1, map);
                DAthinkick(Porder, PB, K1, Pirho, pthick, map);
		if globval.radiation then
		begin
		  DAcad(map[5], 1.0, x1); DAdic(x1, Fdrift2*L, D2);
		end;
                DADrft(L2, D2, map);
                DAthinkick(Porder, PB, K2, Pirho, pthick, map);
		if globval.radiation then
		begin
		  DAcad(map[5], 1.0, x1); DAdic(x1, Fdrift2*L, D2);
		end;
                DADrft(L2, D2, map);
                DAthinkick(Porder, PB, K1, Pirho, pthick, map);
		if globval.radiation then
		begin
		  DAcad(map[5], 1.0, x1); DAdic(x1, Fdrift1*L, D1);
		end;
                DADrft(L1, D1, map);
              end;
              if (PTx2 <> 0) or (Pgap <> 0.0) then
		DAEdgeFocus(Pirho, PTx2, Pgap, map);
            end
	    else
	    begin
              { thin kick }
	      DAthinkick(Porder, PB, 1, 0, pthick, map);
            end;
          end;
        otherwise { Nothing }
      end;
      { Local -> Global }
      DAltoG(map, dS, dT, Pc0, Pc1, Ps1);
    end;
    dadal(D1, 1); dadal(D2, 1); dadal(x1, 1);
  end;{1}

  procedure Wiggler_Alloc;
  var	j	: integer;
  begin
    New(Elem.W);
    with Elem, W^ do
    begin 
      Pmethod := Meth_Linear; PN := 0; PBoBrho := 0d0; Pkx := 0d0;
      for j:=1 to 2 do
      begin
        PdSsys[j] := 0.0; PdSrnd[j] := 0.0;
      end;
      PdTpar := 0.0; PdTsys := 0.0; PdTrnd := 0.0;
      for j:=0 to HOMmax do
      begin
        PBW [j] := 0.0;
      end;
      porder := 0;
    END;
  end;

  procedure Wiggler_Print{ Fnum1 : integer};
  begin
    with ElemFam[Fnum1].ElemF do
    begin
      writeln(f, 'Element[', Fnum1:3, ' ] ');
      writeln(f, '   Name: ', pName, ',  Kind:   wiggler,  L=', PL:15);
      writeln(f);
    end;
  end;

  PROCEDURE Wiggler_SetMatrix{Fnum1, Knum1 : integer;
			      L, lambda, k0, kx : double};
  VAR	t, s, c, k, ky, kz, LL	: double;
	ah, av			: matrix;
  BEGIN
    LL := L/(1+globval.dPparticle);
    if kx = 0d0 then
    begin
      make3by3(ah, 1, LL, 0, 
                   0, 1,  0, 
                   0, 0,  1);
    end
    else
    begin
      kz := 2d0*pi/lambda;
      k := sqrt(sqr(kx/kz)*k0);
      t := LL*k;
      c := cosh_(t);  s := sinh_(t);
      make3by3(ah,  c,  s/k, 0, 
                   k*s,  c,  0, 
                    0,   0,  1);
    end;
    if k0 = 0d0 then
    begin
      make3by3(av, 1, LL, 0, 
                   0, 1,  0, 
                   0, 0,  1);
    end
    else
    begin
      ky := sqrt(sqr(kx)+sqr(kz));
      k := sqrt(sqr(ky/kz)*k0);
      t := LL*k;
      c := cos(t);  s := sin(t);
      make3by3(av,  c,  s/k, 0, 
                  -k*s,  c,  0, 
                    0,   0,  1);
    end;
    with Cell[ ElemFam[Fnum1].Kidlist[Knum1] ].Elem.W^ do
    mergeto4by5(W55, AH, AV);      
  END;

  PROCEDURE Wiggler_SetPB{Fnum1, Knum1, Order : integer};
  begin
    with Cell[ ElemFam[Fnum1].KidList[Knum1] ], Elem, W^ do
    begin
      if abs(Order) > Porder then Porder := abs(Order);
      if (Pmethod = Meth_linear) and (Order = 2) then
        Wiggler_SetMatrix(Fnum1, Knum1, PL, Plperiod, PBW[Order], Pkx);
      cellconcat := false;
    end;
  end;

  PROCEDURE Wiggler_Init{Fnum1 : integer};
  const	order=2;
  var	i	: integer;
	x	: double;
  begin
    with ElemFam[Fnum1] do
    begin
{      ElemF.M^.PB := ElemF.M^.PBpar;}
      ElemF.W^.Porder := order;
      x := ElemF.W^.PBW[quad];
      for i:=1 to nKid do
	with Cell[KidList[i]] do
	begin
	  Wiggler_Alloc(Elem);
	  Elem.Pname := Elemf.Pname;
	  Elem.PL := Elemf.PL;
	  Elem.Pkind := Elemf.Pkind;
	  Elem.W^ := Elemf.W^;
	  dT[1] := 1; dT[2] := 0; dS[1] := 0; dS[2] := 0;
          Wiggler_SetMatrix(Fnum1, i, Elem.PL, Elem.W^.Plperiod, x,
			    Elem.W^.Pkx);
        end;
    end;
  end;  

  procedure wiggler_Pass{var Cell : CEllType; var X : vector};
  const	eps = 1d-18; kx = 0d0;
  var	seg, pthlen, radia		: integer;
	L, L1, L2, D1, D2, K1, K2	: double;
  begin
    with Cell, Elem, W^ do
    begin
      { Global -> Local }
      GtoL(x, dS, dT, 0.0, 0.0, 0.0);
      case Pmethod of
        Meth_Linear:
	  begin
	    LinTrans(5, W55, X); CopyVec(5, x, BeamPos);
	  end;
        Meth_First:
	  begin
	    if globval.pathlength then
	      pthlen := 1
	    else
	      pthlen := 0;
	    if globval.radiation then
	      radia := 1
	    else
	      radia := 0;
	    etwigg(PN, PL, Plperiod, PBoBrho, Pkx, x, crad, pthlen, radia);
	    CopyVec(6, x, BeamPos);
	  end;
        Meth_Second:
	  begin
	    { second order integrator }
            L := PL/PN; L1 := 0.5*L; D1 := 0.5*L/(1+x[5]);
            for seg := 1 to PN do
            begin 
	      Drft(L1, D1, x);
              thinkick(Porder, PBW, L, 0.0, thick, x);
              Drft(L1, D1, x);
            end;
            { Save Beam Pos }
            CopyVec(6, x, BeamPos);
          end;
        Meth_Fourth:
	  begin { 4-th order integrator }
            L  := PL/PN; L1 := Fdrift1*L; L2 := Fdrift2*L;
            D1 := Fdrift1*L/(1+x[5]); D2 := Fdrift2*L/(1+x[5]);
	    K1 := Fkick1*L; K2 := Fkick2*L;
            for seg := 1 to PN do
            begin 
              Drft(L1, D1, x);
              thinkick(Porder, PBW, K1, 0.0, thick, x);
              Drft(L2, D2, x);
              thinkick(Porder, PBW, K2, 0.0, thick, x);
              Drft(L2, D2, x);
              thinkick(Porder, PBW, K1, 0.0, thick, x);
              Drft(L1, D1, x);
            end;
            { Save Beam Pos }
            CopyVec(6, x, BeamPos);
          end;
        Meth_genfun:
	  begin
	    gfwigg(PN, eps, PL, Plperiod, PBoBrho, Pkx, x);
	    CopyVec(5, x, BeamPos);
	  end;
        otherwise { Nothing }
      end;
      { Local -> Global }
      LtoG(x, dS, dT, 0.0, 0.0, 0.0);
    end;
  end;

  procedure Wiggler_Pass_M{var Cell : CEllType; var xref : vector;
			   var x : matrix};
  begin
    with Cell, Elem, W^ do
    begin
      { Global -> Local }
      GtoL_M(x, dT); GtoL(xref, dS, dT, 0.0, 0.0, 0.0);
      case Pmethod of
        Meth_Linear:
	  begin { Tracy integrator }
	    MulLmat(5, W55, x); LinTrans(5, W55, xref);
	  end;
        otherwise { Nothing }
      end;
      { Local -> Global }
      LtoG_M(x, dT); LtoG(xref, dS, dT, 0.0, 0.0, 0.0);
    end;
  end;

  procedure wiggler_DApass{var Cell : CEllType; var map : DAmap};
  var	seg			: integer;
	L, L1, L2, K1, K2	: double;
	D1, D2, x1		: DAvect;
  begin
    daall(D1, 1, 'D1        ', 6, 6);
    daall(D2, 1, 'D2        ', 6, 6);
    daall(x1, 1, 'x1        ', 6, 6);
    with Cell, Elem, W^ do
    begin
      { Global -> Local }
      DAgtoL(map, dS, dT, 0.0, 0.0, 0.0);
      case Pmethod of
        Meth_Second:
	  begin { second order integrator }
            L := PL/PN; L1 := 0.5*L;
	    DAcad(map[5], 1.0, x1); DAdic(x1, 0.5*L, D1);
            for seg := 1 to PN do
            begin 
              DAdrft(L1, D1, map);
	      DAthinkick(Porder, PBW, L, 0.0, thick, map);
	      DAdrft(L1, D1, map);
            end;
          end;
        Meth_Fourth:
	  begin { 4-th order integrator }
            L := PL/PN; L1 := Fdrift1*L; L2 := Fdrift2*L;
            DAcad(map[5], 1.0, x1);
	    DAdic(x1, Fdrift1*L, D1); K1 := Fkick1*L;
            DAdic(x1, Fdrift2*L, D2); K2 := Fkick2*L;
            for seg := 1 to PN do
            begin 
              DADrft(L1, D1, map);
              DAthinkick(Porder, PBW, K1, 0.0, thick, map);
              DADrft(L2, D2, map);
              DAthinkick(Porder, PBW, K2, 0.0, thick, map);
              DADrft(L2, D2, map);
              DAthinkick(Porder, PBW, K1, 0.0, thick, map);
              DADrft(L1, D1, map);
            end;
          end;
        otherwise { Nothing }
      end;
      { Local -> Global }
      DAltoG(map, dS, dT, 0.0, 0.0, 0.0);
    end;
    dadal(D1, 1); dadal(D2, 1); dadal(x1, 1);
  end;

  procedure Cav_Alloc;
  begin
    New(Elem.C);
    with Elem, C^ do
    begin 
      Pvolt:=0.0; Pfreq:=0.0; Ph := 0;
    end;
  end;

  PROCEDURE Cav_Set{Fnum1, Knum1 : integer; Freq, Volt : double; h : integer};
  begin
    with Cell[ElemFam[Fnum1].KidList[Knum1]].Elem.C^ do
    begin
      Pfreq := Freq; Pvolt := Volt; Ph := h;
    end;
  end;

  PROCEDURE Cav_Get{Fnum1, Knum1 : integer; var Freq, Volt : double;
		    var h : integer};
  begin
    with Cell[ElemFam[Fnum1].KidList[Knum1]].Elem.C^ do
    begin
      Freq := Pfreq; Volt := Pvolt; h := Ph;
    end;
  end;

  procedure Cav_Init{Fnum1 : integer};
  var	i	: integer;
  begin
    with ElemFam[Fnum1], ElemF do
      for i:=1 to nKid do
	with Cell[KidList[i]] do
	  Elem := ElemF;
  end;

  procedure Cav_Pass{var Cell : CEllType; var X  : vector};
  const	c0 = 2.99792458d8;
  begin
    with Cell, Elem, C^ do
    begin
      if globval.Cavity_on and (Pvolt <> 0.0) then
      begin
        x[5] := x[5] - Pvolt/(GlobVal.Energy*1d9)*sin(2*Pi*Pfreq/c0*x[6]);
        if globval.pathlength then x[6] := x[6] - Ph/Pfreq*c0;
      end;
      { Save Beam Pos }
      CopyVec(6, x, BeamPos);
    end;
  end;

  procedure Cav_DApass{var Cell : CEllType; var map  : DAmap};
  const	c0 = 2.99792458d8;
  var	x1, x2	: Davect;
  begin
    daall(x1, 1, 'x1        ', 6, 6);
    daall(x2, 1, 'x2        ', 6, 6);
    with Cell, Elem, C^ do
    begin
      if globval.Cavity_on and (Pvolt <> 0.0) then
      begin
	DAcmu(map[6], 2*Pi*Pfreq/c0, x1);
	DAfun('SIN ', x1, x2);
	DAcmu(x2, -Pvolt/(GlobVal.Energy*1d9), x1);
        DAadd(map[5], x1, map[5]);
	if globval.radiation then globval.dE := globval.dE - x1[0];
        if globval.pathlength then DAcsu(map[6], Ph/Pfreq*c0, map[6]);
      end;
    end;
    dadal(x1, 1); dadal(x2, 1);
  end;

  procedure Marker_Init{Fnum1 : integer};
  var	i	: integer;
  begin
    with ElemFam[Fnum1], ElemF do
      for i:=1 to nKid do
	with Cell[KidList[i]] do
        begin
	  Elem := ElemF;
	  dT[1] := 1e0; dT[2] := 0.0; dS[1] := 0.0; dS[2] := 0.0;
        end;
  end;

  Procedure Marker_Pass{var Cell : CEllType; var X : vector};
  begin
    With Cell, Elem do
    begin
      { Global -> Local }
      GtoL(x, dS, dT, 0.0, 0.0, 0.0);
      CopyVec(6, x, BeamPos);
      { Local -> Global }
      LtoG(x, dS, dT, 0.0, 0.0, 0.0);
    end;
  end;

  end.
