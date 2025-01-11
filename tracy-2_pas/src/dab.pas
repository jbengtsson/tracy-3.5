  module dab(input, output);

  %include 'mathlib.def'
  %include 'dabo.def'

  { Interface }

  var	DAdim, DAord, ndim2	: [global] integer;
	DAeps1			: [global] double;


  [global] procedure DAeps(eps : double); forward;

  [global] procedure DAini(o, n, u : integer); forward;

  [global] procedure lieini(no, nv, nd2i : integer); forward;

  [global] procedure DAall(var x : DAvect; nd2 : integer;
                  DAname : DAnambuf; no, nv : integer); forward;

  [global] procedure DAdal(var x : Davect; DAdim : integer); forward;

  [global] procedure DAmapall(var x : DAmap; nd2 : integer;
                     DAname : DAnambuf; no, nv : integer); forward;

  [global] procedure DAmapdal(var x : DAmap; nd2 : integer); forward;

  [global] procedure DAvar(var x : DAvect; r : double; i : integer); forward;

  [global] procedure DAcon(var x : DAvect; r : double); forward;

  [global] procedure ETini(var map : DAmap); forward;

  [global] procedure DApek(var x : DAvect; var jj : ivector; var r : double); forward;

  [global] procedure DApok(var x : DAvect; var jj : ivector; r : double); forward;

  [global] function getmat(var map : DAmap; i, j : integer) : double; forward;

  [global] procedure putmat(var map : DAmap; i, j : integer; r : double); forward;

  [global] procedure getlinmat(nv : integer; var map : DAmap; var mat : matrix); forward;

  [global] procedure putlinmat(nv : integer; var mat : matrix; var map : DAmap); forward;

  [global] procedure DAcop(var x, z : DAvect); forward;

  [global] procedure CopyMap(var map1, map2 : DAmap); forward;

  [global] procedure DAadd(var x, y, z : DAvect); forward;

  [global] procedure DAsub(var x, y, z : DAvect); forward;

  [global] procedure DAmul(var x, y, z : DAvect); forward;

  [global] procedure DAdiv(var x, y, z : DAvect); forward;

  [global] procedure DAcad(var x : DAvect; y : double; var z : DAvect); forward;

  [global] procedure DAcsu(var x : DAvect; y : double; var z : DAvect); forward;

  [global] procedure DAcmu(var x : DAvect; y : double; var z : DAvect); forward;

  [global] procedure DAsuc(var x : DAvect; y : double; var z : DAvect); forward;

  [global] procedure DAcdi(var x : DAvect; y : double; var z : DAvect); forward;

  [global] procedure DAdic(var x : DAvect; y : double; var z : DAvect); forward;

  [global] procedure DApos(var x, z : DAvect); forward;

  [global] procedure DAsqr(var x, z : DAvect); forward;

  [global] procedure DAcma(var x, y : DAvect; rb : double; var z : DAvect); forward;

  [global] procedure DAlin(var x : DAvect; ra : double; var y : DAvect; rb : double;
		  var z : DAvect); forward;

  [global] procedure DAfun(fun : funnambuf; var x, z : DAvect); forward;

  [global] procedure dacct(var x : DAmap; i : integer; var y : DAmap; j : integer;
		  var z : DAmap; k : integer); forward;

  [global] procedure dainv(var x : DAmap; i : integer; var z : DAmap; k : integer); forward;

  [global] procedure Rotmap(n : integer; var map : DAmap; var R : matrix); forward;

  [global] procedure DAwrite(var x : DAvect; var mapfil : text); forward;

  [global] procedure prtmap(var map : DAmap; var mapfil : text); forward;

  [global] procedure prtdamat(n : integer; var map : DAmap); forward;

  { Implementation }

  procedure DAeps{eps : double};

  begin
    DAeps1 := eps;
  end;


  procedure DAini{o, n, u : integer};

  begin
    DAeps1 := 1d-14;
    if (0 < n) and (n <= nvmax) then
      DAdim:=n
    else
      writeln('** DAdimension overflow **');
  end;


  procedure lieini{no, nv, nd2i : integer};

  begin
    if (0 < nv) and (nv <= nvmax) then
      ndim2 := nd2i
    else
      writeln('** DAdimension overflow **');
  end;


  procedure DAall{var x : DAvect; nd2 : integer};

  var	j	: integer;

  begin
    for j:=0 to DAdim do
      x[j] := 0.0;
  end;


  procedure DAdal{var x : Davect; DAdim : integer};

  begin
  end;


  procedure DAmapall{var x : DAmap; nd2 : integer;
                     DAname : DAnambuf; no, nv : integer};

  var	j	: integer;

  begin
    for j:=1 to nd2 do
      DAall(x[j], 1, DAname, no, nv);
  end;


  procedure DAmapdal{var x : DAmap; nd2 : integer};

  var	j	: integer;

  begin
    for j:=1 to nd2 do
      DAdal(x[j], 1);
  end;


  procedure DAvar{var x : DAvect; r : double; i : integer};

  var	j	: integer;

  begin
    x[0] := r;
    for j:=1 to DAdim do
      x[j] := 0.0;
    x[i] := 1.0;
  end;


  procedure DAcon{var x : DAvect; r : double};

  var	j	: integer;

  begin
    x[0] := r;
    for j:=1 to DAdim do
      x[j] := 0.0;
  end;


  procedure ETini{var map : DAmap};

  var	i	: integer;

  begin
    for i:=1 to ndim2 do
      DAvar(map[i], 0.0, i);
  end;


  procedure DApek{var x : DAvect; var jj : ivector; var r : double};

  var	i, nzero	: integer;

  begin
    nzero := 0;
    for i:=1 to DAdim do
    begin
      if jj[i] = 0 then
        nzero := succ(nzero)
      else if jj[i] = 1 then
        r := x[i]
      else
        writeln('Invalid jj in DApek');
    end;
    if nzero = DAdim then
      r := x[0]
    else if nzero < DAdim-1 then
      writeln('Invalid jj in DApek');
  end;


  procedure DApok{var x : DAvect; var jj : ivector; r : double};

  var	i, nzero	: integer;

  begin
    nzero := 0;
    for i:=1 to DAdim do
    begin
      if jj[i] = 0 then
        nzero := succ(nzero)
      else if jj[i] = 1 then
        x[i] := r
      else
        writeln('Invalid jj in DApek');
    end;
    if nzero = DAdim then
      x[0] := r
    else if nzero < DAdim-1 then
      writeln('Invalid jj in DApek');
  end;


  function getmat{var map : DAmap; i, j : integer}{ : double};

  var	k	: integer;
	r	: double;
	jj	: ivector;

  begin
    for k:=1 to DAdim do
      jj[k] := 0;
    if j <> 0 then jj[j] := 1;
    DApek(map[i], jj, r); getmat := r;
  end;


  procedure putmat{var map : DAmap; i, j : integer; r : double};

  var	k	: integer;
	jj	: ivector;

  begin
    for k:=1 to DAdim do
      jj[k] := 0;
    if j <> 0 then jj[j] := 1;
    DApok(map[i], jj, r);
  end;


  procedure getlinmat{nv : integer; var map : DAmap; var mat : matrix};

  var	j, k	: integer;

  begin
    for j:=1 to nv do
      for k:=1 to nv do
	mat[j, k] := getmat(map, j, k);
  end;


  procedure putlinmat{nv : integer; var mat : matrix; var map : DAmap};

  { Puts zeroes in constant part of DA map }

  var	j, k	: integer;

  begin
    for j:=1 to nv do
      for k:=0 to nv do
        if k = 0 then
          putmat(map, j, k, 0.0)
	else
          putmat(map, j, k, mat[j, k]);
  end;


  procedure DAcop{var x, z : DAvect};

  var	i	: integer;

  begin
    for i := 0 to DAdim do
      z[i] := x[i];
  end;


  procedure CopyMap{var map1, map2 : DAmap};

  var	i	: integer;

  begin
    for i := 1 to DAdim do
      DAcop(map1[i], map2[i]);
  end;


  procedure DAadd{var x, y, z : DAvect};

  var	i	: integer;

  begin
    for i := 0 to DAdim do
      z[i] := x[i] + y[i];
  end;


  procedure DAsub{var x, y, z : DAvect};

  var	i	: integer;

  begin
    for i := 0 to DAdim do
      z[i] := x[i] - y[i];
  end;


  procedure DAmul{var x, y, z : DAvect};

  var	i	: integer;

  begin
    z[0] := x[0] * y[0];
    for i := 1 to DAdim do
      z[i] := x[0] * y[i] + x[i] * y[0];
  end;


  procedure DAdiv{var x, y, z : DAvect};

  var	yinv	: DAvect;

  begin
    DAfun('INV ', y, yinv); DAmul(x, yinv, z);
  end;


  procedure DAcad{var x : DAvect; y : double; var z : DAvect};

  var	i	: integer;

  begin
    z[0] := x[0] + y;
    for i := 1 to DAdim do
      z[i] := x[i];
  end;


  procedure DAcsu{var x : DAvect; y : double; var z : DAvect};

  begin
    DAcad(x, -y, z)
  end;


  procedure DAcmu{var x : DAvect; y : double; var z : DAvect};

  var	i	: integer;

  begin
    for i := 0 to DAdim do
      z[i] := x[i] * y;
  end;


  procedure DAsuc{var x : DAvect; y : double; var z : DAvect};

  var	negx	: DAvect;

  begin
    DAcmu(x, -1.0, negx); DAcad(negx, y, z)
  end;


  procedure DAcdi{var x : DAvect; y : double; var z : DAvect};

  begin
    DAcmu(x, 1/y, z);
  end;


  procedure DAdic{var x : DAvect; y : double; var z : DAvect};

  var	xinv	: DAvect;

  begin
    DAfun('INV ', x, xinv); DAcmu(xinv, y, z);
  end;


  procedure DApos{var x, z : DAvect};

  begin
    if x[0] < 0 then DAcmu(x, -1.0, z);
  end;


  procedure DAsqr{var x, z : DAvect};

  begin
    DAmul(x, x, z);
  end;


  procedure DAcma{var x, y : DAvect; rb : double; var z : DAvect};

  var	x1	: DAvect;

  begin
    DAcmu(y, rb, x1); DAadd(x, x1, z);
  end;


  procedure DAlin{var x : DAvect; ra : double; var y : DAvect; rb : double;
		  var z : DAvect};

  var	x1, x2	: DAvect;

  begin
    DAcmu(x, ra, x1); DAcmu(y, rb, x2); DAadd(x1, x2, z);
  end;


  procedure DAfun{fun : funnambuf; var x, z : DAvect};

  procedure DAinv(var x, z : DAvect);

  var	i	: integer;
	a	: double;

  begin
    z[0] := 1/x[0]; a := -1/sqr(x[0]);
    for i := 1 to DAdim do
      z[i] := a * x[i];
  end;

  procedure DAsqrt(var x, z : DAvect);

  var	i	: integer;
	a	: double;

  begin
    a := sqrt(x[0]); z[0] := a; a := 0.5/a;
    for i := 1 to DAdim do
      z[i] := a * x[i];
  end;

  procedure DAexp(var x, z : DAvect);

  var	i	: integer;
	a	: double;

  begin
    a:= exp(x[0]); z[0] := a;
    for i := 1 to DAdim do
      z[i] := a*x[i];
  end;

  procedure DAlog(var x, z : DAvect);

  var	i	: integer;

  begin
    z[0] := ln(x[0]);
    for i := 1 to DAdim do
      z[i] := x[i] / x[0];
  end;

  procedure DAsin(var x, z : DAvect);

  var	i	: integer;
	a	: double;

  begin
    z[0] := sin(x[0]); a := cos(x[0]);
    for i := 1 to DAdim do
      z[i] := a * x[i];
  end;

  procedure DAcos(var x, z : DAvect);

  var	i	: integer;
	a	: double;

  begin
    z[0] := cos(x[0]); a := -sin(x[0]);
    for i := 1 to DAdim do
      z[i] := a * x[i];
  end;

  procedure DAtan(var x, z : DAvect);

  var	c, s	: DAvect;

  begin
    DAcos(x, c); DAsin(x, s); DAdiv(s, c, z);
  end;

  procedure DAarctan(var x, z : DAvect);

  var	i	: integer;
	a	: double;

  begin
    a := x[0]; z[0] := ArcTan(a); a := 1/(1+sqr(a));
    for i := 1 to DAdim do
      z[i] := a * x[i];
  end;

  begin
    if fun = 'INV ' then DAinv(x, z)
    else if fun = 'SQRT' then DAsqrt(x, z)
    else if fun = 'EXP ' then DAexp(x, z)
    else if fun = 'LOG ' then DAlog(x, z)
    else if fun = 'SIN ' then DAsin(x, z)
    else if fun = 'COS ' then DAcos(x, z)
    else if fun = 'TAN ' then DAtan(x, z)
    else if fun = 'ATAN' then DAarctan(x, z)
    else
      writeln('illegal function');
  end;


  procedure dacct{var x : DAmap; i : integer; var y : DAmap; j : integer;
		  var z : DAmap; k : integer};

  var	l, m, n	: integer;

  begin
    for l:=1 to k do
    begin
      z[l, 0] := x[l, 0];
      for m:=1 to j do
        z[l, m] := 0.0;
      for m:=0 to j do
        for n:=1 to i do
          z[l, m] := z[l, m] + x[l, n]*y[n, m];
    end;
  end;


  procedure dainv{var x : DAmap; i : integer; var z : DAmap;
			   k : integer};

  var	mat	: matrix;

  begin
    getlinmat(i, x, mat);
    if invmat(i, mat) then
      putlinmat(k, mat, z)
    else
      writeln('map is singular');
  end;


  procedure Rotmap{n : integer; var map : DAmap; var R : matrix};

  var	mat	: matrix;

  begin
    getlinmat(n, map, mat); MulRmat(n, mat, R); putlinmat(n, mat, map);
  end;


  procedure DAwrite{var x : DAvect; var mapfil : text};

  var	i, j	: integer;
	jj	: ivector;
	zero	: boolean;

  begin
    zero := true;
    writeln(mapfil);
    writeln(mapfil, ' map  ', 1:5, ', NO =', DAord:5, ', NV =', DAdim:5,
	    ', INA =', 1:5);
    writeln(mapfil, ' *********************************************');
    writeln(mapfil);
    writeln(mapfil, '    I  COEFFICIENT          ORDER   EXPONENTS');
    for i:=1 to DAdim do
      jj[i] := 0;
    for i:=0 to DAdim do
    begin
      if abs(x[i]) > DAeps1 then
      begin
        zero := false;
        write(mapfil, i:6, ' ', x[i]:21, DAord:5, '   ');
	if i <> 0 then jj[i] := 1;
        for j:=1 to DAdim div 2 do
        begin
          write(mapfil, jj[2*j-1]:3, jj[2*j]:2);
        end;
	if i <> 0 then jj[i] := 0;
        writeln(mapfil);
        writeln(mapfil, ' ', x[i]);
      end;
    end;
    if zero then writeln(mapfil, '   ALL COMPONENTS ZERO');
    writeln(mapfil);
  end;


  procedure prtmap{var map : DAmap; var mapfil : text};

  var	i	: integer;

  begin
    for i:=1 to ndim2 do
      DAwrite(map[i], mapfil);
  end;


  procedure prtdamat{n : integer; var map : DAmap};

  var	i, j	: integer;

  begin
    writeln('DA map:');
    for i:=1 to n do
    begin
      for j:=1 to n do
        write(map[i, j]);
      writeln;
    end;
  end;

  end.
