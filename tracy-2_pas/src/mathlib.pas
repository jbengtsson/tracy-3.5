  module mathlib(input, output);

  %include 'mathlibo.def'

  { Interface }

  var	pi		: [global] double;
	rseed0, rseed	: [global] integer;
	normcut_	: [global] double;

  { Note that pi has to be initialized }

  { Extensions to the standard functions }

{  [global] FUNCTION dble(x : real) : double; forward;}

{  [global] FUNCTION sngl(x : double) : real; forward;}

  [global] FUNCTION min_(x1, x2 : double) : double; forward;

  [global] FUNCTION max_(x1, x2 : double) : double; forward;

  [global] FUNCTION pwr(x, y : double) : double; forward;
  { power <- x^y }

  [global] FUNCTION tan_(x : double) : double; forward;

  [global] FUNCTION cosh_(x : double) : double; forward;

  [global] FUNCTION sinh_(x : double) : double; forward;

  [global] FUNCTION tanh_(x : double) : double; forward;

  [global] PROCEDURE iniranf(i : integer); forward;

  [global] PROCEDURE newseed; forward;

  [global] FUNCTION ranf : double; forward;

  [global] PROCEDURE setrancut(cut : double); forward;

  [global] FUNCTION normranf : double; forward;

  { Conversion routines }

  [global] FUNCTION dtor (d : double) : double; forward;

  [global] FUNCTION sign(x : double) : integer; forward;

  [global] Function GetAngle(x, y : double) : double; forward;
  { Get the angle phi from x=cos(phi), y=sin(phi) where -pi <= phi <= pi }

  { Matrix routines }

  [global] PROCEDURE UnitMat(n : integer; VAR a : matrix); forward; 

  [global] PROCEDURE CopyVec(n : integer; VAR a, b : vector); forward;

  [global] PROCEDURE CopyMat(n : integer; VAR a, b : matrix); forward;

  [global] PROCEDURE AddVec(n : integer; VAR a, b : vector); forward;

  [global] PROCEDURE SubVec(n : integer; VAR a, b : vector); forward;

  [global] PROCEDURE AddMat(n : integer; VAR a, b : matrix); forward;

  [global] PROCEDURE SubMat(n : integer; VAR a, b : matrix); forward;

  [global] PROCEDURE LinTrans(n : integer; VAR a : matrix; VAR x : vector); forward;

  [global] PROCEDURE MulLMat(n : integer; VAR a, b : matrix); forward;

  [global] PROCEDURE MulRMat(n : integer; VAR a, b : matrix); forward;

  [global] FUNCTION TrMat(n : integer; VAR a : matrix) : double; forward;

  [global] PROCEDURE TpMat(n : integer; VAR a : matrix); forward;

  [global] FUNCTION DetMat(n : integer; VAR a : matrix) : double; forward;

  [global] function InvMat(n : integer; VAR a : matrix) : boolean; forward;

  [global] procedure prtmat(n : integer; var a : matrix); forward;

  { Implementation }

{  FUNCTION sngl}{x : double}{ : real}{;

  begin
    sngl := x;
  end;}


{  FUNCTION dble}{x : real}{ : double}{;

  begin
    dble := x;
  end;}


  FUNCTION min_{x1, x2 : double}{ : double};

  begin
    if x1 < x2 then
      min_ := x1
    else
      min_ := x2;
  end;


  FUNCTION max_{x1, x2 : double}{ : double};

  begin
    if x1 > x2 then
      max_ := x1
    else
      max_ := x2;
  end;


  function pwr{x, y : double) : double};

  { pwr <- x^y  }

  begin
    if x = 0d0 then
    begin
      if y = 0d0 then
      begin
	pwr := 1d0;
	writeln('  *** 0^0 = 1');
      end
      else
	pwr := 0d0;
    end
    else if x > 0d0 then
      pwr := exp(y*ln(x))
    else
    begin
      if y = round(y) then
      begin
	if not odd(round(y)) then
          pwr := exp(y*ln(-x))
	else
          pwr := -exp(y*ln(-x))
      end
      else
	writeln('  *** undefined exponentiation');
    end;
  end;


  FUNCTION tan_{x : double}{ : double};

  BEGIN
    tan_:=sin(x)/cos(x)
  END;


  FUNCTION cosh_{x : double}{ : double};

  VAR  y : double  ;

  BEGIN
    y := exp(x); cosh_ := (y + 1/y)/2;
  END;


  FUNCTION sinh_{x : double}{ : double};

  VAR  y : double  ;

  BEGIN
    y := exp(x); sinh_ := (y - 1/y)/2
  END;


  FUNCTION tanh_{x : double}{ : double};

  VAR  y, z : double  ;

  BEGIN
    y := exp(x); z := 1/y; tanh_ := (y - z)/(y + z)
  END;


  PROCEDURE iniranf{i : integer};

  { Initialize random number generator with seed i }

  BEGIN 
     rseed0 := i; rseed := i;
  END;


  PROCEDURE newseed;

  { Get new seed }

  const k=19; c=656329; m=100000001;

  BEGIN 
     rseed0 := ( k*rseed0 + c ) mod m;
     rseed := ( rseed0 + 54321 ) mod m;
  END;


  FUNCTION ranf {: double};

  { Generate random number with rectangular distribution }

  const k=19; c=656329; m=100000001;

  BEGIN
    rseed := ( k*rseed + c ) mod m;
    ranf := rseed/1e8
  END;


  PROCEDURE setrancut{cut : double};

  { Set cut for normal distribution }

  begin 
    normcut_ := cut;
  end;


  FUNCTION normranf {: double};

  { Generate random numbers with normal distribution (m=0, sigma=1)
    and cut normcut_ }

  const	maxiter = 100;

  var	i, j	: Integer;
	f, w	: double;

  begin
    j := 0;
    repeat
      j := succ(j);
      w := 0d0;
      for i:=1 to 12 do
	w := w + ranf;
      f := w - 6.0;
    until (abs(f) <= abs(normcut_)) or (j > maxiter) ;
    if j > maxiter then writeln('  *** fatal error in normranf');
    normranf := f;
  end;


  FUNCTION dtor {d : double}{ : double};
  
  BEGIN
    if pi = 0d0 then writeln('** pi not initialized');
    dtor := d * pi/180.0;
  END;


  FUNCTION sign{x : double}{ : integer};

  BEGIN
    if x >= 0d0 then
      sign := 1
    else
      sign := -1;
  END;


  Function GetAngle{x, y : double}{ : double};

  { Get the angle phi from x=cos(phi), y=sin(phi) where -pi <= phi <= pi }

  var	z	: double;

  begin
    if pi = 0d0 then writeln('** pi not initialized');
    if x <> 0 then
      z := ArcTan(y/x)
    else
      z := sign(y)*pi/2;
    if x < 0 then
    begin
      if y >= 0 then
	z := z + pi
      else
	z := z - pi;
    end;
    GetAngle := z;
  end;


  PROCEDURE UnitMat{n : integer; VAR a : matrix}; 

  var	i,j	:integer;

  begin
    for i:=1 to n  do
      for j:=1 to n do
	if i = j then
          a[i,j] := 1
	else
	  a[i,j] := 0d0;
  end;


  PROCEDURE CopyVec{n : integer; VAR a, b : vector};

  VAR	i	: integer; 

  BEGIN
    FOR i:=1 TO n DO 
      b[i] := a[i];
  END;


  PROCEDURE CopyMat{n : integer; VAR a, b : matrix};

  VAR	i	: integer;

  BEGIN
    FOR i:=1 TO n DO 
      CopyVec(n, a[i], b[i]);
  END;


  PROCEDURE AddVec{n : integer; VAR a, b : vector};

  VAR	i	: integer;

  BEGIN
    FOR i:=1 TO n DO 
      b[i] := a[i] + b[i];
  END;


  PROCEDURE SubVec{n : integer; VAR a, b : vector};

  VAR	i	: integer;

  BEGIN
    FOR i:=1 TO n DO 
      b[i] := b[i] - a[i];
  END;


  PROCEDURE AddMat{n : integer; VAR a, b : matrix};

  VAR	i	: integer;

  BEGIN
    FOR i:=1 TO n DO
      AddVec(n, a[i], b[i]);
  END;


  PROCEDURE SubMat{n : integer; VAR a, b : matrix};

  VAR	i	: integer;

  BEGIN
    FOR i:=1 TO n DO
      SubVec(n, a[i], b[i]);
  END;


  PROCEDURE LinTrans{n : integer; VAR a : matrix; VAR x : vector};

    VAR	i, j	: integer;
	y	: vector;

    BEGIN
      FOR i:=1 TO n DO
      BEGIN
        y[i] := 0d0;
        FOR j:=1 TO n DO
          y[i] := y[i] + a[i, j] * x[j];
      END;
      CopyVec(n, y, x);
    END;


  PROCEDURE MulLMat{n : integer; VAR a, b : matrix};

  VAR	i, j, k	: integer;
	x 	: double;
	c	: matrix;

  BEGIN
    FOR i:=1 TO n DO
      FOR j:=1 TO n DO
      BEGIN
        x := 0d0;
        FOR k := 1 TO n DO 
          x := x + a[i, k] * b[k, j];
        c[i, j] := x;
      END;
    CopyMat(n, c, b);
  END;


  PROCEDURE MulRMat{n : integer; VAR a, b : matrix};

  VAR	i, j, k	: integer;
	x 	: double;
	c	: matrix;

  BEGIN
    FOR i:=1 TO n DO
      FOR j:=1 TO n DO
      BEGIN
        x := 0d0;
        FOR k := 1 TO n DO 
          x := x + a[i, k] * b[k, j];
        c[i, j] := x;
      END;
    CopyMat(n, c, a);
  END;


  FUNCTION TrMat{n : integer; VAR a : matrix}{ : double};

    VAR	i	: integer;
	x	: double;

    BEGIN
      x := 0d0;
      FOR i:=1 TO n DO
        x := x + a[i, i];
      TrMat := x;
    END;


  PROCEDURE TpMat{n : integer; VAR a : matrix};

    VAR	i, j	: integer;
	x	: double;

    BEGIN
      FOR i:=1 TO n DO
        FOR j:=1 TO i-1 DO
	begin
          x := a[i, j]; a[i, j] := a[j, i]; a[j, i] := x;
	end;
    END;


  FUNCTION DetMat{n : integer; VAR a : matrix}{ : double};

    var		j	: integer;
		d	: double;
		cross	: array [1..matdim] of boolean;


    FUNCTION GdetMat(n : integer) : double;

      VAR	k, sign	: integer;
		det	: double;

      BEGIN
        if n > 1 then
        begin
          det := 0d0;
	  if n mod 2 = 1 then
	    sign := 1
	  else
	    sign := -1;
          for k:=1 to matdim do
	  begin
 	    if not cross[k] then
	    begin
	      cross[k] := true;
	      det := det + sign*a[n, k]*GdetMat(n-1);
	      cross[k] := false;
	      sign := -sign;
	    end;
	  end;
	  GdetMat := det;
        end
        else
          for k:=1 to matdim do
            if not cross[k] then GdetMat := a[n, k];
      END;

    begin
      if n = 2 then
        DetMat := a[1, 1]*a[2, 2] - a[1, 2]*a[2, 1]
      else if n = 3 then
      begin
        d :=     a[1, 1]*(a[2, 2]*a[3, 3] - a[2, 3]*a[3, 2]);
        d := d + a[1, 2]*(a[2, 3]*a[3, 1] - a[2, 1]*a[3, 3]);
        d := d + a[1, 3]*(a[2, 1]*a[3, 2] - a[2, 2]*a[3, 1]);
        DetMat := d;
      end
      else
      begin
        for j:=1 to matdim do
	  if j <= n then
	    cross[j] := false
	  else
	    cross[j] := true;
        DetMat := GdetMat(n);
      end;
    end;


  function InvMat{n : integer; VAR a : matrix}{ : boolean};

  { Gauss-Jordan elimination with full pivoting,
    gaussj from Numerical Recipes }

    LABEL      1;

    TYPE	iv1	= ARRAY[1..matdim] OF integer;
		iv2	= ARRAY[1..matdim, 1..2] OF integer;
		v1	= ARRAY[1..matdim] OF double;

    VAR	i, j, k, l, l1, row, column	: integer;
	amax, t, determ, d		: double;
	b				: matrix;
	ipivot				: iv1;
	index				: iv2;
	pivot				: v1;


    PROCEDURE swap (VAR x, y : double);

      VAR      d : double;

      BEGIN
	d := x; x := y; y := d
      END;


    PROCEDURE interchange;

      VAR        l : integer;

    BEGIN
      IF row <> column THEN
      BEGIN
        determ := -determ;
        FOR l := 1 TO n DO 
          swap(a[row, l], a[column, l]);
      END
    END;  {interchange}

  BEGIN {matinv}
    if n = 2 then
    begin
      d := detMat(2, a);
      IF d <> 0d0 THEN
      BEGIN
        InvMat := true;
        b[1, 1] :=  a[2, 2]/d;
        b[1, 2] := -a[1, 2]/d;
        b[2, 1] := -a[2, 1]/d;
        b[2, 2] :=  a[1, 1]/d;
	CopyMat(n, b, a);
      END
      ELSE
        InvMat := false;
    end
    else
    begin
      determ := 1.0;
      FOR j := 1 TO n DO ipivot[j] := 0;
      i := 1;
      WHILE (i <= n) AND (determ <> 0d0) DO
      BEGIN
        amax := 0d0;
        FOR j := 1 TO n DO
          IF ipivot[j] <> 1 THEN
            FOR k := 1 TO n DO
            BEGIN
              IF ipivot[k] > 1 THEN GOTO 1;
              IF ipivot[k] < 1 THEN
              BEGIN
                IF abs(amax) < abs(a[j, k]) THEN
                BEGIN
		  row := j; column := k; amax := a[j, k]
	        END
              END
	    END;
            IF amax = 0d0 THEN
            begin
	      InvMat := false;
	      determ := 0d0;
	    end
            ELSE
            BEGIN
	      InvMat := true;
              ipivot[column] := ipivot[column] + 1;
              interchange;
              index[i, 1] := row; index[i, 2] := column;
              pivot[i] := a[column, column];
              determ := determ * pivot[i];
              a[column, column] := 1.0;
              FOR l := 1 TO n DO 
                a[column, l] := a[column, l] / pivot[i];
              FOR l1 := 1 TO n DO
              BEGIN
                IF l1 <> column THEN
                BEGIN
                  t := a[l1, column];
                  a[l1, column] := 0d0;
                  FOR l := 1 TO n DO 
                    a[l1, l] := a[l1, l] - a[column, l] * t;
                END
              END
            END; {else }
            i := i + 1;
          END; {while}
          IF determ <> 0d0 THEN
            FOR i := 1 TO n DO
            BEGIN
              l := n + 1 - i;
              IF index[l, 1] <> index[l, 2] THEN
              BEGIN
		row := index[l, 1]; column := index[l, 2];
                FOR k := 1 TO n DO
		  swap(a[k, row], a[k, column])
              END
           END;
    end;
1 :
  END;


  procedure prtmat{n : integer; var a : matrix};

  var	i, j	: integer;

  begin
    writeln('matrix:');
    for i:=1 to n do
    begin
      for j:=1 to n do
        write(a[i, j]);
      writeln;
    end;
  end;

  end.
