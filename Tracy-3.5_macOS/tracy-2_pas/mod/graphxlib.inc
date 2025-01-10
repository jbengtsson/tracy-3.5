  { General routines for plotting }

  const	loch=999999; floch=1e36;

  procedure getminmax(n : integer; var x : graphvect; var xmin, xmax : double);

  { Procedure to find min and max values in an array }

  var	i	: integer;

  begin
    xmin := 1d30; xmax := -1d30;
    for i:=1 to n do
    begin
      xmin := min_(x[i], xmin); xmax := max_(x[i], xmax);
    end;
  end;


  procedure plotbars(n : integer; var px, py : graphvect);

  { Plot data as bars }

  var	i		: integer;
{	px1, py1	: array [1..2] of double; for compability}
	px1, py1	: graphvect;

  begin
    for i:=1 to n do
    begin
      px1[1] := px[i]; py1[1] := 0;
      px1[2] := px[i]; py1[2] := py[i];
      graph1(px1, py1, 2, glsoli, 0, 0.0);
    end;
  end;


  procedure plotfft( n : integer; var x : graphvect);

  var	i, m	: integer;
	pos	: graphvect;
        xplt,
        yplt    : array [1..2] of double;

  begin
    m := n div 2 + 1;
    for i:=1 to m do
      pos[i] := (i-1)/n;
    graph1(pos, x, m, 0, 0, 0.0);
    for i:= 1 to m do
    begin
      xplt[1]:=pos[i];
      xplt[2]:=pos[i];
      yplt[1]:=0.0;
      yplt[2]:=x[i];
      graph1(xplt, yplt, 2,  glsoli, 0, 0.0)
    end;
  end;
