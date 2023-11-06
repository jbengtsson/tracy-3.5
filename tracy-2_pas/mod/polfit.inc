  procedure polfit(n : integer; var x, y : graphvect; m : integer;
		   var b : vector; var sigma : double);

  { Polynomial fit of data by linear chi square }

  const	sigma0 = 1; chi2 = 4;

  var	i, j, k		: integer;
	T		: matrix;

  function gsigma : double;

  { Calculate sigma }

  var	i, j	: integer;
	sum, z	: double;

  begin
    sum := 0d0;
    for i:=1 to n do
    begin
      z := 0d0;
      for j:=1 to m do
	z := z + b[j]*pwr(x[i], j-1);
      sum := sum + sqr(y[i]-z);
    end;
    gsigma := sqrt(sum/n);
  end;

  begin
    writeln('**************** Polfit-diag *****************');
    writeln('**  n,m  ',n:10,m:10);

    for i:=1 to m do
    begin
      b[i] := 0d0;
      for j:=1 to m do
        T[i, j] := 0d0;
    end;

    for i:=1 to n do
    begin
      writeln('** i,x,y:',i:5,' ',x[i]:20,' ',y[i]:20);
      for j:=1 to m do
      begin
	b[j] := b[j] + y[i]*pwr(x[i], j-1)/sqr(sigma0);
	for k:=1 to m do
	  T[j, k] := T[j, k] + pwr(x[i], j-1+k-1)/sqr(sigma0);
      end;
    end;

    if not invmat(m, T) then
      writeln('  matrix is singular');

    writeln('** lintrans in , m=',m:5);
    lintrans(m, T, b);
    writeln('** lintrans out, vector: ');
    for i:=1 to m do writeln('** i,b   :',I:5,'  ',b[i]:20);

    sigma := gsigma;
    { b[i]:11, ' +-', sqrt(chi2*T[i, i]):11 }
  end;
