  module svdcmpmod(input, output);

  { Interface }

  %include 'svd.def'
  %include 'dsvdc.def'

  [global] PROCEDURE svdcmp(VAR A: svdmat; m, n: integer;
			    VAR w: svdarray; VAR V: svdmat);

  var	i, j, job, mn, info	: integer;
	work, e			: svdarray;
	At, Ut, Vt		: svdmat;

  begin
    mn := mnp;
    { return n left singular vectors in U and the right
      singular vectors in V }
    job := 11;
    { transpose matrix for Fortran call }
    for i:=1 to m do
      for j:=1 to n do
	At[j, i] := A[i, j];
    { dsvcmp from numerical recipies has been replaced by
      dsvdc from linpack, due to numerical problems }
    dsvdc(At, mn, m, n, w, e, Ut, mn, Vt, mn, work, job, info);
    for i:=1 to m do
      for j:=1 to n do
      begin
	A[j, i] := Ut[i, j]; V[j, i] := Vt[i, j];
      end;
  end;

  end.
