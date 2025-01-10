  module fft(input, output);

  %include 'mathlib.def'
  %include 'pascommon.def'
  %include 'fft.def'

  { Interface }

  [global] PROCEDURE FFT(n : integer; var xr, xi : bigvect); forward;

  { Implementation }

  PROCEDURE FFT{n : integer; var xr, xi : bigvect};

  { DFT using FFT algorithm }

  var   m, i, j, j1, j2, jj1, jj2, jr, jf, k, l, mh	: integer;
	ifact1, ifact2, nfact1, nfact2, nfact3		: integer;
        lcoef1, lcoef2, ncoef1, ncoef2			: integer;
	tcoef, theat, scoef, ccoef, s, c, r1, r2	: double;
	a, b, sb, cb					: double;

  begin
    m := trunc(ln(n)/ln(2));
    for l:=1 to m do
    begin
      lcoef1 := round(pwr(2, l-1)); lcoef2 := lcoef1*2;
      ncoef1 := n div lcoef2; ncoef2 := ncoef1*2;
      tcoef := (2*pi)/n; theat := lcoef1*tcoef;
      scoef := sin(theat); ccoef := -2*sqr(sin(theat*0.5));
      for k:=1 to lcoef1 do
      begin
	j1 := (k-1)*ncoef2; j2 := j1 + ncoef1;
	s := 0; c := 1;
	for j := 1 to ncoef1 do
	begin
	  jj1 := j + j1; jj2 := j + j2; r1 := xr[jj1]; r2 := xi[jj1];
	  xr[jj1] := r1 + xr[jj2]; xi[jj1] := r2 + xi[jj2];
	  a := r1 - xr[jj2]; b := r2 - xi[jj2];
	  xr[jj2] := a*c + b*s; xi[jj2] := b*c - a*s;
	  sb := ccoef*s + scoef*c + s; cb := ccoef*c - scoef*s + c;
	  s := sb; c := cb;
	end ;
      end;
    end;
    mh := m div 2;
    for l := 1 to mh do
    begin
      ifact1 := round(pwr(2, l-1)); ifact2 := ifact1*4;
      nfact1 := n div ifact2; nfact2 := nfact1*2; nfact3 := nfact2*2;
      for k:=1 to ifact1 do
      begin
        jr := (k-1)*nfact3; j1 := jr + ifact1; j2 := jr + nfact2;
        for j:=1 to nfact1 do
        begin
	  jf := (j-1) div ifact1;  jf := jf*ifact1 + j;
	  jj1 := jf + j1; jj2 := jf + j2;
	  r1 := xr[jj1]; r2 := xi[jj1];
	  xr[jj1] := xr[jj2]; xi[jj1] := xi[jj2];
	  xr[jj2] := r1; xi[jj2] := r2;
        end;
      end;
    end;
    for i:=1 to n do
    begin
      xr[i] := 2/n*xr[i]; xi[i] := 2/n*xi[i];
    end;
  end;

  end.
