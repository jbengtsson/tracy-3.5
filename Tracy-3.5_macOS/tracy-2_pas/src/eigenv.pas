  module eigenv(input, output);

  %include 'mathlib.def'
  %include 'ety.def'
  %include 't2common.def'
  %include 'eigenv.def'

  { Interface }

  [global] procedure gdiag(n : integer; C : double; var A, Ainv, R, M : matrix;
		           var Omega, alphac : double); forward;

  { Implementation }

  [hidden] procedure geigen(n : integer; var fm : matrix; var Vre, Vim : matrix;
		   	    var wr, wi : vector);

  { This routine finds the eigenvalues and eigenvectors of the full matrix fm }

  label 999;

  var	info, i, j, k, c		: integer;
	ort				: vector;
	cosfm, sinfm, nufm1, nufm2	: vector;
	nu1, nu2			: vector;
	aa, vv				: matrix;

  function closest(x, x1, x2, x3 : double) : integer;
  var	dx1, dx2, dx3	: double;
  begin
    dx1 := abs(x-x1); dx2 := abs(x-x2); dx3 := abs(x-x3);
    if (dx1 < dx2) and (dx1 < dx3) then
      closest := 1
    else if (dx2 < dx1) and (dx2 < dx3) then
      closest := 2
    else
      closest := 3
  end;

  procedure SwapMat(n : integer; var A : matrix; i, j : integer);
  {   A[*, i] <==> A[*, j]
      A[i, *] <==> A[j, *]  }
  var	x	: double;
	k	: integer;
  begin
    for k:=1 to n do
    begin
      x := a[k, i]; a[k, i] := a[k, j]; a[k, j] := x;
    end;
  end;

  procedure Swap(var x, y : double);
  var	z	: double;
  begin
    z := x; x := y; y := z;
  end;

  begin
    { copy matrix to temporary storage [the matrix aa is destroyed]}
    for i:=1 to n do
      for j:=1 to n do
        aa[i, j]:=fm[i, j];
    { Compute eigenvalues and eigenvectors using double
    precision Eispack routines: }
    ety(n, 1, n, aa, ort); etyt(n, 1, n, aa, ort, vv);
    ety2(n, 1, n, aa, wr, wi, vv, info);
    if info <> 0 then
    begin
      writeln( '  Error in eig'); goto 999;
    end;
    for i:=1 to n div 2 do
      for j:=1 to n do
      begin
        Vre[j, 2*i-1]:= vv[j, 2*i-1]; Vim[j, 2*i-1] := vv[j, 2*i];
        Vre[j, 2*i] := vv[j, 2*i-1]; Vim[j, 2*i] := -vv[j, 2*i];
      end;
      for i:=1 to n div 2 do
      begin
        j:=2*i-1;
        cosfm[i]:=(fm[j, j] + fm[j+1, j+1])/2;
        if abs(cosfm[i]) > dble(1.0) then
	begin
          writeln('cosfm[', i:1, ']=', cosfm[i]:20, ' > 1.0!');
	  goto 999;
	end;
        sinfm[i] := sign(fm[j, j+1])*SQRT(1-sqr(cosfm[i]));
        nufm1[i] := GetAngle(cosfm[i], sinfm[i])/(2*pi);
	if nufm1[i] < 0 then nufm1[i] := nufm1[i] + 1;
        if nufm1[i] <= 0.5 then
	  nufm2[i] := nufm1[i]
        else
	  nufm2[i] := 1 - nufm1[i];
        nu1[i] := GetAngle(wr[j], wi[j])/(2*pi);
	if nu1[i] < 0 then nu1[i] := nu1[i] + 1;
        if nu1[i] <= 0.5 then
	  nu2[i] := nu1[i]
        else
	  nu2[i] := 1 - nu1[i];
{        writeln(' est. nu   :', nufm2[i], ', eigenvalue:', nu2[i]);}
      end;
      for i:=1 to n div 2 do
      begin
	c := closest(nufm2[i], nu2[1], nu2[2], nu2[3]);
        if c <> i then
        begin
          j := 2*c - 1; k := 2*i - 1;
{          writeln(' swapping ', j:0, ' with ', k:0);}
          SwapMat(n, Vre, j, k); SwapMat(n, Vim, j, k);
          SwapMat(n, Vre, j+1, k+1); SwapMat(n, Vim, j+1, k+1);
          Swap(wr[j], wr[k]); Swap(wi[j], wi[k]);
          Swap(wr[j+1], wr[k+1]); Swap(wi[j+1], wi[k+1]);
          Swap(nu1[i], nu1[c]); Swap(nu2[i], nu2[c]);
        end;
      end;
      for i:=1 to n div 2 do
      if (0.5-nufm1[i])*(0.5-nu1[i]) < 0 then
      begin
        j:=2*i-1;
        SwapMat(n, Vim, j, j+1); Swap(wi[j], wi[j+1]);
      end;
999:
  end;


  procedure gdiag{n : integer; var A, Ainv, R, M : matrix; var Omega : double};
  {  Input M:     One turn transfer matrix

     Output A

                                [ c  s  0  0 ]
                 -1             [-s  c  0  0 ]
          M  -> A   M A = R =   [ 0  0  c  s ]
                                [ 0  0 -s  c ]

     The eigenvectors are normalized so that the real and
     imaginary part of vectors 1 and 3 have +1 antisymmetric
     product:

       Vre1 J aivec1 := 1 ; Vre3 J aivec3 := 1 ;

     the eigenvectors 2 and 4 have the opposite normalization. }
  var	j  				: integer;
	x1, x2				: double;
	wr, wi, eta			: vector;
	JJ, fm, Vre, Vim, B, Binv	: matrix;

  procedure InitJJ;
 
  var	i, j	: integer;

  begin
    for i:=1 to n do
      for j:=1 to n do
        JJ[i, j] := 0.0;
    for j:=1 to n do
      if odd(j) then
      begin
        JJ[j, j+1] := 1.0; JJ[j+1, j] := -1.0;
      end;
  end;

  procedure GetAinv;
  var	i, j, k, sgn	: integer;
	z		: double;
  begin
    for i:=1 to n do
      if odd(i) then
      begin
        z := 0;
        for j:=1 to n do
          for k:=1 to n do
            z := z + Vre[j, i]*JJ[j, k]*Vim[k, i];
	sgn := sign(z);
        z := sqrt(abs(1/z));
        for j:=1 to n do
	begin
          Vre[j, i] := Vre[j, i]*z;
	  Vim[j, i] := sgn*Vim[j, i]*z;
	end;
      end;
    for i:=1 to n do
    begin
      Ainv[1, i] := sign(Vre[1, 1])*Vre[i, 1];
      Ainv[2, i] := sign(Vre[1, 1])*Vim[i, 1];
      Ainv[3, i] := sign(Vre[3, 3])*Vre[i, 3];
      Ainv[4, i] := sign(Vre[3, 3])*Vim[i, 3];
      if n > 4 then 
      begin
        Ainv[5, i] := sign(Vre[5, 5])*Vre[i, 5];
        Ainv[6, i] := sign(Vre[5, 5])*Vim[i, 5];
      end;
    end;
  end;

  procedure GetEta(k : integer; var M : matrix; var Eta: vector);
  { Find dispersion and momentum compaction of the map M where
    M is a 4x4 or a 6x6 Matrix which has this form

	[ N11  N12  m1  0 ]		[ N11  N12  N13  N14  m1  0 ]
	[ N21  N22  m2  0 ]		[ N21  N22  N23  N24  m2  0 ]
	[ 0    0    1   0 ]		[ N31  N32  N33  N34  m3  0 ]
	[ n1   n2   a   1 ]		[ N41  N42  N43  N44  m4  0 ]
					[ 0    0    0    0    1   0 ]
					[ n1   n2   n3   n4   a   1 ]

	Dispersion:  		Eta   = (Inv(I - N)) x m
	Momentum Compaction:	Alpha = -(a + n x Eta)  
  }
  var	 SmallM			: vector;
	 IMinNInv, I		: matrix;
         j			: integer;
  begin
    { k = 6 in tracy }
    UnitMat(k, I); CopyMat(k, I, IMinNInv); SubMat(k, M, IMinNInv);
    If not InvMat(k-2,IminNinv) then writeln('(I-N) is singular');
    for j:=1 to k-2 do smallm[j]:=M[j,k-1];
    smallm[5]:=0.0d0;smallm[6]:=0.0d0;
    CopyVec(k,SmallM,Eta);
    LinTrans(k-2,IMinNInv,Eta);
{   Alpha := 0.0d0;
    for j:= 1 TO (k-2) DO
      Alpha := Alpha + M[6,j]*Eta[j];
    Alpha := Alpha + M[6,5]; }
  end;

  procedure GenB(k : integer; var B, BInv : matrix; var Eta : vector);
  var	j	: integer;
  begin
    UnitMat(k,B);
    for j:= 1 TO (k-2) DO
      B[j, k-1] := Eta[j];
    B[6, 1] := Eta[2]; B[6, 2] := -Eta[1];
    B[6, 3] := Eta[4]; B[6, 4] := -Eta[3];
    CopyMat(k, B, BInv);
    If not InvMat(k, BInv) then writeln('B is singular');
  end;

  begin
    InitJJ;
    CopyMat(n, M, fm); TpMat(n, fm); geigen(n, fm, Vre, Vim, wr, wi);
    for j:=1 to n div 2 do
    begin
      x1 := sqrt(sqr(wr[2*j-1])+sqr(wi[2*j-1]));
      x2 := sqrt(sqr(wr[2*j])+sqr(wi[2*j]));
      globval.rad[j] := ln(sqrt(x1*x2));
    end;
    unitmat(6, Ainv); GetAinv; CopyMat(6, Ainv, A);
    if not InvMat(6, A) then writeln('A^-1 script is singular');
    if n = 4 then
    begin
      getEta(6, M, eta); GenB(6,B, BInv, Eta);
      mullmat(6, B, A); mullmat(6, Ainv, Binv);
      copymat(6, Binv, Ainv);
    end;
    CopyMat(6, A, R); MulLMat(6, M, R); MulLMat(6, Ainv, R);
    if n = 4 then
    begin
      Omega := 0.0; alphac := R[6, 5]/C;
    end;
    if n = 6 then
    begin
      if globval.cavity_on then
      begin
	Omega := GetAngle(R[5, 5], R[5, 6])/(2*pi); alphac := 0.0;
      end
      else
      begin
	Omega := 0.0; alphac := R[6, 5]/C;
      end;
    end;
  end;

  end.
