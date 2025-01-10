  module t2bump(input, output);

  %include 'mathlib.def'
  %include 'dab.def'
  %include 'pascommon.def'
  %include 't2common.def'
  %include 'pascalio.def'
  %include 't2elem.def'
  %include 't2cell.def'
  %include 't2bumpo.def'

  { Closed orbit correction with local bump method. Local bump implies

	theta1 = free parameter

		     ___________
	theta2 = - \/beta1/beta2 sin(phi13)/sin(phi23) theta1

		   ___________
	theta3 = \/beta1/beta3 sin(phi12)/sin(phi23) theta1

    and the new orbit is given by

		 _____________
	y(s) = \/beta1 beta(s)	theta1 sin[phi(s)-phi1], s1 <= s <= s2

		 _____________
	y(s) = \/beta1 beta(s)	theta1 sin[phi(s)-phi1]

		 _____________
	     + \/beta2 beta(s)	theta2 sin[phi(s)-phi2], s2 <= s <= s3

    The rms orbit is minimized

			====
			\		    2
			 >   (x  + a  theta)
			/      i    i
			====

    so that

				   ====
				   \
				    >   a  x
				   /     i  i
				   ====	
			theta = - ------------
				    ====
				    \      2
				     >   a
				    /     i
				    ====	

  }

  { interface }

  var	Nbumps	: [global] ivector2;
	bumps	: [global] array [horizontal..vertical] of bumparray;
	bumpf	: [global] text;
	cod	: [global] array [1..bpmmax, horizontal..vertical] of double;


  [global] procedure getbumprec(ncorr : integer; var corr : fitvect;
				plane, corr1, corr2, corr3 : integer;
				var bump : bumprec); forward;

  [global] procedure SetUpBump(ncorr : integer; var corr : fitvect;
			       dnumin : double; plane : integer); forward;

  [global] procedure setlocbump(plane : integer; maxkick : double;
			        theta : double; var bump : bumprec;
				lastpos : integer); forward;

  [global] procedure execbump(MaxKick : double; lastpos : integer); forward;

  [global] procedure InitBUMP(var ncorr : ivector2; var hcorr, vcorr : fitvect;
			      dnuhmin, dnuvmin : double); forward;

  [global] procedure getbpm(Fnum, Knum : integer; var x, y : double;
			    lastpos : integer); forward;

  { implementation }

  procedure getbpm{Fnum, Knum : integer; var x, y : double; lastpos : integer};

  { Get Knum:th bpm reading }

  begin
    with Cell[ElemFam[Fnum].Kidlist[Knum]] do
    begin
      if Elem_Getpos(Fnum, Knum) < lastpos then
      begin
	x := beampos[1]; y := beampos[3];
      end
      else
      begin
	x := 0d0; y := 0d0;
      end
    end;
  end;


  procedure setcorr(i, plane : integer; strength : double);

  { Set strength of i:th corrector }

  begin
    with Cell[i], Elem, M^ do
      if plane = horizontal then
      begin
	PBpar[dip] := -strength;
	mpole_setPB(Fnum, Knum, dip);
      end
      else if plane = vertical then
      begin
	PBpar[-dip] := strength;
	mpole_setPB(Fnum, Knum, -dip);
      end;
  end;


  function getcorr(i, plane : integer) : double;

  { Get strength of i:th corrector }

  begin
    if plane = horizontal then
      getcorr := -Elem_Getkval(Cell[i].Fnum, Cell[i].Knum, dip)
    else if plane = vertical then
      getcorr := Elem_Getkval(Cell[i].Fnum, Cell[i].Knum, -dip);
  end;


  procedure ClearCorr(ncorr : ivector2; var hcorr, vcorr : fitvect);

  { Clear Corrector }

  var	i	: integer;

  begin
    for i:=1 to ncorr[1] do
      setcorr(hcorr[i], horizontal, 0.0);
    for i:=1 to ncorr[2] do
      setcorr(vcorr[i], vertical, 0.0);
  end;


  procedure getbumprec{plane, corr1, corr2, corr3 : integer;
		       var bump : bumprec};

  var	k, kmod, n2, n3, ind	: integer;
	beta1, beta2, beta3	: double;
	nu1, nu2, nu3		: double;
	sin12, sin13, sin23	: double;
	nubpm			: double;


  procedure getBN(i : integer; var beta1, nu1 : double);

  begin
    with Cell[i] do
    begin
      beta1 := beta[plane]; nu1 := nu[plane];
    end;
  end;


  begin
    with bump do
    begin
      if trace then writeln(bumpf, '  Nbumps =', Nbumps[plane]:4);
      getBN(corr1, beta1, nu1);
      getBN(corr2, beta2, nu2);
      if corr2 < corr1 then nu2 := nu2 + globval.totaltune[plane];
      getBN(corr3, beta3, nu3);
      if corr3 < corr1 then nu3 := nu3 + globval.totaltune[plane];
      sin12 := sin( (nu2 - nu1)*2*Pi );
      sin13 := sin( (nu3 - nu1)*2*Pi );
      sin23 := sin( (nu3 - nu2)*2*Pi );
      corrcoeff[1] := 1.0;
      corrcoeff[2] := -sqrt( beta1/beta2 )*sin13/sin23;
      corrcoeff[3] := sqrt( beta1/beta3 )*sin12/sin23;

      Nbpm := 0; bpmcoeffsum := 0.0; k := 0;
      with ElemFam[globval.bpm] do
      begin
	repeat
	  k := succ(k);
        until (Kidlist[k] > corr1) or (k = nKid) ;
        if Kidlist[k] <= corr1 then
        begin
          kmod := 1; k := 1 + nKid;
	  ind := Kidlist[kmod] + globval.Cell_nLoc;
        end
	else
        begin
          kmod := k; ind := Kidlist[kmod];
        end;
	n2 := corr2;
	if corr2 < corr1 then n2 := n2 + globval.Cell_nLoc;
	n3 := corr3;
	if corr3 < corr1 then n3 := n3 + globval.Cell_nLoc;
        while ind < n3 do
        begin
	  with Cell[Kidlist[kmod]] do
	  begin
	    nbpm := succ(nbpm);
	    bpmind[nbpm] := kmod;
	    nubpm := nu[plane];
	    if k > nKid then nubpm := nubpm + globval.totaltune[plane];
	    bpmcoeff[nbpm] := sqrt( beta1*beta[plane] )*sin( (nubpm-nu1)*2*Pi );
	    if ind >= n2 then
	      bpmcoeff[nbpm] := bpmcoeff[nbpm] + sqrt( beta2*beta[plane] )
				*sin( (nubpm-nu2)*2*Pi )*corrcoeff[2];
	    if trace then
	      writeln(bumpf, '  nbpm =', nbpm:4, 
		      ', bpmind =', bpmind[nbpm]:4,
		      ', bpmcoeff =', bpmcoeff[nbpm]:10:6);
	    bpmcoeffsum := bpmcoeffsum + sqr( bpmcoeff[nbpm] );
	  end;
	  k := succ(k); kmod := (k-1) mod nKid + 1;
	  ind := Kidlist[kmod]; if k > nKid then ind := ind + globval.Cell_nLoc;
	end;
      end;
      if trace then writeln(bumpf, '  bpmcoeffsum =', bpmcoeffsum:10:6);
      corrind[1] := corr1;
      corrind[2] := corr2;
      corrind[3] := corr3;
      if trace then
      begin
	writeln(bumpf, '  corrind[1] =', corrind[1]:4,
		', corrcoeff[1] =', corrcoeff[1]:10:6);
	writeln(bumpf, '  corrind[2] =', corrind[2]:4,
		', corrcoeff[2] =', corrcoeff[2]:10:6);
	writeln(bumpf, '  corrind[3] =', corrind[3]:4,
		', corrcoeff[3] =', corrcoeff[3]:10:6);
      end;
    end;
  end;


  procedure SetUpBump{ncorr : integer; var corr : fitvect;
		      dnumin : double; plane : integer};

  { Set up coefficients for closed orbit correction with local bump method. }

  var	i, j			: integer;
	corr1, corr2, corr3	: integer;


  function getnextcorr(var k : integer) : integer;

  var	nu1, nu2	: double;

  begin
    nu1 := Cell[k].nu[plane]; nu2 := nu1;
    while nu2-nu1 < dnuMin do
    begin
      k := succ(k);
      with Cell[corr[k]] do
      begin
        if k <= ncorr then
          nu2 := nu[plane]
        else
        begin
	  k := k mod ncorr;
	  nu2 := nu[plane] + globval.totaltune[plane];
        end;
      end;
    end;
    getnextcorr := corr[k];
  end;


  begin
    if trace then rewrite_(bumpf, 'bumpfile.dat');
    Nbumps[plane] := 0;
    for i:=1 to ncorr do
    begin
      Nbumps[plane] := succ(Nbumps[plane]);
      corr1 := corr[i];
      j := i;
      corr2 := getnextcorr(j); corr3 := getnextcorr(j);
      getbumprec(ncorr, corr, plane, corr1, corr2, corr3,
		 bumps[plane, Nbumps[plane]]);
    end;
    if trace then close(bumpf);
  end;


  function getlocbump(plane : integer; var bump : bumprec) : double;

  { Get local bump }

  var	i	: integer;
	sum	: double;

  begin
    with bump do
    begin
      sum := 0.0;
      for i:=1 to Nbpm do
        sum := sum + bpmcoeff[i]*cod[bpmind[i], plane];
      if bpmcoeffsum <> 0.0 then
        getlocbump := -sum/bpmcoeffsum
      else
      begin
	getlocbump := 0.0;
	writeln('  ** bpmcoeffsum = 0.0');
      end;
    end;
  end;


  procedure setlocbump{plane : integer; maxkick : double;
		       theta : double; var bump : bumprec; lastpos};

  { Set local bump }

  const	bumpred=0.9;

  var	i	: integer;
	thetal	: double;
	corrstr	: array [1..3] of double;

  begin
    thetal := theta;
    with bump do
    begin
      { Get old corrector strengths }
      for i:=1 to 3 do
        corrstr[i] := getcorr(corrind[i], plane);
      { Reduce strength until not bigger than MaxKick }
      while ( abs(corrstr[1]+thetal*corrcoeff[1]) > MaxKick ) or
	    ( abs(corrstr[2]+thetal*corrcoeff[2]) > MaxKick ) or
	    ( abs(corrstr[3]+thetal*corrcoeff[3]) > MaxKick ) do
        thetal := thetal*bumpred;

      { Set new corrector strengths }
      for i:=1 to 3 do
        setcorr(corrind[i], plane, corrstr[i]+thetal*corrcoeff[i]);

      if corrind[3] > lastpos then
      begin
{	writeln('  plane:  ', plane:3);
        writeln('  corr:   ', corrstr[1]+thetal*corrcoeff[1], ' ',
	        corrstr[2]+thetal*corrcoeff[2], ' ',
	        corrstr[3]+thetal*corrcoeff[3]);
        writeln('  lastpos:', lastpos:4);}
      end;
      { Calculate new cod }
      for i:=1 to Nbpm do
      begin
        if corrind[3] > lastpos then
        begin
{          writeln('  bpm:    ', elem_getpos(globval.bpm, bpmind[i]):4,
		  bpmind[i]);
          writeln('  cod0:   ', cod[bpmind[i], plane]);}
	end;
	if Elem_Getpos(globval.bpm, bpmind[i]) < lastpos then
	  cod[bpmind[i], plane] := cod[bpmind[i], plane] + bpmcoeff[i]*thetal
	else
	  cod[bpmind[i], plane] := 0d0;
        if corrind[3] > lastpos then
{          writeln('  cod1:   ', cod[bpmind[i], plane]);}
      end;
    end;
  end;


  procedure execbump{MaxKick : double; lastpos : integer};

  { COD correction with the local bump method, MaxKick in rad }

  var	i	: integer;
	theta	: double;

  begin
    for i:=1 to ElemFam[globval.bpm].nKid do
      getbpm(globval.bpm, i, cod[i, horizontal], cod[i, vertical], lastpos);

    for i:=1 to Nbumps[horizontal] do
    begin
      theta := getlocbump(horizontal, bumps[horizontal, i]);
      with bumps[horizontal, i] do
        if elem_getpos(globval.bpm, bpmind[1]) < lastpos then
	  setlocbump(horizontal, maxkick, theta, bumps[horizontal, i], lastpos);
    end;
    for i:=1 to Nbumps[vertical] do
    begin
      theta := getlocbump(vertical, bumps[vertical, i]);
      with bumps[vertical, i] do
        if elem_getpos(globval.bpm, bpmind[1]) < lastpos then
          setlocbump(vertical, maxkick, theta, bumps[vertical, i], lastpos);
    end;
  end;


  procedure InitBUMP{var ncorr : ivector2; var hcorr, vcorr : fitvect;
		     dnuhmin, dnuvmin : double};

  { Initialize BUMP }

  begin
    ClearCorr(ncorr, hcorr, vcorr);
    SetUpBump(ncorr[1], hcorr, dnuhmin, horizontal);
    SetUpBump(ncorr[2], vcorr, dnuvmin, vertical);
  end;

  end.
