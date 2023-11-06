
  { physlib.inc }

  { General routines }

  const	   nueps = 1d-6;    nudkL = 0.00002;   nuimax = 10;
	  ksieps = 1d-5;  ksidkpL = 0.01;     ksiimax = 10;
	 dispeps = 1d-10; dispdkL = 0.001;   dispimax = 10;
	npeakmax = 10;
        writetrack = true; {protocol from tracking}


  type	ipeakbuf 	= array [1..npeakmax] of integer;
	peakbuf 	= array [1..npeakmax] of double;


  { Routines for printing }

  procedure printglob;
  begin
    writeln;
    writeln('  dPcommon =     ', globval.dPcommon:7,
	    ', dPparticle =     ', globval.dPparticle:7,
	    ', Energy =     ', globval.Energy:7);
    writeln('  maxAmplx =', globval.maxampl[1]:13:5,
	    ', maxAmply   =', globval.maxampl[2]:13:5);
    write('  MatMeth  =       ', globval.MatMeth);
    write(', Cavity_On  =       ', globval.Cavity_On);
    writeln(', bpm    =          ', globval.bpm:3);
    writeln;
    write('  nux      =', globval.totaltune[1]:13:8,
	  ', nuy        =', globval.totaltune[2]:13:8);
    if globval.Cavity_on then
      writeln(', omega  =', globval.omega:13:9)
    else
    begin
      writeln;
      writeln('  ksix     =', globval.chrom[1]:13:8,
	      ', ksiy       =', globval.chrom[2]:13:8);
    end;
  end;

  procedure printmat(n : integer; var A : matrix);
  var	i, j	: integer;
  begin
    writeln;
    for i:=1 to n do
    begin
      for j:=1 to n do
	write(' ', A[i, j]:12);
      writeln;
    end;
  end;

  procedure printlatt;

  var	i, j	: integer;
	Cell	: Celltype;
	str	: tstring;
	outf	: text;

  begin
    getstr(str, 'linlat.dat');
    rewrite_(outf, str.str);
    writeln(outf, '  #', ' name', '     s  ',
	    '  alphax', '  betax', '   nux', '   etax', '   etapx',
	    '  alphay', '  betay', '   nuy', '   etay', '   etapy');
    writeln(outf, '   ', '     ', '    [m] ',
	    '        ', '   [m] ', '      ', '    [m]', '        ',
	    '        ', '   [m] ', '      ', '    [m]', '        ');
    writeln(outf);
    for i:=0 to globval.Cell_nLoc do
    begin
      getElem(i, Cell);
      write(outf, i:3,':');
      for j:=1 to 6 do
      begin
        if (i <> 0) and (Cell.Elem.Pname[j] <> ' ') then
	  write(outf, Cell.Elem.Pname[j])
	else
	  write(outf, ' ');
      end;
      writeln(outf, Cell.s:6:2, 
        Cell.Alpha[1]:8:3, Cell.Beta[1]:7:3, Cell.Nu[1]:7:3,
        Cell.Eta[1]:7:3, Cell.Etap[1]:7:3,
	Cell.Alpha[2]:8:3, Cell.Beta[2]:7:3, Cell.Nu[2]:7:3,
        Cell.Eta[2]:7:3, Cell.Etap[2]:7:3);
    end;
    close(outf);
  end;

  procedure Printcod;
  var	i, j	: integer;
        Cell	: CellType;
	str	: tstring;
	outf	: text;
  begin
    getstr(str, 'cod.dat');
    rewrite_(outf, str.str);
    for i:=1 to globval.Cell_nLoc do
    begin
      getElem(i, Cell);
{      if Cell.Fnum = globval.BPM then}
      begin
	{ COD is in local coordinates }
        write(outf, i:3,':');
        for j:=1 to 10 do
        begin
          if Cell.Elem.Pname[j] <> ' ' then
	    write(outf, Cell.Elem.Pname[j])
	else
	  write(outf, ' ');
        end;
        writeln(outf, i:5,': ', Cell.beampos[1], ', ', Cell.beampos[3],
		', ', Cell.ds[1], ', ', Cell.ds[2]);
      end;
    end;
    close(outf);
  end;

  procedure getmean(n : integer; var x : graphvect);
  var	i	: integer;
	mean	: double;
  begin
    mean := 0.0;
    for i:=1 to n do
      mean := mean + x[i];
    mean := mean/n;
    for i:=1 to n do
      x[i] := sngl(x[i] - mean);
  end;

  procedure getcod(dP : double; var lastpos : integer);
  { Get closed orbit for a particle with the mumentum deviation dP }
  begin
    Cell_getCOD(globval.CODimax, globval.CODeps, dP, lastpos);
  end;

  procedure getabn(var alpha, beta, nu : vector2);
  { Get alpha, beta nu from one turn matrix }
  var  gamma      : vector2;
  begin
    Cell_GetABGN(globval.oneturnmat, alpha, beta, gamma, nu);
  end;
 
  procedure TraceABN(i0, i1 : integer; alpha, beta, eta, etap : Vector2;
		     dP : double);
  const	n = 4;
  var	i, j, lastpos	: integer;
	sb		: double;
	Ascr		: matrix;
  begin
    UnitMat(6, Ascr);
    for i:=1 to 2 do
    begin
      sb := sqrt(beta[i]);
      j := 2*i-1;
      Ascr[j, j]   := sb;           Ascr[j, j+1]   := 0;
      Ascr[j+1, j] := -alpha[i]/sb; Ascr[j+1, j+1] := 1/sb;
    end;
    Ascr[1, 5] := eta[1]; Ascr[2, 5] := etap[1];
    Ascr[3, 5] := eta[2]; Ascr[4, 5] := etap[2];
    globval.CODvect[1] := 0.0; globval.CODvect[2] := 0.0;
    globval.CODvect[3] := 0.0; globval.CODvect[4] := 0.0;
    globval.CODvect[5] := dP; globval.CODvect[6] := 0.0;
    putglobv_(globval);
    Cell_MatTwiss(i0, i1, Ascr, false, false, dP);
    unitmat(n+1, globval.oneturnmat);
    globval.CODvect[1] := 0.0; globval.CODvect[2] := 0.0;
    globval.CODvect[3] := 0.0; globval.CODvect[4] := 0.0;
    globval.CODvect[5] := dP; globval.CODvect[6] := 0.0;
    Cell_Pass_M(0, globval.Cell_nLoc, globval.CODvect, globval.oneturnmat,
		lastpos);
  end;

  procedure ttwiss(alpha, beta, eta, etap : vector2; dP : double);
  { Calculate Twiss functions for transport line }
  begin
    traceABN(0, globval.Cell_nLoc, alpha, beta, eta, etap, dP);
  end;

  { Routines for fitting }

  PROCEDURE FitTune(qf, qd : integer; nux, nuy : double);
  { Fit tune }
  var	i		: integer;
	nq		: ivector2;
	nu		: vector2;
	qfbuf, qdbuf	: fitvect;
  begin
    nq[1] := GetnKid(qf);
    for i:=1 to nq[1] do
      qfbuf[i] := Elem_Getpos(qf, i);
    nq[2] := GetnKid(qd);
    for i:=1 to nq[2] do
      qdbuf[i] := Elem_Getpos(qd, i);
    nu[1] := nux; nu[2] := nuy;
    Ring_FitTune(nu, nueps, nq, qfbuf, qdbuf, nudkL, nuimax);
  end;  

  PROCEDURE FitChrom(sf, sd : integer; ksix, ksiy : double);
  { Fit chromaticity }
  var	i		: integer;
	ns		: ivector2;
	sfbuf, sdbuf	: fitvect;
	ksi		: vector2;
  begin
    ns[1] := GetnKid(sf);
    for i:=1 to ns[1] do
      sfbuf[i] := Elem_Getpos(sf, i);
    ns[2] := GetnKid(sd);
    for i:=1 to ns[2] do
      sdbuf[i] := Elem_Getpos(sd, i);
    ksi[1] := ksix; ksi[2] := ksiy;
{    Ring_FitChrom(ksi, ksieps, ns, sfbuf, sdbuf, 1.0, 1);}
    Ring_FitChrom(ksi, ksieps, ns, sfbuf, sdbuf, ksidkpL, ksiimax);
  end;

  PROCEDURE FitDisp(q, pos : integer; eta : double);
  { Fit dispersion }
  var	i, nq	: integer;
	qbuf	: fitvect;
  begin
    nq := GetnKid(q);
    for i:=1 to nq do
      qbuf[i] := Elem_Getpos(q, i);
    Ring_FitDisp(pos, eta, dispeps, nq, qbuf, dispdkL, dispimax);
  end;

  procedure inibump(coh, cov : integer);
  const	dnux = 0.02; dnuy = 0.01;
  var	i		: integer;
	ncorr		: ivector2;
	cohbuf, covbuf	: fitvect;
  begin
    ncorr[1] := GetnKid(coh);
    for i:=1 to ncorr[1] do
      cohbuf[i] := Elem_Getpos(coh, i);
    ncorr[2] := GetnKid(cov);
    for i:=1 to ncorr[2] do
      covbuf[i] := Elem_Getpos(cov, i);
    initbump(ncorr, cohbuf, covbuf, dnux, dnuy);
  end;

  { Routines for tracking }

  procedure getfloqs(var x : vector);
  { Get coordinates in Floquet space }
  const	n = 4;
  begin
    LinTrans(n+1, globval.AscrInv, x);
  end;

  procedure track(x, px, y, py, dp : double; nmax : integer;
                  var lastn, lastpos : integer; floq : boolean);
  { Single particle tracking around closed orbit }
  const	n = 4;
  var   i,k 		: integer;
        x1, x2, xf	: vector;
	globval		: globvalrec;
	str		: tstring;
	outf		: text;
  begin
    getcod(dp, lastpos); getglobv_(globval);
    if writetrack then begin
      getstr(str, 'DISK_APPLI:[SCRATCH.STREUN]track.dat');
      rewrite_(outf, str.str);
      writeln(outf, ' version 3.11 tracking 10/09/90 13:21:32 0');
      writeln(outf, 'Tracking with TRACY');
    end;
    if floq then
    begin
      Ring_GetTwiss(false, dp); getglobv_(globval);
      if writetrack then writeln(outf, ' (Floquet space)');
    end;
    if writetrack then 
      writeln(outf, 1:3, nmax:6, 1d0:8, 1d0:8, 0d0:8, 0d0:8,
              globval.totaltune[1], globval.totaltune[2]);
    xf[1] := x; xf[2] := px; xf[3] := y; xf[4] := py;
    if floq then getfloqs(xf);
    if writetrack then
      writeln(outf, 0:6, '  ', xf[1]:12, '  ', xf[2]:12,
    	    '  ', xf[3]:12, '  ', xf[4]:12, '  ', dp:12);
    x2[1] := x + globval.CODvect[1]; x2[2] := px + globval.CODvect[2];
    x2[3] := y + globval.CODvect[3]; x2[4] := py + globval.CODvect[4];
    x2[5] := dp; x2[6] := 0.0;
    lastn := 0; x1[5] := dp; x1[6] := 0.0;
    if globval.MatMeth then Cell_Concat(dp);
    repeat
      lastn := succ(lastn);
      for i:=1 to n do
	x1[i] := x2[i];
      if globval.MatMeth then
        Cell_fPass(x2, lastpos)
      else
        Cell_Pass(0, globval.Cell_nLoc, x2, lastpos);
      for i:=1 to n do
	xf[i] := x2[i] - globval.CODvect[i];
      if floq then getfloqs(xf);
      if writetrack then
        writeln(outf, lastn:6, '  ', xf[1]:12, '  ', xf[2]:12,
	      '  ', xf[3]:12, '  ', xf[4]:12, '  ', xf[5]:12, '  ', xf[6]:12);
    until (lastn = nmax) or (lastpos <> globval.Cell_nLoc) ;
    if globval.MatMeth then Cell_Pass(0, globval.Cell_nLoc, x1, lastpos);
    if writetrack then close(outf);
  end;

  procedure getdynap(var r : double; phi, delta, eps : double;
		     nturn : integer; floqs : boolean);
  { Binary search for dynamical aperture in Floquet space. }
  const	step = 0.1; px = 0.0; py = 0.0;
  var	dr, rmin, rmax	: double;
	lost		: boolean;

  procedure track(r : double);
  var	i, lastn, lastpos	: integer;
	x			: vector;
  begin
    x[1] := r*cos(phi); x[2] := px; x[3] := r*sin(phi); x[4] := py;
    x[5] := delta; x[6] := 0.0;
    { transform to phase space }
    if floqs then lintrans(5, globval.Ascr, x);
    for i:=1 to 4 do
      x[i] := x[i] + globval.codvect[i];
    lastn := 0;
    repeat
      lastn := succ(lastn);
      if globval.MatMeth then
        Cell_fPass(x, lastpos)
      else
        Cell_Pass(0, globval.Cell_nLoc, x, lastpos);
    until (lastn = nturn) or (lastpos <> globval.Cell_nLoc) ;
    lost := lastn <> nturn;
  end;

  begin
    if globval.MatMeth then Cell_Concat(delta);
    { find upper and lower limit }
    track(r);
    dr := r*step;
    if r <= 0.0 then dr := eps;
    if not lost then
    begin
      rmax := r;
      repeat
	dr := 2d0*dr; rmax := rmax + dr;
	track(rmax);
      until lost ;
      rmin := rmax - dr;
    end
    else
    begin
      rmin := r;
      repeat
	dr := 2d0*dr; rmin := rmin - dr;
	if rmin < 0d0 then rmin := 0d0;
	track(rmin);
      until not lost ;
      rmax := rmin + dr;
    end;
    { binary search }
    repeat
      r := rmin + (rmax-rmin)/2d0;
      track(r);
      if not lost then
	rmin := r
      else
	rmax := r;
    until rmax - rmin < eps ;
    r := rmin;
  end;

  procedure getcsAscr;
  { Get Courant-Snyder Ascr }
  var	i, j	: integer;
	phi	: double;
	R	: matrix;
  begin
    unitmat(6, R);
    for i:=1 to 2 do
    begin
      phi := -arctan(globval.Ascr[2*i-1, 2*i]/globval.Ascr[2*i-1, 2*i-1]);
      R[2*i-1, 2*i-1] := cos(phi); R[2*i, 2*i] := R[2*i-1, 2*i-1];
      R[2*i-1, 2*i] := sin(phi); R[2*i, 2*i-1] := -R[2*i-1, 2*i];
    end;
    mulrmat(6, globval.Ascr, R);
    for i:=1 to 2 do
    begin
      if globval.Ascr[2*i-1, 2*i-1] < 0.0 then
      begin
	for j:=1 to 6 do
	begin
	  globval.Ascr[j, 2*i-1] := -globval.Ascr[j, 2*i-1];
	  globval.Ascr[j, 2*i] := -globval.Ascr[j, 2*i];
	end;
      end;
    end;
    if not invmat(6, globval.Ascrinv) then
      writeln('  *** Ascr is singular');
  end;

  procedure dynap(r, delta, eps : double; npoint, nturn : integer;
		  var x, y : graphvect; floqs, print : boolean);
  const	phimax = 360.0; px = 0.0; py = 0.0;
  var	i, lastpos	: integer;
	phi		: double;
  begin
    getcod(delta, lastpos); getglobv_(globval);
    if print then
    begin
      writeln;
      if floqs then
      begin
	Ring_GetTwiss(false, delta); getglobv_(globval);
	writeln('  Dynamical Aperture (Floquet space):');
	writeln('     x^         y^');
      end
      else
      begin
	writeln('  Dynamical Aperture:');
	writeln('   x [mm]     y [mm]');
      end;
    end;
    for i:=1 to npoint do
    begin
      phi := (i-1)*dtor(phimax)/npoint + 1d-6;
      getdynap(r, phi, delta, eps, nturn, floqs);
      x[i] := sngl(r*cos(phi)); y[i] := sngl(r*sin(phi));
      if print then writeln(' ', x[i]:10, ' ', y[i]:10);
    end;
  end;

  function getaper(n : integer; var x, y : graphvect) : double;
  var	i		: integer;
	sum		: double;
  begin
    sum := 0.0;
    for i:=2 to n do
      sum := sum + x[i-1]*y[i] - x[i]*y[i-1];
    sum := sum + x[n]*y[1] - x[1]*y[n];
    sum := abs(sum/2.0);
    writeln; writeln('  Dynamical aperture [mxm]: ', sum:9);
    getaper := sum;
  end;

  procedure gettrack(var n : integer; var x, px, y, py : graphvect);
  var	i	: integer;
	str	: tstring;
	inf	: text;
  begin
    getstr(str, 'DISK_APPLI:[SCRATCH.STREUN]track.dat');
    reset_(inf, str.str);
    readln(inf); readln(inf); readln(inf);
    { Skip initial conditions }
    readln(inf);
    n := 0;
    while not eof(inf) do
    begin
      n := succ(n);
      read(inf, i); readln(inf, x[n], px[n], y[n], py[n]);
    end;
    close(inf);
  end;

  procedure getj(n : integer; var x, px, y, py : graphvect);
  { Calculates the linear invariant }
  var	i	: integer;
  begin
    for i:=1 to n do
    begin
      x[i] := (sqr(x[i]) + sqr(px[i]))/2; y[i] := (sqr(y[i]) + sqr(py[i]))/2;
    end;
  end;

  function fract(x : double) : double;
  begin
    fract := x - trunc(x);
  end;

  procedure getphi(n : integer; var x, px, y, py : graphvect);
  { Calculates the linear phase }
  var	i	: integer;

  function getarg(x, px, nu : double) : double;
  var	phi, val	: double;
  begin
    phi := getangle(x, px);
    if phi < 0.0 then phi := phi + 2*pi;
    val := phi + fract( nu )*2*pi;
    if val < 0.0 then val := val + 2*pi;
    getarg := val;
  end;

  begin
    for i:=1 to n do
    begin
      x[i] := sngl( getarg( x[i], px[i], i*globval.totaltune[1]) );
      y[i] := sngl( getarg( y[i], py[i], i*globval.totaltune[2]) );
    end;
  end;

  { Routines for Fourier analysis }

  procedure sinfft(n : integer; var xr : graphvect);
  { DFT with sine window }
  var	i	: integer;
	xi	: graphvect;
  begin
    for i:=1 to n do
    begin
      xr[i] := sngl(sin((i-1)/n*pi)*xr[i]); xi[i] := 0.0;
    end;
    fft(n, xr, xi);
    for i:=1 to n do
      xr[i] := sqrt(sqr(xr[i]) + sqr(xi[i]));
  end;

  procedure getind(n : integer; k : integer; var ind1, ind3 : integer);
  begin
    if k = 1 then
    begin
      ind1 := 2; ind3 := 2;
    end
    else if k = n div 2 + 1 then
    begin 
      ind1 := n div 2; ind3 := n div 2;
    end
    else
    begin
      ind1 := k - 1; ind3 := k + 1;
    end;
  end;

  procedure getpeak(n : integer; var x : graphvect; var k : integer);
  { Locate peak in DFT spectrum }
  var	ind1, ind2, ind3	: integer;
	peak			: double;
  begin
    peak := 0.0;
    k := 1;
    for ind2:=1 to n div 2 + 1 do
    begin
      getind(n, ind2, ind1, ind3);
      if (x[ind2] > peak) and (x[ind1] < x[ind2]) and (x[ind3] < x[ind2]) then
      begin
	peak := x[ind2]; k := ind2;
      end;
    end;
  end;

  function int2snu(n : integer; var x : graphvect; k : integer) : double;
  { Get frequency by nonlinear interpolation with two samples
    for sine window. The interpolation is:

             1              2 A(k)       1
	nu = - [ k - 1 + ------------- - - ] ,      k-1 <= N nu <= k
	     N           A(k-1) + A(k)   2
  }
  var	ind, ind1, ind3	: integer;
	ampl1, ampl2		: double;
  begin
    getind(n, k, ind1, ind3);
    if x[ind3] > x[ind1] then
    begin
      ampl1 := x[k]; ampl2 := x[ind3]; ind := k;
    end
    else
    begin
      ampl1 := x[ind1]; ampl2 := x[k];
      { Interpolate in right direction for 0 frequency }
      if k <> 1 then
	ind := ind1
      else
	ind := 0;
    end;
    { Avoid division by zero }
    if ampl1+ampl2 <> 0.0 then
      int2snu := ((ind-1)+2*ampl2/(ampl1+ampl2)-0.5)/n
    else
      int2snu := 0.0;
  end;

  function intsampl(n : integer; var x : graphvect; nu : double; 
		     k : integer) : double;
  { Get amplitude by nonlinear interpolation for sine window. The
    distribution is given by:

	          1    sin pi ( k + 1/2 )     sin pi ( k - 1/2 )
	  F(k) =  - ( -------------------- + -------------------- )
	          2      pi ( k + 1/2 )          pi ( k - 1/2 )
  }
  var	corr	: double;

  function sinomega(omega : double) : double;
  {  Function to calculate:

			sin( omega )
			------------
			   omega
  }
  begin
    if omega <> 0.0 then
      sinomega := sin(omega)/omega
    else
      sinomega := 1.0;
  end;

  begin
    corr := ( sinomega( pi*((k-1) + 0.5 - nu*n) )
	    + sinomega( pi*((k-1) - 0.5 - nu*n) ) )/2;
    intsampl := x[k]/corr;
  end;

  function linint(n, k : integer; nu : double; var x : graphvect) : double;
  { Get phase by linear interpolation for rectangular window
    with -pi <= phi <= pi }
  var	i	: integer;
	phi	: double;
	xr, xi	: graphvect;
  begin
    for i:=1 to n do
    begin
      xr[i] := x[i]; xi[i] := 0.0;
    end;
    fft(n, xr, xi);
    phi := getangle(xr[k], xi[k]) - (n*nu-k+1)*pi;
    if phi > pi then
      phi := phi - 2.0*pi
    else if phi < -pi then
      phi := phi + 2.0*pi;
    linint := phi;
  end;

  procedure findres(n : integer; nux, nuy, f : double; var nx, ny : integer);
  { Match f by a linear combination of nux and nuy }
  var	eps	: double;
	found	: boolean;

  procedure fndres;
  var	i, j	: integer;
	delta	: double;
  begin
    for i:=0 to n do
    begin
      for j:=-n to n do
      begin
	delta := abs(i*nux+j*nuy);
	delta := delta - trunc(delta);
	if delta > 0.5 then delta := 1 - delta;
	delta := abs(delta-f);
	delta := delta - trunc(delta);
	if delta > 0.5 then delta := 1 - delta;
	if delta < eps then
	  if (abs(i)+abs(j) < n) and not ((i = 0) and (j < 0)) then
	  begin
	    found := true; nx := i; ny := j;
	  end;
      end;
    end;
  end;

  begin
    found := false; eps := 0.5e-6;
    repeat
      eps := 10*eps; fndres;
    until found ;
  end;

  procedure getpeaks(n : integer; var x : graphvect; nf : integer;
		     var nu, A : peakbuf);
  var	i, k, ind1, ind3	: integer;
  begin
    for i:=1 to nf do
    begin
      getpeak(n, x, k);
      nu[i] := int2snu(n, x, k); A[i] := intsampl(n, x, nu[i], k);
      { Make peak flat to allow for new call }
      getind(n, k, ind1, ind3);
      if x[ind1] > x[ind3] then
        x[k] := x[ind1]
      else
        x[k] := x[ind3];
    end;
  end;

  procedure settol(Fnum : integer; dxrms, dyrms, drrms : double);
  { Set random multipole displacement and roll errors }
  var	Knum	: integer;
	Cell	: celltype;
  begin
    for Knum:=1 to GetnKid(Fnum) do
    begin
      getElem(Elem_GetPos(Fnum, Knum), Cell);
      Cell.Elem.M.PdSrms[1] := dxrms; Cell.Elem.M.PdSrnd[1] := normranf;
      Cell.Elem.M.PdSrms[2] := dyrms; Cell.Elem.M.PdSrnd[2] := normranf;
      Cell.Elem.M.PdTrms := drrms;    Cell.Elem.M.PdTrnd := normranf;
      putElem(Elem_GetPos(Fnum, Knum), Cell);
      Mpole_SetdS(Fnum, Knum); Mpole_SetdT(Fnum, Knum);
    end;
  end;

  procedure setatol(Fnum, Knum : integer; dx, dy, dr : double);
  { Set a known random multipole displacement error }
  var	Cell	: celltype;
  begin
    getElem(Elem_GetPos(Fnum, Knum), Cell);
    Cell.Elem.M.PdSrms[1] := dx; Cell.Elem.M.PdSrnd[1] := 1d0;
    Cell.Elem.M.PdSrms[2] := dy; Cell.Elem.M.PdSrnd[2] := 1d0;
    Cell.Elem.M.PdTrms := dr;    Cell.Elem.M.PdTrnd := 1d0;
    putElem(Elem_GetPos(Fnum, Knum), Cell);
    Mpole_SetdS(Fnum, Knum); Mpole_SetdT(Fnum, Knum);
  end;

  procedure setkpar(Fnum, Knum, Order : integer; k : double);
  { Set multipole strength }
  var	Cell	: Celltype;
  begin
    getElem(Elem_GetPos(Fnum, Knum), Cell);
    Cell.Elem.M.PBpar[order] := k;
    putElem(Elem_GetPos(Fnum, Knum), Cell);
    Mpole_SetPB(Fnum, Knum, Order);
  end;

  procedure setkLpar(Fnum, Knum, Order : integer; kL : double);
  { Set integrated multipole strength }
  var	Cell	: Celltype;
  begin
    getElem(Elem_GetPos(Fnum, Knum), Cell);
    if Cell.Elem.PL <> 0d0 then
      Cell.Elem.M.PBpar[order] := kL/Cell.Elem.PL
    else
      Cell.Elem.M.PBpar[order] := kL;
    putElem(Elem_GetPos(Fnum, Knum), Cell);
    Mpole_SetPB(Fnum, Knum, Order);
  end;

  procedure setdkLpar(Fnum, Knum, Order : integer; dkL : double);
  { Increment integrated multipole strength }
  var	Cell	: Celltype;
  begin
    getElem(Elem_GetPos(Fnum, Knum), Cell);
    if Cell.Elem.PL <> 0d0 then
      Cell.Elem.M.PBpar[order] := Cell.Elem.M.PBpar[order] + dkL/Cell.Elem.PL
    else
      Cell.Elem.M.PBpar[order] := Cell.Elem.M.PBpar[order] + dkL;
    putElem(Elem_GetPos(Fnum, Knum), Cell);
    Mpole_SetPB(Fnum, Knum, Order);
  end;

  procedure setdkrpar(Fnum, Knum, Order : integer; dkrel : double);
  { Relative increment of multipole strength }
  var	Cell	: Celltype;
  begin
    getElem(Elem_GetPos(Fnum, Knum), Cell);
    if (order = dip) and (cell.elem.M.pthick = thick) then
      Cell.Elem.M.PBpar[dip] := Cell.Elem.M.PBpar[dip]
				+ dkrel*Cell.Elem.M.Pirho
    else
      Cell.Elem.M.PBpar[order] := Cell.Elem.M.PBpar[order]
			          + dkrel*Cell.Elem.M.PBpar[order];
    putElem(Elem_GetPos(Fnum, Knum), Cell);
    Mpole_SetPB(Fnum, Knum, Order);
  end;

  procedure setbn(fnum, order : integer; bn : double);
  var	i	: integer;
  begin
    for i:=1 to getnkid(fnum) do
      setkpar(fnum, i, order, bn);
  end;

  function getkpar(Fnum, Knum, Order : integer) : double;
  { Get multipole strength }
  var	Cell	: Celltype;
  begin
    getElem(Elem_GetPos(Fnum, Knum), Cell);
    getkpar := Cell.Elem.M.PBpar[order];
  end;

  function getkLpar(Fnum, Knum, Order : integer) : double;
  { Get integrated multipole strength }
  var	Cell	: Celltype;
  begin
    getElem(Elem_GetPos(Fnum, Knum), Cell);
    if Cell.Elem.PL <> 0d0 then
      getkLpar := Cell.Elem.M.PBpar[order]*Cell.Elem.PL
    else
      getkLpar := Cell.Elem.M.PBpar[order];
  end;

  procedure setdkLsys(Fnum, Order : integer; dkLsys : double);
  { Set systematic multipole errors }
  var	Knum	: integer;
	Cell	: Celltype;
  begin
    for Knum:=1 to GetnKid(Fnum) do
    begin
      getElem(Elem_GetPos(Fnum, Knum), Cell);
      if Cell.Elem.PL <> 0d0 then
        Cell.Elem.M.PBsys[order] := dkLsys/Cell.Elem.PL
      else
        Cell.Elem.M.PBsys[order] := dkLsys;
      putElem(Elem_GetPos(Fnum, Knum), Cell);
      Mpole_SetPB(Fnum, Knum, Order);
    end;
  end;

  procedure setdkLrms(Fnum, Order : integer; dkLrms : double);
  { Set random multipole errors }
  var	Knum	: integer;
	Cell	: Celltype;
  begin
    for Knum:=1 to GetnKid(Fnum) do
    begin
      getElem(Elem_GetPos(Fnum, Knum), Cell);
      if Cell.Elem.PL <> 0d0 then
        Cell.Elem.M.PBrms[order] := dkLrms/Cell.Elem.PL
      else
        Cell.Elem.M.PBrms[order] := dkLrms;
      Cell.Elem.M.PBrnd[order] := normranf;
      putElem(Elem_GetPos(Fnum, Knum), Cell);
      Mpole_SetPB(Fnum, Knum, Order);
    end;
  end;

  procedure setdkrrms(Fnum, Order : integer; dkrrms : double);
  { Set random multipole errors }
  var	Knum	: integer;
	Cell	: Celltype;
  begin
    for Knum:=1 to GetnKid(Fnum) do
    begin
      getElem(Elem_GetPos(Fnum, Knum), Cell);
      if (order = dip) and (cell.elem.M.pthick = thick) then
        Cell.Elem.M.PBrms[dip] := dkrrms*Cell.Elem.M.Pirho
      else
        Cell.Elem.M.PBrms[order] := dkrrms*Cell.Elem.M.PBpar[order];
      Cell.Elem.M.PBrnd[order] := normranf;
      putElem(Elem_GetPos(Fnum, Knum), Cell);
      Mpole_SetPB(Fnum, Knum, Order);
    end;
  end;

  procedure setkL(Fnum, Order : integer);
  { Sum up multipole strengths }
  var	Knum	: integer;
  begin
    for Knum:=1 to GetnKid(Fnum) do
      Mpole_SetPB(Fnum, Knum, Order);
  end;

  procedure setbpmdS(Fnum : integer; dxrms, dyrms : double);
  { Set bpm displacement errors }
  var	Knum	: integer;
	Cell	: Celltype;
  begin
    for Knum:=1 to GetnKid(Fnum) do
    begin
      getElem(Elem_GetPos(Fnum, Knum), Cell);
      Cell.dS[1] := normranf*dxrms; Cell.dS[2] := normranf*dyrms;
      putElem(Elem_GetPos(Fnum, Knum), Cell);
    end;
  end;

  { Routines for closed orbit correction }

  procedure codstat(var mean, sigma, xmax : vector2; lastpos : integer;
		    all : boolean);
  { Get statistics for  closed orbit }
  var	i, j, n        	: integer;
	sum, sum2	: vector2;
        Cell          	: celltype;
  begin
    n := 0;
    for j:=1 to 2 do
    begin
      sum[j] := 0.0; sum2[j] := 0.0; xmax[j] := 0.0;
    end;
    for i:=0 to lastpos do
    begin
      getElem(i, cell);
      if all or (cell.Fnum = globval.BPM) then
      begin
        n := succ(n);
	for j:=1 to 2 do
	begin
	  sum[j] := sum[j] + cell.beampos[2*j-1];
	  sum2[j] := sum2[j] + sqr(cell.beampos[2*j-1]);
	  xmax[j] := max_(xmax[j], abs(cell.beampos[2*j-1]));
	end;
      end;
    end;
    for j:=1 to 2 do
    begin
      if n <> 0 then
	mean[j] := sum[j]/n
      else
	mean[j] := 0.0;
      if (n <> 0) and (n <> 1) then
	sigma[j] := (n*sum2[j]-sqr(sum[j]))/(n*dble(n-1))
      else
	sigma[j] := 0.0;
      if sigma[j] >= 0.0 then
	sigma[j] := sqrt(sigma[j])
      else
	sigma[j] := 0.0;
    end;
  end;

  procedure gemit;
  {  Get equilibrium emittance }
  var	i, j, lastpos	: integer;
	U0		: double;
	eps, partnumb	: vector3;
	Ascr		: DAmap;
	emit, rad, cav	: boolean;
  begin
    getglobv_(globval);
    rad := globval.radiation; emit := globval.emittance;
    cav := globval.cavity_on;
    globval.radiation := true; globval.emittance := true;
    globval.cavity_on := true;
    putglobv_(globval);
    Ring_gettwiss(true, 0.0); getglobv_(globval);
    U0 := globval.dE*globval.energy;
    writeln; writeln('  energy loss per turn [keV]:    U0 =', abs(1d6*U0):10);
    for i:=1 to 6 do
    begin
      Ascr[i, 0] := globval.codvect[i];
      for j:=1 to 6 do
	Ascr[i, j] := globval.Ascr[i, j];
    end;
    cell_dapass(0, globval.cell_nloc, Ascr, lastpos); getglobv_(globval);
    for i:=1 to 3 do
    begin
      eps[i] := -globval.qfluct[i]/(4.0*globval.rad[i]);
      partnumb[i] := 2*globval.energy*(1d0+globval.codvect[5])*globval.rad[i]
		     /U0;
    end;
    globval.radiation := rad; globval.emittance := emit;
    globval.cavity_on := cav;
    putglobv_(globval);
    writeln('  partition numbers:             nx =', partnumb[1]:10,
	    ', ny =', partnumb[2]:10, ', nz =', partnumb[3]:10);
    writeln('  equilibrium emittance [m-rad]: ex =', eps[1]:10,
	    ', ey =', eps[2]:10, ', ez =', eps[3]:10);
  end;
