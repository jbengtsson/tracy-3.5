  { graphxphys.inc }

  { Routines for plotting }

  procedure printpeaks(xpos, ypos : double; nf : integer; var nu, A : peakbuf);
  var	i 	: integer;
	xe, ye	: double;
	str	: tstring;
  begin
    getstr(str, '#8#f         A');
    gratexs(xpos+1.4, ypos, labsmall, str.str, str.len, 0.0, xe, ye, 1);
    for i:=1 to nf do
    begin
      { do not supress leading blanks }
      str.str[1] := chr(esc); str.len := 1;
      concat(str, '#8#'); getreal(str, nu[i], 9, 5); getreale(str, A[i], 11, 2);
      gratexs(xpos, ypos-0.2-i*1.4*labsmall, labsmall, str.str, str.len,
	      0.0, xe, ye, 1);
    end;
  end;

  procedure plotdynap(r, delta, eps : double; npoint, n : integer;
		      floqs : boolean);
  var	x, y		: graphvect;
	str, xstr, ystr	: tstring;
  begin
    dynap(r, delta, eps, npoint, n, x, y, floqs, true);
    graope;
    gravwp(4.0, 2.0, 25.0, 20.0);
    str.str[1] := chr(esc); str.len := 1;
    if floqs then
      concat(str, '#8#  Dynamical aperture (Floquet space)')
    else
      concat(str, '#8#  Dynamical Aperture [m#6##L2#"#8##L1#m]');
    graheds(str.str, str.len, hedmedi);
    graxi1(loch, loch, -3, floch); grayi1(loch, loch, -3, floch);
    getstr(xstr, '#8#x^'); getstr(ystr, '#8#y^');
    gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len, 1);
    graph1(x, y, npoint, 0, gsquare, markbig); gracls;
  end;

  procedure plotps(n : integer; var x, px, y, py : graphvect);
  { Plot phase space }
  const	nf = 4; maxord = 5; prtpeaks = false;
  var	i			: integer;
	xe, ye			: double;
	nx, ny			: ipeakbuf;
	nux, nuy, A		: peakbuf;
	x1, px1, y1, py1	: graphvect;
	str, xstr, ystr		: tstring;
  begin
    for i:=1 to n do
    begin
      x1[i] := 1d3*x[i]; px1[i] := 1d3*px[i];
      y1[i] := 1d3*y[i]; py1[i] := 1d3*py[i];
    end;
    graope; gravwp(4.0, 14.5, 12.0, 8.0);
    getstr(str, '#8#Horizontal Phase Space');
    graheds(str.str, str.len, hedsmall);
    graxi1(loch, loch, loch, floch); grayi1(loch, loch, loch, floch);
    getstr(xstr, '#8#x [mm]'); getstr(ystr, '#8#p#IL#x');
    gracoms(labsmall, xstr.str, xstr.len, loch,  ystr.str, ystr.len, 1);
    graph1(x1, px1, n, 0, gcros, marksmall); graclx(0, 0, 0, 1);
    graope; gravwp(20.5, 14.5, 12.0, 8.0);
    getstr(str, '#8#FFT{x}');
    graheds(str.str, str.len , hedsmall);
    grayi1(loch, loch, loch, floch);
    getstr(xstr, '#6#n#8##IL#x'); getstr(ystr, '#8#A#IL#x'); 
    gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len , 1);
{    getmean(n, x1); }sinfft(n, x1);
    if prtpeaks then
    begin
      getpeaks(n, x1, nf, nux, A); printpeaks(21.5, 21.5, nf, nux, A);
    end;
    plotfft(n, x1); graclx(0, 0, 0, 1);
    graope; gravwp(4.0, 2.0, 12.0, 8.0);
    getstr(str, '#8#Vertical Phase Space');
    graheds(str.str, str.len, hedsmall );
    graxi1(loch, loch, loch, floch); grayi1(loch, loch, loch, floch);
    getstr(xstr, '#8#y [mm]'); getstr(ystr, '#8#p#IL#y');
    gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len, 1);
    graph1(y1, py1, n, 0, gcros, marksmall); graclx(0, 0, 0, 1);
    graope; gravwp(20.5, 2.0, 12.0, 8.0);
    getstr(str, '#8#FFT{y}');
    graheds(str.str, str.len, hedsmall );
    grayi1(loch, loch, loch, floch);
    getstr(xstr, '#6#n#8##IL#y'); getstr(ystr, '#8#A#IL#y'); 
    gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len, 1);
{    getmean(n, y1); }sinfft(n, y1);
    if prtpeaks then
    begin
      getpeaks(n, y1, nf, nuy, A); printpeaks(21.5, 9.0, nf, nuy, A);
    end;
    plotfft(n, y1); {graclx(0, 0, 0, 1);}
    if prtpeaks then
    begin
      for i:=2 to nf do
        findres(maxord, nux[1], nuy[1], nux[i], nx[i], ny[i]);
      getstr(str, '#8#n#IL#x#CL#  n#IL#y#CL#');
      gratexs(28.8, 21.5, labsmall, str.str, str.len, 0.0, xe, ye, 1);
      for i:=2 to nf do
      begin
        str.str[1] := chr(esc); str.len := 1;
        concat(str, '#8#'); getint(str, nx[i], 3); getint(str, ny[i], 3);
        gratexs(28.5+0.1, 21.5-0.2-i*1.4*labsmall, labsmall,
	        str.str, str.len, 0.0, xe, ye, 1);
      end;
      for i:=2 to nf do
        findres(maxord, nux[1], nuy[1], nuy[i], nx[i], ny[i]);
      getstr(str, '#8#n#IL#x#CL#  n#IL#y#CL#');
      gratexs(28.8, 9.0, labsmall, str.str, str.len, 0.0, xe, ye, 1);
      for i:=2 to nf do
      begin
        str.str[1] := chr(esc); str.len := 1;
        concat(str, '#8#'); getint(str, nx[i], 3); getint(str, ny[i], 3);
        gratexs(28.5+0.1, 9.0-0.2-i*1.4*labsmall, labsmall,
	        str.str, str.len, 0.0, xe, ye, 1);
      end;
    end;
    gracls;
  end;

  procedure plotj(n : integer; var x, y : graphvect);
  { Plot the variation of the linear invariant }
  const	nf = 4; prtpeaks = false;
  var	k		: integer;
	nux, nuy, A	: peakbuf;
	pos, Jx, Jy	: graphvect;
	str, xstr, ystr	: tstring;
  begin
    { Include 0 for plot }
    pos[n+1] := n+1; Jx[n+1] := 0.0; Jy[n+1] := 0.0;
    for k:=1 to n do
    begin
      pos[k] := k; Jx[k] := 1d6*x[k]; Jy[k] := 1d6*y[k];
    end;
    graope; gravwp(4.0, 14.5, 12.0, 8.0);
    getstr(str, '#8#J#IL#x#CL#(N) #6#[m#8#m#2##L2#4#4##L1#rad#6#]');
    graheds(str.str, str.len, hedsmall );
    grayi1(loch, loch, loch, floch);
    getstr(xstr, '#8#Turn no N');    getstr(ystr, '');
    gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len, 1);
    graph1(pos, Jx, n+1, 0, 0 , 0.0); plotbars(n, pos, Jx); graclx(0, 0, 0, 1);
    graope; gravwp(20.5, 14.5, 12.0, 8.0);
    str.str[1] := chr(esc); str.len := 1;
    concat(str, '  #8#FFT{J#IL#x#CL#}');
    graheds(str.str, str.len , hedsmall );
    grayi1(loch, loch, -3, floch);
    getstr(xstr, '#6#n#8##IL#x'); getstr(ystr, '#8#A#IL#x'); 
    gracoms(labsmall, xstr.str, xstr.len , loch, ystr.str, ystr.len , 1);
    getmean(n, Jx); sinfft(n, Jx);
    if prtpeaks then
    begin
      getpeaks(n, Jx, nf, nux, A); printpeaks(21.5, 21.5, nf, nux, A);
    end;
    plotfft( n, Jx); graclx(0, 0, 0, 1);
    graope; gravwp(4.0, 2.0, 12.0, 8.0);
    getstr(str, '#8#J#IL#y#CL#(N) #6#[m#8#m#2##L2#4#4##L1#rad#6#]');
    graheds(str.str, str.len, hedsmall );
    grayi1(loch, loch, loch, floch);
    getstr(xstr, '#8#Turn no N'); getstr(ystr, '');
    gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len, 1);
    graph1(pos, Jy, n+1, 0, 0 , 0.0); plotbars(n, pos, Jy); graclx(0, 0, 0, 1);
    graope; gravwp(20.5, 2.0, 12.0, 8.0);
    str.str[1] := chr(esc); str.len := 1;
    concat(str, '  #8#FFT{J#IL#y#CL#}');
    graheds(str.str, str.len , hedsmall );
    grayi1(loch, loch, -3, floch); 
    getstr(xstr, '#6#n#8##IL#y'); getstr(ystr, '#8#A#IL#y'); 
    gracoms(labsmall, xstr.str, xstr.len , loch, ystr.str, ystr.len , 1);
    getmean(n, Jy); sinfft(n, Jy);
    if prtpeaks then
    begin
      getpeaks(n, Jy, nf, nuy, A); printpeaks(21.5, 9.0, nf, nuy, A);
    end;
    plotfft(n, Jy); gracls;
  end;

  procedure plotphi(n : integer; var x, px, y, py : graphvect);
  { Plot the variation of the linear phase }
  const	nf = 4;
  var	k		: integer;
	pos		: graphvect;
	nux, nuy, A	: peakbuf;
	str, xstr, ystr	: tstring;
  begin
    getphi(n, x, px, y, py);
    { Include 0 for plot }
    pos[n+1] := n+1; x[n+1] := 0.0; y[n+1] := 0.0;
    for k:=1 to n do
      pos[k] := k;
    graope; gravwp(2.0, 14.0, 9.5, 9.0);
    str.str[1] := chr(esc); str.len := 1;
    concat(str, '#6##L2#"'); graheds(str.str, str.len, hedsmall );
    getstr(xstr, '#8#turn nb'); getstr(ystr, '#8#[rad]'); 
    gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len, 1);
    graph1(pos, x, n+1, 0, 0 , 0.0); plotbars(n, pos, x); graclx(0, 0, 0, 1);
    graope; gravwp(15.5, 14.0, 9.5, 9.0);
    str.str[1] := chr(esc); str.len := 1;
    concat(str, '#8#  FFT of #6##L2#"#8##IL#x#CL#');
    graheds(str.str, str.len, hedsmall );
    getstr(xstr, '#6#n'); getstr(ystr, '#8#A'); 
    gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len, 1);
    getmean(n, x); sinfft(n, x);
    getpeaks(n, x, nf, nux, A); printpeaks(21.5, 21.5, nf, nux, A);
    plotfft(n, x); graclx(0, 0, 0, 1);
    graope; gravwp(2.0, 1.5, 9.5, 9.0);
    str.str[1] := chr(esc); str.len := 1;
    concat(str, '#8#  #6##IL#y#CL#'); graheds(str.str, str.len, hedsmall );
    getstr(xstr, '#8#turn nb'); getstr(ystr, '#8#[rad]'); 
    gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len, 1);
    graph1(pos, y, n+1, 0, 0 , 0.0); plotbars(n, pos, y); graclx(0, 0, 0, 1);
    graope; gravwp(15.5, 1.5, 9.5, 9.0);
    str.str[1] := chr(esc); str.len := 1;
    getstr(str, '#8#  FFT of #6##IL#y#CL#');
    graheds(str.str, str.len, hedsmall );
    getstr(xstr, '#6#n'); getstr(ystr, '#8#A');
    gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len, 1);
    getmean(n, y); sinfft(n, y);
    getpeaks(n, y, nf, nuy, A); printpeaks(21.5, 9.0, nf, nuy, A);
    plotfft(n, y); gracls;
  end;

  procedure plotpos(lastpos : integer);
  { Plot a particle's position along lattice }
  var  i, ind		: integer;
       xstr, ystr, str	: tstring;
       s, px, py	: graphvect;
       Cell		: CellType;
  begin
    ind := 0;
    for i:=0 to lastpos do
    begin
      getelem(i, Cell);
      if Cell.Fnum = globval.BPM then
      begin
        ind := succ(ind); s[ind] := sngl(Cell.s);
        px[ind] := sngl(Cell.beampos[1]); py[ind] := sngl(Cell.beampos[3]);
      end;
    end;
    graope;
    gravwp(3.5, 14.0, 14.0, 9.0);
    str.str[1] := chr(esc); str.len := 1;
    concat(str, '#8#  Horizontal Position'); 
    graheds(str.str, str.len, hedsmall );
{    grayi1(loch, -1, loch, floch);}
    getstr(xstr, '#8#s [m]'); getstr(ystr, '#8#x [m]'); 
     gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len, 1);
    graph1(s, px, ind, 0, 0, 0.0); plotbars(ind, s, px); graclx(0, 0, 0, 1);
    graope;
    gravwp(3.5, 1.5, 14.0, 9.0);
    str.str[1] := chr(esc); str.len := 1;
    concat(str, '#8#  Vertical Position'); 
    graheds(str.str, str.len, hedsmall );
{    grayi1(loch, -1, loch, floch);}
    getstr(xstr, '#8#s [m]'); getstr(ystr, '#8#y [m]'); 
     gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len, 1);
    graph1(s, py, ind, 0, 0, 0.0); plotbars(ind, s, py); gracls;
  end;

  procedure plotcell(symfac : integer);
  { Plot alpha, beta and eta for one cell }
  const	etascl=10;
  var   i, n, loc, Fnum, Knum			: integer;
	ymin, ymax, ymin1, ymax1, ypos		: double;
	Tx, Tx1, Tx2, mark, xend, yend, b2, b3	: double;
        x, etax, px, py				: graphvect;
        str, xstr, ystr				: tstring;
        Cell   					: CellType;
  begin
    n := globval.Cell_nLoc div symfac + 1;
    for i:=1 to n do
    begin
      getelem(i-1, Cell);
      x[i] := sngl(Cell.s);
      px[i] := sngl(Cell.Beta[1]); py[i] := sngl(Cell.Beta[2]);
      etax[i] := sngl(etascl*Cell.Eta[1]);
    end;
    getminmax(n, px, ymin, ymax); getminmax(n, py, ymin1, ymax1);
    ymin := min_(ymin, ymin1); ymax := max_(ymax, ymax1);
    getminmax(n, etax, ymin1, ymax1);
    ymin := min_(ymin, ymin1); ymax := max_(ymax, ymax1);
    mark := (ymax-ymin)/10.0;
    ypos := 1.1*ymax+1.5*mark; ymax := 1.1*ymax + 3*mark;
    graope; gravwp(2.0, 2.0, 28.0, 18.0);
    getstr(str, '#8#Lattice Functions'); graheds(str.str, str.len, hedmedi );
    grawnd(x[1], x[i], ymin, ymax);
    getstr(xstr, '#8#s [m]'); getstr(ystr, ''); 
    gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len, 1);
    graph1(x, px, n, glsoli, 0 , 0.0); graph1(x, py, n, gldot, 0 , 0.0);
    graph1(x, etax, n, gldash, 0 , 0.0);
    getstr(str, '#6#n#8##IL#x#CL# ='); getreal(str, globval.totaltune[1], 9, 5);
    concat(str, ', #6#o#8##IL#x#CL# ='); getreal(str, globval.chrom[1], 9, 5);
    gratexs(12.0, 22.4, hedsmall, str.str, str.len, 0.0, xend, yend, 1);
    getstr(str, '#6#n#8##IL#y#CL# ='); getreal(str, globval.totaltune[2], 9, 5);
    concat(str, ', #6#o#8##IL#y#CL# ='); getreal(str, globval.chrom[2], 9, 5);
    gratexs(12.0, 21.0, hedsmall, str.str, str.len, 0.0, xend, yend, 1);
    for i:=0 to n-1 do
    begin
      getelem(i+1, Cell);
      px[4*i+1] := sngl(Cell.s - Cell.Elem.PL); py[4*i+1] := sngl(ypos);
      px[4*i+2] := sngl(px[4*i+1]); 
      if Cell.Elem.Pkind = Mpole then
      begin
        Elem_Getbend(Cell.Fnum, Cell.Knum, Tx, Tx1, Tx2);
	{ sign(0) = 0 for Fortran compability! }
        py[4*i+2] := ypos;
	if Tx <> 0.0 then
          py[4*i+2] := py[4*i+2] + sign(Tx)*mark/2.0;
	b2 := getkpar(cell.Fnum, cell.Knum, quad);
	if b2 <> 0.0 then
          py[4*i+2] := py[4*i+2] + sign(b2)*mark;
	b3 := getkpar(cell.Fnum, cell.Knum, sext);
	if b3 <> 0.0 then
          py[4*i+2] := py[4*i+2] + 1.5*sign(b3)*mark;
      end
      else
	py[4*i+2] := sngl(ypos);
      px[4*i+3] := px[4*i+1] + sngl(Cell.Elem.PL); py[4*i+3] := py[4*i+2];
      px[4*i+4] := px[4*i+3]; py[4*i+4] := py[4*i+1];
    end;
    graph1(px, py, n*4, glsoli, 0, 0.0); gracls;
  end;


  procedure plotcod(high : boolean; lastpos : integer; all, absolute : boolean);
  { Plot closed orbit }
  var	i, n			: integer;
	str, xstr, ystr		: tstring;
	nux, nuy, codx, cody	: graphvect;
        cell			: celltype;
        x, y                    : double;
  begin
    n := 0;
    for i:=0 to lastpos do
    begin
      getelem(i, cell);
      if all or (cell.Fnum = globval.BPM) then
      begin
        n := succ(n);
	nux[n] := cell.nu[1];      nuy[n] := cell.nu[2];
        x:=cell.beampos[1];        y:=cell.beampos[2];
        if absolute then begin
          x:=x+cell.ds[1];         y:=y+cell.ds[2];
        end;
	codx[n] := sngl(1d3*x);    cody[n] := sngl(1d3*y);
      end;
    end;
    graope;
    if high then
      gravwp(4.0, 15.0, 10.0, 8.5)
    else
      gravwp(4.0, 2.0, 10.0, 8.5);
    str.str[1] := chr(esc); str.len := 1;
    concat(str, '#8#Hor. Closed Orbit [mm]');
    graheds(str.str, str.len, hedsmall );
    getstr(xstr, '#6#n#8##IL#x#CL#'); getstr(ystr, '#8#x'); 
    gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len, 1);
    graph1(nux, codx, n, 0, 0, 0.0); plotbars(n, nux, codx); graclx(0, 0, 0, 1);
    graope;
    if high then
      gravwp(20.5, 15.0, 10.0, 8.5)
    else
      gravwp(20.5, 2.0, 10.0, 8.5);
    str.str[1] := chr(esc); str.len := 1;
    concat(str, '#8#Ver. Closed Orbit [mm]');
    graheds(str.str, str.len, hedsmall );
    getstr(xstr, '#6#n#8##IL#y#CL#'); getstr(ystr, '#8#y'); 
    gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len, 1);
    graph1(nuy, cody, n, 0, 0, 0.0); plotbars(n, nuy, cody);
    if high then
      graclx(0, 0, 0, 1)
    else
      gracls;
  end;

  procedure codcorrect(bumpimax : integer; thetamax : double);
  { Closed orbit correction using local bump method }
  const	all = false;
  var	i, lastpos		: integer;
        mean, sigma, xmax	: vector2;
        str, xstr, ystr		: tstring;
        n, px, py		: graphvect;
  begin
    getcod(0.0, lastpos);
    plotcod(true, lastpos, false, true);
    codstat(mean, sigma, xmax, lastpos, all);
    n[1] := 0; px[1] := sngl(1d3*sigma[1]); py[1] := sngl(1d3*sigma[2]);
    for i:=1 to bumpimax do
    begin
      execbump(thetamax, lastpos); getcod(0.0, lastpos);
      codstat(mean, sigma, xmax, lastpos, all);
      n[i+1] := i; px[i+1] := sngl(sigma[1]); py[i+1] := sngl(sigma[2]);
    end;
    graope; gravwp(4.0, 1.5, 28.0, 9.0);
    str.str[1] := chr(esc); str.len := 1;
    concat(str, '#8#  Horizontal and Vertical rms orbit');
    graheds(str.str, str.len, hedsmall);
    getstr(xstr, '#8#Iteration'); getstr(ystr, '#8#x, y [mm]'); 
    grawnd(0.0, sngl(bumpimax), 0.0, max_(px[1], py[1]));
    gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len, 1);
    graph1(n, px, bumpimax+1, 0, gcros, markbig);
    graph1(n, py, bumpimax+1, 0, gcirc, markbig);
    gracls;
  end;
