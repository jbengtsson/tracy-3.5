  { Routines for plotting }

  const	esc = 27;

  procedure printpeaks(xpos, ypos : double; nf : integer; var nu, A : peakbuf);
  var	i 	: integer;
	xe, ye	: double;
	str	: tstring;
  begin
    getstr(str, 'f         A');
    gratexs(xpos+1.4, ypos, labsmall, str.str, str.len, 0.0, xe, ye, 1);
    for i:=1 to nf do
    begin
      getstr(str, ''); getreal(str, nu[i], 9, 5); getreale(str, A[i], 11, 2);
      { do not supress leading blanks }
      str.str[1] := chr(esc);
      gratexs(xpos, ypos-0.2-i*1.4*labsmall, labsmall, str.str, str.len,
	      0.0, xe, ye, 1);
    end;
  end;

  procedure plotdynap(r, delta, eps : double; npoint, n : integer);
  var	x, y		: graphvect;
	str, xstr, ystr	: tstring;
  begin
    dynap(r, delta, eps, npoint, n, x, y, true);
    graope;
    gravwp(4.5, 2.5, 25.0, 20.0);
    getstr(str, 'Dynamical aperture (Floquet space)');
    graheds(str.str, str.len, hedmedi);
    graxi1(loch, loch, -3, floch); grayi1(loch, loch, -3, floch);
    getstr(xstr, 'x^'); getstr(ystr, 'y^');
    gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len, 1);
    graph1(x, y, npoint, 0, gsquare, markbig); gracls;
  end;

  procedure plotps(n : integer; var x, px, y, py : graphvect);
  { Plot phase space }
  const	nf = 4; maxord = 5;
  var	i		: integer;
	xe, ye		: double;
	nx, ny		: ipeakbuf;
	nux, nuy, A	: peakbuf;
	str, xstr, ystr	: tstring;
  begin
    graope; gravwp(4.5, 15.0, 8.5, 7.5);
    getstr(str, 'Horizontal phase space'); graheds(str.str, str.len, hedsmall);
    graxi1(loch, loch, -3, floch); grayi1(loch, loch, -3, floch);
    getstr(xstr, 'x [m]'); getstr(ystr, 'px [rad]'); 
    gracoms(labsmall, xstr.str, xstr.len, loch,  ystr.str, ystr.len, 1);
    graph1(x, px, n, 0, gcros, marksmall); graclx(0, 0, 0, 1);
    graope; gravwp(21.0, 15.0, 8.5, 7.5);
    getstr(str, 'FFT of horizontal position');
    graheds(str.str, str.len , hedsmall);
    grayi1(loch, loch, -3, floch);
    getstr(xstr, 'nu'); getstr(ystr, 'A [m]'); 
    gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len , 1);
{    getmean(n, x); }sinfft(n, x);
    getpeaks(n, x, nf, nux, A); printpeaks(25.0, 22.5, nf, nux, A);
    plotfft(n, x); graclx(0, 0, 0, 1);
    graope; gravwp(4.5, 2.5, 11.5, 7.5);
    getstr(str, 'Vertical phase space'); graheds(str.str, str.len, hedsmall );
    graxi1(loch, loch, -3, floch); grayi1(loch, loch, -3, floch);
    getstr(xstr, 'y [m]'); getstr(ystr, 'py [rad]');
    gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len, 1);
    graph1(y, py, n, 0, gcros, marksmall); graclx(0, 0, 0, 1);
    graope; gravwp(21.5, 2.5, 11.5, 7.5);
    getstr(str, 'FFT of vertical position');
    graheds(str.str, str.len, hedsmall );
    grayi1(loch, loch, -3, floch);
    getstr(xstr, 'nu'); getstr(ystr, 'A [m]'); 
    gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len, 1);
{    getmean(n, y); }sinfft(n, y);
    getpeaks(n, y, nf, nuy, A); printpeaks(25.0, 10.0, nf, nuy, A);
    plotfft(n, y); {graclx(0, 0, 0, 1);}
    for i:=2 to nf do
      findres(maxord, nux[1], nuy[1], nux[i], nx[i], ny[i]);
    getstr(str, 'n#IL#x#CL#  n#IL#y#CL#');
    gratexs(30.8, 22.5, labsmall, str.str, str.len, 0.0, xe, ye, 1);
    for i:=2 to nf do
    begin
      getstr(str, ''); getint(str, nx[i], 3); getint(str, ny[i], 3);
      { do not supress leading blanks }
      str.str[1] := chr(esc);
      gratexs(30.5+0.1, 22.5-0.2-i*1.4*labsmall, labsmall,
	      str.str, str.len, 0.0, xe, ye, 1);
    end;
    for i:=2 to nf do
      findres(maxord, nux[1], nuy[1], nuy[i], nx[i], ny[i]);
    getstr(str, 'n#IL#x#CL#  n#IL#y#CL#');
    gratexs(30.8, 10.0, labsmall, str.str, str.len, 0.0, xe, ye, 1);
    for i:=2 to nf do
    begin
      getstr(str, ''); getint(str, nx[i], 3); getint(str, ny[i], 3);
      { do not supress leading blanks }
      str.str[1] := chr(esc);
      gratexs(30.5+0.1, 10.0-0.2-i*1.4*labsmall, labsmall,
	      str.str, str.len, 0.0, xe, ye, 1);
    end;
    gracls;
  end;

  procedure plotj(n : integer; var x, px, y, py : graphvect);
  { Plot the variation of the linear invariant }
  const	nf = 4;
  var	k		: integer;
	nux, nuy, A	: peakbuf;
	pos		: graphvect;
	str, xstr, ystr	: tstring;
  begin
    getj(n, x, px, y, py);
    { Include 0 for plot }
    pos[n+1] := n+1; x[n+1] := 0.0; y[n+1] := 0.0;
    for k:=1 to n do
      pos[k] := k;
    graope; gravwp(2.5, 14.0, 9.5, 9.0);
    getstr(str, 'Ix [m rad]'); graheds(str.str, str.len, hedsmall );
    grayi1(loch, loch, -6, floch);
    getstr(xstr, 'turn nb');    getstr(ystr, ' ');
    gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len, 1);
    graph1(pos, x, n+1, 0, 0 , 0.0); plotbars(n, pos, x); graclx(0, 0, 0, 1);
    graope; gravwp(15.5, 14.0, 9.5, 9.0);
    getstr(str, 'FFT of Ix'); graheds(str.str, str.len , hedsmall );
    grayi1(loch, loch, -6, floch);
    getstr(xstr, 'nu'); getstr(ystr, 'A [m rad]'); 
    gracoms(labsmall, xstr.str, xstr.len , loch, ystr.str, ystr.len , 1);
    getmean(n, x); sinfft(n, x);
    getpeaks(n, x, nf, nux, A); printpeaks(25.0, 22.5, nf, nux, A);
    plotfft( n, x); graclx(0, 0, 0, 1);
    graope; gravwp(2.5, 1.5, 9.5, 9.0);
    getstr(str, 'Iy [m rad]'); graheds(str.str, str.len, hedsmall );
    grayi1(loch, loch, -6, floch);
    getstr(xstr, 'turn nb'); getstr(ystr, ' ');
    gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len, 1);
    graph1(pos, y, n+1, 0, 0 , 0.0); plotbars(n, pos, y); graclx(0, 0, 0, 1);
    graope; gravwp(15.5, 1.5, 9.5, 9.0);
    getstr(str, 'FFT of Iy'); graheds(str.str, str.len , hedsmall );
    grayi1(loch, loch, -6, floch); 
    getstr(xstr, 'nu'); getstr(ystr, 'A [m rad]'); 
    gracoms(labsmall, xstr.str, xstr.len , loch, ystr.str, ystr.len , 1);
    getmean(n, y); sinfft(n, y);
    getpeaks(n, y, nf, nuy, A); printpeaks(25.0, 10.0, nf, nuy, A);
    plotfft(n, y); gracls;
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
    graope; gravwp(2.5, 14.0, 9.5, 9.0);
    getstr(str, 'phix'); graheds(str.str, str.len, hedsmall );
    getstr(xstr, 'turn nb'); getstr(ystr, '[rad]'); 
    gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len, 1);
    graph1(pos, x, n+1, 0, 0 , 0.0); plotbars(n, pos, x); graclx(0, 0, 0, 1);
    graope; gravwp(15.5, 14.0, 9.5, 9.0);
    getstr(str, 'FFT of phix'); graheds(str.str, str.len, hedsmall );
    getstr(xstr, 'nu'); getstr(ystr, 'A'); 
    gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len, 1);
    getmean(n, x); sinfft(n, x);
    getpeaks(n, x, nf, nux, A); printpeaks(25.0, 22.5, nf, nux, A);
    plotfft(n, x); graclx(0, 0, 0, 1);
    graope; gravwp(2.5, 1.5, 9.5, 9.0);
    getstr(str, 'phiy'); graheds(str.str, str.len, hedsmall );
    getstr(xstr, 'turn nb'); getstr(ystr, '[rad]'); 
    gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len, 1);
    graph1(pos, y, n+1, 0, 0 , 0.0); plotbars(n, pos, y); graclx(0, 0, 0, 1);
    graope; gravwp(15.5, 1.5, 9.5, 9.0);
    getstr(str, 'FFT of phiy'); graheds(str.str, str.len, hedsmall );
    getstr(xstr, 'nu'); getstr(ystr, 'A'); 
    gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len, 1);
    getmean(n, y); sinfft(n, y);
    getpeaks(n, y, nf, nuy, A); printpeaks(25.0, 10.0, nf, nuy, A);
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
    getstr(str, 'Hor. position'); 
    graheds(str.str, str.len, hedsmall );
{    grayi1(loch, -1, loch, floch);}
    getstr(xstr, 's [m]'); getstr(ystr, 'x [m]'); 
     gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len, 1);
    graph1(s, px, ind, 0, 0, 0.0); plotbars(ind, s, px); graclx(0, 0, 0, 1);
    graope;
    gravwp(3.5, 1.5, 14.0, 9.0);
    getstr(str, 'Ver. position'); 
    graheds(str.str, str.len, hedsmall );
{    grayi1(loch, -1, loch, floch);}
    getstr(xstr, 's [m]'); getstr(ystr, 'y [m]'); 
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
    n := globval.Cell_nLoc div symfac;
    for i:=1 to n do
    begin
      getelem(i, Cell);
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
    getstr(str, 'Lattice functions'); graheds(str.str, str.len, hedmedi );
    grawnd(x[1], x[i], ymin, ymax);
    getstr(xstr, 's [m]'); getstr(ystr, ' '); 
    gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len, 1);
    graph1(x, px, n, 1, 0 , 0.0); graph1(x, py, n, 2, 0 , 0.0);
    graph1(x, etax, n, 3, 0 , 0.0);
    getstr(str, '#2#n#1##IL#x#CL# ='); getreal(str, globval.totaltune[1], 9, 5);
    concat(str, ', #2#o#1##IL#x#CL# ='); getreal(str, globval.chrom[1], 9, 5);
    gratexs(12.0, 22.4, hedmedi, str.str, str.len, 0.0, xend, yend, 1);
    getstr(str, '#2#n#1##IL#y#CL# ='); getreal(str, globval.totaltune[2], 9, 5);
    concat(str, ', #2#o#1##IL#y#CL# ='); getreal(str, globval.chrom[2], 9, 5);
    gratexs(12.0, 21.0, hedmedi, str.str, str.len, 0.0, xend, yend, 1);
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


  procedure plotcod(lastpos : integer);
  { Plot closed orbit }
  var	i, n		: integer;
	str, xstr, ystr	: tstring;
	x, codx, cody	: graphvect;
        Cell		: CellType;
  begin
    n := getnkid(globval.bpm);
    for i:=1 to n do
    begin
      getElem(elem_getpos(globval.bpm, i), Cell);
      x[i] := i;
      codx[i] := sngl(1d3*Cell.beampos[1]+Cell.ds[1]);
      cody[i] := sngl(1d3*Cell.beampos[3]+Cell.ds[2]);
    end;
    graope; gravwp(4.5, 14.0, 8.5, 9.0);
    getstr(str, 'Horizontal position');
    graheds(str.str, str.len, hedsmall );
    getstr(xstr, 's [m]'); getstr(ystr, 'x [mm]'); 
    gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len, 1);
    graph1(x, codx, n, 0, 0, 0.0); plotbars(n, x, codx); graclx(0, 0, 0, 1);
    graope; gravwp(21.0, 14.0, 8.5, 9.0);
    getstr(str, 'Vertical position'); graheds(str.str, str.len, hedsmall );
    getstr(xstr, 's [m]'); getstr(ystr, 'y [mm]'); 
    gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len, 1);
    graph1(x, cody, n, 0, 0, 0.0); plotbars(n, x, cody); graclx(0, 0, 0, 1);
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
    plotcod(lastpos);
    codstat(mean, sigma, xmax, lastpos, all);
    n[1] := 0; px[1] := sngl(1d3*sigma[1]); py[1] := sngl(1d3*sigma[2]);
    for i:=1 to bumpimax do
    begin
      execbump(thetamax, lastpos); getcod(0.0, lastpos);
      codstat(mean, sigma, xmax, lastpos, all);
      n[i+1] := i; px[i+1] := sngl(sigma[1]); py[i+1] := sngl(sigma[2]);
    end;
    graope; gravwp(4.5, 1.5, 28.0, 9.0);
    getstr(str, 'Horizontal and Vertical rms orbit');
    graheds(str.str, str.len, hedsmall);
    getstr(xstr, 'Iteration'); getstr(ystr, 'x, y [mm]'); 
    grawnd(0.0, sngl(bumpimax), 0.0, max_(px[1], py[1]));
    gracoms(labsmall, xstr.str, xstr.len, loch, ystr.str, ystr.len, 1);
    graph1(n, px, bumpimax+1, 0, gcros, markbig);
    graph1(n, py, bumpimax+1, 0, gcirc, markbig);
    gracls;
  end;
