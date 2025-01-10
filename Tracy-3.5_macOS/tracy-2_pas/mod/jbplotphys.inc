  { plotphys.inc }

  { Routines for plotting }

  procedure plotfft(wn, n : integer; var x : graphvect);

  var	i, m	: integer;
	pos	: graphvect;

  begin
    m := n div 2 + 1;
    for i:=1 to m do
      pos[i] := (i-1)/n;
    getscale(m, pos, x, wn); plotbars(m, pos, x);
  end;


  procedure getnpeak(n : integer; x : graphvect; nq : integer);

  { Get frequency and amplitude of a peak by exact interpolation in
    the case of a sine window. }

  const	xpos=0.15; ypos=0.77;

  var	i, k, ind1, ind3	: integer;
	nu			: double;
	str			: tstring;

  begin
    gschh(0.03);
    getstr(str, '      f          A');
    gtxs(xpos, ypos, str.len, str.str);
    for i:=1 to nq do
    begin
      getpeak(n, x, k); nu := int2snu(n, x, k);
      str.len := 0; getreal(str, nu, 10, 5);
      getreale(str, intsampl(n, x, nu, k), 11, 2);
      gtxs(xpos, ypos-i*0.06, str.len, str.str);
      { Make peak flat to allow for new call }
      getind(n, k, ind1, ind3);
      if x[ind1] > x[ind3] then
        x[k] := x[ind1]
      else
        x[k] := x[ind3];
    end;
  end;


  procedure plotdynap(r, delta, eps : double; npoint, n : integer);

  var	i		: integer;
	phi		: double;
	x, y		: graphvect;
	str, xstr, ystr	: tstring;

  begin
    for i:=1 to npoint do
    begin
      phi := (i-1)*90d0/(npoint-1);
      getdynap(r, phi, 0d0, 0d0, delta, eps, n);
      x[i] := sngl(r*cos(dtor(phi))); y[i] := sngl(r*sin(dtor(phi)));
    end;

    gsvp(1, 0.0, 1.0, 0.1, 1.0); gswn(1, 0.0, 1.0, 0.0, 1.0);
    gselnt(1);
    getstr(str, 'Dynamical aperture'); gethead(str);
    getstr(xstr, 'x [m]'); getstr(ystr, 'y [m]'); getlabel(xstr, ystr);
    getscale(npoint, x, y, 1); gsmk(gsquare); gpm(npoint, x, y);
  end;


  procedure plotps;

  { Plot phase space }

  const	nq=5;

  var	n		: integer;
	x, px, y, py	: graphvect;
	str, xstr, ystr	: tstring;

  begin
    gettrack(n, x, px, y, py);

    gsmk(gplus);

    gsvp(1, 0.0, 0.5, 0.55, 1.0); gswn(1, 0.0, 1.0, 0.0, 1.0); gselnt(1);
    getstr(str, 'Hor. phase space'); gethead(str);
    getstr(xstr, 'x'); getstr(ystr, 'px'); getlabel(xstr, ystr);
    getscale(n, x, px, 1); gpm(n, x, px);

    gsvp(2, 0.5, 1.0, 0.55, 1.0); gswn(2, 0.0, 1.0, 0.0, 1.0); gselnt(2);
    getstr(str, 'FFT of hor. pos.'); gethead(str);
    getstr(xstr, 'nu'); getstr(ystr, 'A'); getlabel(xstr, ystr);
    getmean(n, x); sinfft(n, x); getnpeak(n, x, nq); plotfft(2, n, x);

    gsvp(3, 0.0, 0.5, 0.1, 0.55); gswn(3, 0.0, 1.0, 0.0, 1.0); gselnt(3);
    getstr(str, 'Ver. phase space'); gethead(str);
    getstr(xstr, 'y'); getstr(ystr, 'py'); getlabel(xstr, ystr);
    getscale(n, y, py, 3); gpm(n, y, py);

    gsvp(4, 0.5, 1.0, 0.1, 0.55); gswn(4, 0.0, 1.0, 0.0, 1.0); gselnt(4);
    getstr(str, 'FFT of ver. position'); gethead(str);
    getstr(xstr, 'nu'); getstr(ystr, 'A'); getlabel(xstr, ystr);
    getmean(n, y); sinfft(n, y); getnpeak(n, y, nq); plotfft(4, n, y);
  end;


  procedure plotj;

  { Plot the variation of the linear invariant }

  const	nq=5;

  var	n, k			: integer;
	x, px, y, py, pos	: graphvect;
	str, xstr, ystr		: tstring;

  begin
    gettrack(n, x, px, y, py); getj(n, x, px, y, py);

    { Include 0 for plot }
    pos[n+1] := n+1; x[n+1] := 0.0; y[n+1] := 0.0;

    for k:=1 to n do
      pos[k] := k;

    gsvp(1, 0.0, 0.5, 0.55, 1.0); gswn(1, 0.0, 1.0, 0.0, 1.0); gselnt(1);
    getstr(str, 'Ix'); gethead(str);
    getstr(xstr, 'turn nb'); ystr.len := 0; getlabel(xstr, ystr);
    getscale(n+1, pos, x, 1); plotbars(n, pos, x);

    gsvp(2, 0.5, 1.0, 0.55, 1.0); gswn(2, 0.0, 1.0, 0.0, 1.0); gselnt(2);
    getstr(str, 'FFT of Ix'); gethead(str);
    getstr(xstr, 'nu'); getstr(ystr, 'A'); getlabel(xstr, ystr);
    getmean(n, x); sinfft(n, x); getnpeak(n, x, nq); plotfft(2, n, x);

    gsvp(3, 0.0, 0.5, 0.1, 0.55); gswn(3, 0.0, 1.0, 0.0, 1.0); gselnt(3);
    getstr(str, 'Iy'); gethead(str);
    getstr(xstr, 'turn nb'); ystr.len := 0; getlabel(xstr, ystr);
    getscale(n+1, pos, y, 3); plotbars(n, pos, y);

    gsvp(4, 0.5, 1.0, 0.1, 0.55); gswn(4, 0.0, 1.0, 0.0, 1.0); gselnt(4);
    getstr(str, 'FFT of Iy'); gethead(str);
    getstr(xstr, 'nu'); getstr(ystr, 'A'); getlabel(xstr, ystr);
    getmean(n, y); sinfft(n, y); getnpeak(n, y, nq); plotfft(4, n, y);
  end;


  procedure plotphi;

  { Plot the variation of the linear phase }

  const	nq=5;

  var	n, k			: integer;
	x, px, y, py, pos	: graphvect;
	str, xstr, ystr		: tstring;

  begin
    gettrack(n, x, px, y, py); getphi(n, x, px, y, py);

    { Include 0 for plot }
    pos[n+1] := n+1; x[n+1] := 0.0; y[n+1] := 0.0;

    for k:=1 to n do
      pos[k] := k;

    gsvp(1, 0.0, 0.5, 0.55, 1.0); gswn(1, 0.0, 1.0, 0.0, 1.0); gselnt(1);
    getstr(str, 'phix'); gethead(str);
    getstr(xstr, 'turn nb'); getstr(ystr, '[rad]'); getlabel(xstr, ystr);
    getscale(n+1, pos, x, 1); plotbars(n, pos, x);

    gsvp(2, 0.5, 1.0, 0.55, 1.0); gswn(2, 0.0, 1.0, 0.0, 1.0); gselnt(2);
    getstr(str, 'FFT of phix'); gethead(str);
    getstr(xstr, 'nu'); getstr(ystr, 'A'); getlabel(xstr, ystr);
    getmean(n, x); sinfft(n, x); getnpeak(n, x, nq); plotfft(2, n, x);

    gsvp(3, 0.0, 0.5, 0.1, 0.55); gswn(3, 0.0, 1.0, 0.0, 1.0); gselnt(3);
    getstr(str, 'phiy'); gethead(str);
    getstr(xstr, 'turn nb'); getstr(ystr, '[rad]'); getlabel(xstr, ystr);
    getscale(n+1, pos, y, 3); plotbars(n, pos, y);

    gsvp(4, 0.5, 1.0, 0.1, 0.55); gswn(4, 0.0, 1.0, 0.0, 1.0); gselnt(4);
    getstr(str, 'FFT of phiy'); gethead(str);
    getstr(xstr, 'nu'); getstr(ystr, 'A'); getlabel(xstr, ystr);
    getmean(n, y); sinfft(n, y); getnpeak(n, y, nq); plotfft(4, n, y);
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
      getElem(i, Cell);
      if Cell.Fnum = globval.BPM then
      begin
        ind := succ(ind);
        s[ind] := sngl(Cell.s);
        px[ind] := sngl(Cell.beampos[1]);
	py[ind] := sngl(Cell.beampos[3]);
      end;
    end;

    gsvp(1, 0.0, 0.5, 0.55, 1.0); gswn(1, 0.0, 1.0, 0.0, 1.0); gselnt(1);
    getstr(str, 'Hor. position'); gethead(str);
    getstr(xstr, 's [m]'); getstr(ystr, 'x [m]'); getlabel(xstr, ystr);
    getscale(ind, s, px, 1); plotbars(ind, s, px);

    gsvp(3, 0.0, 0.5, 0.1, 0.55); gswn(3, 0.0, 1.0, 0.0, 1.0); gselnt(3);
    getstr(str, 'Ver. position'); gethead(str);
    getstr(xstr, 's [m]'); getstr(ystr, 'y [m]'); getlabel(xstr, ystr);
    getscale(ind, s, py, 3); plotbars(ind, s, py);
  end;


  procedure plotcell(symfac : integer);

  { Plot alpha, beta and eta for one cell }

  const	eps=1d-3; etascl=100;

  var   i, Fnum, Knum	 			: integer;
	circ, ymin, ymax, ymin1, ymax1, ypos	: double;
	Tx, Tx1, Tx2, mark			: double;
	rxmin, rxmax, rymin, rymax		: double;
        x, etax, px, py				: graphvect;
        str, xstr, ystr				: tstring;
        Cell   					: CellType;

  FUNCTION sign(x : double) : integer;

  { special implementation of sign }

  BEGIN
    if x = 0d0 then
      sign := 0
    else if x > 0d0 then
      sign := 1
    else
      sign := -1;
  END;

  begin
    getElem(globval.Cell_nLoc, Cell);
    circ := Cell.s; i := 0; getElem(i, Cell);
    while (Cell.s < circ/symfac+eps) and (i <= globval.Cell_nLoc) do
    begin
      x[i+1] := sngl( Cell.s );
      px[i+1] := sngl( Cell.Beta[1] );
      py[i+1] := sngl( Cell.Beta[2] );
      etax[i+1] := sngl( etascl*Cell.Eta[1] );
      i := succ(i);
      getElem(i, Cell);
    end;

    getminmax(i, px, ymin, ymax);
    getminmax(i, py, ymin1, ymax1);
    ymin := min_(ymin, ymin1); ymax := max_(ymax, ymax1);
    getminmax(i, etax, ymin1, ymax1);
    ymin := min_(ymin, ymin1); ymax := max_(ymax, ymax1);
    mark := (ymax-ymin)/10.0;
    ypos := 1.1*ymax+1.5*mark; ymax := 1.1*ymax + 3*mark;

    gsvp(1, 0.0, 1.0, 0.1, 1.0); gswn(1, 0.0, 1.0, 0.0, 1.0); gselnt(1);
{    getstr(str, 'Lattice functions'); gethead(str);
    getstr(xstr, 's [m]'); ystr.len := 0; getlabel(xstr, ystr);}
    rxmin := sngl(x[1]); rxmax := sngl(x[i]);
    rymin := sngl(ymin); rymax := sngl(ymax);
    cframe(1, rxmin, rxmax, rymin, rymax);
    gsln(gldash); gpl(i, x, px);
    gsln(gldot); gpl(i, x, py);
    gsln(gldasd);  gpl(i, x, etax);

    i := 0;
    for Fnum:=1 to globval.Elem_nFam do
      for Knum:=1 to GetnKid(Fnum) do
      begin
        getElem(Elem_Getpos(Fnum, Knum), Cell);
	if (Cell.s < circ/symfac+eps) and (i <= globval.Cell_nLoc) then
	begin
	  px[4*i+1] := sngl(Cell.s - Cell.Elem.PL); py[4*i+1] := sngl(ypos);
          px[4*i+2] := sngl(px[4*i+1]); 
	  if (Cell.Elem.Pkind = Mpole) or (Cell.Elem.Pkind = Wigl) then
	  begin
	    Elem_Getbend(Fnum, Knum, Tx, Tx1, Tx2);
            py[4*i+2] := sngl(ypos
	              + sign(Tx)*mark/2.0
	              + sign(Elem_Getkval(Fnum, Knum, Quad))*mark
	              + sign(Elem_Getkval(Fnum, Knum, Sext))*mark);
	  end
	  else
	    py[4*i+2] := sngl(ypos);
          px[4*i+3] := px[4*i+1] + sngl(Cell.Elem.PL); py[4*i+3] := py[4*i+2];
          px[4*i+4] := px[4*i+3]; py[4*i+4] := py[4*i+1];
          i := succ(i);
	end;
      end;
    gsln(glsoli); gpl(4*i, px, py);
  end;


  procedure plotcell1(symfac : integer);

  { Plot alpha, beta and eta for one cell }

  const	eps = 1d-3;

  var   n	 	: integer;
	circ		: double;
        x, etax, px, py	: graphvect;
        str, xstr, ystr	: tstring;
        Cell   		: CellType;

  begin
    getElem(globval.Cell_nLoc, Cell);
    circ := Cell.s; n := 0; getElem(n, Cell);
    while (Cell.s < circ/symfac+eps) and (n <= globval.Cell_nLoc) do
    begin
      x[n+1] := sngl( Cell.s );
      px[n+1] := sngl( Cell.Beta[1] );
      py[n+1] := sngl( Cell.Beta[2] );
      etax[n+1] := sngl( Cell.Eta[1] );
      n := succ(n);
      getElem(n, Cell);
    end;

    gsvp(1, 0.0, 1.0, 0.55, 1.0); gswn(1, 0.0, 1.0, 0.0, 1.0); gselnt(1);
{    getstr(str, 'Lattice functions'); gethead(str);
    getstr(xstr, 's [m]'); ystr.len := 0; getlabel(xstr, ystr);}
    getscale(n, x, px, 1);
    gsln(gldash); gpl(n, x, px);
    gsvp(2, 0.0, 1.0, 0.1, 0.55); gswn(2, 0.0, 1.0, 0.0, 1.0); gselnt(2);
    getscale(n, x, py, 2);
    gsln(gldot); gpl(n, x, py);
{    gsln(gldasd);  gpl(n, x, etax);}
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
      codx[i] := sngl(Cell.beampos[1] + Cell.ds[1]);
      cody[i] := sngl(Cell.beampos[3] + Cell.ds[2]);
    end;

    gsvp(1, 0.0, 0.5, 0.55, 1.0); gswn(1, 0.0, 1.0, 0.0, 1.0); gselnt(1);
    getstr(str, 'Hor. orbit'); gethead(str);
    getstr(xstr, 's [m]'); getstr(ystr, 'x [m]'); getlabel(xstr, ystr);
    getscale(n, x, codx, 1); plotbars(n, x, codx);

    gsvp(2, 0.5, 1.0, 0.55, 1.0); gswn(2, 0.0, 1.0, 0.0, 1.0); gselnt(2);
    getstr(str, 'Ver. orbit'); gethead(str);
    getstr(xstr, 's [m]'); getstr(ystr, 'y [m]'); getlabel(xstr, ystr);
    getscale(n, x, cody, 2); plotbars(n, x, cody);
  end;


  procedure codcorrect(bumpimax : integer; thetamax : double);

  { Closed orbit correction using local bump method }

  var	i, lastpos			: integer;
        mean, sigma, xmax		: double;
	rxmin, rxmax, rymin, rymax	: double;
        str, xstr, ystr			: tstring;
        px, py				: graphvect;

  begin
    gsmk(gast);

    getcod(0.0, lastpos); plotcod(lastpos);

    codstat(1, mean, sigma, xmax, lastpos);

    gsvp(3, 0.0, 0.5, 0.1, 0.55); gswn(3, 0.0, 1.0, 0.0, 1.0); gselnt(3);
    getstr(str, 'Hor. rms orbit'); gethead(str);
    getstr(xstr, 'iteration'); getstr(ystr, 'x rms [m]'); getlabel(xstr, ystr);

    rxmin := sngl(0.0); rxmax := sngl(bumpimax);
    rymin := sngl(0.0); rymax := sngl(sigma);
    cframe(3, rxmin, rxmax, rymin, rymax);
    px[1] := 0; py[1] := sngl(sigma); gpm(1, px, py);

    codstat(2, mean, sigma, xmax, lastpos);

    gsvp(4, 0.5, 1.0, 0.1, 0.55); gswn(4, 0.0, 1.0, 0.0, 1.0); gselnt(4);
    getstr(str, 'Ver. rms orbit'); gethead(str);
    getstr(xstr, 'iteration'); getstr(ystr, 'y rms [m]'); getlabel(xstr, ystr);

    rxmin := sngl(0.0); rxmax := sngl(bumpimax);
    rymin := sngl(0.0); rymax := sngl(sigma);
    cframe(4, rxmin, rxmax, rymin, rymax);
    px[1] := 0; py[1] := sngl(sigma); gpm(1, px, py);

    for i:=1 to bumpimax do
    begin
      execbump(thetamax, lastpos); getcod(0.0, lastpos);

      codstat(1, mean, sigma, xmax, lastpos);
      gselnt(3); px[1] := i; py[1] := sngl(sigma); gpm(1, px, py);
      codstat(2, mean, sigma, xmax, lastpos);
      gselnt(4); px[1] := i; py[1] := sngl(sigma); gpm(1, px, py);
    end;
  end;
