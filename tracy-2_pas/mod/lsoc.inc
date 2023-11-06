  { lsoc.inc }

  { Least-square closed orbit correction by singular value decomposition }

  procedure gcmat(bpm, corr, plane : integer; var A, V : svdmat;
		  var w : svdarray);

  { Get correlation matrix

		-----------
	      \/beta  beta
		    i     j
	A   = ------------- cos(nu pi - 2 pi|nu  - nu |)
	 ij   2 sin(pi nu)                     i     j         }

  var	i, j, m, n				: integer;
	nu, betai, betaj, nui, nuj, spiq	: double;
	cell					: celltype;
  begin
    m := getnkid(bpm); n := getnkid(corr);
    writeln;
    writeln(  '  no of BPMs:                  ', m:4);
    if plane = 1 then
      writeln('  no of horizontal correctors: ', n:4)
    else
      writeln('  no of vertical correctors:   ', n:4);
    nu := globval.totaltune[plane]; spiq := sin(pi*nu);
    for i:=1 to m do
    begin
      getelem(elem_getpos(bpm, i), cell);
      betai := cell.beta[plane]; nui := cell.nu[plane];
      for j:=1 to n do
      begin
	getelem(elem_getpos(corr, j), cell);
	betaj := cell.beta[plane]; nuj := cell.nu[plane];
	A[i, j] := sqrt(betai*betaj)/(2*spiq)*cos(nu*pi-abs(2*pi*(nui-nuj)));
      end;
    end;
    svdcmp(A, m, n, w, V);
{    writeln;
    for i:=1 to n do
    begin
      write('  w[', i:2, ']: ', w[i]:9);
      if i mod 4 = 0 then writeln;
    end;
    if n mod 4 <> 0 then writeln;}
  end;

  procedure lsoc(niter, bpm, corr, plane : integer; var A, V : svdmat;
		 var w : svdarray);
  const	all = false;
  var	i, j, m, n, lastpos	: integer;
	mean, xmax, sigma	: vector2;
	b, x			: svdarray;
	cell			: celltype;
  begin
    m := getnkid(bpm); n := getnkid(corr);
{    for i:=1 to n do
      if plane = 1 then
	setkpar(corr, i, dip, 0d0)
      else
	setkpar(corr, i, -dip, 0d0);}
    getcod(0d0, lastpos); codstat(mean, sigma, xmax, lastpos, all);
    writeln;
    if plane = 1 then
      writeln('  Horizontal orbit:           max [mm] =', 1d3*xmax[1]:9:5,
	      ', rms [mm] =', 1d3*sigma[1]:9:5)
    else
      writeln('  Vertical orbit:             max [mm] =', 1d3*xmax[2]:9:5,
	      ', rms [mm] =', 1d3*sigma[2]:9:5);
    for i:=1 to niter do
    begin
      for j:=1 to m do
      begin
        getelem(elem_getpos(bpm, j), cell); b[j] := -cell.beampos[2*plane-1];
      end;
      svbksb(A, w, V, m, n, b, x);
      for j:=1 to n do
	if plane = 1 then
	  setkpar(corr, j, dip, getkpar(corr, j, dip)-x[j])
	else
	  setkpar(corr, j, -dip, getkpar(corr, j, -dip)+x[j]);
      getcod(0d0, lastpos); codstat(mean, sigma, xmax, lastpos, all);
      if plane = 1 then
	writeln('  Corrected horizontal orbit: max [mm] =', 1d3*xmax[1]:9:5,
		', rms [mm] =', 1d3*sigma[1]:9:5)
      else
	writeln('  Corrected vertical orbit:   max [mm] =', 1d3*xmax[2]:9:5,
		', rms [mm] =', 1d3*sigma[2]:9:5);
    end;
  end;
