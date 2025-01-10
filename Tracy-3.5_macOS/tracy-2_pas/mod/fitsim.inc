  const	nmax = 20;

  type	intbuf	= array [1..nmax] of integer;
	realbuf	= array [1..nmax] of double;


  procedure fitlat(ftol : double; niter, nprint : integer;
		   nknob : integer; var Fknob, Kknob : intbuf;
		   var knobs : realbuf;
		   nset : integer; var Floc, Kloc : intbuf; 
		   var alphax, betax, etax, etaxp,
		   alphay, betay : realbuf);

  { Fitt lattice functions by multidimensional Newton search with
    singular value decomposition }

  const	nparam = 6;

  var	i, j, iter, cut	: integer;
	eps		: double;
	Cell		: CellType;
	first		: boolean;

  procedure nr;

  const	h = 1d-4;

  var	i, j, n		: integer;
	M, V		: svdmat;
	w, b, x		: svdarray;

  begin
    iter := succ(iter);
    { compute the deviations }
    if first then
    begin
      writeln;
      writeln('  iteration ', 0:1);
      writeln;
    end;
    for i:=1 to nset do
    begin
      getelem(elem_getpos(Floc[i], Kloc[i]), cell);
      b[i] := -(cell.alpha[1]-alphax[i]);
      b[i+nset] := -(cell.beta[1]-betax[i]);
      b[i+2*nset] := -(cell.eta[1]-etax[i]);
      b[i+3*nset] := -(cell.etap[1]-etaxp[i]);
      b[i+4*nset] := -(cell.alpha[2]-alphay[i]);
      b[i+5*nset] := -(cell.beta[2]-betay[i]);

      if first then
      begin
        writeln('  alphax: ', cell.alpha[1]:6:3, ' (', alphax[i]:6:3, ')',
		', betax:  ', cell.beta[1]:6:3,  ' (', betax[i]:6:3, ')');
        writeln('  etax:   ', cell.eta[1]:6:3,   ' (', etax[i]:6:3, ')',
		', etaxp:  ', cell.etap[1]:6:3,  ' (', etaxp[i]:6:3, ')');
        writeln('  alphay: ', cell.alpha[2]:6:3, ' (', alphay[i]:6:3, ')',
	        ', betay:  ', cell.beta[2]:6:3,  ' (', betay[i]:6:3, ')');
      end;
    end;

    { compute the Jacobian }
    if nknob > nset*nparam then
      n := nknob
    else
      n := nset*nparam;
    for i:=1 to n do
      for j:=1 to n do
	M[i, j] := 0d0;

    Ring_gettwiss(true, 0);

    for i:=1 to nknob do
    begin
      for j:=1 to nset do
      begin
        getelem(elem_getpos(Floc[j], Kloc[j]), cell);
	M[j, i] := cell.alpha[1];
	M[nset+j, i] := cell.beta[1];
	M[2*nset+j, i] := cell.eta[1];
	M[3*nset+j, i] := cell.etap[1];
	M[4*nset+j, i] := cell.alpha[2];
	M[5*nset+j, i] := cell.beta[2];
      end;
    end;

    for i:=1 to nknob do
    begin
      setkpar(Fknob[i], Kknob[i], quad, knobs[i]+h);
      { set the same strenght for all kids }
      for j:=2 to getnkid(Fknob[i]) do
        setkpar(Fknob[i], j, quad, knobs[i]+h);
      Ring_gettwiss(true, 0);
      for j:=1 to nset do
      begin
        getelem(elem_getpos(Floc[j], Kloc[j]), cell);
	M[j, i] := (cell.alpha[1] - M[j, i])/h;
	M[nset+j, i] := (cell.beta[1] - M[nset+j, i])/h;
	M[2*nset+j, i] := (cell.eta[1] - M[2*nset+j, i])/h;
	M[3*nset+j, i] := (cell.etap[1] - M[3*nset+j, i])/h;
	M[4*nset+j, i] := (cell.alpha[2] - M[4*nset+j, i])/h;
	M[5*nset+j, i] := (cell.beta[2] - M[5*nset+j, i])/h;
      end;
      setkpar(Fknob[i], Kknob[i], quad, knobs[i]);
      { set the same strenght for all kids }
      for j:=2 to getnkid(Fknob[i]) do
        setkpar(Fknob[i], j, quad, knobs[i]);
    end;

    writeln;
    for i:=1 to nset*nparam do
    begin
      for j:=1 to nknob do
      begin
	write(' ', M[i, j]:9);
	if j mod 8 = 0 then writeln;
      end;
      if j mod 8 <> 0 then writeln;
      writeln;
    end;

    svdcmp(M, nset*nparam, nknob, w, V);

    writeln;
    writeln('  singular value decomposition');
    writeln;
    for i:=1 to nknob do
    begin
      write('  w[', i:2, ']: ', w[i]:9);
      if i mod 4 = 0 then writeln;
    end;
    if i mod 4 <> 0 then writeln;

    if first then
    begin
      writeln;
      write('  cut from> '); readln(cut); writeln;
      first := false;
    end;

    for i:=cut to nset*nparam do
      w[i] := 0d0;

    svbksb(M, w, V, nset*nparam, nknob, b, x);

    writeln;
    writeln('  increments');
    writeln;
    for i:=1 to nknob do
    begin
      knobs[i] := knobs[i] + x[i];
      write(' ', x[i]:9);
      if i mod 4 = 0 then writeln;
    end;
    if i mod 4 <> 0 then writeln;
    writeln;

    for i:=1 to nknob do
      setkpar(Fknob[i], Kknob[i], quad, knobs[i]);
    { set the same strenght for all kids }
    for i:=1 to nknob do
      for j:=2 to getnkid(Fknob[i]) do
        setkpar(Fknob[i], j, quad, knobs[i]);

    Ring_gettwiss(true, 0);

    eps := 0d0;
    if iter mod nprint = 0 then
    begin
      writeln;
      writeln('  iteration ', iter:1);
      writeln;
    end;
    for i:=1 to nset do
    begin
      getelem(elem_getpos(Floc[i], Kloc[i]), cell);
      b[i] := -(cell.alpha[1]-alphax[i]);
      b[i+nset] := -(cell.beta[1]-betax[i]);
      b[i+2*nset] := -(cell.eta[1]-etax[i]);
      b[i+3*nset] := -(cell.etap[1]-etaxp[i]);
      b[i+4*nset] := -(cell.alpha[2]-alphay[i]);
      b[i+5*nset] := -(cell.beta[2]-betay[i]);

      eps := eps + abs(b[i]) + abs(b[i+nset]) + abs(b[i+2*nset])
	  + abs(b[i+3*nset]);

      if iter mod nprint = 0 then
      begin
        writeln('  alphax: ', cell.alpha[1]:6:3, ' (', alphax[i]:6:3, ')',
		', betax:  ', cell.beta[1]:6:3,  ' (', betax[i]:6:3, ')');
        writeln('  etax:   ', cell.eta[1]:6:3,   ' (', etax[i]:6:3, ')',
		', etaxp:  ', cell.etap[1]:6:3,  ' (', etaxp[i]:6:3, ')');
        writeln('  alphay: ', cell.alpha[2]:6:3, ' (', alphay[i]:6:3, ')',
	        ', betay:  ', cell.beta[2]:6:3,  ' (', betay[i]:6:3, ')');
      end;
    end;
  end;


  begin
    Ring_gettwiss(true, 0);

    first := true;
    repeat
      nr;
    until (iter = niter) or (eps < ftol) ;

  end;
