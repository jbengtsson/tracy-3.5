  procedure sext_anal(scl : double; mx, my, md, nx, ny : integer;
		      var n : integer; var x, y : graphvect; var r : double);
  var	i		: integer;
	phi, sx, sy 	: double;
	cell		: celltype;
  begin
    Ring_gettwiss(false, 0.0);
    n := 0; sx := 0.0; sy := 0.0;
    for i:=1 to globval.cell_nloc do
    begin
      getelem(i, cell);
      if (cell.elem.pkind = mpole) and (cell.elem.m.porder = sext) then
      begin
	n := succ(n); x[n] := 0.0; y[n] := 0.0;
	r := scl*cell.elem.m.pbpar[sext]
	     *pwr(2.0*cell.beta[1], mx/2.0)*pwr(2.0*cell.beta[2], my/2.0);
	if (md <> 0) and (cell.eta[1] <> 0.0) then r := r*pwr(cell.eta[1], md);
	phi := 2.0*pi*(nx*cell.nu[1]+ny*cell.nu[2]);
	n := succ(n); x[n] := r*cos(phi); y[n] := r*sin(phi);
	sx := sx + x[n]; sy := sy + y[n];
      end;
    end;
    r := sqrt(sqr(sx)+sqr(sy));
  end;
