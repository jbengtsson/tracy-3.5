  var	mfile1	: text;

  procedure prtmfile;

  { Routine to generate machine file.

    Type codes:			marker					-1
				drift					 0
				multipole				 1
				cavity					 2
				thin kick				 3

    Integration methods:	not applicable				 0
				2:nd order symplectic integrator	 2
				4:th order symplectic integrator	 4


    Format:	name, family no, kid no, element no
		type code, integration method, no of integration steps


    The following lines follows depending on element type.


    type

    drift:	L

    multipole:	hor., ver. displacement, roll angle (design), roll angle (error)
  		L, 1/rho, entrance angle, exit angle
		no of nonzero multipole coeff.
		n, b , a
		    n   n
		    .
		    .
		    .

    cavity:	cavity voltage/beam energy [eV], omega/c, beam energy [eV]

    thin kick:	hor., ver. displacement, roll angle (total)
		no of nonzero multipole coeff.
		n, b , a
		    n   n
		    .
		    .
		    .
  }

  const	snamelen = 10;

	c0 = 2.99792458d8;

	{ type codes }
	marker_ = -1; drift_ = 0; mpole_ = 1; cavity_ = 2; thinkick_ = 3;


  var	i, h						: integer;
	Pmethod, PN, nmpole				: integer;
	pi, fkick1, fdrift1, fkick2, fdrift2, frf, Vrf	: double;
	Cell						: Celltype;


  procedure prtName(Name : Partsname; Fnum, Knum : integer);

  var 	j	: integer;

  begin
    for j:=1 to snamelen do
      write(mfile1, Name[j]);
    write(mfile1, '       ', Fnum, Knum);
  end;


  procedure prtDrift(L : double);

  begin
    writeln(mfile1, drift_:3, 0:3, 0:3);
    writeln(mfile1, ' ', L:23);
  end;


  function getnmpole(Fnum, Knum, Order : integer; L : double) : integer;

  var	i, k	: integer;

  begin
    k := 0;
    for i:=1 to Order do
      if (cell.elem.M.PB[-i] <> 0d0) or (cell.elem.M.PB[i] <> 0d0) then
	k := succ(k);
    getnmpole := k;
  end;


  procedure prtThick(L, irho, Tx1, Tx2 : double;
		       Pmethod, PN, nmpole : integer;
		       dS : vector2; dTpar, dTsys, dTrms, dTrnd : double);

  begin
    writeln(mfile1, mpole_:3, Pmethod:3, '  ', PN:1);
    writeln(mfile1, ' ', dS[1]:23, ' ', dS[2]:23,
	    ' ', dTpar:23, ' ', dTsys+dTrms*dTrnd:23);
    writeln(mfile1, ' ', L:23, ' ', irho:23, ' ', Tx1:23, ' ', Tx2:23);
    writeln(mfile1, '  ', nmpole:1);
  end;


  procedure prtThin(nmpole : integer; dS : vector2;
		      dTpar, dTsys, dTrms, dTrnd : double);

  begin
    writeln(mfile1, thinkick_:3, 0:3, 0:3);
    writeln(mfile1, ' ', dS[1]:23, ' ', dS[2]:23,
	    ' ', dTpar+dTsys+dTrms*dTrnd:23);
    writeln(mfile1, '  ', nmpole:1);
  end;


  procedure prtHOM(Order : integer);

  var	i	: integer;

  begin
    for i:=1 to Order do
    begin
      if (cell.elem.M.PB[-i] <> 0d0) or (cell.elem.M.PB[i] <> 0d0) then
	writeln(mfile1, i:3, ' ', cell.elem.M.PB[i]:23, ' ',
		cell.elem.M.PB[-i]:23);
    end;
  end;


  procedure prtCavity(frf, Vrf : double; h : integer);

  begin
    writeln(mfile1, cavity_:3, 0:3, 0:3);
    writeln(mfile1, ' ', Vrf/(globval.Energy*1d9):23,
            ' ', 2*pi*frf/c0:23, ' ', h:1, ' ', (globval.Energy*1d9):23);
  end;


  begin
    pi := 4*arctan(1.0);

    getglobv_(globval);

    for i:=1 to globval.Cell_nLoc do
    begin
      getelem(i, Cell);
      case Cell.Elem.Pkind of
	drift:	begin
		  prtName(Cell.Elem.Pname, cell.Fnum, cell.Knum);
		  writeln(mfile1,  '   ', i:4, ' ');
		  prtDrift(Cell.Elem.PL);
		end;

        mpole:	begin
		  prtName(Cell.Elem.Pname, cell.Fnum, cell.Knum);
		  writeln(mfile1, '   ', i:4, ' ');
		  Mpole_GetPmeth(cell.Fnum, cell.Knum, Pmethod, PN,
				 fdrift1, fkick1, fdrift2, fkick2);
		  if cell.elem.PL <> 0d0 then
		  begin
		    prtThick(cell.elem.PL, cell.elem.M.Pirho,
			       cell.elem.M.PTx1, cell.elem.M.PTx2,
			       cell.elem.M.Pmethod, cell.elem.M.PN,
			       getnmpole(cell.Fnum, cell.Knum,
					 cell.elem.M.Porder, cell.elem.PL),
			       Cell.dS,
			       Cell.Elem.M.PdTpar, Cell.Elem.M.PdTsys,
			       Cell.Elem.M.PdTrms, Cell.Elem.M.PdTrnd);
		    prtHOM(cell.elem.M.Porder);
		  end
		  else
		  begin
		    prtThin(getnmpole(cell.Fnum, cell.Knum,
					cell.elem.M.Porder, cell.elem.PL),
			      Cell.dS,
			      Cell.Elem.M.PdTpar, Cell.Elem.M.PdTsys,
			      Cell.Elem.M.PdTrms, Cell.Elem.M.PdTrnd);
		    prtHOM(cell.elem.M.Porder);
		  end;
                end;

        cavity:	begin
		  prtName(Cell.Elem.Pname, cell.Fnum, cell.Knum);
		  writeln(mfile1, '   ', i:4, ' ');
		  Cav_Get(cell.Fnum, cell.Knum, frf, Vrf, h);
		  prtCavity(frf, Vrf, h);
		end;

        marker:	begin
		  prtName(Cell.Elem.Pname, cell.Fnum, cell.Knum);
		  writeln(mfile1,  '   ', i:4, ' ');
		  writeln(mfile1, marker_:3, 0:3, 0:3);
		end;

	undef:	writeln(mfile1, ' ?? ');
      end;
    end;
  end;
