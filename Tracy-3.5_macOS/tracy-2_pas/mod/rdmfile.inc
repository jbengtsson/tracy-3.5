  var   mfile2  : text;

  procedure rdmfile;

  { Routine to read a machine file.

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

  const snamelen = 10;

  var   i, j, nmpole, n		: integer;
        frf, Vrf, E0		: double;
        harnum, s, dTerror	: double;
        Cell			: Celltype;


  procedure rdname(var Name : Partsname);

  var   j       : integer;

  begin
    for j:=1 to snamelen do
      read(mfile2, Name[j]);
  end;


  procedure rdkind(var elem : elemtype);

  var	kind	: integer;

  begin
    read(mfile2, kind);
    case Kind of
      -1: elem.Pkind := marker;
       0: elem.Pkind := drift;
       1: begin
	    elem.Pkind := mpole; elem.M.Pthick := thick;
	  end;
       2: elem.Pkind := cavity;
       3: begin
	    elem.Pkind := mpole; elem.M.Pthick := thin;
	  end;
    end;
    if (kind <> 1) and (kind <> 3) then readln(mfile2);
  end;


  Procedure clrMpole(var M : mpoletype);

  begin
    M.Pmethod := 0; M.PN := 0;
    M.PdSsys[1] := 0d0; M.PdSsys[2] := 0d0;
    M.PdSrms[1] := 0d0; M.PdSrms[2] := 0d0;
    M.PdSrnd[1] := 0d0; M.PdSrnd[2] := 0d0;
    M.Pdtpar := 0d0; M.Pdtsys := 0d0; M.Pdtrms := 0d0; M.Pdtrnd := 0d0;
    for j := -HOMmax to HOMmax do
    begin
      M.Pbpar[j] := 0d0; M.Pbsys[j] := 0d0;
      M.Pbrms[j] := 0d0; M.Pbrnd[j] := 0d0;
      M.Pb[j] := 0d0;
    end;
    M.Porder := 0;
    M.Pthick := thin;
    M.PTx1 := 0d0; M.PTx2 := 0d0; M.Pgap := 0.0; M.Pirho := 0d0;
    M.Pc0 := 0d0; M.Pc1 := 0d0; M.Ps1 := 0d0;
  end;


  procedure clrelem(var elem : elemtype);

  begin
    Elem.PL := 0d0;
    Elem.Pkind := undef;
  end;


  procedure clrcell(var cell : celltype);

  begin
    Cell.S := 0d0;
    Cell.ds[1] := 0d0; Cell.ds[2] := 0d0;
    Cell.dt[1] := 0d0; Cell.dt[2] := 0d0;
  end;


  procedure clrglobval(var globval : globvalrec);

  var	i	: integer;

  begin
    globval.dPcommon := 0d0;
    globval.dPparticle := 0d0;
    globval.maxampl[1] := 0d0; globval.maxampl[2] := 0d0;
    globval.TotalTune[1] := 0d0; globval.TotalTune[2] := 0d0;
    globval.Omega := 0d0;
    globval.Chrom[1] := 0d0; globval.Chrom[2] := 0d0;
    globval.Energy := 0d0;
    globval.Cell_nLoc := 0;
    globval.elem_nFam := 0;
    globval.CODimax := 0;
    globval.CODeps := 0d0;
    globval.bpm := 0;
  end;


  function imax(x1, x2 : integer) : integer;

  begin
    if x1 > x2 then
      imax := x1
    else
      imax := x2;
  end;


  begin
{    clrglobval(globval);}
    s := 0d0;
    while not eof(mfile2) do
    begin
      clrcell(Cell); clrelem(Cell.Elem); clrMpole(cell.elem.M);
      rdname(Cell.Elem.Pname); readln(mfile2, Cell.Fnum, Cell.Knum, i);
      rdkind(cell.elem);
      case Cell.Elem.Pkind of
        marker: begin
                  Cell.dT[1] := 1.0; Cell.dT[2] := 0.0;
                end;
        drift:  begin
                  readln(mfile2, Cell.Elem.PL);
                end;
        cavity: begin
		  getelem(i, cell);
                  readln(mfile2, Vrf, frf, harnum, E0);
                end;
        mpole:  begin
		  readln(mfile2, cell.elem.M.Pmethod, cell.elem.M.PN);
                  if cell.elem.M.Pthick = thick then
                  begin
		    { Thick Multipole }
                    readln(mfile2, Cell.dS[1], Cell.dS[2], Cell.Elem.M.Pdtpar,
                           dTerror);
                    Cell.dT[1] := cos(dtor(dTerror+Cell.Elem.M.Pdtpar));
                    Cell.dT[2] := sin(dtor(dTerror+Cell.Elem.M.Pdtpar));
                    readln(mfile2, Cell.Elem.PL, Cell.Elem.m.Pirho,
                           Cell.Elem.M.PTx1, Cell.Elem.M.PTx2);
		    if cell.elem.M.Pirho <> 0.0 then cell.elem.M.Porder := 1;
                  end
                  else
                  begin
		    { Thin Multipole }
		    cell.elem.M.Pmethod := 2;
                    readln(mfile2, Cell.dS[1], Cell.dS[2], dTerror);  
                    Cell.dT[1] := cos(dtor(dTerror));
		    Cell.dT[2] := sin(dtor(dTerror));
                  end;
		  cell.elem.M.Pc0 := sin(cell.elem.PL*cell.elem.M.Pirho/2d0);
		  cell.elem.M.Pc1 := cos(dtor(cell.elem.M.PdTpar))
				     *cell.elem.M.Pc0;
		  cell.elem.M.Ps1 := sin(dtor(cell.elem.M.PdTpar))
				     *cell.elem.M.Pc0;
                  readln(mfile2, nmpole); { number of non zero multipoles }
                  for j:=1 to nmpole do
                  begin
                    read(mfile2, n);
		    readln(mfile2, Cell.Elem.M.PB[n], Cell.Elem.M.PB[-n]);
		    cell.elem.M.Porder := imax(n, cell.elem.M.Porder);
                  end;
                end;
      end;
      s := s + Cell.Elem.PL; Cell.S := s;
      putelem(i, cell);
    end;

    globval.cell_nloc := i; putglobv_(globval);
  end;
