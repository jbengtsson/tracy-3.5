  module t2cell(input, output);

  %include 'mathlib.def'
  %include 'dab.def'

  %include 'pascommon.def'
  %include 't2common.def'
  %include 't2elem.def'
  %include 't2cello.def'

  { interface }

  var	ntransfmat	: [global] integer;
	transfmat	: [global] array [1..maxtransfmat] of matrix;
	kicks		: [global] array [1..maxtransfmat, 1..maxkicks] of integer;

  [global] procedure Cell_SetdP(dP:double); forward;

  [global] procedure Cell_Init; forward;

  [global] procedure Elem_Pass(i : integer; var x : vector); forward;

  [global] procedure Elem_Pass_M(i : integer; var xref : vector;
				    var x : matrix); forward;

  [global] procedure Elem_DApass(i : integer; var map : DAmap); forward;

  [global] procedure Cell_Pass(i0, i1 : integer; var x : vector;
		      var lastpos : integer); forward;
     { Track particle from i0 to i1, set x[5]:=dP }

  [global] procedure Cell_Pass_M(i0, i1 : integer; var xref : vector;
		      var mat : matrix; var lastpos : integer); forward;
     { Track matrix from i0 to i1 around ref. orbit }

  [global] procedure Cell_DApass(i0, i1 : integer; var map : DAmap;
				 var lastpos : integer); forward;
     { Track matrix from i0 to i1 around ref. orbit }

  [global] procedure Cell_Concat(dP : double); forward;

  [global] procedure Cell_fPass(var x : vector; var lastpos : integer); forward;

  [global] procedure Cell_fPass_M(var xref : vector; var mat : matrix;
				  var lastpos : integer); forward;

  [global] procedure Cell_GetCOD(imax : integer; eps, dP : double;
				 var lastpos : integer); forward;

  { implementation }

  procedure Cell_SetdP{dP:double};
  var i, j	: integer;
  begin
    globval.dPparticle := dP;
    for i:=1 to globval.Elem_NFam do
      with ElemFam[i], ElemF do
      begin
        case Pkind of 
          drift  : begin
		     for j:=1 to nKid do
		       Drift_SetMatrix(i, j);
		   end;
          mpole  : begin
		     for j:=1 to nKid do
		       Mpole_SetPB(i, j, 2);
		   end;
          wigl   : begin
		     for j:=1 to nKid do
		       Wiggler_SetPB(i, j, 2);
		   end;
	  cavity : { nothing };
	  marker : { nothing };
	  otherwise writeln('** undefined type');
        end;
      end;
  end;

  procedure Cell_Init;
  var i		: integer;
      Stotal	: double;
  begin
    if debug then writeln('**  Cell_Init');
    InitFkick;
    for i:=1 to globval.Elem_NFam do
      with ElemFam[i], ElemF do
      begin
        if debug then writeln('**  Cell_Init,   i:=',i:3,':   ',Pname);
        case Pkind of 
          drift : Drift_Init(i);
          mpole : Mpole_Init(i);
          wigl  : Wiggler_Init(i);
	  cavity: Cav_Init(i);
	  marker: Marker_Init(i);
	  otherwise writeln('** undefined type');
        end;
      end;
    Stotal:=0;
    for i:=0 to globval.Cell_nLoc do
    begin
      Stotal:=Stotal+Cell[i].Elem.PL;
      Cell[i].S:=Stotal;
    end;   
  end;

  function CheckAmpl(var x : vector; loc : integer;
		     var lastpos : integer) : boolean;
  begin
    CheckAmpl := true;
    if not (abs(x[1]) < globval.maxampl[1]) or
       not (abs(x[3]) < globval.maxampl[2]) then
    begin
      lastpos := loc;
      if trace then
        writeln('Particle lost at element nb:', loc:5,
   		'  s =', Cell[loc].s:10:5);
      CheckAmpl := false;
    end;
  end;

  procedure Elem_Pass{i : integer; var x : vector};

  begin
    with Cell[i], Elem do
    begin
      case Pkind of 
        drift  : Drift_Pass(Cell[i], x);
        mpole  : Mpole_Pass(Cell[i], x);
        wigl   : Wiggler_Pass(Cell[i], x);
        cavity : Cav_Pass(Cell[i], x);
        marker : Marker_Pass(Cell[i], x);
	otherwise writeln('** undefined type');
      end;
    end;
  end;

  procedure Elem_Pass_M{i : integer; var xref : vector; var x : matrix};
  begin
    with Cell[i], Elem do
    begin
      case Pkind of 
        drift  : Drift_Pass_M(Cell[i], xref, x);
        mpole  : Mpole_Pass_M(Cell[i], xref, x);
        wigl   : Wiggler_Pass_M(Cell[i], xref, x);
	cavity : { nothing };
	marker : { nothing };
	otherwise writeln('** undefined type');
      end;
    end;
  end;

  procedure Elem_DApass{i : integer; var map : DAmap};
  begin
    with Cell[i], Elem do
    begin
      case Pkind of 
        drift  : Drift_DApass(Cell[i], map);
        mpole  : Mpole_DApass(Cell[i], map);
        wigl   : Wiggler_DApass(Cell[i], map);
        cavity : cav_DApass(Cell[i], map);
	marker : { nothing };
	otherwise writeln('** undefined type');
      end;
    end;
  end;

  procedure Cell_Pass{i0, i1 : integer; var x : vector;
		      var lastpos : integer};
  label	999;
  var	i	: integer;
  begin
    if x[5] <> globval.dPparticle then Cell_SetDP(x[5]);
    lastpos := i1;
    if not CheckAmpl(x, i0, lastpos) then goto 999;
    cell[i0].beampos := x;
    for i:=i0+1 to i1 do
    begin
      Elem_Pass(i, x);
      if not CheckAmpl(x, i, lastpos) then goto 999;
    end;
999:
  end;

  procedure Cell_Pass_M{i0, i1 : integer; var xref : vector; var mat : matrix;
			var lastpos : integer};
  label	999;
  var	i	: integer;
  begin
    if xref[5] <> globval.dPparticle then Cell_SetDP(xref[5]);
    lastpos := i1;
    if not CheckAmpl(xref, i0, lastpos) then goto 999;
    for i:=i0+1 to i1 do
    begin
      Elem_Pass_M(i, xref, mat);
      if not CheckAmpl(xref, i, lastpos) then goto 999;
    end;
999:
  end;

  procedure Cell_DApass{i0, i1 : integer; var map : DAmap;
			var lastpos : integer};
  label	999;
  var	i	: integer;

  function CheckAmpl(var map : DAmap; loc : integer;
		     var lastpos : integer) : boolean;
  begin
    CheckAmpl := true;
    if not (abs(map[1, 0]) < globval.maxampl[1]) or
       not (abs(map[3, 0]) < globval.maxampl[2]) then
    begin
      lastpos := loc;
      if trace then
        writeln('Particle lost at element nb:', loc:5,
		'  s =', Cell[loc].s:10:5);
      CheckAmpl := false;
    end;
  end;

  begin
    lastpos := i1;
    if not CheckAmpl(map, i0, lastpos) then goto 999;
    if globval.radiation then globval.dE := 0.0;
    if globval.emittance then
      for i:=1 to 3 do
	globval.qfluct[i] := 0d0;
    for i:=i0+1 to i1 do
    begin
      Elem_DApass(i, map);
      if not CheckAmpl(map, i, lastpos) then goto 999;
    end;
999:
  end;

  procedure Cell_Concat{dP : double};
  const	n=4;
  var	j, i1	: integer;
	PB2	: double;

  function linearel(i : integer) : boolean;
  begin
    with Cell[i], Elem do
      if Pkind = Drift then
	linearel := true
      else if Pkind = Mpole then
      begin
	if (M^.Porder <= quad) and (M^.PB[-quad] = 0d0) then
	  linearel := true
	else
	  linearel := false;
      end
      else if Pkind = Wigl then
	linearel := true
      else if Pkind = Cavity then
	linearel := true
      else if Pkind = marker then
	linearel := true;
  end;

  procedure GtoL_dP(var mat : matrix; var dT : vector2);
  var	k	: integer;
	dS0	: vector2;
	x	: vector;
  begin
    dS0[1] := 0.0; dS0[2] := 0.0;
    for k:=1 to n do
      x[k] := mat[k, n+1];
    GtoL(x, dS0, dT, 0.0, 0.0, 0.0);
    for k:=1 to n do
      mat[k, n+1] := x[k];
  end;

  procedure LtoG_dP(var mat : matrix; var dT : vector2);
  var	k	: integer;
	dS0	: vector2;
	x	: vector;
  begin
    dS0[1] := 0.0; dS0[2] := 0.0;
    for k:=1 to n do
      x[k] := mat[k, n+1];
    LtoG(x, dS0, dT, 0.0, 0.0, 0.0);
    for k:=1 to n do
      mat[k, n+1] := x[k];
  end;

  begin
    if dP <> globval.dPparticle then
    begin
      Cell_SetDP(dP); cellconcat := false;
    end;
    if not cellconcat then
    begin
      if trace then writeln('concatenating');
      cellconcat := true;
      i1 := 0; ntransfmat := 1;
      UnitMat(n+1, transfmat[ntransfmat]); kicks[ntransfmat, 1] := 0;
      repeat
        while linearel(i1+1) and (i1+1 <= globval.Cell_nLoc) do
        begin
	  i1 := succ(i1);
	  with Cell[i1], Elem do
	  begin
	    case Pkind of
	      drift : MulLsMat(D^.D55, transfmat[ntransfmat]);
	      mpole : with M^ do 
		      begin
			GtoL_M(transfmat[ntransfmat], dT);
			GtoL_dP(transfmat[ntransfmat], dT);
			GtoL(transfmat[ntransfmat, n+1], dS, dT, Pc0, Pc1, Ps1);
	                if Pthick = thick then
		        begin
		          MulLsMat(AU55, transfmat[ntransfmat]);
			  { Assuming there is no quadrupole kick }
			  PB2 := PB[quad]; PB[quad] := 0.0;
			  thinkick(Porder, PB, PL, 0, pthick, transfmat[ntransfmat, n+1]);
			  PB[quad] := PB2;
		          MulLsMat(AD55, transfmat[ntransfmat]);
		        end
	                else
		          { Dipole kick }
		          thinkick(Porder, PB, 1, 0, pthick, transfmat[ntransfmat, n+1]);
			LtoG_M(transfmat[ntransfmat], dT);
			LtoG_dP(transfmat[ntransfmat], dT);
			LtoG(transfmat[ntransfmat, n+1], dS, dT, Pc0, Pc1, Ps1);
		      end;
	      wigl  : MulLsMat(W^.W55, transfmat[ntransfmat]);
	      cavity: { nothing };
	      marker: { nothing };
	      otherwise writeln('** undefined type');
	    end;
	  end;
        end;
        j := 0;
        while not linearel(i1+j+1) and ((i1+j+1) <= globval.Cell_nLoc) do
        begin
	  j := succ(j);
	  with Cell[i1+j], Elem do
	  begin
	    if Pkind = mpole then
	    begin
	      with M^ do
	      begin
	        GtoL_M(transfmat[ntransfmat], dT);
		GtoL_dP(transfmat[ntransfmat], dT);
	        GtoL(transfmat[ntransfmat, n+1], dS, dT, Pc0, Pc1, Ps1);
	        if Pthick = thick then MulLsMat(AU55, transfmat[ntransfmat]);
	        kicks[ntransfmat, j] := i1+j; kicks[ntransfmat, j+1] := 0;
	        ntransfmat := succ(ntransfmat);
	        UnitMat(n+1, transfmat[ntransfmat]); kicks[ntransfmat, 1] := 0;
	        if Pthick = thick then MulLsMat(AD55, transfmat[ntransfmat]);
	        LtoG_M(transfmat[ntransfmat], dT);
		LtoG_dP(transfmat[ntransfmat], dT);
	        LtoG(transfmat[ntransfmat, n+1], dS, dT, Pc0, Pc1, Ps1);
	      end;
	    end;
	  end;
        end;
        i1 := i1 + j;
      until i1 = globval.Cell_nLoc ;
    end;
  end;

  procedure Cell_fPass{var x : vector; var lastpos : integer};
  label	999;
  const	n=4;
  var	i, j	: integer;
	PB2	: double;
  BEGIN
    lastpos := globval.Cell_nLoc;
    if not CheckAmpl(x, 1, lastpos) then goto 999;
    for i:=1 to ntransfmat do
    begin
      LinsTrans(transfmat[i], x);
      j := 0;
      while kicks[i, j+1] <> 0 do
      begin
	j := succ(j);
	with Cell[kicks[i, j]], Elem, M^ do
	begin
	  CopyVec(n, x, Cell[kicks[i, j]].BeamPos);
          if Pthick = thick then
	  begin
	    PB2 := PB[quad]; PB[quad] := 0.0;
	    thinkick(Porder, PB, PL, 0, pthick, x);
            PB[quad] := PB2;
	  end
	  else
	    thinkick(Porder, PB, 1, 0, pthick, x);
	end;
	if not CheckAmpl(x, kicks[i, j], lastpos) then goto 999;
      end;
    end;
999:
  end;

  procedure Cell_fPass_M{var xref : vector; var mat : matrix;
			 var lastpos : integer};
  label	999;
  const	n=4;
  var	i, j	: integer;
	PB2	: double;
  BEGIN
    lastpos := globval.Cell_nLoc;
    if not CheckAmpl(xref, 1, lastpos) then goto 999;
    for i:=1 to ntransfmat do
    begin
      MulLsMat(transfmat[i], mat); LinsTrans(transfmat[i], xref);
      j := 0;
      while kicks[i, j+1] <> 0 do
      begin
	j := succ(j);
	with Cell[kicks[i, j]], Elem, M^ do
	begin
          if Pthick = thick then
	  begin
	    PB2 := PB[quad]; PB[quad] := 0.0;
	    thinkick_M(Porder, PB, PL, 0, pthick, xref, mat);
	    thinkick(Porder, PB, PL, 0, pthick, xref);
            PB[quad] := PB2;
	  end
	  else
	  begin
	    thinkick_M(Porder, PB, 1, 0, pthick, xref, mat);
	    thinkick(Porder, PB, 1, 0, pthick, xref);
	  end;
	end;
	if not CheckAmpl(xref, kicks[i, j], lastpos) then goto 999;
      end;
    end;
999:
  end;

  [hidden] function xabs(n : integer; var x : vector) : double;
  var	i	: integer;
	sum	: double;
  begin
    sum := 0.0;
    for i:=1 to n do
      sum := sum + sqr(x[i]);
    xabs := sqrt(sum);
  end;

  procedure Cell_MatGetCOD(imax : integer; eps, dP : double;
			   var lastpos : integer);
  { Calculate closed orbit }
  const	n=4;
  var	i, j		: integer;
	dxabs		: double;
	x0, x1, dx	: vector;
	A		: matrix;
  begin
    if globval.MatMeth then Cell_Concat(dP);
    CopyVec(n, globval.CODvect, x0); x0[5] := dP; x0[6] := 0.0;
    i := 0;
    repeat
      i := succ(i);
      UnitMat(n+2, globval.OneturnMat);
      CopyVec(n+2, x0, x1);
      Cell_fPass_M(x1, globval.OneturnMat, lastpos);
{      Cell_Pass_M(0, globval.Cell_nLoc, x1, globval.OneturnMat, lastpos);}
      if lastpos = globval.Cell_nLoc then CopyVec(n+2, x0, globval.CODvect);
      CopyVec(n, x0, dx); SubVec(n, x1, dx);
      CopyMat(n, globval.OneturnMat, A);
      for j:=1 to n do
	A[j, j] := A[j, j] - 1;
      if InvMat(n, A) then
      begin
	if lastpos = globval.Cell_nLoc then
        begin
          LinTrans(n, A, dx);
          for j:=1 to n do
            x0[j] := x0[j] + dx[j];
	end;
      end
      else
        writeln('  *** A is singular');
      dxabs := xabs(4, dx);
      if trace then
      begin
        writeln('--- CODLOOP', i:3 , ', Err=', dxabs:10, '/', eps:10); 
        writeln(x0[1], ' ', x0[2]); writeln(x0[3], ' ', x0[4]);
	writeln(x0[5], ' ', x0[6]);
      end;
    until (i >= imax) or ( dxabs <= eps ) or (lastpos <> globval.Cell_nLoc);
    CopyVec(n+2, globval.CODvect, x0);
    Cell_Pass(0, globval.Cell_nLoc, x0, lastpos);
    if (dxabs <= eps) and (lastpos = globval.Cell_nLoc) then
      status.CODflag := true
    else
    begin
      writeln(' GetCOD failed ... '); status.CODflag := false;
    end;
  end;

  procedure Cell_DAgetCOD(imax : integer; eps, dP : double;
			  var lastpos : integer);
  { Calculate closed orbit }
  var	i, j, n, k			: integer;
	dxabs				: double;
	jj				: ivector;
	x0, x1, dx			: vector;
	dd, DAdx, map, map1, map2	: DAmap;
  begin
    damapall(dd, 6, 'dd        ', 6, 6);  damapall(DAdx, 6, 'DAdx      ', 6, 6);
    damapall(map, 6, 'map       ', 6, 6); damapall(map1, 6, 'map1      ', 6, 6);
    damapall(map2, 6, 'map2      ', 6, 6);
    if globval.Cavity_on then
      n := 6
    else
      n := 4;
    globval.dPparticle := dP;
    CopyVec(4, globval.CODvect, x0); x0[5] := dP; x0[6] := 0.0;
    i := 0;
    for j:=1 to ndim2 do
    begin
      k := j; jj[k] := 0; DAvar(dd[k], 0.0, k);
    end;
    repeat
      i := succ(i);
      for j:=1 to ndim2 do
      begin
        k := j; DAvar(map[k], x0[k], k);
      end;
      Cell_DApass(0, globval.Cell_nLoc, map, lastpos);
      if lastpos = globval.Cell_nLoc then CopyVec(6, x0, globval.CODvect);
      for j:=1 to ndim2 do
      begin
        DAsub(map[j], dd[j], map1[j]); DApek(map[j], jj, x1[j]);
        DApok(map1[j], jj, 0.0); DAcon(DAdx[j], x0[j]-x1[j]);
      end;
      DAinv(map1, n, map2, n); DAcct(map2, n, DAdx, n, map1, n);
      if lastpos = globval.Cell_nLoc then
      begin
        for j:=1 to n do
        begin
          DApek(map1[j], jj, dx[j]); x0[j] := x0[j] + dx[j];
        end;
      end;
      dxabs := xabs(n, dx);
      if trace then
      begin
        writeln('--- CODLOOP', i:3 , ', Err=', dxabs:10, '/', eps:10); 
        writeln(x0[1], ' ', x0[2]); writeln(x0[3], ' ', x0[4]);
	writeln(x0[5], ' ', x0[6]);
      end;
    until (i >= imax) or ( dxabs <= eps ) or (lastpos <> globval.Cell_nLoc);
    CopyVec(6, globval.CODvect, x0);
    getlinmat(6, map, globval.Oneturnmat);
    Cell_Pass(0, globval.Cell_nLoc, x0, lastpos);
    if (dxabs <= eps) and (lastpos = globval.Cell_nLoc) then
      status.CODflag := true
    else
    begin
      writeln(' GetCOD failed ... '); status.CODflag := false;
    end;
    damapdal(dd, 6); damapdal(DAdx, 6); damapdal(map, 6); damapdal(map1, 6);
    damapdal(map2, 6);
  end;

  procedure Cell_GetCOD{imax : integer; eps, dP : double;
			var lastpos : integer};
  begin
    if globval.MatMeth then
      Cell_MatGetCOD(imax, eps, dP, lastpos)
    else
      Cell_DAGetCOD(imax, eps, dP, lastpos);
  end;

  end.
