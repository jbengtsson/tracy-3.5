  module t2ring(input, output);

  %include 'mathlib.def'
  %include 'dab.def'
  %include 'eigenv.def'

  %include 'pascommon.def'
  %include 't2common.def'
  %include 't2elem.def'

  %include 't2cell.def'
  %include 't2ringo.def'

  { interface }

  var	stable	: [global] boolean;

  [global] PROCEDURE Cell_GetABGN(var M : matrix;
				  var alpha, beta, gamma, nu : vector2);
				  forward;

  [global] procedure Cell_MatTwiss(i0, i1 : integer; var Ascr : matrix;
				   chroma, ring : boolean; dP : double); forward;

  [global] procedure Cell_DATwiss(i0, i1 : integer; var Ascr : DAmap;
				  chroma, ring : boolean; dP : double); forward;

  [global] procedure Ring_Getchrom(dP : double); forward;

  [global] procedure Ring_GetTwiss(chroma : boolean; dP : double); forward;

  [global] PROCEDURE Ring_Fittune(var nu : vector2; eps : double;
				  var nq : ivector2; var qf, qd : fitvect;
				  dkL : double; imax : integer); forward;

  [global] PROCEDURE Ring_Fitchrom(var ksi : vector2; eps : double;
				   var ns : ivector2; var sf, sd : fitvect;
				   dkpL : double; imax : integer); forward;

  [global] PROCEDURE Ring_FitDisp(pos : integer; eta, eps : double;
				  nq : integer; var q : fitvect;
				  dkL : double; imax : integer); forward;

  { implementation }

  PROCEDURE Cell_GetABGN{var M : matrix;
			 var alpha, beta, gamma, nu : vector2};

    { Get Twiss parameters from the transfer matrix M

			    [Nx   ]
			M = [     ]
			    [   Ny]

     where, in the case of mid-plane symmetry

	    [ cos(mu) + alpha sin(mu)        beta sin(mu)      ]
	N = [                                                  ]
	    [    -gamma sin(mu)        cos(mu) - alpha sin(mu) ]

    }

  var	i, j	: integer;
	c, s	: double;


  function arccos(x : double) : double;

  { 0 <= phi <= pi }

  begin
    if x = 0d0 then
      arccos := 1d0
    else if x >= 0.0 then
      arccos := arctan( sqrt(1/sqr(x)-1) )
    else
      arccos := pi - arctan( sqrt(1/sqr(x)-1) );
  end;


  procedure getnu(var nu : vector2);

  { Get nu from any symplectic matrix }

  const	n = 4;

  var	i, j					: integer;
	sgn, detp, detm, b, c, th, tv, b2mc	: double;
	M1					: matrix;

  begin
    CopyMat(n, M, M1);

    for i:=1 to n do
      M1[i, i] := M1[i, i] - 1.0;
    detp := DetMat(n, M1);

    for i:=1 to n do
      M1[i, i] := M1[i, i] + 2.0;
    detm := DetMat(n, M1);

    for i:=1 to n do
      M1[i, i] := M1[i, i] - 1.0;

    b := (detp-detm)/16.0; c := (detp+detm)/8.0 - 1.0;

    th := (M1[1, 1]+M1[2, 2])/2.0; tv := (M1[3, 3]+M1[4, 4])/2.0;

    b2mc := sqr(b) - c;
    if b2mc < 0d0 then
    begin
      nu[1] := -1.0; nu[2] := -1.0;
      writeln('*** Getnu: unstable in tune');
    end
    else
    begin
      if th > tv then
	sgn := 1.0
      else
	sgn := -1.0;

      nu[1] := arccos(-b+sgn*sqrt(b2mc))/(2.0*Pi);
      nu[2] := arccos(-b-sgn*sqrt(b2mc))/(2.0*Pi);

      for i:=1 to 2 do
      begin
        j := 2*i - 1;
        if M1[j, j+1] < 0d0 then nu[i] := 1d0 - nu[i];
      end;
    end;
  end;


  begin
    stable := true;
    for i:=1 to 2 do
    begin
      j := 2*i - 1;
      c := (M[j, j] + M[j+1, j+1])/2;
      stable := abs(c) < 1.0;
      if stable then
      begin
        s := SQRT(1-sqr(c))*sign(M[j, j+1]);
        alpha[i] := (M[j, j] - M[j+1, j+1])/(2*s);
        beta[i] := M[j, j+1]/s;
        gamma[i] := -M[j+1, j]/s;            
	getnu(nu);
      end;
    end;
  end;


  procedure Cell_Geteta(i0, i1 : integer; ring : boolean; dP : double);

  { Calculate eta and eta' by numerical differentiation }

  const	n = 4;

  var	i, j, k, lastpos	: integer;
	xref			: vector;
	codbuf			: array [0..Cell_nLocmax] of vector;

  begin
    if ring then
      Cell_GetCOD(globval.CODimax, globval.CODeps,
		  dP-globval.dPcommon/2d0, lastpos)
    else
    begin
      CopyVec(n+2, globval.CODvect, xref);
      xref[5] := dP - globval.dPcommon/2d0;
      Cell_Pass(i0, i1, xref, lastpos);
    end;
    { Store reference orbit }
    for i:=i0 to i1 do
      CopyVec(n+2, Cell[i].beampos, codbuf[i]);
    if ring then
      Cell_GetCOD(globval.CODimax, globval.CODeps,
		  dP+globval.dPcommon/2d0, lastpos)
    else
    begin
      CopyVec(n+2, globval.CODvect, xref);
      xref[5] := dP + globval.dPcommon/2d0;
      Cell_Pass(i0, i1, xref, lastpos);
    end;
    for i:=i0 to i1 do
    begin
      with Cell[i] do
      begin
	for j:=1 to 2 do
	begin
	  k := 2*j - 1;
	  Eta[j] := (beampos[k]-codbuf[i, k])/globval.dPcommon;
	  Etap[j] := (beampos[k+1]-codbuf[i, k+1])/globval.dPcommon;
	end;
      end;
    end;
  end;


  procedure Cell_MatTwiss{i0, i1 : integer; var Ascr : matrix;
			  chroma, ring : boolean; dP : double};

  { Calculate Twiss parameters from element i0 to element i1 }

  const	n=4;

  var	i, j, k		: integer;
        nu1, dnu	: vector2;
	xref		: vector;
	Ascr0, Ascr1	: matrix;

  procedure getprm(var Ascr : matrix; var alpha, beta : vector2);

  var 	i, j	: integer;

  begin
    for i:=1 to 2 do
    begin
      j := 2*i - 1;
      alpha[i] := -( Ascr[j, j]*Ascr[j+1, j] + Ascr[j, j+1]*Ascr[j+1, j+1] );
      beta[i] := sqr(Ascr[j, j]) + sqr(Ascr[j, j+1]);
    end;
  end;


  begin
    if dP <> globval.dPparticle then Cell_SetDP(dP);
    for j:=1 to 2 do
      nu1[j] := 0.0;
    with Cell[i0] do
    begin
      getprm(Ascr, alpha, beta); Nu := nu1;
    end;
    CopyMat(n+1, Ascr, Ascr0); CopyVec(n+2, globval.CODvect, xref);
    for i:=i0+1 to i1 do
    begin
      CopyMat(n+1, Ascr0, Ascr1); Elem_Pass_M(i, xref, Ascr1);
      with Cell[i] do
      begin
        getprm(Ascr1, alpha, beta);
        for j:=1 to 2 do
        begin
          k := 2*j - 1;
	  dnu[j] := ( getangle( Ascr1[k, k], Ascr1[k, k+1] )
	         - getangle( Ascr0[k, k], Ascr0[k, k+1] ) )/(2*pi);
          if dnu[j] < -1e-16 then dnu[j] := dnu[j] + 1.0;
          nu1[j] := nu1[j] + dnu[j]; Nu[j] := nu1[j];
	  { Only correct for bare lattice }
	  Eta[j] := Ascr1[k, 5]; Etap[j] := Ascr1[k+1, 5];
        end;
      end;
      CopyMat(n+1, Ascr1, Ascr0);
    end;
    if chroma then Cell_Geteta(i0, i1, ring, dP);
  end;


  procedure Cell_DATwiss{i0, i1 : integer; Ascr : DAmap;
			 chroma, ring : boolean; dP : double};

  { Calculate Twiss parameters from element i0 to element i1 }

  const	n=4;

  var	i, j, k		: integer;
        nu1, dnu	: vector2;
	Ascr0, Ascr1	: DAmap;

  procedure DAgetprm(var Ascr : DAmap; var alpha, beta : vector2);

  var	i, j	: integer;

  begin
    for i:=1 to 2 do
    begin
      j := 2*i - 1;
      alpha[i] := -( getmat(Ascr, j, j)*getmat(Ascr, j+1, j)
	       + getmat(Ascr, j, j+1)*getmat(Ascr, j+1, j+1) );
      beta[i] := sqr(getmat(Ascr, j, j)) + sqr(getmat(Ascr, j, j+1));
    end;
  end;


  begin
    for j:=1 to 2 do
      nu1[j] := 0.0;
    with Cell[i0] do
    begin
      DAgetprm(Ascr, alpha, beta); Nu := nu1;
    end;
    CopyMap(Ascr, Ascr0);
    for j:=1 to n+2 do
      Ascr0[j, 0] := globval.CODvect[j];
    if globval.radiation then globval.dE := 0.0;
    for i:=i0+1 to i1 do
    begin
      CopyMap(Ascr0, Ascr1);
      Elem_DApass(i, Ascr1);
      with Cell[i] do
      begin
        DAgetprm(Ascr1, alpha, beta);
        for j:=1 to 2 do
        begin
	  k := 2*j - 1;
	  dnu[j] := ( getangle( getmat(Ascr1, k, k), getmat(Ascr1, k, k+1) )
	         - getangle( getmat(Ascr0, k, k), getmat(Ascr0, k, k+1) ) )/(2*pi);
          if dnu[j] < -1e-16 then dnu[j] := dnu[j] + 1.0;
          nu1[j] := nu1[j] + dnu[j];
          Nu[j] := nu1[j];
	  eta[j] := getmat(Ascr1, k, 5)*getmat(Ascr1, 6, 6)
		    -  getmat(Ascr1, k, 6)*getmat(Ascr1, 6, 5);
	  etap[j] := getmat(Ascr1, k+1, 5)*getmat(Ascr1, 6, 6)
		     -  getmat(Ascr1, k+1, 6)*getmat(Ascr1, 6, 5);
        end;
      end;
      CopyMap(Ascr1, Ascr0);
    end;
    if chroma and not globval.Cavity_on then Cell_Geteta(i0, i1, ring, dP);
  end;


  procedure Ring_Getchrom{dP : double};

  const	n=4;

  var	j, lastpos			: integer;
	alpha, beta, gamma, nu, nu0	: vector2;

  begin
    { Get tune }
    Cell_GetCOD(globval.CODimax, globval.CODeps, dP-globval.dPcommon/2d0,
	        lastpos);
    if status.codflag then
    begin
      Cell_GetABGN(globval.OneTurnMat, alpha, beta, gamma, nu0);
      Cell_GetCOD(globval.CODimax, globval.CODeps, dP+globval.dPcommon/2d0,
		  lastpos);
      if status.codflag then
      begin
        Cell_GetABGN(globval.OneTurnMat, alpha, beta, gamma, nu);
        if stable then
        begin
          { Get chromaticities }
          for j:=1 to 2 do
            globval.chrom[j] := (nu[j]-nu0[j])/globval.dPcommon;
          status.chromflag := true;
        end
        else
	  writeln('  Lattice is unstable');
      end;
    end
    else
      writeln('  Lattice is unstable');
  end;


  procedure Ring_MatTwiss(chroma : boolean; dP : double);

  { M : One turn transfer matrix

                                  [ cx  sx  .   . ]
		   		  [               ]
                   -1             [-sx  cx  .   . ]
	    M  -> A   M A = R =   [               ]
				  [ .   .   cy  sy]
				  [               ]
                                  [ .   .  -sy  cy]

  }

  const	n=4;

  var	j, lastpos		: integer;
	alpha, beta, gamma, nu	: vector2;
	eta0			: vector;
	R, C, Ascr		: matrix;

  begin
    { Get closed orbit }
    Cell_GetCOD(globval.CODimax, globval.CODeps, dP, lastpos);
    if status.codflag then
    begin
      { Check if stable }
      Cell_GetABGN(globval.OneTurnMat, alpha, beta, gamma, nu);

      { Get eigenvalues and eigenvectors for the one turn transfer matrix }
      gdiag(n, cell[globval.cell_nloc].s, globval.Ascr, globval.Ascrinv, R,
	    globval.OneTurnMat, globval.Omega, globval.alphac);

      { Only correct for bare lattice }
      for j:=1 to n+1 do
        eta0[j] := globval.oneturnmat[j, n+1];
      unitmat(n, C); submat(n, globval.oneturnmat, C);
      if not invmat(n, C) then
        writeln('** matrix is singular');
      lintrans(n, C, eta0);
      for j:=1 to n+1 do
      begin
	globval.Ascr[n+1, j] := 0.0; globval.Ascr[j, n+1] := eta0[j];
      end;

      CopyMat(n+1, globval.Ascr, Ascr);
      Cell_MatTwiss(0, globval.Cell_nLoc, Ascr, chroma, true, dP);
      globval.TotalTune := Cell[globval.Cell_nLoc].Nu;
      status.tuneflag := true;
      if chroma then
      begin
	Ring_Getchrom(dP);
        Cell_GetCOD(globval.CODimax, globval.CODeps, dP, lastpos);
      end;
    end;
  end;


  procedure Ring_DATwiss(chroma : boolean; dP : double);

  var	j, n, lastpos		: integer;
	alpha, beta, gamma, nu	: vector2;
	R			: matrix;
	DAAScr			: DAmap;

  begin
    if globval.Cavity_on then
      n := 6
    else
      n := 4;
    { Get closed orbit }
    Cell_GetCOD(globval.CODimax, globval.CODeps, dP, lastpos);
    if status.codflag then
    begin
      { Check if stable }
      Cell_GetABGN(globval.OneTurnMat, alpha, beta, gamma, nu);
      { Get eigenvalues and eigenvectors for the one turn transfer matrix }
      gdiag(n, cell[globval.cell_nloc].s, globval.Ascr, globval.Ascrinv, R,
	    globval.OneTurnMat, globval.Omega, globval.alphac);
      putlinmat(n, globval.Ascr, DAAscr);
      if not globval.CAvity_on then
      begin
        for j:=1 to 4 do
	begin
          DAAscr[j, 5] := 0.0; DAAscr[j, 6] := 0.0;
        end;
        DAcon(DAAscr[5], 0.0); DAcon(DAAscr[6], 0.0);
      end;
      Cell_DATwiss(0, globval.Cell_nLoc, DAAscr, chroma, true, dP);
      globval.TotalTune := Cell[globval.Cell_nLoc].Nu;
      status.tuneflag := true;
      if chroma and not globval.Cavity_on then
      begin
	Ring_Getchrom(dP);
        Cell_GetCOD(globval.CODimax, globval.CODeps, dP, lastpos);
      end;
    end;
  end;


  procedure Ring_GetTwiss{chroma : boolean; dP : double};

  begin
    if trace then writeln('  enter ring_gettwiss');
    if globval.MatMeth then
      Ring_MatTwiss(chroma, dP)
    else
      Ring_DATwiss(chroma, dP);
    if trace then writeln('  exit ring_gettwiss');
  end;


PROCEDURE Ring_Fittune{var nu : vector2; eps : double;
		       var nq : ivector2; var qf, qd : fitvect;
		       dkL : double; imax : integer};

  { Fit tune }

  label 999;

  const	dP = 0.0;

  var	i, j, k, lastpos	: integer;
	nu0, nu1		: vector2;
	alpha, beta, gamma	: vector2;
	dkL1, dnu		: Vector;
	A			: Matrix;

  procedure shiftk(Elnum : integer; dk : double);

  begin
    with Cell[Elnum], Elem, M^ do
    begin
      PBpar[quad] := PBpar[quad] + dk;
      mpole_setPB(Fnum, Knum, quad);
    end;
  end;


  procedure checkifstable;

  begin
    if not stable then
    begin
      writeln('  lattice is unstable');
      goto 999;
    end;
  end;

  begin
    if trace then
      writeln('  Tune fit, nux =', nu[1]:10:5, ', nuy =', nu[2]:10:5,
	      ', eps =', eps:10, ', imax =', imax:4, ', dkL = ', dkL);
    Ring_GetTwiss(false, dP);
    checkifstable;
    nu0 := globval.TotalTune;
    i := 0;
    repeat
      i := succ(i);
      { First vary kf then kd }
      for j:=1 to 2 do
      begin
	for k:=1 to nq[j] do
        begin
	  if j = 1 then
	    Shiftk(qf[k], dkL)
	  else
	    Shiftk(qd[k], dkL);
        end;
        Cell_GetCOD(globval.CODimax, globval.CODeps, dP, lastpos);
	Cell_GetABGN(globval.OneTurnmat, alpha, beta, gamma, nu1);
	checkifstable;
	for k:=1 to 2 do
	begin
          dnu[k] := (nu1[k] - trunc(nu1[k])) - (nu0[k] - trunc(nu0[k]));
	  if abs(dnu[k]) > 0.5 then dnu[k] := 1 - abs(dnu[k]);
          A[k, j] := dnu[k]/dkL;
	end;
	{ Restore strength }
	for k:=1 to nq[j] do
	  if j = 1 then
	    Shiftk(qf[k], -dkL)
	  else
	    Shiftk(qd[k], -dkL);
      end;
      if not InvMat(2, A) then
      begin
	writeln('  A is singular');
	goto 999;
      end;
      for j:=1 to 2 do
        dkL1[j] := nu[j] - nu0[j];
      LinTrans(2, A, dkL1);
      for j:=1 to 2 do
	for k:=1 to nq[j] do
	  if j = 1 then
	    Shiftk(qf[k], dkL1[j])
	  else
	    Shiftk(qd[k], dkL1[j]);
      Ring_GetTwiss(false, dP);
      checkifstable;
      nu0 := globval.TotalTune;
      if trace then
        writeln('  Nux = ', nu0[1]:10:6, nu1[1]:10:6,
	        ', Nuy = ', nu0[2]:10:6, nu1[2]:10:6, 
      	        ', QF*L = ', Elem_getkval(Cell[qf[1]].Fnum, 1, quad),
		', QD*L = ', Elem_getkval(Cell[qd[1]].Fnum, 1, quad),
		' @', i:3);
    until ( sqrt( sqr(nu[1]-nu0[1]) + sqr(nu[2]-nu0[2]) ) < eps ) or
	  ( i = imax ) ;
999:
  end;


PROCEDURE Ring_Fitchrom{var ksi : vector2; eps : double;
			var ns : ivector2; var sf, sd : fitvect;
		        dkpL : double; imax : integer};

  label 999;

  const	dP = 0.0;

  var	i, j, k, lastpos	: integer;
	ksi0 			: vector2;
	dkpL1, dksi		: Vector;
	A			: Matrix;
	rad			: boolean;

  procedure shiftkp(Elnum : integer; dkp : double);

  begin
    with Cell[Elnum], Elem, M^ do
    begin
      PBpar[sext] := PBpar[sext] + dkp;
      mpole_setPB(Fnum, Knum, sext);
    end;
  end;


  begin
    if trace then
      writeln('  Chromaticity fit, ksix =', ksi[1]:10:5, ', ksiy =', ksi[2]:10:5,
	      ', eps =', eps:10, ', imax =', imax:4, ', dkpL =', dkpL:10:5);
    { Turn off radiation }
    rad := globval.radiation; globval.radiation := false;
    Cell_GetCOD(globval.CODimax, globval.CODeps, dP, lastpos);
    Ring_GetChrom(dP);
    for j:=1 to 2 do
      ksi0[j] := globval.Chrom[j];
    i := 0;
    repeat
      i := succ(i);
      { First vary sf then sd }
      for j:=1 to 2 do
      begin
	for k:=1 to ns[j] do
	  if j = 1 then
	    Shiftkp(sf[k], dkpL)
	  else
	    Shiftkp(sd[k], dkpL);
        Cell_GetCOD(globval.CODimax, globval.CODeps, dP, lastpos);
        Ring_GetChrom(dP);
        for k:=1 to 2 do
        begin
          dksi[k] := globval.Chrom[k] - ksi0[k];
          A[k, j] := dksi[k]/dkpL;
        end;
	{ Restore strength }
	for k:=1 to ns[j] do
	  if j = 1 then
	    Shiftkp(sf[k], -dkpL)
	  else
	    Shiftkp(sd[k], -dkpL);
      end;
      if not InvMat(2, A) then
      begin
	writeln('  A is singular');
	goto 999;
      end;
      for j:=1 to 2 do
	dkpL1[j] := ksi[j] - ksi0[j];
      LinTrans(2, A, dkpL1);
      for j:=1 to 2 do
	for k:=1 to ns[j] do
	  if j = 1 then
	    Shiftkp(sf[k], dkpL1[j])
	  else
	    Shiftkp(sd[k], dkpL1[j]);
      Cell_GetCOD(globval.CODimax, globval.CODeps, dP, lastpos);
      Ring_GetChrom(dP);
      for j:=1 to 2 do
        ksi0[j] := globval.Chrom[j];
      if trace then
        writeln('  ksix =', ksi0[1]:10:6, ', ksiy =', ksi0[2]:10:6, 
      	        ', SF = ', Elem_getkval(Cell[sf[1]].Fnum, 1, sext),
		', SD = ', Elem_getkval(Cell[sd[1]].Fnum, 1, sext), ' @', i:3);
    until ( sqrt( sqr(ksi[1]-ksi0[1]) + sqr(ksi[2]-ksi0[2]) ) < eps ) or
	  ( i = imax ) ;
999:
    { Restore radiation }
    globval.radiation := rad;
  end;


PROCEDURE Ring_FitDisp{pos : integer; eta, eps : double;
		       nq : integer; var q : fitvect;
	 	       dkL : double; imax : integer};

  label 999;

  const	dP = 0.0;

  var	i, j			: integer;
	dkL1, Eta0, deta	: double;
	rad			: boolean;

  procedure shiftk(Elnum : integer; dk : double);

  begin
    with Cell[Elnum], Elem, M^ do
    begin
      PBpar[quad] := PBpar[quad] + dk;
      mpole_setPB(Fnum, Knum, quad);
    end;
  end;


  procedure checkifstable;

  begin
    if not stable then
    begin
      writeln('  lattice is unstable');
      goto 999;
    end;
  end;

  begin
    if trace then
      writeln('  Dispersion fit, etax =', eta:10:5, 
	      ', eps =', eps:10, ', imax =', imax:4, ', dkL =', dkL:10:5);
    { Turn off radiation }
    rad := globval.radiation; globval.radiation := false;
    Ring_GetTwiss(true, dP);
    checkifstable;
    Eta0 := Cell[pos].Eta[1];
    i := 0;
    while (abs(eta-eta0) > eps) and (i < imax) do
    begin
      i := succ(i);
      for j:=1 to nq do
	Shiftk(q[j], dkL);
      Ring_GetTwiss(true, dP);
      checkifstable;
      deta := Cell[pos].Eta[1] - Eta0;
      if deta <> 0.0 then
        dkL1:= (Eta-Eta0)*dkL/deta - dkL
      else
      begin
	writeln('  deta is 0');
	goto 999;
      end;
      for j:=1 to nq do
	Shiftk(q[j], dkL1);
      Ring_GetTwiss(true, dP);
      checkifstable;
      eta0 := Cell[pos].Eta[1];
      if trace then
        writeln('  Dispersion = ', eta0,
	        ', kL =', Elem_getkval(Cell[q[1]].Fnum, 1, quad), ' @', i:3);
    end;
999: 
    { Restore radiation }
    globval.radiation := rad;
  end;

  end.
