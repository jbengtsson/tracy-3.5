  module t2lat(input, output);

  %include 'mathlib.def'
  %include 'dab.def'

  %include 'pascommon.def'
  %include 't2common.def'
  %include 't2elem.def'
  %include 't2lat.def'

{ interface }

[global] Function Lattice_Read(var fi : text; { lattice input   file  }
			       var fo : text  { lattice message file  }
			       ):boolean; forward;

{ implementation }

Function Lattice_Read;

label 9999;

const	NoBmax  = 60;		{ maximum number of blocks (NoB) }
	NoBEmax = 10000;	{ maximum number of elements in a block (Elem_NFam) }
	UDImax  = 300;
	LatLLng = 81;
	Lat_nkw_max = 100;	{no. of key words}

	{ tables }

	emax = 322;
	emin = -292;
	kmax = 15;		{max no. of signIFicant digitis}
	smax = 600;		{size   of string table}
	xmax = maxint;		{2**8 - 1}
	nmax = maxint;

type	Latlinetype	= packed array [1..latllng] of char;
	Lat_symbol	= ( bndsym, defsym, dfcsym, drfsym, elmsym, exesym, fcssym, horsym,
              monsym, qdsym, sexsym, versym, plus, minus, lparent, rparent,
              eql, comma, lbrack, rbrack, neq, andsy, semicolon,
              times, rdiv, intcon, realcon, becomes, colon, leq, pwrsym,
              lss, geo, gtr, period, charcon, stringcon, ident,
              geq, lsym, bobrhosym, kxsym, ksym, tsym, t1sym, t2sym,
	      gapsym, thksym, invsym,
              thnsym, endsym, extsym, tsksym, bemsym, corsym, prnsym,
              tblsym, possym, prmsym, udisym, squote, linsym, mthsym,
              celsym, matsym, cavsym, symsym, chmsym, cctsym,
              usesym, andsym, dspsym, kicksym, wglsym, nsym,
              mrksym, b0sym, nbdsym, frgsym, latsym, mpsym, dbnsym,
              kssym, homsym, lmdsym, dtsym, xytsym,
              vrfsym, harnumsym, frqsym, gstsym, typsym, tltsym ) ;{\}

	symset		= set of Lat_symbol;

	Lat_keytype	= array[1..Lat_nkw_max] of alfa;
	Lat_ksytype	= array[1..Lat_nkw_max] of Lat_symbol;
	Lat_spstype	= array[char] of Lat_symbol;

	BlockStype	= array[1..NoBmax] of RECORD
						Bname		: partsname; { name of a beam line }
						BSTART, BOWARI	: integer;
					      END;      

	index		=    - xmax .. +  xmax;
	types		=   (notyp, ints, reals, bools, chars, arrays, records);
	typset		=   set of types;

  var	Symmetry	: integer; 
	Ring		: boolean;	{ true is CELL is a ring }

	NoB	 	: integer;	{ Number of defined Blocks }
	BlockS		: BlockStype;

	Bstack		: array[1..NoBEmax] of integer;
	Bpointer	: integer;

	UDIC		: integer;	{ Number of user defined constants }
	UDItable	: array[1..udimax] of RECORD
						Uname	: partsname;
						Uvalue	: double;
					      END;

	nkw		:  integer;       { number of key word }
	sym		:  Lat_symbol;        {last   symbol read by GetSym}
	id		:  alfa;          {identIFier freom GetSym}
	inum		:  integer;       {integer from GetSym}
	rnum		:  double;          {double   number from GetSym}
	chin		:  char;          {last   character read from source program}
	cc		:  integer;       {character counter}
	lc		:  integer;       {program location counter}
	ll		:  integer;       {length of current line}
	errpos		:  integer;
	line		: LatLineType;

	key		: Lat_Keytype;
	ksy		: Lat_Ksytype;
	sps		: Lat_spstype;

	defbegsys	: symset;
	elmbegsys	: symset;
	skipflag	: boolean;
	rsvwd		: boolean;


   FUNCTION rtot(x:double) : double; forward;

   FUNCTION CheckElementtable(name : partsname) : integer; forward;

   function CheckBLOCKStable(name : partsname) : integer; forward;


   PROCEDURE InitUDItable; forward;

   FUNCTION CheckUDItable(name : partsname) : integer; forward;

   PROCEDURE EnterUDItable(name : partsname; X : double); forward;

   PROCEDURE ModUDItable(N : integer; X : double); forward;

   PROCEDURE RefUDItable(name : partsname; var X : double); forward;

   procedure PrintUpname(var name:partsname);forward;

   procedure PrintUpname1(var name:partsname;var pos:integer);forward;

   procedure PrintUpname2(var name:partsname);forward;

   PROCEDURE abort;forward;

   PROCEDURE ENDskip(var fo:text; var errpos, cc:integer;
                     var skipflag:boolean);forward;

   PROCEDURE Lat_Error(n:integer; var fo:text; var cc, errpos: integer);forward;

   PROCEDURE Lat_Nextch(var fi, fo:text;
                        var cc, ll, errpos, lc:integer;
                        var chin:char;
                        var  skipflag:boolean;
                        var line:latlinetype);forward;

   PROCEDURE Lat_errorm(cmnt:str80; var fi, fo:text;
              var cc, ll, errpos, lc:integer;
              var chin:char;
              var  skipflag:boolean;
              var line:latlinetype);forward;

  PROCEDURE Lat_GetSym(var fi, fo:text;
              var cc, ll, errpos, lc, nkw, inum:integer;
                  emax, emin, kmax, nmax:integer;
              var chin:char;
              var id:alfa;
              var rnum:double;
              var  skipflag, rsvwd:boolean;
              var line:latlinetype;
              var sym:Lat_symbol;
              var key:Lat_keytype;
              var ksy:Lat_ksytype;
              var sps:Lat_spstype);forward;

  function Lat_EVAL(var fi, fo:text;
              var cc, ll, errpos, lc, nkw, inum:integer;
                  emax, emin, kmax, nmax:integer;
              var chin:char;
              var id:alfa;
              var rnum:double;
              var  skipflag, rsvwd:boolean;
              var line:latlinetype;
              var sym:Lat_symbol;
              var key:Lat_keytype;
              var ksy:Lat_ksytype;
              var sps:Lat_spstype):double;forward;

 procedure Lat_ProcessBlockInput(var fi, fo:text;
              var cc, ll, errpos, lc, nkw, inum:integer;
                  emax, emin, kmax, nmax:integer;
              var chin:char;
              var id:alfa;
              var BlockName:partsname;
              var rnum:double;
              var  skipflag, rsvwd:boolean;
              var line:latlinetype;
              var sym:Lat_symbol;
              var key:Lat_keytype;
              var ksy:Lat_ksytype;
              var sps:Lat_spstype);forward;

  function Lat_CheckWiggler(var fo:text;i:integer):boolean;forward;

  function Lat_DealElement(
              var fi, fo:text;
              var cc, ll, errpos, lc, nkw, inum:integer;
                  emax, emin, kmax, nmax:integer;
              var chin:char;
              var id:alfa;
              var ElementName, BlockName:partsname;
              var rnum:double;
              var  skipflag, rsvwd:boolean;
              var line:latlinetype;
              var sym:Lat_symbol;
              var key:Lat_keytype;
              var ksy:Lat_ksytype;
              var sps:Lat_spstype):boolean;forward;


  FUNCTION rtot{(x:double  ):double } ;

  begin
    rtot := x/6.2831853
  end;


  FUNCTION CheckElementtable{name : partsname}{ : integer};

  var	i, j	: integer;

  BEGIN { globval.Elem_nFam = Number of parts in a Element }
     j := 0;
     if globval.Elem_nFam <= Elem_nFamMax then
     begin
       for i:=1 to globval.Elem_nFam do
         if ElemFam[i].ElemF.Pname = name then j:=i;
     end
     else
       writeln('Elem_nFamMax exceeded: ', globval.Elem_nFam:1, 
	       '(', Elem_nFamMax:1, ')');
     CheckElementtable:=j
  END;


  function CheckBLOCKStable{name : partsname}{ : integer};

  var i, j	: integer;

  BEGIN { NoB = Number of Block defs }
    j := 0;
    if NoB <= NoBmax then
    begin
      for i:=1 to NoB do
        if BLOCKS[i].Bname = name then j:=i;
    end
    else
      writeln('** NoBmax exhausted: ', NoB:1, '(', NoBmax:1, ')');
    CheckBLOCKStable := j
  END;


  PROCEDURE InitUDItable;

  BEGIN
    udic := 0
  END;


  FUNCTION CheckUDItable{name : partsname}{ : integer};

  var	i, j	: integer;

  BEGIN
    j := 0;
    if udic <= UDImax then
    begin
      for i:=1 to udic do
        if UDItable[i].Uname = name then j := i;
    end
    else
      writeln('** UDImax exhausted: ', udic:1, '(', UDImax:1, ')'); 
    CheckUDItable := j;
  END;


  PROCEDURE EnterUDItable{name : partsname; X : double};

  BEGIN
    udic := udic + 1;
    if udic <= UDImax then
    begin
      with UDItable[udic] do
      BEGIN
        Uname := name; Uvalue := x;
      END;
    end
    else
      writeln('** UDImax exhausted: ', udic:1, '(', UDImax:1, ')'); 
   END;


   PROCEDURE ModUDItable{ (N : integer; X : double};

   BEGIN
     with UDItable[N] do
       Uvalue:=x;
   END;


   PROCEDURE RefUDItable{(name:partsname;var X:double)};

   var	k	: integer;

   BEGIN
     k := CheckUDItable(name); X := UDItable[k].Uvalue;
   END;  


   procedure PrintUpname{(var name:partsname)};

   var i	: integer;
       ch	: char;

   begin
     for i:=1 to NameLength do
     begin
      ch := name[i];
      if (ord('a') <= ord(ch)) and (ord(ch) <= ord('z')) then
        ch := chr(ord(ch) - ord('a') + ord('A'));
      write(fo, ch:1);
     end;
   end;


   procedure PrintUpname1{var name : partsname; var pos : integer};

   var	i	: integer;
	ch	: char;

   begin{1}
     pos:=0;
     for i:=1 to NameLength do
     begin{2}
      ch:=name[i];
      if (ord('a')<=ord(ch)) and (ord(ch)<=ord('z')) then
                           ch:=chr(ord(ch)-ord('a')+ord('A'));
      if ch<>' ' then begin write(fo, ch:1);pos:=pos+1 end;
     end;{2}
   end;{1}

  procedure PrintUpname2{(var name:partsname)};

  var	i	: integer;
	ch	: char;

  begin
    for i:=1 to NameLength do
    begin
      ch:=name[i];
      if (ord('a') <= ord(ch)) and (ord(ch) <= ord('z')) then
	ch:=chr(ord(ch)-ord('a')+ord('A'));
      write(fo, ch:1);
    end;
  end;

  PROCEDURE abort;

  const	n = 3;

  var	i	: integer;

  BEGIN
    for i:=1 to n do write(chr(7));
    writeln;
    writeln('>>>>> error detected in the lattice file <<<<<<');
    writeln;
    ErrFlag:=true;
    {goto 9999}
    writeln(SQRT(1-2));
  END;

  PROCEDURE ENDskip{var fo:text; var errpos, cc:integer;
                    var skipflag:boolean};

  BEGIN {underline skips part of input}
    while errpos < cc do
    BEGIN
      write(fo, '-'); 
      errpos := errpos + 1
    END;
    skipflag  := false
  END;
                          

  PROCEDURE Lat_Error{ n:integer; var fo:text;
                       var cc, errpos: integer};
  BEGIN
    IF errpos = 0 THEN {write(fo, ' ****')}
    IF cc > errpos THEN
    BEGIN
      write(fo, ' ':cc-errpos, '^', ord(n):2);
      errpos := cc + 3;
    END
  END;


   PROCEDURE Lat_Nextch{var fi, fo:text;
                        var cc, ll, errpos, lc:integer;
                        var chin:char;
                        var  skipflag:boolean;
                        var line:latlinetype};
  {read   next character; process line END}

   BEGIN { Lat_Nextch }
     IF cc = ll THEN
     BEGIN
       IF eof(fi) THEN
       BEGIN
         writeln(fo);
         writeln(fo, 'program incomplete');
         {errormsg;}
         abort;
       END;
       IF errpos <> 0 THEN
       BEGIN
         IF skipflag THEN ENDskip(fo, errpos, cc, skipflag); 
         writeln(fo);
         errpos := 0
       END;
       { write(fo, }{lc: 5, }{ '  ');}
       lc:=lc+1;

       ll := 0; cc := 0;
       while not eoln(fi) do
       begin
         ll := succ(ll);
         read(fi, chin); write(fo, chin);
	 line[ll] := chin;
       end;
       ll := succ(ll);
       readln(fi); line[ll] := ' ';
       {read(fi, line[ll]);} writeln(fo);
     end;
     cc := cc + 1;
     chin := line[cc];
     { upper case to lower case }
     if chin in ['A'..'Z'] then chin := chr(ord(chin) - ord('A') + ord('a'));
     { tab }
     if chin = chr(9) then chin := ' ';
   END; { Lat_Nextch }

   PROCEDURE Lat_errorm;
   BEGIN
      IF errpos = 0 THEN
      BEGIN  {write(fo, ' ****')} END;
      IF cc > errpos THEN
          BEGIN
             write(fo, ' ':cc-errpos, '^', cmnt);
             errpos := cc + 3;
          END;
      while not eof(fi) do 
          Lat_Nextch(fi, fo, cc, ll, errpos, lc, chin, skipflag, line);
      ErrFlag:=true;
      Abort;
    END;


  PROCEDURE Lat_GetSym;

  label  1;

  var	i, j, k, e, mysign	: integer;
	parsename		: boolean;

  procedure NextCh;

  begin
    Lat_Nextch(fi, fo, cc, ll, errpos, lc, chin, skipflag, line);
  end;

  PROCEDURE readscale;

  var	s, sign	: integer;

  BEGIN 
    Nextch;
    while chin=' ' do Nextch;
    sign:=1; s:=0;
    IF chin='+' THEN Nextch
    ELSE  IF chin='-' THEN
    BEGIN
      Nextch; sign:=-1;
    END;
    IF not ((chin>='0') and (chin<='9')) THEN
      Lat_Error(40, fo, cc, errpos)
    ELSE
      REPEAT
        s:=10*s+ord(chin)-ord('0');
        Nextch;
      until not ((chin>='0') and (chin<='9'));
    e:=s*sign+e
  END  { readscale  };


  PROCEDURE adjustscale;

  var	s	: integer;
	d, t	: double;

  BEGIN
    IF k+e > emax THEN
      Lat_Error(21, fo, cc, errpos)
    ELSE IF k+e < emin THEN
      rnum := 0.0
    ELSE
    BEGIN
      s := abs(e); t := 1.0; d := 10.0;
      REPEAT
        while not odd(s)  do
        BEGIN
          s := s div 2; d := sqr(d)
        END;
        s := s-1; t := d*t
      until s=0;

      IF e>=0 THEN rnum:=rnum*t ELSE
	rnum:=rnum/t
    END
  END  {  adjustscale  };


   BEGIN {GetSym}
      rsvwd := false; mysign := 1; parsename := false;
   1: while chin = ' ' do Nextch;
      case chin of
         'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 
         'j', 'k', 'l', 'm', 
         'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 
         'w', 'x', 'y', 'z', '"':
             BEGIN {identIFier or wordsymbol}
                k := 0; id := '          ';
                REPEAT
		   if chin = '"' then parsename := not parsename;
                   IF k < NameLength THEN
                   BEGIN
                     k := k + 1; id[k] := chin
                   END;
                   Nextch;
                until not parsename and
		      not (chin in ['a' .. 'z', '0' .. '9', '_']) ;

                {writeln(fo, 'GetSym detected: id=', id);}
                
                i := 1;
                j := nkw;               {binary search}
                REPEAT
                   k := (i + j) div 2;
                   IF id <= key[k] THEN
                      j := k - 1;
                   IF id >= key[k] THEN
                      i := k + 1;
                 {  writeln(fo, '  bunary: id=', id, '  key[', k:3, ']=', key[k], 
                           'i=', i:4, ' j=', j:4, ' k=', k:4, ' i-1-j=', (i-1-j):4);}
                until i > j;
                IF i - 1 > j
                THEN 
                  BEGIN
                    sym := ksy[k];
                    rsvwd := true;
                  {  writeln(fo, 'GetSym detected reserved word: id=', id, 
                                '  k=', k:4, '  key[', k:4, ']=', key[k]);}
                  END
                ELSE
                BEGIN
                  IF id = 't         ' THEN
		    sym := tsym
		  ELSE IF id = 'gap       ' THEN
		    sym := gapsym
		  ELSE IF id = 'l         ' THEN
		    sym := lsym
		  ELSE IF id = 'n         ' THEN
		    sym := nsym
		  ELSE IF id = 'bobrho    ' THEN
		    sym := bobrhosym
		  ELSE IF id = 'kx        ' THEN
		    sym := kxsym
		  ELSE IF id = 'k         ' THEN
		    sym := ksym
		  ELSE IF id = 'harnum    ' THEN
		    sym := harnumsym
		  ELSE
                    sym := ident
                END
             END;
          '0', '1', '2', '3', '4', '5', '6', '7', '8', '9':
             BEGIN {number}
                k := 0;
                inum := 0;
                sym := intcon;
                REPEAT
                   inum := inum * 10 + ord(chin) - ord('0');
                   k := k + 1;
                   Nextch;
                until not (chin in ['0' .. '9']);
                IF (k > kmax) or (inum > nmax) THEN
                   BEGIN
                      Lat_Error(21, fo, cc, errpos);
                      inum := 0;
                      k := 0
                   END;
                 IF chin='.'
                 THEN BEGIN Nextch;
                    IF chin='.'
                    THEN chin:=':'
                    ELSE BEGIN
                      sym:=realcon; rnum:=inum; e:=0;
                      while (chin >= '0') and (chin <= '9' ) do
                      BEGIN
                        e:=e-1;
                        rnum := 10.0*rnum + (ord(chin)-ord('0'));
                        Nextch;
                      END;
                      while chin=' ' do Nextch;
                      IF e = 0 THEN Lat_Error(40, fo, cc, errpos);
                      if (chin = 'd') or (chin = 'D') or
			 (chin = 'e') or (chin = 'E') then readscale;
                      IF e <> 0 THEN adjustscale;
                    END
                  END ELSE
                    if (chin = 'd') or (chin = 'D') or
		       (chin = 'e') or (chin = 'E') then
		    begin
                      sym := realcon; rnum := inum; e := 0;
                      readscale;
                      IF e <> 0 THEN adjustscale
                    END;
                  IF sym=intcon THEN inum:=mysign*inum ELSE
                  IF sym=realcon THEN rnum:=mysign*rnum;
             END;

          ':' {, col}:
             BEGIN
                Nextch;
                IF chin = '='
                THEN
                   BEGIN
                      sym := becomes;
                      Nextch;
                   END
                ELSE
                   sym := colon
             END;
          '<':
             BEGIN
                Nextch;
                IF chin = '='
                THEN
                   BEGIN
                      sym := leq;
                      Nextch;
                   END
                ELSE
                   IF chin = '>'
                   THEN
                      BEGIN
                         sym := neq;
                         Nextch;
                      END
                   ELSE
                      sym := lss
             END;
          '>':
             BEGIN
                Nextch;
                IF chin = '='
                THEN
                   BEGIN
                      sym := geq;
                      Nextch;
                   END
                ELSE
                   sym := gtr
             END;
          '.':
             BEGIN
                Nextch;
                sym := period
             END;
          '*':
             BEGIN
                Nextch;
                IF chin='*' THEN
                   BEGIN
                     sym:=pwrsym;
                     Nextch;
                   END
                ELSE sym := times;
             END;
	  '{' : begin
                  repeat
		    nextch;
                  until chin = '}';
                  nextch;
                  goto 1
		end;
          '+', '-', '/', '(', ')', '=', ',', ';', '[', ']', '''' :
             BEGIN
                sym := sps[chin];
               { IF chin='+' THEN BEGIN nextch; goto 1 END ELSE
               IF chin='-' THEN BEGIN nextch; mysign:=-1; goto 1 END ELSE}
                Nextch;
             END;
          '$', '!', '@', '?', '_', '&', '\', '^' :
             BEGIN
                Lat_Error(24, fo, cc, errpos);
                Nextch;
                goto 1
             END;
      otherwise
         BEGIN
           Lat_Error(24, fo, cc, errpos); Nextch;
           goto 1
         END;
      END
   END  {GetSym};



{*****************************
 *                           *
 *        E V A L            *
 *                           *
 *****************************}
    
  function Lat_EVAL{
              var fi, fo:text;
              var cc, ll, errpos, lc, nkw, inum:integer;
                  emax, emin, kmax, nmax:integer;
              var chin:char;
              var id:alfa;
              var rnum:double;
              var  skipflag, rsvwd:boolean;
              var line:latlinetype;
              var sym:Lat_symbol;
              var key:Lat_keytype;
              var ksy:Lat_ksytype;
              var sps:Lat_spstype):double};

    label 999, 888;
    const tmax=100;
         
    var S:array[0..tmax] of double;
        t:integer;
        facbegsys:symset;


  PROCEDURE GetSym;     { reads next symbol  }

  BEGIN
    Lat_GetSym(fi, fo, cc, ll, errpos, lc, nkw, inum, emax, emin, kmax, nmax, 
               chin, id, rnum, skipflag, rsvwd, line, 
               sym, key, ksy, sps);
  END;

 
  PROCEDURE errorm(cmnt:str80);

  BEGIN
    Lat_errorm(cmnt, fi, fo, cc, ll, errpos, lc, chin, skipflag, line);
  END;


  PROCEDURE test(s1:symset; cmnt: str80);

  BEGIN
    IF not (sym in s1) THEN errorm(cmnt);
  END;


  PROCEDURE getest(s1:symset; cmnt: str80);
 
  BEGIN
    GetSym;
    IF not (sym in s1) THEN errorm(cmnt);
  END;


  function ArcSin(x:double):double;

  BEGIN
    IF abs(x)>1 THEN goto 999;
    IF x=1 THEN
      ArcSin:=Pi/2
    ELSE IF x=-1 THEN
      ArcSin:=-Pi/2
    ELSE
      ArcSin:=ArcTan(x/sqrt(1-x*x));
  END;


  function ArcCos(x:double):double;

  BEGIN
    IF abs(x)>1 THEN goto 999;
    IF x=1 THEN
      ArcCos:=0
    ELSE IF x=-1 THEN
      ArcCos:=Pi
    ELSE
      ArcCos:=ArcTan(sqrt(1-x*x)/x);
  END;


  procedure writes;

  BEGIN
    {writeln('PUSH:  s[', t:3, ']=', s[t]);}
  END;


  procedure PUSH(x:double);

  BEGIN
    t:=t+1;
    IF t=tmax THEN
    BEGIN
      writeln('** Lat_Eval: stack overflow');
      goto 999;
    END;
    S[t]:=x;
    writes;
  END;


    function BlockLength(ii:integer):double;
     var k1, k2, k3:integer;
         S:double;
     begin
       S:=0.0;
       IF ii<>0 THEN  
         BEGIN
           k2:=Blocks[ii].Bstart ;
           k3:=Blocks[ii].Bowari ;
           FOR k1:=k2 to k3 do 
           BEGIN
             S:=S+ElemFam[Bstack[k1]].ElemF.PL;
           END;
         END;
       BlockLength:=S;
     end;

    procedure Expression;
    var addop: Lat_symbol;

      procedure Term;
      var mulop: Lat_symbol;

        procedure Factor;
        var i: integer; x:double;
            fname: partsname;

           function GetKparm(direction:integer):double;
             begin
               Getest([lbrack], '<[> expected');
               Getsym;
               Expression;
               test(  [rbrack], '<]> expected');
               if direction =1 
                  then GetKparm:=ElemFam[i].ElemF.M^.PBpar[trunc(S[t])]
                  else GetKparm:=ElemFam[i].ElemF.M^.PBpar[-trunc(S[t])];
               t:=t-1;
              { GetSym;}
             end;

           BEGIN {  factor  }
             IF not (sym in facbegsys) THEN goto 999;
             BEGIN
               {while sym in facbegsys do}
               IF sym = ident THEN {***************************}
               BEGIN{1:  of ident } 
                 i := CheckUDItable(id);
                 IF i <> 0 THEN { UDI }
                 BEGIN
                   x := UDItable[i].Uvalue; PUSH(x); GetSym;
                 END
                 ELSE
                 BEGIN
                   i:=CheckElementTable(id);
                   IF i<>0 THEN
                   BEGIN
                     getest([period], '<.> expected');
                     {--> new }
                     getsym;
                     fname:=id;
		     with ElemFam[i].ElemF, M^ do
		     begin
                       IF fname='l         ' THEN 
                         x := PL
		       ELSE IF fname='t         ' THEN 
		       begin
		         if PL <> 0.0 then
			   x := Pirho*PL*180/pi
		         else
			   x := Pirho*180/pi;
		       end
		       ELSE IF fname='t1        ' THEN 
                         x := PTx1
		       ELSE IF fname='t2        ' THEN 
                         x := PTx2
		       ELSE IF fname='gap       ' THEN 
                         x := Pgap
		       ELSE IF fname='tilt      ' THEN 
                       begin
                         if Pkind=mpole then
			   x := PdTpar
                         else if Pkind=wigl then
			   x := PdTpar;
                       end
		       ELSE IF fname='n         ' THEN 
                       begin
                         if Pkind in [mpole, wigl] then x := PN;
                       end
		       ELSE IF fname='b         ' THEN 
                         x := GetKparm(1)
		       ELSE IF fname='a         ' THEN 
                         x := GetKparm(2)
		       ELSE
                         { error detected }
                         Getest([], '  illegal extension...');
		     end;
                     PUSH(x); GetSym; 
                   END
		   else
                   begin
                     i:=CheckBLOCKStable(id);
                     IF i <> 0 THEN
                     BEGIN
                       getest([period], '<.> expected');
                       getest([lsym], 'illegal component');
                       x := BlockLength(i); 
                       PUSH(x);
                       GetSym; 
                     END
                     ELSE { function ? }
                     BEGIN{4: function ?}
                       fname:=id;                 
                       GetSym;
                       case sym of{5}
                         semicolon: GetSym;
                         lparent  : BEGIN{6: of lparent}
	                              GetSym;
        	                      Expression;
                             	      IF fname='sin       ' THEN
					S[t]:=sin(S[t])
				      ELSE IF fname='cos       ' THEN
					S[t]:=cos(S[t])   
				      ELSE IF fname='tan       ' THEN
					S[t]:=tan_(S[t])  
				      ELSE IF fname='arcsin    ' THEN
					S[t]:=arcsin(S[t])
				      ELSE IF fname='arccos    ' THEN
					S[t]:=arccos(S[t])
				      ELSE IF fname='arctan    ' THEN
					S[t]:=arctan(S[t])
				      ELSE IF fname='sqrt      ' THEN
					S[t]:=sqrt(S[t])  
				      ELSE IF fname='log       ' THEN
					S[t]:=ln(S[t])    
				      ELSE IF fname='exp       ' THEN
					S[t]:=exp(S[t]);
	                              writes;
                                      IF sym=rparent THEN GetSym ELSE goto 999;
                                    END;{6:of lparent}
                       END;{5: of case }
                     END{4: of function?}
                   end
                 END
               END{1: of ident}
               ELSE IF sym = realcon THEN 
               BEGIN
                 PUSH(rnum); GetSym;
               END

                ELSE IF sym = intcon THEN
                        BEGIN
                          x:=inum;
                          PUSH(X);
                          GetSym;
                        END

                ELSE IF sym = lparent THEN
                  BEGIN 
                     GetSym; 
                     expression;
                     IF sym = rparent THEN GetSym 
                                      ELSE goto 999;
                  END ELSE goto 999;

                IF sym=pwrsym THEN
                  BEGIN
                      GetSym;
                      IF sym<>intcon THEN goto 999;
                      S[t]:=PWR(S[t], inum);
                      GetSym;
                  END;
     
             END
          END { factor };

      BEGIN { term }
          factor;
          while sym in [times, rdiv] do
           BEGIN mulop:=sym; 
                 GetSym;
                 factor;
                 IF mulop = times THEN 
                    BEGIN
                       s[t-1]:=s[t-1]*s[t];t:=t-1;
                        writes;
                    END ELSE
                 IF mulop = rdiv THEN
                    BEGIN
                       s[t-1]:=s[t-1]/s[t];t:=t-1;
                        writes;
                    END;
           END
      END { term };

    BEGIN { Expression }

       IF sym in [plus, minus] THEN
          BEGIN 
            addop:=sym; 
            GetSym;
            term;
            IF addop = minus THEN S[t]:=-S[t];
          END ELSE term;

       while sym in [plus, minus] do
         BEGIN 
            addop:=sym;
            GetSym;
            term;
            IF addop = plus  THEN 
                    BEGIN
                      s[t-1]:=s[t-1]+s[t];
                      t:=t-1;
                      writes;
                    END ELSE
            IF addop = minus THEN
                    BEGIN
                      s[t-1]:=s[t-1]-s[t];
                      t:=t-1;
                      writes;
                    END;
         END

    END { Expression };

  BEGIN { eval }
    facbegsys:=[intcon, realcon, ident, lparent];
    GetSym;
    t:=0;
    Expression;
    IF t=1 THEN Lat_EVAL:=s[1]
           ELSE goto 999;
    goto 888;

  999: BEGIN
          ErrFlag:=true;
          TEST([], '** Lat_Eval: error');
       END;
  888: { exit }
    
  END;


 procedure Lat_ProcessBlockInput;

   var i:integer;
   

   PROCEDURE errorm(cmnt:str80);
     BEGIN
      Lat_errorm(cmnt, fi, fo, cc, ll, errpos, lc, chin, skipflag, line);
     END;

   procedure GetSym;
    begin
      Lat_GetSym(fi, fo, cc, ll, errpos, lc, nkw, inum, emax, emin, kmax, nmax, 
               chin, id, rnum, skipflag, rsvwd, line, 
               sym, key, ksy, sps);
    end;

   PROCEDURE test(s1:symset; cmnt: str80);
      BEGIN
          IF not (sym in s1) THEN errorm(cmnt);
      END {test};

   PROCEDURE getest(s1:symset; cmnt: str80);
      BEGIN  GetSym;
          IF not (sym in s1) THEN errorm(cmnt);
      END {test};

  function EVAL:double;
   begin
    Eval:=  Lat_EVAL(fi, fo, cc, ll, errpos, lc, nkw, inum, emax, emin, kmax, nmax, 
               chin, id, rnum, skipflag, rsvwd, 
               line, sym, key, ksy, sps);
   end;

   procedure DeBlock(ii, k4:integer);
     var k1, k2, k3, k5:integer;
     begin
       k2:=Blocks[ii].Bstart ;
       k3:=Blocks[ii].Bowari ;
       for k5:=1  to k4 do
       FOR k1:=k2 to k3 do 
       BEGIN{11}
         Bpointer:=Bpointer+1;
         Bstack[Bpointer]:=Bstack[k1];
       END;{11}
       GetSym;IF sym=comma THEN GetSym;
     end;
  
   procedure GetBlock;
      var i, ii, k1, k4:integer;

        procedure InsideParent(k4:integer);
          var b, b1, b2, k1:integer;
          begin
              b1:=bpointer+1;
              GetSym;
	      GetBlock;
              b2:=bpointer;
              if k4>=2 then 
              begin
                for k1:=2 to k4 do
                for b:=b1 to b2 do
                begin
                  Bpointer:=Bpointer+1;
                  Bstack[Bpointer]:=Bstack[b];
                end;
              end;
	      test([rparent], '<)> expected');
	      Getest([comma, semicolon, rparent], '<, > or <;> expected');
              if sym=comma then GetSym;
          end;

        procedure Doinverse;
          var b, b1, b2, b3, k1, k4:integer;
          begin
              Getest([lparent], '<(> expected after INV');
              b1:=bpointer+1;
              GetSym;
	      GetBlock;
              b2:=bpointer;
              k4:=b2-b1;
              if k4>=2 then 
              begin
                k4:=k4 div 2;
                b3:=b1+k4;
                k1:=0;
                for b:=b1 to b3 do
                begin
                  b3:=Bstack[b];
                  Bstack[b]:=Bstack[b2-k1];
                  Bstack[b2-k1]:=b3;
                  k1:=k1+1;
                end;
              end;
	      test([rparent], '<)> expected');
	      Getest([comma, semicolon, rparent], '<, > or <;> expected');
              if sym=comma then GetSym;
          end;

      begin
          REPEAT{7}
           test([ident, intcon, lparent, invsym], 
	        '<Element/Block name>, <integer>, <INV> or <(> expected');
	   if sym=lparent then
	   begin
             InsideParent(1);
	   end else
           if sym=invsym then
           begin
             Doinverse;
           end else
           IF sym=ident THEN 
           BEGIN{8}
             i:=CheckElementtable(id);
             IF i<>0 THEN
             BEGIN{9}
               Bpointer:=Bpointer+1; Bstack[Bpointer]:=i;
               GetSym;IF sym=comma THEN GetSym;
             END {9} ELSE
             BEGIN {9}
               ii:=CheckBlockstable(id);
               IF ii<>0 THEN  
               BEGIN{10}
                 DeBlock(ii, 1);
               END{10} ELSE
               BEGIN{10}
                 ii := CheckUDItable(id);
                 IF ii <> 0 THEN { UDI }
                 BEGIN{11}
                   K4 := Round(UDItable[ii].Uvalue); GetSym;
                 END{11} ELSE test([], 'invalid identIFier');
                 test([times], '<*> expected');
                 IF sym=times THEN
                 BEGIN{11}
                   Getest([ident, lparent, invsym], 
                        '<element/Block name>, <INV> or <(> expected');
		   if sym=lparent then
	           begin
	             InsideParent(k4);
	           end else
                   if sym=invsym then
                   begin
                     Doinverse;
                   end else
                   IF sym=ident THEN
                   BEGIN{12}
                     i:=CheckElementtable(id);
                     IF i<>0 THEN
                     BEGIN{13}
                       FOR k1:=1 to K4 do
                       BEGIN{14}
                         Bpointer:=Bpointer+1; Bstack[Bpointer]:=i;
                       END;{14}
                       GetSym;IF sym=comma THEN GetSym;
                     END {13} ELSE
                     BEGIN {13}
                       ii:=CheckBlockstable(id);
                       IF ii=0 THEN test([], 'invalid name');
                       DeBlock(ii, k4);
                     END;{13}
                    END;{12}
                  END;{11} 
                END;{10}
              END;{9}
            END{8} ELSE
            IF sym=intcon THEN 
            BEGIN{8}
              K4:=inum;
              GetSym;test([times], '<*> expected');
              IF sym=times THEN
              BEGIN{9}
                GetSym;test([ident, lparent, invsym], 
                           '<element/Block name>, <INV> or <(> expected');
    	        if sym=lparent then
	        begin
	          InsideParent(k4);
	        end else
                if sym=invsym then
                begin
                  Doinverse;
                end else
                IF sym=ident THEN
                BEGIN{10}
                  i:=CheckElementtable(id);
                  IF i<>0 THEN
                  BEGIN{11}
                    FOR k1:=1 to K4 do
                    BEGIN{12}
                      Bpointer:=Bpointer+1; Bstack[Bpointer]:=i;
                    END;{12}
                    GetSym;IF sym=comma THEN GetSym;
                  END {11} ELSE
                  BEGIN {11}
                    ii:=CheckBlockstable(id);
                    IF ii=0 THEN test([], 'invalid name');
                    DeBlock(ii, k4);
                  END;{11}
                END;{10}
              END;{9} 
            END{8} ELSE
           IF sym=minus THEN 
           BEGIN{8}
             Getsym;
             i:=CheckElementtable(id);
             IF i<>0 THEN
             BEGIN{9}
               Bpointer:=Bpointer+1; Bstack[Bpointer]:=-i;
               GetSym;IF sym=comma THEN GetSym;
             END {9} ELSE test([], '<element name> expected.');
            END;{8}
          UNTIL{7} not (sym in [ident, intcon, minus, invsym]);      
      end;
   
  BEGIN 
    i:=CheckElementtable(BlockName);
    IF i <> 0 THEN
      test([], '<Block name>: conflict with Element name')
    ELSE 
    BEGIN
      { Increment number of defined blocks }                 
      NoB := NoB + 1;
      if NoB <= NoBmax then
      begin
        with Blocks[NoB] do 
        BEGIN
          Bname := BlockName; Bstart := Bpointer + 1;
          GetBlock;
          test([semicolon], '<;> expected');
          GetSym;
          Bowari := Bpointer;                    
        END;
      end
      else
        writeln('** NoBmax exhausted: ', NoB:1, '(', NoBmax:1, ')');
    END;
  END;{ ProcessBlockInput }


  function Lat_CheckWiggler{i : integer}{ : boolean};
 
  var	a, Lambda, L, diff	: double;
	NN			: integer;

  begin
    Lat_CheckWiggler:=false;
    with ElemFam[i], ElemF, W^ do
    begin
      Lambda := Plperiod;
      L := PL;
      a := L/Lambda;
      NN := round(a+0.01);
      diff := abs((L-NN*Lambda)/L);
      if (diff >= 1e-5) then
      begin
        Writeln('>>> Incorrect definition of ', Pname);
        writeln;
        writeln('    L      ( total length ) =', L:20:12, ' [m]');
        writeln('    Lambda ( wave  length ) =', Lambda:20:12, ' [m]');
        writeln('    # of Period = L/Lambda  =', (L/Lambda):20:12, ' ?????');
        writeln;
      end;
    end;
    Lat_CheckWiggler := true;
  end;


function Lat_DealElement;

  label 9999;
  var	t, t1, t2, gap, QL, QK, QKx	: double;
	QKS, dt, Frf, Vrf		: double;
	k1, k2, harnum			: integer;
	sym1				: Lat_symbol;
	mysys				: symset;
	B				: array[-HOMmax..HOMmax] of double;
	BA				: array[-HOMmax..HOMmax] of boolean;
	DBNsavemax			: integer;
	DBNsave				: array[1..nKidMax] of DBNameType;


   PROCEDURE errorm(cmnt:str80);
     BEGIN
       Lat_errorm(cmnt, fi, fo, cc, ll, errpos, lc, chin, skipflag, line);
     END;

   procedure GetSym;
    begin
      Lat_GetSym(fi, fo, cc, ll, errpos, lc, nkw, inum, emax, emin, kmax, nmax, 
               chin, id, rnum, skipflag, rsvwd, line, 
               sym, key, ksy, sps);
    end;

   PROCEDURE test(s1:symset; cmnt: str80);
      BEGIN
          IF not (sym in s1) THEN errorm(cmnt);
      END {test};

   PROCEDURE getest(s1:symset; cmnt: str80);
      BEGIN  GetSym;
          IF not (sym in s1) THEN errorm(cmnt);
      END {test};
    
  function EVAL:double;
   begin
    Eval:=  Lat_EVAL(fi, fo, cc, ll, errpos, lc, nkw, inum, emax, emin, kmax, nmax, 
               chin, id, rnum, skipflag, rsvwd, 
               line, sym, key, ksy, sps);
   end;

 procedure ProcessBlockInput;
    Begin
       Lat_ProcessBlockInput(fi, fo, 
               cc, ll, errpos, lc, nkw, inum, emax, emin, kmax, nmax, 
               chin, id, BlockName, rnum, skipflag, rsvwd, 
               line, sym, key, ksy, sps)
    END;{ ProcessBlockInput }


  procedure CheckWiggler(i:integer);
  begin
    if not Lat_CheckWiggler(fo, i) then Goto 9999;
  end;

  procedure GetHOM;
  var i:integer;
      x, y:double;
  begin
    getest([lparent], '<(> expected');
    repeat    
      i:=round(Eval);
      if (i<1) or (HOMmax<i) then getest([], 'invalid value detected');
      test([comma], '<, > expected');
      x:=Eval;
      test([comma], '<, > expected');
      y:=Eval;
      B[i]:=x; B[-i]:=y; BA[i]:=true;
      test([comma, rparent], '<, > or <)> expected');
    until sym=rparent;
    getsym;
  end;


  procedure ClearHOMandDBN;

  var	i	: integer;

  begin
    for i:=-HOMmax to HOMmax do
    begin
      B[i] := 0.0; BA[i] := false;
    end;
    DBNsavemax := 0;
  end;


  procedure AssignHOM(elem:integer);

  var i	: integer;

  begin
    with ElemFam[elem].ElemF.M^ do
    begin
      for i:=-HOMmax to HOMmax do
      if BA[i] then
      begin
        PBpar[i] := B[i]; Porder := abs(i);
      end;
    end;
  end;


  procedure GetDBN;

  var	i	: integer;

  begin
    getest([lparent], '<(> expected:GetDBN');
    repeat
      getest([squote], '<''> expected:GetDBN');
      DBNsavemax:=DBNsavemax+1;
      for i:=1 to DBNameLen do
        DBNsave[DBNsavemax, i]:=' ';
      i := 0;
      while (chin <> '''') and (i < DBNameLen) do
      begin
        i := i + 1;
        DBNsave[DBNsavemax, i] := chin;
        Lat_Nextch(fi, fo, cc, ll, errpos, lc, chin, skipflag, line);
      end;
      getest([squote], '<''> expected:GetDBN');
      getest([comma, rparent], '<, > or <)> expected:GetDBN');
    until sym=rparent;
    getsym;
  end;


  procedure SetDBN;

  var	i	: integer;

  procedure adjdbname(var DBname1, DBname2 : DBnametype);

  var	i, j, k, offset	: integer;
	first, blank	: boolean;

  begin
    blank := true;
    for j:=1 to DBnamelen do
      if DBname1[j] <> ' ' then blank := false;
    first := true; j := 0; offset := 0;
    if not blank then
    begin
      repeat
	j := succ(j);
        DBname2[j+offset] := DBname1[j];
        if first and (DBname1[j] = ' ') then
        begin
          first := false; k := -1;
          repeat
            k := succ(k);
          until (DBname1[j+k+1] <> ' ') or (j+k+1 = DBnamelen);
          for i:=j+k+1 to 8 do
	    DBname2[i] := ' ';
	  offset := 8 - j;
        end;
      until j >= DBnamelen-offset ;
    end;
  end;

  begin
    if globval.Elem_nFam <= Elem_nFamMax then
    begin
      with ElemFam[globval.Elem_nFam] do
      begin
        NoDBN := DBNsavemax;
        if NoDBN > 0 then
          for i:=1 to NoDBN do
	    adjdbname(DBNsave[i], DBNlist[i]);
      end;
    end
    else
      writeln('Elem_nFamMax exceeded: ', globval.Elem_nFam:1, 
	      '(', Elem_nFamMax:1, ')');
  end;


begin
  Lat_DealElement:=false;

  case sym of{3.5: ****************************************}

  {****************************************************************************
       Drift 
   ****************************************************************************

 <name> : Drift, 
             L=<length>; [m]

 * Example

 L1 : Drift, L=0.30;

*** Implementation **********************************************************}

drfsym: BEGIN
          getest([comma], '<comma> expected');
          getest([lsym], '<L> expected');
          getest([eql], '<=> expected');
          rnum := Eval;
          test([semicolon], '<;> expected');
          GetSym;
          globval.Elem_nFam := globval.Elem_nFam + 1;
          if globval.Elem_nFam <= Elem_nFamMax then
          begin
            with ElemFam[globval.Elem_nFam], ElemF do 
            begin
              Pname := ElementName; PL := rnum; Pkind := drift;
              Drift_Alloc(ElemF);
            end;
          end
          else
            writeln('Elem_nFamMax exceeded: ', globval.Elem_nFam:1, 
		    '(', Elem_nFamMax:1, ')');
        end;

{****************************************************************************
     bending 
 ****************************************************************************

 <name> : Bending, 
             L      = <length>, ( [m] )
             T      = <bending angle>, ( [degree] )
             T1     = <entrance angle>, ( [degree] )
             T2     = <exit angle>, ( [degree] )
             gap    = <total magnet gap>, ( [m] )
             K      = <K-value>, ( [m^-2] ) 
	              ( K > 0 : focusing in horizontal plane )
	              ( K < 0 : defocusing in horizontal plane )
             N      = <# of kicks>, 
             method = <method>, ( 2 or 4. The method to divide Q into slices.)
                      ( The detail of <method> will be discussed later.)
             Default value is 2.
             Tilt   = <Rotation>, ( [deg], Tilt angle, designed.  )
             HOM    = (i, <Bi>, <Ai>, ( higher order component in USA notation )
                       j, <Bj>, <Aj>, ( Systematic error Only )
                       ............    ( Random errors are assigned )
                       n, <Bn>, <An>); ( in a Program File using procedures )

 * Example

 B  : bending, L=0.70, T=10.0, T1:=5.0, T2:=5.0, K=-1.0, N=8, Method=2;

*** Implementation **********************************************************}

BndSym :BEGIN{4}
        getest([comma], '<, > expected');
        getsym;
        QL:=0.0;{ L }
        QK:=0.0;{ K }
        k1:=0  ;{ N }
        t :=0.0;{ T }
        t1:=0.0;{ T1 }
        t2:=0.0;{ T2 }
        gap:=0.0;{ gap }
        k2:=Meth_Linear  ;{ method }
        dt:=0.0;
        ClearHOMandDBN;
        mysys:=[lsym, ksym, nsym, mthsym, tsym, tltsym, t1sym, t2sym, gapsym, 
		HOMsym, dbnsym];
        REPEAT
          test(mysys, 'illegal parameter');
          sym1:=sym;
          getest([eql], '<=> expected');
          case sym1 of
            lsym  : QL:=Eval;
            ksym  : QK:=Eval;
            nsym  : K1:=Round(Eval);
            tsym  : t :=Eval;
            tltsym: dt:=Eval;
            t1sym : t1:=Eval;
            t2sym : t2:=Eval;
            gapsym: gap:=Eval;
            mthsym: begin
		      k2:=Round(Eval);
                      if not (k2 in [Meth_Linear, Meth_Second, Meth_Fourth ]) then
                        getest([], 'Check integrator..');
                    end; 
            HOMsym: GetHOM;
            DBNsym: GetDBN;
          END;
          test([comma, semicolon], '<, > or <;> expected');
          IF sym=comma THEN getsym;
        until not ( sym in mysys );
        test([semicolon], '<;> expected.');
        GetSym;
        globval.Elem_nFam := globval.Elem_nFam + 1;
        if globval.Elem_nFam <= Elem_nFamMax then
        begin
          with ElemFam[globval.Elem_nFam], ElemF do 
          begin
            Pname := ElementName; PL := QL; Pkind := mpole;
            Mpole_Alloc(ElemF);
            with M^ do
            begin
              Pmethod:=k2; PN:=K1;
	      if PL <> 0.0 then
		Pirho:=T*pi/180/PL
	      else
		Pirho:=T*pi/180;
	      ptx1 := t1; ptx2 := t2;
	      Pgap:=gap;
              AssignHOM(globval.Elem_nFam);
              SetDBN;
              PBpar[2]:=QK;
              PdTpar:=dt;
            end;
          end;
        end
        else
          writeln('Elem_nFamMax exceeded: ', globval.Elem_nFam:1, 
                  '(', Elem_nFamMax:1, ')');
       end;

{****************************************************************************
     Quadrupole 
 ****************************************************************************

 <name> : Quadrupole, 
             L=<length>, ( [m] )
             K =<K-value>, ( [m-2] )
             N =<# of kicks>, 
             Method=<method>, 
             Tilt=<Rotation>, ( [deg], Tilt angle, designed.  )
             HOM=(i, <Bi>, <Ai>, ( higher order component in USA notation )
                  j, <Bj>, <Aj>, ( Systematic error Only )
                  ............    ( Random errors are assigned )
                  n, <Bn>, <An>); ( in a Program File using procedures )

 * Example

 Q1 : Quadrupole, L=0.5, K=2.213455, N=4, Method=4;

*** Implementation **********************************************************}

QdSym: BEGIN{4}
         getest([comma], '<, > expected');
         getsym;
         QL:=0.0;{ L }
         QK:=0.0;{ K }
         k1:=0  ;{ N }
         k2:=Meth_Linear  ;{ method }
         dt:=0.0;
         ClearHOMandDBN;
         mysys:=[lsym, ksym, nsym, mthsym, tltsym, HOMsym, DBNsym];
         REPEAT {5: read L, K, N, T, T1, T2 }
           test(mysys, 'illegal parameter');
           sym1:=sym;
           getest([eql], '<=> expected');
           case sym1 of{6}
             lsym:   QL:=Eval;
             ksym:   QK:=Eval;
             nsym:   K1:=Round(Eval);
             mthsym: begin
		       k2:=Round(Eval);
                       if not (k2 in [Meth_Linear, Meth_First, Meth_Second,
			 Meth_Fourth ]) then
			 getest([], 'Check integrator..');
                     end;
             dtsym:  dT:=Eval;
             homsym: GetHOM;
             DBNsym: GetDBN;
           END;
           test([comma, semicolon], '<, > or <;> expected');
           IF sym=comma THEN getsym;
         until not ( sym in mysys);{5}
         test([semicolon], '<;> expected.');
         GetSym;
         globval.Elem_nFam := globval.Elem_nFam + 1;
         if globval.Elem_nFam <= Elem_nFamMax then
         begin
           with ElemFam[globval.Elem_nFam], ElemF do 
           begin
             Pname := ElementName; PL := QL; Pkind := Mpole;
             Mpole_Alloc(ElemF);
             with M^ do
             begin
               Pmethod := k2; PN := K1;
               PdTpar:=dT;
               SetDBN;
               PBpar[2]:=QK;
             end;
           END;
         end
         else
           writeln('Elem_nFamMax exceeded: ', globval.Elem_nFam:1, 
                   '(', Elem_nFamMax:1, ')');
       END;

{****************************************************************************
     Sextupole 
 ****************************************************************************

 <name> : Sextupole, 
             L=<length>, ( [m] )
             K =<K-value>, ( [m-3] )
             Tilt=<Rotation>, ( [degree], Tilt angle, designed.  )
             HOM=(i, <Bi>, <Ai>, ( higher order component in USA notation )
                  j, <Bj>, <Aj>, ( Systematic error Only )
                  ............    ( Random errors are assigned )
                  n, <Bn>, <An>); ( in a Program File using procedures )

 * Example

 SF : Sextupole, K=-10.236345;

*** Implementation **********************************************************}

SexSym: BEGIN{4}
          QL := 0.0;{ L }
          QK:=0.0;{ K }
          k1:=0  ;{ N }
          k2:=Meth_Linear  ;{ method }
          dt:=0.0;
          ClearHOMandDBN;
          getest([comma, semicolon], '<, > or <;> expected');
          if sym=comma then 
          begin
            getsym;
            mysys:=[lsym, ksym, nsym, mthsym, tltsym, HOMsym, DBNsym];
            REPEAT {5: read L, K, N, T, T1, T2 }
              test(mysys, 'illegal parameter');
              sym1:=sym;
              getest([eql], '<=> expected');
              case sym1 of{6}
                lsym: QL := Eval;
                ksym: QK:=Eval;
                nsym: K1:=Round(Eval);
                mthsym: begin
			  k2:=Round(Eval);
                          if not (k2 in [Meth_Linear, Meth_Second, Meth_Fourth ]) then
			    getest([], 'Check integrator..');
                        end; 
                tltsym: dT:=Eval;
                homsym: GetHOM;
                DBNsym: GetDBN;
              END;
              test([comma, semicolon], '<, > or <;> expected');
              IF sym=comma THEN getsym;
            until not ( sym in mysys);{5}
            test([semicolon], '<;> expected.');
          end;
          GetSym;
          globval.Elem_nFam := globval.Elem_nFam + 1;
          if globval.Elem_nFam <= Elem_nFamMax then
          begin
            with ElemFam[globval.Elem_nFam], ElemF do 
            begin
              Pname := ElementName; PL := QL; Pkind := Mpole;
              Mpole_Alloc(ElemF);
              with M^ do
              begin
                Pmethod := k2; PN := K1;
                if PL <> 0.0 then
	          Pthick := thick
	        else
	          Pthick := thin;
                PdTpar:=dT;
                AssignHOM(globval.Elem_nFam);
                SetDBN;
                PBpar[3]:=QK;
              end;
            end;
          end
          else
            writeln('Elem_nFamMax exceeded: ', globval.Elem_nFam:1, 
                    '(', Elem_nFamMax:1, ')');
        end;

{****************************************************************************
     Cavity
 ****************************************************************************

 <name> : Cavity, 
             Frequency = <Frf>, ( [Hz] )
             Voltage   = <Vrf>, ( [V]  )
	     harnum    = <h>

 * Example

 CAV : Cavity, Frequency = 499.95e6, Voltage=1.22e6, harnum=328;

*** Implementation **********************************************************}

CavSym: begin
          ClearHOMandDBN;
          getest([comma], '<, > expected');
          getsym;
          Frf := 0.0;{ Frf }
          Vrf := 0.0;{ Vrf }
          harnum := 0;{ Voff }
          mysys := [Frqsym, Vrfsym, harnumsym, DBNsym];
          repeat 
            test(mysys, 'illegal parameter');
            sym1 := sym;
            getest([eql], '<=> expected');
            case sym1 of
              Frqsym:    Frf := Eval;
              Vrfsym:    Vrf := Eval;
              harnumsym: harnum := round(Eval);
              DBNsym:GetDBN;
            end;
            test([comma, semicolon], '<, > or <;> expected');
            IF sym = comma THEN getsym;
          until not ( sym in mysys);
          test([semicolon], '<;> expected.');
          GetSym;
          globval.Elem_nFam := globval.Elem_nFam + 1;
          if globval.Elem_nFam <= Elem_nFamMax then
          begin
            with ElemFam[globval.Elem_nFam], ElemF do 
            begin
              Pname := ElementName; Pkind := Cavity;
              Cav_Alloc(ElemF);
              with C^ do
              begin
                Pvolt := Vrf; { Voltage [V] }
                Pfreq := Frf; { Frequency in Hz }
		Ph := harnum;
                SetDBN;
              end;
            end;
          end
          else
            writeln('Elem_nFamMax exceeded: ', globval.Elem_nFam:1, 
                    '(', Elem_nFamMax:1, ')');
        end;


{****************************************************************************
     Orbit Corrector 
 ****************************************************************************

  <name> : Corrector, <direction>, L=<length>, ;

  <name> :== Alphanumeric string. Up to NameLength character length.
             BEGIN with an alphabet.
  <direction> :== 'horizontal'|'vertical'

 * Example

 COH : Corrector, horizontal;

*** Implementation **********************************************************}


CorSym: BEGIN{4}
          QL := 0.0;{ L }
          QK := 0.0;{ K }
          k1 := 0;{ N }
          k2 := Meth_Linear  ;{ method }
          dt := 0.0;
          ClearHOMandDBN;
          getest([comma], '<, > expected');
          if sym=comma then 
          begin
            getsym;
            mysys:=[lsym, nsym, mthsym, horsym, versym, tltsym, DBNsym];
            REPEAT {5: read L, K, N, T, T1, T2 }
              test(mysys, 'illegal parameter');
              sym1:=sym;
              if sym in [lsym, nsym, mthsym, tltsym, DBNsym] then
              begin
                getest([eql], '<=> expected');
                if sym=eql then
                case sym1 of{6}
                  lsym:   QL:=Eval;
                  nsym:   K1:=Round(Eval);
                  mthsym: begin
			    k2:=Round(Eval);
                            if not (k2 in [Meth_Linear, Meth_Second, Meth_Fourth ]) then
			      getest([], 'Check integrator..');
                          end; 
                  tltsym: dT:=Eval;
                  DBNsym: GetDBN;
                END;
              end
	      else 
              begin
                if sym1=horsym then
		  qk:= 1 
                else
		  qk:=-1;
                Getsym;
              end;
              test([comma, semicolon], '<, > or <;> expected');
              IF sym=comma THEN getsym;
            until not ( sym in mysys);{5}
            test([semicolon], '<;> expected.');
          end;
          GetSym;
          globval.Elem_nFam := globval.Elem_nFam + 1;
          if globval.Elem_nFam <= Elem_nFamMax then
          begin
            with ElemFam[globval.Elem_nFam], ElemF do 
            begin
              Pname := ElementName; PL := QL; Pkind := mpole;
              Mpole_Alloc(ElemF);
              with M^ do
              begin
                SetDBN;
                if PL <> 0.0 then
		  Pthick := thick
	        else
		  Pthick := thin;
                Pmethod := k2; PN := K1;
                PdTpar := dT;
              end;
            end;
          end
          else
            writeln('Elem_nFamMax exceeded: ', globval.Elem_nFam:1, 
                    '(', Elem_nFamMax:1, ')');
        end;

{****************************************************************************
     Beam Position Monitor
 ****************************************************************************

 <name> : Beam Position Monitor;

  <name>:== Alphanumeric string. Up to NameLength character length.
            BEGIN with an alphabet.

 * Example

 BPM : Beam Position Monitor;

*** Implementation **********************************************************}

BemSym: BEGIN
          ClearHOMandDBN;
          getest([possym], '<position> expected');
          getest([monsym], '<monitor> expected');
          getest([comma, semicolon], '<, > or <;> expected');
          if sym = comma then
          begin
            getest([DBNsym], 'illegal parameter'); sym1 := sym;
            getest([eql], '<=> expected');
            if sym1 = DBNsym then GetDBN;
            test([semicolon], '<;> expected');
          end;
          GetSym;
          globval.Elem_nFam := globval.Elem_nFam + 1;
          if globval.Elem_nFam <= Elem_nFamMax then
          begin
            with ElemFam[globval.Elem_nFam], ElemF do 
            BEGIN
              Pname := ElementName; Pkind := Mpole;
              Mpole_Alloc(ElemF);
              with M^ do
              begin
                Pthick:=thin; SetDBN;
              end;
            END;
          end
          else
            writeln('Elem_nFamMax exceeded: ', globval.Elem_nFam:1, 
                    '(', Elem_nFamMax:1, ')');
        END;


{****************************************************************************
     Marker
 ****************************************************************************

 <name> : Marker;

  <name>:== Alphanumeric string. Up to NameLength character length.
            BEGIN with an alphabet.

 * Example

 SYM : Marker;

*** Implementation **********************************************************}

MrkSym: begin
          ClearHOMandDBN;
          getest([comma, semicolon], '<, > or <;> expected');
          if sym = comma then
          begin
            getest([DBNsym], 'illegal parameter'); sym1 := sym;
            getest([eql], '<=> expected');
            if sym1 = DBNsym then GetDBN;
            test([semicolon], '<;> expected');
          end;
          GetSym;
          globval.Elem_nFam := globval.Elem_nFam + 1;
          if globval.Elem_nFam <= Elem_nFamMax then
          begin
            with ElemFam[globval.Elem_nFam], ElemF do 
            begin
              Pname := ElementName; PL := 0.0; Pkind := marker;
              SetDBN;
            end;
          end
          else
            writeln('Elem_nFamMax exceeded: ', globval.Elem_nFam:1, 
                    '(', Elem_nFamMax:1, ')');
        end;


{****************************************************************************
     Ghost
 ****************************************************************************

 <name> : Ghost;

  <name>:== Alphanumeric string. Up to NameLength character length.
            BEGIN with an alphabet.

 * Example

OBAKE : Ghost;

*** Implementation **********************************************************}
{----------->>>
GstSym:BEGIN
        getest([comma], '<, > expexted');
        getest([typsym], '<type> expected');
        getest([eql], '<=> expected');
        QL:=Eval;
        test([semicolon], '<;> expected');
        getsym;
        if sym=DBNsym then GetDBN;
        globval.Elem_nFam := globval.Elem_nFam + 1;
        if globval.Elem_nFam <= Elem_nFamMax then
        begin
          with ElemFam[globval.Elem_nFam].ElemF do 
            with ElementT[globval.Elem_nFam] do 
            BEGIN
              Pname:=ElementName; Pkind:=Ghost; PN:=round(QL);
              SetDBN;
            END;
        end
        else
          writeln('Elem_nFamMax exceeded: ', globval.Elem_nFam:1, 
                  '(', Elem_nFamMax:1, ')');
       END;
<<-----------------------------}


{****************************************************************************
     Multipole 
 ****************************************************************************

 <name> : Multipole, 
             L=<length>, ( [m] )
             T =<bending angle>, ( [degree] )
             T1=<entrance angle>, ( [degree] )
             T2=<exit angle>, ( [degree] )
             Tilt=<Rotation>, ( [deg], Tilt angle, designed.  )
             N =<# of kicks>, 
             method=<method>, ( 2 or 4. The method to divide Q into slices.)
                     ( The detail of <method> will be discussed later.)
             Default value is 2.
             HOM=(i, <Bi>, <Ai>, ( higher order component in USA notation )
                  j, <Bj>, <Aj>, ( Systematic error Only )
                  ............    ( Random errors are assigned )
                  n, <Bn>, <An>); ( in a Program File using procedures )

 * Example

 B  : multipole, L=0.70, T=10.0, T1=5.0, T2=5.0, 
                 HOM=(2, -1.0, 0), N=8, Method=2;


 QF  : multipole, L=0.70, 
                 HOM=(2, 2.50, 0.0, 
                      4, 1.01e7, 0.0), 
                  N=8, Method=2;

*** Implementation **********************************************************}

MPsym:BEGIN{4}
        getest([comma], '<, > expected');
        getsym;
        QL:=0.0;{ L }
        QK:=0.0;{ K }
        k1:=0  ;{ N }
        t :=0.0;{ T }
        t1:=0.0;{ T1 }
        t2:=0.0;{ T2 }
        k2:=Meth_Linear  ;{ method }
        dt:=0.0;
        ClearHOMandDBN;
        mysys:=[lsym, nsym, mthsym, tsym, t1sym, t2sym, 
                tltsym, HOMsym, DBNsym];
        REPEAT { read L, K, N }
          test(mysys, 'illegal parameter');
          sym1:=sym;
          getest([eql], '<=> expected');
          case sym1 of
            lsym :  QL:=Eval;
            nsym :  K1:=Round(Eval);
            tsym :  t :=Eval;
            tltsym: dt :=Eval;
            t1sym:  t1:=Eval;
            t2sym:  t2:=Eval;
            mthsym: begin
		      k2:=Round(Eval);
                      if not (k2 in [Meth_Linear, Meth_Second, Meth_Fourth ]) then
			getest([], 'Check integrator..');
                    end; 
            HOMsym: GetHOM;
            DBNsym: GetDBN;
          END;
          test([comma, semicolon], '<, > or <;> expected');
          IF sym=comma THEN getsym;
        until not (sym in mysys) ;
        test([semicolon], '<;> expected.');
        GetSym;
        globval.Elem_nFam := globval.Elem_nFam + 1;
        if globval.Elem_nFam <= Elem_nFamMax then
        begin
          with ElemFam[globval.Elem_nFam], ElemF do 
          begin
            Mpole_Alloc(ElemF);
            with M^ do
            begin
              Pname:=ElementName;
              Pkind:=mpole;
              PL:=QL;
              if PL <> 0d0 then
	      begin
		Pthick := thick; Pirho:=T*pi/180/PL
	      end
	      else
	      begin
		Pthick := thin; Pirho:=T*pi/180;
	      end;
              PN:=K1; Pmethod:=k2;
	      ptx1 := t1; ptx2 := t2;
              PdTpar:=dT;
              AssignHOM(globval.Elem_nFam);
              SetDBN;
            end;                    
          END;
        end
        else
          writeln('Elem_nFamMax exceeded: ', globval.Elem_nFam:1, 
                  '(', Elem_nFamMax:1, ')');
      END;

{****************************************************************************
     Wiggler 
 ****************************************************************************

 <name> : Wiggler, 
             L      = <length [m]>,
             BoBrho = <B/Brho [m-1]>,
             Lambda = <period [m]>,
	     kx	    = <[m]>,
             N      = <no of integration steps>,
             Method = <method>, 

 * Example

 U143: wiggler, L=4.80, K=0.5, Lambda=0.15, N=20, Method=0;

*** Implementation **********************************************************}

wglsym: BEGIN
	  getest([comma], '<, > expected');
          getsym;
	  QL := 0d0; QK := 0d0; QKx := 0d0; QKS := 0d0;
	  k1 := 0; k2 := Meth_Linear;
          dt := 0d0;
	  ClearHOMandDBN;
          mysys:=[lsym, lmdsym, bobrhosym, kxsym, nsym, mthsym, tltsym, DBNsym];
          REPEAT
            test(mysys, 'illegal parameter');
            sym1:=sym;
            getest([eql], '<=> expected');
            case sym1 of{6}
	      lsym      : QL := Eval;
              bobrhosym : QK := Eval;
	      kxsym     : Qkx := Eval;
              nsym      : k1 := round(Eval);
              mthsym    : begin
		            k2 := Round(Eval);
                            if not (k2 in [Meth_Linear, Meth_First, Meth_Second,
			       Meth_Fourth, Meth_genfun]) then
			       getest([], 'Check integrator..');
                          end;
              lmdsym    : QKS := Eval;
              dtsym     : dT := Eval;
	      DBNsym    : GetDBN;
	    end;
            test([comma, semicolon], '<, > or <;> expected');
            IF sym=comma THEN getsym;
          until not ( sym in mysys);{5}
	  test([semicolon], '<;> expected');
	  getsym;
          globval.Elem_nFam := globval.Elem_nFam + 1;
          if globval.Elem_nFam <= Elem_nFamMax then
          begin
	    with ElemFam[globval.Elem_nFam], ElemF do 
	    begin
	      Pname := ElementName; PL := QL; Pkind := wigl;
	      Wiggler_Alloc(ElemF);
	      with W^ do
	      begin
                Pmethod := k2; PN := K1;
                PdTpar := dT;
                SetDBN;
		PBoBrho := QK;
		Pkx := QKx;
		{ Equivalent vertical gradient }
                PBW[2] := sqr(QK)/2d0;
	        Plperiod := QKS;
	      end;
	    end;
	    CheckWiggler(globval.Elem_nFam);
          end
          else
            writeln('Elem_nFamMax exceeded: ', globval.Elem_nFam:1, 
                    '(', Elem_nFamMax:1, ')');
	end;

{***************************************************************************
        BLOCK DEFINITION        
 ***************************************************************************}

    ident, intcon, invsym: { Block Definition }
			   ProcessBlockInput;
  END;{3.5:of CASE}
  Lat_DealElement:=true;
9999:
end;{of procedure Lat_DealElement}


   PROCEDURE errorm(cmnt:str80);

   BEGIN
      IF errpos = 0 THEN
      BEGIN  {write(fo, ' ****')} END;
      IF cc > errpos THEN
          BEGIN
             write(fo, ' ':cc-errpos, '^', cmnt);
             errpos := cc + 3;
          END;
      while not eof(fi) do 
         Lat_Nextch(fi, fo, cc, ll, errpos, lc, chin, skipflag, line);
      ErrFlag:=true;
      goto 9999;      
    END  {error};

  PROCEDURE GetSym;     { reads next symbol  }
   BEGIN

    Lat_GetSym(fi, fo, cc, ll, errpos, lc, nkw, inum, emax, emin, kmax, nmax, 
               chin, id, rnum, skipflag, rsvwd, line, 
               sym, key, ksy, sps);
   END  {GetSym};

   PROCEDURE test(s1:symset; cmnt: str80);
      BEGIN
          IF not (sym in s1) THEN errorm(cmnt);
      END {test};

   PROCEDURE getest(s1:symset; cmnt: str80);
      BEGIN  GetSym;
          IF not (sym in s1) THEN errorm(cmnt);
      END {test};


  PROCEDURE init_reserved_words;

    PROCEDURE Reg(name:alfa; ks:Lat_symbol);

    BEGIN
      nkw:=nkw+1;  key[nkw]:=name; ksy[nkw]:=ks
    END;

  BEGIN
    nkw:=0;
    Reg('and            ', andsym); Reg('b0             ', b0sym);
    Reg('beam           ', bemsym); Reg('bending        ', bndsym);
    Reg('cavity         ', cavsym); Reg('cell           ', celsym);     
    Reg('chromaticity   ', chmsym); Reg('corrector      ', corsym);
    Reg('dbname         ', dbnsym); Reg('define         ', defsym);   
    Reg('dispersion     ', dspsym); Reg('drift          ', drfsym);
    Reg('dt             ', dtsym);  Reg('end            ', endsym);    
    Reg('execute        ', exesym); Reg('exit           ', extsym);
    Reg('focusing       ', fcssym);
    Reg('frequency      ', frqsym); Reg('fringe         ', frgsym);
    Reg('galilean       ', xytsym); Reg('ghost          ', gstsym);
    Reg('hom            ', homsym); Reg('horizontal     ', horsym);
    Reg('inv            ', invsym); Reg('kicker         ', kicksym);
    Reg('ks             ', kssym);  Reg('lambda         ', lmdsym);
    Reg('lattice        ', latsym); Reg('marker         ', mrksym);
    Reg('matrix         ', matsym); Reg('method         ', mthsym);
    Reg('monitor        ', monsym); Reg('multipole      ', mpsym);
    Reg('nonlinear      ', nbdsym); Reg('parameter      ', prmsym);
    Reg('position       ', possym); Reg('print          ', prnsym);	
    Reg('quadrupole     ', qdsym);  Reg('sextupole      ', sexsym);	
    Reg('symmetry       ', symsym); Reg('t1             ', t1sym);     
    Reg('t2             ', t2sym);  Reg('gap            ', gapsym);
    Reg('table          ', tblsym);     
    Reg('task           ', tsksym); Reg('tilt           ', tltsym);
    Reg('type           ', typsym);
    Reg('use            ', usesym); Reg('vertical       ', versym);
    Reg('harnum         ', harnumsym); Reg('voltage        ', vrfsym);
    Reg('wiggler        ', wglsym);

    sps['+'] := plus;      sps['-'] := minus;
    sps['('] := lparent;   sps[')'] := rparent;
    sps['='] := eql;       sps[','] := comma;
    sps['[']:=lbrack;      sps[']']:=rbrack;
    sps[''''] := squote;   sps['&'] := andsy;
    sps[';'] := semicolon; sps['/'] := rdiv;
    sps[':'] := colon;

    lc := 0; { reset line counter } 
    ll := 0; { reset line length  }
    cc := 0; { reset char counter }
    errpos := 0; { reset error position }
    chin := ' ';   { reset current char   }
    skipflag:=false; { reset skip flag  }
    defbegsys:=[ident];
    elmbegsys:=[qdsym, sexsym, corsym, bemsym, gstsym, 
                mrksym, nbdsym, frgsym, xytsym, 
                drfsym, bndsym, wglsym, MPSYM, cavsym];
  END;


{*****************************************************
 *                                                   *
 *                  P A R S E R                      *
 *                                                   *
 *****************************************************}
  
  PROCEDURE DealWithDefns;

  var	idsave, ElementName, BlockName, IdentName	: partsname;
	i, j, k, k1					: integer;
    
  function EVAL:double;

  begin
    Eval:=  Lat_EVAL(fi, fo, cc, ll, errpos, lc, nkw, inum, emax, emin, kmax, nmax, 
            chin, id, rnum, skipflag, rsvwd, line, sym, key, ksy, sps);
  end;

{************* DEAL WITH DEFINITIONS ********************************}

BEGIN{0}
  GetSym;
  IF sym=latsym THEN
  {***** The first word must be 'lattice' *********}
  BEGIN{1}
    getest([semicolon], '<;> expected');
    getest([ident], '<identIFiler> expected');

 {**************************************************************************}
  IF sym=ident THEN
  REPEAT{2}
    IF sym=ident THEN 
    BEGIN{2.5:-----------------------}
      idsave:=id;
      ElementName:=id;
      BlockName:=id;
      Getest([colon, eql, period], '<colon>, <=> or <.> expected');
      IF sym=colon THEN

{*****************************
 *                           *
 *     DEAL WITH ELEMENT     *
 *                           *
 ***************************** }    

 BEGIN{3}
   getest(elmbegsys+[ident, intcon, invsym], 
          '<identifier>, <integer> or <INV> expected');
   IF sym in (elmbegsys+[ident, intcon, invsym]) THEN
      if not Lat_DealElement(fi, fo, 
               cc, ll, errpos, lc, nkw, inum, emax, emin, kmax, nmax, 
               chin, id, ElementName, BlockName, 
               rnum, skipflag, rsvwd, line, sym, key, ksy, sps)   
      then goto 9999;

END{3:of IF sym=colon} ELSE 

{*************************************************************
 *                                                           *
 *         P A R A M E T E R  A S S I G N M E N T            *
 *                                                           *
 *************************************************************}

IF sym=eql THEN
BEGIN{3:of parameter}
  IdentName:=idsave;
  i := CheckUDItable(IdentName);
  IF i = 0 THEN
    EnterUDItable(IdentName, Eval)
  ELSE
    ModUDItable(i, Eval);
  test([semicolon], '<;> expected');
  GetSym;
END{3:of parameter }
else {8888888888}
IF sym=period THEN
BEGIN{3:of parameter}
{-----
  IdentName:=idsave;
  i:=CheckElementtable(IdentName);
  IF i=0 THEN Test([], '<element name> expected');
  getest([lsym, tsym, t1sym, t2sym, gapsym, ksym], 'illegal component');
  sym1:=sym;
  getest([eql], '<=> expected');
  case sym1 of
    lsym:  ElemFam[i].ElemF.PL :=Eval; 
    ksym:  ElemFam[i].ElemF.Pk :=Eval; 
    tsym:  ElemFam[i].ElemF.Pt :=Eval; 
    t1sym: ElemFam[i].ElemF.Pt1:=Eval; 
    t2sym: ElemFam[i].ElemF.Pt2:=Eval;
    gapsym: ElemFam[i].ElemF.Pgap:=Eval;
  END;
  test([semicolon], '<;> expected');
  GetSym;
  -----}
END;{3:of parameter }

END;{2.5}

{*************************************************************
 *                                                           *
 *         C E L L    D E F I N I T I O N                    *
 *                                                           *
 *************************************************************
  
 CELL : <block name>, SYMMETRY=<symmetry>;

    <block name>:== name of a block.
    <symmetry>:== number of supersymmetry:== number of the block/ring

 * Example
  
   CELL : BL1, Symmetry=12;

*** Implementation *******************************************}

IF sym = celsym THEN
BEGIN{3}
  getest([colon], '<colon> expected');
  getest([ident], '<Block name> expected');
  i:=CheckBlockstable(id);
  IF i=0 THEN test([], '<Block name> expected');
  k:=0;
  IF i<>0 THEN
  BEGIN{4}
    with Blocks[i] do 
    BEGIN{5}
      FOR j:=Bstart to Bowari do 
      BEGIN{6} 
        k:=k+1;
        k1:=Bstack[j];
        if k <= Cell_nLocmax then
          CELL[k].Fnum := k1
	else
          writeln('** Cell_nLocMax exhausted: ', k:1, '(', Cell_nLocMax:1, ')'); 
      END{6};
    END{5};
  END{4};
  globval.Cell_NLoc:=k;{number of Elements in a cell}
  getest([comma], '<, > expected');
  getest([symsym], '<symmetry> expented');
  getest([eql], '<=> expected');
  symmetry:=Round(Eval);
  if symmetry>=1 then ring:=true
     else begin symmetry:=1; ring:=false end;
  test([semicolon], '<;> expected');
  GetSym;
END;{3: of celsym}


case{3} sym of

{*****************************************

     PRINT element-name
     PRINT block_name
     PRINT parameter

 *****************************************}    

 prnsym:
   BEGIN{4}
    getest([ident], '<identifiler> expected');
    IdentName:=id;
    i:=CheckElementTable(Identname);
    if i<>0 then {PrintElementParam(i)} else
    begin{5}
      i:=CheckBlocksTable(Identname);
      if i<>0 then {PrintBlockParam(i)} else
      begin{6}
        i := CheckUDItable(IdentName);
        if i <> 0 then
	  {PrintUDIParam(i)}
	else
          getest([], ' invalid expression');
      end;{6}
    end;{5}
    if i<>0 then
    begin
      GeTest([semicolon], '<;> expected');
      GetSym;
    end;
   END;{4}

   otherwise;

 END;{3:of CASE}

UNTIL{2} not (sym in [ident, chmsym, dspsym, celsym, prnsym]);

test([endsym], '<END> expected');
getest([semicolon], '<;> expexted');

END{1} ELSE 
BEGIN{1}
        test([], '<illegal operand> detected');
END;{1}
END;{0}


  PROCEDURE GetEnergy;

  var	k			: integer;
        cgam, cu, xmc2, hbc, con	: double;

  begin
     k := CheckUDItable('energy         ');
     if k = 0 then 
     begin
       writeln('> Beam energy is not defined.');
       write('  Input beam energy in [GeV] := '); readln(globval.energy);
       EnterUDItable('energy         ', globval.energy);
     end
     else
     begin
       RefUDItable('energy         ', globval.energy);
     end;

     CGAM := 8.846056192D-05; CU := 55.0d0/24d0/SQRT(3d0);
     XMC2 := 0.5110034D-03; HBC := 1.9732858D-16;
     CON := 3d0*CU*CGAM*HBC/4.0d0/PI/pwr(XMC2, 3d0);
     CRAD := CGAM*pwr(globval.energy, 3d0)/2d0/PI;
     CFLUC := CON*pwr(globval.energy, 5d0);
   end;


  PROCEDURE GetDP;

  var k:integer;

  begin
    k := CheckUDItable('dp             ');
        if k=0 then 
          begin
            writeln('> dP/P is not defined.');
            write('  Input dP/P := ');readln(globval.dpcommon);
            EnterUDItable('dp             ', globval.dpcommon);
          end else
          begin
            RefUDItable('dp             ', globval.dpcommon);
          end;
  end;

PROCEDURE GetCODEPS;
  var k:integer;
  begin k:=CheckUDItable('codeps         ');
        if k=0 then 
          begin
            writeln('> CODEPS is not defined.');
            write('  Input CODEPS := ');readln(globval.codeps);
            EnterUDItable('codeps         ', globval.energy);
          end else
          begin
            RefUDItable('codeps         ', globval.codeps);
          end;
  end;

  function Circumference:double;
  var i:integer;
      S:double;
  begin
    S:=0.0;
    for i:=1 to globval.Cell_NLoc do
      S:=S+ElemFam[Cell[i].Fnum].ElemF.PL;
    Circumference:=S;
  end;

  procedure RegisterKids;

  var	i	: integer;

  begin
    if globval.Elem_nFam <= Elem_nFamMax then
    begin
      for i:=1 to globval.Elem_nFam do
        ElemFam[i].nKid := 0;
    end
    else
      writeln('Elem_nFamMax exceeded: ', globval.Elem_nFam:1, 
              '(', Elem_nFamMax:1, ')');
  
    for i:=1 to globval.Cell_nLoc do
    begin
      with ElemFam[Cell[i].Fnum] do
      begin
        nKid := nKid + 1;
	if nKid <= nKidMax then
	begin
          KidList[nKid] := i;
          Cell[i].Knum := nKid;
	end
	else
          writeln('nKidMax exceeded: ', nKid:1, '(', nKidMax:1, ')');
      end;
    end;
  end;


  procedure PrintResult;

  var	j, nKid	: integer;

  begin
    writeln;
    writeln('  LAT Statistics:');
    writeln('  Number of constants: UDImax               =',
	    udic:5,                 ', UDImax          =', UDImax:5);
    writeln('  Number of Families : globval.Elem_nFam    =',
	    globval.Elem_nFam:5,    ', Elem_nFamMax    =', Elem_nFamMax:5);
    nKid := 0;
    for j:=1 to globval.Elem_nFam do
      if ElemFam[j].nKid > nKid then
	nKid := ElemFam[j].nKid;
    writeln('  Max number of Kids : nKidMax              =',
	    nKid:5,                 ', nKidMax         =', nKidMax:5);
    writeln('  Number of Blocks   : NoB                  =',
	    NoB:5,                  ', NoBmax          =', NoBmax:5);
    writeln('  Number of Elements : globval.Cell_nLoc    =',
	    globval.Cell_nLoc:5,    ', Cell_nLocmax    =', Cell_nLocmax:5);
    writeln('  Circumference      : ',Circumference:12:7,' [m]');
    writeln;
  end;


BEGIN
  udic:=0;
  globval.Cell_nLoc:=0; globval.Elem_nFam:=0;
  NoB:=0;
  Symmetry:=0;

  globval.CODeps:= 0.0; globval.dPcommon := 0.0; globval.Energy := 0.0;

  ErrFlag:=false;
  init_reserved_words;
  GetSym;
  IF sym=defsym THEN  DealWithDefns;

  IF symmetry <> 0 then
  begin
    GetEnergy; GetCODeps; GetdP;
  end;

  Close(fi); Close(fo);  
  RegisterKids;
  PrintResult;
9999:Lattice_Read:=not ErrFlag;
END;{ of ReadLattice }

  end.
