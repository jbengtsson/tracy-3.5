program tracy2(input, output);

{ J. Bengtsson, Lawrence Berkeley Lab. }

{ PSI version, modified for Swiss quality graphics:  graphx }

%include 'mathlib.def'
%include 'dab.def'
%include 'eigenv.def'
%include 'stringlib.def'

%include 'pascommon.def'
%include 't2common.def'

%include 'graphx.def'

%include 't2elem.def'

%include 't2lat.def'

%include 't2cell.def'
%include 't2ring.def'
%include 't2bump.def'

%include 'fft.def'

%include 'pascalio.def'

%include 'amoeba.def'

%include 'svd.def'
%include 'dsvdc.def'
%include 'svbksb.def'

procedure  Pascals;

{ Pascal-s, by N. Wirth, ETH, Zurich }

label 9999;

const
   nkwmax = 100;        {max no. of key words}
   llng = 132;          {inputline lenght}
   emax = 322;
   emin = -292;
   kmax = 15;           {max no. of significant digits}
   tmax = 5000;         {size of table}
   bmax = 1000;         {size of block-table}
   amax = 1000;         {size of array table}
   c2max= 1000;         {size of real constant table}
   csmax= 30;           {max no. of cases}
   cmax = 15000;        {size of code}
   lmax = 10;           {maximum level}
   smax = 5000;         {size of string table}
   ermax= 58;
   omax = 300;		{highest order code}
   xmax = maxint;       {2**16 - 1, array index}
   nmax = maxint;	{2**16 - 1, biggest integer constant}

   lineleng  = 132;     {output line lenght}
   linelimit = 400;     {max lines to print}
   stacksize = 1500000;  {stacksize}

   tabsym = 9;

   falseconst = 0;	{implementation dependent, quel bordel}
   trueconst = 1; { VAX }
{   trueconst = 16777216; IBM PC, SUN, RSX 6000 }
 
type

   symbol =
     (intcon,     realcon,   charcon,   stringcon,   notsy, 
      plus,       minus,     times,     idiv,        rdiv, 
      imod,       andsy,     orsy,      eql,         neq, 
      geq,        gtr,       lss,       leq,         lparent, 
      rparent,    lbrack,    rbrack,    comma,       semicolon, 
      period,     colon,     becomes,   constsy,     typesy, 
      varsy,      funcsy,    procsy,    arraysy,     recordsy, 
      programsy,  ident,     beginsy,   ifsy,        casesy,
      repeatsy,   whilesy,   forsy,     endsy,       elsesy,
      untilsy,    ofsy,      dosy,      tosy,        downtosy, 
      thensy,     incsy  );

   index	= - xmax .. +  xmax;
   object	= (konstant, vvariable, type1, prozedure, funktion);
   types	= (notyp, ints, reals, bools, chars, arrays, records, texts);
   symset	= set of symbol;
   typset	= set of types;

   item		= record
		    typ		: types;
		    ref1	: index;
		  end;

   order	= packed record
		    f	: - omax .. + omax;
		    x	: - lmax .. + lmax;
		    y	: - nmax .. + nmax;
		  end;

var

   nkw    : integer;                { number of key words }
   i	  : integer;
   sy     : symbol;                 {last symbol read by insymbol}
   errs   : set of 0..ermax;
   id     : alfa;                   {identifier freom insymbol}
   inum   : integer;                {integer from insymbol}
   rnum   : double;		    {real number from insymbol}
   sleng  : integer;                {string length}
   chin			: char;	    {last character read from source program}
   cc			: integer;  {character counter}
   lc			: integer;  {program location counter}
   linecounter		: integer;
   ll			: integer;  {length of current line}
   errpos, incllev	: integer;
   progname		: alfa;
   constbegsys, typebegsys, blockbegsys,
   facbegsys, statbegsys			: symset;

   line			: packed array [1..llng] of char;
   key			: array [1..nkwmax] of alfa;
   ksy			: array [1..nkwmax] of symbol;
   sps			: array [char] of symbol;
   t, a, b, sx, c1, c2	: integer;        {indices to tables}
   iflag, oflag, skipflag, prtsource, stackdump, prtables	: boolean;

   stantyps	: typset;
   display	: array[0..lmax] of integer;

   tab		: array[0..tmax] of packed record
					     name	: alfa;
					     link	: index;
					     obj	: object;
					     typ	: types;
					     ref1	: index;
					     normal	: boolean;
					     lev	: 0.. lmax;
					     adr	: integer;
					   end;

   atab		: array[1..amax] of packed record
					     inxtyp, eltyp	: types;
					     elref, low, high,
					     elsize, size	: index;
					   end;

   btab		: array[1..bmax] of packed record
					     last, lastpar,
					     psize, vsize	: index;
					   end;

   stab		: packed array[0..smax] of char; { string table }

   rconst	: array[1..c2max] of double;

   code		: array[0..cmax] of order;

   linerror	: array[0..cmax] of integer;
   prtlin	: boolean;

   inf, outf		: str40;
   filnb, fnamelen	: integer;


PROCEDURE abort;

begin
  ErrFlag:=true; goto 9999
end;

   PROCEDURE errormsg;

   var	k	: integer;
	msg	: array[0..ermax] of alfa;

   begin
      msg[ 0] := 'undef id       '; msg[ 1] := 'multi def      ';
      msg[ 2] := 'identifier     '; msg[ 3] := 'program        ';
      msg[ 4] := ') expected     '; msg[ 5] := ': expected     ';
      msg[ 6] := 'syntax         '; msg[ 7] := 'ident, var     ';
      msg[ 8] := 'of             '; msg[ 9] := '(expected     ';
      msg[10] := 'id, array      '; msg[11] := '[ expected     ';
      msg[12] := '] expected     '; msg[13] := '..             ';
      msg[14] := '; expected     '; msg[15] := 'func. type     ';
      msg[16] := '=              '; msg[17] := 'boolean        ';
      msg[18] := 'convar typ     '; msg[19] := 'type           ';
      msg[20] := 'prog.param     '; msg[21] := 'too big        ';
      msg[22] := '.              '; msg[23] := 'typ (case)     ';
      msg[24] := 'character      '; msg[25] := 'const id       ';
      msg[26] := 'index type     '; msg[27] := 'indexbound     ';
      msg[28] := 'no array       '; msg[29] := 'type id        ';
      msg[30] := 'undef type     '; msg[31] := 'no record      ';
      msg[32] := 'boole type     '; msg[33] := 'arith type     ';
      msg[34] := 'integer        '; msg[35] := 'types          ';
      msg[36] := 'param type     '; msg[37] := 'varib id       ';
      msg[38] := 'string         '; msg[39] := 'no.of pars     ';
      msg[40] := 'real numbr     '; msg[41] := 'type           ';
      msg[42] := 'real type      '; msg[43] := 'integer        ';
      msg[44] := 'var, const     '; msg[45] := 'var, proc      ';
      msg[46] := 'types (:=)     '; msg[47] := 'typ (case)     ';
      msg[48] := 'type           '; msg[49] := 'store ovfl     ';
      msg[50] := 'constant       '; msg[51] := ':= expectd     ';
      msg[52] := 'then           '; msg[53] := 'until          ';
      msg[54] := 'do             '; msg[55] := 'to downto      ';
      msg[56] := 'begin          '; msg[57] := 'end            ';
      msg[58] := 'factor         ';

      writeln(psout); writeln;
      writeln(psout,' Key words'); writeln(' Keywords');
      k:=0;
      while errs <> [] do
      begin
         while not (k in errs ) do  k:=k+1;
         writeln(psout,k:2,'  ', msg[k]); writeln(k:2,'  ', msg[k]);
         errs := errs - [k];
      end  {  while errs  }
   end  {errormsg};

  PROCEDURE endskip;

  begin {underline skips part of input}
    while errpos < cc do
    begin
      write(psout,'-'); write('-');
      errpos := errpos + 1
    end;
    skipflag  := false
  end  {endskip};

  PROCEDURE nextch;
  {read   next character; process line end}

  begin { nextch }
    if cc = ll then
    begin
      if eof(psin[incllev]) then
      begin
	if incllev = 0 then
        begin
          writeln(psout); writeln;
          writeln(psout, 'program incomplete'); writeln('program incomplete');
          errormsg;
          abort;
        end
        else
        begin
	  close(psin[incllev]); incllev := pred(incllev);
        end;
      end;

      if errpos <> 0 then
      begin
        if skipflag then endskip;
        writeln(psout); writeln;
        errpos := 0;
      end;

      linecounter := linecounter + 1;
      prtlin := false;

      if prtsource {and (incllev = 0)} then
	write(psout, linecounter:4, ' ', lc:5, chr(tabsym));

      ll := 0; cc := 0;
      while not eoln(psin[incllev]) do
      begin
	ll := succ(ll);
        read(psin[incllev], chin);
        if prtsource {and (incllev = 0)} then write(psout, chin);
	line[ll] := chin;
      end;
      ll := succ(ll);
      readln(psin[incllev]); line[ll] := ' ';
      {read(psin[incllev], line[ll]);}
      if prtsource {and (incllev = 0)} then writeln(psout);
    end;
    cc := cc + 1;
    chin := line[cc];
  end; { nextch }


  PROCEDURE error(n: integer);

  var	i, cc0	: integer;

  begin
    if not prtlin then
    begin
      if not prtsource {or (incllev <> 0)} then
      begin
        write(psout, linecounter:4, ' ', lc:5, chr(tabsym));
        for i:=1 to ll do
          write(psout, line[i]);
        writeln(psout);
      end;
      write(linecounter:4, ' ', lc:5, chr(tabsym));
      for i:=1 to ll do
        write(line[i]);
      writeln;
      prtlin := true;
    end;

    if errpos = 0 then
    begin
      write(psout, ' *********', chr(tabsym));
      write(' *********', chr(tabsym));
      cc0 := cc - 2;
    end
    else
      cc0 := cc;
    if cc > errpos then
    begin
      for i:=1 to cc0-errpos do
      begin
        if line[errpos+i] = chr(tabsym) then
        begin
          write(psout, chr(tabsym)); write(chr(tabsym));
        end
	else
        begin
          write(psout, ' '); write(' ');
        end;
      end;
      write(psout, '^', n:2); write('^', n:2);
      errpos := cc + 3;
      errs := errs + [n]
    end
  end  {error};

  PROCEDURE fatal(n: integer);

  var	msg: array[1..7] of alfa;

  begin
    writeln(psout); writeln;
    errormsg;
    msg[1]    := 'identifier     ';
    msg[2]    := 'PROCEDUREs     ';
    msg[3]    := 'reals          ';
    msg[4]    := 'arrays         ';
    msg[5]    := 'levels         ';
    msg[6]    := 'code           ';
    msg[7]    := 'strings        ';

    writeln(psout,'Compiler table for ', msg[n], ' is too small');
    writeln('Compiler table for ', msg[n], ' is too small');
    abort;
  end  {fatal};

PROCEDURE insymbol;     { reads next symbol  }

{ Scanner }

label	1, 2, 3;

var	i, j, k, e	: integer;
	str1		: alfa;
	parsename	: boolean;


  PROCEDURE readscale;

  var   s, sign : integer;

  begin
    nextch;
    sign:=1; s:=0;
    if chin='+' then nextch
    else  if chin='-' then begin
              nextch;  sign:=-1;
              end;
    if not ((chin>='0') and (chin<='9'))
    then  error(40)
    else  repeat
            s:=10*s+ord(chin)-ord('0');
            nextch
          until not ((chin>='0') and (chin<='9'));
    e:=s*sign+e
  end  { readscale  };

  PROCEDURE adjustscale;

  var   s : integer;
        d,t : double;

  begin
    if k+e > emax  then error(21)
    else if k+e < emin
         then rnum:=0.0
         else begin
           s:=abs(e); t:=1.0; d:=10.0;
           repeat
             while not odd(s)  do
             begin
               s:=s div 2;  d:=sqr(d)
             end;
               s:=s-1; t:=d*t
           until s=0;

           if e>=0 then rnum:=rnum*t
                   else rnum:=rnum/t
         end
  end  {  adjustscale  };

  PROCEDURE options;

    PROCEDURE switch(var b:boolean);

      begin
        b:=chin='+';
        if not b
        then if not (chin='-')
             then begin { print error message  }
               while (chin<>'*') and (chin<>',') do nextch;
             end
             else nextch
        else nextch
     end  { switch };

     begin  { options  }
       repeat
         nextch;
         if chin <> '*' then
	 begin
           if chin = 't' then
	   begin
             nextch; switch(prtables)
           end
	   else if chin = 's' then
	   begin
             nextch; switch(stackdump)
           end
	   else if chin = 'p' then
	   begin
             nextch; switch(prtsource)
           end;
         end
       until chin<>','
   end { options  };

   begin {insymbol}
     parsename := false;
1:   while (chin = ' ') or (chin = chr(tabsym)) do nextch;
      case chin of
         'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i',
         'j', 'k', 'l', 'm',
         'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v',
         'w', 'x', 'y', 'z',
         'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I',
         'J', 'K', 'L', 'M',
         'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V',
         'W', 'X', 'Y', 'Z', '_', '"':
             begin {identifier  or wordsymbol}
                k := 0;
                id := blankname;
                repeat
		   if chin = '"' then parsename := not parsename;
                   if k < NameLength then
                      begin
                         k := k + 1;
                         { upper case to lower case }
                         if chin in ['A'..'Z'] then
                           chin := chr(ord(chin) - ord('A') + ord('a'));
                         id[k] := chin
                      end;
                   nextch
                until not parsename and
                      not (chin in ['a'..'z', 'A'..'Z', '0'..'9', '_']);
                i := 1;
                j := nkw;               {binary search}
                repeat
                   k := (i + j) div 2;
                   if id <= key[k] then
                      j := k - 1;
                   if id >= key[k] then
                      i := k + 1
                until i > j;
                if i - 1 > j
                then
                   sy := ksy[k]
                else
                   sy := ident
             end;
          '0', '1', '2', '3', '4', '5', '6', '7', '8', '9':
             begin {number}
                k := 0;
                inum := 0;
                sy := intcon;
                repeat
                   inum := inum * 10 + ord(chin) - ord('0');
                   k := k + 1;
                   nextch
                until not (chin in ['0' .. '9']);
                if (k > kmax) or (inum > nmax) then
                   begin
                      error(21);
                      inum := 0;
                      k := 0
                   end;
                 if chin='.'
                 then begin nextch;
                    if chin='.'
                    then chin:=':'
                    else begin
                      sy:=realcon; rnum:=inum; e:=0;
                      while (chin >= '0') and (chin <= '9' ) do
                      begin
                        e:=e-1;
                        rnum := 10.0*rnum + (ord(chin)-ord('0'));
                        nextch
                      end;
                      if e = 0 then error(40);
                      if (chin = 'd') or (chin = 'D') or
			 (chin = 'e') or (chin = 'E') then readscale;
                      if e <> 0 then adjustscale
                    end
                  end else
                    if (chin = 'd') or (chin = 'D') or
		       (chin = 'e') or (chin = 'E') then
		    begin
                      sy := realcon; rnum := inum; e := 0;
                      readscale;
                      if e <> 0 then adjustscale
                    end;
             end;

          ':' {, col}:
             begin
                nextch;         {mod mh}
                if chin = '='
                then
                   begin
                      sy := becomes;
                      nextch
                   end
                else
                   sy := colon
             end;
          '<':
             begin
                nextch;
                if chin = '='
                then
                   begin
                      sy := leq;
                      nextch
                   end
                else
                   if chin = '>'
                   then
                      begin
                         sy := neq;
                         nextch
                      end
                   else
                      sy := lss
             end;
          '>':
             begin
                nextch;
                if chin = '='
                then
                   begin
                      sy := geq;
                      nextch
                   end
                else
                   sy := gtr
             end;
          '.':
             begin
                nextch;
                if chin = '.'
                then
                   begin
                      sy := colon;
                      nextch
                   end
                else
                   sy := period
                end;
          '''':
             begin
                k := 0;
             2: nextch;
                if chin = '''' then
                   begin
                      nextch;
                      if chin <> '''' then  goto 3
                   end;
                if sx + k = smax then
                   fatal(7);
                stab[sx+k] := chin;
                k := k + 1;
                if cc = 1
                then
                   begin {end of line}
                      k := 0;
                   end
                else
                   goto 2;
             3: if k = 1 then
                begin
                  sy := charcon;
                  inum := ord(stab[sx])
                end
                else if k = 0 then
                begin
{                 error(38);}
                  sy := charcon;
                  inum := 0
                end
                else
                begin
                  sy := stringcon;
                  inum := sx;
                  sleng := k;
                  sx := sx + k
                end
             end;
          '(' : begin
                  nextch;
                  if chin <> '*' then
                    sy := lparent
                  else
                  begin {comment}
                    nextch;
                    if chin = '$' then options;
                    repeat
                      while chin <> '*' do
                        nextch;
                      nextch
                    until chin = ')';
                    nextch;
                    goto 1
                  end
                end;
	  '{' : begin
                  repeat
		    nextch;
                  until chin = '}';
                  nextch;
                  goto 1
		end;
          '+', '-', '/', '*', ')', '=', ',', ';','[',']' :
             begin
                sy := sps[chin];
                nextch
             end;
	  '%' : begin{ include statement }
		  k := 0; str1 := blankname;
		  while chin <> ' ' do
		  begin
		    nextch;
		    k := succ(k); str1[k] := chin;
		  end;
		  if str1 = 'include        ' then
		    write('  include:')
		  else
		    writeln('%include expected');
		  while chin = ' ' do
		    nextch;
		  if chin <> '''' then writeln(''' expected');
		  k := 0;
		  while chin <> ' ' do
		  begin
		    nextch;
		    k := succ(k); fname[k] := chin;
		  end;
		  if fname[k-1] = '''' then
		  begin
		    fname[k-1] := ' '; k := k - 2;
		    for i:=k+1 to 80 do
		      fname[i] := ' ';
		    writeln(' ', fname:k);
		    incllev := succ(incllev);
		    reset_(psin[incllev], fname);
		  end
		  else
		    writeln(''' expected');
		  goto 1;
		end;
          '$', '!', '@', '^', '?','&', '\' :
             begin
                error(24);
                nextch;
                goto 1
              end;
      otherwise
         begin
           error(24); nextch; goto 1
         end;
      end
   end  {insymbol};

PROCEDURE enterstandardids(x0: alfa;  x1: object;
                           x2: types; x3: integer);

begin
  t := t+1;
  with tab[t] do
  begin
    name := x0;  link := t-1;  obj := x1;
    typ  := x2;  ref1  := 0;    normal := true;
    lev := 0;  adr := x3;
  end;
end; { enterstandardids }

PROCEDURE enterarray(tp: types; l, h: integer);

   begin
      if l > h  then  error(27);
      if (abs(l) > xmax) or (abs(h) > xmax) then
          begin
             error(27);
             l := 0;
             h := 0;
          end;
      if a = amax
      then
          fatal(4)
      else
          begin
             a := a + 1;
             with atab[a] do
                begin
                   inxtyp := tp;
                   low := l;
                   high := h
                end
          end {enterarray};
   end  {enterarray};

PROCEDURE enterblock;

   begin
      if b = bmax
      then
          fatal(2)
      else
          begin
             b := b + 1;
             btab[b].last := 0;
             btab[b].lastpar := 0
          end
   end  {enterblock};

PROCEDURE enterreal(x : double);

begin
  if c2 = c2max-1  then  fatal(3)
  else  begin
    rconst[c2+1] := x;  c1 := 1;
    while rconst[c1] <> x do c1 := c1 + 1;
    if c1 > c2 then c2 := c1
  end
end  {  enterreal  };

PROCEDURE emit(fct: integer);

  begin
    if lc = cmax then fatal(6);

    linerror[lc] := linecounter;

    code[lc].f := fct;
    lc := lc  + 1;
  end  {emit};

PROCEDURE emit1(fct, b: integer);

  begin
    if lc = cmax then fatal(6);

    linerror[lc] := linecounter;

    with code[lc] do
    begin
      f := fct;
      y := b;
    end;
    lc := lc  + 1;
  end  {emit1};

PROCEDURE emit2(fct, a, b: integer);

  begin
    if lc = cmax then fatal(6);

    linerror[lc] := linecounter;

    with code[lc] do
    begin
      f := fct;
      x := a;
      y := b
    end;
    lc := lc  + 1;
  end  {emit2};

PROCEDURE printtables;

type
        alfa1 = packed array [1..10] of char;
        alfa2 = packed array [1..4]  of char;

var     i		:  integer;
        currentinstr	: order;
	instr		: array [0..63] of packed array [1..6] of char;
        idt1		: alfa1;
        idt2		: alfa2;


  procedure initinstructions;

  begin
    instr[0] :='lodadd'; instr[1] :='lodval'; instr[2] :='lodind';
    instr[3] :='upddis'; instr[4] :='nop   '; instr[5] :='nop   ';
    instr[6] :='nop   '; instr[7] :='nop   '; instr[8] :='stfunc';
    instr[9] :='offset'; instr[10]:='jump  '; instr[11]:='jmpfal';
    instr[12]:='switch'; instr[13]:='nop   '; instr[14]:='forfup';
    instr[15]:='forsup'; instr[16]:='forfdn'; instr[17]:='forsdn';
    instr[18]:='markst'; instr[19]:='call  '; instr[20]:='index1';
    instr[21]:='indexm'; instr[22]:='loadbl'; instr[23]:='copybl';
    instr[24]:='ldint '; instr[25]:='ldreal'; instr[26]:='float ';
    instr[27]:='read  '; instr[28]:='wrstrg'; instr[29]:='write ';
    instr[30]:='wr1fld'; instr[31]:='halt  '; instr[32]:='leavep';
    instr[33]:='leavef'; instr[34]:='lodins'; instr[35]:='not   ';
    instr[36]:='negate'; instr[37]:='wr2fld'; instr[38]:='store ';
    instr[39]:='equalr'; instr[40]:='noteqr'; instr[41]:='lessr ';
    instr[42]:='lesser'; instr[43]:='greal '; instr[44]:='gereal';
    instr[45]:='equal '; instr[46]:='notequ'; instr[47]:='lessth';
    instr[48]:='lesseq'; instr[49]:='greatt'; instr[50]:='greate';
    instr[51]:='or    '; instr[52]:='add   '; instr[53]:='sub   ';
    instr[54]:='addrea'; instr[55]:='subrea'; instr[56]:='and   ';
    instr[57]:='mult  '; instr[58]:='div   '; instr[59]:='mod   ';
    instr[60]:='multr '; instr[61]:='rdiv  '; instr[62]:='readln';
    instr[63]:='writln';
  end;


begin
  initinstructions;
  writeln(psout);  writeln(psout);  writeln(psout);
  writeln(psout,'   Identifiers         Link  Obj   Typ  Ref  Nrm  Lev  Adr   Remarks');
  writeln(psout);
  for i:=btab[1].last to t do
     with tab[i] do
       begin
         case obj of
           konstant:  idt1 := '  constant';
           vvariable: idt1 := '  variable';
           type1:     idt1 := ' user type';
           prozedure: idt1 := ' PROCEDURE';
           funktion:  idt1 := '  FUNCTION';
           otherwise    idt1 := ' undefined';
         end;
         case typ of
           notyp:   idt2 := 'proc';
           ints:    idt2 := ' int';
           reals:   idt2 := 'real';
           bools:   idt2 := 'bool';
           chars:   idt2 := 'char';
           arrays:  idt2 := 'arry';
           records: idt2 := ' rec';
	   texts:   idt2 := 'text';
           otherwise  idt2 := 'undf';
         end;
         writeln(psout, i:5, ' ', name, link:5, ord(obj):5, 
                 '   ', idt2:4, ref1:4, 
                 ord(normal):5, lev:5, adr:5, ' ', idt1:10);
       end;

  writeln(psout);  writeln(psout);  writeln(psout);
  writeln(psout,'blocks    Last Lpar Psze Vsze');
  writeln(psout);
  for i:=1 to b do
    with btab[i] do
      writeln(psout,i:4,last:9,lastpar:5,psize:5,vsize:5);

  writeln(psout);  writeln(psout);  writeln(psout);
  writeln(psout,'Arrays    Xtyp Etyp Eref   Low High Elsz Size');
  writeln(psout);

  for i:=1 to a do
    with atab[i] do
      writeln(psout,i:4,ord(inxtyp):9,ord(eltyp):5,
              elref:5,low:5,high:5,elsize:5,size:5);

  writeln(psout);  writeln(psout);  writeln(psout);

  writeln(psout, ' code:');
  for i:=0 to lc-1 do
  begin
    write(psout, i:4, ' ');
    currentinstr := code[i];
    if currentinstr.f <= 63 then
    begin
      write(psout, instr[currentinstr.f], ' ');
      if currentinstr.f < 31 then
        if currentinstr.f <= 3 then
 	  write(psout, currentinstr.x:2, currentinstr.y:5)
        else
	  write(psout, currentinstr.y:7)
      else
        write(psout, '      ');
      writeln(psout);
    end
    else
      writeln(psout, 'user def.');
  end;
end { printtables };

PROCEDURE block(fsys: symset; isfun: boolean; level: integer);

{ Parser }

   type
      conrec =  record
		  case tp: types of
                    ints, chars, bools:		(i : integer);
                    reals:			(r : double);
		    notyp, arrays, records:	();
                end;
   var
      dx: integer;      {data   allocation index}
      prt: integer;     {t-index of this PROCEDURE}
      prb: integer;     {b-index of this PROCEDURE}
      x  : integer;

   PROCEDURE skip(fsys: symset; n: integer);
      begin
          error(n);
          skipflag := true;
          while not (sy in fsys) do
             insymbol;
          if skipflag then
             endskip;
      end {skip};

   PROCEDURE test(s1, s2: symset; n: integer);

   begin
     if not (sy in s1) then skip(s1 + s2, n)
   end {test};


   PROCEDURE testsemicolon;

      begin
          if sy = semicolon
          then
             insymbol
          else begin
             error(14);
             if sy in [comma,colon] then insymbol
          end;
          test([ident,incsy] + blockbegsys,   fsys, 6);
      end {testsemicolon};

    PROCEDURE  enter(id: alfa; k: object);
      var
          j, l: integer;

      begin
          if t = tmax
          then
             fatal(1)
          else
             begin
                tab[0].name := id;
                j := btab[display[level]].last;
                l := j;
                while tab[j].name <> id do
                   j := tab[j].link;
                if j <> 0
                then
                   error(1)
                else
                   begin
                      t := t + 1;
                      with tab[t] do
                         begin
                            name := id;
                            link := l;
                            obj := k;
                            typ := notyp;
                            ref1 := 0;
                            lev := level;
                            adr := 0
                         end;
                      btab[display[level]].last := t
                   end
             end
      end {enter};


   FUNCTION loc(id: alfa): integer;

      var
          i, j: integer;        {locate id in table}

      begin
          i := level;
          tab[0].name   := id;  {sentinel}
          repeat
             j := btab[display[i]].last;
             while tab[j].name <> id do
                j := tab[j].link;
             i := i - 1;
          until (i < 0) or (j <> 0);
          if j = 0 then
             error(0);
          loc := j
      end {loc};

   PROCEDURE entervariable;

      begin
          if sy = ident
          then
             begin
                enter(id, vvariable);
                insymbol
             end
          else
             error(2)
      end {entervariable};

   PROCEDURE constant(fsys: symset; var c: conrec);

      var
          x, sign: integer;

      begin
          c.tp := notyp;
          c.i := 0;
          test(constbegsys, fsys, 50);
          if sy in constbegsys
          then
             begin
                if sy = charcon
                then
                   begin
                      c.tp := chars;
                      c.i := inum;
                      insymbol
                   end
                else
                   begin
                      sign := 1;
                      if sy in [plus, minus] then
                         begin
                            if sy = minus then
                               sign := - 1;
                            insymbol
                         end;
                      if sy = ident
                      then
                         begin
                            x := loc(id);
                            if x <> 0 then
                               if tab[x].obj <> konstant
                               then
                                  error(25)
                               else
                                  begin
                                     c.tp := tab[x].typ;
                                     if c.tp = reals
                                     then c.r:=sign*rconst[tab[x].adr]
                                     else c.i:=sign*tab[x].adr
                                  end;
                            insymbol
                         end
                      else
                         if sy = intcon
                         then
                            begin
                               c.tp := ints;
                               c.i := sign * inum;
                               insymbol
                            end else if sy = realcon
                                     then begin
                                       c.tp := reals; c.r := sign*rnum;
                                       insymbol
                                     end else skip(fsys,50)
                   end;
                test(fsys, [], 6);
             end
      end {constant};


   PROCEDURE typ(fsys:  symset; var tp: types; var rf, sz: integer);

      var
          x: integer;
          eltp: types;
          elrf: integer;
          elsz, offset, t0, t1: integer;


      PROCEDURE arraytyp(var aref, arsz: integer);

          var
             eltp: types;
             low, high: conrec;
             elrf, elsz: integer;

          begin
             constant([colon, rbrack, rparent, ofsy] +  fsys, low);
             if low.tp = reals
             then begin
               error(27);
               low.tp := ints;  low.i := 0
             end;
             if sy = colon
             then
                insymbol
             else
                error(13);
             constant([rbrack,  comma, rparent, ofsy] +  fsys, high);
             if high.tp <> low.tp then
                begin
                   error(27);
                   high.i := low.i
                end;
             enterarray(low.tp, low.i, high.i);
             aref := a;
             if sy = comma
             then
                begin
                   insymbol;
                   eltp := arrays;
                   arraytyp(elrf, elsz)
                end
             else
                begin
                   if sy = rbrack
                   then
                      insymbol
                   else begin
                      error(2);
                      if sy = rparent then insymbol
                   end;
                   if sy = ofsy
                   then
                      insymbol
                   else
                      error(8);
                   typ(fsys, eltp, elrf, elsz)
                end;
             with atab[aref] do
                begin
                   arsz := (high - low + 1) * elsz;
                   size := arsz;
                   eltyp := eltp;
                   elref := elrf;
                   elsize := elsz
                end;
          end {arraytyp};


      begin {typ}
          tp := notyp;
          rf := 0;
          sz := 0;
          test(typebegsys, fsys, 10);
          if sy in typebegsys
          then
             begin
                if sy = ident
                then
                   begin
                      x := loc(id);
                      if x <> 0 then
                         with tab[x] do
                            if obj <> type1
                            then
                               error(29)
                            else
                               begin
                                  tp := typ;
                                  rf := ref1;
                                  sz := adr;
                                  if tp = notyp then
                                     error(30)
                               end;
                      insymbol
                   end
                else
                   if sy = arraysy
                   then
                      begin
                         insymbol;
                         if sy = lbrack
                         then
                            insymbol
                         else begin
                            error(11);
                            if sy = lparent
                            then insymbol
                         end;
                         tp := arrays;
                         arraytyp(rf, sz)
                      end
                   else begin  { records }
                     insymbol;
                     enterblock;
                     tp := records; rf := b;
                     if level = lmax then fatal(5);
                     level := level+1; display[level] := b; offset := 0;
                     while not (sy in fsys-[semicolon,comma,ident]+[endsy]) do
                     begin { field section }
                       if sy = ident
                       then begin
                         t0 := t; entervariable;
                         while sy = comma do
                         begin
                           insymbol; entervariable
                         end;
                         if sy = colon then insymbol else error(5);
                         t1 := t;
                       typ(fsys+[semicolon,endsy,comma,ident],eltp,elrf,elsz);
                      while t0 < t1 do
                      begin
                        t0 := t0+1;
                        with tab[t0] do
                        begin
                          typ := eltp;
                          ref1 := elrf;   normal := true;
                          adr := offset; offset := offset + elsz
                        end
                      end
                    end;  { sy = ident }
                    if sy <> endsy
                    then begin
                      if sy = semicolon
                      then insymbol
                      else begin
                        error(14);
                        if sy = comma then insymbol
                      end;
                      test([ident,endsy,semicolon],fsys,6)
                    end
                  end;  { field section }

                  btab[rf].vsize := offset; sz := offset;
                  btab[rf].psize := 0;
                  insymbol; level := level-1
               end;  { records }
               test(fsys, [], 6);
             end
      end {typ};


   PROCEDURE parameterlist;  {formal parameter list}

      var
          tp: types;
          rf, sz, x, t0: integer;
          valpar: boolean;

      begin
          insymbol;
          tp := notyp;
          rf := 0;
          sz := 0;
          test([ident,  varsy], fsys +  [rparent], 7);
          while sy in [ident, varsy] do
             begin
                if sy <> varsy
                then
                   valpar := true
                else
                   begin
                      insymbol;
                      valpar := false
                   end;
                t0 := t;
                entervariable;
                while sy = comma do
                   begin
                      insymbol;
                      entervariable;
                   end;
                if sy = colon
                then
                   begin
                      insymbol;
                      if sy <> ident
                      then
                         error(2)
                      else
                         begin
                            x := loc(id);
                            insymbol;
                            if x <> 0 then
                               with tab[x] do
                                  if obj <> type1
                                  then
                                     error(29)
                                  else
                                     begin
                                        tp := typ;
                                        rf := ref1;
                                        if valpar
                                        then
                                           sz := adr
                                        else
                                           sz := 1
                                     end;
                         end;
                      test([semicolon,  rparent], [comma, ident] +
                         fsys, 14)
                   end
                else
                   error(5);
                while t0 < t do
                   begin
                      t0 := t0 + 1;
                      with tab[t0] do
                         begin
                            typ := tp;
                            ref1 := rf;
                            normal := valpar;
                            adr := dx;
                            lev := level;
                            dx := dx + sz
                         end
                   end;
                if sy <> rparent then
                   begin
                      if sy = semicolon
                      then
                         insymbol
                      else begin
                         error(14);
                         if sy = comma then insymbol
                      end;
                      test([ident, varsy], [rparent] + fsys, 6);
                   end
             end;  { while }
          if sy = rparent
          then
             begin
                insymbol;
                test([semicolon, colon], fsys, 6);
             end
          else
             error(4)
      end {parameterlist};


   PROCEDURE constdec0;

     var NopTemp:integer;

     begin
       for NopTemp:=1 to globval.Elem_nFam do
       begin
         id:=ElemFam[NopTemp].ElemF.PName;
         enter(id, konstant);
         with tab[t] do
         begin
           typ := ints;
           ref1 := 0;
           lev := level;
           adr := NopTemp;
           normal := true;
         end
       end;
     end;


   PROCEDURE constdec;

      var
          c: conrec;

      begin
          insymbol;
          test([ident], blockbegsys, 2);
          while sy = ident do
             begin
                enter(id, konstant);
                insymbol;
                if sy = eql
                then
                   insymbol
                else begin
                   error(16);
                   if sy = becomes then insymbol
                end;
                constant([semicolon, comma, ident] + fsys, c);
                tab[t].typ :=   c.tp;
                tab[t].ref1 :=   0;
                if c.tp = reals
                then begin
                  enterreal(c.r); tab[t].adr := c1
                end else tab[t].adr := c.i;
                testsemicolon
              end
      end { constdeclaration };

   PROCEDURE typedeclaration;

      var
          tp: types;
          rf, sz, t1: integer;

      begin
          insymbol;
          test([ident], blockbegsys, 2);
          while sy = ident do
             begin
                enter(id, type1);
                t1 := t;
                insymbol;
                if sy = eql
                then
                   insymbol
                else begin
                   error(16);
                   if sy = becomes then insymbol
                end;
                typ([semicolon, comma,  ident]  + fsys, tp, rf, sz);
                with tab[t1] do
                   begin
                      typ := tp;
                      ref1 := rf;
                      adr := sz
                   end;
                testsemicolon
             end
      end {typedeclaration};


   PROCEDURE vardeclaration;

      var
          t0, t1, rf, sz: integer;
          tp: types;

      begin
        insymbol;
        while sy = ident do
        begin
          t0 := t;
          entervariable;
          while sy = comma do
          begin
            insymbol;
            entervariable;
          end;
          if sy = colon then
            insymbol
          else
            error(5);
          t1 := t;
          typ([semicolon, comma,  ident]  + fsys, tp, rf, sz);
          while t0 < t1 do
          begin
            t0 := t0 + 1;
            with tab[t0] do
            begin
              typ := tp;
	      if typ = texts then
	      begin
		filnb := succ(filnb);
		ref1 := filnb;
	      end
	      else
                ref1 := rf;
              lev := level;
              adr := dx;
              normal := true;
              dx := dx + sz
            end;
          end;
          testsemicolon
        end;
      end {variab|edeclaration};

   PROCEDURE procdeclaration;

      var
          isfun: boolean;

      begin
          isfun := sy = funcsy;
          insymbol;
          if sy <> ident then
             begin
                error(2);
                id := blankname;
             end;
          if isfun
          then
             enter(id, funktion)
          else
             enter(id, prozedure);
          tab[t].normal := true;
          insymbol;
          block([semicolon] +   fsys, isfun, level + 1);
          if sy = semicolon
          then
             insymbol
          else
             error(14);
          emit(32 + ord(isfun)) {exit}
      end {PROCEDUREdeclaration};

{**************************************************}
{*                                                *}
{*      S T A T E M E N T                         *}
{*                                                *}
{**************************************************}

   PROCEDURE statement(fsys: symset);

   var    i: integer;
          x: item;

   PROCEDURE expression(fsys: symset; var x: item);
   forward;

      PROCEDURE selector_type(fsys:  symset; var v: item; mytyp:types);

          var
             x: item;
             a, j: integer;

          begin
             repeat
               if sy = period
               then begin
                 insymbol;
                 if sy <> ident
                 then error(2)
                 else begin
                    if v.typ <> records
                    then error(31)
                    else begin
                      j:=btab[v.ref1].last;
                      tab[0].name := id;
                      while tab[j].name <> id do j:=tab[j].link;
                      if j = 0 then error(0);
                      v.typ := tab[j].typ;
                      v.ref1 := tab[j].ref1;
                      a := tab[j].adr;
                      if a <> 0 then emit1(9,a)
                    end;
                    insymbol
                  end
             end else begin { array selector }
             if sy <> lbrack then
                error(11);
             repeat
                insymbol;
                expression(fsys + [comma, rbrack], x);
                if v.typ <> mytyp { arrays }
                then
                   error(28)
                else
                   begin
                      a := v.ref1;
                      if atab[a].inxtyp <> x.typ
                      then
                         error(26)
                      else if atab[a].elsize = 1
                           then emit1(20,a)
                           else emit1(21,a);
                      v.typ := atab[a].eltyp;
                      v.ref1 := atab[a].elref
                   end
             until sy <> comma;
             if sy = rbrack
             then
                insymbol
             else begin
                error(12);
                if sy = rparent then insymbol
             end
           end
          until not (sy in [lbrack,lparent,period]);
          test(fsys, [], 6);
          end {selector};

      PROCEDURE selector(fsys:  symset; var v: item);
          begin
            selector_type(fsys,v,arrays);
          end {selector};


      PROCEDURE call(fsys: symset; i: integer);

          var
             x: item;
             lastp, cp, k: integer;

          begin
            emit1(18, i);              {markstack}
            lastp := btab[tab[i].ref1].lastpar;
            cp := i;
            if sy = lparent then
            begin {actual parameter list}
              repeat
                insymbol;
                if cp >= lastp then
                  error(39)
                else
                begin
                  cp := cp + 1;
                  if tab[cp].normal then
                  begin {value parameter}
                    expression(fsys + [comma, colon, rparent], x);
                    if x.typ = tab[cp].typ then
                    begin
                      if x.ref1 <> tab[cp].ref1 then
                        error(36)
                      else if x.typ = arrays then
			emit1(22, atab[x.ref1].size)
                      else if x.typ=records then
			emit1(22,btab[x.ref1].vsize)
                      end
		      else if (x.typ=ints) and (tab[cp].typ=reals) then
			emit1(26,0)
                      else if x.typ <> notyp then
			error(36);
                    end
                    else
                    begin     {variable      parameter}
                      if sy <> ident then
                        error(2)
                      else
                      begin
                        k := loc(id);
                        insymbol;
                        if k <> 0 then
                        begin
                          if tab[k].obj <> vvariable then
                            error(37);
                          x.typ := tab[k].typ;
                          x.ref1 := tab[k].ref1;
                          if tab[k].normal then
                            emit2(0, tab[k].lev, tab[k].adr)
                          else
                            emit2(1, tab[k].lev, tab[k].adr);
                          if sy in [lbrack,lparent,period] then
                            selector(fsys + [comma, colon, rparent], x);
                          if (x.typ <> tab[cp].typ) or
			     ((x.ref1 <> tab[cp].ref1) and (x.typ <> texts)) then
                             error(36)
                        end
                      end
                    end
                  end;
                  test([comma, rparent], fsys, 6);
              until sy <> comma;
              if sy = rparent then
                insymbol
              else
                error(4)
            end;
            if cp < lastp then
              error(39);           {too few actual parameters}
            emit1(19, btab[tab[i].ref1].psize - 1);
            if tab[i].lev < level then
              emit2(3, tab[i].lev, level)
          end {call};



	{ Modified }
	procedure getlparent;

	begin
	  if sy = lparent then insymbol else error(9);
	end;

	procedure getrparent;

	begin
	  if sy = rparent then insymbol else error(4);
	end;

	procedure getstrcon(komma : boolean);

	begin
	  if komma then
	  begin
	    if sy <> comma then error(6);
	    insymbol;
	  end;
          if sy = stringcon then
          begin
            emit1(24, sleng);
	    insymbol;
          end
	  else if sy = charcon then
	  begin
	    if inum <> 0 then
              emit1(24, 1)
	    else
              emit1(24, 0);
	    insymbol;
	  end
          else
	    error(38);
	end;

	procedure getparam(komma, valpar : boolean;
			   typ : types; var x : item );

	var	k	: integer;

	begin
	  if komma then
	  begin
	    if sy <> comma then error(6);
	    insymbol;
	  end;
	  if valpar then
	  begin {value parameter}
	    expression(fsys + [comma,colon,rparent],  x);
	    if x.typ = typ then
	    begin
	      if x.typ = arrays then
	        emit1(22, atab[x.ref1].size)
	      else if x.typ = records then
	        emit1(22, btab[x.ref1].vsize);
	    end
	    else if (x.typ = ints) and (typ = reals) then
	      emit1(26,0)
	    else if x.typ <> notyp then
	      error(36);
	  end
	  else
	  begin {variable parameter}
	    if sy <> ident then
	      error(2)
	    else
	    begin
	      k := loc(id);
	      insymbol;
	      if k <> 0 then
	      begin
		if tab[k].obj <> vvariable then error(37);
		x.typ := tab[k].typ;
		x.ref1 := tab[k].ref1;
		if tab[k].normal then
	          emit2(0, tab[k].lev, tab[k].adr)
		else
		  emit2(1, tab[k].lev, tab[k].adr);
		if sy in [lbrack,lparent,period] then
		  selector(fsys + [comma, colon, rparent], x);
		if x.typ <> typ then error(36);
	      end;
	    end;
	  end;
	end;


{**************************************************}
{*                                                *}
{*        F U N C T I O N                         *}
{*                                                *}
{**************************************************}

      FUNCTION  resulttype(a, b: types): types;

          begin
             if (a > reals) or (b > reals)
             then
                begin
                   error(33);
                   resulttype := notyp
                end
             else
                if (a = notyp) or (b = notyp)
                then
                   resulttype := notyp
                else if a = ints
                     then if b = ints
                          then resulttype := ints
                          else begin
                            resulttype := reals; emit1(26,1)
                          end
                     else begin
                       resulttype := reals;
                       if b=ints then emit1(26,0)
                     end
          end {resulttyp};

{**************************************************}
{*                                                *}
{*        E X P R E S S I O N                     *}
{*                                                *}
{**************************************************}

      PROCEDURE expression;

          var
             y: item;
             op: symbol;

{**************************************************}
{*                                                *}
{*        S I M P L E  E X P R E S S I O N        *}
{*                                                *}
{**************************************************}

          PROCEDURE simpleexpression(fsys: symset; var x: item);

             var
                y: item;
                op: symbol;



{**************************************************}
{*                                                *}
{*        T E R M                                 *}
{*                                                *}
{**************************************************}

             PROCEDURE term(fsys: symset; var x: item);

                var
                   y: item;
                   op: symbol;
                   ts: typset;



{**************************************************}
{*                                                *}
{*          F A C T O R                           *}
{*                                                *}
{**************************************************}

                PROCEDURE factor(fsys: symset; var x: item);

                   var	i, f, fnb	: integer;


                PROCEDURE standfct(n : integer);

		procedure gencode(ts : typset);

		begin
		  if x.typ in ts then
		    emit1(8, n)
		  else if x.typ <> notyp then
		    error(48);
		end;


                begin  { standard FUNCTION no. n }

                   case n of
                   0, 2:{ abs, sqr }
                       begin
	                  getlparent;
                          expression(fsys+[rparent], x);
                          tab[i].typ := x.typ;
                          if x.typ = reals then n := n+1;
			  getrparent;
			  gencode([ints, reals]);
                       end;
                   4, 5:{ odd, chr }
                       begin
                          getlparent;
 		          getparam(false, true, ints, x);
			  getrparent;
			  gencode([ints]);
                       end;
                   6:{ ord }
                       begin
	                  getlparent;
                          expression(fsys+[rparent], x);
			  getrparent;
			  gencode([ints, bools, chars]);
                       end;
                   7, 8:{ succ, pred }
                       begin
	                  getlparent;
                          expression(fsys+[rparent], x);
                          tab[i].typ := x.typ;
			  getrparent;
			  gencode([ints, bools, chars]);
                       end;
                   9, 10, 11, 12, 13, 14, 15, 16:
                       begin
	                  getlparent;
 		          getparam(false, true, reals, x);
			  getrparent;
			  gencode([ints, reals]);
                       end;
		   17, 18: begin { eof, eoln }
			     getlparent;
			     getparam(false, true, texts, x);
			     fnb := x.ref1;
			     getrparent;
			     emit2(8, fnb, n);
			   end;


      { Modified }
 19 : begin { pwr }
        getlparent;
        getparam(false, true, reals, x);
        getparam(true, true, reals, x);
        getrparent;
        emit1(8,  19);
      end;
 20 : begin { tan }
        getlparent;
        getparam(false, true, reals, x);
        getrparent;
        emit1(8,  20);
      end;
 21 : begin { cosh_ }
        getlparent;
        getparam(false, true, reals, x);
        getrparent;
        emit1(8,  21);
      end;
 22 : begin { sinh_ }
        getlparent;
        getparam(false, true, reals, x);
        getrparent;
        emit1(8,  22);
      end;
 23 : begin { tanh_ }
        getlparent;
        getparam(false, true, reals, x);
        getrparent;
        emit1(8,  23);
      end;
 24 : begin { ranf }
        emit1(8,  24);
      end;
 25 : begin { dtor }
        getlparent;
        getparam(false, true, reals, x);
        getrparent;
        emit1(8,  25);
      end;
 26 : begin { }
        emit1(8,  26);
      end;
 27 : begin { sign }
        getlparent;
        getparam(false, true, reals, x);
        getrparent;
        emit1(8,  27);
      end;
 28 : begin { getangle }
        getlparent;
        getparam(false, true, reals, x);
        getparam(true, true, reals, x);
        getrparent;
        emit1(8,  28);
      end;
 29 : begin { trmat }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, false, arrays, x);
        getrparent;
        emit1(8,  29);
      end;
 30 : begin { detmat }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, false, arrays, x);
        getrparent;
        emit1(8,  30);
      end;
 31 : begin { invmat }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, false, arrays, x);
        getrparent;
        emit1(8,  31);
      end;
 32 : begin { strlen_ }
        getlparent;
        getparam(false, false, records, x);
        getrparent;
        emit1(8,  32);
      end;
 33 : begin { strind }
        getlparent;
        getparam(false, false, records, x);
	getstrcon(true);
        getrparent;
        emit1(4, inum);
      end;
 34 : begin { Elem_getkval }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, true, ints, x);
        getparam(true, true, ints, x);
        getrparent;
        emit1(8,  34);
      end;
 35 : begin { min_ }
        getlparent;
        getparam(false, true, reals, x);
        getparam(true, true, reals, x);
        getrparent;
        emit1(8,  35);
      end;
 36 : begin { max_ }
        getlparent;
        getparam(false, true, reals, x);
        getparam(true, true, reals, x);
        getrparent;
        emit1(8,  36);
      end;
 37 : begin { normranf }
        emit1(8, 37);
      end;
 38 : begin { dble }
        getlparent;
        getparam(false, true, ints, x);
        getrparent;
        emit1(8, 38);
      end;
 39 : begin { getnKid }
        getlparent;
        getparam(false, true, ints, x);
        getrparent;
        emit1(8,  39);
      end;
 40 : begin { Elem_getpos }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, true, ints, x);
        getrparent;
        emit1(8,  40);
      end;
 41 : begin { }
        emit1(8,  41);
      end;
 42 : begin { Elem_getord }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, true, ints, x);
        getrparent;
        emit1(8,  42);
      end;
 43 : begin { clock }
        emit1(8,  43);
      end;
 44 : begin { }
        emit1(8,  44);
      end;
 45 : begin { }
        emit1(8,  45);
      end;
 46 : begin { }
        emit1(8,  46);
      end;


{*begin   GEM.NECK ***}

                end;{ of case}

		x.typ := tab[i].typ;
              end { standfct };



                 begin {factor}
                      x.typ := notyp;
                      x.ref1 := 0;
                      test(facbegsys, fsys, 58);
                      while sy in facbegsys do
                         begin
                            if sy = ident
                            then
                               begin
                                  i := loc(id);
                                  insymbol;
                                  with tab[i]   do
                                     case obj of
                                     konstant:
                                           begin
                                              x.typ := typ;
                                              x.ref1 := 0;
                                              if x.typ = reals
                                              then emit1(25,adr)
                                              else emit1(24,adr)
                                           end;
                                     vvariable:
                                           begin
                                              x.typ := typ;
                                              x.ref1 := ref1;
                                              if sy in [lbrack,lparent,period]
                                              then
                                                 begin
                                                    if normal
                                                    then
                                                       f := 0
                                                    else
                                                       f := 1;
                                                    emit2(f, lev, adr);
                                                    selector(fsys, x);
                                                    if x.typ in stantyps
                                                    then
                                                       emit(34)
                                                 end
                                              else
                                                 begin
                                                    if x.typ in stantyps
                                                    then
                                                       if normal
                                                       then
                                                          f := 1
                                                       else
                                                          f := 2
                                                    else
                                                       if normal
                                                       then
                                                          f := 0
                                                       else
                                                          f := 1;
                                                    emit2(f, lev, adr)
                                                 end
                                           end;
                              type1, prozedure:
                                           error(44);
                                      funktion:
                                           begin
                                              x.typ := typ;
                                              if lev <> 0
                                              then
                                                 call(fsys, i)
                                              else
                                                 standfct(adr)
                                           end
                                     end {case, with}
                               end
                            else
                               if sy in [charcon, intcon, realcon]
                               then begin
                                 if sy = realcon
                                 then begin
                                   x.typ := reals;
                                   enterreal(rnum);
                                   emit1(25,c1)
                                 end else begin
                                     if sy = charcon
                                     then
                                        x.typ := chars
                                     else
                                        x.typ := ints;
                                     emit1(24, inum);
                                 end;
                                 x.ref1 := 0;
                                 insymbol
                               end
                               else
                                  if sy = lparent
                                  then
                                     begin
                                        insymbol;
                                        expression(fsys + [rparent], x);
                                        if sy = rparent
                                        then
                                           insymbol
                                        else
                                           error(4)
                                     end
                                  else
                                     if sy = notsy then
                                        begin
                                           insymbol;
                                           factor(fsys, x);
                                           if x.typ = bools
                                           then
                                              emit(35)
                                           else
                                              if x.typ <> notyp then
                                                 error(32)
                                        end;
                            test(fsys, facbegsys, 6);
                         end {while}
                   end {factor};


                begin {term}
                   factor(fsys + [times, rdiv, idiv, imod, andsy], x);
                   while sy in [times, rdiv, idiv, imod, andsy] do
                      begin
                         op := sy;
                         insymbol;
                         factor(fsys + [times, rdiv, idiv, imod, andsy], y);
                         if op = times
                         then
                            begin
                               x.typ := resulttype(x.typ, y.typ);

                         case x.typ of
                          notyp: ;
                          ints:  emit(57);
                          reals: emit(60);
                         end
                     end else if op = rdiv
                         then begin
                            if x.typ = ints
                            then begin
                               emit1(26,1);
                               x.typ := reals
                            end;
                            if y.typ = ints
                            then begin
                               emit1(26,0);
                               y.typ := reals
                            end;
                            if (x.typ=reals) and (y.typ=reals)
                            then emit(61)
                            else begin
                              if (x.typ<>notyp) and (y.typ<>notyp)
                              then error(33);
                              x.typ := notyp
                            end
                         end else if op = andsy
                            then
                               begin
                                  if (x.typ = bools) and (y.typ = bools)
                                  then
                                     emit(56)
                                  else
                                     begin
                                        if (x.typ <> notyp) and (y.typ <>
                                           notyp)
                                        then
                                           error(32);
                                        x.typ := notyp
                                     end
                               end
                            else
                               begin {op in[idiv, imod]}
                                  if (x.typ = ints) and (y.typ = ints)
                                  then
                                     if op = idiv
                                     then
                                        emit(58)
                                     else
                                        emit(59)
                                  else
                                     begin
                                        if (x.typ <> notyp) and (y.typ <>
                                           notyp)
                                        then
                                           error(34);
                                        x.typ := notyp
                                     end
                               end
                      end { while }
                end {term};


             begin {simpleexpression}
                if sy in [plus, minus]
                then
                   begin
                      op := sy;
                      insymbol;
                      term(fsys + [plus, minus], x);
                      if x.typ > reals
                      then error(33)
                      else if op = minus then
		      { Modified }
			if x.typ = ints then
			  emit1(36, 0)
			else
			  emit1(36, 1);
                 end else term(fsys+[plus,minus,orsy], x);
                while sy in [plus, minus, orsy] do
                   begin
                      op := sy;
                      insymbol;
                      term(fsys + [plus, minus, orsy], y);
                      if op = orsy
                      then
                         begin
                            if (x.typ = bools) and (y.typ = bools)
                            then
                               emit(51)
                            else
                               begin
                                  if (x.typ <> notyp) and (y.typ <> notyp)
                                  then
                                     error(32);
                                  x.typ := notyp
                               end
                         end else begin
                           x.typ := resulttype(x.typ,y.typ);

                            case x.typ of
                     notyp:  ;
                     ints:   if op = plus
                             then emit(52)
                             else emit(53);
                     reals:  if op = plus
                             then emit(54)
                             else emit(55)
                            end { case }
                         end
                   end
             end {simpleexpression};


          begin {expression};
             simpleexpression(fsys + [eql, neq, lss, leq, gtr,  geq], x);
             if sy in [eql, neq, lss, leq, gtr, geq]
             then
                begin
                   op := sy;
                   insymbol;
                   simpleexpression(fsys, y);
                   if (x.typ in [notyp, ints, bools, chars]) and (x.typ
                      = y.typ)
                   then
                      case op of
                         eql:
                            emit(45);
                         neq:
                            emit(46);
                         lss:
                            emit(47);
                         leq:
                            emit(48);
                         gtr:
                            emit(49);
                         geq:
                            emit(50);
                      end
                   else begin
                     if x.typ = ints
                     then begin
                       x.typ := reals;
                       emit1(26,1)
                     end else if y.typ = ints
                              then begin
                                 y.typ := reals;
                                 emit1(26,0)
                              end;
                      if (x.typ=reals) and (y.typ=reals)
                      then case op of

                           eql:  emit(39);
                           neq:  emit(40);
                           lss:  emit(41);
                           leq:  emit(42);
                           gtr:  emit(43);
                           geq:  emit(44);

                           end
                      else error(35);
                   end;
                   x.typ := bools
                end
          end {expression};


      PROCEDURE assignment(lv,  ad: integer);

          var
             x, y: item;
             f: integer;         {tab[i]. obj in [variable,prozedure]}

          begin
             x.typ := tab[i].typ;
             x.ref1 := tab[i].ref1;
             if tab[i].normal
             then
                f := 0
             else
                f := 1;
             emit2(f, lv, ad);
             if sy in [lbrack,lparent,period] then
                selector([becomes, eql] + fsys, x);
             if sy = becomes
             then
                insymbol
             else begin
                error(51);
                if sy = eql then insymbol
             end;
             expression(fsys, y);
             if x.typ = y.typ
             then
                if x.typ in stantyps
                then
                   emit(38)
                else
                   if x.ref1 <> y.ref1
                   then
                      error(46)
                   else
                      if x.typ = arrays
                      then
                         emit1(23, atab[x.ref1].size)
                      else
                         emit1(23, btab[x.ref1].vsize)
             else if (x.typ=reals) and (y.typ=ints)
             then begin
                emit1(26,0);
                emit(38)
             end else if (x.typ<>notyp) and (y.typ<>notyp)
                      then error(46)
          end {assignment};


      PROCEDURE compoundstatement;

          begin
             insymbol;
             statement([semicolon, endsy] + fsys);
             while sy in [semicolon] + statbegsys do
                begin
                   if sy = semicolon
                   then
                      insymbol
                   else
                      error(14);
                   statement([semicolon, endsy] + fsys)
                end;
             if sy = endsy
             then
                insymbol
             else
                error(57)
          end {compoundstatement};


      PROCEDURE ifstatement;

          var
             x: item;
             lc1, lc2: integer;

          begin
             insymbol;
             expression(fsys + [thensy, dosy], x);
             if not (x.typ in [bools, notyp])   then
                error(17);
             lc1 := lc;
             emit(11);          {jmpc}
             if sy = thensy
             then
                insymbol
             else begin
                error(52);
                if sy = dosy
                then insymbol
             end;

             statement(fsys + [elsesy]);

             if sy = elsesy
             then
                begin
                   insymbol;
                   lc2 := lc;
                   emit(10);
                   code[lc1].y := lc;
                   statement(fsys);
                   code[lc2].y := lc
                end
             else
                code[lc1].y := lc
          end {ifstatement};

      PROCEDURE casestatement;
      var     x : item;
              i,j,k,lc1 : integer;
              casetab : array[1..csmax] of
                         packed record
                           val, lc: index
                         end;
              exittab : array[1..cmax] of integer;

      PROCEDURE caselabel;
      var     lab : conrec;
              k   : integer;
      begin
        constant(fsys+[comma,colon], lab);
        if lab.tp <> x.typ
        then error(47)
        else if i = csmax
             then fatal(6)
             else begin
                i := i+1;     k := 0;
                casetab[i].val := lab.i;
                casetab[i].lc  := lc;
                repeat
                  k := k+1
                until casetab[k].val = lab.i;
                if k < i then error(1);  { multiple definition  }
             end
      end { caselabel };

      PROCEDURE onecase;
      begin
        if sy in constbegsys
        then begin
          caselabel;
          while sy = comma do
          begin
            insymbol;  caselabel
          end;
          if sy = colon
          then insymbol else error(5);
          statement([semicolon,endsy] + fsys);
          j := j+1;
          exittab[j] := lc; emit(10)
       end
      end { onecase };

   begin  { casestatement  }
     insymbol;
     i := 0;  j := 0;
     expression(fsys+[ofsy,comma,colon],x);
     if not (x.typ in [ints,bools,chars,notyp])
     then error(23);
     lc1 := lc;  emit(12);  { jmx  }
     if sy = ofsy then insymbol else error(8);
     onecase;
     while sy = semicolon do
     begin
       insymbol;
       onecase
     end;
     code[lc1].y := lc;
     for k:=1 to i do
     begin
       emit1(13,casetab[k].val);
       emit1(13,casetab[k].lc)
     end;
     emit1(10,0);
     for k:=1 to j do code[exittab[k]].y := lc;
     if sy = endsy then insymbol else error(57)
  end  { casestatement  };

      PROCEDURE repeatstatement;

          var
             x: item;
             lc1: integer;

          begin
             lc1 := lc;
             insymbol;
             statement([semicolon, untilsy] +   fsys);
             while sy in [semicolon] + statbegsys do
                begin
                   if sy = semicolon
                   then
                      insymbol
                   else
                      error(14);
                   statement([semicolon, untilsy] + fsys)
                end;
             if sy = untilsy
             then
                begin
                   insymbol;
                   expression(fsys, x);
                   if not (x.typ in [bools, notyp]) then
                      error(17);
                   emit1(11, lc1)
                end
             else
                error(53)
          end {repeatstatement};


      PROCEDURE whilestatement;

          var
             x: item;
             lc1, lc2: integer;

          begin
             insymbol;
             lc1 := lc;
             expression(fsys + [dosy], x);
             if not (x.typ in [bools, notyp])   then
                error(17);
             lc2 := lc;
             emit(11);
             if sy = dosy
             then
                insymbol
             else
                error(54);
             statement(fsys);
             emit1(10, lc1);
             code[lc2].y := lc
          end {whilestatement};


      PROCEDURE forstatement;

          var
             cvt: types;
             x: item;
             i, f, lc1, lc2: integer;

          begin
             insymbol;
             if sy = ident
             then
                begin
                   i := loc(id);
                   insymbol;
                   if i = 0 then
                      cvt := ints
                   else if tab[i].obj = vvariable
                   then
                      begin
                         cvt := tab[i].typ;
                         if not tab[i].normal
                         then
                            error(37)
                         else
                            emit2(0, tab[i].lev, tab[i].adr);
                         if not (cvt in [notyp, ints, bools, chars])
                         then
                            error(18)
                      end
                   else
                      begin
                         error(37);
                         cvt := ints
                      end
                end
             else
                skip([becomes,tosy,downtosy,dosy] + fsys, 2);
             if sy = becomes
             then
                begin
                   insymbol;
                   expression([tosy,downtosy,dosy] + fsys, x);
                   if x.typ <> cvt then
                      error(19);
                end
             else
                skip([tosy,downtosy,dosy] + fsys, 51);
                f := 14;
             if sy in [tosy,downtosy]
             then
                begin
                   if sy = downtosy then f := 16;
                   insymbol;
                   expression([dosy] + fsys, x);
                   if x.typ <> cvt then
                      error(19)
                end
             else
                skip([dosy] + fsys, 55);
             lc1 := lc;
             emit(f);
             if sy = dosy
             then
                insymbol
             else
                error(54);
             lc2 := lc;
             statement(fsys);
             emit1(f+1, lc2);
             code[lc1].y := lc
          end {forstatement};


{**************************************************}
{*                                                *}
{*        S T A N D A R D  P R O C E D U R E      *}
{*                                                *}
{**************************************************}


      PROCEDURE standproc(n: integer);

          var
             i, f, fnb	: integer;
             x, y	: item;
	     first	: boolean;

	  begin
            case n of
              1, 2: begin { read }
                      if not iflag then
		      begin
                        error(20); iflag := true
                      end;
		      fnb := 1;
                      if sy = lparent then
                      begin
			first := true;
                        repeat
                          insymbol;
                          if sy <> ident then
                            error(2)
                          else
                          begin
                            i := loc(id);
                            insymbol;
			    if first and (tab[i].typ = texts) then
			      fnb := tab[i].ref1
			    else
			    begin
                              if i <> 0 then
                                if tab[i].obj <> vvariable then
                                  error(37)
                                else
                                begin
                                  x.typ := tab[i].typ;
                                  x.ref1 := tab[i].ref1;
                                  if tab[i].normal then
                                    f := 0
                                  else
                                    f := 1;
                                  emit2(f, tab[i].lev, tab[i].adr);
                                  if sy in [lbrack,lparent,period] then
                                    selector(fsys + [comma, rparent], x);
                                  if x.typ in [ints,reals,chars, notyp] then
                                    emit2(27, fnb, ord(x.typ))
                                  else
                                    error(41);
                              end;
			    end;
                          end;
                          test([comma, rparent],   fsys, 6);
			  first := false;
                        until sy <> comma;
			getrparent;
                      end;
                      if n = 2 then
			emit1(62, fnb);
                    end;
              3, 4: begin { write }
		      fnb := 2;
                      if sy = lparent then
                      begin
			first := true;
                        repeat
                          insymbol;
                          if first and (sy = ident) then
			  begin
			    i := loc(id);
			    if tab[i].typ = texts then
			    begin
			      fnb := tab[i].ref1;
		              insymbol;
			    end;
			  end;
			  if not first or (first and (fnb = 2)) then
			  begin
                            if sy = stringcon then
                            begin
                              emit1(24, sleng);
			      emit2(28, fnb, inum);
                              insymbol;
                            end
                            else
                            begin
                              expression(fsys+[comma,colon, rparent], x);
                              if not (x.typ in stantyps) then
                                error(41);
                              if sy = colon then
			      begin
                                insymbol;
                                expression(fsys+[comma,colon,rparent],y);
                                if y.typ <> ints then error(43);
                                if sy = colon then
			        begin
                                  if x.typ <> reals then error(42);
                                  insymbol;
                                  expression(fsys+[comma,rparent],y);
                                  if y.typ <> ints then error(43);
                                  emit1(37, fnb);
                                end
			        else
				  emit2(30, fnb, ord(x.typ));
                              end
			      else
			        emit2(29, fnb, ord(x.typ));
			    end;
			  end;
			  first := false;
                        until sy <> comma;
		        getrparent;
                      end;
                      if n = 4 then
			emit1(63, fnb);
                    end;


  5 : begin { iniranf }
        getlparent;
        getparam(false, true, ints, x);
        getrparent;
        emit(64);
      end;
  6 : begin { newseed }
        emit(65);
      end;
  7 : begin { unitmat }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, false, arrays, x);
        getrparent;
        emit(66);
      end;
  8 : begin { copymat }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, false, arrays, x);
        getparam(true, false, arrays, x);
        getrparent;
        emit(67);
      end;
  9 : begin { copyvec }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, false, arrays, x);
        getparam(true, false, arrays, x);
        getrparent;
        emit(68);
      end;
 10 : begin { addvec }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, false, arrays, x);
        getparam(true, false, arrays, x);
        getrparent;
        emit(69);
      end;
 11 : begin { subvec }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, false, arrays, x);
        getparam(true, false, arrays, x);
        getrparent;
        emit(70);
      end;
 12 : begin { addmat }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, false, arrays, x);
        getparam(true, false, arrays, x);
        getrparent;
        emit(71);
      end;
 13 : begin { submat }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, false, arrays, x);
        getparam(true, false, arrays, x);
        getrparent;
        emit(72);
      end;
 14 : begin { mullmat }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, false, arrays, x);
        getparam(true, false, arrays, x);
        getrparent;
        emit(73);
      end;
 15 : begin { lintrans }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, false, arrays, x);
        getparam(true, false, arrays, x);
        getrparent;
        emit(74);
      end;
 16 : begin { tpmat }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, false, arrays, x);
        getrparent;
        emit(75);
      end;
 17 : begin { graope }
        emit( 76);
      end;
 18 : begin { graph1 }
        getlparent;
        getparam(false, false, arrays, x);
        getparam(true, false, arrays, x);
        getparam(true, true, ints, x);
        getparam(true, true, ints, x);
        getparam(true, true, ints, x);
        getparam(true, true, reals, x);
        getrparent;
        emit( 77);
      end;
 19 : begin { gracls }
        emit( 78);
      end;
 20 : begin { graclx }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, true, ints, x);
        getparam(true, true, ints, x);
        getparam(true, true, ints, x);
        getrparent;
        emit( 79);
      end;
 21 : begin { graheds }
        getlparent;
        getparam(false, false, arrays, x);
        getparam(true, true, ints, x);
        getparam(true, true, reals, x);
        getrparent;
        emit( 80);
      end;
 22 : begin { gravwp }
        getlparent;
        getparam(false, true, reals, x);
        getparam(true, true, reals, x);
        getparam(true, true, reals, x);
        getparam(true, true, reals, x);
        getrparent;
        emit( 81);
      end;
 23 : begin { grawnd }
        getlparent;
        getparam(false, true, reals, x);
        getparam(true, true, reals, x);
        getparam(true, true, reals, x);
        getparam(true, true, reals, x);
        getrparent;
        emit( 82);
      end;
 24 : begin { gracoms }
        getlparent;
        getparam(false, true, reals, x);
        getparam(true, false, arrays, x);
        getparam(true, true, ints, x);
        getparam(true, true, ints, x);
        getparam(true, false, arrays, x);
        getparam(true, true, ints, x);
        getparam(true, true, ints, x);
        getrparent;
        emit( 83);
      end;
 25 : begin { graxi1 }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, true, ints, x);
        getparam(true, true, ints, x);
        getparam(true, true, reals, x);
        getrparent;
        emit( 84);
      end;
 26 : begin { grayi1 }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, true, ints, x);
        getparam(true, true, ints, x);
        getparam(true, true, reals, x);
        getrparent;
        emit( 85);
      end;
 27 : begin { graxi2 }
        getlparent;
        getparam(false, true, reals, x);
        getparam(true, true, ints, x);
        getparam(true, true, ints, x);
        getrparent;
        emit( 86);
      end;
 28 : begin { grayi2 }
        getlparent;
        getparam(false, true, reals, x);
        getparam(true, true, ints, x);
        getparam(true, true, ints, x);
        getrparent;
        emit(87);
      end;
 29 : begin { gratexs }
        getlparent;
        getparam(false, true, reals, x);
        getparam(true, true, reals, x);
        getparam(true, true, reals, x);
        getparam(true, false, arrays, x);
        getparam(true, true, ints, x);
        getparam(true, true, reals, x);
        getparam(true, false, reals, x);
        getparam(true, false, reals, x);
        getparam(true, true, ints, x);
        getrparent;
        emit( 88);
      end;
 30 : begin { gravec }
        getlparent;
        getparam(false, true, reals, x);
        getparam(true, true, reals, x);
        getparam(true, true, reals, x);
        getparam(true, true, reals, x);
        getrparent;
        emit(89);
      end;
 31 : begin { gratyp }
        getlparent;
        getparam(false, true, ints, x);
        getrparent;
        emit(90);
      end;
 32 : begin { Cell_DATwiss }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, true, ints, x);
        getparam(true, false, arrays, x);
        getparam(true, true, bools, x);
        getparam(true, true, bools, x);
        getparam(true, true, reals, x);
        getrparent;
        emit(91);
      end;
 33 : begin { gramar }
        getlparent;
        getparam(false, true, reals, x);
        getparam(true, true, reals, x);
        getparam(true, true, ints, x);
        getparam(true, true, reals, x);
        getrparent;
        emit(92);
      end;
 34 : begin { gdiag }
        getlparent;
        getparam(false, true, ints, x);
        getparam(false, true, reals, x);
        getparam(true, false, arrays, x);
        getparam(true, false, arrays, x);
        getparam(true, false, arrays, x);
        getparam(true, false, arrays, x);
        getparam(true, false, reals, x);
        getparam(true, false, reals, x);
        getrparent;
        emit(93);
      end;
 35 : begin { svbksb }
        getlparent;
        getparam(false, false, arrays, x);
        getparam(true, false, arrays, x);
        getparam(true, false, arrays, x);
        getparam(true, true, ints, x);
        getparam(true, true, ints, x);
        getparam(true, false, arrays, x);
        getparam(true, false, arrays, x);
        getrparent;
        emit(94);
      end;
 36 : begin { dsvdc }
        getlparent;
        getparam(false, false, arrays, x);
        getparam(true, true, ints, x);
        getparam(true, true, ints, x);
        getparam(true, false, arrays, x);
        getparam(true, false, arrays, x);
        getrparent;
        emit(95);
      end;
 37 : begin { break }
        getlparent;
        getparam(false, true, bools, x);
        getrparent;
        emit(96);
      end;
 38 : begin { getstr }
        getlparent;
        getparam(false, false, records, x);
	getstrcon(true);
        getrparent;
        emit1(97, inum);
      end;
 39 : begin { copystr }
        getlparent;
        getparam(false, false, records, x);
        getparam(true, false, records, x);
        getrparent;
        emit(98);
      end;
 40 : begin { concat }
        getlparent;
        getparam(false, false, records, x);
	getstrcon(true);
        getrparent;
        emit1(99, inum);
      end;
 41 : begin { getint }
        getlparent;
        getparam(false, false, records, x);
        getparam(true, true, ints, x);
        getparam(true, true, ints, x);
        getrparent;
        emit(100);
      end;
 42 : begin { getreal }
        getlparent;
        getparam(false, false, records, x);
        getparam(true, true, reals, x);
        getparam(true, true, ints, x);
        getparam(true, true, ints, x);
        getrparent;
        emit(101);
      end;
 43 : begin { getreale }
        getlparent;
        getparam(false, false, records, x);
        getparam(true, true, reals, x);
        getparam(true, true, ints, x);
        getparam(true, true, ints, x);
        getrparent;
        emit(102);
      end;
 44 : begin { writestr }
        getlparent;
        getparam(false, false, records, x);
        getparam(true, true, ints, x);
        getrparent;
        emit(103);
      end;
 45 : begin { graph0 }
        getlparent;
        getparam(false, false, arrays, x);
        getparam(true, false, arrays, x);
        getparam(true, true, ints, x);
        getrparent;
        emit(104);
      end;
 46 : begin { mulrmat }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, false, arrays, x);
        getparam(true, false, arrays, x);
        getrparent;
        emit(105);
      end;
 47 : begin { Cell_Pass }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, true, ints, x);
        getparam(true, false, arrays, x);
        getparam(true, false, ints, x);
        getrparent;
        emit(106);
      end;
 48 : begin { Cell_Pass_M }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, true, ints, x);
        getparam(true, false, arrays, x);
        getparam(true, false, arrays, x);
        getparam(true, false, ints, x);
        getrparent;
        emit(107);
      end;
 49 : begin { setrancut }
        getlparent;
        getparam(false, true, reals, x);
        getrparent;
        emit(108);
      end;
 50 : begin { Elem_getbend }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, true, ints, x);
        getparam(true, false, reals, x);
        getparam(true, false, reals, x);
        getparam(true, false, reals, x);
        getrparent;
        emit(109);
      end;
 51 : begin { Cell_Getcod }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, true, reals, x);
        getparam(true, true, reals, x);
        getparam(true, false, ints, x);
        getrparent;
        emit(110);
      end;
 52 : begin { Cell_GetABGN }
        getlparent;
        getparam(false, false, arrays, x);
        getparam(true, false, arrays, x);
        getparam(true, false, arrays, x);
        getparam(true, false, arrays, x);
        getparam(true, false, arrays, x);
        getrparent;
        emit(111);
      end;
 53 : begin { trace }
        getlparent;
        getparam(false, true, bools, x);
        getrparent;
        emit(112);
      end;
 54 : begin { Cell_MatTwiss }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, true, ints, x);
        getparam(true, false, arrays, x);
        getparam(true, true, bools, x);
        getparam(true, true, bools, x);
        getparam(true, true, reals, x);
        getrparent;
        emit(113);
      end;
 55 : begin { FFT }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, false, arrays, x);
        getparam(true, false, arrays, x);
        getrparent;
        emit(114);
      end;
 56 : begin { Ring_GetTwiss }
        getlparent;
        getparam(false, true, bools, x);
        getparam(true, true, reals, x);
        getrparent;
        emit(115);
      end;
 57 : begin { Ring_Fittune }
        getlparent;
        getparam(false, false, arrays, x);
        getparam(true, true, reals, x);
        getparam(true, false, arrays, x);
        getparam(true, false, arrays, x);
        getparam(true, false, arrays, x);
        getparam(true, true, reals, x);
        getparam(true, true, ints, x);
        getrparent;
        emit(116);
      end;
 58 : begin { Ring_Fitchrom }
        getlparent;
        getparam(false, false, arrays, x);
        getparam(true, true, reals, x);
        getparam(true, false, arrays, x);
        getparam(true, false, arrays, x);
        getparam(true, false, arrays, x);
        getparam(true, true, reals, x);
        getparam(true, true, ints, x);
        getrparent;
        emit(117);
      end;
 59 : begin { Ring_Fitdisp }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, true, reals, x);
        getparam(true, true, reals, x);
        getparam(true, true, ints, x);
        getparam(true, false, arrays, x);
        getparam(true, true, reals, x);
        getparam(true, true, ints, x);
        getrparent;
        emit(118);
      end;
 60 : begin { }
        emit(119);
      end;
 61 : begin { Cell_Concat }
        getlparent;
        getparam(false, true, reals, x);
        getrparent;
        emit(120);
      end;
 62 : begin { Cell_fPass }
        getlparent;
        getparam(false, false, arrays, x);
        getparam(true, false, ints, x);
        getrparent;
        emit(121);
      end;
 63 : begin { getelem }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, false, records, x);
        getrparent;
        emit(122);
      end;
 64 : begin { Cell_DApass }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, true, ints, x);
        getparam(true, false, arrays, x);
        getparam(true, false, ints, x);
        getrparent;
        emit(123);
      end;
 65 : begin { initbump }
        getlparent;
        getparam(false, false, arrays, x);
        getparam(true, false, arrays, x);
        getparam(true, false, arrays, x);
        getparam(true, true, reals, x);
        getparam(true, true, reals, x);
        getrparent;
        emit(124);
      end;
 66 : begin { execbump }
        getlparent;
        getparam(false, true, reals, x);
        getparam(true, true, ints, x);
	getrparent;
        emit(125);
      end;
 67 : begin { getglobv_ }
        getlparent;
        getparam(false, false, records, x);
	getrparent;
        emit(126);
      end;
 68 : begin { putglobv_ }
        getlparent;
        getparam(false, false, records, x);
	getrparent;
        emit(127);
      end;
 69 : begin { grabeg }
        emit(128);
      end;
 70 : begin { graend }
        emit(129);
      end;
 71 : begin { grafpl }
        getlparent;
        getparam(false, false, arrays, x);
        getparam(true, false, arrays, x);
        getparam(true, true, ints, x);
        getrparent;
        emit(130);
      end;
 72 : begin { gracoi }
        getlparent;
        getparam(false, true, ints, x);
        getrparent;
        emit(131);
      end;
 73 : begin { grapat }
        getlparent;
        getparam(false, true, ints, x);
        getrparent;
        emit(132);
      end;
 74 : begin { gravcs }
        getlparent;
        getparam(false, true, reals, x);
        getparam(true, true, reals, x);
        getparam(true, true, reals, x);
        getparam(true, true, reals, x);
        getrparent;
        emit(133);
      end;
 75 : begin { gravcf }
        getlparent;
        getparam(false, true, ints, x);
        getrparent;
        emit(134);
      end;
 76 : begin { }
        emit(135);
      end;
 77 : begin { }
        emit(136);
      end;
 78 : begin { }
        emit(137);
      end;
 79 : begin { reset_ }
        getlparent;
        getparam(false, false, texts, x);
	fnb := x.ref1;
        getparam(true, false, arrays, x);
        getrparent;
        emit1(138, fnb);
      end;
 80 : begin { rewrite_ }
        getlparent;
        getparam(false, false, texts, x);
	fnb := x.ref1;
        getparam(true, false, arrays, x);
        getrparent;
        emit1(139, fnb);
      end;
 81 : begin { close }
        getlparent;
        getparam(false, false, texts, x);
	fnb := x.ref1;
        getrparent;
        emit1(140, fnb);
      end;
 82 : begin { }
        emit(141);
      end;
 83 : begin { }
        emit(142);
      end;
 84 : begin { }
        emit(143);
      end;
 85 : begin { }
        emit(144);
      end;
 86 : begin { Mpole_setds }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, true, ints, x);
        getrparent;
        emit(145);
      end;
 87 : begin { }
        emit(146);
      end;
 88 : begin { }
        emit(147);
      end;
 89 : begin { }
        emit(148);
      end;
 90 : begin { }
        emit(149);
      end;
 91 : begin { Mpole_setdt }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, true, ints, x);
        getrparent;
        emit(150);
      end;
 92 : begin { }
        emit(151);
      end;
 93 : begin { }
        emit(152);
      end;
 94 : begin { }
        getrparent;
        emit(153);
      end;
 95 : begin { Mpole_setpb }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, true, ints, x);
        getparam(true, true, ints, x);
        getrparent;
        emit(154);
      end;
 96 : begin { Mpole_GetPmeth }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, true, ints, x);
        getparam(true, false, ints, x);
        getparam(true, false, ints, x);
        getparam(true, false, reals, x);
        getparam(true, false, reals, x);
        getparam(true, false, reals, x);
        getparam(true, false, reals, x);
        getrparent;
        emit(155);
      end;
 97 : begin { Cav_Set }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, true, ints, x);
        getparam(true, true, reals, x);
        getparam(true, true, reals, x);
        getparam(true, true, ints, x);
        getrparent;
        emit(156);
      end;
 98 : begin { Cav_Get }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, true, ints, x);
        getparam(true, false, reals, x);
        getparam(true, false, reals, x);
        getparam(true, false, ints, x);
        getrparent;
        emit(157);
      end;
 99 : begin { GetDBN }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, true, ints, x);
        getparam(true, false, arrays, x);
        getrparent;
        emit(158);
      end;
100 : begin { }
        emit(159);
      end;
101 : begin { }
        emit(160);
      end;
102 : begin { }
        emit(161);
      end;
103 : begin { }
        emit(162);
      end;
104 : begin { }
        emit(163);
      end;
105 : begin { }
        emit(164);
      end;
106 : begin { }
        emit(165);
      end;
107 : begin { }
        emit(166);
      end;
108 : begin { }
        emit(167);
      end;
109 : begin { }
        emit(168);
      end;
110 : begin { }
        emit(169);
      end;
111 : begin { amoeba }
	getlparent;
	getparam(false, false, arrays, x);
	getparam(true, false, arrays, x);
	getparam(true, true, ints, x);
	getparam(true, false, reals, x);
	getparam(true, true, reals, x);
	getparam(true, true, ints, x);
	getparam(true, false, ints, x);
	getrparent;
        emit(170);
      end;
112 : begin { Cell_SetdP }
	getlparent;
	getparam(false, true, reals, x);
	getrparent;
        emit(171);
      end;
113 : begin { putelem }
        getlparent;
        getparam(false, true, ints, x);
        getparam(true, false, records, x);
        getrparent;
        emit(172);
      end;


            end {case}
	  end {standproc};



      begin {statement}
          if sy in statbegsys + [ident]
          then
             case sy of
                ident:
                   begin
                      i := loc(id);
                      insymbol;
                      if i <> 0
                      then
                         case tab[i].obj of
                   konstant, type1:
                               error(45);
                         vvariable:
                               assignment(tab[i].lev,   tab[i].adr);
                         prozedure:
                               if tab[i].lev <> 0
                               then
                                  call(fsys, i)
                               else
                                  standproc(tab[i].adr);
                          funktion:
                               if tab[i].ref1 = display[level]
                               then
                                  assignment(tab[i].lev + 1, 0)
                               else
                                  error(45)
                         end { case }
                   end;
                beginsy: compoundstatement;
                   ifsy: ifstatement;
                 casesy: casestatement;
                whilesy: whilestatement;
               repeatsy: repeatstatement;
                  forsy: forstatement;

                    end; { case }
           test(fsys, [], 14)
      end {statement};

   begin {block}
      dx := 5;
      prt := t;
      if level  > lmax then fatal(5);
      test([lparent, colon, semicolon], fsys, 140);
      enterblock;
      display[level]    := b;
      prb := b;
      tab[prt].typ := notyp;
      tab[prt].ref1 := prb;
      if (sy =  lparent) and (level > 1) then
          parameterlist;
      btab[prb].lastpar := t;
      btab[prb].psize := dx;

      if isfun
      then
          if sy = colon
          then
             begin
                insymbol;       {FUNCTION type}
                if sy = ident
                then
                   begin
                      x := loc(id);
                      insymbol;
                      if x <> 0 then
                         if tab[x].obj <> type1
                         then
                            error(29)
                         else
                            if tab[x].typ in stantyps
                            then
                               tab[prt].typ := tab[x].typ
                            else
                               error(15)
                   end
                else
                   skip([semicolon] +   fsys, 213)
             end
          else error(5);{ end of function..}

      if sy = semicolon
        then insymbol
        else error(14);

      if (level = 1) and (finame <> '') then constdec0;

      repeat
        repeat
          if sy = constsy then constdec;
          if sy = typesy  then typedeclaration;
          if sy = varsy   then vardeclaration;
          btab[prb].vsize := dx;
          while sy in [procsy, funcsy] do procdeclaration;
        until not (sy in [constsy,typesy,varsy,procsy,funcsy]);
        test([beginsy], blockbegsys   + statbegsys, 56);
      until sy  in statbegsys;

      tab[prt].adr := lc;
      insymbol;
      statement([semicolon, endsy] +    fsys);
      while sy  in [semicolon] + statbegsys do
          begin
             if sy = semicolon
             then
                insymbol
             else
                error(14);
             statement([semicolon, endsy] + fsys)
          end;
      if sy = endsy
      then
          insymbol
      else
          error(57);
      test(fsys + [period], [], 6);
   end  {block};

PROCEDURE interpret;

{ Pascal-s interpreter }

var	ir		: order;
	pc		: integer;
	t		: integer;
	b		: integer;
	ocnt		: integer;
	h1,h2,h3,h4	: integer;
	lncnt, chrcnt	: integer;

       ps: (run,fin,caschk,divchk,inxchk,stkchk,linchk,lngchk,redchk);

   fld     : array[1..4] of integer;
   display : array[0..lmax] of integer;
   s       : array[1..stacksize] of
          record
            case cn: types of
              ints:			(i : integer);
              reals:			(r : double);
              bools:			(b : boolean);
              chars:			(c : char);
	      notyp, arrays, records:	();
          end;


procedure pmdump(prtarr : boolean);

{ Post mortem dump  }

const	colnb=2;

var	i, j, k, imax, index, h1, h2, h3, blkcnt	: integer;

begin
  h1 := b;  blkcnt := 10;
  repeat
    writeln;
    blkcnt := blkcnt - 1;
    if blkcnt = 0 then h1 := 0; h2 := s[h1+4].i;
    if h1 <> 0 then 
    writeln(' ', tab[h2].name, ' called at ', linerror[s[h1+1].i-1]:1 );
    h2 := btab[tab[h2].ref1].last;
    while h2 <> 0 do
      with tab[h2] do
      begin
        if obj = vvariable then
	begin
	  if typ in stantyps+[arrays] then
            if normal then h3 := h1+adr else h3:=s[h1+adr].i;

	  if typ in stantyps then
	  begin
            write('    ', name, ' = ');
            case typ of
              ints: writeln(' ', s[h3].i:1);
              reals: writeln(s[h3].r:13);
              bools: writeln(s[h3].b);
	      chars: writeln(s[h3].c);
	      otherwise writeln;
            end;
          end
	  else if prtarr and (typ = arrays) then
	  begin
	    index := ref1;
	    imax := ord(atab[ref1].high-atab[ref1].low) + 1;
	    while atab[index].elref <> 0 do
	    begin
	      index := atab[index].elref;
	      imax := imax*(ord(atab[index].high-atab[index].low)+1);
	    end;
	    for i:=1 to imax do
	    begin
	        if (i-1) mod colnb = 0 then
		  write('  ')
	        else
		 write(', ');
		j := 1;
		while name[j] <> ' ' do
		begin
		  write(name[j]);
		  j := succ(j);
		end;
		write('[');
	        k := ref1; j := i;
		repeat
		  if k = ref1 then
		    write((j-1) div atab[k].elsize+atab[k].low:1)
		  else
		    write(',', (j-1) div atab[k].elsize+atab[k].low:1);
		  j := j - ((j-1) div atab[k].elsize)*atab[k].elsize;
	          k := atab[k].elref;
	        until k = 0 ;
	        write('] = ');

	      case atab[index].eltyp of
		ints:  write(s[h3+i-1].i:1);
		reals: write(s[h3+i-1].r:13);
		bools: write(s[h3+i-1].b);
		chars: write(s[h3+i-1].c);
	      end;
	      if ((i-1) mod colnb = colnb-1) or (i = imax) then writeln;
	    end;
	  end;

	end;
        h2 := link;
      end;
      h1 := s[h1+3].i
  until h1 < 0;
end;


PROCEDURE dump;

var   p,h3 : integer;

begin
  h3 := tab[h2].lev;
  writeln(psout); writeln(psout);
  writeln(psout,'        calling ',tab[h2].name);
  writeln(psout,'          level ',h3:4);
  writeln(psout,'  start of code ',pc:4);
  writeln(psout); writeln(psout);
  writeln(psout,'  contents of display '); writeln(psout);

  for p:=h3 downto 0 do writeln(psout,p:6,display[p]:12);
  writeln(psout); writeln(psout);
  writeln(psout,' top of stack   ',t:4,' frame base ':14,b:4);
  writeln(psout); writeln(psout);
  writeln(psout,'stack contents':22); writeln(psout);

  for p:=t downto 1 do writeln(psout,p:14,s[p].i:8);

  writeln(psout,'< = = = >':22)
end; { dump }


PROCEDURE ExecPcode;

   var 	k	: integer;

	job, mn, m, n, info		: integer;
	snglx, sngly			: real;
	ch				: char;
	simpvec1, simpvec2		: simpvect;
	vec1, vec2			: vector;
	mat1, mat2, mat3, mat4		: matrix;
	map				: DAmap;
	ivec21				: ivector2;
	vec21, vec22, vec23, vec24	: vector2;
	px, py				: graphvect;
	bx, by				: bigvect;
	strb1, strb2			: strbuf;
	str1, str2			: tstring;
	DBname				: DBnametype;

	xf, xd				: fitvect;

	a1, u1, v1			: svdmat;
	b1, w1, x1, work, e1		: svdarray;


	procedure popstrcon(t : integer; var str : tstring);

	var	i	: integer;

	begin
	  h1 := s[t].i; h2 := ir.y;
          str.len := h1;
          for i:=str.len+1 to strlmax do
            str.str[i] := ' ';
	  if str.len = 1 then
	    { Character constant }
	    str.str[1] := chr(h2)
	  else if str.len > 1 then
	    { String constant }
	    repeat 
	      str.str[str.len-h1+1] := stab[h2];
	      h1 := h1 - 1; h2 := h2 + 1;
	    until h1 = 0 ;
	end;


	procedure pushstrbuf(t, len : integer;
		   var str : packed array [low..high : integer] of char);

	var	j	: integer;

	begin
	  for j:=low to low-1+len do
	    s[s[t].i+j-low].i := ord(str[j]);
	end;


	procedure popstrbuf(t, len : integer;
		   var str : packed array [low..high : integer] of char);

	var	j	: integer;

	begin
	  for j:=low to low-1+len do
	    str[j] := chr(s[s[t].i+j-low].i); 
	end;


	procedure pushstr(t : integer; var str : tstring);

	var	j, len	: integer;

	begin
	  len := strlen_(str);
	  { To make sure assignment of an empty string }
	  s[s[t].i].i := len;
	  for j:=1 to len do
	    s[s[t].i+j].i := ord(str.str[j]);
	  for j:=len+1 to strlmax do
	    s[s[t].i+j].i := ord(' ');
	end;


	procedure popstr(t : integer; var str : tstring);

	var	j, len	: integer;

	begin
	  { To make sure assignment of an empty string }
	  str.len := s[s[t].i].i;
	  len := strlen_(str);
	  for j:=1 to len do
	    str.str[j] := chr(s[s[t].i+j].i); 
	  for j:=len+1 to strlmax do
	    str.str[j] := ' ';
	end;


	procedure pushDBname(t : integer; var DBname : DBnametype);

	var	j	: integer;

	begin
	  for j:=1 to DBnamelen do
	    s[s[t].i+j-1].i := ord(DBname[j]); 
	end;


	procedure pushivec(t, n : integer;
			   var ivec : array [low..high : integer] of integer);

	var	j	: integer;

	begin
	  for j:=low to n-low+1 do
	    s[s[t].i+j-low].i := ivec[j];
	end;


	procedure popivec(t, n : integer;
			 var ivec : array [low..high : integer] of integer);

	var	j	: integer;

	begin
	  for j:=low to n-low+1 do
	    ivec[j] := s[s[t].i+j-low].i; 
	end;


	procedure pushvec(t, n : integer;
			  var vec : array [low..high : integer] of double);

	var	j	: integer;

	begin
	  for j:=low to n-low+1 do
	    s[s[t].i+j-low].r := vec[j];
	end;


	procedure popvec(t, n : integer;
			 var vec : array [low..high : integer] of double);

	var	j	: integer;

	begin
	  for j:=low to n-low+1 do
	    vec[j] := s[s[t].i+j-low].r; 
	end;


	procedure pushmat(t, n : integer; var mat : matrix);

	var	i, j	: integer;

	begin
	  for i:=1 to n do
	    for j:=1 to n do
	      s[s[t].i+(i-1)*matdim+j-1].r := mat[i, j];
	end;


	procedure popmat(t, n : integer; var mat : matrix);

	var	i, j	: integer;

	begin
	  for i:=1 to n do
	    for j:=1 to n do
	      mat[i, j] := s[s[t].i+(i-1)*matdim+j-1].r;
	end;


	procedure pushsvdmat(t : integer; tp : boolean; m, n : integer;
			     var u : svdmat);

	var	i, j	: integer;

	begin
	  for i:=1 to m do
	    for j:=1 to n do
	      if tp then
		{ Fortran call, pass transposed matrices }
		s[s[t].i+(i-1)*mnp+j-1].r := u[j, i]
	      else
		s[s[t].i+(i-1)*mnp+j-1].r := u[i, j];
	end;


	procedure popsvdmat(t : integer; tp : boolean; m, n : integer;
			    var u : svdmat);

	var	i, j	: integer;

	begin
	  for i:=1 to m do
	    for j:=1 to n do
	      if tp then
		{ Fortran call, pass transposed matrices }
		u[j, i] := s[s[t].i+(i-1)*mnp+j-1].r
	      else
		u[i, j] := s[s[t].i+(i-1)*mnp+j-1].r;
	end;


	procedure pushDAmap(t : integer; var map : DAmap);

	var	i, j	: integer;

	begin
	  for i:=1 to nvmax do
	    for j:=0 to nvmax do
	      s[s[t].i+(i-1)*(nvmax+1)+j].r := map[i, j]; 
	end;


	procedure popDAmap(t : integer; var map : DAmap);

	var	i, j	: integer;

	begin
	  for i:=1 to nvmax do
	    for j:=0 to nvmax do
	      map[i, j] := s[s[t].i+(i-1)*(nvmax+1)+j].r;
	end;


	procedure pushgvec(t, n : integer; var vec : graphvect);

	var	j	: integer;

	begin
	  for j:=1 to n do
	    s[s[t].i+j-1].r := vec[j];
	end;


	procedure popgvec(t, n : integer; var vec : graphvect);

	var	j	: integer;

	begin
	  for j:=1 to n do
	    vec[j] := sngl(s[s[t].i+j-1].r); 
	end;


	procedure pushbvec(t, n : integer; var vec : bigvect);

	var	j	: integer;

	begin
	  for j:=1 to n do
	    s[s[t].i+j-1].r := vec[j];
	end;


	procedure popbvec(t, n : integer; var vec : bigvect);

	var	j	: integer;

	begin
	  for j:=1 to n do
	    vec[j] := sngl(s[s[t].i+j-1].r); 
	end;


	procedure pushglobval(t : integer);

	var	i, j, ind	: integer;

	begin
	  s[s[t].i].r := globval.dPcommon;
	  s[s[t].i+1].r := globval.dPparticle;
	  s[s[t].i+2].r := globval.maxampl[1];
	  s[s[t].i+3].r := globval.maxampl[2];

	  s[s[t].i+4].r := globval.Totaltune[1];
	  s[s[t].i+5].r := globval.Totaltune[2];
	  s[s[t].i+6].r := globval.Omega;
	  s[s[t].i+7].r := globval.alphac;
	  s[s[t].i+8].r := globval.chrom[1];
	  s[s[t].i+9].r := globval.chrom[2];

	  s[s[t].i+10].r := globval.Energy;

	  s[s[t].i+11].i := globval.Cell_nloc;
	  s[s[t].i+12].i := globval.Elem_nFam;

	  s[s[t].i+13].i := globval.CODimax;
	  s[s[t].i+14].r := globval.CODeps;

	  for i:=1 to matdim do
	    s[s[t].i+14+i].r := globval.CODvect[i];

	  ind := matdim + 15;
	  s[s[t].i+ind].i := globval.bpm;

	  for i:=1 to matdim do
	    for j:=1 to matdim do
	    begin
	      s[s[t].i+ind+matdim*(i-1)+j].r := globval.oneturnmat[i, j];
	      s[s[t].i+ind+sqr(matdim)+matdim*(i-1)+j].r := globval.Ascr[i, j];
	      s[s[t].i+ind+2*sqr(matdim)+matdim*(i-1)+j].r := globval.Ascrinv[i, j];
	    end;
          ind := ind + 3*sqr(matdim) + 1;
	  s[s[t].i+ind].b := globval.MatMeth;
	  s[s[t].i+ind+1].b := globval.Cavity_on;
	  s[s[t].i+ind+2].b := globval.radiation;
	  s[s[t].i+ind+3].b := globval.emittance;
	  s[s[t].i+ind+4].r := globval.dE;
	  ind := ind + 4;
	  for i:=1 to 3 do
	    s[s[t].i+ind+i].r := globval.rad[i];
	  ind := ind + 3;
	  for i:=1 to 3 do
	    s[s[t].i+ind+i].r := globval.qfluct[i];
	  ind := ind + 4;
	  s[s[t].i+ind].b := globval.pathlength;
	end;


	procedure popglobval(t : integer);

	var	i, ind	: integer;

	begin
	  globval.dPcommon := s[s[t].i].r;
	  globval.dPparticle := s[s[t].i+1].r;
	  globval.maxampl[1] := s[s[t].i+2].r;
	  globval.maxampl[2] := s[s[t].i+3].r;

{	  globval.Energy := s[s[t].i+10].r;}

	  globval.CODimax := s[s[t].i+13].i;
	  globval.CODeps := s[s[t].i+14].r;

	  for i:=1 to matdim do
	    globval.CODvect[i] := s[s[t].i+14+i].r;

	  ind := matdim + 15;
	  globval.bpm := s[s[t].i+ind].i;
	  ind := ind + 3*sqr(matdim) + 1;
          globval.MatMeth := s[s[t].i+ind].b;
          globval.Cavity_on := s[s[t].i+ind+1].b;
	  globval.radiation := s[s[t].i+ind+2].b;
	  globval.emittance := s[s[t].i+ind+3].b;
	  globval.dE := s[s[t].i+ind+4].r;
	  ind := ind + 4;
	  for i:=1 to 3 do
	    globval.rad[i] := s[s[t].i+ind+i].r;
	  ind := ind + 3;
	  for i:=1 to 3 do
	    globval.qfluct[i] := s[s[t].i+ind+i].r;
	  ind := ind + 4;
	  globval.pathlength := s[s[t].i+ind].b;
	end;


	procedure pushelem(t, i : integer);

	var	j, k, ind	: integer;

	begin
	  { Can not use with due to Cell.s }
	  s[s[t].i+0].i := Cell[i].Fnum;  s[s[t].i+1].i := Cell[i].Knum;
	  s[s[t].i+2].r := Cell[i].S;
	  s[s[t].i+3].r := Cell[i].dS[1]; s[s[t].i+4].r := Cell[i].dS[2];
	  s[s[t].i+5].r := Cell[i].dT[1]; s[s[t].i+6].r := Cell[i].dT[2];

          for j:=1 to Namelength do
            s[s[t].i+6+j].i := ord(Cell[i].Elem.PName[j]);
	  ind := s[t].i + Namelength + 7;
	  s[ind+0].r := Cell[i].Elem.PL; s[ind+1].i := ord(Cell[i].Elem.Pkind);
	  if Cell[i].Elem.Pkind = Mpole then
	  begin
	    with Cell[i].Elem.M^ do
	    begin
	      s[ind+2].i := Pmethod;
	      s[ind+3].i := PN;
	      s[ind+4].r := PdSsys[1]; s[ind+5].r := PdSsys[2];
	      s[ind+6].r := PdSrms[1]; s[ind+7].r := PdSrms[2];
	      s[ind+8].r := PdSrnd[1]; s[ind+9].r := PdSrnd[2];

	      s[ind+10].r := PdTpar; s[ind+11].r := PdTsys;
	      s[ind+12].r := PdTrms; s[ind+13].r := PdTrnd;

	      ind := ind + 14 + HOMmax;
	      for j:=-HOMmax to HOMmax do
	        s[ind+j].r := PBpar[j];
	      ind := ind + 2*HOMmax + 1;
	      for j:=-HOMmax to HOMmax do
	        s[ind+j].r := PBsys[j];
	      ind := ind + 2*HOMmax + 1;
	      for j:=-HOMmax to HOMmax do
	        s[ind+j].r := PBrms[j];
	      ind := ind + 2*HOMmax + 1;
	      for j:=-HOMmax to HOMmax do
	        s[ind+j].r := PBrnd[j];
	      ind := ind + 2*HOMmax + 1;
	      for j:=-HOMmax to HOMmax do
	        s[ind+j].r := PB[j];

	      ind := ind + HOMmax + 1;

	      s[ind+0].i := Porder; s[ind+1].i := ord(Pthick);
	      if Pthick = Thick then
	      begin
	        s[ind+2].r := PTx1; s[ind+3].r := PTx2; s[ind+4].r := Pgap;
	        s[ind+5].r := Pirho;
	        s[ind+6].r := Pc0; s[ind+7].r := Pc1; s[ind+8].r := Ps1;
		ind := ind + 8;
		  for j:=1 to matdim do
		    for k:=1 to matdim do
		    begin
		      s[ind+matdim*(j-1)+k].r := AU55[j, k];
		      s[ind+sqr(matdim)+matdim*(j-1)+k].r := AD55[j, k];
		    end;
	      end
	      else
		ind := ind + 8;
	    end;
	    ind := ind + 2*sqr(matdim) + 1;
	  end
	  else
	    ind := ind + 14 + 10*HOMmax + 5 + 8 + 2*sqr(matdim) + 1;

	  s[ind+0].r := Cell[i].Nu[1];    s[ind+1].r := Cell[i].Nu[2];
	  s[ind+2].r := Cell[i].Alpha[1]; s[ind+3].r := Cell[i].Alpha[2];
	  s[ind+4].r := Cell[i].Beta[1];  s[ind+5].r := Cell[i].Beta[2];
	  s[ind+6].r := Cell[i].Eta[1];   s[ind+7].r := Cell[i].Eta[2];
	  s[ind+8].r := Cell[i].Etap[1];  s[ind+9].r := Cell[i].Etap[2];
	  ind := ind + 9;
	  for j:=1 to matdim do
	    s[ind+j].r := Cell[i].BeamPos[j];
	end;


	procedure popelem(t, i : integer);

	var	j, k, ind	: integer;

	begin
	  { Can not use with due to Cell.s }
	  Cell[i].Fnum := s[s[t].i+0].i;  Cell[i].Knum := s[s[t].i+1].i;
	  Cell[i].S := s[s[t].i+2].r;
	  Cell[i].dS[1] := s[s[t].i+3].r; Cell[i].dS[2] := s[s[t].i+4].r;
	  Cell[i].dT[1] := s[s[t].i+5].r; Cell[i].dT[2] := s[s[t].i+6].r;

          for j:=1 to Namelength do
            Cell[i].Elem.PName[j] := chr(s[s[t].i+6+j].i);
	  ind := s[t].i + Namelength + 7;
	  Cell[i].Elem.PL := s[ind+0].r;
	  case s[ind+1].i of
	    0: Cell[i].Elem.Pkind := drift;
	    1: Cell[i].Elem.Pkind := Wigl;
	    2: Cell[i].Elem.Pkind := Mpole;
	    3: Cell[i].Elem.Pkind := Cavity;
	    4: Cell[i].Elem.Pkind := marker;
	    5: Cell[i].Elem.Pkind := undef;
	    otherwise writeln('** undefined type');
	  end;
	  if Cell[i].Elem.Pkind = Mpole then
	  begin
	    with Cell[i].Elem.M^ do
	    begin
	      Pmethod := s[ind+2].i;
	      PN := s[ind+3].i;
	      PdSsys[1] := s[ind+4].r; PdSsys[2] := s[ind+5].r;
	      PdSrms[1] := s[ind+6].r; PdSrms[2] := s[ind+7].r;
	      PdSrnd[1] := s[ind+8].r; PdSrnd[2] := s[ind+9].r;

	      PdTpar := s[ind+10].r; PdTsys := s[ind+11].r;
	      PdTrms := s[ind+12].r; PdTrnd := s[ind+13].r;

	      ind := ind + 14 + HOMmax;
	      for j:=-HOMmax to HOMmax do
	        PBpar[j] := s[ind+j].r;
	      ind := ind + 2*HOMmax + 1;
	      for j:=-HOMmax to HOMmax do
	        PBsys[j] := s[ind+j].r;
	      ind := ind + 2*HOMmax + 1;
	      for j:=-HOMmax to HOMmax do
	        PBrms[j] := s[ind+j].r;
	      ind := ind + 2*HOMmax + 1;
	      for j:=-HOMmax to HOMmax do
	        PBrnd[j] := s[ind+j].r;
	      ind := ind + 2*HOMmax + 1;
	      for j:=-HOMmax to HOMmax do
	        PB[j] := s[ind+j].r;

	      ind := ind + HOMmax + 1;
	      Porder := s[ind+0].i;
	      case s[ind+1].i of
		0: Pthick := thick;
		1: Pthick := thin;
		otherwise writeln('** undefined type');
	      end;
	      if Pthick = Thick then
	      begin
	        PTx1 := s[ind+2].r; PTx2 := s[ind+3].r; Pgap := s[ind+4].r;
	        Pirho := s[ind+5].r;
	        Pc0 := s[ind+6].r; Pc1 := s[ind+7].r; Ps1 := s[ind+8].r;
		ind := ind + 8;
		  for j:=1 to matdim do
		    for k:=1 to matdim do
		    begin
		      AU55[j, k] := s[ind+matdim*(j-1)+k].r;
		      AD55[j, k] := s[ind+sqr(matdim)+matdim*(j-1)+k].r;
		    end;
	      end
	      else
		ind := ind + 8;
	    end;
	    ind := ind + 2*sqr(matdim) + 1;
	  end
	  else
	    ind := ind + 14 + 10*HOMmax + 5 + 8 + 2*sqr(matdim) + 1;

	  Cell[i].Nu[1] := s[ind+0].r;    Cell[i].Nu[2] := s[ind+1].r;
	  Cell[i].Alpha[1] := s[ind+2].r; Cell[i].Alpha[2] := s[ind+3].r;
	  Cell[i].Beta[1] := s[ind+4].r;  Cell[i].Beta[2] := s[ind+5].r;
	  Cell[i].Eta[1] := s[ind+6].r;   Cell[i].Eta[2] := s[ind+7].r;
	  Cell[i].Etap[1] := s[ind+8].r;  Cell[i].Etap[2] := s[ind+9].r;
	  ind := ind + 9;
	  for j:=1 to matdim do
	    Cell[i].BeamPos[j] := s[ind+j].r;
	end;


        { Local implementation of some of the routines for string handling
          with modified parameters }

	procedure getstr(var str2, str1 : tstring);

	{ Put str1 in str2 }

	var	i, len1	: integer;

	begin
	  len1 := strlen_(str1);
	  for i:=1 to len1 do
	    str2.str[i] := str1.str[i];
	  str2.len := len1;
	  for i:=len1+1 to strlmax do
	    str2.str[i] := ' ';
	end;


	procedure concat(var str2, str1 : tstring);

	{ Concatenate str1 to str2 }

	var	i, len1, len2	: integer;

	begin
	  len1 := strlen_(str1);
	  len2 := strlen_(str2);
	  for i:=1 to len1 do
	    str2.str[len2+i] := str1.str[i];
	  str2.len := len1 + len2;
	end;


	function strind(var object, pattern : tstring) : integer;

	{ Locates the first occurrence of a pattern tstring within an object
	  tstring. Index is set to:

		location of the first occurrence
		0 if pattern is not found
		1 if both pattern and object are empty tstrings
		0 if pattern, but not object, is an empty tstring

	}

	label	10;

	var	i, j, olen, plen, ind	: integer;

	begin
	  olen := strlen_(object);
	  plen := strlen_(pattern);
	  if plen = 0 then
	    if olen = 0 then
	      ind := 1
	    else
	      ind := 0
	  else
	  begin
	    ind := 0;
	    i := 1;
	    while i <= olen do
	    begin
	      if pattern.str[1] = object.str[i] then
	      begin
	        ind := i;
	        for j:=1 to plen do
		  if i+j-1 <= olen then
		  begin
		    if pattern.str[j] <> object.str[i+j-1] then ind := 0;
		  end
		  else
		    ind := 0;
	        if ind <> 0 then
		  goto 10
	        else
		  i := i + plen - 1;
	      end;
	      i := succ(i);
	    end;
	  end;
10:	  strind := ind;
	end;


   begin
      case ir.f of
             0:
                begin {load address}
                   t := t + 1;
                   if t > stacksize
                   then
                      ps := stkchk
                   else
                      s[t].i :=   display[ir.x] + ir.y
                end;
             1:
                begin {load value}
                   t := t + 1;
                   if t > stacksize
                   then
                      ps := stkchk
                   else
                      s[t] :=   s[display[ir.x] + ir.y]
                end;
             2:
                begin {load indirect}
                   t := t + 1;
                   if t > stacksize
                   then
                      ps := stkchk
                   else
                      s[t] :=   s[s[display[ir.x] + ir.y].i]
                end;
             3:
                begin {update display}
                   h1 := ir.y;
                   h2 := ir.x;
                   h3 := b;
                   repeat
                      display[h1] := h3;
                      h1 := h1 - 1;
                      h3 := s[h3 + 2].i
                   until h1 = h2
                end;

  4 : begin
	{ Has to be put here since inum has to be passed }
	popstrcon(t, str1);
	popstr(t-1, str2);
        s[t-1].i := strind(str2, str1);
        t := t - 1;
      end;

             8:
                case ir.y of
                    0: s[t].i := abs(s[t].i);
                    1: s[t].r := abs(s[t].r);
                    2: s[t].i := sqr(s[t].i);
                    3: s[t].r := sqr(s[t].r);
                    4: s[t].b := odd(s[t].i);
                    5: s[t].i := ord(chr(s[t].i));

		    { Modified }
                    6: s[t].i := ord(s[t].i);
                    7: s[t].i := succ(s[t].i);
                    8: s[t].i := pred(s[t].i);

                    9: s[t].i := round(s[t].r);
                   10: s[t].i := trunc(s[t].r);
                   11: s[t].r := sin(s[t].r);
                   12: s[t].r := cos(s[t].r);
                   13: s[t].r := exp(s[t].r);
                   14: s[t].r := ln(s[t].r);
                   15: s[t].r := sqrt(s[t].r);
                   16: s[t].r := ArcTan(s[t].r);
                   17: begin
			 if ir.x = 1 then
			   s[t].b := eof
			 else
			   s[t].b := eof(prr[ir.x]);
		       end;
                   18: begin
			 if ir.x = 1 then
			   s[t].b := eoln
			 else
			   s[t].b := eoln(prr[ir.x]);
		       end;

		   { Modified }
 19 : begin
	s[t-1].r := pwr(s[t-1].r, s[t].r);
	t := t - 1;
      end;
 20 : s[t].r := tan_(s[t].r);
 21 : s[t].r := cosh_(s[t].r);
 22 : s[t].r := sinh_(s[t].r);
 23 : s[t].r := tanh_(s[t].r);
 24 : begin
	s[t+1].r := ranf;
	t := t + 1;
      end;
 25 : s[t].r := dtor(s[t].r);
 26 : ;
 27 : s[t].i := sign(s[t].r);
 28 : begin
        s[t-1].r := getangle(s[t-1].r, s[t].r);
        t := t - 1;
      end;
 29 : begin
	n := s[t-1].i;
	popmat(t, n, mat1);
        s[t-1].r := trmat(n, mat1);
        t := t - 1;
      end;
 30 : begin
	n := s[t-1].i;
	popmat(t, n, mat1);
        s[t-1].r := detmat(n, mat1);
        t := t - 1;
      end;
 31 : begin
	n := s[t-1].i;
	popmat(t, n, mat1);
        s[t-1].b := invmat(n, mat1);
	pushmat(t, n, mat1);
        t := t - 1;
      end;
 32 : begin
	popstr(t, str1);
        s[t].i := strlen_(str1);
      end;
{33 : strind has been put elsewhere }
 34 : begin
        s[t-2].r := Elem_getkval(s[t-2].i, s[t-1].i, s[t].i);
        t := t - 2;
      end;
 35 : begin
	s[t-1].r := min_(s[t-1].r, s[t].r);
	t := t - 1;
      end;
 36 : begin
        s[t-1].r := max_(s[t-1].r, s[t].r);
        t := t - 1;
      end;
 37 : begin
        t := t + 1;
        if t > stacksize then
          ps := stkchk
        else
          s[t].r := normranf;
      end;
 38 : begin
        s[t].r := dble(s[t].i);
      end;
 39 : begin { getnKid }
	s[t].i := ElemFam[s[t].i].nKid;
      end;
 40 : begin { Elem_getpos }
        s[t-1].i := ElemFam[s[t-1].i].KidList[s[t].i];
        t := t - 1;
      end;
 41 : begin { }
      end;
 42 : begin { Elem_getord }
        s[t-1].i := Elem_getord(s[t-1].i, s[t].i);
        t := t - 1;
      end;
 43 : begin { clock }
        s[t+1].i := clock;
	t := t + 1;
      end;
 44 : begin { }
      end;
 45 : begin { }
      end;
 46 : begin { }
      end;

                end;

              9: s[t].i := s[t].i + ir.y;  { offset }

             10: pc := ir.y;     {jump}
             11:
                begin {conditional jump}
                   if not s[t].b  then
                      pc := ir.y;
                   t := t - 1
                end;
             12:
                begin { switch }
                  h1 := s[t].i;   t := t-1;
                  h2 := ir.y;    h3 := 0;
                  repeat
                    if code[h2].f <> 13
                    then begin
                    h3 := 1;
                    ps := caschk;
                  end else if code[h2].y = h1
                           then begin
                             h3 := 1;
                             pc := code[h2+1].y
                           end else h2 := h2 + 2
                  until h3 <> 0
                end;
             14:
                begin {for1up}
                   h1 := s[t -  1].i;
                   if h1 <= s[t].i
                   then
                      s[s[t -   2].i].i := h1
                   else
                      begin
                         t := t - 3;
                         pc := ir.y
                      end
                end;
             15:
                begin {for2up}
                  h2 := s[t-2].i;
                  h1 := s[h2].i + 1;
                  if h1 <= s[t].i
                  then begin
                     s[h2].i := h1;  pc := ir.y
                  end else t := t-3;
                end;
             16:
                begin  {for1down }
                  h1 := s[t-1].i;
                  if h1 >= s[t].i
                  then s[s[t-2].i].i := h1
                  else begin
                    pc := ir.y;  t := t-3
                  end
                end;
             17:
                begin  {for2down}
                  h2 := s[t-2].i;
                  h1 := s[h2].i - 1;
                  if h1 >= s[t].i
                  then begin
                    s[h2].i := h1;  pc := ir.y
                  end else t := t-3;
                end;
             18:
                begin
                   h1 := btab[tab[ir.y].ref1].vsize;
                   if t + h1 > stacksize
                   then
                      ps := stkchk
                   else
                      begin
                         t := t + 5;
                         s[t -  1].i := h1 - 1;
                         s[t].i      := ir.y
                      end;
                end;
             19:
                begin
                   h1 := t - ir.y;
                   h2 := s[h1 + 4].i;     {h2 points to   tab}
                   h3 := tab[h2].lev;
                   display[h3 + 1] :=   h1;
                   h4 := s[h1 + 3].i + h1;
                   s[h1 + 1].i := pc;
                   s[h1 + 2].i := display[h3];
                   s[h1 + 3].i := b;
                   for k:=t+1 to h4 do
		     s[k].i := 0;
                   b := h1;
                   t := h4;
                   pc := tab[h2].adr;
                   if stackdump then dump
                end;
              20:
                begin  { index1 }
                   h1 := ir.y;
                   h2 := atab[h1].low;
                   h3 := s[t].i;
                   if h3 < h2
                   then ps := inxchk
                   else if h3 > atab[h1].high
                        then ps := inxchk
                        else begin
                          t := t - 1;
                          s[t].i := s[t].i + (h3-h2)
                        end
                end;
              21:
                begin {index}
                   h1 := ir.y;  {h1 points to   atab}
                   h2 := atab[h1].low;
                   h3 := s[t].i;
                   if h3 < h2
                   then
                      ps := inxchk
                   else
                      if h3 > atab[h1].high
                      then
                         ps := inxchk
                      else
                         begin
                            t := t - 1;
                            s[t].i := s[t].i+(h3-h2)*atab[h1].elsize
                         end
                end;
            22:             
                begin {load block}
                   h1 := s[t].i;
                   t := t - 1;
                   h2 := ir.y + t;
                   if h2 > stacksize
                   then
                      ps := stkchk
                   else
                      while t < h2 do
                         begin
                            t := t + 1;
                            s[t] := s[h1];
                            h1 := h1 + 1
                         end
                end;
             23:
                begin {copy block}
                   h1 := s[t -  1].i;
                   h2 := s[t].i;
                   h3 := h1 + ir.y;
                   while h1 < h3 do
                      begin
                         s[h1] := s[h2];
                         h1 := h1 + 1;
                         h2 := h2 + 1
                      end;
                   t := t - 2
                end;
             24:
                begin  { literal }
                  t := t+1;
                  if t > stacksize
                  then ps := stkchk 
                  else 
                  begin
                    s[t].i := ir.y;
{                    s[t].b := s[t].i = trueconst; for IBM PC }
                  end;
                end;
             25:
                begin  {load real}
                  t := t+1;
                  if t > stacksize
                  then ps := stkchk else s[t].r := rconst[ir.y];
                end;
             26:
                begin  { float }
                  h1 := t - ir.y;
                  s[h1].r := s[h1].i
                end;
             27:
                begin  { read  }
		  if ir.x = 1 then
		  begin 
                    if eof then
		      ps := redchk
                    else
		      case ir.y of
                        1: read(s[s[t].i].i);
                        2: read(s[s[t].i].r);
                        4: begin
			     read(ch);
			     s[s[t].i].i := ord(ch);
			   end;
                      end;
 		  end
		  else
		  begin
                    if eof(prr[ir.x]) then
		      ps := redchk
                    else
		      case ir.y of
                        1: read(prr[ir.x], s[s[t].i].i);
                        2: read(prr[ir.x], s[s[t].i].r);
                        4: begin
			     read(prr[ir.x], ch);
			     s[s[t].i].i := ord(ch);
			   end;
                      end;
		  end;
                  t := t-1
                end;
             28:
                begin  { write string }
                  h1 := s[t].i;
                  h2 := ir.y;  t := t-1;
                  chrcnt := chrcnt + h1;
                  if chrcnt > lineleng then ps := lngchk;
                  repeat
		    if ir.x = 2 then
		    begin
                      write(stab[h2]); write(psout, stab[h2]);
		    end
		    else
                      write(prr[ir.x], stab[h2]);
                    h1 := h1-1;
                    h2 := h2+1
                  until h1 = 0
                end;
             29:
                begin  { write1 }
                  chrcnt := chrcnt + fld[ir.y];
                  if chrcnt > lineleng
                  then ps := lngchk
                  else case ir.y of
                    1: begin
			 if ir.x = 2 then
			 begin 
                           write(s[t].i:fld[1]); write(psout, s[t].i:fld[1]);
			 end
			 else
                           write(prr[ir.x], s[t].i:fld[1]);
                       end;
                    2: begin
			 if ir.x = 2 then 
			 begin
                           write(s[t].r:fld[2]); write(psout, s[t].r:fld[2]);
			 end
			 else
                           write(prr[ir.x], s[t].r:fld[2]);
                       end;
                    3: if s[t].b then 
                       begin
			 if ir.x = 2 then 
			 begin
                           write('  true'); write(psout, '  true');
			 end
			 else
                           write(prr[ir.x], '  true');
                       end
		       else 
                       begin
			 if ir.x = 2 then 
			 begin
                           write(' false'); write(psout, ' false');
			 end
			 else
                           write(prr[ir.x], ' false');
                       end;
                    4: begin
			 if ir.x = 2 then 
			 begin
                           write(chr(s[t].i)); write(psout, chr(s[t].i));
			 end
			 else
                           write(prr[ir.x], chr(s[t].i));
                       end;
                    end;
                  t := t-1
                 end;
             30:
                begin  { write2  }
                  chrcnt := chrcnt + s[t].i;
                  if chrcnt > lineleng
                  then ps := lngchk
                  else case ir.y of
                     1: begin
			  if ir.x = 2 then 
			  begin
                            write(s[t-1].i:s[t].i); write(psout, s[t-1].i:s[t].i);
			  end
			  else
                            write(prr[ir.x], s[t-1].i:s[t].i);
                        end;
                     2: begin
			  if ir.x = 2 then 
			  begin
                            write(s[t-1].r:s[t].i); write(psout, s[t-1].r:s[t].i);
			  end
			  else
                            write(prr[ir.x], s[t-1].r:s[t].i);
                        end;
                    3: if s[t].b then 
                       begin
			 if ir.x = 2 then 
			 begin
                           write('  true'); write(psout, '  true');
			 end
			 else
                           write(prr[ir.x], '  true');
                       end else 
                       begin
			 if ir.x = 2 then 
			 begin
                           write(' false'); write(psout, ' false');
			 end
			 else
                           write(prr[ir.x], ' false');
                       end;
                    4: begin
			 if ir.x = 2 then 
			 begin
                           write(chr(s[t-1].i):s[t].i);
			   write(psout, chr(s[t-1].i):s[t].i);
			 end
			 else
                           write(prr[ir.x], chr(s[t-1].i):s[t].i);
                       end;
                    end;
                  t := t-2
                end;
             31: ps := fin;
             32:
                begin
                  t := b-1;
                  pc := s[b+1].i;   b := s[b+3].i
                end;
             33:
                begin
                  t := b;
                  pc := s[b+1].i;   b := s[b+3].i
                end;
             34:begin
                   s[t] := s[s[t].i];
                end;
             35: s[t].b := not s[t].b;
             36: case ir.y of
		   0: s[t].i := -s[t].i;
		   1: s[t].r := -s[t].r;
		 end;
             37:begin
                  chrcnt := chrcnt + s[t-1].i;
                  if chrcnt > lineleng
                  then ps := lngchk
                  else 
                  begin
		    if ir.y = 2 then 
		    begin
                      write(s[t-2].r:s[t-1].i:s[t].i);
                      write(psout, s[t-2].r:s[t-1].i:s[t].i);
		    end
		    else
                      write(prr[ir.y], s[t-2].r:s[t-1].i:s[t].i);
                  end;
                  t := t-3
                end;
             38:begin {store}
                   s[s[t - 1].i] := s[t];
                   t := t - 2
                end;
             39:
                begin
                   t := t-1;
                   s[t].b := s[t].r = s[t+1].r
                end;
             40:
                begin
                  t := t-1;
                  s[t].b := s[t].r <> s[t+1].r
                end;
             41:
                begin
                  t := t-1;
                  s[t].b := s[t].r < s[t+1].r
                end;
             42:
                begin
                  t := t-1;
                  s[t].b := s[t].r <= s[t+1].r
                end;
             43:
                begin
                  t := t-1;
                  s[t].b := s[t].r > s[t+1].r
                end;
             44:
                begin
                   t := t-1;
                   s[t].b := s[t].r >= s[t+1].r
                end;
             45:
                begin
                   t := t - 1;
                   s[t].b := s[t].i = s[t+ 1].i
                end;
             46:
                begin
                   t := t - 1;
                   s[t].b := s[t].i <> s[t + 1].i
                end;
             47:
                begin
                   t := t - 1;
                   s[t].b := s[t].i < s[t + 1].i
                end;
             48:
                begin
                   t := t - 1;
                   s[t].b := s[t].i <= s[t + 1].i
                end;
             49:
                begin
                   t := t - 1;
                   s[t].b := s[t].i > s[t + 1].i
                end;
             50:
                begin
                   t := t - 1;
                   s[t].b := s[t].i >= s[t + 1].i
                end;
             51:
                begin
                   t := t - 1;
                   s[t].b := s[t].b or s[t + 1].b
                end;
             52:
                begin
                   t := t - 1;
                   s[t].i := s[t].i + s[t + 1].i
                end;
             53:
                begin
                   t := t - 1;
                   s[t].i := s[t].i - s[t + 1].i
                end;
             54:
                begin
                   t := t - 1;
                   s[t].r := s[t].r + s[t+1].r
                end;
             55:
                begin
                   t := t - 1;
                   s[t].r := s[t].r - s[t+1].r
                end;
             56:
                begin
                   t := t - 1;
                   s[t].b := s[t].b and s[t + 1].b
                end;
             57:
                begin
                   t := t - 1;
                   s[t].i := s[t].i * s[t + 1].i
                end;
             58:
                begin
                   t := t - 1;
                   if s[t + 1].i = 0
                   then
                      ps := divchk
                   else
                      s[t].i :=   s[t].i div s[t + 1].i
                end;
             59:
                begin
                   t := t - 1;
                   if s[t + 1].i = 0
                   then
                      ps := divchk
                   else
                      s[t].i :=   s[t].i mod s[t + 1].i
                end;
             60:
                begin
                  t := t-1;
                  s[t].r := s[t].r * s[t+1].r;
                end;
             61:
                begin
                  t := t-1;
                  s[t].r := s[t].r / s[t+1].r;
                end;
             62: begin
		    if ir.y = 1 then 
		    begin
                      if eof
                      then
                        ps := redchk
                      else
                        readln;
		    end
		    else
		    begin
                      if eof(prr[ir.y])
                      then
                        ps := redchk
                      else
                        readln(prr[ir.y]);
		    end;
		 end;
             63:
                begin
		  if ir.y = 2 then 
		  begin
                    writeln;
                    writeln(psout);
		  end
		  else
                    writeln(prr[ir.y]);
                  lncnt := lncnt + 1;
                  chrcnt := 0;
                {  if lncnt > linelimit then ps := linchk;}
                end;

	    { Modified }
 64 : begin
        iniranf(s[t].i);
        t := t - 1;
      end;
 65 : newseed;
 66 : begin
	n := s[t-1].i;
	popmat(t, n, mat1);
        unitmat(n, mat1);
	pushmat(t, n, mat1);
        t := t - 2;
      end;
 67 : begin
	n := s[t-2].i;
	popmat(t-1, n, mat1);
        copymat(n, mat1, mat2);
	pushmat(t, n, mat2);
        t := t - 3;
      end;
 68 : begin
	n := s[t-2].i;
	popvec(t-1, n, vec1);
        copyvec(n, vec1, vec2);
	pushvec(t, n, vec2);
        t := t - 3;
      end;
 69 : begin
	n := s[t-2].i;
	popvec(t-1, n, vec1);
	popvec(t, n, vec2);
        addvec(n, vec1, vec2);
	pushvec(t, n, vec2);
        t := t - 3;
      end;
 70 : begin
	n := s[t-2].i;
	popvec(t-1, n, vec1);
	popvec(t, n, vec2);
        subvec(n, vec1, vec2);
	pushvec(t, n, vec2);
        t := t - 3;
      end;
 71 : begin
	n := s[t-2].i;
	popmat(t-1, n, mat1);
	popmat(t, n, mat2);
        addmat(n, mat1, mat2);
	pushmat(t, n, mat2);
        t := t - 3;
      end;
 72 : begin
	n := s[t-2].i;
	popmat(t-1, n, mat1);
	popmat(t, n, mat2);
        submat(n, mat1, mat2);
	pushmat(t, n, mat2);
        t := t - 3;
      end;
 73 : begin
	n := s[t-2].i;
	popmat(t-1, n, mat1);
	popmat(t, n, mat2);
        mullmat(n, mat1, mat2);
	pushmat(t, n, mat2);
        t := t - 3;
      end;
 74 : begin
	n := s[t-2].i;
	popmat(t-1, n, mat1);
	popvec(t, n, vec1);
        lintrans(n, mat1, vec1);
	pushvec(t, n, vec1);
        t := t - 3;
      end;
 75 : begin
	n := s[t-1].i;
	popmat(t, n, mat1);
        tpmat(n, mat1);
	pushmat(t, n, mat1);
        t := t - 2;
      end;
 76 : begin
        graope;
      end;
 77 : begin
	popgvec(t-5, s[t-3].i, px);
	popgvec(t-4, s[t-3].i, py);
        graph1(px, py, s[t-3].i, s[t-2].i, s[t-1].i, sngl(s[t].r));
        t := t - 6;
      end;
 78 : gracls;
 79 : begin
	graclx(s[t-3].i, s[t-2].i, s[t-1].i, s[t].i);
        t := t - 4;
      end;
 80 : begin
	popstrbuf(t-2, s[t-1].i, strb1);
	graheds(strb1, s[t-1].i, sngl(s[t].r));
	t := t - 3;
      end;
 81 : begin
	gravwp(sngl(s[t-3].r), sngl(s[t-2].r), sngl(s[t-1].r), sngl(s[t].r));
        t := t - 4;
      end;
 82 : begin
	grawnd(sngl(s[t-3].r), sngl(s[t-2].r), sngl(s[t-1].r), sngl(s[t].r));
        t := t - 4;
      end;
 83 : begin
	popstrbuf(t-5, s[t-4].i, strb1);
	popstrbuf(t-2, s[t-1].i, strb2);
	gracoms(sngl(s[t-6].r), strb1, s[t-4].i, s[t-3].i, strb2, s[t-1].i,
	        s[t].i);
        t := t - 7;
      end;
 84 : begin
	graxi1(s[t-3].i, s[t-2].i, s[t-1].i, sngl(s[t].r));
	t := t - 4;
      end;
 85 : begin
	grayi1(s[t-3].i, s[t-2].i, s[t-1].i, sngl(s[t].r));
	t := t - 4;
      end;
 86 : begin
	graxi2(sngl(s[t-2].r), s[t-1].i, s[t].i);
        t := t - 3;
      end;
 87 : begin
	grayi2(sngl(s[t-2].r), s[t-1].i, s[t].i);
        t := t - 3;
      end;
 88 : begin
	popstrbuf(t-5, s[t-4].i, strb1);
	gratexs(sngl(s[t-8].r), sngl(s[t-7].r), sngl(s[t-6].r), strb1,
	        s[t-4].i, sngl(s[t-3].r), snglx, sngly, s[t].i);
	s[s[t-2].i].r := snglx; s[s[t-1].i].r := sngly;
	t := t - 9;
      end;
 89 : begin
	gravec(sngl(s[t-3].r), sngl(s[t-2].r), sngl(s[t-1].r), sngl(s[t].r));
	t := t - 4;
      end;
 90 : begin
	gratyp(s[t].i);
	t := t - 1;
      end;
 91 : begin
	popDAmap(t-3, map);
        Cell_DATwiss(s[t-5].i, s[t-4].i, map, s[t-2].b, s[t-1].b, s[t].r);
        t := t - 6;
      end;
 92 : begin
	gramar(sngl(s[t-3].r), sngl(s[t-2].r), s[t-1].i, sngl(s[t].r));
	t := t - 4;
      end;
 93 : begin
	popmat(t-2, 6, mat4);
	gdiag(s[t-7].i, s[t-6].r, mat1, mat2, mat3, mat4, s[s[t-1].i].r,
	      s[s[t].i].r);
	pushmat(t-5, 6, mat1);
	pushmat(t-4, 6, mat2);
	pushmat(t-3, 6, mat3);
	t := t - 8;
      end;
 94 : begin
	m := s[t-3].i; n := s[t-2].i;
	popsvdmat(t-6, false, m, n, u1);
	popvec(t-5, n, w1);
	popsvdmat(t-4, false, n, n, v1);
	popvec(t-1, m, b1);
	svbksb(u1, w1, v1, m, n, b1, x1);
	pushvec(t, n, x1);
	t := t - 7;
      end;
 95 : begin
	m := s[t-3].i; n := s[t-2].i;
	popsvdmat(t-4, true, m, n, a1);
	mn := mnp; job := 11;
	{ dsvcmp from numerical recipies has been replaced by
	  dsvdc from linpack, due to numerical problems }
	dsvdc(a1, mn, m, n, w1, e1, u1, mn, v1, mn, work, job, info);
	pushsvdmat(t-4, true, m, n, u1);
	pushvec(t-1, n, w1);
	pushsvdmat(t, true, n, n, v1);
	t := t - 5;
      end;
 96 : begin{ break }
	write('break at ', linerror[pc-1]:1);
	pmdump(s[t].b);
      end;
 97 : begin
 	popstrcon(t, str1);
        getstr(str2, str1);
	pushstr(t-1, str2);
        t := t - 2;
      end;
 98 : begin
	popstr(t, str1);
        copystr(str2, str1);
	pushstr(t-1, str2);
        t := t - 2;
      end;
 99 : begin
	popstrcon(t, str1);
	popstr(t-1, str2);
        concat(str2, str1);
	pushstr(t-1, str2);
        t := t - 2;
      end;
100 : begin
	popstr(t-2, str1);
        getint(str1, s[t-1].i, s[t].i);
	pushstr(t-2, str1);
        t := t - 3;
      end;
101 : begin
	popstr(t-3, str1);
        getreal(str1, s[t-2].r, s[t-1].i, s[t].i);
	pushstr(t-3, str1);
        t := t - 4;
      end;
102 : begin
	popstr(t-3, str1);
        getreale(str1, s[t-2].r, s[t-1].i, s[t].i);
	pushstr(t-3, str1);
        t := t - 4;
      end;
103 : begin
	popstr(t-1, str1);
        writestr(str1, s[t].i, output);
        writestr(str1, s[t].i, psout);
        t := t - 2;
      end;
104 : begin
	popgvec(t-2, s[t-3].i, px);
	popgvec(t-1, s[t-3].i, py);
        graph0(px, py, s[t].i);
        t := t - 3;
      end;
105 : begin
	n := s[t-2].i;
	popmat(t-1, n, mat1);
	popmat(t, n, mat2);
        mulrmat(n, mat1, mat2);
	pushmat(t-1, n, mat1);
        t := t - 3;
      end;
106 : begin
	popvec(t-1, 6, vec1);
        Cell_Pass(s[t-3].i, s[t-2].i, vec1, s[s[t].i].i);
	pushvec(t-1, 6, vec1);
        t := t - 4;
      end;
107 : begin
	popvec(t-2, 6, vec1);
	popmat(t-1, 6, mat1);
        Cell_Pass_M(s[t-4].i, s[t-3].i, vec1, mat1, s[s[t].i].i);
	pushvec(t-2, 6, vec1);
	pushmat(t-1, 6, mat1);
        t := t - 5;
      end;
108 : begin
	setrancut(s[t].r);
	t := t - 1;
      end;
109 : begin { Elem_getbend }
	Elem_GetBend(s[t-4].i, s[t-3].i, s[s[t-2].i].r, s[s[t-1].i].r,
		        s[s[t].i].r);
	t := t - 5;
      end;
110 : begin
        Cell_Getcod(s[t-3].i, s[t-2].r, s[t-1].r, s[s[t].i].i);
        t := t - 4;
      end;
111 : begin
	popmat(t-4, 6, mat1);
        Cell_GetABGN(mat1, vec21, vec22, vec23, vec24);
	pushvec(t-3, 2, vec21);
	pushvec(t-2, 2, vec22);
	pushvec(t-1, 2, vec23);
	pushvec(t, 2, vec24);
        t := t - 5;
      end;
112 : begin
	trace := s[t].b;
	t := t - 1;
      end;
113 : begin
	popmat(t-3, 6, mat1);
        Cell_MatTwiss(s[t-5].i, s[t-4].i, mat1, s[t-2].b, s[t-1].b, s[t].r);
        t := t - 6;
      end;
114 : begin
	n := s[t-2].i;
	popbvec(t-1, n, bx);
	popbvec(t, n, by);
	fft(n, bx, by);
	pushbvec(t-1, n, bx);
	pushbvec(t, n, by);
        t := t - 3;
      end;
115 : begin
        Ring_GetTwiss(s[t-1].b, s[t].r);
        t := t - 2;
      end;
116 : begin
	popvec(t-6, 2, vec21);
	popivec(t-4, 2, ivec21);
	popivec(t-3, ivec21[1], xf);
	popivec(t-2, ivec21[2], xd);
        Ring_Fittune(vec21, s[t-5].r, ivec21, xf, xd, s[t-1].r, s[t].i);
        t := t - 7;
      end;
117 : begin
	popvec(t-6, 2, vec21);
	popivec(t-4, 2, ivec21);
	popivec(t-3, ivec21[1], xf);
	popivec(t-2, ivec21[2], xd);
        Ring_Fitchrom(vec21, s[t-5].r, ivec21, xf, xd, s[t-1].r, s[t].i);
        t := t - 7;
      end;
118 : begin
	popivec(t-2, s[t-3].i, xf);
        Ring_Fitdisp(s[t-6].i, s[t-5].r, s[t-4].r, s[t-3].i, xf,
		     s[t-1].r, s[t].i);
        t := t - 7;
      end;
119 : begin
      end;
120 : begin
        Cell_Concat(s[t].r);
        t := t - 1;
      end;
121 : begin
	popvec(t-1, 5, vec1);
	Cell_fPass(vec1, s[s[t].i].i);
	pushvec(t-1, 5, vec1);
        t := t - 2;
      end;
122 : begin { getelem }
	pushelem(t, s[t-1].i);
        t := t - 2;
      end;
123 : begin { Cell_DApass }
	popDAmap(t-1, map);
        Cell_DApass(s[t-3].i, s[t-2].i, map, s[s[t].i].i);
	pushDAmap(t-1, map);
        t := t - 4;
      end;
124 : begin
	popivec(t-4, 2, ivec21);
	popivec(t-3, ivec21[1], xf);
	popivec(t-2, ivec21[2], xd);
	initbump(ivec21, xf, xd, s[t-1].r, s[t].r);
        t := t - 5;
      end;
125 : begin
	execbump(s[t-1].r, s[t].i);
        t := t - 2;
      end;
126 : begin { getglobv_ }
	pushglobval(t);
        t := t - 1;
      end;
127 : begin { putglobv_ }
	popglobval(t);
        t := t - 1;
      end;
128 : begin
        grabeg;
      end;
129 : begin
        graend;
      end;
130 : begin
	popgvec(t-2, s[t].i, px);
	popgvec(t-1, s[t].i, py);
        grafpl(px, py, s[t].i);
        t := t - 3;
      end;
131 : begin
	gracoi(s[t].i);
	t := t - 1;
      end;
132 : begin
	grapat(s[t].i);
	t := t - 1;
      end;
133 : begin
	gravcs(sngl(s[t-3].r), sngl(s[t-2].r), sngl(s[t-1].r), sngl(s[t].r));
	t := t - 4;
      end;
134 : begin
	gravcf(s[t].i);
	t := t - 1;
      end;
135 : begin
      end;
136 : begin
      end;
137 : begin
      end;
138 : begin{ reset_ }
	popstrbuf(t, strlmax, strb1);
	k := ir.y; { for SUN compability }
        reset_(prr[k], strb1);
        t := t - 2;
      end;
139 : begin{ rewrite_ }
	popstrbuf(t, strlmax, strb1);
	k := ir.y; { for SUN compability }
        rewrite_(prr[k], strb1);
        t := t - 2;
      end;
140 : begin
	close(prr[ir.y]);
        t := t - 1;
      end;
141 : begin
      end;
142 : begin
      end;
143 : begin
      end;
144 : begin
      end;
145 : begin
        Mpole_setds(s[t-1].i, s[t].i);
        t := t - 2;
      end;
146 : begin
      end;
147 : begin
      end;
148 : begin
      end;
149 : begin
      end;
150 : begin
        Mpole_setdt(s[t-1].i, s[t].i);
        t := t - 2;
      end;
151 : begin
      end;
152 : begin
      end;
153 : begin
      end;
154 : begin
        Mpole_setpb(s[t-2].i, s[t-1].i, s[t].i);
        t := t - 3;
      end;
155 : begin
        Mpole_GetPmeth(s[t-7].i, s[t-6].i, s[s[t-5].i].i, s[s[t-4].i].i,
		       s[s[t-3].i].r, s[s[t-2].i].r, s[s[t-1].i].r,
		       s[s[t].i].r);
        t := t - 8;
      end;
156 : begin
        Cav_Set(s[t-4].i, s[t-3].i, s[t-2].r, s[t-1].r, s[t].i);
        t := t - 5;
      end;
157 : begin
        Cav_Get(s[t-4].i, s[t-3].i, s[s[t-2].i].r, s[s[t-1].i].r, s[s[t].i].i);
        t := t - 5;
      end;
158 : begin
        GetDBN(s[t-2].i, s[t-1].i, DBname);
        pushDBname(t, DBname);
        t := t - 3;
      end;
159 : begin
      end;
160 : begin
      end;
161 : begin
      end;
162 : begin
      end;
163 : begin
      end;
164 : begin
      end;
165 : begin
      end;
166 : begin
      end;
167 : begin
      end;
168 : begin
      end;
169 : begin
      end;
170 : begin { amoeba }
	popvec(t-6, s[t-4].i, simpvec1);
	popvec(t-5, s[t-4].i, simpvec2);
	amoeba(simpvec1, simpvec2, s[t-4].i, s[s[t-3].i].r, s[t-2].r, s[t-1].i,
	       s[s[t].i].i);
	pushvec(t-6, s[t-4].i, simpvec1);
        t := t - 7;
      end;
171 : begin { Cell_SetdP }
	Cell_SetdP(s[t].r);
        t := t - 1;
      end;
172 : begin { putelem }
	popelem(t, s[t-1].i);
        t := t - 2;
      end;

          end;{ of Case }
   end{ ExecPcode};

   begin {interpret}
          s[1].i := 0;      s[2].i := 0;
          s[3].i := -1;     s[4].i := btab[1].last;
      display[1] := 0;           t := btab[2].vsize - 1;
               b := 0;          pc := tab[s[4].i].adr;
           lncnt := 0;        ocnt := 0;
          chrcnt := 0;          ps := run;

          fld[1] := 10;     fld[2] := 23;
          fld[3] := 10;     fld[4] := 1;
      
      writeln(psout);
      writeln(psout,'  Total of ', lc:1, ' steps (max ', cmax:1, ' steps)');

      repeat
          ir := code[pc];
          pc := pc + 1;    ocnt := ocnt + 1;
          ExecPcode;
      until ps <> run;

      if ps <>  fin
      then
          begin
            writeln;
            write(' Halt at line ', linerror[pc-1]:1,
	      ' of source because of ');
             case ps of
                caschk:
                   writeln('Undefined case');
                divchk:
                   writeln('Division by 0');
                inxchk:
                   writeln('Invalid index');
                stkchk:
                   writeln('Storage overflow');
                linchk:
                   writeln('Too much output');
                lngchk:
                   writeln('Line too long');
                redchk:
                   writeln('Reading past end of file');
             end;

	    pmdump(false);
          end;
end; { interpret }


PROCEDURE init_reserved_words;

  PROCEDURE Reg(nam:alfa; sy:symbol);

  begin
    nkw:=nkw+1;
    Key[nkw]:=nam; ksy[nkw]:=sy
  end;

begin
   nkw:=0;
   Reg('and            ', andsy   );   Reg('array          ', arraysy );
   Reg('begin          ', beginsy );   Reg('case           ', casesy  );
   Reg('const          ', constsy );   Reg('div            ', idiv    );
   Reg('do             ', dosy    );   Reg('downto         ', downtosy);
   Reg('else           ', elsesy  );   Reg('end            ', endsy   );
   Reg('for            ', forsy   );   Reg('function       ', funcsy  );
   Reg('if             ', ifsy    );   Reg('include        ', incsy   );
   Reg('mod            ', imod    );
   Reg('not            ', notsy   );   Reg('of             ', ofsy    );
   Reg('or             ', orsy    );   Reg('procedure      ', procsy  );
   Reg('program        ', programsy);  Reg('record         ', recordsy);
   Reg('repeat         ', repeatsy);   Reg('then           ', thensy  );
   Reg('to             ', tosy    );   Reg('type           ', typesy  );
   Reg('until          ', untilsy );   Reg('var            ', varsy   );
   Reg('while          ', whilesy );

   sps['+'] := plus;      sps['-'] := minus;
   sps['('] := lparent;   sps[')'] := rparent;
   sps['='] := eql;       sps[','] := comma;
   sps['[']:=lbrack;      sps[']'] := rbrack;
   sps[''''] := neq;      sps['&'] := andsy;
   sps[';'] := semicolon; sps['*'] := times;
   sps['/'] := rdiv;
end;


PROCEDURE enterdis;

begin
   enterstandardids('               ',vvariable, notyp, 0);   {sentinel}
   enterstandardids('false          ',konstant, bools, falseconst);
   enterstandardids('true           ',konstant, bools, trueconst);
   enterstandardids('text           ',type1, texts, 1);
   enterstandardids('double         ',type1, reals, 1);
   enterstandardids('char           ',type1, chars, 1);
   enterstandardids('boolean        ',type1, bools, 1);
   enterstandardids('integer        ',type1, ints, 1);
   enterstandardids('abs            ',funktion, reals, 0);
   enterstandardids('sqr            ',funktion, reals, 2);
   enterstandardids('odd            ',funktion, bools, 4);
   enterstandardids('chr            ',funktion, chars, 5);
   enterstandardids('ord            ',funktion, ints,  6);
   enterstandardids('succ           ',funktion, chars, 7);
   enterstandardids('pred           ',funktion, chars, 8);
   enterstandardids('round          ',funktion, ints , 9);
   enterstandardids('trunc          ',funktion, ints ,10);
   enterstandardids('sin            ',funktion, reals,11);
   enterstandardids('cos            ',funktion, reals,12);
   enterstandardids('exp            ',funktion, reals,13);
   enterstandardids('ln             ',funktion, reals,14);
   enterstandardids('sqrt           ',funktion, reals,15);
   enterstandardids('arctan         ',funktion, reals,16);
   enterstandardids('eof            ',funktion, bools,17);
   enterstandardids('eoln           ',funktion, bools,18);

   { Modified }
   enterstandardids('pwr            ', funktion, reals, 19);
   enterstandardids('tan_           ', funktion, reals, 20);
   enterstandardids('cosh_          ', funktion, reals, 21);
   enterstandardids('sinh_          ', funktion, reals, 22);
   enterstandardids('tanh_          ', funktion, reals, 23);
   enterstandardids('ranf           ', funktion, reals, 24);
   enterstandardids('dtor           ', funktion, reals, 25);
   enterstandardids('               ', funktion, reals, 26);
   enterstandardids('sign           ', funktion, ints,  27);
   enterstandardids('getangle       ', funktion, reals, 28);
   enterstandardids('trmat          ', funktion, reals, 29);
   enterstandardids('detmat         ', funktion, reals, 30);
   enterstandardids('invmat         ', funktion, bools, 31);
   enterstandardids('strlen_        ', funktion, ints,  32);
   enterstandardids('strind         ', funktion, ints,  33);
   enterstandardids('elem_getkval   ', funktion, reals, 34);
   enterstandardids('min_           ', funktion, reals, 35);
   enterstandardids('max_           ', funktion, reals, 36);
   enterstandardids('normranf       ', funktion, reals, 37);
   enterstandardids('dble           ', funktion, reals, 38);
   enterstandardids('getnkid        ', funktion, ints,  39);
   enterstandardids('elem_getpos    ', funktion, ints,  40);
   enterstandardids('               ', funktion, reals, 41);
   enterstandardids('elem_getord    ', funktion, ints,  42);
   enterstandardids('clock          ', funktion, ints,  43);
   enterstandardids('               ', funktion, reals,  44);
   enterstandardids('               ', funktion, reals,  45);
   enterstandardids('               ', funktion, reals,  46);

   enterstandardids('read           ',prozedure, notyp, 1);
   enterstandardids('readln         ',prozedure, notyp, 2);
   enterstandardids('write          ',prozedure, notyp, 3);
   enterstandardids('writeln        ',prozedure, notyp, 4);

   { Modified }
   enterstandardids('iniranf        ', prozedure, notyp,   5);
   enterstandardids('newseed        ', prozedure, notyp,   6);
   enterstandardids('unitmat        ', prozedure, notyp,   7);
   enterstandardids('copymat        ', prozedure, notyp,   8);
   enterstandardids('copyvec        ', prozedure, notyp,   9);
   enterstandardids('addvec         ', prozedure, notyp,  10);
   enterstandardids('subvec         ', prozedure, notyp,  11);
   enterstandardids('addmat         ', prozedure, notyp,  12);
   enterstandardids('submat         ', prozedure, notyp,  13);
   enterstandardids('mullmat        ', prozedure, notyp,  14);
   enterstandardids('lintrans       ', prozedure, notyp,  15);
   enterstandardids('tpmat          ', prozedure, notyp,  16);

   enterstandardids('graope         ', prozedure, notyp,  17);
   enterstandardids('graph1         ', prozedure, notyp,  18);
   enterstandardids('gracls         ', prozedure, notyp,  19);
   enterstandardids('graclx         ', prozedure, notyp,  20);
   enterstandardids('graheds        ', prozedure, notyp,  21);
   enterstandardids('gravwp         ', prozedure, notyp,  22);
   enterstandardids('grawnd         ', prozedure, notyp,  23);
   enterstandardids('gracoms        ', prozedure, notyp,  24);
   enterstandardids('graxi1         ', prozedure, notyp,  25);
   enterstandardids('grayi1         ', prozedure, notyp,  26);
   enterstandardids('graxi2         ', prozedure, notyp,  27);
   enterstandardids('grayi2         ', prozedure, notyp,  28);
   enterstandardids('gratexs        ', prozedure, notyp,  29);
   enterstandardids('gravec         ', prozedure, notyp,  30);
   enterstandardids('gratyp         ', prozedure, notyp,  31);
   enterstandardids('cell_datwiss   ', prozedure, notyp,  32);
   enterstandardids('gramar         ', prozedure, notyp,  33);
   enterstandardids('gdiag          ', prozedure, notyp,  34);
   enterstandardids('svbksb         ', prozedure, notyp,  35);
   enterstandardids('svdcmp         ', prozedure, notyp,  36);
   enterstandardids('break          ', prozedure, notyp,  37);
   enterstandardids('getstr         ', prozedure, notyp,  38);
   enterstandardids('copystr        ', prozedure, notyp,  39);
   enterstandardids('concat         ', prozedure, notyp,  40);
   enterstandardids('getint         ', prozedure, notyp,  41);
   enterstandardids('getreal        ', prozedure, notyp,  42);
   enterstandardids('getreale       ', prozedure, notyp,  43);
   enterstandardids('writestr       ', prozedure, notyp,  44);
   enterstandardids('graph0         ', prozedure, notyp,  45);
   enterstandardids('mulrmat        ', prozedure, notyp,  46);
   enterstandardids('cell_pass      ', prozedure, notyp,  47);
   enterstandardids('cell_pass_m    ', prozedure, notyp,  48);
   enterstandardids('setrancut      ', prozedure, notyp,  49);
   enterstandardids('elem_getbend   ', prozedure, notyp,  50);
   enterstandardids('cell_getcod    ', prozedure, notyp,  51);
   enterstandardids('cell_getabgn   ', prozedure, notyp,  52);
   enterstandardids('trace          ', prozedure, notyp,  53);
   enterstandardids('cell_mattwiss  ', prozedure, notyp,  54);
   enterstandardids('fft            ', prozedure, notyp,  55);
   enterstandardids('ring_gettwiss  ', prozedure, notyp,  56);
   enterstandardids('ring_fittune   ', prozedure, notyp,  57);
   enterstandardids('ring_fitchrom  ', prozedure, notyp,  58);
   enterstandardids('ring_fitdisp   ', prozedure, notyp,  59);
   enterstandardids('               ', prozedure, notyp,  60);
   enterstandardids('cell_concat    ', prozedure, notyp,  61);
   enterstandardids('cell_fpass     ', prozedure, notyp,  62);
   enterstandardids('getelem        ', prozedure, notyp,  63);
   enterstandardids('cell_dapass    ', prozedure, notyp,  64);
   enterstandardids('initbump       ', prozedure, notyp,  65);
   enterstandardids('execbump       ', prozedure, notyp,  66);
   enterstandardids('getglobv_      ', prozedure, notyp,  67);
   enterstandardids('putglobv_      ', prozedure, notyp,  68);
   enterstandardids('grabeg         ', prozedure, notyp,  69);
   enterstandardids('graend         ', prozedure, notyp,  70);
   enterstandardids('grafpl         ', prozedure, notyp,  71);
   enterstandardids('gracoi         ', prozedure, notyp,  72);
   enterstandardids('grapat         ', prozedure, notyp,  73);
   enterstandardids('gravcs         ', prozedure, notyp,  74);
   enterstandardids('gravcf         ', prozedure, notyp,  75);
   enterstandardids('               ', prozedure, notyp,  76);
   enterstandardids('               ', prozedure, notyp,  77);
   enterstandardids('               ', prozedure, notyp,  78);
   enterstandardids('reset_         ', prozedure, notyp,  79);
   enterstandardids('rewrite_       ', prozedure, notyp,  80);
   enterstandardids('close          ', prozedure, notyp,  81);
   enterstandardids('               ', prozedure, notyp,  82);
   enterstandardids('               ', prozedure, notyp,  83);
   enterstandardids('               ', prozedure, notyp,  84);
   enterstandardids('               ', prozedure, notyp,  85);
   enterstandardids('mpole_setds    ', prozedure, notyp,  86);
   enterstandardids('               ', prozedure, notyp,  87);
   enterstandardids('               ', prozedure, notyp,  88);
   enterstandardids('               ', prozedure, notyp,  89);
   enterstandardids('               ', prozedure, notyp,  90);
   enterstandardids('mpole_setdt    ', prozedure, notyp,  91);
   enterstandardids('               ', prozedure, notyp,  92);
   enterstandardids('               ', prozedure, notyp,  93);
   enterstandardids('               ', prozedure, notyp,  94);
   enterstandardids('mpole_setpb    ', prozedure, notyp,  95);
   enterstandardids('mpole_getpmeth ', prozedure, notyp,  96);
   enterstandardids('cav_set        ', prozedure, notyp,  97);
   enterstandardids('cav_get        ', prozedure, notyp,  98);
   enterstandardids('getdbn         ', prozedure, notyp,  99);

   enterstandardids('               ', prozedure, notyp, 100);
   enterstandardids('               ', prozedure, notyp, 101);
   enterstandardids('               ', prozedure, notyp, 102);
   enterstandardids('               ', prozedure, notyp, 103);
   enterstandardids('               ', prozedure, notyp, 104);
   enterstandardids('               ', prozedure, notyp, 105);
   enterstandardids('               ', prozedure, notyp, 106);
   enterstandardids('               ', prozedure, notyp, 107);
   enterstandardids('               ', prozedure, notyp, 108);
   enterstandardids('               ', prozedure, notyp, 109);
   enterstandardids('               ', prozedure, notyp, 110);

   enterstandardids('amoeba         ', prozedure, notyp, 111);
   enterstandardids('cell_setdp     ', prozedure, notyp, 112);
   enterstandardids('putelem        ', prozedure, notyp, 113);

   enterstandardids('               ',prozedure, notyp, 0);
end;


  procedure  getlatfilnam;

  var	i, len	: integer;

  begin
    insymbol;
    if sy = ident then
    begin
      len := 0;
      while id[len+1] <> ' ' do
      begin
        len := succ(len); fname[len] := id[len];
      end;
      insymbol;
      if sy = period then
      begin
        len := succ(len);
        fname[len] := '.';
        insymbol;
        if sy = ident then
        begin
          i := 0;
          while id[i+1] <> ' ' do
          begin
            i := succ(i); fname[len+i] := id[i];
          end;
          len := len + i;
          for i:=len+1 to 40 do
            fname[i] := ' ';

          reset_(fi, fname);
          foname := fname; foname[len] := 'x';
          rewrite_(fo, foname);

          insymbol;
          if sy = rparent then
            insymbol
          else
            error(4);
        end
	  else error(4);
      end
      else
        error(4);
    end
    else
      error(2);
    if sy <> semicolon then error(14);
  end;


begin { Pascals }
   init_reserved_words;

   constbegsys  := [plus, minus, intcon, realcon, charcon, ident];
   typebegsys   := [ident, arraysy, recordsy];
   blockbegsys  := [constsy, typesy, varsy, procsy, funcsy, beginsy];
   facbegsys    := [intcon, realcon, charcon, ident, lparent, notsy];
   statbegsys   := [beginsy, ifsy, whilesy, repeatsy, forsy, casesy];
   stantyps     := [notyp, ints, reals, bools, chars, texts];

   t:=-1;
   enterdis;

   linecounter := 0;

   lc := 0;  ll := 0;  cc := 0;   errpos := 0;  errs := [];  chin := ' ';
   t := - 1;  a := 0;   b := 1;       sx := 0;    c2 := 0;
   prtlin := false; incllev := 0; filnb := 2;

   write('_File: ' );
   fnamelen := 0;
   while not eoln do
   begin
     fnamelen := succ(fnamelen);
     read(fname[fnamelen]);
   end;
   writeln;

   for i:=1 to 40 do
   begin
     if i <= fnamelen then
     begin
       inf[i]  := fname[i]; outf[i] := fname[i];
     end
     else
     begin
       inf[i]  := ' '; outf[i] := ' ';
     end;
   end;

   inf[fnamelen+1] := '.'; inf[fnamelen+2] := 'i';
   inf[fnamelen+3] := 'n'; inf[fnamelen+4] := 'p';
   outf[fnamelen+1] := '.'; outf[fnamelen+2] := 'o';
   outf[fnamelen+3] := 'u'; outf[fnamelen+4] := 't';

   reset_(psin[incllev], inf);
   rewrite_(psout, outf);

   writeln;
   writeln('  Source input  file : ', inf);
   writeln('  Source output file : ', outf);

   prtsource := false; prtables := false;  stackdump := false;

   insymbol;

   display[0] := 1;
   skipflag := false;  iflag := false;  oflag := false;

   if sy <> programsy
   then
     error(3)
   else
   begin
     insymbol;{ Get Program Name } 
     if sy <> ident then
       error(2)
     else
     begin
       progname := id;
       iflag := true;
       oflag := true;
       insymbol;{ Get Semicolon? }
     end;

     if sy = lparent then
     begin
       getlatfilnam;
       if Lattice_Read(fi, fo) then
	 Cell_Init
       else
       begin
         writeln('** Error in Lattice_Read');
	 goto 9999;
       end;
     end;
   end;

   enterdis;

   with btab[1] do
   begin
     last := t;
     lastpar := 1;
     psize := 0;
     vsize := 0
   end;

   block(blockbegsys +  statbegsys, false, 1);

   close(psin[incllev]);

   if sy <> period then error(22);
     emit(31);  { halt }

   if prtables then printtables;

   if errs = [] then
   begin
     writeln;
     interpret;
   end
   else
   begin
     writeln(psout);
     writeln(psout,'Compiled with errors');
     errormsg;
   end;
9999:
   close(psout);
end; { Pascals }


begin { Main program }
  t2init;
  Pascals;
end.
