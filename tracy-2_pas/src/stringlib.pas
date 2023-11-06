	module stringlib(input, output);

	%include 'mathlib.def'
	%include 'stringlib.def'

	{ Interface }

	[global] function strlen_(var str : tstring) : integer; forward;

	[global] procedure getstr(var str : tstring;
                         vstr : packed array [low..high : integer] of char);
                         forward;

	{ probably not needed, can be done by concat }
	[global] procedure copystr(var outstr, instr : tstring); forward;

	[global] procedure concat(var str2 : tstring;
                         str1 : packed array [low..high : integer] of char);
                         forward;

	[global] procedure getint(var str : tstring; i, blanks : integer); forward;

	[global] procedure getreal(var str : tstring; x : double;
			  blanks, ndec : integer); forward;

	[global] procedure getreale(var str : tstring; x : double;
			   blanks, ndec : integer); forward;

	[global] function strind(var object: tstring;
                        pattern : packed array [low..high : integer] of char)
                        : integer; forward;

	[global] procedure writestr(var str : tstring; blanks : integer;
			   var outf : text); forward;

	{ Implementation }

	function strlen_{var str : tstring}{ : integer};

	{ Get length of str }

	begin
	  strlen_ := str.len;
	end;


	procedure getstr{var str : tstring;
                         vstr : packed array [low..high : integer] of char};

	{ Put vstr in str }

	var	i	: integer;

	begin
	  for i:=low to high do
	    str.str[i] := vstr[i];
	  str.len := high - low + 1;
	  for i:=high+1 to strlmax do
	    str.str[i] := ' ';
	end;


	procedure copystr{var outstr, instr : tstring};

	{ Copy instr to outstr }

	var	i, len	: integer;

	begin
	  len := strlen_(instr);
	  for i:=1 to len do
	    outstr.str[i] := instr.str[i];
	  outstr.len := len;
	end;


	procedure concat{var str2 : tstring;
                         str1 : packed array [low..high : integer] of char};

	{ Concatenate str1 to str2 }

	var	i, len2	: integer;

	begin
	  len2 := strlen_(str2);
	  for i:=low to high do
	    str2.str[len2+i] := str1[i];
	  str2.len := high - low + 1 + len2;
	end;


	procedure getint{var str : tstring; i, blanks : integer};

	{ Concatenate integer i to str. Will fill with leading zeroes
	  if blanks is negative. }

	var	iabs, j, k, n, len, x	: integer;

	begin
	  iabs := abs(i);
	  x := 1; n := 0;
	  while iabs div x >= 10 do
	  begin
	    x := x*10; n := succ(n);
	  end;

	  len := strlen_(str); k := abs(blanks)-(n+1);
	  if i < 0 then
	  begin
	    k := k - 1;
	    if blanks < 0 then
	    begin
	      len := len + 1;
	      str.str[len] := '-';
	    end;
	  end;

	  for j:=1 to k do
	    if blanks >= 0 then
	      str.str[len+j] := ' '
	    else
	      str.str[len+j] := '0';
	  if k > 0 then len := len + k;

	  if ( blanks >= 0 ) and ( i < 0 ) then
	  begin
	    len := len + 1; str.str[len] := '-';
	  end;

	  for j:=n downto 0 do
	  begin
	    k := iabs div x;
	    str.str[len+n-j+1] := chr(k+ord('0'));
	    iabs := iabs - k*x; x := x div 10;
	  end;

	  str.len := len+n+1;
	end;


	procedure getreal{var str : tstring; x : double; blanks, ndec : integer};

	{ Concatenate real x to str, fix notation }

	var	ix, xdec, i		: integer;
		x1, potens, fracx	: double;
		str1			: packed array [1..2] of char;

	begin
	  potens := pwr(10, ndec);
	  ix := trunc(x);
	  fracx := abs(x-ix)*potens;
	  xdec := trunc(fracx);
	  if fracx-xdec > 0.5 then
	  begin
	    if x >= 0 then
	      x1 := x + 0.1/potens
	    else
	      x1 := x - 0.1/potens;
	    ix := trunc(x1);
	    fracx := abs(x1-ix)*potens;
	    xdec := trunc(fracx);
	  end;
	  i := blanks-1-ndec; if i < 0 then i := 0;
	  { Assure correct sign for -1 < x < 0 }
	  if (ix = 0) and (x < 0) and (xdec <>0) then
	  begin
	    getint(str, ix, i-1);
	    i := strlen_(str);
	    str.str[i] := '-'; str.str[i+1] := '0';
	    str.len := i + 1;
	  end
	  else
	    getint(str, ix, i);
          str1 := '. ';
	  concat(str, str1);
          str.len := str.len - 1;
	  getint(str, xdec, -ndec);
	end;


	procedure getreale{var str : tstring; x : double; blanks, ndec : integer};

	{ Concatenate real x to str, scientific notation }

	var	iexp	: integer;
		absx	: double;
		str1	: packed array [1..2] of char;

	begin
	  absx := abs(x); if (absx < 1) and (absx <> 0) then absx := 1/absx;
	  iexp := 0;
	  while absx/pwr(10, iexp) >= 10 do
	    iexp := succ(iexp);
	  absx := abs(x);
	  if (absx < 1) and (absx <> 0) then iexp := -iexp-1;
	  getreal(str, x/pwr(10, iexp), blanks-4, ndec);
	  if iexp >= 0 then
	  begin
            str1 := 'E+';
	    concat(str, str1);
	    getint(str, iexp, -2);
	  end
	  else
	  begin
            str1 := 'E ';
	    concat(str, str1);
            str.len := str.len - 1;
	    getint(str, iexp, -3);
	  end;
	end;


	function strind{var object: tstring;
                        pattern : packed array [low..high : integer] of char}
                        { : integer};

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
	  plen := high - low + 1;
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
	      if pattern[low] = object.str[i] then
	      begin
	        ind := i;
	        for j:=low+1 to plen do
		  if i+j-1 <= olen then
		  begin
		    if pattern[j] <> object.str[i+j-1] then ind := 0;
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


	procedure writestr{var str : tstring; blanks : integer; var outf : text};

	{ Write tstring. If blanks is positive then right adjusted else
	  left adjusted }

	var	i, len	: integer;

	begin
	  len := strlen_(str);
	  if blanks > 0 then
	    for i:=1 to blanks-len do
	      write(outf, ' ');
	  for i:=1 to len do
	    write(outf, str.str[i]);
	  if blanks < 0 then
	    for i:=1 to -blanks-len do
	      write(outf, ' ');
	end;

	end.
