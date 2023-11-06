	const	strlmax=80;

	type	strbuf = packed array [1..strlmax] of char;
		tstring		= record
				    len : integer;
				    str : strbuf;
				  end;

	function strlen_(var str : tstring) : integer; external;

	procedure getstr(var str : tstring;
                         vstr : packed array [low..high : integer] of char);
                         external;

	procedure copystr(var outstr, instr : tstring); external;

	procedure concat(var str2 : tstring;
                         str1 : packed array [low..high : integer] of char);
                         external;

	procedure getint(var str : tstring; i, blanks : integer); external;

	procedure getreal(var str : tstring; x : double;
			  blanks, ndec : integer); external;

	procedure getreale(var str : tstring; x : double;
			   blanks, ndec : integer); external;

	function strind(var object: tstring;
                        pattern : packed array [low..high : integer] of char)
                        : integer; external;

	procedure writestr(var str : tstring; blanks : integer;
			   var outf : text); external;
